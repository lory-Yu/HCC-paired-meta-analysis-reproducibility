#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(metafor)
  library(ggplot2)
  library(pheatmap)
  library(optparse)
})

option_list <- list(make_option("--project-root", dest = "project_root", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)
out_root <- file.path(root, "results", "analysis")
ann_root <- file.path(out_root, "annotation")
gene_root <- file.path(out_root, "gene_expression")
de_root <- file.path(out_root, "differential")
meta_root <- file.path(out_root, "meta")
figure_root <- file.path(out_root, "figures")
for (path in c(out_root, ann_root, gene_root, de_root, meta_root, figure_root)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

cohorts <- c("GSE121248", "GSE45114", "GSE57555", "GSE57957", "GSE76427", "GSE84402")
platform_by_cohort <- c(
  GSE121248 = "GPL570", GSE45114 = "GPL5918", GSE57555 = "GPL16699",
  GSE57957 = "GPL10558", GSE76427 = "GPL10558", GSE84402 = "GPL570"
)

manifest <- fread(file.path(root, "00_protocol", "sample_manifest_v1.0.tsv"), na.strings = "")
manifest <- manifest[status == "include"]
exclusions <- fread(file.path(root, "results", "qc", "qc_exclusions.tsv"), na.strings = "")
excluded_gsm <- if (nrow(exclusions)) gsub('"', '', exclusions$gsm) else character()

single_numeric_id <- function(x) {
  x <- as.character(x)
  vapply(x, function(value) {
    if (is.na(value) || !nzchar(trimws(value))) return(NA_character_)
    ids <- unique(regmatches(value, gregexpr("[0-9]+", value))[[1]])
    ids <- ids[nzchar(ids)]
    if (length(ids) == 1L) ids else NA_character_
  }, character(1))
}

single_symbol <- function(x) {
  x <- as.character(x)
  vapply(x, function(value) {
    if (is.na(value) || !nzchar(trimws(value))) return(NA_character_)
    pieces <- unique(trimws(unlist(strsplit(value, "///|//|;|,|\\|"))))
    pieces <- pieces[nzchar(pieces) & pieces != "---"]
    if (length(pieces) == 1L) pieces else NA_character_
  }, character(1))
}

map_unique_orgdb <- function(keys, keytype, column) {
  keys <- as.character(keys)
  out <- rep(NA_character_, length(keys))
  ok <- !is.na(keys) & nzchar(keys)
  if (!any(ok)) return(out)
  mapped <- mapIds(
    org.Hs.eg.db, keys = unique(keys[ok]), keytype = keytype, column = column,
    multiVals = function(z) { z <- unique(z); if (length(z) == 1L) z else NA_character_ }
  )
  out[ok] <- unname(mapped[keys[ok]])
  out
}

read_geo_annotation <- function(path) {
  tab <- fread(cmd = paste("gzip -cd", shQuote(path)), skip = "ID\tGene title", sep = "\t",
               check.names = FALSE, fill = TRUE, data.table = TRUE)
  tab[, .(
    row_id = as.character(ID),
    gene_id = single_numeric_id(`Gene ID`),
    submitted_symbol = single_symbol(`Gene symbol`)
  )]
}

read_platform_table <- function(path, header) {
  tab <- fread(path, skip = header, sep = "\t", check.names = FALSE, fill = TRUE, data.table = TRUE)
  tab <- tab[grepl("^[0-9]+$", as.character(ID))]
  tab
}

message("Building frozen platform annotations")
annotations <- list()
annotations$GPL570 <- read_geo_annotation(file.path(root, "data", "downloads", "GPL570", "GPL570.annot.gz"))
annotations$GPL10558 <- read_geo_annotation(file.path(root, "data", "downloads", "GPL10558", "GPL10558.annot.gz"))

gpl5918 <- read_platform_table(
  file.path(root, "data", "downloads", "GPL5918", "GPL5918.txt"),
  "ID\tBlock\tRow\tColumn"
)
gpl5918[, submitted_symbol := single_symbol(HUGOname)]
gpl5918[, gene_id := map_unique_orgdb(submitted_symbol, "SYMBOL", "ENTREZID")]
annotations$GPL5918 <- gpl5918[, .(
  row_id = as.character(ID), gene_id = as.character(gene_id), submitted_symbol
)]

gpl16699 <- read_platform_table(
  file.path(root, "data", "downloads", "GPL16699", "GPL16699.txt"),
  "ID\tCOL\tROW\tNAME"
)
annotations$GPL16699 <- gpl16699[, .(
  row_id = as.character(ID), gene_id = single_numeric_id(LOCUSLINK_ID),
  submitted_symbol = single_symbol(GENE_SYMBOL)
)]

annotation_summary <- list()
for (platform in names(annotations)) {
  ann <- annotations[[platform]]
  ann[, gene_symbol := map_unique_orgdb(gene_id, "ENTREZID", "SYMBOL")]
  ann[, mapping_status := fifelse(is.na(gene_id), "excluded_nonunique_or_unmapped", "unique_entrez")]
  if (anyDuplicated(ann$row_id)) stop(platform, " platform row IDs are duplicated")
  fwrite(ann, file.path(ann_root, paste0(platform, "_probe_annotation.tsv.gz")),
         sep = "\t", na = "NA", compress = "gzip")
  annotation_summary[[platform]] <- data.table(
    platform = platform,
    platform_rows = nrow(ann),
    unique_entrez_probes = sum(!is.na(ann$gene_id)),
    unique_entrez_genes = uniqueN(ann$gene_id[!is.na(ann$gene_id)]),
    mapping_fraction = mean(!is.na(ann$gene_id))
  )
}
fwrite(rbindlist(annotation_summary), file.path(ann_root, "annotation_summary.tsv"), sep = "\t")

probe_ids_for_cohort <- function(cohort, expr) {
  platform <- platform_by_cohort[[cohort]]
  if (platform == "GPL570") return(rownames(expr))
  if (platform == "GPL10558") {
    path <- file.path(root, "data", "downloads", cohort, paste0(cohort, "_non-normalized.txt.gz"))
    ids <- fread(cmd = paste("gzip -cd", shQuote(path)), select = 1L, header = TRUE)[[1]]
    if (length(ids) != nrow(expr)) stop(cohort, " Illumina row count differs from QC matrix")
    return(as.character(ids))
  }
  ann <- annotations[[platform]]
  if (nrow(ann) != nrow(expr)) stop(cohort, " platform row count differs from QC matrix")
  ann$row_id
}

collapse_to_gene <- function(expr, annotation, probe_ids) {
  index <- match(probe_ids, annotation$row_id)
  gene_id <- annotation$gene_id[index]
  keep <- !is.na(gene_id) & nzchar(gene_id)
  if (!any(keep)) stop("no uniquely Entrez-mapped probes remain")
  summed <- rowsum(expr[keep, , drop = FALSE], group = gene_id[keep], reorder = TRUE, na.rm = TRUE)
  counts <- as.numeric(table(factor(gene_id[keep], levels = rownames(summed))))
  gene_expr <- summed / counts
  rownames(gene_expr) <- rownames(summed)
  gene_expr
}

run_paired_limma <- function(gene_expr, cohort_manifest, cohort) {
  primary <- cohort_manifest[complete_pair == "yes" & !gsm %in% excluded_gsm]
  valid_patients <- primary[, .(
    n = .N, groups = uniqueN(group), ntumor = sum(group == "tumor"), ncontrol = sum(group == "control")
  ), by = patient_id][n == 2L & groups == 2L & ntumor == 1L & ncontrol == 1L, patient_id]
  primary <- primary[patient_id %in% valid_patients]
  primary <- primary[order(patient_id, factor(group, levels = c("control", "tumor")))]
  if (nrow(primary) != 2L * length(valid_patients)) stop(cohort, " paired sample construction failed")
  missing <- setdiff(primary$gsm, colnames(gene_expr))
  if (length(missing)) stop(cohort, " gene matrix lacks primary GSMs: ", paste(missing, collapse = ","))
  y <- gene_expr[, primary$gsm, drop = FALSE]
  group <- factor(primary$group, levels = c("control", "tumor"))
  patient <- factor(primary$patient_id)
  design <- model.matrix(~ patient + group)
  coefficient <- "grouptumor"
  if (!coefficient %in% colnames(design) || qr(design)$rank != ncol(design)) {
    stop(cohort, " paired fixed-effect design is not full rank")
  }
  fit <- eBayes(lmFit(y, design))
  j <- match(coefficient, colnames(design))
  se <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
  result <- data.table(
    gene_id = rownames(y),
    gene_symbol = mapIds(org.Hs.eg.db, keys = rownames(y), keytype = "ENTREZID",
                         column = "SYMBOL", multiVals = "first"),
    log2FC = fit$coefficients[, j],
    SE = se,
    t = fit$t[, j],
    df_total = fit$df.total,
    p_value = fit$p.value[, j],
    FDR = p.adjust(fit$p.value[, j], method = "BH"),
    average_expression = fit$Amean,
    n_pairs = length(valid_patients),
    contrast = "tumor-control"
  )
  result <- result[order(FDR, -abs(log2FC))]
  list(table = result, n_pairs = length(valid_patients), samples = primary)
}

de_results <- list()
cohort_run_summary <- list()
for (cohort in cohorts) {
  message("Annotating and modelling ", cohort)
  platform <- platform_by_cohort[[cohort]]
  expr <- readRDS(file.path(root, "results", "qc", cohort, "qc_normalized_expression.rds"))
  probe_ids <- probe_ids_for_cohort(cohort, expr)
  gene_expr <- collapse_to_gene(expr, annotations[[platform]], probe_ids)
  saveRDS(gene_expr, file.path(gene_root, paste0(cohort, "_gene_expression.rds")), compress = "xz")
  fwrite(data.table(gene_id = rownames(gene_expr), gene_expr),
         file.path(gene_root, paste0(cohort, "_gene_expression.tsv.gz")),
         sep = "\t", na = "NA", compress = "gzip")
  model <- run_paired_limma(gene_expr, manifest[gse == cohort], cohort)
  model$table[, cohort := cohort]
  setcolorder(model$table, c("cohort", setdiff(names(model$table), "cohort")))
  fwrite(model$table, file.path(de_root, paste0(cohort, "_paired_limma.tsv.gz")),
         sep = "\t", na = "NA", compress = "gzip")
  fwrite(model$samples, file.path(de_root, paste0(cohort, "_primary_samples.tsv")), sep = "\t", na = "NA")
  de_results[[cohort]] <- model$table
  cohort_run_summary[[cohort]] <- data.table(
    cohort = cohort, platform = platform, n_pairs = model$n_pairs,
    genes_tested = nrow(model$table),
    cohort_FDR_lt_0_05 = sum(model$table$FDR < 0.05, na.rm = TRUE)
  )
}
cohort_summary <- rbindlist(cohort_run_summary)
fwrite(cohort_summary, file.path(de_root, "cohort_differential_summary.tsv"), sep = "\t")

meta_input <- rbindlist(lapply(names(de_results), function(cohort) {
  de_results[[cohort]][is.finite(log2FC) & is.finite(SE) & SE > 0,
                        .(cohort, gene_id, gene_symbol, log2FC, SE, n_pairs)]
}))
fwrite(meta_input, file.path(meta_root, "meta_input_long.tsv.gz"),
       sep = "\t", na = "NA", compress = "gzip")

genes <- meta_input[, .N, by = gene_id][N >= 4L, gene_id]
message("Running REML/Hartung-Knapp meta-analysis for ", length(genes), " genes")
meta_rows <- vector("list", length(genes))
meta_failures <- list()
meta_warnings <- list()
meta_retries <- list()

fit_meta_model <- function(z, control = NULL) {
  fit_warnings <- character()
  error_message <- NA_character_
  fitted <- tryCatch(
    withCallingHandlers({
      args <- list(yi = z$log2FC, sei = z$SE, method = "REML", test = "knha")
      if (!is.null(control)) args$control <- control
      fit <- do.call(rma.uni, args)
      list(fit = fit, prediction = predict(fit))
    }, warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  list(fitted = fitted, warnings = unique(fit_warnings), error = error_message)
}

for (i in seq_along(genes)) {
  gene <- genes[[i]]
  z <- meta_input[gene_id == gene]
  initial <- fit_meta_model(z)
  chosen <- initial
  optimizer_used <- "default"
  if (is.null(initial$fitted) || length(initial$warnings)) {
    retry <- fit_meta_model(
      z,
      control = list(optimizer = "optim", optmethod = "Nelder-Mead", maxiter = 1000)
    )
    meta_retries[[length(meta_retries) + 1L]] <- data.table(
      gene_id = gene,
      k = nrow(z),
      initial_error = initial$error,
      initial_warning = if (length(initial$warnings)) paste(initial$warnings, collapse = " | ") else NA_character_,
      retry_error = retry$error,
      retry_warning = if (length(retry$warnings)) paste(retry$warnings, collapse = " | ") else NA_character_,
      retry_succeeded = !is.null(retry$fitted)
    )
    if (!is.null(retry$fitted) && (is.null(initial$fitted) || !length(retry$warnings))) {
      chosen <- retry
      optimizer_used <- "optim_Nelder-Mead_retry"
    }
  }
  fitted <- chosen$fitted
  if (length(chosen$warnings)) {
    meta_warnings[[length(meta_warnings) + 1L]] <- data.table(
      gene_id = gene, optimizer_used = optimizer_used,
      warning = paste(chosen$warnings, collapse = " | ")
    )
  }
  if (is.null(fitted)) {
    meta_failures[[length(meta_failures) + 1L]] <- data.table(
      gene_id = gene, k = nrow(z), error = chosen$error
    )
  } else {
    fit <- fitted$fit
    prediction <- fitted$prediction
    positive_fraction <- mean(z$log2FC > 0)
    negative_fraction <- mean(z$log2FC < 0)
    meta_rows[[i]] <- data.table(
      gene_id = gene,
      gene_symbol = { s <- z$gene_symbol[!is.na(z$gene_symbol) & nzchar(z$gene_symbol)];
                      if (length(s)) s[[1]] else NA_character_ },
      k = nrow(z), total_pairs = sum(z$n_pairs),
      optimizer_used = optimizer_used,
      pooled_log2FC = as.numeric(fit$b[1]), SE = fit$se,
      CI_low = fit$ci.lb, CI_high = fit$ci.ub,
      prediction_low = prediction$pi.lb, prediction_high = prediction$pi.ub,
      p_value = fit$pval, tau2 = fit$tau2, I2 = fit$I2,
      Q = fit$QE, Q_p_value = fit$QEp,
      positive_fraction = positive_fraction,
      negative_fraction = negative_fraction,
      direction_consistency = max(positive_fraction, negative_fraction),
      majority_direction = ifelse(positive_fraction >= negative_fraction, "up", "down")
    )
  }
  if (i %% 2000L == 0L) message("  meta genes completed: ", i)
}
meta <- rbindlist(meta_rows, fill = TRUE)
failure_table <- if (length(meta_failures)) rbindlist(meta_failures, fill = TRUE) else
  data.table(gene_id = character(), k = integer(), error = character())
warning_table <- if (length(meta_warnings)) rbindlist(meta_warnings, fill = TRUE) else
  data.table(gene_id = character(), optimizer_used = character(), warning = character())
retry_table <- if (length(meta_retries)) rbindlist(meta_retries, fill = TRUE) else
  data.table(gene_id = character(), k = integer(), initial_error = character(),
             initial_warning = character(), retry_error = character(),
             retry_warning = character(), retry_succeeded = logical())
for (column in intersect(c("error", "warning"), names(failure_table))) {
  set(failure_table, j = column, value = gsub("[\r\n]+", " ", failure_table[[column]]))
}
for (column in intersect(c("warning"), names(warning_table))) {
  set(warning_table, j = column, value = gsub("[\r\n]+", " ", warning_table[[column]]))
}
for (column in intersect(c("initial_error", "initial_warning", "retry_error", "retry_warning"), names(retry_table))) {
  set(retry_table, j = column, value = gsub("[\r\n]+", " ", retry_table[[column]]))
}
fwrite(failure_table, file.path(meta_root, "meta_analysis_failures.tsv"), sep = "\t", na = "NA")
fwrite(warning_table, file.path(meta_root, "meta_analysis_warnings.tsv"), sep = "\t", na = "NA")
fwrite(retry_table, file.path(meta_root, "meta_analysis_retries.tsv"), sep = "\t", na = "NA")
meta[, FDR := p.adjust(p_value, method = "BH")]
meta[, passes_primary_rule := FDR < 0.05 & abs(pooled_log2FC) >= log2(1.5) & direction_consistency >= 0.75]
meta <- meta[order(FDR, -abs(pooled_log2FC))]
fwrite(meta, file.path(meta_root, "meta_all_genes.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")
fwrite(meta[passes_primary_rule == TRUE], file.path(meta_root, "meta_primary_hits.tsv"), sep = "\t", na = "NA")
fwrite(head(meta[passes_primary_rule == TRUE], 20L), file.path(meta_root, "meta_top20_primary_hits.tsv"), sep = "\t", na = "NA")

wide <- dcast(meta_input, gene_id + gene_symbol ~ cohort, value.var = "log2FC")
fwrite(wide, file.path(meta_root, "cohort_log2FC_wide.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")

plot_meta <- copy(meta)
plot_meta[, neg_log10_FDR := -log10(pmax(FDR, .Machine$double.xmin))]
plot_meta[, class := fifelse(passes_primary_rule & pooled_log2FC > 0, "Primary up",
                      fifelse(passes_primary_rule & pooled_log2FC < 0, "Primary down", "Other"))]
volcano <- ggplot(plot_meta, aes(pooled_log2FC, neg_log10_FDR, color = class)) +
  geom_point(size = 0.75, alpha = 0.65) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 2, linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.35) +
  scale_color_manual(values = c("Primary up" = "#B33A3A", "Primary down" = "#2878A5", "Other" = "#BDBDBD")) +
  labs(x = "Pooled log2 fold change (tumour - control)", y = "-log10 meta-analysis FDR", color = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(figure_root, "meta_volcano.png"), volcano, width = 7.2, height = 5.6, dpi = 320)

forest_ids <- head(meta[passes_primary_rule == TRUE, gene_id], 12L)
if (length(forest_ids)) {
  forest_cohort <- meta_input[gene_id %in% forest_ids,
    .(gene_id, gene_symbol, source = cohort, effect = log2FC,
      CI_low = log2FC - 1.96 * SE, CI_high = log2FC + 1.96 * SE,
      estimate_type = "Cohort")]
  forest_pooled <- meta[gene_id %in% forest_ids,
    .(gene_id, gene_symbol, source = "Pooled random effects", effect = pooled_log2FC,
      CI_low, CI_high, estimate_type = "Pooled")]
  forest_data <- rbind(forest_cohort, forest_pooled, fill = TRUE)
  forest_data[, gene_label := ifelse(is.na(gene_symbol) | gene_symbol == "", gene_id,
                                     paste0(gene_symbol, " (", gene_id, ")"))]
  forest_data[, source := factor(source, levels = rev(c(cohorts, "Pooled random effects")))]
  forest_plot <- ggplot(forest_data, aes(effect, source)) +
    geom_vline(xintercept = 0, color = "#777777", linewidth = 0.35, linetype = 2) +
    geom_segment(aes(x = CI_low, xend = CI_high, yend = source, color = estimate_type), linewidth = 0.55) +
    geom_point(aes(color = estimate_type, shape = estimate_type), size = 1.8) +
    facet_wrap(~gene_label, ncol = 3, scales = "free_x") +
    scale_color_manual(values = c(Cohort = "#4D4D4D", Pooled = "#B2182B")) +
    scale_shape_manual(values = c(Cohort = 16, Pooled = 18)) +
    labs(x = "log2 fold change (tumour - control), 95% CI", y = NULL,
         color = NULL, shape = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "top", strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  ggsave(file.path(figure_root, "top12_forest_plot.png"), forest_plot,
         width = 11, height = 10, dpi = 320)
}

top_ids <- head(meta[passes_primary_rule == TRUE, gene_id], 50L)
if (length(top_ids) >= 2L) {
  heat <- as.data.frame(wide[gene_id %in% top_ids])
  rownames(heat) <- ifelse(is.na(heat$gene_symbol) | heat$gene_symbol == "",
                          heat$gene_id, paste0(heat$gene_symbol, " (", heat$gene_id, ")"))
  heat <- as.matrix(heat[, cohorts, drop = FALSE])
  png(file.path(figure_root, "top50_cohort_effect_heatmap.png"), width = 1900, height = 2400, res = 240)
  pheatmap(heat, cluster_cols = FALSE, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
           main = "Top primary-rule genes: cohort log2 fold changes", border_color = NA,
           fontsize_row = 6.5, fontsize_col = 9)
  dev.off()
}

analysis_summary <- data.table(
  metric = c("cohorts_completed", "primary_pairs", "genes_meta_requested", "genes_meta_analyzed",
             "genes_meta_failed", "genes_retried", "genes_with_meta_warnings", "primary_rule_hits",
             "primary_rule_hits_up", "primary_rule_hits_down",
             "primary_hits_prediction_interval_excludes_zero", "primary_hits_I2_gt_50"),
  value = c(nrow(cohort_summary), sum(cohort_summary$n_pairs), length(genes), nrow(meta),
            nrow(failure_table), nrow(retry_table), nrow(warning_table), sum(meta$passes_primary_rule),
            sum(meta$passes_primary_rule & meta$majority_direction == "up"),
            sum(meta$passes_primary_rule & meta$majority_direction == "down"),
            sum(meta$passes_primary_rule & (meta$prediction_low > 0 | meta$prediction_high < 0)),
            sum(meta$passes_primary_rule & meta$I2 > 50))
)
fwrite(analysis_summary, file.path(out_root, "analysis_summary.tsv"), sep = "\t")

writeLines(c(
  "# Paired differential expression and meta-analysis run",
  "",
  paste0("Generated UTC: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "",
  paste0("Six paired limma models completed with ", sum(cohort_summary$n_pairs), " patient pairs."),
  paste0(length(genes), " genes had valid effects in at least four cohorts; ",
         nrow(meta), " models completed and ", nrow(failure_table), " failed."),
  paste0(nrow(retry_table), " genes required a prespecified numerical retry; ",
         nrow(warning_table), " final fits retained warnings."),
  paste0(sum(meta$passes_primary_rule), " genes met the frozen primary rule."),
  "",
  "The contrast is tumour minus control. Cohort models use patient fixed effects.",
  "Meta-analysis uses REML random effects and Hartung-Knapp inference.",
  "No cohort, sample, threshold, or effect direction was changed after inspecting differential results."
), file.path(out_root, "run_summary.md"), useBytes = TRUE)
