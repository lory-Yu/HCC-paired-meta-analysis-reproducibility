suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(metafor)
  library(optparse)
})

option_list <- list(make_option("--project-root", dest = "project_root", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)

out_root <- file.path(root, "results", "analysis", "sensitivity")
loo_root <- file.path(out_root, "leave_one_cohort_out")
full_root <- file.path(out_root, "all_samples")
for (path in c(out_root, loo_root, full_root)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

existing_artifacts <- c(
  file.path(out_root, "fixed_effect_all_genes.tsv.gz"),
  file.path(out_root, "fixed_effect_candidate_comparison.tsv"),
  file.path(loo_root, "leave_one_cohort_out_all_genes.tsv.gz"),
  file.path(loo_root, "leave_one_cohort_out_failures.tsv"),
  file.path(loo_root, "candidate_leave_one_out_long.tsv.gz"),
  file.path(loo_root, "candidate_leave_one_out_summary.tsv"),
  file.path(full_root, "all_samples_cohort_summary.tsv"),
  file.path(full_root, "all_samples_meta_all_genes.tsv.gz"),
  file.path(full_root, "all_samples_candidate_comparison.tsv"),
  file.path(out_root, "sensitivity_summary.tsv")
)
recompute <- identical(Sys.getenv("RECOMPUTE_SENSITIVITY", "0"), "1")
if (!recompute && all(file.exists(existing_artifacts)) && all(file.info(existing_artifacts)$size > 0)) {
  message("Reusing existing full-universe sensitivity artifacts; set RECOMPUTE_SENSITIVITY=1 to recompute.")
  quit(save = "no", status = 0)
}

cohorts <- c("GSE121248", "GSE45114", "GSE57555", "GSE57957", "GSE76427", "GSE84402")
primary_meta <- fread(file.path(root, "results", "analysis", "meta", "meta_all_genes.tsv.gz"))
meta_input <- fread(file.path(root, "results", "analysis", "meta", "meta_input_long.tsv.gz"))
candidates <- fread(file.path(root, "results", "analysis", "candidate_lock", "candidate_set_v1.0.tsv"))
candidate_ids <- as.character(candidates$gene_id)
manifest <- fread(file.path(root, "00_protocol", "sample_manifest_v1.0.tsv"), na.strings = "")
manifest <- manifest[status == "include"]
exclusions <- fread(file.path(root, "results", "qc", "qc_exclusions.tsv"), na.strings = "")
excluded_gsm <- if (nrow(exclusions)) as.character(exclusions$gsm) else character()

fit_reml <- function(z) {
  fit_once <- function(control = NULL) {
    warnings <- character()
    error <- NA_character_
    value <- tryCatch(
      withCallingHandlers({
        args <- list(yi = z$log2FC, sei = z$SE, method = "REML", test = "knha")
        if (!is.null(control)) args$control <- control
        fit <- do.call(rma.uni, args)
        pred <- predict(fit)
        list(fit = fit, pred = pred)
      }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) {
        error <<- conditionMessage(e)
        NULL
      }
    )
    list(value = value, warnings = unique(warnings), error = error)
  }
  initial <- fit_once()
  chosen <- initial
  optimizer <- "default"
  retried <- FALSE
  if (is.null(initial$value) || length(initial$warnings)) {
    retried <- TRUE
    retry <- fit_once(list(optimizer = "optim", optmethod = "Nelder-Mead", maxiter = 1000))
    if (!is.null(retry$value) && (is.null(initial$value) || !length(retry$warnings))) {
      chosen <- retry
      optimizer <- "optim_Nelder-Mead_retry"
    }
  }
  list(value = chosen$value, optimizer = optimizer, retried = retried,
       warning = if (length(chosen$warnings)) paste(unique(chosen$warnings), collapse = " | ") else NA_character_,
       error = chosen$error)
}

meta_row <- function(z, fit_result, extra = list()) {
  if (is.null(fit_result$value)) return(NULL)
  fit <- fit_result$value$fit
  pred <- fit_result$value$pred
  positive_fraction <- mean(z$log2FC > 0)
  negative_fraction <- mean(z$log2FC < 0)
  base <- data.table(
    gene_id = as.character(z$gene_id[[1]]),
    gene_symbol = { s <- z$gene_symbol[!is.na(z$gene_symbol) & nzchar(z$gene_symbol)]; if (length(s)) s[[1]] else NA_character_ },
    k = nrow(z), total_samples_or_pairs = sum(z$n_units),
    optimizer_used = fit_result$optimizer, retried = fit_result$retried,
    pooled_log2FC = as.numeric(fit$b[1]), SE = fit$se,
    CI_low = fit$ci.lb, CI_high = fit$ci.ub,
    prediction_low = pred$pi.lb, prediction_high = pred$pi.ub,
    p_value = fit$pval, tau2 = fit$tau2, I2 = fit$I2,
    Q = fit$QE, Q_p_value = fit$QEp,
    positive_fraction = positive_fraction, negative_fraction = negative_fraction,
    direction_consistency = max(positive_fraction, negative_fraction),
    majority_direction = ifelse(positive_fraction >= negative_fraction, "up", "down"),
    fit_warning = gsub("[\r\n]+", " ", fit_result$warning)
  )
  if (length(extra)) for (nm in names(extra)) base[, (nm) := extra[[nm]]]
  base
}

message("Fixed-effect comparison for the primary meta-analysis universe")
fixed <- meta_input[, {
  w <- 1 / SE^2
  estimate <- sum(w * log2FC) / sum(w)
  standard_error <- sqrt(1 / sum(w))
  q <- sum(w * (log2FC - estimate)^2)
  q_df <- .N - 1L
  positive_fraction <- mean(log2FC > 0)
  negative_fraction <- mean(log2FC < 0)
  list(
    gene_symbol = { s <- gene_symbol[!is.na(gene_symbol) & nzchar(gene_symbol)]; if (length(s)) s[[1]] else NA_character_ },
    k = .N, total_pairs = sum(n_pairs), pooled_log2FC = estimate, SE = standard_error,
    CI_low = estimate - qnorm(0.975) * standard_error,
    CI_high = estimate + qnorm(0.975) * standard_error,
    p_value = 2 * pnorm(-abs(estimate / standard_error)),
    Q = q, Q_p_value = pchisq(q, df = q_df, lower.tail = FALSE),
    I2 = if (q > 0) max(0, (q - q_df) / q) * 100 else 0,
    positive_fraction = positive_fraction, negative_fraction = negative_fraction,
    direction_consistency = max(positive_fraction, negative_fraction),
    majority_direction = ifelse(positive_fraction >= negative_fraction, "up", "down")
  )
}, by = gene_id][k >= 4L]
fixed[, FDR := p.adjust(p_value, method = "BH")]
fixed[, passes_primary_rule := FDR < 0.05 & abs(pooled_log2FC) >= log2(1.5) & direction_consistency >= 0.75]
fixed <- fixed[order(FDR, -abs(pooled_log2FC))]
fwrite(fixed, file.path(out_root, "fixed_effect_all_genes.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")

fixed_candidates <- merge(
  candidates[, .(gene_id = as.character(gene_id), gene_symbol,
                 random_log2FC = pooled_log2FC, random_FDR = FDR, random_I2 = I2)],
  fixed[, .(gene_id = as.character(gene_id), fixed_log2FC = pooled_log2FC,
            fixed_CI_low = CI_low, fixed_CI_high = CI_high, fixed_FDR = FDR,
            fixed_I2 = I2, fixed_passes_primary_rule = passes_primary_rule)],
  by = "gene_id", all.x = TRUE
)
fixed_candidates[, same_direction := sign(random_log2FC) == sign(fixed_log2FC)]
fwrite(fixed_candidates, file.path(out_root, "fixed_effect_candidate_comparison.tsv"), sep = "\t", na = "NA")

message("Leave-one-cohort-out REML/Hartung-Knapp analysis")
loo_rows <- list()
loo_failures <- list()
row_index <- 0L
failure_index <- 0L
primary_gene_ids <- as.character(primary_meta$gene_id)
for (omitted in cohorts) {
  z_omitted <- meta_input[cohort != omitted]
  eligible <- z_omitted[, .N, by = gene_id][N >= 3L, as.character(gene_id)]
  eligible <- intersect(primary_gene_ids, eligible)
  message("  omitting ", omitted, ": ", length(eligible), " genes")
  for (i in seq_along(eligible)) {
    gene <- eligible[[i]]
    z <- z_omitted[as.character(gene_id) == gene,
                   .(cohort, gene_id = as.character(gene_id), gene_symbol, log2FC, SE, n_units = n_pairs)]
    fitted <- fit_reml(z)
    if (is.null(fitted$value)) {
      failure_index <- failure_index + 1L
      loo_failures[[failure_index]] <- data.table(omitted_cohort = omitted, gene_id = gene,
                                                   error = gsub("[\r\n]+", " ", fitted$error))
    } else {
      row_index <- row_index + 1L
      loo_rows[[row_index]] <- meta_row(z, fitted, list(omitted_cohort = omitted))
    }
    if (i %% 4000L == 0L) message("    genes completed: ", i)
  }
}
loo <- rbindlist(loo_rows, fill = TRUE)
loo[, FDR := p.adjust(p_value, method = "BH"), by = omitted_cohort]
loo[, passes_primary_rule := FDR < 0.05 & abs(pooled_log2FC) >= log2(1.5) & direction_consistency >= 0.75]
setcolorder(loo, c("omitted_cohort", setdiff(names(loo), "omitted_cohort")))
loo <- loo[order(omitted_cohort, FDR, -abs(pooled_log2FC))]
fwrite(loo, file.path(loo_root, "leave_one_cohort_out_all_genes.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")
loo_failure_table <- if (length(loo_failures)) rbindlist(loo_failures) else
  data.table(omitted_cohort = character(), gene_id = character(), error = character())
fwrite(loo_failure_table, file.path(loo_root, "leave_one_cohort_out_failures.tsv"), sep = "\t", na = "NA")

candidate_primary <- candidates[, .(gene_id = as.character(gene_id), primary_log2FC = pooled_log2FC,
                                    primary_FDR = FDR, primary_I2 = I2)]
candidate_loo_long <- merge(loo[as.character(gene_id) %in% candidate_ids], candidate_primary,
                            by = "gene_id", all.x = TRUE)
candidate_loo_long[, same_direction_as_primary := sign(pooled_log2FC) == sign(primary_log2FC)]
fwrite(candidate_loo_long, file.path(loo_root, "candidate_leave_one_out_long.tsv.gz"),
       sep = "\t", na = "NA", compress = "gzip")
candidate_loo_summary <- candidate_loo_long[, .(
  scenarios_completed = .N,
  same_direction_all = .N == length(cohorts) && all(same_direction_as_primary),
  passes_primary_rule_all = .N == length(cohorts) && all(passes_primary_rule),
  min_log2FC = min(pooled_log2FC), max_log2FC = max(pooled_log2FC),
  max_FDR = max(FDR), min_I2 = min(I2), max_I2 = max(I2),
  any_prediction_interval_crosses_zero = any(prediction_low <= 0 & prediction_high >= 0)
), by = gene_id]
candidate_loo_summary <- merge(candidate_primary, candidate_loo_summary, by = "gene_id", all.x = TRUE)
fwrite(candidate_loo_summary, file.path(loo_root, "candidate_leave_one_out_summary.tsv"), sep = "\t", na = "NA")

message("All retained samples with patient blocking")
all_sample_effects <- list()
all_sample_summary <- list()
for (cohort in cohorts) {
  samples <- manifest[gse == cohort & !gsm %in% excluded_gsm]
  expr <- readRDS(file.path(root, "results", "analysis", "gene_expression", paste0(cohort, "_gene_expression.rds")))
  missing <- setdiff(samples$gsm, colnames(expr))
  if (length(missing)) stop(cohort, " all-sample matrix lacks GSMs: ", paste(missing, collapse = ","))
  y <- expr[, samples$gsm, drop = FALSE]
  group <- factor(samples$group, levels = c("control", "tumor"))
  design <- model.matrix(~group)
  block <- factor(samples$patient_id)
  correlation_fit <- duplicateCorrelation(y, design, block = block)
  if (!is.finite(correlation_fit$consensus)) stop(cohort, " duplicateCorrelation did not return a finite consensus")
  fit <- eBayes(lmFit(y, design, block = block, correlation = correlation_fit$consensus))
  coefficient <- "grouptumor"
  j <- match(coefficient, colnames(design))
  standard_error <- fit$stdev.unscaled[, j] * sqrt(fit$s2.post)
  patient_counts <- samples[, .(n = .N, groups = uniqueN(group)), by = patient_id]
  result <- data.table(
    cohort = cohort, gene_id = rownames(y),
    gene_symbol = meta_input[match(rownames(y), gene_id), gene_symbol],
    log2FC = fit$coefficients[, j], SE = standard_error,
    t = fit$t[, j], df_total = fit$df.total, p_value = fit$p.value[, j],
    FDR = p.adjust(fit$p.value[, j], method = "BH"), average_expression = fit$Amean,
    n_samples = ncol(y), n_patients = uniqueN(samples$patient_id),
    consensus_correlation = correlation_fit$consensus,
    contrast = "tumor-control"
  )
  result <- result[order(FDR, -abs(log2FC))]
  fwrite(result, file.path(full_root, paste0(cohort, "_all_samples_blocked_limma.tsv.gz")),
         sep = "\t", na = "NA", compress = "gzip")
  all_sample_effects[[cohort]] <- result
  all_sample_summary[[cohort]] <- data.table(
    cohort = cohort, n_samples = ncol(y), n_tumor = sum(group == "tumor"),
    n_control = sum(group == "control"), n_patients = uniqueN(samples$patient_id),
    n_repeated_patients = sum(patient_counts$n > 1L),
    n_complete_pair_patients = sum(patient_counts$n == 2L & patient_counts$groups == 2L),
    consensus_correlation = correlation_fit$consensus,
    genes_tested = nrow(result), cohort_FDR_lt_0_05 = sum(result$FDR < 0.05, na.rm = TRUE)
  )
}
all_sample_summary <- rbindlist(all_sample_summary)
fwrite(all_sample_summary, file.path(full_root, "all_samples_cohort_summary.tsv"), sep = "\t", na = "NA")

all_meta_input <- rbindlist(lapply(all_sample_effects, function(x) {
  x[is.finite(log2FC) & is.finite(SE) & SE > 0,
    .(cohort, gene_id = as.character(gene_id), gene_symbol, log2FC, SE, n_units = n_samples)]
}))
fwrite(all_meta_input, file.path(full_root, "all_samples_meta_input.tsv.gz"),
       sep = "\t", na = "NA", compress = "gzip")
all_genes <- all_meta_input[, .N, by = gene_id][N >= 4L, gene_id]
all_meta_rows <- vector("list", length(all_genes))
all_meta_failures <- list()
for (i in seq_along(all_genes)) {
  gene <- all_genes[[i]]
  z <- all_meta_input[gene_id == gene]
  fitted <- fit_reml(z)
  if (is.null(fitted$value)) {
    all_meta_failures[[length(all_meta_failures) + 1L]] <- data.table(
      gene_id = gene, error = gsub("[\r\n]+", " ", fitted$error))
  } else {
    all_meta_rows[[i]] <- meta_row(z, fitted)
  }
  if (i %% 4000L == 0L) message("  all-sample meta genes completed: ", i)
}
all_meta <- rbindlist(all_meta_rows, fill = TRUE)
all_meta[, FDR := p.adjust(p_value, method = "BH")]
all_meta[, passes_primary_rule := FDR < 0.05 & abs(pooled_log2FC) >= log2(1.5) & direction_consistency >= 0.75]
all_meta <- all_meta[order(FDR, -abs(pooled_log2FC))]
fwrite(all_meta, file.path(full_root, "all_samples_meta_all_genes.tsv.gz"), sep = "\t", na = "NA", compress = "gzip")
all_failure_table <- if (length(all_meta_failures)) rbindlist(all_meta_failures) else
  data.table(gene_id = character(), error = character())
fwrite(all_failure_table, file.path(full_root, "all_samples_meta_failures.tsv"), sep = "\t", na = "NA")

all_candidate_comparison <- merge(
  candidates[, .(gene_id = as.character(gene_id), gene_symbol,
                 primary_log2FC = pooled_log2FC, primary_FDR = FDR, primary_I2 = I2)],
  all_meta[, .(gene_id, all_samples_log2FC = pooled_log2FC,
               all_samples_CI_low = CI_low, all_samples_CI_high = CI_high,
               all_samples_FDR = FDR, all_samples_I2 = I2,
               all_samples_passes_primary_rule = passes_primary_rule)],
  by = "gene_id", all.x = TRUE
)
all_candidate_comparison[, same_direction := sign(primary_log2FC) == sign(all_samples_log2FC)]
fwrite(all_candidate_comparison, file.path(full_root, "all_samples_candidate_comparison.tsv"), sep = "\t", na = "NA")

summary <- data.table(
  metric = c(
    "locked_candidates", "leave_one_out_models", "leave_one_out_failures",
    "candidates_same_direction_all_leave_one_out", "candidates_pass_rule_all_leave_one_out",
    "fixed_effect_genes", "fixed_effect_rule_hits", "candidates_same_direction_fixed_effect",
    "candidates_pass_rule_fixed_effect", "all_samples_retained", "all_samples_meta_genes",
    "all_samples_meta_failures", "candidates_same_direction_all_samples",
    "candidates_pass_rule_all_samples"
  ),
  value = c(
    nrow(candidates), nrow(loo), nrow(loo_failure_table),
    sum(candidate_loo_summary$same_direction_all, na.rm = TRUE),
    sum(candidate_loo_summary$passes_primary_rule_all, na.rm = TRUE),
    nrow(fixed), sum(fixed$passes_primary_rule), sum(fixed_candidates$same_direction, na.rm = TRUE),
    sum(fixed_candidates$fixed_passes_primary_rule, na.rm = TRUE), sum(all_sample_summary$n_samples),
    nrow(all_meta), nrow(all_failure_table), sum(all_candidate_comparison$same_direction, na.rm = TRUE),
    sum(all_candidate_comparison$all_samples_passes_primary_rule, na.rm = TRUE)
  )
)
fwrite(summary, file.path(out_root, "sensitivity_summary.tsv"), sep = "\t")

writeLines(c(
  "# Sensitivity analysis run",
  "",
  paste0("Locked discovery candidates: ", nrow(candidates), "."),
  paste0("Leave-one-cohort-out models completed: ", nrow(loo), "; failures: ", nrow(loo_failure_table), "."),
  paste0("All-sample blocked models used ", sum(all_sample_summary$n_samples), " QC-retained samples."),
  "Fixed-effect and all-sample results are sensitivity analyses and do not redefine the discovery candidate lock."
), file.path(out_root, "run_summary.md"), useBytes = TRUE)
