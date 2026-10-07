#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(Biobase)
  library(optparse)
})

option_list <- list(
  make_option("--project-root", dest = "project_root", type = "character", help = "Project root")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)
qc_root <- file.path(root, "results", "qc")
dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)
manifest <- fread(file.path(root, "00_protocol", "sample_manifest_v1.0.tsv"), na.strings = c(""))
manifest <- manifest[status == "include"]

extract_gsm <- function(x) {
  names_only <- basename(x)
  ok <- grepl("GSM[0-9]+", names_only, ignore.case = TRUE)
  out <- rep(NA_character_, length(names_only))
  out[ok] <- toupper(sub(".*?(GSM[0-9]+).*", "\\1", names_only[ok],
                           ignore.case = TRUE, perl = TRUE))
  out
}

robust_z <- function(x) {
  center <- median(x, na.rm = TRUE)
  spread <- mad(x, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
  (x - center) / spread
}

order_by_manifest <- function(expr, cohort_manifest) {
  if (anyDuplicated(colnames(expr))) stop("duplicate GSM columns in expression matrix")
  missing <- setdiff(cohort_manifest$gsm, colnames(expr))
  extra <- setdiff(colnames(expr), cohort_manifest$gsm)
  if (length(missing)) stop("missing GSM columns: ", paste(missing, collapse = ","))
  if (length(extra)) expr <- expr[, setdiff(colnames(expr), extra), drop = FALSE]
  expr[, cohort_manifest$gsm, drop = FALSE]
}

read_affy_cohort <- function(cohort, cohort_manifest) {
  suppressPackageStartupMessages(library(affy))
  raw_dir <- file.path(root, "data", "work", cohort, "raw")
  files <- list.files(raw_dir, pattern = "[.]CEL$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("no decompressed CEL files")
  gsms <- extract_gsm(files)
  files <- files[match(cohort_manifest$gsm, gsms)]
  if (anyNA(files)) stop("not all eligible GSMs have CEL files")
  raw <- ReadAffy(filenames = files)
  normalized <- rma(raw, verbose = FALSE)
  expr <- exprs(normalized)
  colnames(expr) <- cohort_manifest$gsm
  list(expr = expr, detected = NULL, distribution_flags = FALSE,
       method = "Affymetrix RMA from raw CEL", features = nrow(expr))
}

read_gpr_cohort <- function(cohort, cohort_manifest) {
  raw_dir <- file.path(root, "data", "work", cohort, "raw")
  files <- list.files(raw_dir, pattern = "[.]gpr$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("no decompressed GPR files")
  gsms <- extract_gsm(files)
  files <- files[match(cohort_manifest$gsm, gsms)]
  if (anyNA(files)) stop("not all eligible GSMs have GPR files")
  # These legacy GenePix headers contain ISO-8859-1/GBK scanner paths and, in
  # two files, quoted header metadata that base read.table interprets as an
  # unmatched quote. The tabular FEATURES blocks themselves are structurally
  # complete. Read only those blocks, verify row identity, and build an RGList.
  read_features <- function(path) {
    prefix <- readLines(path, n = 100L, encoding = "latin1", warn = FALSE)
    header_line <- grep('^"?Block"?\\t"?Column"?\\t"?Row"?\\t"?Name"?', prefix)[1]
    if (is.na(header_line)) stop("GenePix FEATURES header not found in ", basename(path))
    fread(path, skip = header_line - 1L, encoding = "Latin-1",
          check.names = FALSE, data.table = TRUE)
  }
  tabs <- lapply(files, read_features)
  required <- c("Block", "Column", "Row", "ID", "F635 Median", "F532 Median",
                "B635 Median", "B532 Median")
  if (any(vapply(tabs, function(z) !all(required %in% names(z)), logical(1)))) {
    stop("one or more GPR feature blocks lack required GenePix columns")
  }
  # The later printing batch renames two control-spot IDs; physical feature
  # coordinates remain identical across all arrays, so align on those stable
  # coordinates rather than on control labels.
  keys <- lapply(tabs, function(z) paste(z$Block, z$Column, z$Row, sep = ":"))
  if (any(vapply(keys[-1], function(z) !identical(z, keys[[1]]), logical(1)))) {
    stop("GPR feature rows are not aligned across eligible samples")
  }
  make_matrix <- function(column) do.call(cbind, lapply(tabs, function(z) as.numeric(z[[column]])))
  rg <- list(
    R = make_matrix("F635 Median"), G = make_matrix("F532 Median"),
    Rb = make_matrix("B635 Median"), Gb = make_matrix("B532 Median"),
    genes = data.frame(ID = tabs[[1]]$ID, Name = tabs[[1]]$Name,
                       stringsAsFactors = FALSE)
  )
  class(rg) <- "RGList"
  rg <- backgroundCorrect(rg, method = "normexp", offset = 50)
  ma <- normalizeWithinArrays(rg, method = "loess")
  expr <- ma$M
  colnames(expr) <- cohort_manifest$gsm
  list(expr = expr, detected = NULL, distribution_flags = TRUE,
       method = "GenePix normexp background plus within-array loess M-values", features = nrow(expr))
}

read_agilent_cohort <- function(cohort, cohort_manifest) {
  raw_dir <- file.path(root, "data", "work", cohort, "raw")
  files <- list.files(raw_dir, pattern = "^GSM[0-9]+.*[.]txt$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  gsms <- extract_gsm(files)
  files <- files[match(cohort_manifest$gsm, gsms)]
  if (anyNA(files)) stop("not all eligible GSMs have Agilent TXT files")
  ag <- read.maimages(files, source = "agilent", green.only = TRUE)
  ag <- backgroundCorrect(ag, method = "normexp", offset = 16)
  expr <- normalizeBetweenArrays(log2(pmax(ag$E, 0) + 1), method = "quantile")
  colnames(expr) <- cohort_manifest$gsm
  list(expr = expr, detected = NULL, distribution_flags = FALSE,
       method = "Agilent normexp background, log2, between-array quantile normalization", features = nrow(expr))
}

read_illumina_cohort <- function(cohort, cohort_manifest) {
  path <- file.path(root, "data", "downloads", cohort, paste0(cohort, "_non-normalized.txt.gz"))
  if (!file.exists(path)) stop("non-normalized Illumina series file missing")
  tab <- fread(cmd = paste("gzip -cd", shQuote(path)), check.names = FALSE)
  signal_columns <- grep("[.]AVG_Signal$", names(tab), value = TRUE)
  if (!length(signal_columns)) stop("AVG_Signal columns not found")
  labels <- sub("[.]AVG_Signal$", "", signal_columns)
  if (cohort == "GSE57957") {
    expected_labels <- sub("^HCC_([0-9]+)([TN])$", "\\2\\1", cohort_manifest$source_name)
  } else if (cohort == "GSE76427") {
    number <- sub("^HCC", "", cohort_manifest$patient_id)
    expected_labels <- ifelse(cohort_manifest$group == "tumor", paste0("PT", number), paste0("ANTT", number))
  } else {
    stop("unsupported Illumina cohort")
  }
  index <- match(expected_labels, labels)
  if (anyNA(index)) {
    stop("Illumina matrix labels missing: ", paste(expected_labels[is.na(index)], collapse = ","))
  }
  signal <- as.matrix(tab[, ..signal_columns])
  storage.mode(signal) <- "double"
  signal <- signal[, index, drop = FALSE]
  colnames(signal) <- cohort_manifest$gsm
  det_columns <- paste0(labels, ".Detection Pval")
  detected <- NULL
  if (all(det_columns %in% names(tab))) {
    det <- as.matrix(tab[, ..det_columns])
    storage.mode(det) <- "double"
    det <- det[, index, drop = FALSE]
    colnames(det) <- cohort_manifest$gsm
    detected <- colMeans(det < 0.05, na.rm = TRUE)
  }
  expr <- normalizeBetweenArrays(log2(pmax(signal, 0) + 1), method = "quantile")
  list(expr = expr, detected = detected, distribution_flags = FALSE,
       method = "Illumina probe-summary log2 plus between-array quantile normalization for QC", features = nrow(expr))
}

compute_qc <- function(expr, cohort_manifest, detected = NULL, distribution_flags = TRUE) {
  keep <- rowSums(is.finite(expr)) >= ceiling(0.8 * ncol(expr))
  expr <- expr[keep, , drop = FALSE]
  sample_median <- apply(expr, 2, median, na.rm = TRUE)
  sample_iqr <- apply(expr, 2, IQR, na.rm = TRUE)
  missing_fraction <- colMeans(!is.finite(expr))
  row_variance <- apply(expr, 1, var, na.rm = TRUE)
  top <- order(row_variance, decreasing = TRUE)[seq_len(min(5000L, length(row_variance)))]
  correlations <- cor(expr[top, , drop = FALSE], use = "pairwise.complete.obs", method = "spearman")
  diag(correlations) <- NA_real_
  median_correlation <- apply(correlations, 2, median, na.rm = TRUE)
  pc <- prcomp(t(expr[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pc1 <- pc$x[, 1]
  pc2 <- if (ncol(pc$x) >= 2) pc$x[, 2] else rep(0, ncol(expr))

  z_median <- robust_z(sample_median)
  z_iqr <- robust_z(sample_iqr)
  z_corr <- robust_z(median_correlation)
  z_detected <- if (is.null(detected)) rep(NA_real_, ncol(expr)) else robust_z(detected)
  # Quantile-normalized distributions are intentionally equalized; tiny
  # residual rank/tie differences must not become exclusion evidence.
  flag_median <- distribution_flags & abs(z_median) > 3.5
  flag_iqr <- distribution_flags & abs(z_iqr) > 3.5
  flag_missing <- missing_fraction > 0.05
  flag_correlation <- z_corr < -3.5
  flag_detection <- if (is.null(detected)) rep(FALSE, ncol(expr)) else z_detected < -3.5
  flag_count <- flag_median + flag_iqr + flag_missing + flag_correlation + flag_detection

  data.table(
    gsm = colnames(expr),
    group = cohort_manifest$group[match(colnames(expr), cohort_manifest$gsm)],
    patient_id = cohort_manifest$patient_id[match(colnames(expr), cohort_manifest$gsm)],
    complete_pair = cohort_manifest$complete_pair[match(colnames(expr), cohort_manifest$gsm)],
    median = sample_median,
    iqr = sample_iqr,
    missing_fraction = missing_fraction,
    detected_fraction = if (is.null(detected)) NA_real_ else detected[colnames(expr)],
    median_spearman_correlation = median_correlation,
    pc1 = pc1,
    pc2 = pc2,
    robust_z_median = z_median,
    robust_z_iqr = z_iqr,
    robust_z_correlation = z_corr,
    robust_z_detected_fraction = z_detected,
    distribution_flags_enabled = distribution_flags,
    flag_median = flag_median,
    flag_iqr = flag_iqr,
    flag_missing = flag_missing,
    flag_correlation = flag_correlation,
    flag_detection = flag_detection,
    independent_qc_flags = flag_count,
    qc_decision = ifelse(flag_count >= 2, "exclude_by_prespecified_rule", ifelse(flag_count == 1, "review_single_flag", "retain"))
  )
}

plot_qc <- function(metrics, expr, cohort, out_dir) {
  cols <- ifelse(metrics$group == "tumor", "#B33A3A", "#2878A5")
  png(file.path(out_dir, "PCA.png"), width = 1800, height = 1500, res = 220)
  plot(metrics$pc1, metrics$pc2, col = cols, pch = 19,
       xlab = "PC1", ylab = "PC2", main = paste(cohort, "QC PCA"))
  legend("topright", legend = c("tumor", "control"), col = c("#B33A3A", "#2878A5"), pch = 19, bty = "n")
  flagged <- metrics$qc_decision != "retain"
  if (any(flagged)) text(metrics$pc1[flagged], metrics$pc2[flagged], labels = metrics$gsm[flagged], pos = 3, cex = 0.55)
  dev.off()

  png(file.path(out_dir, "sample_distributions.png"), width = 2200, height = 1500, res = 220)
  boxplot(expr, outline = FALSE, las = 2, cex.axis = 0.35, col = cols,
          ylab = "QC-normalized expression", main = paste(cohort, "sample distributions"))
  dev.off()
}

cohorts <- c("GSE121248", "GSE45114", "GSE57555", "GSE57957", "GSE76427", "GSE84402")
summary_rows <- list()
all_metrics <- list()
errors <- character()

for (cohort in cohorts) {
  message("QC: ", cohort)
  cohort_manifest <- manifest[gse == cohort]
  cohort_dir <- file.path(qc_root, cohort)
  dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)
  result <- tryCatch({
    x <- switch(
      cohort,
      GSE121248 = read_affy_cohort(cohort, cohort_manifest),
      GSE45114 = read_gpr_cohort(cohort, cohort_manifest),
      GSE57555 = read_agilent_cohort(cohort, cohort_manifest),
      GSE57957 = read_illumina_cohort(cohort, cohort_manifest),
      GSE76427 = read_illumina_cohort(cohort, cohort_manifest),
      GSE84402 = read_affy_cohort(cohort, cohort_manifest)
    )
    x$expr <- order_by_manifest(x$expr, cohort_manifest)
    metrics <- compute_qc(x$expr, cohort_manifest, x$detected, x$distribution_flags)
    metrics[, cohort := cohort]
    setcolorder(metrics, c("cohort", setdiff(names(metrics), "cohort")))
    fwrite(metrics, file.path(cohort_dir, "sample_qc_metrics.tsv"), sep = "\t", na = "NA")
    saveRDS(x$expr, file.path(cohort_dir, "qc_normalized_expression.rds"), compress = "xz")
    plot_qc(metrics, x$expr, cohort, cohort_dir)
    all_metrics[[cohort]] <- metrics
    summary_rows[[cohort]] <- data.table(
      cohort = cohort,
      status = "completed",
      eligible_samples = nrow(cohort_manifest),
      complete_pair_samples = sum(cohort_manifest$complete_pair == "yes"),
      features_after_platform_processing = x$features,
      no_qc_flag = sum(metrics$qc_decision == "retain"),
      retained_for_analysis = sum(metrics$qc_decision != "exclude_by_prespecified_rule"),
      single_flag_review = sum(metrics$qc_decision == "review_single_flag"),
      excluded_by_prespecified_rule = sum(metrics$qc_decision == "exclude_by_prespecified_rule"),
      method = x$method,
      error = ""
    )
    TRUE
  }, error = function(e) {
    errors <<- c(errors, paste0(cohort, ": ", conditionMessage(e)))
    summary_rows[[cohort]] <<- data.table(
      cohort = cohort, status = "failed", eligible_samples = nrow(cohort_manifest),
      complete_pair_samples = sum(cohort_manifest$complete_pair == "yes"),
      features_after_platform_processing = NA_integer_, no_qc_flag = NA_integer_,
      retained_for_analysis = NA_integer_,
      single_flag_review = NA_integer_, excluded_by_prespecified_rule = NA_integer_,
      method = "", error = conditionMessage(e)
    )
    FALSE
  })
}

summary_table <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_table, file.path(qc_root, "qc_summary.tsv"), sep = "\t", na = "NA")
metrics_table <- rbindlist(all_metrics, fill = TRUE)
if (nrow(metrics_table)) fwrite(metrics_table, file.path(qc_root, "all_sample_qc_metrics.tsv"), sep = "\t", na = "NA")

decisions <- metrics_table[qc_decision != "retain", .(
  cohort, gsm,
  decision = qc_decision,
  rule = ifelse(qc_decision == "exclude_by_prespecified_rule", "at least two independent QC flags", "one QC flag; review only"),
  evidence_1 = paste0("flag_count=", independent_qc_flags),
  evidence_2 = paste0("median_z=", round(robust_z_median, 3), ";iqr_z=", round(robust_z_iqr, 3),
                      ";corr_z=", round(robust_z_correlation, 3), ";detection_z=",
                      round(robust_z_detected_fraction, 3), ";missing=", round(missing_fraction, 5)),
  decision_time = format(Sys.time(), tz = "UTC", usetz = TRUE),
  notes = "QC decision uses only pre-specified technical metrics; PCA position is not an exclusion rule"
)]
fwrite(decisions, file.path(qc_root, "qc_decisions.tsv"), sep = "\t", na = "NA")
fwrite(decisions[decision == "exclude_by_prespecified_rule"],
       file.path(qc_root, "qc_exclusions.tsv"), sep = "\t", na = "NA")

review_lines <- c(
  "# Six-cohort GEO QC review",
  "",
  paste0("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "",
  "No differential-expression model, candidate-gene test, or outcome-guided exclusion was run.",
  "The independent biological unit remains the patient. PCA is descriptive and is not an exclusion criterion.",
  "",
  "## Cohort status",
  "",
  "| Cohort | Status | Eligible samples | Retained for analysis | Single-flag review (retained) | Excluded by frozen rule |",
  "|---|---|---:|---:|---:|---:|",
  apply(summary_table, 1, function(row) paste0("| ", row[["cohort"]], " | ", row[["status"]], " | ",
    row[["eligible_samples"]], " | ", row[["retained_for_analysis"]], " | ", row[["single_flag_review"]], " | ",
    row[["excluded_by_prespecified_rule"]], " |")),
  "",
  "## Interpretation boundary",
  "",
  "QC completion establishes file identity and technical usability only. It does not establish differential expression, biological replication, therapeutic efficacy, or direct compound-target evidence."
)
if (length(errors)) review_lines <- c(review_lines, "", "## Failures", "", paste0("- ", errors))
writeLines(review_lines, file.path(qc_root, "qc_review.md"), useBytes = TRUE)

if (length(errors)) stop("one or more cohorts failed QC: ", paste(errors, collapse = " | "))
