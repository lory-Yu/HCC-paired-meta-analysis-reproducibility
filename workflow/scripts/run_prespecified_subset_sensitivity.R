suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
  library(optparse)
})

option_list <- list(make_option("--project-root", dest = "project_root", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)

out_dir <- file.path(root, "results", "analysis", "sensitivity", "prespecified_subsets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta_input <- fread(file.path(root, "results", "analysis", "meta", "meta_input_long.tsv.gz"))
primary <- fread(file.path(root, "results", "analysis", "candidate_lock", "candidate_set_v1.0.tsv"))
primary[, gene_id := as.character(gene_id)]
meta_input[, gene_id := as.character(gene_id)]

scenarios <- list(
  affymetrix_only = c("GSE121248", "GSE84402"),
  illumina_only = c("GSE57957", "GSE76427"),
  adjacent_non_tumour_only = c("GSE121248", "GSE57957", "GSE76427")
)

fit_gene <- function(z) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(
      rma.uni(yi = z$log2FC, sei = z$SE, method = "REML", test = "knha"),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )
  list(value = value, warning = paste(unique(warnings), collapse = " | "))
}

all_results <- list()
scenario_summaries <- list()

for (scenario in names(scenarios)) {
  included <- scenarios[[scenario]]
  z_all <- meta_input[cohort %in% included]
  minimum_k <- length(included)
  genes <- z_all[, .N, by = gene_id][N == minimum_k, gene_id]
  message(scenario, ": fitting ", length(genes), " genes across ", minimum_k, " cohorts")

  rows <- vector("list", length(genes))
  for (i in seq_along(genes)) {
    gene <- genes[[i]]
    z <- z_all[gene_id == gene]
    fitted <- fit_gene(z)
    if (is.null(fitted$value)) next
    fit <- fitted$value
    pred <- suppressWarnings(predict(fit))
    positive_fraction <- mean(z$log2FC > 0)
    negative_fraction <- mean(z$log2FC < 0)
    rows[[i]] <- data.table(
      scenario = scenario,
      included_cohorts = paste(included, collapse = ";"),
      gene_id = gene,
      gene_symbol = { s <- z$gene_symbol[!is.na(z$gene_symbol) & nzchar(z$gene_symbol)]; if (length(s)) s[[1]] else NA_character_ },
      k = nrow(z),
      total_pairs = sum(z$n_pairs),
      pooled_log2FC = as.numeric(fit$b[1]),
      SE = fit$se,
      CI_low = fit$ci.lb,
      CI_high = fit$ci.ub,
      prediction_low = pred$pi.lb,
      prediction_high = pred$pi.ub,
      p_value = fit$pval,
      tau2 = fit$tau2,
      I2 = fit$I2,
      positive_fraction = positive_fraction,
      negative_fraction = negative_fraction,
      direction_consistency = max(positive_fraction, negative_fraction),
      majority_direction = ifelse(positive_fraction >= negative_fraction, "up", "down"),
      fit_warning = ifelse(nzchar(fitted$warning), fitted$warning, NA_character_)
    )
  }
  result <- rbindlist(rows, fill = TRUE)
  result[, FDR := p.adjust(p_value, method = "BH")]
  result[, formal_primary_rule_evaluable := FALSE]
  result[, formal_primary_rule_reason := paste0("subset has k=", minimum_k, "; frozen primary rule requires >=4 cohorts")]
  result <- result[order(FDR, -abs(pooled_log2FC))]
  fwrite(result, file.path(out_dir, paste0(scenario, "_all_genes.tsv.gz")), sep = "\t", na = "NA", compress = "gzip")

  candidate <- merge(
    primary[, .(gene_id, gene_symbol, primary_log2FC = pooled_log2FC, primary_FDR = FDR)],
    result[, .(gene_id, subset_log2FC = pooled_log2FC, subset_CI_low = CI_low,
               subset_CI_high = CI_high, subset_FDR = FDR, subset_I2 = I2,
               subset_direction_consistency = direction_consistency,
               subset_warning = fit_warning)],
    by = "gene_id", all.x = TRUE, suffixes = c("", "_subset")
  )
  candidate[, same_direction := sign(primary_log2FC) == sign(subset_log2FC)]
  fwrite(candidate, file.path(out_dir, paste0(scenario, "_candidate_comparison.tsv")), sep = "\t", na = "NA")

  estimable <- candidate[is.finite(subset_log2FC)]
  scenario_summaries[[scenario]] <- data.table(
    scenario = scenario,
    included_cohorts = paste(included, collapse = ";"),
    n_cohorts = minimum_k,
    all_genes_estimable = nrow(result),
    candidates_estimable = nrow(estimable),
    candidates_same_direction = sum(estimable$same_direction),
    same_direction_fraction = mean(estimable$same_direction),
    effect_pearson_r = cor(estimable$primary_log2FC, estimable$subset_log2FC, method = "pearson"),
    formal_primary_rule_evaluable = FALSE,
    interpretation = "direction-and-effect sensitivity only; not a new discovery screen"
  )
  all_results[[scenario]] <- result
}

# The two cohort-exclusion analyses are exactly represented by the corresponding
# full-universe leave-one-cohort-out runs, so summarize rather than recompute them.
loo <- fread(file.path(root, "results", "analysis", "sensitivity", "leave_one_cohort_out", "candidate_leave_one_out_long.tsv.gz"))
loo[, gene_id := as.character(gene_id)]
exclusion_summary <- loo[omitted_cohort %in% c("GSE57555", "GSE45114"), .(
  candidates_estimable = .N,
  candidates_same_direction = sum(same_direction_as_primary),
  same_direction_fraction = mean(same_direction_as_primary),
  effect_pearson_r = cor(primary_log2FC, pooled_log2FC, method = "pearson"),
  candidates_meeting_original_effect_direction_fdr_rule = sum(passes_primary_rule)
), by = omitted_cohort]
exclusion_summary[, scenario := fifelse(omitted_cohort == "GSE57555", "exclude_small_GSE57555", "exclude_two_colour_GSE45114")]
setcolorder(exclusion_summary, c("scenario", "omitted_cohort", setdiff(names(exclusion_summary), c("scenario", "omitted_cohort"))))
fwrite(exclusion_summary, file.path(out_dir, "cohort_exclusion_summary.tsv"), sep = "\t", na = "NA")

summary <- rbindlist(scenario_summaries)
fwrite(summary, file.path(out_dir, "subset_sensitivity_summary.tsv"), sep = "\t", na = "NA")

writeLines(c(
  "Takeaway: All remaining prespecified platform/control-definition sensitivity analyses were completed without redefining the frozen discovery set.",
  "",
  "# Prespecified subset sensitivity analyses",
  "",
  "The exclusion of GSE57555 and GSE45114 is already represented by the corresponding full-universe leave-one-cohort-out runs. Affymetrix-only, Illumina-only and adjacent-non-tumour-only analyses were fitted separately across their complete eligible gene universes.",
  "",
  "Because these subsets contain only two, two and three cohorts, respectively, they cannot satisfy the frozen primary requirement of at least four contributing cohorts. They are therefore direction-and-effect sensitivity analyses, not new discovery screens. No candidate lock, threshold or confidence stratum was changed.",
  "",
  paste0("- Affymetrix-only: ", summary[scenario == "affymetrix_only", candidates_same_direction], "/", summary[scenario == "affymetrix_only", candidates_estimable], " discovery candidates retained direction; Pearson r=", sprintf("%.3f", summary[scenario == "affymetrix_only", effect_pearson_r]), "."),
  paste0("- Illumina-only: ", summary[scenario == "illumina_only", candidates_same_direction], "/", summary[scenario == "illumina_only", candidates_estimable], " retained direction; Pearson r=", sprintf("%.3f", summary[scenario == "illumina_only", effect_pearson_r]), "."),
  paste0("- Adjacent-non-tumour-only: ", summary[scenario == "adjacent_non_tumour_only", candidates_same_direction], "/", summary[scenario == "adjacent_non_tumour_only", candidates_estimable], " retained direction; Pearson r=", sprintf("%.3f", summary[scenario == "adjacent_non_tumour_only", effect_pearson_r]), "."),
  "",
  "Full-universe results, candidate comparisons and the two cohort-exclusion summaries are stored in this directory."
), file.path(out_dir, "subset_sensitivity_review.md"), useBytes = TRUE)

