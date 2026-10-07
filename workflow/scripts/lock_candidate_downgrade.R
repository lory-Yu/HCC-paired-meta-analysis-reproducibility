suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

option_list <- list(make_option("--project-root", dest = "project_root", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)

lock_dir <- file.path(root, "results", "analysis", "candidate_lock")
sens_dir <- file.path(root, "results", "analysis", "sensitivity")
supp_dir <- file.path(sens_dir, "exploratory")
dir.create(lock_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

v1_path <- file.path(lock_dir, "candidate_set_v1.0.tsv")
v2_meta_path <- file.path(lock_dir, "candidate_set_v2_metadata.txt")
high_path <- file.path(lock_dir, "candidate_set_v2_high_confidence.tsv")
int_path <- file.path(lock_dir, "candidate_set_v2_intermediate.tsv")
if (file.exists(v2_meta_path) && file.exists(high_path) && file.exists(int_path)) {
  message("Candidate set v2 already locked; refusing to overwrite.")
  quit(save = "no", status = 0)
}

candidates <- fread(v1_path)
candidates[, gene_id := as.character(gene_id)]
if (nrow(candidates) != 446L) stop("Expected 446 locked discovery candidates")

loo <- fread(file.path(sens_dir, "leave_one_cohort_out", "candidate_leave_one_out_summary.tsv"))
loo[, gene_id := as.character(gene_id)]
loo_long <- fread(file.path(sens_dir, "leave_one_cohort_out", "candidate_leave_one_out_long.tsv.gz"))
loo_long[, gene_id := as.character(gene_id)]
fixed <- fread(file.path(sens_dir, "fixed_effect_candidate_comparison.tsv"))
fixed[, gene_id := as.character(gene_id)]
all_samples <- fread(file.path(sens_dir, "all_samples", "all_samples_candidate_comparison.tsv"))
all_samples[, gene_id := as.character(gene_id)]
warnings <- fread(file.path(root, "results", "analysis", "meta", "meta_analysis_warnings.tsv"))
warnings[, gene_id := as.character(gene_id)]
retries <- fread(file.path(root, "results", "analysis", "meta", "meta_analysis_retries.tsv"))
retries[, gene_id := as.character(gene_id)]
primary_failures <- fread(file.path(root, "results", "analysis", "meta", "meta_analysis_failures.tsv"))
primary_failures[, gene_id := as.character(gene_id)]
all_failures <- fread(file.path(sens_dir, "all_samples", "all_samples_meta_failures.tsv"))
all_failures[, gene_id := as.character(gene_id)]
loo_failures <- fread(file.path(sens_dir, "leave_one_cohort_out", "leave_one_cohort_out_failures.tsv"))
if (nrow(loo_failures)) loo_failures[, gene_id := as.character(gene_id)]

loo_warn <- loo_long[, .(
  loco_any_warning = any(!is.na(fit_warning) & nzchar(fit_warning)),
  loco_any_retry = any(retried %in% TRUE)
), by = gene_id]

audit <- merge(candidates, loo[, .(
  gene_id, loco_scenarios_completed = scenarios_completed,
  loco_same_direction_all = same_direction_all,
  loco_passes_primary_rule_all = passes_primary_rule_all,
  loco_min_log2FC = min_log2FC, loco_max_log2FC = max_log2FC,
  loco_max_FDR = max_FDR, loco_any_prediction_interval_crosses_zero = any_prediction_interval_crosses_zero
)], by = "gene_id", all.x = TRUE)
audit <- merge(audit, loo_warn, by = "gene_id", all.x = TRUE)
audit <- merge(audit, fixed[, .(
  gene_id, fixed_log2FC, fixed_FDR, fixed_I2,
  fixed_passes_primary_rule, fixed_same_direction = same_direction
)], by = "gene_id", all.x = TRUE)
audit <- merge(audit, all_samples[, .(
  gene_id, all_samples_log2FC, all_samples_FDR, all_samples_I2,
  all_samples_passes_primary_rule, all_samples_same_direction = same_direction
)], by = "gene_id", all.x = TRUE)

audit[, primary_prediction_interval_excludes_zero := prediction_low > 0 | prediction_high < 0]
audit[, primary_I2_lt_50 := I2 < 50]
audit[, primary_I2_ge_50 := I2 >= 50]
audit[, primary_I2_ge_75 := I2 >= 75]
audit[, unresolved_boundary_warning := gene_id %in% warnings$gene_id]
audit[, used_alternative_optimizer := optimizer_used != "default"]
audit[, primary_retried := gene_id %in% retries$gene_id]
audit[, primary_model_failed := gene_id %in% primary_failures$gene_id]
audit[, all_samples_model_failed := gene_id %in% all_failures$gene_id]
audit[, loco_model_failed := if (nrow(loo_failures)) gene_id %in% loo_failures$gene_id else FALSE]
audit[, any_model_failed := primary_model_failed | all_samples_model_failed | loco_model_failed]
audit[is.na(loco_any_warning), loco_any_warning := FALSE]
audit[is.na(loco_any_retry), loco_any_retry := FALSE]
audit[is.na(loco_same_direction_all), loco_same_direction_all := FALSE]
audit[is.na(fixed_same_direction), fixed_same_direction := FALSE]
audit[is.na(all_samples_same_direction), all_samples_same_direction := FALSE]
audit[is.na(loco_scenarios_completed), loco_scenarios_completed := 0L]

audit[, exploratory_reasons := ""]
add_reason <- function(flag, label) {
  audit[flag == TRUE, exploratory_reasons := ifelse(
    exploratory_reasons == "", label, paste(exploratory_reasons, label, sep = ";")
  )]
}
add_reason(!audit$loco_same_direction_all, "LOCO_direction_flip")
add_reason(audit$loco_scenarios_completed < 6L, "LOCO_scenarios_lt_6")
add_reason(audit$primary_I2_ge_75, "I2_ge_75")
add_reason(audit$unresolved_boundary_warning, "unresolved_boundary_warning")
add_reason(audit$used_alternative_optimizer, "alternative_optimizer")
add_reason(audit$any_model_failed, "model_failure")

audit[, high_confidence :=
        passes_primary_rule == TRUE &
        loco_same_direction_all == TRUE &
        loco_scenarios_completed == 6L &
        fixed_same_direction == TRUE &
        all_samples_same_direction == TRUE &
        primary_I2_lt_50 == TRUE &
        primary_prediction_interval_excludes_zero == TRUE &
        unresolved_boundary_warning == FALSE &
        used_alternative_optimizer == FALSE &
        any_model_failed == FALSE]
audit[, exploratory_only := nzchar(exploratory_reasons)]
audit[, intermediate_confidence :=
        passes_primary_rule == TRUE &
        loco_same_direction_all == TRUE &
        loco_scenarios_completed == 6L &
        fixed_same_direction == TRUE &
        all_samples_same_direction == TRUE &
        !high_confidence &
        !exploratory_only]
audit[, leftover := !high_confidence & !intermediate_confidence & !exploratory_only]
if (any(audit$leftover)) stop("Unclassified discovery candidates remain")
if (sum(audit$high_confidence) + sum(audit$intermediate_confidence) + sum(audit$exploratory_only) != 446L) {
  stop("Confidence strata do not sum to 446")
}

audit[, confidence_status := fifelse(
  high_confidence, "high_confidence",
  fifelse(intermediate_confidence, "intermediate_confidence", "exploratory_only")
)]
audit[, downgrade_is_not_negative := TRUE]
audit[, threshold_unchanged := TRUE]

setcolorder(audit, c(
  "gene_id", "gene_symbol", "confidence_status", "exploratory_reasons",
  "passes_primary_rule", "k", "total_pairs", "pooled_log2FC", "CI_low", "CI_high",
  "prediction_low", "prediction_high", "FDR", "I2", "tau2", "optimizer_used"
))
fwrite(audit, file.path(sens_dir, "candidate_downgrade_audit.tsv"), sep = "\t", na = "NA")

keep_cols <- names(candidates)
high <- audit[high_confidence == TRUE, ..keep_cols]
intermediate <- audit[intermediate_confidence == TRUE, ..keep_cols]
exploratory <- audit[exploratory_only == TRUE, ..keep_cols]
high <- high[order(FDR, -abs(pooled_log2FC))]
intermediate <- intermediate[order(FDR, -abs(pooled_log2FC))]
exploratory <- exploratory[order(FDR, -abs(pooled_log2FC))]

high_path <- file.path(lock_dir, "candidate_set_v2_high_confidence.tsv")
int_path <- file.path(lock_dir, "candidate_set_v2_intermediate.tsv")
expl_path <- file.path(supp_dir, "candidate_set_v2_exploratory_only.tsv")
fwrite(high, high_path, sep = "\t", na = "NA")
fwrite(intermediate, int_path, sep = "\t", na = "NA")
fwrite(exploratory, expl_path, sep = "\t", na = "NA")

file_sha <- function(path) strsplit(system2("sha256sum", path, stdout = TRUE), " +")[[1]][1]
high_sha <- file_sha(high_path)
int_sha <- file_sha(int_path)
expl_sha <- file_sha(expl_path)
v1_sha <- file_sha(v1_path)
audit_sha <- file_sha(file.path(sens_dir, "candidate_downgrade_audit.tsv"))

if (nrow(high) == 0L) {
  validation_note <- "high_confidence_count=0; do_not_expand_candidate_set; use_intermediate_as_labelled_secondary_validation_set"
} else {
  validation_note <- "primary_external_validation_set=high_confidence; secondary_labelled_set=intermediate_confidence; exploratory_supplement_only"
}

metadata <- c(
  "candidate_lock_version=v2",
  "source=results/analysis/candidate_lock/candidate_set_v1.0.tsv plus prespecified sensitivity strata",
  "discovery_candidate_count=446",
  paste0("v1_sha256=", v1_sha),
  paste0("high_confidence_count=", nrow(high)),
  paste0("intermediate_confidence_count=", nrow(intermediate)),
  paste0("exploratory_only_count=", nrow(exploratory)),
  paste0("high_confidence_sha256=", high_sha),
  paste0("intermediate_sha256=", int_sha),
  paste0("exploratory_sha256=", expl_sha),
  paste0("downgrade_audit_sha256=", audit_sha),
  "downgrade_does_not_delete_candidates=true",
  "downgrade_is_not_a_negative_result=true",
  "thresholds_not_retuned_after_downgrade=true",
  "generated_before_tcga_icgc_scrna_protein_chembl_lincs_depmap=true",
  validation_note
)
writeLines(metadata, file.path(lock_dir, "candidate_set_v2_metadata.txt"), useBytes = TRUE)

counts <- data.table(
  metric = c(
    "discovery_candidates",
    "loco_same_direction_all_six",
    "fixed_effect_same_direction",
    "all_438_samples_same_direction",
    "prediction_interval_excludes_zero",
    "I2_lt_50",
    "I2_ge_75",
    "unresolved_boundary_warning",
    "used_alternative_optimizer",
    "high_confidence",
    "intermediate_confidence",
    "exploratory_only"
  ),
  value = c(
    446L,
    sum(audit$loco_same_direction_all),
    sum(audit$fixed_same_direction),
    sum(audit$all_samples_same_direction),
    sum(audit$primary_prediction_interval_excludes_zero),
    sum(audit$primary_I2_lt_50),
    sum(audit$primary_I2_ge_75),
    sum(audit$unresolved_boundary_warning),
    sum(audit$used_alternative_optimizer),
    nrow(high),
    nrow(intermediate),
    nrow(exploratory)
  )
)
fwrite(counts, file.path(sens_dir, "confidence_stratum_counts.tsv"), sep = "\t")

message("Candidate set v2 locked: high=", nrow(high),
        " intermediate=", nrow(intermediate),
        " exploratory=", nrow(exploratory))
