#!/usr/bin/env Rscript
# Post-rejection multiplicity reconstruction for TCGA-LIHC and ICGC-LIRI-JP.
# Old validation outputs are read-only. All outputs are written to a new CBC directory.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(optparse)
})

option_list <- list(make_option("--project-root", dest = "project_root", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)
out <- file.path(root, "results", "validation", "cbc_external_validation_v1")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

high_path <- file.path(root, "results", "analysis", "candidate_lock", "candidate_set_v2_high_confidence.tsv")
int_path <- file.path(root, "results", "analysis", "candidate_lock", "candidate_set_v2_intermediate.tsv")
icgc_count_path <- file.path(root, "results", "validation", "rnaseq_normalized_v2", "icgc_liri_jp_paired_full_counts.tsv.gz")
tcga_count_path <- file.path(root, "data", "external", "tcga_lihc", "TCGA-LIHC.star_counts.tsv.gz")

sha_of <- function(path) strsplit(system2("sha256sum", path, stdout = TRUE), " +")[[1]][1]
stopifnot(
  identical(sha_of(high_path), "55110e282bf9d26619de70d6e8c1e5972d050766b034a6fa06009021ae9fd3cd"),
  identical(sha_of(int_path), "23812a47c8df31301aca8a9061ecd7da56dd41afa5f214a3fc6e4a585494aa59")
)

high <- fread(high_path)
inter <- fread(int_path)
high[, `:=`(gene_id = as.character(gene_id), candidate_order = .I,
            evidence_tier = "discovery_robust")]
inter[, `:=`(gene_id = as.character(gene_id), candidate_order = .I,
             evidence_tier = "intermediate")]
candidates <- rbind(
  high[, .(gene_id, gene_symbol, discovery_log2FC = pooled_log2FC,
           discovery_FDR = FDR, evidence_tier, candidate_order)],
  inter[, .(gene_id, gene_symbol, discovery_log2FC = pooled_log2FC,
            discovery_FDR = FDR, evidence_tier, candidate_order)]
)
stopifnot(nrow(high) == 173L, nrow(inter) == 119L, nrow(candidates) == 292L)
stopifnot(!anyDuplicated(high$gene_symbol), !anyDuplicated(inter$gene_symbol))

map_unique <- function(keys, keytype, column) {
  keys <- as.character(keys)
  out <- rep(NA_character_, length(keys))
  ok <- !is.na(keys) & nzchar(keys)
  if (!any(ok)) return(out)
  mapped <- mapIds(
    org.Hs.eg.db, keys = unique(keys[ok]), keytype = keytype, column = column,
    multiVals = function(z) {
      z <- unique(z[!is.na(z) & nzchar(z)])
      if (length(z) == 1L) z else NA_character_
    }
  )
  out[ok] <- unname(mapped[keys[ok]])
  out
}

fit_voom_all <- function(counts, pair_ids, group, key_name, cohort) {
  stopifnot(ncol(counts) == length(pair_ids), length(group) == length(pair_ids))
  stopifnot(!anyDuplicated(rownames(counts)))
  patient <- factor(pair_ids)
  group <- factor(group, levels = c("normal", "tumor"))
  design <- model.matrix(~ patient + group)

  dge0 <- DGEList(counts = counts)
  keep <- filterByExpr(dge0, group = group)
  dge <- dge0[keep, , keep.lib.sizes = FALSE]
  rm(dge0)
  dge <- calcNormFactors(dge, method = "TMM")
  lib_size_min <- min(dge$samples$lib.size)
  lib_size_median <- median(dge$samples$lib.size)
  lib_size_max <- max(dge$samples$lib.size)
  norm_factor_min <- min(dge$samples$norm.factors)
  norm_factor_median <- median(dge$samples$norm.factors)
  norm_factor_max <- max(dge$samples$norm.factors)
  v <- voom(dge, design, plot = FALSE)
  rm(dge)
  invisible(gc())
  feature_names <- rownames(v$E)
  mean_log_cpm <- rowMeans(v$E)
  fit <- eBayes(lmFit(v, design))
  rm(v)
  invisible(gc())
  idx <- match("grouptumor", colnames(fit$coefficients))
  if (is.na(idx)) stop("Missing group coefficient for ", cohort)

  # `key` is a reserved formal argument of data.table(); use a neutral
  # temporary column and rename it after construction.
  ans <- data.table(
    feature_key = feature_names,
    validation_log2FC = as.numeric(fit$coefficients[, idx]),
    SE = as.numeric(fit$stdev.unscaled[, idx] * sqrt(fit$s2.post)),
    nominal_P = as.numeric(fit$p.value[, idx]),
    df_total = as.numeric(fit$df.total),
    mean_logCPM = mean_log_cpm
  )
  ans[, genomewide_FDR := p.adjust(nominal_P, method = "BH")]
  ans[, validation_CI_low := validation_log2FC - qt(0.975, df_total) * SE]
  ans[, validation_CI_high := validation_log2FC + qt(0.975, df_total) * SE]
  setnames(ans, "feature_key", key_name)
  ans[, cohort := cohort]

  audit <- data.table(
    cohort = cohort,
    samples = ncol(counts),
    pairs = uniqueN(pair_ids),
    genes_before_filter = nrow(counts),
    genes_after_filter = nrow(ans),
    genomewide_BH_family_n = sum(is.finite(ans$nominal_P)),
    library_size_min = lib_size_min,
    library_size_median = lib_size_median,
    library_size_max = lib_size_max,
    norm_factor_min = norm_factor_min,
    norm_factor_median = norm_factor_median,
    norm_factor_max = norm_factor_max,
    residual_df_min = min(fit$df.residual),
    residual_df_max = max(fit$df.residual),
    method = "edgeR TMM + limma-voom; patient fixed effect; tumour-minus-normal",
    multiplicity = "BH across the cohort-specific complete filtered transcriptome"
  )
  list(full = ans, before_keys = rownames(counts), audit = audit)
}

candidate_results <- function(base, fit_obj, key_name, cohort) {
  fitted <- copy(fit_obj$full)
  model_cols <- c(
    key_name, "validation_log2FC", "SE", "nominal_P", "df_total", "mean_logCPM",
    "genomewide_FDR", "validation_CI_low", "validation_CI_high", "cohort"
  )
  fitted <- fitted[, ..model_cols]
  ans <- merge(base, fitted, by = key_name, all.x = TRUE, sort = FALSE)
  setorder(ans, evidence_tier, candidate_order)
  ans[, present_before_filter := as.character(get(key_name)) %in% fit_obj$before_keys]
  ans[, retained_after_filter := is.finite(nominal_P)]
  ans[, estimable := present_before_filter & retained_after_filter]
  ans[, same_direction := fifelse(
    estimable,
    sign(validation_log2FC) == sign(discovery_log2FC),
    NA
  )]
  ans[, candidate_family_FDR := NA_real_]
  for (tier in c("discovery_robust", "intermediate")) {
    idx <- ans$evidence_tier == tier & is.finite(ans$nominal_P)
    ans[idx, candidate_family_FDR := p.adjust(nominal_P, method = "BH")]
  }
  ans[, nominal_directional := estimable & same_direction == TRUE & nominal_P < 0.05]
  ans[, candidate_family_directional := estimable & same_direction == TRUE & candidate_family_FDR < 0.05]
  ans[, genomewide_directional := estimable & same_direction == TRUE & genomewide_FDR < 0.05]
  ans[, estimation_note := fifelse(
    present_before_filter == FALSE,
    paste0("absent_from_", cohort, "_count_matrix"),
    fifelse(retained_after_filter == FALSE, "filtered_low_expression_before_voom", "")
  )]
  ans[, `:=`(
    n_pairs = ifelse(cohort == "TCGA-LIHC", 50L, 199L),
    design = "paired_TMM_voom_patient_fixed_effect",
    post_rejection_analysis = TRUE
  )]
  ans
}

message("TCGA-LIHC: reconstructing integer counts from Xena log2(STAR count + 1)")
tcga_pairs <- fread(file.path(root, "results", "validation", "tcga_lihc", "paired_patients.tsv"))
tcga_samples <- c(tcga_pairs$normal_sample, tcga_pairs$tumor_sample)
tcga <- fread(tcga_count_path, select = c("Ensembl_ID", tcga_samples))
ens_versioned <- tcga$Ensembl_ID
ens <- sub("\\..*$", "", ens_versioned)
tcga_log <- as.matrix(tcga[, -1, with = FALSE])
storage.mode(tcga_log) <- "numeric"
tcga_counts0 <- round(pmax(2^tcga_log - 1, 0))
rownames(tcga_counts0) <- ens_versioned
entrez <- map_unique(ens, "ENSEMBL", "ENTREZID")
ok <- !is.na(entrez)
tcga_counts <- rowsum(tcga_counts0[ok, , drop = FALSE], group = entrez[ok], reorder = FALSE)
storage.mode(tcga_counts) <- "numeric"
tcga_fit <- fit_voom_all(
  tcga_counts,
  pair_ids = c(tcga_pairs$patient_id, tcga_pairs$patient_id),
  group = c(rep("normal", nrow(tcga_pairs)), rep("tumor", nrow(tcga_pairs))),
  key_name = "gene_id",
  cohort = "TCGA-LIHC"
)
tcga_fit$full[, gene_symbol := map_unique(gene_id, "ENTREZID", "SYMBOL")]
setcolorder(tcga_fit$full, c("cohort", "gene_id", "gene_symbol"))
tcga_candidates <- candidate_results(candidates, tcga_fit, "gene_id", "TCGA-LIHC")
tcga_audit <- copy(tcga_fit$audit)
# The full TCGA fit is already represented by the candidate table and will be
# written immediately below; release the large count/model objects before the
# ICGC fit so both complete-universe analyses can run in constrained workers.
fwrite(tcga_fit$full, file.path(out, "tcga_complete_filtered_transcriptome.tsv.gz"), sep = "\t", na = "NA")
rm(tcga_counts, tcga, tcga_fit)
invisible(gc())

message("ICGC-LIRI-JP: reading the frozen full paired-sample raw-count matrix")
icgc_pairs <- fread(file.path(root, "results", "validation", "icgc_liri_jp", "paired_donors.tsv"))
icgc <- fread(icgc_count_path, integer64 = "double")
symbols <- icgc$gene_symbol
icgc_counts <- as.matrix(icgc[, -1, with = FALSE])
storage.mode(icgc_counts) <- "numeric"
rownames(icgc_counts) <- symbols
expected_icgc <- c(icgc_pairs$normal_sample, icgc_pairs$tumor_sample)
stopifnot(identical(colnames(icgc_counts), expected_icgc))
icgc_fit <- fit_voom_all(
  icgc_counts,
  pair_ids = c(icgc_pairs$icgc_donor_id, icgc_pairs$icgc_donor_id),
  group = c(rep("normal", nrow(icgc_pairs)), rep("tumor", nrow(icgc_pairs))),
  key_name = "gene_symbol",
  cohort = "ICGC-LIRI-JP"
)
icgc_fit$full[, gene_id := map_unique(gene_symbol, "SYMBOL", "ENTREZID")]
setcolorder(icgc_fit$full, c("cohort", "gene_symbol", "gene_id"))
icgc_candidates <- candidate_results(candidates, icgc_fit, "gene_symbol", "ICGC-LIRI-JP")

fwrite(icgc_fit$full, file.path(out, "icgc_complete_filtered_transcriptome.tsv.gz"), sep = "\t", na = "NA")
fwrite(tcga_candidates, file.path(out, "tcga_locked_candidate_multiplicity.tsv"), sep = "\t", na = "NA")
fwrite(icgc_candidates, file.path(out, "icgc_locked_candidate_multiplicity.tsv"), sep = "\t", na = "NA")

cross <- merge(
  tcga_candidates[, .(
    gene_id, gene_symbol, evidence_tier, candidate_order, discovery_log2FC, discovery_FDR,
    tcga_estimable = estimable, tcga_log2FC = validation_log2FC,
    tcga_nominal_P = nominal_P, tcga_candidate_family_FDR = candidate_family_FDR,
    tcga_genomewide_FDR = genomewide_FDR, tcga_same_direction = same_direction,
    tcga_nominal_directional = nominal_directional,
    tcga_candidate_family_directional = candidate_family_directional,
    tcga_genomewide_directional = genomewide_directional
  )],
  icgc_candidates[, .(
    gene_symbol, evidence_tier, candidate_order,
    icgc_estimable = estimable, icgc_log2FC = validation_log2FC,
    icgc_nominal_P = nominal_P, icgc_candidate_family_FDR = candidate_family_FDR,
    icgc_genomewide_FDR = genomewide_FDR, icgc_same_direction = same_direction,
    icgc_nominal_directional = nominal_directional,
    icgc_candidate_family_directional = candidate_family_directional,
    icgc_genomewide_directional = genomewide_directional
  )],
  by = c("gene_symbol", "evidence_tier", "candidate_order"), all.x = TRUE, sort = FALSE
)
setorder(cross, evidence_tier, candidate_order)
cross[, jointly_estimable := tcga_estimable == TRUE & icgc_estimable == TRUE]
cross[, dual_same_direction := jointly_estimable & tcga_same_direction == TRUE & icgc_same_direction == TRUE]
cross[, dual_nominal_directional := jointly_estimable & tcga_nominal_directional == TRUE & icgc_nominal_directional == TRUE]
cross[, dual_candidate_family_directional := jointly_estimable & tcga_candidate_family_directional == TRUE & icgc_candidate_family_directional == TRUE]
cross[, dual_genomewide_directional := jointly_estimable & tcga_genomewide_directional == TRUE & icgc_genomewide_directional == TRUE]
fwrite(cross, file.path(out, "tcga_icgc_locked_candidate_cross_validation.tsv"), sep = "\t", na = "NA")

count_rows <- list()
add_count <- function(tier, scope, metric, numerator, denominator, denominator_definition) {
  count_rows[[length(count_rows) + 1L]] <<- data.table(
    evidence_tier = tier, scope = scope, metric = metric,
    numerator = as.integer(numerator), denominator = as.integer(denominator),
    proportion = ifelse(denominator > 0, numerator / denominator, NA_real_),
    denominator_definition = denominator_definition
  )
}

for (tier in c("discovery_robust", "intermediate")) {
  total <- candidates[evidence_tier == tier, .N]
  for (cohort_name in c("TCGA-LIHC", "ICGC-LIRI-JP")) {
    z <- if (cohort_name == "TCGA-LIHC") tcga_candidates[evidence_tier == tier] else icgc_candidates[evidence_tier == tier]
    est_n <- z[estimable == TRUE, .N]
    add_count(tier, cohort_name, "estimable", est_n, total, "locked tier total")
    add_count(tier, cohort_name, "same_direction", z[estimable == TRUE & same_direction == TRUE, .N], est_n, "estimable in this cohort")
    add_count(tier, cohort_name, "nominal_P_lt_0.05_and_same_direction", z[nominal_directional == TRUE, .N], est_n, "estimable in this cohort")
    add_count(tier, cohort_name, "candidate_family_FDR_lt_0.05_and_same_direction", z[candidate_family_directional == TRUE, .N], est_n, "estimable in this cohort")
    add_count(tier, cohort_name, "genomewide_FDR_lt_0.05_and_same_direction", z[genomewide_directional == TRUE, .N], est_n, "estimable in this cohort")
  }
  zc <- cross[evidence_tier == tier]
  joint_n <- zc[jointly_estimable == TRUE, .N]
  add_count(tier, "TCGA-and-ICGC", "jointly_estimable", joint_n, total, "locked tier total")
  add_count(tier, "TCGA-and-ICGC", "dual_same_direction", zc[dual_same_direction == TRUE, .N], joint_n, "jointly estimable in both cohorts")
  add_count(tier, "TCGA-and-ICGC", "dual_nominal_P_lt_0.05_and_same_direction", zc[dual_nominal_directional == TRUE, .N], joint_n, "jointly estimable in both cohorts")
  add_count(tier, "TCGA-and-ICGC", "dual_candidate_family_FDR_lt_0.05_and_same_direction", zc[dual_candidate_family_directional == TRUE, .N], joint_n, "jointly estimable in both cohorts")
  add_count(tier, "TCGA-and-ICGC", "dual_genomewide_FDR_lt_0.05_and_same_direction", zc[dual_genomewide_directional == TRUE, .N], joint_n, "jointly estimable in both cohorts")
}
counts <- rbindlist(count_rows)
fwrite(counts, file.path(out, "multiplicity_counts_with_denominators.tsv"), sep = "\t", na = "NA")

audit <- rbind(tcga_audit, icgc_fit$audit, fill = TRUE)
fwrite(audit, file.path(out, "normalization_and_multiplicity_audit.tsv"), sep = "\t", na = "NA")

manifest <- data.table(
  role = c("discovery_robust_lock", "intermediate_lock", "TCGA_source_counts", "ICGC_frozen_full_counts"),
  path = c(high_path, int_path, tcga_count_path, icgc_count_path),
  sha256 = vapply(c(high_path, int_path, tcga_count_path, icgc_count_path), sha_of, character(1))
)
fwrite(manifest, file.path(out, "locked_input_manifest.tsv"), sep = "\t")

# Avoid data.table non-standard evaluation ambiguity in narrative lookup.
lookup <- function(tier_value, scope_value, metric_value) {
  z <- counts[evidence_tier == tier_value & scope == scope_value & metric == metric_value]
  stopifnot(nrow(z) == 1L)
  paste0(z$numerator, "/", z$denominator)
}

review <- c(
  "# CBC 外部验证多重检验口径重构结果",
  "",
  "本分析由拒稿编辑意见触发，属于冻结后纠正分析。发现阶段锁定的 173 个基因在本报告中称为 discovery robust set；候选锁、样本和模型未改变。",
  "",
  "两个外部队列均在各自完整过滤后转录组上执行 TMM + limma-voom 患者固定效应模型。每个队列同时报告未经校正 P 值、锁定候选家族 BH-FDR 和全转录组 BH-FDR，三者不得互换。",
  "",
  paste0("- TCGA 完整过滤后检验宇宙：", audit[cohort == "TCGA-LIHC", genes_after_filter], " 个基因。"),
  paste0("- ICGC 完整过滤后检验宇宙：", audit[cohort == "ICGC-LIRI-JP", genes_after_filter], " 个基因。"),
  paste0("- discovery robust set 在 TCGA 可估计：", lookup("discovery_robust", "TCGA-LIHC", "estimable"), "；同方向：", lookup("discovery_robust", "TCGA-LIHC", "same_direction"), "；候选家族 FDR：", lookup("discovery_robust", "TCGA-LIHC", "candidate_family_FDR_lt_0.05_and_same_direction"), "；全转录组 FDR：", lookup("discovery_robust", "TCGA-LIHC", "genomewide_FDR_lt_0.05_and_same_direction"), "。"),
  paste0("- discovery robust set 在 ICGC 可估计：", lookup("discovery_robust", "ICGC-LIRI-JP", "estimable"), "；同方向：", lookup("discovery_robust", "ICGC-LIRI-JP", "same_direction"), "；候选家族 FDR：", lookup("discovery_robust", "ICGC-LIRI-JP", "candidate_family_FDR_lt_0.05_and_same_direction"), "；全转录组 FDR：", lookup("discovery_robust", "ICGC-LIRI-JP", "genomewide_FDR_lt_0.05_and_same_direction"), "。"),
  paste0("- 两库共同可估计：", lookup("discovery_robust", "TCGA-and-ICGC", "jointly_estimable"), "；双库同方向：", lookup("discovery_robust", "TCGA-and-ICGC", "dual_same_direction"), "；双库候选家族 FDR：", lookup("discovery_robust", "TCGA-and-ICGC", "dual_candidate_family_FDR_lt_0.05_and_same_direction"), "；双库全转录组 FDR：", lookup("discovery_robust", "TCGA-and-ICGC", "dual_genomewide_FDR_lt_0.05_and_same_direction"), "。"),
  "",
  "解释边界：候选家族校正回答锁定候选的定向复核问题；全转录组校正是更严格敏感性分析。任何一种表达重复均不等同于机制、诊断效能、预后效应或治疗获益。"
)
writeLines(review, file.path(out, "multiplicity_reconstruction_review_zh.md"), useBytes = TRUE)

message("CBC external validation multiplicity reconstruction complete")
print(audit[, .(cohort, genes_after_filter, genomewide_BH_family_n)])
print(counts[evidence_tier == "discovery_robust"])
