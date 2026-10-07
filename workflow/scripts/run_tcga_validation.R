#!/usr/bin/env Rscript
# TCGA-LIHC paired validation for locked candidate set v2.
# Does not modify discovery meta tables or candidate locks.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(optparse)
})

option_list <- list(make_option("--project-root", dest = "project_root", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)

data_dir <- file.path(root, "data", "external", "tcga_lihc")
out_dir <- file.path(root, "results", "validation", "tcga_lihc")
fig_dir <- file.path(out_dir, "validation_figures")
for (p in c(out_dir, fig_dir)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

high_path <- file.path(root, "results", "analysis", "candidate_lock", "candidate_set_v2_high_confidence.tsv")
int_path <- file.path(root, "results", "analysis", "candidate_lock", "candidate_set_v2_intermediate.tsv")
expected_high <- "55110e282bf9d26619de70d6e8c1e5972d050766b034a6fa06009021ae9fd3cd"
expected_int <- "23812a47c8df31301aca8a9061ecd7da56dd41afa5f214a3fc6e4a585494aa59"
sha_of <- function(path) strsplit(system2("sha256sum", path, stdout = TRUE), " +")[[1]][1]
stopifnot(identical(sha_of(high_path), expected_high))
stopifnot(identical(sha_of(int_path), expected_int))

high <- fread(high_path); high[, gene_id := as.character(gene_id)]
inter <- fread(int_path); inter[, gene_id := as.character(gene_id)]
candidates <- rbind(
  high[, .(gene_id, gene_symbol, discovery_log2FC = pooled_log2FC, discovery_FDR = FDR,
           discovery_I2 = I2, evidence_tier = "high_confidence")],
  inter[, .(gene_id, gene_symbol, discovery_log2FC = pooled_log2FC, discovery_FDR = FDR,
            discovery_I2 = I2, evidence_tier = "intermediate_confidence")]
)
stopifnot(nrow(high) == 173L, nrow(inter) == 119L, nrow(candidates) == 292L)

expr_path <- file.path(data_dir, "TCGA-LIHC.star_counts.tsv.gz")
clin_path <- file.path(data_dir, "TCGA-LIHC.clinical.tsv.gz")
stopifnot(file.exists(expr_path), file.exists(clin_path))

message("Reading Xena GDC TCGA-LIHC STAR counts matrix")
expr <- fread(expr_path)
ensembl_raw <- expr[[1]]
ensembl <- sub("\\..*$", "", ensembl_raw)
mat <- as.matrix(expr[, -1, with = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- ensembl_raw
sample_ids <- colnames(mat)

# Xena GDC ETL stores star_counts as log2(count + 1).
expression_unit <- "xena_gdc_log2_star_count_plus_1"

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

gene_map <- data.table(
  ensembl_versioned = ensembl_raw,
  ensembl_id = ensembl,
  gene_id = map_unique(ensembl, "ENSEMBL", "ENTREZID"),
  mapped_symbol = map_unique(ensembl, "ENSEMBL", "SYMBOL")
)
# Keep unique Entrez only; ambiguous Ensembl→Entrez already NA.
# Aggregate duplicate Entrez by arithmetic mean across Ensembl rows (prespecified).
gene_map_ok <- gene_map[!is.na(gene_id)]
mat_entrez <- mat[gene_map_ok$ensembl_versioned, , drop = FALSE]
agg <- lapply(split(seq_len(nrow(gene_map_ok)), gene_map_ok$gene_id), function(idx) {
  if (length(idx) == 1L) mat_entrez[idx, , drop = FALSE] else {
    m <- colMeans(mat_entrez[idx, , drop = FALSE], na.rm = TRUE)
    matrix(m, nrow = 1, dimnames = list(gene_map_ok$gene_id[idx][1], colnames(mat_entrez)))
  }
})
expr_entrez <- do.call(rbind, agg)
rownames(expr_entrez) <- names(agg)

clin <- fread(clin_path, na.strings = c("", "NA", "--", "not reported", "Not Reported"))
# clinical rows may exceed expression columns; keep expression barcodes.
clin[, sample_barcode := as.character(sample)]
clin <- clin[sample_barcode %in% sample_ids]
clin[, patient_id := substr(sample_barcode, 1, 12)]
clin[, sample_type := as.character(`sample_type.samples`)]
clin[, sample_type_id := as.character(`sample_type_id.samples`)]
# Prefer official sample_type; also parse barcode code as audit.
clin[, barcode_type_code := substr(sample_barcode, 14, 15)]
clin[, group := fcase(
  sample_type == "Primary Tumor" | barcode_type_code == "01", "tumor",
  sample_type == "Solid Tissue Normal" | barcode_type_code == "11", "normal",
  default = "other"
)]

# Deterministic duplicate resolution: among same patient+group, keep lexicographically
# first full sample barcode. Never choose by expression value.
setorder(clin, patient_id, group, sample_barcode)
clin[, keep_primary := FALSE]
clin[group %in% c("tumor", "normal"), keep_primary := seq_len(.N) == 1L, by = .(patient_id, group)]
excluded <- clin[group %in% c("tumor", "normal") & !keep_primary,
                 .(sample_barcode, patient_id, group, sample_type, reason = "duplicate_patient_group_lexicographic_first_kept")]
excluded <- rbind(
  excluded,
  clin[group == "other",
       .(sample_barcode, patient_id, group, sample_type, reason = "not_primary_tumor_or_solid_tissue_normal")],
  fill = TRUE
)

kept <- clin[keep_primary == TRUE]
tumor <- kept[group == "tumor", .(patient_id, tumor_sample = sample_barcode, tumor_type = sample_type)]
normal <- kept[group == "normal", .(patient_id, normal_sample = sample_barcode, normal_type = sample_type)]
pairs <- merge(tumor, normal, by = "patient_id")
setorder(pairs, patient_id)
n_pairs <- nrow(pairs)
message("TCGA-LIHC complete tumor-normal pairs: ", n_pairs)

sample_manifest <- clin[, .(
  sample_barcode, patient_id, sample_type, sample_type_id, barcode_type_code, group,
  in_expression_matrix = sample_barcode %in% sample_ids,
  used_in_paired_primary = sample_barcode %in% c(pairs$tumor_sample, pairs$normal_sample),
  used_in_unpaired_secondary = group %in% c("tumor", "normal") & keep_primary
)]
fwrite(sample_manifest, file.path(out_dir, "sample_manifest.tsv"), sep = "\t", na = "NA")
fwrite(excluded, file.path(out_dir, "excluded_samples.tsv"), sep = "\t", na = "NA")
fwrite(pairs, file.path(out_dir, "paired_patients.tsv"), sep = "\t", na = "NA")

cand_ids <- candidates$gene_id
present <- cand_ids[cand_ids %in% rownames(expr_entrez)]
missing_genes <- cand_ids[!cand_ids %in% rownames(expr_entrez)]
cand_mat <- expr_entrez[present, , drop = FALSE]
fwrite(
  data.table(gene_id = rownames(cand_mat), cand_mat, keep.rownames = FALSE),
  file.path(out_dir, "candidate_expression_matrix.tsv.gz"),
  sep = "\t", compress = "gzip"
)

message("Running paired limma jointly on locked candidates")
paired_samples <- c(pairs$normal_sample, pairs$tumor_sample)
patient <- factor(c(pairs$patient_id, pairs$patient_id))
group <- factor(c(rep("normal", n_pairs), rep("tumor", n_pairs)), levels = c("normal", "tumor"))
design_paired <- model.matrix(~ patient + group)
present_ids <- candidates$gene_id[candidates$gene_id %in% rownames(cand_mat)]
y_paired <- cand_mat[present_ids, paired_samples, drop = FALSE]
fit_paired_all <- eBayes(lmFit(y_paired, design_paired))
coef <- "grouptumor"
paired_fit <- data.table(
  gene_id = present_ids,
  validation_log2FC = as.numeric(fit_paired_all$coefficients[, coef]),
  SE = fit_paired_all$stdev.unscaled[, coef] * sqrt(fit_paired_all$s2.post),
  validation_P = as.numeric(fit_paired_all$p.value[, coef]),
  df_total = fit_paired_all$df.total
)
paired_fit[, validation_CI_low := validation_log2FC - qt(0.975, df_total) * SE]
paired_fit[, validation_CI_high := validation_log2FC + qt(0.975, df_total) * SE]
paired_fit[, n := n_pairs]
paired_fit[, design := "paired"]

assign_status <- function(dt) {
  dt <- merge(candidates, dt, by = "gene_id", all.x = TRUE)
  dt[is.na(design), `:=`(design = "paired", n = n_pairs, validation_status = "not_estimable",
                         note = "gene_absent_from_TCGA_matrix_after_unique_Entrez_mapping")]
  dt[!is.na(validation_P), same_direction := sign(validation_log2FC) == sign(discovery_log2FC)]
  for (tier in c("high_confidence", "intermediate_confidence")) {
    idx <- dt$evidence_tier == tier & !is.na(dt$validation_P)
    dt[idx, validation_FDR := p.adjust(validation_P, method = "BH")]
  }
  dt[, validation_status := fcase(
    !is.na(validation_status) & validation_status == "not_estimable", "not_estimable",
    evidence_tier == "high_confidence" & same_direction == TRUE & validation_FDR < 0.05, "validated_primary",
    evidence_tier == "intermediate_confidence" & same_direction == TRUE & validation_FDR < 0.05, "validated_secondary",
    same_direction == TRUE, "direction_consistent",
    same_direction == FALSE, "direction_inconsistent",
    default = "not_estimable"
  )]
  dt
}

paired_all <- assign_status(paired_fit)
paired_high <- paired_all[evidence_tier == "high_confidence"]
paired_int <- paired_all[evidence_tier == "intermediate_confidence"]
fwrite(paired_all, file.path(out_dir, "paired_or_unpaired_model_results.tsv"), sep = "\t", na = "NA")
fwrite(paired_all, file.path(out_dir, "candidate_validation_summary.tsv"), sep = "\t", na = "NA")
fwrite(paired_high, file.path(out_dir, "paired_validation.tsv"), sep = "\t", na = "NA")

# Secondary unpaired: all primary tumor vs solid tissue normal (one sample/patient/group).
unpaired_samples <- kept[group %in% c("tumor", "normal")]
g_unpaired <- factor(unpaired_samples$group, levels = c("normal", "tumor"))
design_unpaired <- model.matrix(~ g_unpaired)
y_unpaired <- cand_mat[present_ids, unpaired_samples$sample_barcode, drop = FALSE]
fit_unpaired_all <- eBayes(lmFit(y_unpaired, design_unpaired))
coef_u <- "g_unpairedtumor"
unpaired <- data.table(
  gene_id = present_ids,
  unpaired_log2FC = as.numeric(fit_unpaired_all$coefficients[, coef_u]),
  unpaired_SE = fit_unpaired_all$stdev.unscaled[, coef_u] * sqrt(fit_unpaired_all$s2.post),
  unpaired_P = as.numeric(fit_unpaired_all$p.value[, coef_u]),
  unpaired_df = fit_unpaired_all$df.total,
  n_tumor = sum(g_unpaired == "tumor"),
  n_normal = sum(g_unpaired == "normal")
)
unpaired[, unpaired_CI_low := unpaired_log2FC - qt(0.975, unpaired_df) * unpaired_SE]
unpaired[, unpaired_CI_high := unpaired_log2FC + qt(0.975, unpaired_df) * unpaired_SE]
unpaired <- merge(candidates, unpaired, by = "gene_id", all.x = TRUE)
unpaired[!is.na(unpaired_P), same_direction := sign(unpaired_log2FC) == sign(discovery_log2FC)]
for (tier in c("high_confidence", "intermediate_confidence")) {
  idx <- unpaired$evidence_tier == tier & !is.na(unpaired$unpaired_P)
  unpaired[idx, unpaired_FDR := p.adjust(unpaired_P, method = "BH")]
}
fwrite(unpaired, file.path(out_dir, "unpaired_secondary.tsv"), sep = "\t", na = "NA")

# Optional clinical association on tumor samples only for high_confidence (secondary).
tumor_only <- kept[group == "tumor"]
clin_fields <- clin[, .(
  sample_barcode,
  pathologic_stage = `ajcc_pathologic_stage.diagnoses`,
  tumor_grade = `tumor_grade.diagnoses`,
  vital_status = `vital_status.demographic`
)]
tumor_only <- merge(tumor_only, clin_fields, by = "sample_barcode", all.x = TRUE)
clinical_rows <- list()
for (gid in high$gene_id) {
  if (!gid %in% rownames(cand_mat)) next
  dt <- copy(tumor_only)
  dt[, expr := as.numeric(cand_mat[gid, sample_barcode])]
  dt <- dt[is.finite(expr)]
  stage <- dt[!is.na(pathologic_stage) & nzchar(pathologic_stage)]
  stage_p <- if (uniqueN(stage$pathologic_stage) >= 2L && nrow(stage) >= 10L) {
    tryCatch(anova(lm(expr ~ factor(pathologic_stage), data = stage))$`Pr(>F)`[1], error = function(e) NA_real_)
  } else NA_real_
  clinical_rows[[length(clinical_rows) + 1L]] <- data.table(
    gene_id = gid,
    gene_symbol = high[gene_id == gid, gene_symbol][1],
    n_tumor_expression = nrow(dt),
    stage_anova_p = stage_p,
    analysis = "secondary_clinical_association_not_efficacy"
  )
}
clinical_tab <- if (length(clinical_rows)) rbindlist(clinical_rows) else data.table()
if (nrow(clinical_tab)) {
  clinical_tab[, stage_anova_FDR := p.adjust(stage_anova_p, method = "BH")]
  fwrite(clinical_tab, file.path(out_dir, "clinical_association_secondary.tsv"), sep = "\t", na = "NA")
}

# Figures
plot_dt <- paired_high[!is.na(validation_log2FC)]
p1 <- ggplot(plot_dt, aes(discovery_log2FC, validation_log2FC, color = validation_status)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#777777") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#777777") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4, color = "#555555") +
  geom_point(size = 1.6, alpha = 0.85) +
  labs(title = "Discovery vs TCGA-LIHC paired effects (high_confidence)",
       x = "GEO meta pooled log2FC", y = "TCGA paired log2FC", color = NULL,
       caption = paste0("n_pairs=", n_pairs, "; unit=", expression_unit)) +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "discovery_vs_tcga_effect_concordance.png"), p1, width = 7.2, height = 6.2, dpi = 320)

dir_counts <- paired_high[, .N, by = .(same_direction, validation_status)]
p2 <- ggplot(paired_high[!is.na(same_direction), .N, by = same_direction], aes(same_direction, N, fill = same_direction)) +
  geom_col(width = 0.65) + geom_text(aes(label = N), vjust = -0.3, size = 3.5) +
  labs(title = "TCGA paired direction concordance (high_confidence)", x = "Same direction as discovery", y = "Genes") +
  theme_classic(base_size = 11) + theme(legend.position = "none") +
  expand_limits(y = max(1, paired_high[!is.na(same_direction), .N, by = same_direction][, max(N)]) * 1.1)
ggsave(file.path(fig_dir, "high_confidence_direction_concordance.png"), p2, width = 5.8, height = 4.8, dpi = 320)

top_forest <- paired_high[validation_status == "validated_primary"][order(validation_FDR, -abs(validation_log2FC))][1:min(20L, .N)]
if (nrow(top_forest)) {
  top_forest[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  p3 <- ggplot(top_forest, aes(validation_log2FC, gene_symbol)) +
    geom_vline(xintercept = 0, color = "#777777", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = validation_CI_low, xmax = validation_CI_high), height = 0.25, color = "#4C78A8") +
    geom_point(size = 2, color = "#4C78A8") +
    labs(title = "TCGA paired validation forest (top validated_primary)",
         x = "TCGA paired log2FC (95% CI)", y = NULL) +
    theme_classic(base_size = 10)
  ggsave(file.path(fig_dir, "tcga_candidate_validation_forest.png"), p3, width = 7.0, height = 6.5, dpi = 320)
}

miss_dt <- paired_high[, .(N = .N), by = .(category = fifelse(is.na(validation_log2FC), "not_estimable", "estimable"))]
p4 <- ggplot(miss_dt, aes(category, N, fill = category)) +
  geom_col(width = 0.6) + geom_text(aes(label = N), vjust = -0.3) +
  labs(title = "TCGA high_confidence estimability audit", x = NULL, y = "Genes") +
  theme_classic(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(fig_dir, "tcga_missing_not_estimable_audit.png"), p4, width = 5.6, height = 4.6, dpi = 320)

# Data manifest
dl_date <- as.character(Sys.Date())
file_info <- function(path) {
  info <- file.info(path)
  data.table(
    file_name = basename(path),
    local_path = path,
    bytes = as.numeric(info$size),
    sha256 = sha_of(path)
  )
}
manifest <- rbind(
  file_info(expr_path)[, `:=`(
    database = "UCSC Xena GDC Hub / AWS gdc-hub",
    project_id = "TCGA-LIHC",
    release_or_version = "Xena GDC Hub matrices dated 2024-10-11 (hub documents GDC Data Release 41.0); GDC API status at download-time checked separately",
    download_date = dl_date,
    url = "https://gdc.xenahubs.net/download/TCGA-LIHC.star_counts.tsv.gz",
    expression_unit = expression_unit,
    sample_type = "RNA-seq STAR counts matrix columns are sample barcodes",
    paired = "pairing derived from patient barcode + Primary Tumor / Solid Tissue Normal",
    missing_values = "finite-value filter per gene/pair",
    sample_filter_rule = "Primary Tumor (01) and Solid Tissue Normal (11); exclude recurrent; lexicographic first aliquot per patient/group",
    official_citation = "Grossman et al. N Engl J Med 2016 (GDC); Goldman et al. Nat Biotechnol 2020 (UCSC Xena)"
  )],
  file_info(clin_path)[, `:=`(
    database = "UCSC Xena GDC Hub / AWS gdc-hub",
    project_id = "TCGA-LIHC",
    release_or_version = "Xena GDC Hub clinical matrix dated 2024-10-11",
    download_date = dl_date,
    url = "https://gdc.xenahubs.net/download/TCGA-LIHC.clinical.tsv.gz",
    expression_unit = "NA",
    sample_type = "GDC phenotype/clinical fields including sample_type.samples",
    paired = "supports pairing via sample barcode / case",
    missing_values = "clinical fields may be missing",
    sample_filter_rule = "intersect with expression barcodes",
    official_citation = "Grossman et al. N Engl J Med 2016 (GDC); Goldman et al. Nat Biotechnol 2020 (UCSC Xena)"
  )]
)
fwrite(manifest, file.path(out_dir, "data_manifest.tsv"), sep = "\t", na = "NA")

# GDC API cross-check note
gdc_status_path <- file.path(out_dir, "gdc_api_status.json")
if (file.exists(file.path(data_dir, "gdc_api_status.json"))) {
  file.copy(file.path(data_dir, "gdc_api_status.json"), gdc_status_path, overwrite = TRUE)
}

summary_counts <- data.table(
  metric = c(
    "n_expression_samples", "n_primary_tumor_in_matrix", "n_solid_tissue_normal_in_matrix",
    "n_paired_patients", "high_confidence_n", "high_same_direction", "high_FDR_lt_0_05",
    "high_validated_primary", "intermediate_n", "intermediate_same_direction",
    "intermediate_FDR_lt_0_05", "missing_genes_in_tcga"
  ),
  value = c(
    length(sample_ids),
    sum(substr(sample_ids, 14, 15) == "01"),
    sum(substr(sample_ids, 14, 15) == "11"),
    n_pairs,
    nrow(paired_high),
    sum(paired_high$same_direction == TRUE, na.rm = TRUE),
    sum(paired_high$validation_FDR < 0.05, na.rm = TRUE),
    sum(paired_high$validation_status == "validated_primary", na.rm = TRUE),
    nrow(paired_int),
    sum(paired_int$same_direction == TRUE, na.rm = TRUE),
    sum(paired_int$validation_FDR < 0.05, na.rm = TRUE),
    length(missing_genes)
  )
)
fwrite(summary_counts, file.path(out_dir, "validation_counts.tsv"), sep = "\t")

review <- c(
  "Takeaway: TCGA-LIHC 配对验证只检验已锁定的候选集 v2，不回改发现阈值。",
  "",
  "# TCGA-LIHC 独立验证审阅",
  "",
  paste0("- 表达矩阵：UCSC Xena GDC Hub `TCGA-LIHC.star_counts.tsv.gz`（", expression_unit, "）。"),
  "- 表型：`TCGA-LIHC.clinical.tsv.gz`，使用官方 `sample_type.samples`。",
  paste0("- 完整肿瘤—正常配对数：", n_pairs, "。"),
  "- 主验证集：173 个 high_confidence；次级：119 个 intermediate。",
  "- FDR：high_confidence 与 intermediate 分别在各自候选集内做 BH，不把 17,900 基因组重新筛选。",
  paste0("- high_confidence 方向一致：", sum(paired_high$same_direction == TRUE, na.rm = TRUE),
         "；候选集内 FDR<0.05：", sum(paired_high$validation_FDR < 0.05, na.rm = TRUE),
         "；validated_primary：", sum(paired_high$validation_status == "validated_primary", na.rm = TRUE), "。"),
  paste0("- intermediate 方向一致：", sum(paired_int$same_direction == TRUE, na.rm = TRUE),
         "；候选集内 FDR<0.05：", sum(paired_int$validation_FDR < 0.05, na.rm = TRUE), "。"),
  paste0("- 映射后缺失基因数：", length(missing_genes), "。"),
  "- 未配对分析和临床关联仅为次级，不能替代配对主验证，也不能写成疗效。",
  "- 本阶段支持组织层面的表达复现，不支持治疗获益、直接靶点结合或因果机制表述。"
)
writeLines(review, file.path(out_dir, "validation_review_zh.md"), useBytes = TRUE)

message("TCGA-LIHC validation complete: pairs=", n_pairs,
        " high_validated_primary=", sum(paired_high$validation_status == "validated_primary", na.rm = TRUE))
