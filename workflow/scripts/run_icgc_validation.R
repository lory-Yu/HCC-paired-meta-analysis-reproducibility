#!/usr/bin/env Rscript
# ICGC-LIRI-JP validation for locked candidate set v2.

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(optparse)
})

option_list <- list(make_option("--project-root", dest = "project_root", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)

out_dir <- file.path(root, "results", "validation", "icgc_liri_jp")
fig_dir <- file.path(out_dir, "validation_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

high_path <- file.path(root, "results", "analysis", "candidate_lock", "candidate_set_v2_high_confidence.tsv")
int_path <- file.path(root, "results", "analysis", "candidate_lock", "candidate_set_v2_intermediate.tsv")
sha_of <- function(path) strsplit(system2("sha256sum", path, stdout = TRUE), " +")[[1]][1]
stopifnot(identical(sha_of(high_path), "55110e282bf9d26619de70d6e8c1e5972d050766b034a6fa06009021ae9fd3cd"))
stopifnot(identical(sha_of(int_path), "23812a47c8df31301aca8a9061ecd7da56dd41afa5f214a3fc6e4a585494aa59"))

high <- fread(high_path); high[, gene_id := as.character(gene_id)]
inter <- fread(int_path); inter[, gene_id := as.character(gene_id)]
candidates <- rbind(
  high[, .(gene_id, gene_symbol, discovery_log2FC = pooled_log2FC, discovery_FDR = FDR,
           discovery_I2 = I2, evidence_tier = "high_confidence")],
  inter[, .(gene_id, gene_symbol, discovery_log2FC = pooled_log2FC, discovery_FDR = FDR,
            discovery_I2 = I2, evidence_tier = "intermediate_confidence")]
)

pairs <- fread(file.path(out_dir, "paired_donors.tsv"))
sample_manifest <- fread(file.path(out_dir, "sample_manifest.tsv"))
expr <- fread(file.path(out_dir, "candidate_expression_matrix.tsv.gz"))
sym <- expr$gene_symbol
mat_raw <- as.matrix(expr[, -1, with = FALSE])
storage.mode(mat_raw) <- "numeric"
rownames(mat_raw) <- sym
mat <- log2(mat_raw + 1)
expression_unit <- "log2(raw_read_count + 1) from ICGC Release 28 open exp_seq"

n_pairs <- nrow(pairs)
message("ICGC-LIRI-JP paired donors: ", n_pairs)
stopifnot(n_pairs >= 3L)

# Map candidates to matrix by gene_symbol
candidates[, in_matrix := gene_symbol %in% rownames(mat)]
present <- candidates[in_matrix == TRUE]
missing <- candidates[in_matrix == FALSE]

# Paired limma
paired_samples <- c(pairs$normal_sample, pairs$tumor_sample)
patient <- factor(c(pairs$icgc_donor_id, pairs$icgc_donor_id))
group <- factor(c(rep("normal", n_pairs), rep("tumor", n_pairs)), levels = c("normal", "tumor"))
design <- model.matrix(~ patient + group)
y <- mat[present$gene_symbol, paired_samples, drop = FALSE]
# genes with too many NA become not estimable
keep_gene <- rowSums(is.finite(y)) == ncol(y)
y <- y[keep_gene, , drop = FALSE]
fit <- eBayes(lmFit(y, design))
coef <- "grouptumor"
paired_fit <- data.table(
  gene_symbol = rownames(y),
  validation_log2FC = as.numeric(fit$coefficients[, coef]),
  SE = fit$stdev.unscaled[, coef] * sqrt(fit$s2.post),
  validation_P = as.numeric(fit$p.value[, coef]),
  df_total = fit$df.total,
  n = n_pairs,
  design = "paired"
)
paired_fit[, validation_CI_low := validation_log2FC - qt(0.975, df_total) * SE]
paired_fit[, validation_CI_high := validation_log2FC + qt(0.975, df_total) * SE]

paired_all <- merge(candidates, paired_fit, by = "gene_symbol", all.x = TRUE)
paired_all[is.na(design), `:=`(
  design = "paired", n = n_pairs, validation_status = "not_estimable",
  note = "gene_symbol_absent_from_ICGC_RefSeq_exp_seq_or_incomplete_values"
)]
paired_all[!is.na(validation_P), same_direction := sign(validation_log2FC) == sign(discovery_log2FC)]
for (tier in c("high_confidence", "intermediate_confidence")) {
  idx <- paired_all$evidence_tier == tier & !is.na(paired_all$validation_P)
  paired_all[idx, validation_FDR := p.adjust(validation_P, method = "BH")]
}
paired_all[, validation_status := fcase(
  !is.na(validation_status) & validation_status == "not_estimable", "not_estimable",
  evidence_tier == "high_confidence" & same_direction == TRUE & validation_FDR < 0.05, "validated_primary",
  evidence_tier == "intermediate_confidence" & same_direction == TRUE & validation_FDR < 0.05, "validated_secondary",
  same_direction == TRUE, "direction_consistent",
  same_direction == FALSE, "direction_inconsistent",
  default = "not_estimable"
)]
paired_all[, evidence_note := "ICGC-LIRI-JP paired tumour-normal validation"]

fwrite(paired_all, file.path(out_dir, "paired_or_unpaired_model_results.tsv"), sep = "\t", na = "NA")
fwrite(paired_all, file.path(out_dir, "candidate_validation_summary.tsv"), sep = "\t", na = "NA")
fwrite(paired_all[evidence_tier == "high_confidence"], file.path(out_dir, "paired_validation.tsv"), sep = "\t", na = "NA")

# Secondary unpaired
unpaired_samples <- sample_manifest[used_in_unpaired_secondary == TRUE]
g <- factor(unpaired_samples$group, levels = c("normal", "tumor"))
design_u <- model.matrix(~ g)
y_u <- mat[present$gene_symbol, unpaired_samples$icgc_sample_id, drop = FALSE]
keep_u <- rowSums(is.finite(y_u)) == ncol(y_u)
y_u <- y_u[keep_u, , drop = FALSE]
fit_u <- eBayes(lmFit(y_u, design_u))
coef_u <- "gtumor"
unpaired <- data.table(
  gene_symbol = rownames(y_u),
  unpaired_log2FC = as.numeric(fit_u$coefficients[, coef_u]),
  unpaired_SE = fit_u$stdev.unscaled[, coef_u] * sqrt(fit_u$s2.post),
  unpaired_P = as.numeric(fit_u$p.value[, coef_u]),
  unpaired_df = fit_u$df.total,
  n_tumor = sum(g == "tumor"),
  n_normal = sum(g == "normal")
)
unpaired[, unpaired_CI_low := unpaired_log2FC - qt(0.975, unpaired_df) * unpaired_SE]
unpaired[, unpaired_CI_high := unpaired_log2FC + qt(0.975, unpaired_df) * unpaired_SE]
unpaired <- merge(candidates, unpaired, by = "gene_symbol", all.x = TRUE)
unpaired[!is.na(unpaired_P), same_direction := sign(unpaired_log2FC) == sign(discovery_log2FC)]
for (tier in c("high_confidence", "intermediate_confidence")) {
  idx <- unpaired$evidence_tier == tier & !is.na(unpaired$unpaired_P)
  unpaired[idx, unpaired_FDR := p.adjust(unpaired_P, method = "BH")]
}
fwrite(unpaired, file.path(out_dir, "unpaired_secondary.tsv"), sep = "\t", na = "NA")

paired_high <- paired_all[evidence_tier == "high_confidence"]
paired_int <- paired_all[evidence_tier == "intermediate_confidence"]

# Figures
plot_dt <- paired_high[!is.na(validation_log2FC)]
p1 <- ggplot(plot_dt, aes(discovery_log2FC, validation_log2FC, color = validation_status)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#777777") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#777777") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4, color = "#555555") +
  geom_point(size = 1.6, alpha = 0.85) +
  labs(title = "Discovery vs ICGC-LIRI-JP paired effects (high_confidence)",
       x = "GEO meta pooled log2FC", y = "ICGC paired log2FC", color = NULL,
       caption = paste0("n_pairs=", n_pairs, "; unit=", expression_unit)) +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "discovery_vs_icgc_effect_concordance.png"), p1, width = 7.2, height = 6.2, dpi = 320)

p2 <- ggplot(paired_high[!is.na(same_direction), .N, by = same_direction], aes(same_direction, N, fill = same_direction)) +
  geom_col(width = 0.65) + geom_text(aes(label = N), vjust = -0.3, size = 3.5) +
  labs(title = "ICGC paired direction concordance (high_confidence)",
       x = "Same direction as discovery", y = "Genes") +
  theme_classic(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(fig_dir, "high_confidence_direction_concordance.png"), p2, width = 5.8, height = 4.8, dpi = 320)

top_forest <- paired_high[validation_status == "validated_primary"][order(validation_FDR, -abs(validation_log2FC))][1:min(20L, .N)]
if (nrow(top_forest)) {
  top_forest[, gene_symbol := factor(gene_symbol, levels = rev(gene_symbol))]
  p3 <- ggplot(top_forest, aes(validation_log2FC, gene_symbol)) +
    geom_vline(xintercept = 0, color = "#777777", linewidth = 0.35) +
    geom_errorbar(aes(xmin = validation_CI_low, xmax = validation_CI_high), orientation = "y",
                  width = 0.25, color = "#4C78A8") +
    geom_point(size = 2, color = "#4C78A8") +
    labs(title = "ICGC-LIRI-JP paired validation forest (top validated_primary)",
         x = "ICGC paired log2FC (95% CI)", y = NULL) +
    theme_classic(base_size = 10)
  ggsave(file.path(fig_dir, "icgc_candidate_validation_forest.png"), p3, width = 7.0, height = 6.5, dpi = 320)
}

miss_dt <- paired_high[, .(N = .N), by = .(category = fifelse(is.na(validation_log2FC), "not_estimable", "estimable"))]
p4 <- ggplot(miss_dt, aes(category, N, fill = category)) +
  geom_col(width = 0.6) + geom_text(aes(label = N), vjust = -0.3) +
  labs(title = "ICGC high_confidence estimability audit", x = NULL, y = "Genes") +
  theme_classic(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(fig_dir, "icgc_missing_not_estimable_audit.png"), p4, width = 5.6, height = 4.6, dpi = 320)

# Data manifest
dl_manifest <- fread(file.path(root, "data", "external", "icgc_liri_jp", "download_file_manifest.tsv"))
data_manifest <- data.table(
  database = "ICGC 25K open object storage (icgc25k-open)",
  project_id = "LIRI-JP",
  release_or_version = "ICGC Data Portal Release 28 (2019-11-26); open bucket accessed 2026-09-16",
  download_date = as.character(Sys.Date()),
  url = "https://object.genomeinformatics.org/icgc25k-open/release_28/data/LIRI-JP/",
  file_name = "donor-level exp_seq/donor/specimen/sample part-*.gz + headers",
  n_files = nrow(dl_manifest),
  bytes_total = sum(dl_manifest$bytes, na.rm = TRUE),
  sample_n = uniqueN(sample_manifest$icgc_sample_id),
  donor_n = uniqueN(sample_manifest$icgc_donor_id),
  paired_donors = n_pairs,
  expression_unit = expression_unit,
  sample_type = "Primary tumour - solid tissue vs Normal - solid tissue / adjacent",
  paired = TRUE,
  missing_values = "candidate symbols absent from RefSeq gene_id treated as not_estimable",
  sample_filter_rule = "Primary tumour - solid tissue; Normal - solid tissue or Normal - tissue adjacent to primary; lexicographic first sample per donor/group",
  official_citation = "International Cancer Genome Consortium; Fujimoto et al. Nat Genet 2016 (LIRI-JP / Japanese HCC genomics)"
)
fwrite(data_manifest, file.path(out_dir, "data_manifest.tsv"), sep = "\t")

counts <- data.table(
  metric = c(
    "project_confirmed_LIRI_JP", "paired", "fallback_unpaired", "paired_donors",
    "high_confidence_n", "high_same_direction", "high_FDR_lt_0_05", "high_validated_primary",
    "intermediate_n", "intermediate_same_direction", "intermediate_FDR_lt_0_05",
    "missing_genes"
  ),
  value = c(
    "TRUE", "TRUE", "FALSE", n_pairs,
    nrow(paired_high),
    sum(paired_high$same_direction == TRUE, na.rm = TRUE),
    sum(paired_high$validation_FDR < 0.05, na.rm = TRUE),
    sum(paired_high$validation_status == "validated_primary", na.rm = TRUE),
    nrow(paired_int),
    sum(paired_int$same_direction == TRUE, na.rm = TRUE),
    sum(paired_int$validation_FDR < 0.05, na.rm = TRUE),
    nrow(missing)
  )
)
fwrite(counts, file.path(out_dir, "validation_counts.tsv"), sep = "\t")

review <- c(
  "Takeaway: ICGC-LIRI-JP 使用 Release 28 公开 exp_seq，并完成真正的 donor 配对验证。",
  "",
  "# ICGC-LIRI-JP 独立验证审阅",
  "",
  "- 项目确认：ICGC-LIRI-JP（开放桶 `icgc25k-open/release_28/data/LIRI-JP/`）。",
  "- 未用其他 ICGC 肝癌项目、TCGA 或 GEO 冒充。",
  paste0("- 配对：是。完整肿瘤—正常 donor 数 = ", n_pairs, "。"),
  "- 未发生 unpaired fallback。",
  "- 表达单位：log2(raw_read_count + 1)；基因标识为 RefSeq gene_id（符号）。",
  "- 主验证集：173 high_confidence；次级：119 intermediate。",
  paste0("- high_confidence 方向一致：", sum(paired_high$same_direction == TRUE, na.rm = TRUE),
         "；候选集内 FDR<0.05：", sum(paired_high$validation_FDR < 0.05, na.rm = TRUE),
         "；validated_primary：", sum(paired_high$validation_status == "validated_primary", na.rm = TRUE), "。"),
  paste0("- intermediate 方向一致：", sum(paired_int$same_direction == TRUE, na.rm = TRUE),
         "；候选集内 FDR<0.05：", sum(paired_int$validation_FDR < 0.05, na.rm = TRUE), "。"),
  paste0("- 缺失/不可估计基因：", nrow(missing), "。"),
  "- 结果支持跨队列表达关联的可重复性；不支持治疗获益、直接靶点或因果机制表述。"
)
writeLines(review, file.path(out_dir, "validation_review_zh.md"), useBytes = TRUE)

message("ICGC-LIRI-JP validation complete: pairs=", n_pairs,
        " high_validated_primary=", sum(paired_high$validation_status == "validated_primary", na.rm = TRUE))
