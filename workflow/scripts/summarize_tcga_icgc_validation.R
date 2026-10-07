#!/usr/bin/env Rscript
# Combine TCGA-LIHC and ICGC-LIRI-JP validation summaries.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(make_option("--project-root", dest = "project_root", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)

val_root <- file.path(root, "results", "validation")
fig_dir <- file.path(val_root, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

tcga <- fread(file.path(val_root, "tcga_lihc", "candidate_validation_summary.tsv"))
icgc <- fread(file.path(val_root, "icgc_liri_jp", "candidate_validation_summary.tsv"))
tcga_counts <- fread(file.path(val_root, "tcga_lihc", "validation_counts.tsv"))
icgc_counts <- fread(file.path(val_root, "icgc_liri_jp", "validation_counts.tsv"))

summary <- merge(
  tcga[, .(gene_id = as.character(gene_id), gene_symbol, evidence_tier,
           discovery_log2FC, discovery_FDR,
           tcga_log2FC = validation_log2FC, tcga_CI_low = validation_CI_low,
           tcga_CI_high = validation_CI_high, tcga_P = validation_P, tcga_FDR = validation_FDR,
           tcga_n = n, tcga_design = design, tcga_same_direction = same_direction,
           tcga_status = validation_status)],
  icgc[, .(gene_id = as.character(gene_id),
           icgc_log2FC = validation_log2FC, icgc_CI_low = validation_CI_low,
           icgc_CI_high = validation_CI_high, icgc_P = validation_P, icgc_FDR = validation_FDR,
           icgc_n = n, icgc_design = design, icgc_same_direction = same_direction,
           icgc_status = validation_status)],
  by = "gene_id", all = TRUE
)
summary[, both_same_direction := tcga_same_direction == TRUE & icgc_same_direction == TRUE]
summary[, both_validated := tcga_status %in% c("validated_primary", "validated_secondary") &
          icgc_status %in% c("validated_primary", "validated_secondary")]
fwrite(summary, file.path(val_root, "tcga_icgc_validation_summary.tsv"), sep = "\t", na = "NA")

concord <- summary[evidence_tier == "high_confidence", .(
  n = .N,
  tcga_same = sum(tcga_same_direction == TRUE, na.rm = TRUE),
  icgc_same = sum(icgc_same_direction == TRUE, na.rm = TRUE),
  both_same = sum(both_same_direction == TRUE, na.rm = TRUE),
  tcga_validated = sum(tcga_status == "validated_primary", na.rm = TRUE),
  icgc_validated = sum(icgc_status == "validated_primary", na.rm = TRUE),
  both_validated = sum(both_validated == TRUE, na.rm = TRUE)
)]
fwrite(concord, file.path(val_root, "tcga_icgc_direction_concordance.tsv"), sep = "\t")

prov <- rbind(
  fread(file.path(val_root, "tcga_lihc", "data_manifest.tsv"))[, source := "TCGA-LIHC"],
  fread(file.path(val_root, "icgc_liri_jp", "data_manifest.tsv"))[, source := "ICGC-LIRI-JP"],
  fill = TRUE
)
fwrite(prov, file.path(val_root, "tcga_icgc_data_provenance.tsv"), sep = "\t", na = "NA")

# Blockers: none for current successful downloads; keep explicit empty table with header note
blockers <- data.table(
  database = character(),
  project = character(),
  blocker = character(),
  impact = character(),
  workaround_or_status = character()
)
# Document version nuance, not a hard blocker
blockers <- rbind(blockers, data.table(
  database = "GDC/UCSC Xena",
  project = "TCGA-LIHC",
  blocker = "Xena GDC Hub matrix dated 2024-10-11 documents Release 41.0; live GDC API status at analysis time was Data Release 46.0",
  impact = "Matrix sample count (424 STAR) matches GDC API STAR-Counts file count; release labels differ",
  workaround_or_status = "Recorded in provenance; no silent cohort substitution"
))
fwrite(blockers, file.path(val_root, "data_blockers_notes.tsv"), sep = "\t")

writeLines(c(
  "# External validation data blockers",
  "",
  "No hard data blockers prevented TCGA-LIHC or ICGC-LIRI-JP analysis.",
  "",
  "## Notes (not silent substitutions)",
  "",
  "1. **TCGA Xena vs live GDC release label**",
  "   - Used UCSC Xena GDC Hub `TCGA-LIHC.star_counts.tsv.gz` and `TCGA-LIHC.clinical.tsv.gz` (Last-Modified 2024-10-11; hub documents GDC Data Release 41.0).",
  "   - Live GDC API `/status` at download/analysis time reported Data Release 46.0 (2026-08-10).",
  "   - GDC API STAR-Counts file total for TCGA-LIHC = 424, matching the Xena matrix column count.",
  "   - No alternate cancer type or unofficial matrix was substituted.",
  "",
  "2. **ICGC portal retirement**",
  "   - Interactive ICGC Data Portal was retired; open Release 28 data were obtained from the official successor open bucket `icgc25k-open`.",
  "   - Project path confirmed as `release_28/data/LIRI-JP/`.",
  "   - Paired analysis completed; no unpaired fallback.",
  "",
  "3. **Not used**",
  "   - No GEO discovery cohorts reused as external validation.",
  "   - No non-LIRI ICGC liver projects.",
  "   - No Canvas artifacts."
), file.path(val_root, "data_blockers.md"), useBytes = TRUE)

# Combined figures
high <- summary[evidence_tier == "high_confidence"]
p1 <- ggplot(high[!is.na(tcga_log2FC)], aes(discovery_log2FC, tcga_log2FC, color = tcga_status)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#777777") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#777777") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.35) +
  geom_point(alpha = 0.85, size = 1.5) +
  labs(title = "Discovery vs TCGA-LIHC paired effects", x = "Discovery log2FC", y = "TCGA log2FC", color = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "discovery_vs_tcga_effect_concordance.png"), p1, width = 7.2, height = 6.2, dpi = 320)

p2 <- ggplot(high[!is.na(icgc_log2FC)], aes(discovery_log2FC, icgc_log2FC, color = icgc_status)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "#777777") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#777777") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.35) +
  geom_point(alpha = 0.85, size = 1.5) +
  labs(title = "Discovery vs ICGC-LIRI-JP paired effects", x = "Discovery log2FC", y = "ICGC log2FC", color = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "discovery_vs_icgc_effect_concordance.png"), p2, width = 7.2, height = 6.2, dpi = 320)

dir_long <- rbind(
  high[, .(cohort = "TCGA-LIHC", same_direction = tcga_same_direction)],
  high[, .(cohort = "ICGC-LIRI-JP", same_direction = icgc_same_direction)]
)[!is.na(same_direction), .N, by = .(cohort, same_direction)]
p3 <- ggplot(dir_long, aes(cohort, N, fill = same_direction)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = N), position = position_dodge(width = 0.7), vjust = -0.3, size = 3.3) +
  labs(title = "High-confidence direction concordance", x = NULL, y = "Genes", fill = "Same direction") +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "high_confidence_direction_concordance.png"), p3, width = 6.8, height = 5.0, dpi = 320)

both_val <- high[both_validated == TRUE][order(-abs(discovery_log2FC))][1:min(15L, .N)]
if (nrow(both_val)) {
  forest <- rbind(
    both_val[, .(gene_symbol, cohort = "TCGA", log2FC = tcga_log2FC, lo = tcga_CI_low, hi = tcga_CI_high)],
    both_val[, .(gene_symbol, cohort = "ICGC", log2FC = icgc_log2FC, lo = icgc_CI_low, hi = icgc_CI_high)]
  )
  forest[, gene_symbol := factor(gene_symbol, levels = rev(unique(both_val$gene_symbol)))]
  p4 <- ggplot(forest, aes(log2FC, gene_symbol, color = cohort)) +
    geom_vline(xintercept = 0, color = "#777777", linewidth = 0.35) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.25,
                  position = position_dodge(width = 0.55)) +
    geom_point(position = position_dodge(width = 0.55), size = 1.8) +
    labs(title = "TCGA/ICGC paired validation forest (both validated, top |effect|)",
         x = "Validation log2FC (95% CI)", y = NULL, color = NULL) +
    theme_classic(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "tcga_icgc_candidate_validation_forest.png"), p4, width = 7.4, height = 7.0, dpi = 320)
}

miss <- rbind(
  high[, .(cohort = "TCGA-LIHC", status = fifelse(is.na(tcga_log2FC), "not_estimable", "estimable"))],
  high[, .(cohort = "ICGC-LIRI-JP", status = fifelse(is.na(icgc_log2FC), "not_estimable", "estimable"))]
)[, .N, by = .(cohort, status)]
p5 <- ggplot(miss, aes(cohort, N, fill = status)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = N), position = position_dodge(width = 0.7), vjust = -0.3, size = 3.3) +
  labs(title = "Missing / not-estimable audit (high_confidence)", x = NULL, y = "Genes", fill = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "missing_not_estimable_audit.png"), p5, width = 6.8, height = 4.8, dpi = 320)

# Copy cohort-specific figures into shared folder when present
file.copy(file.path(val_root, "tcga_lihc", "validation_figures", "discovery_vs_tcga_effect_concordance.png"),
          file.path(fig_dir, "tcga_discovery_vs_effect_concordance.png"), overwrite = TRUE)
file.copy(file.path(val_root, "icgc_liri_jp", "validation_figures", "discovery_vs_icgc_effect_concordance.png"),
          file.path(fig_dir, "icgc_discovery_vs_effect_concordance.png"), overwrite = TRUE)

tcga_n_pairs <- as.integer(tcga_counts[metric == "n_paired_patients", value])
icgc_n_pairs <- as.integer(icgc_counts[metric == "paired_donors", value])

review <- c(
  "Takeaway: TCGA-LIHC（50对）与 ICGC-LIRI-JP（199对）均完成配对独立验证；候选集 v2 未被回改。",
  "",
  "# TCGA / ICGC 外部验证综合审阅",
  "",
  "## 发现阶段口径（未改）",
  "",
  "- intent-to-analyse：439；QC后：438；排除：GSE57957/GSM1398656；主要配对：167对。",
  "- 主模型：tumour − control；Meta：REML + Hartung–Knapp。",
  "- 外部验证仅使用锁定的 173 high_confidence + 119 intermediate；不回筛 17,900 基因。",
  "",
  "## TCGA-LIHC",
  "",
  paste0("- 配对样本数：", tcga_n_pairs, "。"),
  paste0("- high_confidence 方向一致：", concord$tcga_same, " / 173。"),
  paste0("- high_confidence 候选集内 FDR<0.05：", tcga_counts[metric == "high_FDR_lt_0_05", value], "。"),
  paste0("- high_confidence validated_primary：", concord$tcga_validated, "。"),
  "- 表达单位：Xena GDC Hub log2(STAR count + 1)。",
  "",
  "## ICGC-LIRI-JP",
  "",
  "- 确认使用 LIRI-JP（Release 28 公开桶）。",
  paste0("- 是否配对：是；配对数 = ", icgc_n_pairs, "。"),
  "- 是否 fallback：否。",
  paste0("- high_confidence 方向一致：", concord$icgc_same, " / 173。"),
  paste0("- high_confidence 候选集内 FDR<0.05：", icgc_counts[metric == "high_FDR_lt_0_05", value], "。"),
  paste0("- high_confidence validated_primary：", concord$icgc_validated, "。"),
  paste0("- 两库方向均一致：", concord$both_same, "；两库均 validated：", concord$both_validated, "。"),
  "",
  "## 允许与禁止的表述",
  "",
  "- 允许：high_confidence 候选在独立 TCGA-LIHC / ICGC-LIRI-JP 中得到表达方向支持；结果支持跨队列表达关联可重复性。",
  "- 禁止：证明治疗获益、证明直接靶点或因果机制、声称所有候选均在两库验证成功。",
  "",
  "## 下一步",
  "",
  "本阶段在 TCGA/ICGC 后停止，等待审阅；尚未进入单细胞、HPA/CPTAC、ChEMBL/BindingDB、LINCS 或 DepMap。"
)
writeLines(review, file.path(val_root, "tcga_icgc_validation_review_zh.md"), useBytes = TRUE)

message("Combined TCGA/ICGC summary written")
