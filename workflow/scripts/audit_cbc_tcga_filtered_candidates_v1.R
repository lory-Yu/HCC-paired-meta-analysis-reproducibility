#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--project-root", dest = "project_root", type = "character")
)))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)
out <- file.path(root, "results", "validation", "cbc_external_validation_v1")

count_path <- file.path(root, "data", "external", "tcga_lihc", "TCGA-LIHC.star_counts.tsv.gz")
pair_path <- file.path(root, "results", "validation", "tcga_lihc", "paired_patients.tsv")
candidate_path <- file.path(out, "tcga_locked_candidate_multiplicity.tsv")

map_unique <- function(keys, keytype, column) {
  mapped <- mapIds(
    org.Hs.eg.db, keys = unique(as.character(keys)), keytype = keytype,
    column = column,
    multiVals = function(z) {
      z <- unique(z[!is.na(z) & nzchar(z)])
      if (length(z) == 1L) z else NA_character_
    }
  )
  unname(mapped[as.character(keys)])
}

pairs <- fread(pair_path)
samples <- c(pairs$normal_sample, pairs$tumor_sample)
group <- factor(c(rep("normal", nrow(pairs)), rep("tumor", nrow(pairs))),
                levels = c("normal", "tumor"))

raw <- fread(count_path, select = c("Ensembl_ID", samples))
ens_versioned <- raw$Ensembl_ID
ens <- sub("\\..*$", "", ens_versioned)
mat_log <- as.matrix(raw[, -1, with = FALSE])
storage.mode(mat_log) <- "numeric"
mat <- round(pmax(2^mat_log - 1, 0))
entrez <- map_unique(ens, "ENSEMBL", "ENTREZID")
ok <- !is.na(entrez)
counts <- rowsum(mat[ok, , drop = FALSE], group = entrez[ok], reorder = FALSE)
storage.mode(counts) <- "numeric"

dge <- DGEList(counts = counts)
keep <- filterByExpr(dge, group = group)
cpm_mat <- cpm(dge, normalized.lib.sizes = FALSE)

cand <- fread(candidate_path, colClasses = list(character = "gene_id"))
target <- cand[evidence_tier == "discovery_robust" & present_before_filter == TRUE & retained_after_filter == FALSE,
               .(gene_id, gene_symbol, discovery_log2FC)]
stopifnot(nrow(target) == 5L, all(target$gene_id %in% rownames(counts)))

rows <- lapply(seq_len(nrow(target)), function(i) {
  g <- target$gene_id[i]
  x <- counts[g, ]
  y <- cpm_mat[g, ]
  data.table(
    gene_id = g,
    gene_symbol = target$gene_symbol[i],
    discovery_log2FC = target$discovery_log2FC[i],
    tcga_filterByExpr_retained = unname(keep[g]),
    normal_nonzero_n = sum(x[group == "normal"] > 0),
    tumour_nonzero_n = sum(x[group == "tumor"] > 0),
    normal_count_median = median(x[group == "normal"]),
    tumour_count_median = median(x[group == "tumor"]),
    normal_cpm_median = median(y[group == "normal"]),
    tumour_cpm_median = median(y[group == "tumor"]),
    samples_cpm_ge_1 = sum(y >= 1),
    samples_total = length(y),
    audit_interpretation = "present in the released matrix but removed by the unchanged cohort-wide filterByExpr rule"
  )
})

ans <- rbindlist(rows)
fwrite(ans, file.path(out, "tcga_filtered_discovery_robust_candidates_audit.tsv"), sep = "\t")

review <- c(
  "# TCGA 低表达候选审计",
  "",
  "本审计只解释 5 个在 TCGA 原始矩阵中存在、但被预设 `filterByExpr` 规则过滤的发现稳健基因。",
  "未放宽过滤阈值、未重新拟合被过滤基因，也未根据 ICGC 或发现阶段效应恢复候选。",
  "",
  sprintf("- 审计基因数：%d", nrow(ans)),
  sprintf("- filterByExpr 保留数：%d", sum(ans$tcga_filterByExpr_retained)),
  "- 结论：这些基因在 TCGA 中继续记为 not estimable；低表达过滤不是方向不一致。",
  "",
  "详细的非零样本数、中位计数和 CPM 见 `tcga_filtered_discovery_robust_candidates_audit.tsv`。"
)
writeLines(review, file.path(out, "tcga_filtered_candidates_review_zh.md"))
message("TCGA filtered-candidate audit: PASS")
