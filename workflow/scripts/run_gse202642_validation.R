#!/usr/bin/env Rscript
# GSE202642 patient/tissue-level single-cell localisation for locked candidate set v2.
# Cells are used only for localisation; statistical units are tissue samples (GSM).

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(limma)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(irlba)
  library(igraph)
  library(uwot)
  library(ggplot2)
  library(optparse)
})

option_list <- list(make_option("--project-root", dest = "project_root", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)

data_dir <- file.path(root, "data", "external", "gse202642")
out_dir <- file.path(root, "results", "validation", "gse202642")
fig_dir <- file.path(out_dir, "figures")
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
  high[, .(gene_id, gene_symbol, discovery_log2FC = pooled_log2FC, evidence_tier = "high_confidence")],
  inter[, .(gene_id, gene_symbol, discovery_log2FC = pooled_log2FC, evidence_tier = "intermediate_confidence")]
)
stopifnot(nrow(candidates) == 292L)

manifest <- fread(file.path(out_dir, "sample_manifest.tsv"))
manifest[, barcode_suffix := as.character(barcode_suffix)]

mtx_path <- file.path(data_dir, "GSE202642_matrix.mtx.gz")
bc_path <- file.path(data_dir, "GSE202642_barcodes.tsv.gz")
ft_path <- file.path(data_dir, "GSE202642_features.tsv.gz")
stopifnot(file.exists(mtx_path), file.exists(bc_path), file.exists(ft_path))

message("Reading GSE202642 matrix (genes x cells)")
features <- fread(cmd = paste("zcat", shQuote(ft_path)), header = FALSE,
                  col.names = c("ensembl_id", "gene_symbol", "feature_type"))
barcodes <- fread(cmd = paste("zcat", shQuote(bc_path)), header = FALSE, col.names = "barcode")
mat <- readMM(gzfile(mtx_path))
mat <- as(mat, "dgCMatrix")
if (nrow(mat) != nrow(features) || ncol(mat) != nrow(barcodes)) {
  # 10x MTX is often genes x cells; verify orientation
  if (nrow(mat) == nrow(barcodes) && ncol(mat) == nrow(features)) {
    mat <- t(mat)
  } else {
    stop("Matrix dimensions do not match features/barcodes")
  }
}
rownames(mat) <- make.unique(features$ensembl_id)
colnames(mat) <- barcodes$barcode

# Map barcode suffix to GSM sample
suffix <- sub("^.*-", "", barcodes$barcode)
meta <- data.table(
  cell_barcode = barcodes$barcode,
  barcode_suffix = suffix
)
meta <- merge(meta, manifest, by = "barcode_suffix", all.x = TRUE, sort = FALSE)
setkey(meta, cell_barcode)
meta <- meta[colnames(mat)]
stopifnot(!anyNA(meta$gsm_id))

# QC
mito <- grepl("^MT-", features$gene_symbol, ignore.case = TRUE) |
  grepl("^MT-", rownames(mat), ignore.case = TRUE)
# Ensembl MT genes: use symbol
mito <- features$gene_symbol %in% c(
  "MT-ND1","MT-ND2","MT-CO1","MT-CO2","MT-ATP8","MT-ATP6","MT-CO3",
  "MT-ND3","MT-ND4L","MT-ND4","MT-ND5","MT-ND6","MT-CYB"
) | grepl("^MT-", features$gene_symbol)
nCount <- Matrix::colSums(mat)
nFeature <- Matrix::colSums(mat > 0)
pct_mito <- Matrix::colSums(mat[mito, , drop = FALSE]) / pmax(nCount, 1) * 100
keep <- nCount >= 500 & nFeature >= 200 & nFeature <= 6000 & pct_mito <= 20
message(sprintf("QC retain %d / %d cells", sum(keep), length(keep)))
mat <- mat[, keep, drop = FALSE]
meta <- meta[keep]
meta[, `:=`(nCount = nCount[keep], nFeature = nFeature[keep], pct_mito = pct_mito[keep])]

# Gene filter: detected in >= 20 cells
gene_keep <- Matrix::rowSums(mat > 0) >= 20
mat <- mat[gene_keep, , drop = FALSE]
features_f <- features[gene_keep]

# Log-normalize (Seurat-like)
lib <- Matrix::colSums(mat)
scale_factor <- 1e4
mat_norm <- mat
mat_norm@x <- log1p(mat_norm@x / rep(lib, diff(mat_norm@p)) * scale_factor)

# HVG by variance of log-normalised counts
set.seed(20260916)
rs <- Matrix::rowSums(mat_norm)
rss <- Matrix::rowSums(mat_norm^2)
n <- ncol(mat_norm)
vars <- pmax(as.numeric(rss) / n - (as.numeric(rs) / n)^2, 0)
hvg_n <- min(2000L, nrow(mat_norm))
hvg <- order(vars, decreasing = TRUE)[seq_len(hvg_n)]
mat_hvg <- mat_norm[hvg, , drop = FALSE]

n_cells <- ncol(mat_hvg)
message(sprintf("PCA on %d HVG x %d cells", nrow(mat_hvg), n_cells))
X <- as.matrix(mat_hvg)
X <- X - rowMeans(X)
pcs <- irlba(t(X), nv = 30)
pca <- pcs$u %*% diag(pcs$d)
rownames(pca) <- colnames(mat_hvg)
rm(X); gc()

# kNN + Louvain
message("Building kNN / Louvain clusters")
if (requireNamespace("FNN", quietly = TRUE)) {
  nn_idx <- FNN::get.knn(pca, k = 20)$nn.index
} else {
  # Approximate neighbors via uwot
  nn_idx <- uwot::umap(pca, n_neighbors = 20, ret_nn = TRUE, n_epochs = 1)$nn$euclidean$idx[, -1, drop = FALSE]
}
edges <- rbind(
  cbind(rep(seq_len(nrow(nn_idx)), each = ncol(nn_idx)), as.integer(t(nn_idx))),
  cbind(as.integer(t(nn_idx)), rep(seq_len(nrow(nn_idx)), each = ncol(nn_idx)))
)
g <- graph_from_edgelist(edges, directed = FALSE)
g <- simplify(g)
cl <- cluster_louvain(g)
meta$cluster <- as.character(membership(cl))

# Marker-based broad annotation
marker_sets <- list(
  Hepatocyte = c("ALB","APOA1","APOA2","CYP3A4","HP","TTR"),
  Malignant_hepatocyte = c("AFP","GPC3","MDK","SPINK1","AKR1B10"),
  Cholangiocyte = c("KRT19","EPCAM","SOX9","CFTR"),
  Endothelial = c("PECAM1","VWF","CLEC4G","FCN2","STAB2"),
  Fibroblast_CAF = c("COL1A1","COL1A2","DCN","ACTA2","PDGFRA","LUM"),
  Kupffer_macrophage = c("CD68","CD163","MARCO","VSIG4","C1QA"),
  T_cell = c("CD3D","CD3E","CD2","TRAC","IL7R"),
  NK_cell = c("NKG7","GNLY","KLRD1","KLRF1"),
  B_cell = c("MS4A1","CD79A","CD79B","IGKC"),
  Plasma = c("JCHAIN","MZB1","SSR4","XBP1"),
  Dendritic = c("CLEC9A","XCR1","CD1C","FCER1A"),
  Neutrophil = c("FCGR3B","CSF3R","S100A8","S100A9"),
  Cycling = c("MKI67","TOP2A","PCNA","STMN1")
)

# Build symbol -> row index on normalised matrix
sym <- features_f$gene_symbol
names(sym) <- rownames(mat_norm)
# mean expression per cluster for markers
cluster_ids <- sort(unique(meta$cluster))
score_dt <- rbindlist(lapply(names(marker_sets), function(ct) {
  genes <- marker_sets[[ct]]
  rows <- which(features_f$gene_symbol %in% genes)
  if (!length(rows)) return(NULL)
  sub <- mat_norm[rows, , drop = FALSE]
  data.table(
    cluster = meta$cluster,
    score = Matrix::colMeans(sub)
  )[, .(cell_type = ct, mean_score = mean(score)), by = cluster]
}))
# assign each cluster to max score cell type
assign <- score_dt[, .SD[which.max(mean_score)], by = cluster]
setnames(assign, "cell_type", "cell_type_broad")
meta <- merge(meta, assign[, .(cluster, cell_type_broad)], by = "cluster", all.x = TRUE, sort = FALSE)
setkey(meta, cell_barcode)
meta <- meta[colnames(mat)]

# UMAP
message("UMAP")
um <- umap(pca, n_neighbors = 30, min_dist = 0.3, metric = "euclidean", verbose = FALSE)
meta$UMAP1 <- um[, 1]
meta$UMAP2 <- um[, 2]

# Cell-type marker table (top genes by cluster mean vs others - simple)
message("Cluster markers")
# Pseudobulk by sample x cell type for candidates
# Map candidates Ensembl/Entrez to matrix rows via symbol and entrez
map_entrez <- tryCatch(
  mapIds(org.Hs.eg.db, keys = features_f$ensembl_id, keytype = "ENSEMBL", column = "ENTREZID", multiVals = "first"),
  error = function(e) rep(NA_character_, nrow(features_f))
)
features_f$gene_id <- unname(map_entrez)
cand_rows <- match(candidates$gene_id, features_f$gene_id)
# also try symbol
miss <- is.na(cand_rows)
cand_rows[miss] <- match(candidates$gene_symbol[miss], features_f$gene_symbol)

# Localization: mean log-norm expression and pct expressing by cell type (all tissues)
loc_list <- list()
pb_list <- list()
for (i in seq_len(nrow(candidates))) {
  gi <- candidates$gene_id[i]
  gs <- candidates$gene_symbol[i]
  ri <- cand_rows[i]
  if (is.na(ri)) {
    loc_list[[i]] <- data.table(
      gene_id = gi, gene_symbol = gs, evidence_tier = candidates$evidence_tier[i],
      cell_type_broad = NA_character_, n_cells = 0L, mean_lognorm = NA_real_,
      pct_expressing = NA_real_, status = "not_detected_in_matrix"
    )
    next
  }
  expr <- mat_norm[ri, ]
  raw <- mat[ri, ]
  tmp <- data.table(
    cell_type_broad = meta$cell_type_broad,
    gsm_id = meta$gsm_id,
    tissue_class = meta$tissue_class,
    patient_label = meta$patient_label,
    expr = as.numeric(expr),
    det = as.numeric(raw > 0)
  )
  loc <- tmp[, .(
    n_cells = .N,
    mean_lognorm = mean(expr),
    pct_expressing = mean(det) * 100
  ), by = cell_type_broad]
  loc[, `:=`(gene_id = gi, gene_symbol = gs, evidence_tier = candidates$evidence_tier[i], status = "ok")]
  loc_list[[i]] <- loc

  # patient/tissue-level pseudobulk: mean lognorm within sample x cell type
  pb <- tmp[, .(
    n_cells = .N,
    pseudobulk_mean_lognorm = mean(expr),
    pct_expressing = mean(det) * 100
  ), by = .(gsm_id, tissue_class, patient_label, cell_type_broad)]
  pb[, `:=`(gene_id = gi, gene_symbol = gs, evidence_tier = candidates$evidence_tier[i])]
  pb_list[[i]] <- pb
}
localization <- rbindlist(loc_list, fill = TRUE)
pseudobulk <- rbindlist(pb_list, fill = TRUE)

# Dominant cell type per gene (max mean among types with >=50 cells)
dom <- localization[status == "ok" & n_cells >= 50][
  , .SD[which.max(mean_lognorm)], by = .(gene_id, gene_symbol, evidence_tier)]
setnames(dom, c("cell_type_broad", "mean_lognorm", "pct_expressing"),
         c("dominant_cell_type", "dominant_mean_lognorm", "dominant_pct_expressing"))

# Directional tumour vs adjacent support within cell type: sample-level units
# Require >=3 sample units per tissue_class
dir_rows <- list()
for (ct in sort(unique(meta$cell_type_broad))) {
  for (i in seq_len(nrow(candidates))) {
    gi <- candidates$gene_id[i]
    pb_g <- pseudobulk[gene_id == gi & cell_type_broad == ct & n_cells >= 10]
    n_t <- uniqueN(pb_g[tissue_class == "tumor"]$gsm_id)
    n_a <- uniqueN(pb_g[tissue_class == "adjacent"]$gsm_id)
    if (n_t < 3 || n_a < 3) {
      dir_rows[[length(dir_rows) + 1L]] <- data.table(
        gene_id = gi, gene_symbol = candidates$gene_symbol[i],
        evidence_tier = candidates$evidence_tier[i], cell_type_broad = ct,
        n_tumor_samples = n_t, n_adjacent_samples = n_a,
        mean_tumor = NA_real_, mean_adjacent = NA_real_,
        log2FC_tumor_minus_adjacent = NA_real_,
        wilcox_p = NA_real_, status = "insufficient_sample_units"
      )
      next
    }
    t_vals <- pb_g[tissue_class == "tumor", pseudobulk_mean_lognorm]
    a_vals <- pb_g[tissue_class == "adjacent", pseudobulk_mean_lognorm]
    wt <- suppressWarnings(wilcox.test(t_vals, a_vals, exact = FALSE))
    dir_rows[[length(dir_rows) + 1L]] <- data.table(
      gene_id = gi, gene_symbol = candidates$gene_symbol[i],
      evidence_tier = candidates$evidence_tier[i], cell_type_broad = ct,
      n_tumor_samples = n_t, n_adjacent_samples = n_a,
      mean_tumor = mean(t_vals), mean_adjacent = mean(a_vals),
      log2FC_tumor_minus_adjacent = mean(t_vals) - mean(a_vals),
      wilcox_p = wt$p.value, status = "tested"
    )
  }
}
directional <- rbindlist(dir_rows, fill = TRUE)
directional[status == "tested", fdr := p.adjust(wilcox_p, method = "BH")]
directional[status != "tested", fdr := NA_real_]
directional <- merge(
  directional,
  candidates[, .(gene_id, discovery_log2FC)],
  by = "gene_id", all.x = TRUE
)
directional[, same_direction := {
  ifelse(status != "tested", NA,
         sign(log2FC_tumor_minus_adjacent) == sign(discovery_log2FC) &
           abs(log2FC_tumor_minus_adjacent) > 0)
}]

# Cluster marker table (top 5 symbols by mean in cluster)
marker_tbl <- rbindlist(lapply(cluster_ids, function(cl0) {
  cells <- meta$cluster == cl0
  if (sum(cells) < 20) return(NULL)
  mu_in <- Matrix::rowMeans(mat_norm[, cells, drop = FALSE])
  mu_out <- Matrix::rowMeans(mat_norm[, !cells, drop = FALSE])
  diff <- mu_in - mu_out
  top <- order(diff, decreasing = TRUE)[seq_len(min(10L, length(diff)))]
  data.table(
    cluster = cl0,
    cell_type_broad = meta[cluster == cl0, cell_type_broad][1],
    n_cells = sum(cells),
    gene_symbol = features_f$gene_symbol[top],
    ensembl_id = features_f$ensembl_id[top],
    mean_in = mu_in[top],
    mean_out = mu_out[top],
    diff = diff[top]
  )
}))

# Counts summary
counts <- data.table(
  metric = c(
    "n_raw_cells", "n_qc_cells", "n_tumor_samples", "n_adjacent_samples",
    "n_cell_types", "n_candidates", "n_candidates_detected",
    "independent_unit", "min_units_rule"
  ),
  value = c(
    as.character(length(keep)), as.character(sum(keep)),
    as.character(uniqueN(meta[tissue_class == "tumor"]$gsm_id)),
    as.character(uniqueN(meta[tissue_class == "adjacent"]$gsm_id)),
    as.character(uniqueN(meta$cell_type_broad)),
    "292", as.character(sum(!is.na(cand_rows))),
    "tissue_sample_GSM_not_cell", "3_samples_per_tissue_class_within_cell_type"
  )
)

fwrite(meta, file.path(out_dir, "cell_metadata.tsv.gz"), sep = "\t")
fwrite(marker_tbl, file.path(out_dir, "cell_type_marker_table.tsv"), sep = "\t")
fwrite(localization, file.path(out_dir, "candidate_celltype_localization.tsv"), sep = "\t")
fwrite(pseudobulk, file.path(out_dir, "candidate_patient_pseudobulk.tsv.gz"), sep = "\t")
fwrite(dom, file.path(out_dir, "candidate_dominant_celltype.tsv"), sep = "\t")
fwrite(directional, file.path(out_dir, "candidate_celltype_tumor_adjacent_pseudobulk.tsv"), sep = "\t")
fwrite(counts, file.path(out_dir, "validation_counts.tsv"), sep = "\t")
fwrite(assign, file.path(out_dir, "cluster_to_celltype.tsv"), sep = "\t")

# Figures
p1 <- ggplot(meta, aes(UMAP1, UMAP2, color = cell_type_broad)) +
  geom_point(size = 0.1, alpha = 0.5) + theme_bw(base_size = 11) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "GSE202642 UMAP by broad cell type", color = "Cell type")
ggsave(file.path(fig_dir, "umap_celltype.png"), p1, width = 9, height = 6, dpi = 150)

p2 <- ggplot(meta, aes(UMAP1, UMAP2, color = tissue_class)) +
  geom_point(size = 0.1, alpha = 0.5) + theme_bw(base_size = 11) +
  labs(title = "GSE202642 UMAP by tissue class")
ggsave(file.path(fig_dir, "umap_tissue_class.png"), p2, width = 8, height = 6, dpi = 150)

# Feature plots for top few high-confidence dominant genes
top_plot <- head(dom[evidence_tier == "high_confidence"][order(-dominant_mean_lognorm)], 6)
for (j in seq_len(nrow(top_plot))) {
  ri <- cand_rows[match(top_plot$gene_id[j], candidates$gene_id)]
  if (is.na(ri)) next
  df <- data.table(UMAP1 = meta$UMAP1, UMAP2 = meta$UMAP2, expr = as.numeric(mat_norm[ri, ]))
  pj <- ggplot(df, aes(UMAP1, UMAP2, color = expr)) +
    geom_point(size = 0.1) +
    scale_color_gradient(low = "grey90", high = "#b2182b") +
    theme_bw() +
    labs(title = paste0("Feature: ", top_plot$gene_symbol[j]))
  ggsave(file.path(fig_dir, paste0("feature_", top_plot$gene_symbol[j], ".png")),
         pj, width = 7, height = 5, dpi = 120)
}

# Dominant cell type bar for high confidence
dom_high <- dom[evidence_tier == "high_confidence", .N, by = dominant_cell_type][order(-N)]
p3 <- ggplot(dom_high, aes(reorder(dominant_cell_type, N), N)) +
  geom_col(fill = "#2c7fb8") + coord_flip() + theme_bw() +
  labs(title = "High-confidence candidates: dominant cell type", x = NULL, y = "N genes")
ggsave(file.path(fig_dir, "high_confidence_dominant_celltype.png"), p3, width = 8, height = 5, dpi = 150)

# Source data
fwrite(dom_high, file.path(fig_dir, "high_confidence_dominant_celltype_source.tsv"), sep = "\t")

message("Done GSE202642 localisation")
message(paste(capture.output(print(counts)), collapse = "\n"))
