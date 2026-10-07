#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(optparse)
})

option_list <- list(make_option("--project-root", dest = "project_root", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)
out <- file.path(root, "results", "validation", "cbc_external_validation_v1")

high_path <- file.path(root, "results", "analysis", "candidate_lock", "candidate_set_v2_high_confidence.tsv")
int_path <- file.path(root, "results", "analysis", "candidate_lock", "candidate_set_v2_intermediate.tsv")
map_path <- file.path(root, "00_protocol", "cbc_icgc_historical_symbol_map_v1.tsv")
raw_path <- file.path(root, "results", "validation", "rnaseq_normalized_v2", "icgc_liri_jp_paired_full_counts.tsv.gz")
icgc_full_path <- file.path(out, "icgc_complete_filtered_transcriptome.tsv.gz")
tcga_path <- file.path(out, "tcga_locked_candidate_multiplicity.tsv")

sha_of <- function(path) strsplit(system2("sha256sum", path, stdout = TRUE), " +")[[1]][1]
stopifnot(
  identical(sha_of(high_path), "55110e282bf9d26619de70d6e8c1e5972d050766b034a6fa06009021ae9fd3cd"),
  identical(sha_of(int_path), "23812a47c8df31301aca8a9061ecd7da56dd41afa5f214a3fc6e4a585494aa59")
)

high <- fread(high_path)
inter <- fread(int_path)
high[, `:=`(gene_id = as.character(gene_id), candidate_order = .I, evidence_tier = "discovery_robust")]
inter[, `:=`(gene_id = as.character(gene_id), candidate_order = .I, evidence_tier = "intermediate")]
candidates <- rbind(
  high[, .(gene_id, gene_symbol, discovery_log2FC = pooled_log2FC, discovery_FDR = FDR, evidence_tier, candidate_order)],
  inter[, .(gene_id, gene_symbol, discovery_log2FC = pooled_log2FC, discovery_FDR = FDR, evidence_tier, candidate_order)]
)

symbol_map <- fread(map_path, colClasses = list(character = c("gene_id", "current_symbol", "icgc_matrix_symbol")))
stopifnot(!anyDuplicated(symbol_map$current_symbol), !anyDuplicated(symbol_map$icgc_matrix_symbol))
stopifnot(all(symbol_map$current_symbol %in% candidates$gene_symbol))
stopifnot(all(symbol_map$gene_id == candidates$gene_id[match(symbol_map$current_symbol, candidates$gene_symbol)]))

raw_symbols <- fread(raw_path, select = "gene_symbol")$gene_symbol
stopifnot(all(symbol_map$icgc_matrix_symbol %in% raw_symbols))

# Verify that every frozen historical symbol has an alias relationship to the
# locked Entrez ID. Ambiguous aliases are allowed only because the frozen map
# chooses one historical symbol before effect estimates are read.
for (i in seq_len(nrow(symbol_map))) {
  rel <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbol_map$icgc_matrix_symbol[i],
    keytype = "ALIAS",
    columns = c("ENTREZID", "SYMBOL")
  )
  if (!symbol_map$gene_id[i] %in% as.character(rel$ENTREZID)) {
    stop("Alias-to-Entrez verification failed for ", symbol_map$current_symbol[i])
  }
}

candidates[, `:=`(
  icgc_matrix_symbol = gene_symbol,
  mapping_status = fifelse(gene_symbol %in% raw_symbols, "current_symbol", "unmapped"),
  mapping_basis = fifelse(gene_symbol %in% raw_symbols, "current symbol present in frozen matrix", "")
)]
candidates[symbol_map, on = .(gene_symbol = current_symbol), `:=`(
  icgc_matrix_symbol = i.icgc_matrix_symbol,
  mapping_status = "historical_symbol_harmonized",
  mapping_basis = i.mapping_basis
)]

icgc_full <- fread(icgc_full_path)
model_cols <- c(
  "gene_symbol", "validation_log2FC", "SE", "nominal_P", "df_total", "mean_logCPM",
  "genomewide_FDR", "validation_CI_low", "validation_CI_high", "cohort"
)
icgc_model <- icgc_full[, ..model_cols]
setnames(icgc_model, "gene_symbol", "icgc_matrix_symbol")

ans <- merge(candidates, icgc_model, by = "icgc_matrix_symbol", all.x = TRUE, sort = FALSE)
setorder(ans, evidence_tier, candidate_order)
ans[, present_before_filter := icgc_matrix_symbol %in% raw_symbols]
ans[, retained_after_filter := is.finite(nominal_P)]
ans[, estimable := present_before_filter & retained_after_filter]
ans[, same_direction := fifelse(estimable, sign(validation_log2FC) == sign(discovery_log2FC), NA)]
ans[, candidate_family_FDR := NA_real_]
for (tier in c("discovery_robust", "intermediate")) {
  idx <- ans$evidence_tier == tier & is.finite(ans$nominal_P)
  ans[idx, candidate_family_FDR := p.adjust(nominal_P, method = "BH")]
}
ans[, candidate_family_directional := estimable & same_direction == TRUE & candidate_family_FDR < 0.05]
ans[, genomewide_directional := estimable & same_direction == TRUE & genomewide_FDR < 0.05]
ans[, estimation_note := fifelse(
  !present_before_filter, "absent_after_identifier_harmonization",
  fifelse(!retained_after_filter, "filtered_low_expression_before_voom", "")
)]
ans[, post_hoc_identifier_harmonization := TRUE]
setcolorder(ans, c(
  "gene_id", "gene_symbol", "icgc_matrix_symbol", "mapping_status", "mapping_basis",
  "evidence_tier", "candidate_order", "discovery_log2FC", "discovery_FDR",
  "present_before_filter", "retained_after_filter", "estimable", "estimation_note",
  "validation_log2FC", "SE", "validation_CI_low", "validation_CI_high", "nominal_P",
  "candidate_family_FDR", "genomewide_FDR", "same_direction",
  "candidate_family_directional", "genomewide_directional", "df_total", "mean_logCPM",
  "cohort", "post_hoc_identifier_harmonization"
))
fwrite(ans, file.path(out, "icgc_locked_candidate_identifier_harmonized_sensitivity.tsv"), sep = "\t", na = "NA")

tcga <- fread(tcga_path)
tcga[, gene_id := as.character(gene_id)]
cross <- merge(
  tcga[, .(
    gene_id, gene_symbol, evidence_tier, candidate_order, discovery_log2FC,
    tcga_estimable = estimable, tcga_log2FC = validation_log2FC,
    tcga_same_direction = same_direction,
    tcga_candidate_family_directional = candidate_family_directional,
    tcga_genomewide_directional = genomewide_directional
  )],
  ans[, .(
    gene_id, gene_symbol, evidence_tier, candidate_order, icgc_matrix_symbol, mapping_status,
    icgc_estimable = estimable, icgc_log2FC = validation_log2FC,
    icgc_same_direction = same_direction,
    icgc_candidate_family_directional = candidate_family_directional,
    icgc_genomewide_directional = genomewide_directional
  )],
  by = c("gene_id", "gene_symbol", "evidence_tier", "candidate_order"), all.x = TRUE, sort = FALSE
)
setorder(cross, evidence_tier, candidate_order)
cross[, jointly_estimable := tcga_estimable == TRUE & icgc_estimable == TRUE]
cross[, dual_same_direction := jointly_estimable & tcga_same_direction == TRUE & icgc_same_direction == TRUE]
cross[, dual_candidate_family_directional := jointly_estimable & tcga_candidate_family_directional == TRUE & icgc_candidate_family_directional == TRUE]
cross[, dual_genomewide_directional := jointly_estimable & tcga_genomewide_directional == TRUE & icgc_genomewide_directional == TRUE]
fwrite(cross, file.path(out, "tcga_icgc_identifier_harmonized_sensitivity.tsv"), sep = "\t", na = "NA")

counts <- rbindlist(lapply(c("discovery_robust", "intermediate"), function(tier) {
  z <- cross[evidence_tier == tier]
  total <- nrow(z)
  joint <- z[jointly_estimable == TRUE, .N]
  data.table(
    evidence_tier = tier,
    metric = c("locked_total", "icgc_estimable", "jointly_estimable", "dual_same_direction", "dual_candidate_family_directional", "dual_genomewide_directional"),
    numerator = c(total, z[icgc_estimable == TRUE, .N], joint, z[dual_same_direction == TRUE, .N], z[dual_candidate_family_directional == TRUE, .N], z[dual_genomewide_directional == TRUE, .N]),
    denominator = c(total, total, total, joint, joint, joint)
  )
}))
counts[, proportion := numerator / denominator]
fwrite(counts, file.path(out, "identifier_harmonization_counts.tsv"), sep = "\t", na = "NA")

review <- c(
  "# ICGC identifier-harmonization sensitivity analysis",
  "",
  "This is a post hoc technical sensitivity analysis. It reuses the frozen ICGC model and its 19,353-gene genome-wide BH family; no expression model or candidate threshold was changed.",
  "",
  paste0("- Historical symbols frozen before effect lookup: ", nrow(symbol_map), "."),
  paste0("- Discovery robust ICGC estimable: ", counts[evidence_tier == "discovery_robust" & metric == "icgc_estimable", numerator], "/173."),
  paste0("- Discovery robust jointly estimable: ", counts[evidence_tier == "discovery_robust" & metric == "jointly_estimable", numerator], "/173."),
  paste0("- Dual same direction: ", counts[evidence_tier == "discovery_robust" & metric == "dual_same_direction", numerator], "/", counts[evidence_tier == "discovery_robust" & metric == "dual_same_direction", denominator], "."),
  paste0("- Dual candidate-family FDR and direction: ", counts[evidence_tier == "discovery_robust" & metric == "dual_candidate_family_directional", numerator], "/", counts[evidence_tier == "discovery_robust" & metric == "dual_candidate_family_directional", denominator], "."),
  paste0("- Dual genome-wide FDR and direction: ", counts[evidence_tier == "discovery_robust" & metric == "dual_genomewide_directional", numerator], "/", counts[evidence_tier == "discovery_robust" & metric == "dual_genomewide_directional", denominator], "."),
  "",
  "Interpretation boundary: this analysis addresses identifier-version missingness only. It is not a new biological cohort and does not establish mechanism or clinical utility."
)
writeLines(review, file.path(out, "identifier_harmonization_review.md"), useBytes = TRUE)
print(counts)

