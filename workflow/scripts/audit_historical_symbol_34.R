#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
stopifnot(length(file_arg) == 1)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
root <- normalizePath(file.path(script_dir, "..", ".."))
out <- file.path(root, "results/validation/medical_oncology_revision_v1")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

map <- read.delim(
  file.path(root, "00_protocol/cbc_icgc_historical_symbol_map_v1.tsv"),
  stringsAsFactors = FALSE,
  colClasses = c(gene_id = "character", current_symbol = "character", icgc_matrix_symbol = "character", mapping_basis = "character")
)
high <- read.delim(file.path(root, "results/analysis/candidate_lock/candidate_set_v2_high_confidence.tsv"), stringsAsFactors = FALSE)
cross_old <- read.delim(
  file.path(root, "results/validation/cbc_external_validation_v1/tcga_icgc_locked_candidate_cross_validation.tsv"),
  stringsAsFactors = FALSE
)
icgc_h <- read.delim(
  file.path(root, "results/validation/cbc_external_validation_v1/icgc_locked_candidate_identifier_harmonized_sensitivity.tsv"),
  stringsAsFactors = FALSE
)
icgc_filtered <- read.delim(
  gzfile(file.path(root, "results/validation/cbc_external_validation_v1/icgc_complete_filtered_transcriptome.tsv.gz")),
  stringsAsFactors = FALSE, colClasses = "character"
)
raw_syms <- unique(icgc_filtered$gene_symbol)

cat(
  "map n=", nrow(map),
  " dup current=", sum(duplicated(map$current_symbol)),
  " dup hist=", sum(duplicated(map$icgc_matrix_symbol)),
  " dup gene_id=", sum(duplicated(map$gene_id)), "\n"
)

hist_overlap_other <- intersect(map$icgc_matrix_symbol, setdiff(high$gene_symbol, map$current_symbol))
cat("historical symbols that are current symbols of other high genes:", paste(hist_overlap_other, collapse = ","), "\n")

old_high <- subset(cross_old, evidence_tier == "discovery_robust")
already <- old_high$gene_symbol[old_high$gene_symbol %in% map$current_symbol & toupper(as.character(old_high$icgc_estimable)) == "TRUE"]
cat("recovered already ICGC-estimable under current symbol n=", length(already), "\n")

safe_select <- function(...) {
  tryCatch(as.data.frame(AnnotationDbi::select(...), stringsAsFactors = FALSE), error = function(e) {
    data.frame(ENTREZID = character(), SYMBOL = character(), stringsAsFactors = FALSE)
  })
}

rows <- list()
for (i in seq_len(nrow(map))) {
  gid <- map$gene_id[i]
  cur <- map$current_symbol[i]
  hist <- map$icgc_matrix_symbol[i]

  rel_hist <- safe_select(org.Hs.eg.db, keys = hist, keytype = "ALIAS", columns = c("ENTREZID", "SYMBOL"))
  entrez_from_hist_alias <- unique(as.character(rel_hist$ENTREZID))
  entrez_from_hist_alias <- entrez_from_hist_alias[!is.na(entrez_from_hist_alias) & nzchar(entrez_from_hist_alias)]
  official_symbols_for_hist <- unique(as.character(rel_hist$SYMBOL))
  official_symbols_for_hist <- official_symbols_for_hist[!is.na(official_symbols_for_hist) & nzchar(official_symbols_for_hist)]

  hist_as_official <- safe_select(org.Hs.eg.db, keys = hist, keytype = "SYMBOL", columns = c("ENTREZID", "SYMBOL"))
  hist_official_ids <- unique(as.character(hist_as_official$ENTREZID))
  hist_official_ids <- hist_official_ids[!is.na(hist_official_ids) & nzchar(hist_official_ids)]

  cur_official <- safe_select(org.Hs.eg.db, keys = gid, keytype = "ENTREZID", columns = "SYMBOL")
  official_now <- if (nrow(cur_official)) as.character(cur_official$SYMBOL[1]) else NA_character_

  in_matrix <- hist %in% raw_syms
  cur_in_matrix <- cur %in% raw_syms
  contains_target <- gid %in% entrez_from_hist_alias
  multi_entrez <- length(entrez_from_hist_alias) > 1

  ambiguity_flags <- character()
  if (!contains_target) ambiguity_flags <- c(ambiguity_flags, "hist_alias_not_linked_to_locked_entrez")
  if (multi_entrez) {
    ambiguity_flags <- c(ambiguity_flags, paste0("hist_alias_maps_to_multiple_entrez:", paste(entrez_from_hist_alias, collapse = "|")))
  }
  other_official <- setdiff(hist_official_ids, gid)
  if (length(other_official)) {
    ambiguity_flags <- c(ambiguity_flags, paste0("hist_is_official_symbol_of_other_entrez:", paste(other_official, collapse = "|")))
  }
  if (!in_matrix) ambiguity_flags <- c(ambiguity_flags, "hist_absent_from_icgc_matrix")
  if (sum(map$icgc_matrix_symbol == hist) > 1) ambiguity_flags <- c(ambiguity_flags, "many_to_one_candidate_collision")
  if (cur_in_matrix) ambiguity_flags <- c(ambiguity_flags, "current_symbol_already_in_matrix_unexpected")
  if (!is.na(official_now) && official_now != cur) {
    ambiguity_flags <- c(ambiguity_flags, paste0("current_symbol_not_orgdb_official:", official_now))
  }

  if (length(ambiguity_flags) == 0 && contains_target && !multi_entrez) {
    verdict <- "UNAMBIGUOUS_SAME_GENE"
  } else if (contains_target && multi_entrez) {
    verdict <- "AMBIGUOUS_MULTI_ENTREZ_ALIAS"
  } else if (length(ambiguity_flags)) {
    verdict <- "NEEDS_REVIEW"
  } else {
    verdict <- "UNAMBIGUOUS_SAME_GENE"
  }

  rows[[i]] <- data.frame(
    gene_id = gid,
    current_symbol = cur,
    icgc_matrix_symbol = hist,
    mapping_basis = map$mapping_basis[i],
    hist_alias_entrez_ids = paste(entrez_from_hist_alias, collapse = "|"),
    n_hist_alias_entrez = length(entrez_from_hist_alias),
    official_symbols_resolved_from_hist = paste(official_symbols_for_hist, collapse = "|"),
    orgdb_official_for_entrez = official_now,
    hist_in_icgc_matrix = in_matrix,
    current_in_icgc_matrix = cur_in_matrix,
    ambiguity_flags = if (length(ambiguity_flags)) paste(ambiguity_flags, collapse = "; ") else "",
    verdict = verdict,
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, rows)
write.table(audit, file.path(out, "historical_symbol_34_ambiguity_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(table(audit$verdict))
cat("\nNON-UNAMBIGUOUS:\n")
print(subset(audit, verdict != "UNAMBIGUOUS_SAME_GENE")[, c("current_symbol", "icgc_matrix_symbol", "n_hist_alias_entrez", "hist_alias_entrez_ids", "ambiguity_flags", "verdict")])
cat("\nDocumented selected-over-ambiguous-basis rows:\n")
print(subset(audit, grepl("ambiguous", mapping_basis, ignore.case = TRUE))[, c("current_symbol", "icgc_matrix_symbol", "mapping_basis", "n_hist_alias_entrez", "hist_alias_entrez_ids", "verdict", "ambiguity_flags")])

dr <- subset(icgc_h, evidence_tier == "discovery_robust" & toupper(as.character(estimable)) == "TRUE")
cat(
  "\nEstimable discovery_robust:", nrow(dr),
  " unique matrix symbols:", length(unique(dr$icgc_matrix_symbol)),
  " unique gene_id:", length(unique(dr$gene_id)),
  " unique current:", length(unique(dr$gene_symbol)), "\n"
)
cat("dup matrix among estimable:", sum(duplicated(dr$icgc_matrix_symbol)), "\n")
cat("dup gene_id among estimable:", sum(duplicated(dr$gene_id)), "\n")

old_missing <- old_high$gene_symbol[toupper(as.character(old_high$icgc_estimable)) == "FALSE"]
cat("old ICGC non-estimable:", length(old_missing), "\n")
cat("recovered intersection with old missing:", length(intersect(old_missing, map$current_symbol)), "\n")
cat("recovered not in old missing:", paste(setdiff(map$current_symbol, old_missing), collapse = ","), "\n")
cat("old missing not in 34-map:", paste(setdiff(old_missing, map$current_symbol), collapse = ","), "\n")
