suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(make_option("--project-root", dest = "project_root", type = "character"))
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$project_root)) stop("--project-root is required")
root <- normalizePath(opt$project_root, mustWork = TRUE)
sensitivity_root <- file.path(root, "results", "analysis", "sensitivity")
figure_root <- file.path(sensitivity_root, "figures")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

fixed <- fread(file.path(sensitivity_root, "fixed_effect_candidate_comparison.tsv"))
all_samples <- fread(file.path(sensitivity_root, "all_samples", "all_samples_candidate_comparison.tsv"))
loo <- fread(file.path(sensitivity_root, "leave_one_cohort_out", "candidate_leave_one_out_long.tsv.gz"))

fixed_plot <- fixed[, .(gene_id, primary = random_log2FC, sensitivity = fixed_log2FC,
                        passes = fixed_passes_primary_rule, analysis = "Fixed effect")]
all_plot <- all_samples[, .(gene_id, primary = primary_log2FC, sensitivity = all_samples_log2FC,
                            passes = all_samples_passes_primary_rule, analysis = "All 438 samples")]
comparison <- rbind(fixed_plot, all_plot)
labels <- comparison[, .(
  label = sprintf("Pearson r = %.3f\nSame direction = %d/%d",
                  cor(primary, sensitivity, use = "complete.obs"),
                  sum(sign(primary) == sign(sensitivity), na.rm = TRUE), .N)
), by = analysis]

scatter <- ggplot(comparison, aes(primary, sensitivity, color = passes)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.45, color = "#555555") +
  geom_point(size = 1.25, alpha = 0.75) +
  geom_text(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.15, size = 3.2) +
  facet_wrap(~analysis, nrow = 1) +
  scale_color_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#BDBDBD"),
                     labels = c(`TRUE` = "Passes sensitivity rule", `FALSE` = "Does not pass")) +
  labs(x = "Primary random-effects log2FC", y = "Sensitivity-analysis log2FC", color = NULL) +
  coord_equal() + theme_classic(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(figure_root, "candidate_effect_concordance.png"), scatter,
       width = 9.4, height = 4.7, dpi = 320)

loo_by_cohort <- loo[, .(
  `Same direction` = sum(same_direction_as_primary),
  `FDR < 0.05` = sum(FDR < 0.05),
  `Effect threshold` = sum(abs(pooled_log2FC) >= log2(1.5)),
  `Complete rule` = sum(passes_primary_rule)
), by = omitted_cohort]
loo_long <- melt(loo_by_cohort, id.vars = "omitted_cohort", variable.name = "criterion",
                 value.name = "candidates")
loo_long[, omitted_cohort := factor(omitted_cohort, levels = unique(loo_by_cohort$omitted_cohort))]
loo_bar <- ggplot(loo_long, aes(omitted_cohort, candidates, fill = criterion)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_hline(yintercept = 446, color = "#555555", linewidth = 0.35, linetype = 2) +
  scale_fill_manual(values = c("Same direction" = "#2166AC", "FDR < 0.05" = "#FDAE61",
                               "Effect threshold" = "#66BD63", "Complete rule" = "#D73027")) +
  labs(x = "Omitted cohort", y = "Locked candidates (of 446)", fill = NULL) +
  theme_classic(base_size = 10.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")
ggsave(file.path(figure_root, "leave_one_out_candidate_stability.png"), loo_bar,
       width = 9.2, height = 5.4, dpi = 320)

pass_distribution <- loo[, .(scenarios_passing = sum(passes_primary_rule)), by = gene_id][, .N, by = scenarios_passing]
pass_distribution <- merge(data.table(scenarios_passing = 0:6), pass_distribution,
                           by = "scenarios_passing", all.x = TRUE)
pass_distribution[is.na(N), N := 0L]
fwrite(pass_distribution, file.path(sensitivity_root, "leave_one_cohort_out", "candidate_pass_count_distribution.tsv"), sep = "\t")
pass_plot <- ggplot(pass_distribution, aes(factor(scenarios_passing), N)) +
  geom_col(fill = "#4C78A8", width = 0.72) +
  geom_text(aes(label = N), vjust = -0.25, size = 3.4) +
  labs(x = "LOCO scenarios passing the complete rule (of 6)", y = "Locked candidates") +
  theme_classic(base_size = 11) + expand_limits(y = max(pass_distribution$N) * 1.08)
ggsave(file.path(figure_root, "leave_one_out_pass_distribution.png"), pass_plot,
       width = 6.6, height = 4.6, dpi = 320)
