args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep('^--file=', args, value = TRUE)
if (length(file_arg)) {
  script_dir <- dirname(normalizePath(sub('^--file=', '', file_arg[1]), winslash = '/', mustWork = TRUE))
} else {
  script_dir <- normalizePath(getwd(), winslash = '/', mustWork = TRUE)
}
base_dir <- normalizePath(file.path(script_dir, '..'), winslash = '/', mustWork = TRUE)
out_fig <- file.path(base_dir, 'figures')
pipeline_demo_dir <- normalizePath(file.path(base_dir, '..', '2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE', 'EXAMPLE_FILES', 'OUTPUT_DEMO'),
                                   winslash = '/', mustWork = TRUE)
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

infile <- file.path(pipeline_demo_dir, 'summary_terap_empirical_postfit_metrics.csv')
if (!file.exists(infile)) stop('Missing input file: ', infile)
d <- read.csv(infile, stringsAsFactors = FALSE)
ord <- c('chronos_discrete', 'chronos_clock', 'chronos_correlated', 'treePL', 'chronos_relaxed')
d <- d[match(ord, d$candidate), ]
labels <- c('Discrete', 'Clock', 'Correlated', 'treePL', 'Relaxed')
cols <- c('#1b9e77', '#2c7fb8', '#d95f0e', '#6baed6', '#7570b3')
gap_label <- unique(d$gap_mode)
gap_label <- gsub('_', '-', gap_label)

png(file.path(out_fig, 'postfit_metric_family_values.png'), width = 2200, height = 1300, res = 170)
layout(matrix(c(1,2,3,4,5,5), nrow = 2, byrow = TRUE))
par(mar = c(8,5,4,1), oma = c(0.2, 0.2, 2.2, 0.2))
plot_panel <- function(vals, ttl, ylab, note) {
  bp <- barplot(vals, names.arg = labels, col = cols, las = 2,
                ylab = ylab, main = ttl,
                ylim = c(0, max(vals, na.rm = TRUE) * 1.22))
  text(bp, vals, labels = sprintf('%.3f', vals), pos = 3, cex = 0.85)
  mtext(note, side = 1, line = 6.4, cex = 0.9)
}
plot_panel(d$pulse_default_selector_error, 'Pulse preservation (default)', 'Selector error', 'Lower is better')
plot_panel(d$pulse_burst_selector_error, 'Pulse preservation (burst)', 'Selector error', 'Lower is better')
plot_panel(d$fossil_gap_burden, paste0('Gap burden (', gap_label, ')'), 'Gap burden', 'Lower is better')
plot_panel(d$rate_irregularity, 'Rate plausibility', 'Rate irregularity', 'Lower is better')

par(mar = c(8,5,4,1))
bp <- barplot(d$rank_mean_4families, names.arg = labels, col = cols, las = 2,
              ylab = 'Mean rank across 4 families',
              main = 'Overall post-fit rank',
              ylim = c(0, max(d$rank_mean_4families, na.rm = TRUE) * 1.25))
text(bp, d$rank_mean_4families,
     labels = paste0(sprintf('%.2f', d$rank_mean_4families), ' (rank ', d$rank_mean_4families_rank, ')'),
     pos = 3, cex = 0.85)
mtext('Lower is better', side = 1, line = 6.4, cex = 0.9)
mtext('Terapontoid example: post-fit evaluation metrics across chronos models and treePL', side = 3, outer = TRUE, line = 0.5, cex = 1.6, font = 2)
dev.off()
message('Done. Wrote: ', normalizePath(file.path(out_fig, 'postfit_metric_family_values.png'), winslash = '/'))
