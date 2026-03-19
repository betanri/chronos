args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(script_dir)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

method_csv <- "data/pulse_metric_720_method_mean_summary.csv"
wins_csv <- "data/pulse_metric_720_wins_table.csv"
clock_csv <- "data/pulse_metric_720_delta_by_clock_model.csv"
paired_csv <- "data/pulse_metric_720_condition_paired_treepl_vs_chronos.csv"

if (!file.exists(method_csv)) stop("Missing input file: ", method_csv)
if (!file.exists(wins_csv)) stop("Missing input file: ", wins_csv)
if (!file.exists(clock_csv)) stop("Missing input file: ", clock_csv)
if (!file.exists(paired_csv)) stop("Missing input file: ", paired_csv)

meth <- read.csv(method_csv, stringsAsFactors = FALSE)
wins <- read.csv(wins_csv, stringsAsFactors = FALSE)
clk <- read.csv(clock_csv, stringsAsFactors = FALSE)
paired <- read.csv(paired_csv, stringsAsFactors = FALSE)

# Figure 6: win counts across representative-condition comparisons.
png("figures/fig6_pulse_rate_win_counts.png", width = 1700, height = 950, res = 160)
par(mar = c(10, 6, 4, 2))
mat <- rbind(wins$chronos_better_conditions, wins$treePL_better_conditions)
colnames(mat) <- c("pulse score", "burst loss", "tempo composite", "MAE", "rate irregularity")
barplot(mat,
        beside = TRUE,
        col = c("#74c476", "#6baed6"),
        ylim = c(0, max(mat) * 1.2),
        las = 2,
        ylab = "Conditions won (24 representative conditions)",
        main = "720 benchmark: representative-condition wins by metric")
legend("topright", legend = c("chronos", "treePL"), fill = c("#74c476", "#6baed6"), bty = "n")
mtext("Pulse and rate metrics use one representative tree per simulation condition; MAE uses all replicate values.", side = 1, line = 8, cex = 0.95)
dev.off()

# Figure 7: method means for pulse and rate summaries.
png("figures/fig7_pulse_rate_method_means.png", width = 1800, height = 1000, res = 160)
layout(matrix(1:4, nrow = 1))
par(mar = c(7, 5, 4, 1))
plot_pair <- function(col_name, main, ylab, better) {
  vals <- meth[[col_name]]
  names(vals) <- meth$method
  cols <- c(chronos = "#74c476", treePL = "#6baed6")
  bp <- barplot(vals[c("chronos", "treePL")],
                col = cols[c("chronos", "treePL")],
                ylab = ylab,
                main = main,
                ylim = c(0, max(vals, na.rm = TRUE) * 1.25),
                las = 2)
  text(bp, vals[c("chronos", "treePL")], labels = sprintf("%.3f", vals[c("chronos", "treePL")]), pos = 3, cex = 0.9)
  mtext(better, side = 1, line = 5, cex = 0.9)
}
plot_pair("pulse_score", "Pulse score", "Mean value", "Higher is better")
plot_pair("pulse_error", "Pulse error", "Mean value", "Lower is better")
plot_pair("burst_loss", "Burst loss", "Mean value", "Lower is better")
plot_pair("rate_irregularity", "Rate irregularity", "Mean value", "Lower is better")
mtext("720 benchmark: method means under pulse-preservation and rate-plausibility metrics", side = 3, outer = TRUE, line = -1.5, cex = 1.4)
dev.off()

# Figure 8: where treePL still wins by true clock model.
count_metric <- function(x, higher_better = FALSE) {
  if (higher_better) {
    c(treePL = sum(x > 0, na.rm = TRUE),
      chronos = sum(x < 0, na.rm = TRUE),
      ties = sum(x == 0, na.rm = TRUE))
  } else {
    c(treePL = sum(x < 0, na.rm = TRUE),
      chronos = sum(x > 0, na.rm = TRUE),
      ties = sum(x == 0, na.rm = TRUE))
  }
}

metric_specs <- list(
  pulse_score = list(col = "delta_pulse_score_treepl_minus_chronos", higher = TRUE),
  burst_loss = list(col = "delta_burst_loss_treepl_minus_chronos", higher = FALSE),
  tempo_composite = list(col = "delta_tempo_composite_treepl_minus_chronos", higher = FALSE),
  mae_mean30 = list(col = "delta_mae_mean30_treepl_minus_chronos", higher = FALSE),
  rate_irregularity = list(col = "delta_rate_irregularity_treepl_minus_chronos", higher = FALSE)
)

clock_levels <- c("strict", "autocorrelated", "independent", "discrete")
win_rows <- list()
adv_mat <- matrix(NA_real_, nrow = length(metric_specs), ncol = length(clock_levels),
                  dimnames = list(names(metric_specs), clock_levels))

for (m in names(metric_specs)) {
  spec <- metric_specs[[m]]
  for (clk_name in clock_levels) {
    x <- paired[paired$clock_model == clk_name, spec$col]
    cc <- count_metric(x, spec$higher)
    win_rows[[length(win_rows) + 1L]] <- data.frame(
      metric = m,
      clock_model = clk_name,
      treePL_wins = cc["treePL"],
      chronos_wins = cc["chronos"],
      ties = cc["ties"],
      stringsAsFactors = FALSE
    )
  }
}
clock_wins <- do.call(rbind, win_rows)
write.csv(clock_wins, "data/pulse_metric_720_clock_model_wins.csv", row.names = FALSE)

# Positive means chronos better on average; negative means treePL better on average.
adv_mat["pulse_score", ] <- -clk$delta_pulse_score_treepl_minus_chronos[match(clock_levels, clk$clock_model)]
adv_mat["burst_loss", ] <- clk$delta_burst_loss_treepl_minus_chronos[match(clock_levels, clk$clock_model)]
adv_mat["tempo_composite", ] <- clk$delta_tempo_composite_treepl_minus_chronos[match(clock_levels, clk$clock_model)]
adv_mat["mae_mean30", ] <- clk$delta_mae_mean30_treepl_minus_chronos[match(clock_levels, clk$clock_model)]
adv_mat["rate_irregularity", ] <- clk$delta_rate_irregularity_treepl_minus_chronos[match(clock_levels, clk$clock_model)]

win_mat <- matrix(NA_real_, nrow = length(metric_specs), ncol = length(clock_levels),
                  dimnames = list(names(metric_specs), clock_levels))
for (m in names(metric_specs)) {
  for (clk_name in clock_levels) {
    win_mat[m, clk_name] <- clock_wins$treePL_wins[
      clock_wins$metric == m & clock_wins$clock_model == clk_name
    ][1]
  }
}

png("figures/fig8_postfit_by_clock_model.png", width = 2000, height = 1100, res = 160)
layout(matrix(c(1, 2), nrow = 1), widths = c(1.25, 1))
par(mar = c(7, 8, 4, 2))
zmax <- max(abs(adv_mat), na.rm = TRUE)
adv_cols <- colorRampPalette(c("#2b8cbe", "#f7f7f7", "#238b45"))(200)
image(1:ncol(adv_mat), 1:nrow(adv_mat), t(adv_mat[nrow(adv_mat):1, ]),
      col = adv_cols, zlim = c(-zmax, zmax), xaxt = "n", yaxt = "n",
      xlab = "True clock model", ylab = "", main = "Mean advantage by true clock model")
axis(1, at = 1:ncol(adv_mat), labels = colnames(adv_mat), las = 2)
axis(2, at = 1:nrow(adv_mat), labels = rev(c("pulse score", "burst loss", "tempo composite", "MAE", "rate irregularity")), las = 2)
for (i in seq_len(nrow(adv_mat))) {
  for (j in seq_len(ncol(adv_mat))) {
    text(j, nrow(adv_mat) - i + 1, labels = sprintf("%.2f", adv_mat[i, j]), cex = 0.85)
  }
}
mtext("Positive = chronos better on average; negative = treePL better on average", side = 1, line = 5.5, cex = 0.9)

par(mar = c(7, 7, 4, 2))
win_cols <- colorRampPalette(c("#f7fbff", "#9ecae1", "#08519c"))(7)
image(1:ncol(win_mat), 1:nrow(win_mat), t(win_mat[nrow(win_mat):1, ]),
      col = win_cols, zlim = c(0, 6), xaxt = "n", yaxt = "n",
      xlab = "True clock model", ylab = "", main = "treePL wins (out of 6 conditions)")
axis(1, at = 1:ncol(win_mat), labels = colnames(win_mat), las = 2)
axis(2, at = 1:nrow(win_mat), labels = rev(c("pulse score", "burst loss", "tempo composite", "MAE", "rate irregularity")), las = 2)
for (i in seq_len(nrow(win_mat))) {
  for (j in seq_len(ncol(win_mat))) {
    text(j, nrow(win_mat) - i + 1, labels = sprintf("%d", win_mat[i, j]), cex = 0.9)
  }
}
mtext("This panel highlights the specific clock-model strata where treePL still wins.", side = 1, line = 5.5, cex = 0.9)
dev.off()

# Optional compact tables for by-clock interpretation.
write.csv(clk, "data/pulse_metric_720_delta_by_clock_model.csv", row.names = FALSE)
message("Done. Pulse/rate figures written in: ", normalizePath(file.path(getwd(), "figures"), winslash = "/"))
