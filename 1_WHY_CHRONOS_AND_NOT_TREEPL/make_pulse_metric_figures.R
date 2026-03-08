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

if (!file.exists(method_csv)) stop("Missing input file: ", method_csv)
if (!file.exists(wins_csv)) stop("Missing input file: ", wins_csv)
if (!file.exists(clock_csv)) stop("Missing input file: ", clock_csv)

meth <- read.csv(method_csv, stringsAsFactors = FALSE)
wins <- read.csv(wins_csv, stringsAsFactors = FALSE)
clk <- read.csv(clock_csv, stringsAsFactors = FALSE)

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
mtext("720 benchmark: method means under the new pulse-preservation and rate-plausibility metrics", side = 3, outer = TRUE, line = -1.5, cex = 1.4)
dev.off()

# Optional compact table for by-clock deltas used in README interpretation.
write.csv(clk, "data/pulse_metric_720_delta_by_clock_model.csv", row.names = FALSE)
message("Done. Pulse/rate figures written in: ", normalizePath(file.path(getwd(), "figures"), winslash = "/"))
