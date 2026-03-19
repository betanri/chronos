args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(script_dir)

data_dir <- file.path(script_dir, "data")
fig_dir <- file.path(script_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

read_if_populated <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- read.csv(path, stringsAsFactors = FALSE)
  if (!nrow(x)) return(NULL)
  x
}

method_df <- read_if_populated(file.path(data_dir, "benchmark_suite_method_summary.csv"))
recovery_df <- read_if_populated(file.path(data_dir, "benchmark_suite_recovery_summary.csv"))
rank_df <- read_if_populated(file.path(data_dir, "benchmark_suite_rank_summary.csv"))
pcr_df <- read_if_populated(file.path(data_dir, "pcr_summary_table.csv"))
tree_panel_df <- read_if_populated(file.path(data_dir, "representative_tree_panel_manifest.csv"))

if (!is.null(rank_df)) {
  ranks <- rank_df
  meth_levels <- unique(ranks$method_family)
  bench_levels <- unique(ranks$benchmark_label)
  mat <- matrix(NA_real_, nrow = length(meth_levels), ncol = length(bench_levels),
                dimnames = list(meth_levels, bench_levels))
  for (i in seq_len(nrow(ranks))) {
    mat[ranks$method_family[i], ranks$benchmark_label[i]] <- ranks$mean_rank[i]
  }

  png(file.path(fig_dir, "fig1_benchmark_suite_rank_summary.png"), width = 1700, height = 900, res = 150)
  par(mar = c(8, 6, 4, 2))
  barplot(mat,
          beside = TRUE,
          col = c("#2b8cbe", "#31a354", "#e6550d")[seq_len(nrow(mat))],
          ylab = "Mean rank",
          main = "Benchmark-suite rank summary across A-E")
  legend("topright", legend = rownames(mat),
         fill = c("#2b8cbe", "#31a354", "#e6550d")[seq_len(nrow(mat))], bty = "n")
  dev.off()
}

if (!is.null(method_df)) {
  methods <- method_df
  bench_levels <- unique(methods$benchmark_label)
  meth_levels <- c("chronos-clock", "chronos-discrete", "chronos-correlated", "chronos-relaxed",
                   "chronos overall", "RelTime", "treePL")
  methods$method <- factor(methods$method, levels = meth_levels)

  png(file.path(fig_dir, "fig2_benchmark_suite_model_specific_mae.png"), width = 2200, height = 1200, res = 150)
  par(mfrow = c(ceiling(length(bench_levels) / 2), 2), mar = c(10, 5, 4, 1))
  for (b in bench_levels) {
    sub <- methods[methods$benchmark_label == b, ]
    sub <- sub[order(sub$method), ]
    vals <- sub$mae
    names(vals) <- as.character(sub$method)
    vals <- vals[!is.na(vals)]
    cols <- ifelse(grepl("^chronos-", names(vals)), "#4292c6",
                   ifelse(names(vals) == "chronos overall", "#08519c",
                          ifelse(names(vals) == "RelTime", "#31a354", "#e6550d")))
    barplot(vals,
            las = 2,
            col = cols,
            ylab = "MAE",
            main = sprintf("%s (%.2f reps)", b, unique(sub$completed_reps)))
    mtext("Lower is better", side = 1, line = 8, cex = 0.8)
  }
  dev.off()
}

if (!is.null(recovery_df)) {
  rec <- recovery_df
  bench_levels <- unique(rec$benchmark_label)
  model_levels <- c("clock", "discrete", "correlated", "relaxed")

  png(file.path(fig_dir, "fig3_benchmark_suite_model_recovery.png"), width = 1800, height = 1000, res = 150)
  par(mfrow = c(1, length(bench_levels)), mar = c(8, 5, 4, 1))
  for (b in bench_levels) {
    sub <- rec[rec$benchmark_label == b, ]
    vals <- setNames(sub$recovery_rate, sub$true_model)
    vals <- vals[model_levels]
    barplot(vals,
            ylim = c(0, 1),
            las = 2,
            col = c("#08519c", "#2171b5", "#6baed6", "#bdd7e7"),
            ylab = "Recovery rate",
            main = sprintf("%s (%.2f reps)", b, unique(sub$completed_reps)))
    abline(h = c(0.25, 0.5, 0.75, 1.0), col = "grey90", lty = 3)
  }
  dev.off()
}

if (!is.null(pcr_df)) {
  pcr <- pcr_df
  datasets <- unique(pcr$dataset)
  fams <- unique(pcr$metric_family)
  score <- matrix(0, nrow = length(fams), ncol = length(datasets), dimnames = list(fams, datasets))
  for (i in seq_len(nrow(pcr))) {
    score[pcr$metric_family[i], pcr$dataset[i]] <- match(pcr$best_method[i], c("chronos", "RelTime", "treePL"))
  }
  cols <- c("chronos" = "#08519c", "RelTime" = "#31a354", "treePL" = "#e6550d")
  png(file.path(fig_dir, "fig5_pcr_summary.png"), width = 1500, height = 900, res = 150)
  par(mar = c(8, 6, 4, 2))
  barplot(score,
          beside = TRUE,
          col = c("#6baed6", "#74c476", "#fd8d3c")[seq_len(nrow(score))],
          ylab = "Best-method code",
          main = "Compact PCR summary")
  legend("topright", legend = names(cols), fill = unname(cols), bty = "n")
  dev.off()
}

if (!is.null(tree_panel_df)) {
  message("Representative tree panel manifest detected. Fig 4 is intended to be assembled as one row per benchmark (A-E) with fixed columns: reference, chronos, RelTime, treePL.")
}

message("Suite scaffold refresh complete: ", normalizePath(fig_dir, winslash = "/"))
