args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(script_dir)

infile <- "data/mae_treepl_chronos_720_c1_n2_t2_final.csv"
if (!file.exists(infile)) stop("Missing input file: ", infile)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

d <- read.csv(infile)
finite <- is.finite(d$treePL_MAE) & is.finite(d$chronos_MAE)

summ_by <- function(var) {
  lev <- sort(unique(d[[var]]))
  out <- lapply(lev, function(v) {
    s <- d[d[[var]] == v, ]
    both <- is.finite(s$treePL_MAE) & is.finite(s$chronos_MAE)
    data.frame(
      group = as.character(v),
      n_total = nrow(s),
      n_both = sum(both),
      treePL_mean = mean(s$treePL_MAE, na.rm = TRUE),
      chronos_mean = mean(s$chronos_MAE, na.rm = TRUE),
      chronos_win_rate = if (sum(both) > 0) mean(s$chronos_MAE[both] < s$treePL_MAE[both]) else NA_real_
    )
  })
  do.call(rbind, out)
}

by_clock <- summ_by("clock_model")
by_mu <- summ_by("mu")
by_het <- summ_by("heterotachy")

rec_by_true <- do.call(rbind, lapply(sort(unique(d$true_chronos_model)), function(m) {
  s <- d[d$true_chronos_model == m, ]
  data.frame(
    true_model = m,
    n = nrow(s),
    recovered = sum(s$chronos_model == s$true_chronos_model),
    rate = mean(s$chronos_model == s$true_chronos_model)
  )
}))
recovery_tbl <- as.data.frame.matrix(table(d$true_chronos_model, d$chronos_model))

write.csv(by_clock, "by_clock_summary.csv", row.names = FALSE)
write.csv(by_mu, "by_mu_summary.csv", row.names = FALSE)
write.csv(by_het, "by_heterotachy_summary.csv", row.names = FALSE)
write.csv(rec_by_true, "chronos_recovery_summary.csv", row.names = FALSE)
write.csv(recovery_tbl, "chronos_recovery_confusion_matrix.csv", row.names = TRUE)

# Figure 1: overall MAE distribution
png("figures/fig1_overall_mae_boxplot.png", width = 1400, height = 900, res = 150)
par(mar = c(6, 6, 4, 2))
boxplot(
  list(
    treePL = log10(d$treePL_MAE[is.finite(d$treePL_MAE)] + 1e-8),
    chronos = log10(d$chronos_MAE[is.finite(d$chronos_MAE)] + 1e-8)
  ),
  col = c("#9ecae1", "#a1d99b"),
  ylab = "log10(MAE)",
  main = "Overall dating error across 720 tests"
)
mtext("Lower is better; treePL non-finite runs excluded from log-scale panels", side = 1, line = 4, cex = 0.9)
dev.off()

# Figure 2: mean MAE by true clock model
order_clock <- c("strict", "autocorrelated", "independent", "discrete")
bc <- by_clock[match(order_clock, by_clock$group), ]
png("figures/fig2_mae_by_true_clock.png", width = 1500, height = 900, res = 150)
par(mar = c(8, 6, 4, 2))
mat <- rbind(bc$treePL_mean, bc$chronos_mean)
colnames(mat) <- bc$group
barplot(
  mat,
  beside = TRUE,
  col = c("#6baed6", "#74c476"),
  ylim = c(0, max(mat, na.rm = TRUE) * 1.15),
  ylab = "Mean MAE",
  main = "Mean dating error by true clock regime"
)
legend("topright", legend = c("treePL", "chronos"), fill = c("#6baed6", "#74c476"), bty = "n")
dev.off()

# Figure 3: MAE heatmaps for treePL and chronos on the same color scale
mu_vals <- sort(unique(d$mu))
het_vals <- sort(unique(d$heterotachy))

heat_treepl <- outer(mu_vals, het_vals, Vectorize(function(m, h) {
  s <- d[d$mu == m & d$heterotachy == h, ]
  mean(s$treePL_MAE, na.rm = TRUE)
}))
heat_chronos <- outer(mu_vals, het_vals, Vectorize(function(m, h) {
  s <- d[d$mu == m & d$heterotachy == h, ]
  mean(s$chronos_MAE, na.rm = TRUE)
}))

all_heat <- c(heat_treepl, heat_chronos)
zlim <- range(all_heat[is.finite(all_heat)], na.rm = TRUE)
cols <- colorRampPalette(c("#f7fbff", "#9ecae1", "#08519c"))(200)

png("figures/fig3_mae_heatmap_mu_heterotachy.png", width = 1800, height = 900, res = 150)
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(4, 4, 1.4))

par(mar = c(5, 5, 4, 1))
image(het_vals, mu_vals, t(heat_treepl), col = cols, zlim = zlim,
      xlab = "Heterotachy", ylab = "Extinction (mu)",
      main = "treePL mean MAE")
axis(1, at = het_vals, labels = het_vals)
axis(2, at = mu_vals, labels = mu_vals)

par(mar = c(5, 2, 4, 1))
image(het_vals, mu_vals, t(heat_chronos), col = cols, zlim = zlim,
      xlab = "Heterotachy", ylab = "", yaxt = "n",
      main = "chronos mean MAE")
axis(1, at = het_vals, labels = het_vals)

par(mar = c(5, 2, 4, 4))
legend_image <- matrix(seq(zlim[1], zlim[2], length.out = 200), ncol = 1)
image(1, seq(zlim[1], zlim[2], length.out = 200), t(legend_image),
      col = cols, xaxt = "n", xlab = "", ylab = "MAE", main = "Scale")
layout(1)
dev.off()

# Figure 4: interaction plot (mu x heterotachy) with shared Y-axis across methods
comb <- expand.grid(mu = mu_vals, heterotachy = het_vals)
comb$treePL <- mapply(function(m, h) {
  s <- d[d$mu == m & d$heterotachy == h, ]
  mean(s$treePL_MAE, na.rm = TRUE)
}, comb$mu, comb$heterotachy)
comb$chronos <- mapply(function(m, h) {
  s <- d[d$mu == m & d$heterotachy == h, ]
  mean(s$chronos_MAE, na.rm = TRUE)
}, comb$mu, comb$heterotachy)

y_all <- c(comb$treePL, comb$chronos)
y_all <- y_all[is.finite(y_all)]
ylim_shared <- range(y_all, na.rm = TRUE)
if (diff(ylim_shared) == 0) ylim_shared <- ylim_shared + c(-0.01, 0.01)

plot_method <- function(method_col, method_name) {
  sub <- comb[, c("mu", "heterotachy", method_col)]
  names(sub)[3] <- "MAE"
  plot(NA,
       xlim = range(mu_vals),
       ylim = ylim_shared,
       xlab = expression(mu),
       ylab = "Mean MAE",
       main = paste0(method_name, " interaction: heterotachy x extinction"))
  for (i in seq_along(het_vals)) {
    h <- het_vals[i]
    tmp <- sub[sub$heterotachy == h, ]
    tmp <- tmp[order(tmp$mu), ]
    lines(tmp$mu, tmp$MAE,
          type = "b",
          lwd = 2,
          pch = if (i == 1) 16 else 17,
          lty = if (i == 1) 1 else 2,
          col = if (i == 1) "#08519c" else "#cb181d")
  }
  legend("topleft",
         legend = paste("H =", het_vals),
         lty = c(1, 2),
         pch = c(16, 17),
         col = c("#08519c", "#cb181d"),
         lwd = 2,
         bty = "n")
}

png("figures/fig4_interaction_heterotachy_extinction.png", width = 1700, height = 850, res = 150)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
plot_method("treePL", "treePL")
plot_method("chronos", "chronos")
dev.off()

message("Done. Updated tables and figures in: ", normalizePath(getwd(), winslash = "/"))
