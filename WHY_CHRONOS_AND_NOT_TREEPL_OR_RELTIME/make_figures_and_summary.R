args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  script_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(script_dir)

preferred_infile <- "data/mae_treepl_chronos_reltime_720_c1_n2_t2_final.csv"
fallback_infile <- "data/mae_treepl_chronos_720_c1_n2_t2_final.csv"
if (file.exists(preferred_infile)) {
  infile <- preferred_infile
} else if (file.exists(fallback_infile)) {
  infile <- fallback_infile
} else {
  stop("Missing input file. Tried: ", preferred_infile, " and ", fallback_infile)
}

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

d <- read.csv(infile)
has_rt <- "reltime_MAE" %in% names(d)
if (!has_rt) {
  d$reltime_MAE <- NA_real_
  warning("Input file has no reltime_MAE column. RelTime summaries will be NA until the refreshed three-method CSV is provided.")
}
method_cols <- c(treePL = "treePL_MAE", chronos = "chronos_MAE", RelTime = "reltime_MAE")
method_names <- names(method_cols)
method_labels <- c(treePL = "treePL", chronos = "chronos", RelTime = "RelTime")[method_names]
method_colors <- c(treePL = "#6baed6", chronos = "#74c476", RelTime = "#fd8d3c")[method_names]
plot_methods <- method_names[vapply(method_cols, function(col) any(is.finite(d[[col]])), logical(1))]

mean_if_finite <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

summ_by <- function(var) {
  lev <- sort(unique(d[[var]]))
  out <- lapply(lev, function(v) {
    s <- d[d[[var]] == v, ]
    finite_list <- lapply(method_cols, function(col) is.finite(s[[col]]))
    all_methods <- Reduce(`&`, finite_list)
    row <- data.frame(
      group = as.character(v),
      n_total = nrow(s),
      n_all_methods = sum(all_methods)
    )
    for (i in seq_along(method_cols)) {
      method <- names(method_cols)[i]
      col <- method_cols[[i]]
      row[[paste0(method, "_mean")]] <- mean_if_finite(s[[col]])
    }
    if (sum(all_methods) > 0) {
      score_mat <- as.matrix(s[all_methods, unname(method_cols), drop = FALSE])
      row_mins <- apply(score_mat, 1, min, na.rm = TRUE)
      for (i in seq_along(method_cols)) {
        method <- names(method_cols)[i]
        row[[paste0(method, "_win_rate")]] <- mean(score_mat[, i] == row_mins)
      }
    } else {
      for (method in names(method_cols)) {
        row[[paste0(method, "_win_rate")]] <- NA_real_
      }
    }
    row
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

# Figure 1: overall MAE distribution across available methods
png("figures/fig1_overall_mae_boxplot.png", width = 1400, height = 900, res = 150)
par(mar = c(6, 6, 4, 2))
box_data <- lapply(method_cols, function(col) log10(d[[col]][is.finite(d[[col]])] + 1e-8))
names(box_data) <- method_labels
boxplot(
  box_data[plot_methods],
  col = unname(method_colors[plot_methods]),
  ylab = "log10(MAE)",
  main = "Overall dating error across the exact-root benchmark (chronos, treePL, RelTime)"
)
mtext("Lower is better; non-finite runs are excluded from log-scale panels", side = 1, line = 4, cex = 0.9)
dev.off()

# Figure 2: mean MAE by true clock model
order_clock <- c("strict", "autocorrelated", "independent", "discrete")
bc <- by_clock[match(order_clock, by_clock$group), ]
png("figures/fig2_mae_by_true_clock.png", width = 1500, height = 900, res = 150)
par(mar = c(8, 6, 4, 2))
mat <- do.call(rbind, lapply(plot_methods, function(method) bc[[paste0(method, "_mean")]]))
colnames(mat) <- bc$group
rownames(mat) <- method_labels[plot_methods]
barplot(
  mat,
  beside = TRUE,
  col = unname(method_colors[plot_methods]),
  ylim = c(0, max(mat, na.rm = TRUE) * 1.15),
  ylab = "Mean MAE",
  main = "Mean dating error by true clock regime"
)
legend("topright", legend = method_labels[plot_methods], fill = unname(method_colors[plot_methods]), bty = "n")
dev.off()

# Figure 3: MAE heatmaps for the exact-root comparison on the same color scale
mu_vals <- sort(unique(d$mu))
het_vals <- sort(unique(d$heterotachy))

heat_maps <- lapply(method_cols[plot_methods], function(col) {
  outer(mu_vals, het_vals, Vectorize(function(m, h) {
    s <- d[d$mu == m & d$heterotachy == h, ]
    mean_if_finite(s[[col]])
  }))
})
all_heat <- unlist(heat_maps, use.names = FALSE)
zlim <- range(all_heat[is.finite(all_heat)], na.rm = TRUE)
cols <- colorRampPalette(c("#f7fbff", "#9ecae1", "#08519c"))(200)

png("figures/fig3_mae_heatmap_mu_heterotachy.png", width = 2200, height = 900, res = 150)
layout(matrix(seq_len(length(plot_methods) + 1), nrow = 1), widths = c(rep(4, length(plot_methods)), 1.4))

for (i in seq_along(plot_methods)) {
  method <- plot_methods[i]
  heat <- heat_maps[[i]]
  if (i == 1) {
    par(mar = c(5, 5, 4, 1))
  } else {
    par(mar = c(5, 2, 4, 1))
  }
  image(het_vals, mu_vals, t(heat), col = cols, zlim = zlim,
        xlab = "Heterotachy",
        ylab = if (i == 1) "Extinction (mu)" else "",
        yaxt = if (i == 1) "s" else "n",
        main = paste0(method_labels[[method]], " mean MAE"))
  axis(1, at = het_vals, labels = het_vals)
  if (i == 1) axis(2, at = mu_vals, labels = mu_vals)
}

par(mar = c(5, 2, 4, 4))
legend_image <- matrix(seq(zlim[1], zlim[2], length.out = 200), ncol = 1)
image(1, seq(zlim[1], zlim[2], length.out = 200), t(legend_image),
      col = cols, xaxt = "n", xlab = "", ylab = "MAE", main = "Scale")
layout(1)
dev.off()

# Figure 4: interaction plot (mu x heterotachy) with shared Y-axis across methods
comb <- expand.grid(mu = mu_vals, heterotachy = het_vals)
for (method in method_names) {
  col <- method_cols[[method]]
  comb[[method]] <- mapply(function(m, h) {
    s <- d[d$mu == m & d$heterotachy == h, ]
    mean_if_finite(s[[col]])
  }, comb$mu, comb$heterotachy)
}

y_all <- unlist(comb[plot_methods], use.names = FALSE)
y_all <- y_all[is.finite(y_all)]
ylim_shared <- range(y_all, na.rm = TRUE)
if (diff(ylim_shared) == 0) ylim_shared <- ylim_shared + c(-0.01, 0.01)
y_ticks <- pretty(ylim_shared, n = 6)

plot_method <- function(method_col, method_name) {
  sub <- comb[, c("mu", "heterotachy", method_col)]
  names(sub)[3] <- "MAE"
  plot(NA,
       xlim = range(mu_vals),
       ylim = ylim_shared,
       yaxt = "n",
       xlab = expression(mu),
       ylab = "Mean MAE",
       main = paste0(method_name, " interaction"))
  axis(2, at = y_ticks, labels = y_ticks)
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

png("figures/fig4_interaction_heterotachy_extinction.png", width = 2300, height = 850, res = 150)
par(mfrow = c(1, length(plot_methods)), mar = c(5, 5, 4, 2))
for (method in plot_methods) {
  plot_method(method, method_labels[[method]])
}
dev.off()

# Cache-busting duplicate with same content for GitHub web view refreshes
png("figures/fig4_interaction_heterotachy_extinction_same_scale.png", width = 2300, height = 850, res = 150)
par(mfrow = c(1, length(plot_methods)), mar = c(5, 5, 4, 2))
for (method in plot_methods) {
  plot_method(method, method_labels[[method]])
}
dev.off()

message("Done. Updated tables and figures in: ", normalizePath(getwd(), winslash = "/"))
