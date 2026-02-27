suppressPackageStartupMessages({ library(ape) })

base_dir <- "2_CHRONOS_CUSTOM_DATING_TREE_PIPELINE"
out_fig <- file.path(base_dir, "figures")
out_csv <- file.path(base_dir, "EXAMPLE_FILES", "OUTPUT_DEMO")

tempo <- read.csv(file.path(out_csv, "summary_terap_empirical_model_fits.csv"), stringsAsFactors = FALSE)

tempo <- tempo[order(tempo$tempo_composite, tempo$tempo_mae_all), ]
best2 <- tempo$model[1:2]
worst2 <- tempo$model[(nrow(tempo)-1):nrow(tempo)]

read_tree <- function(model) {
  p <- file.path(out_csv, sprintf("TERAP_ML_MAIN_chronos_dated_model%s.tre", model))
  read.tree(p)
}

phy <- read.tree(file.path(out_csv, "TERAP_ML_MAIN_phylogram_used.tree"))
mods <- list(
  discrete = read_tree("discrete"),
  clock = read_tree("clock"),
  correlated = read_tree("correlated"),
  relaxed = read_tree("relaxed")
)

normalize_tree <- function(tr) {
  d <- node.depth.edgelength(tr)
  m <- max(d[1:Ntip(tr)])
  tr$edge.length <- tr$edge.length / m
  tr
}

get_norm_height <- function(tr) {
  d <- node.depth.edgelength(tr)
  m <- max(d[1:Ntip(tr)])
  h <- (m - d) / m
  h
}

plot_one <- function(tr, title_txt, mark_early = TRUE, cex_title = 1.1) {
  plot(tr, show.tip.label = FALSE, no.margin = TRUE, main = title_txt, cex.main = cex_title)
  if (mark_early) {
    lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    h <- get_norm_height(tr)
    n <- Ntip(tr)
    internal_nodes <- (n + 1):(n + tr$Nnode)
    early <- internal_nodes[h[internal_nodes] >= quantile(h[internal_nodes], 0.75, na.rm = TRUE)]
    idx <- early - n
    points(lp$xx[early], lp$yy[idx], pch = 16, cex = 0.55, col = "#d7301f")
  }
}

# 1) Tree panel with best two and worst two
p_phy <- normalize_tree(phy)
p_mods <- lapply(mods, normalize_tree)

order_models <- c(best2, worst2)
order_models <- unique(order_models)
if (!all(c("discrete", "clock", "correlated", "relaxed") %in% order_models)) {
  order_models <- c("discrete", "clock", "correlated", "relaxed")
}

png(file.path(out_fig, "branching_tempo_tree_panel.png"), width = 2200, height = 1300, res = 180)
layout(matrix(c(1,2,3,4,5,6), nrow = 2, byrow = TRUE))
par(mar = c(1,1,3,1))
plot_one(p_phy, "Input phylogram (normalized)", mark_early = TRUE)
for (m in order_models) {
  rec <- tempo[tempo$model == m, ][1, ]
  ttl <- sprintf("%s | tempo_all=%.3f | early=%.3f", m, rec$tempo_mae_all, rec$tempo_mae_early_q75)
  plot_one(p_mods[[m]], ttl, mark_early = TRUE)
}
plot.new()
legend("center",
       legend = c("Red points: early-burst internal nodes (top quartile normalized height)",
                  sprintf("Best two by composite: %s, %s", best2[1], best2[2]),
                  sprintf("Worst two by composite: %s, %s", worst2[1], worst2[2])),
       bty = "n", cex = 1.0)
dev.off()

# 2) Node-height ECDF vs phylogram
get_internal_h <- function(tr) {
  h <- get_norm_height(tr)
  n <- Ntip(tr)
  h[(n + 1):(n + tr$Nnode)]
}

h_phy <- get_internal_h(phy)
hs <- lapply(mods, get_internal_h)

png(file.path(out_fig, "branching_tempo_node_height_ecdf.png"), width = 1400, height = 950, res = 180)
par(mar = c(5,5,3,1))
plot(ecdf(h_phy), do.points = FALSE, verticals = FALSE, lwd = 3, col = "black",
     xlab = "Normalized node height (0 tip-side, 1 root-side)", ylab = "ECDF",
     main = "Branching-tempo shape vs input phylogram")
cols <- c(discrete = "#1b9e77", clock = "#2c7fb8", correlated = "#d95f0e", relaxed = "#7570b3")
for (m in names(hs)) {
  lines(ecdf(hs[[m]]), do.points = FALSE, verticals = FALSE, lwd = 2.5, col = cols[m], lty = ifelse(m %in% best2, 1, 2))
}
legend("bottomright",
       legend = c("phylogram", "discrete", "clock", "correlated", "relaxed"),
       col = c("black", cols), lwd = c(3,2.5,2.5,2.5,2.5), lty = c(1,1,1,2,2), bty = "n")
dev.off()

# 3) Tempo metrics bar chart
show <- tempo[match(c("discrete","clock","correlated","relaxed"), tempo$model), ]

png(file.path(out_fig, "branching_tempo_metric_bars.png"), width = 1400, height = 950, res = 180)
par(mar = c(8,5,3,1))
mat <- rbind(show$tempo_mae_all, show$tempo_mae_early_q75)
colnames(mat) <- show$model
barplot(mat, beside = TRUE, col = c("#6baed6", "#fb6a4a"), las = 2,
        ylab = "Tempo MAE (lower is better)", main = "Branching-tempo metric by model")
legend("topleft", legend = c("overall MAE", "early_q75 MAE"), fill = c("#6baed6", "#fb6a4a"), bty = "n")
dev.off()

write.csv(tempo, file.path(out_csv, "summary_terap_empirical_model_fits_ranked.csv"), row.names = FALSE)
