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
  (m - d) / m
}

clade_key <- function(tips) paste(sort(tips), collapse = "|")

internal_node_keys <- function(tr) {
  n <- Ntip(tr)
  nodes <- (n + 1):(n + tr$Nnode)
  keys <- sapply(nodes, function(nd) clade_key(extract.clade(tr, nd)$tip.label))
  names(keys) <- as.character(nodes)
  keys
}

# Pick two focal "burst" clades from the phylogram:
# early internal nodes with large descendant sets, forcing non-nested clades.
pick_focal_nodes <- function(tr) {
  n <- Ntip(tr)
  h <- get_norm_height(tr)
  nodes <- (n + 1):(n + tr$Nnode)
  early <- nodes[h[nodes] >= quantile(h[nodes], 0.75, na.rm = TRUE)]

  tips_list <- lapply(early, function(nd) sort(extract.clade(tr, nd)$tip.label))
  sizes <- sapply(tips_list, length)
  ord <- order(sizes, decreasing = TRUE)

  selected <- integer(0)
  selected_tips <- list()

  for (k in ord) {
    cand_nd <- early[k]
    cand_tips <- tips_list[[k]]
    nested <- FALSE
    if (length(selected_tips) > 0) {
      for (st in selected_tips) {
        if (all(cand_tips %in% st) || all(st %in% cand_tips)) {
          nested <- TRUE
          break
        }
      }
    }
    if (!nested) {
      selected <- c(selected, cand_nd)
      selected_tips[[length(selected_tips) + 1]] <- cand_tips
    }
    if (length(selected) == 2) break
  }

  if (length(selected) < 2) {
    selected <- early[ord[1:2]]
  }
  selected
}

p_phy <- normalize_tree(phy)
p_mods <- lapply(mods, normalize_tree)
focal_phy_nodes <- pick_focal_nodes(p_phy)
focal_keys <- sapply(focal_phy_nodes, function(nd) clade_key(extract.clade(p_phy, nd)$tip.label))

find_nodes_by_keys <- function(tr, keys) {
  kmap <- internal_node_keys(tr)
  idx <- match(keys, unname(kmap))
  as.integer(names(kmap)[idx])
}

plot_one <- function(tr, title_txt, focal_keys) {
  plot(tr, show.tip.label = FALSE, no.margin = TRUE, main = title_txt, cex.main = 0.95)
  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  n <- Ntip(tr)
  nodes <- find_nodes_by_keys(tr, focal_keys)
  nodes <- nodes[is.finite(nodes)]
  if (length(nodes) == 0) return(invisible(NULL))

  dy <- 0.07 * diff(lp$y.lim)
  for (i in seq_along(nodes)) {
    nd <- nodes[i]
    xt <- lp$xx[nd]
    yt <- lp$yy[nd - n]
    y0 <- yt + ifelse(i %% 2 == 1, dy, -dy)
    x0 <- max(xt - 0.18, lp$x.lim[1] + 0.02)
    arrows(x0, y0, xt - 0.01, yt, length = 0.12, col = "#d7301f", lwd = 2.3)
  }
}

# 1) Tree panel with best two and worst two
order_models <- c(best2, worst2)
order_models <- unique(order_models)
if (!all(c("discrete", "clock", "correlated", "relaxed") %in% order_models)) {
  order_models <- c("discrete", "clock", "correlated", "relaxed")
}

png(file.path(out_fig, "branching_tempo_tree_panel_arrows_v2.png"), width = 2500, height = 1450, res = 180)
layout(matrix(c(1,2,3,4,5,6), nrow = 2, byrow = TRUE))
par(mar = c(1,1,4,1), oma = c(0,0,1,0))
plot_one(p_phy, "input phylogram (normalized)", focal_keys)
for (m in order_models) {
  rec <- tempo[tempo$model == m, ][1, ]
  ttl <- sprintf("%s | all %.3f | early %.3f", m, rec$tempo_mae_all, rec$tempo_mae_early_q75)
  plot_one(p_mods[[m]], ttl, focal_keys)
}
plot.new()
legend("center",
       legend = c("Red arrows: two focal burst clades used for visual comparison",
                  sprintf("Best two by composite: %s, %s", best2[1], best2[2]),
                  sprintf("Worst two by composite: %s, %s", worst2[1], worst2[2])),
       bty = "n", cex = 1.0)
mtext("Figure A: Tree-shape comparison (best vs worst)", side = 3, outer = TRUE, line = -0.5, cex = 1.5, font = 2)
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
