#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(phytools)
})

# ============================================================
# Chronos empirical dating pipeline for Terapontoid trees
# - Option A: build internal calibrations via congruification
# - Option B: use manual calibration CSV
# - Uses robust chronos selector (c1_n2_t2 defaults)
# ============================================================

# -------------------------
# USER SETTINGS
# -------------------------
USE_CONGRUIF <- FALSE

# Set working dir.
setwd("/Users/ricardobetancur/Desktop/Proxy_Misplaced/chronosPL/Terapontoid_Trees")

# Target tree input can be:
# 1) named-Newick file with many trees: "NAME = (....);"
# 2) standard single-tree Newick/NEXUS file
TARGET_TREE_FILE <- "TARGET_PHYLOGRAM_Terapontoidei_ML_MAIN.tre"

# If TARGET_TREE_FILE has multiple trees, choose by name OR index.
TARGET_TREE_NAME <- ""   # e.g., "Terapontoid_ASTRAL_S1" (leave "" to use TARGET_TREE_INDEX)
TARGET_TREE_INDEX <- 1L

# Reference time tree for congruification (stock tree).
REFERENCE_TIME_TREE <- "REFERENCE_TIMETREE_Eupercaria_Congru_yesplectoBEAST_ASTRAL_sub1.tree"

# If not using congruification, provide manual internal calibrations CSV.
# Required columns: taxonA,taxonB,age_min,age_max
MANUAL_CAL_CSV <- "MANUAL_CALIBRATIONS.csv"

# Optional: fix root age to this value. Set NA to skip root constraint.
ROOT_AGE <- NA_real_

# Robust chronos selector settings (c1_n2_t2 defaults).
LAMBDA_GRID <- c(0.1, 1, 10)
CHRONOS_MODELS <- c("clock", "correlated", "relaxed", "discrete")
K_FIT_GRID <- c(2L, 3L, 5L, 10L)
N_RETRIES <- 2L
PLOG_CLOCK_SWITCH_THRESH <- 1
PLOG_NONCLOCK_SWITCH_THRESH <- 2
PLOG_TIE_EPS <- 2
# Empirical model-sensitivity protocol: run only strict and stricter thresholds.
CLOCK_SWITCH_SENSITIVITY <- c(1, 2)
TEMPO_EQ_ABS_TOL <- 0.001
TEMPO_EQ_REL_TOL <- 0.01

# Optional subset mode for large trees (empirical acceleration).
# Keeps calibration taxa, root-to-tip extremes, and diversified spread.
USE_SUBSET <- FALSE
SUBSET_N <- 400L
SUBSET_EXTREME_FRAC <- 0.05
SUBSET_SEED <- 1L
# If TRUE, subset is used only for tuning model/lambda; final selected settings are
# then applied to the full tree for final dated outputs.
SUBSET_TUNE_ON_SUBSET_ONLY <- TRUE

# Output
OUT_BASE_DIR <- "CHRONOS_OUT"
OUT_PREFIX <- "terap_empirical"
CLEAN_PREVIOUS_PREFIX_OUTPUTS <- TRUE

# Optional diagnostic stress mode (empirical-safe).
# Enable with env: DIAGNOSTIC_MODE=TRUE
# Optional env overrides:
#   DIAG_LAMBDA_GRID=0.01,0.1,1,10,100
#   DIAG_N_RETRIES=5
#   DIAG_OUT_SUFFIX=_diagflex
DIAGNOSTIC_MODE <- tolower(Sys.getenv("DIAGNOSTIC_MODE", unset = "false")) %in% c("1", "true", "yes")
if (DIAGNOSTIC_MODE) {
  diag_lam <- strsplit(Sys.getenv("DIAG_LAMBDA_GRID", unset = "0.01,0.1,1,10,100"), ",")[[1]]
  diag_lam <- as.numeric(trimws(diag_lam))
  diag_lam <- diag_lam[is.finite(diag_lam) & diag_lam > 0]
  if (length(diag_lam)) LAMBDA_GRID <- unique(diag_lam)
  diag_retries <- as.integer(Sys.getenv("DIAG_N_RETRIES", unset = "5"))
  if (is.finite(diag_retries) && diag_retries >= 1L) N_RETRIES <- diag_retries
  OUT_PREFIX <- paste0(OUT_PREFIX, Sys.getenv("DIAG_OUT_SUFFIX", unset = "_diagflex"))
}

# -------------------------
# Helpers
# -------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

is_nexus_file <- function(path, max_lines = 50) {
  con <- file(path, "r")
  on.exit(close(con), add = TRUE)
  for (i in seq_len(max_lines)) {
    ln <- readLines(con, n = 1, warn = FALSE)
    if (!length(ln)) break
    s <- trimws(ln)
    if (!nzchar(s)) next
    su <- toupper(s)
    return(grepl("^#?NEXUS\\b", su) || grepl("\\bBEGIN\\s+TREES\\b", su))
  }
  FALSE
}

read_named_newick <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- x[grepl("=", x, fixed = TRUE)]
  if (!length(x)) stop("No named-Newick lines found in: ", path)
  nm <- trimws(sub("\\s*=.*$", "", x))
  nw <- trimws(sub("^.*=\\s*", "", x))
  tr <- lapply(nw, function(s) read.tree(text = s))
  names(tr) <- nm
  class(tr) <- "multiPhylo"
  tr
}

read_any_tree <- function(path) {
  if (!file.exists(path)) stop("Tree file not found: ", path)
  tr <- NULL
  if (is_nexus_file(path)) {
    tr <- try(read.nexus(path), silent = TRUE)
    if (inherits(tr, "try-error")) tr <- try(read.tree(path), silent = TRUE)
  } else {
    tr <- try(read.tree(path), silent = TRUE)
    if (inherits(tr, "try-error")) tr <- try(read.nexus(path), silent = TRUE)
  }
  if (inherits(tr, "try-error")) stop("Could not parse tree file: ", path)
  tr
}

load_target_tree <- function(path, tree_name = "", tree_index = 1L) {
  x <- try(read_named_newick(path), silent = TRUE)
  if (!inherits(x, "try-error")) {
    if (nzchar(tree_name)) {
      if (!(tree_name %in% names(x))) stop("TARGET_TREE_NAME not found: ", tree_name)
      return(list(tree = x[[tree_name]], tree_id = tree_name))
    }
    if (tree_index < 1L || tree_index > length(x)) stop("TARGET_TREE_INDEX out of range.")
    return(list(tree = x[[tree_index]], tree_id = names(x)[tree_index]))
  }

  y <- read_any_tree(path)
  if (inherits(y, "multiPhylo")) {
    if (nzchar(tree_name) && !is.null(names(y)) && tree_name %in% names(y)) {
      return(list(tree = y[[tree_name]], tree_id = tree_name))
    }
    if (tree_index < 1L || tree_index > length(y)) stop("TARGET_TREE_INDEX out of range for multiPhylo.")
    id <- if (!is.null(names(y)) && nzchar(names(y)[tree_index])) names(y)[tree_index] else paste0("tree_", tree_index)
    return(list(tree = y[[tree_index]], tree_id = id))
  }
  list(tree = y, tree_id = tools::file_path_sans_ext(basename(path)))
}

normalize_calib_df <- function(df) {
  names(df) <- tolower(names(df))
  if (!("taxona" %in% names(df)) && "taxon_a" %in% names(df)) df$taxona <- df$taxon_a
  if (!("taxonb" %in% names(df)) && "taxon_b" %in% names(df)) df$taxonb <- df$taxon_b
  if (!("age_min" %in% names(df)) && "minage" %in% names(df)) df$age_min <- df$minage
  if (!("age_max" %in% names(df)) && "maxage" %in% names(df)) df$age_max <- df$maxage
  if (!("age_min" %in% names(df)) && "age" %in% names(df)) df$age_min <- df$age
  if (!("age_max" %in% names(df)) && "age" %in% names(df)) df$age_max <- df$age
  need <- c("taxona", "taxonb", "age_min", "age_max")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Calibration CSV missing: ", paste(miss, collapse = ", "))
  out <- data.frame(
    taxonA = as.character(df$taxona),
    taxonB = as.character(df$taxonb),
    age_min = as.numeric(df$age_min),
    age_max = as.numeric(df$age_max),
    stringsAsFactors = FALSE
  )
  out[is.finite(out$age_min) & is.finite(out$age_max), , drop = FALSE]
}

build_calib_pairs_congruif <- function(target_tree, ref_tree) {
  if (!is.ultrametric(ref_tree)) {
    # Congruification expects a time tree as stock; extend tips only.
    ref_tree <- force.ultrametric(ref_tree, method = "extend")
  }
  cg <- congruify.phylo(ref_tree, target_tree)
  cc <- cg$calibrations
  if (is.null(cc) || nrow(cc) == 0) stop("Congruification returned zero calibrations.")
  data.frame(
    taxonA = as.character(cc$taxonA),
    taxonB = as.character(cc$taxonB),
    age_min = as.numeric(cc$MinAge),
    age_max = as.numeric(cc$MaxAge),
    stringsAsFactors = FALSE
  )
}

build_chronos_calib <- function(phy, pair_df, root_age = NA_real_) {
  nodes <- integer(0)
  age_min <- numeric(0)
  age_max <- numeric(0)

  if (!is.null(pair_df) && nrow(pair_df)) {
    for (i in seq_len(nrow(pair_df))) {
      a <- pair_df$taxonA[i]
      b <- pair_df$taxonB[i]
      if (!(a %in% phy$tip.label) || !(b %in% phy$tip.label)) next
      nd <- getMRCA(phy, c(a, b))
      if (is.null(nd) || !is.finite(nd)) next
      nodes <- c(nodes, nd)
      age_min <- c(age_min, pair_df$age_min[i])
      age_max <- c(age_max, pair_df$age_max[i])
    }
  }

  if (is.finite(root_age)) {
    r <- Ntip(phy) + 1L
    nodes <- c(nodes, r)
    age_min <- c(age_min, root_age)
    age_max <- c(age_max, root_age)
  }
  if (!length(nodes)) stop("No usable calibration nodes on this target tree.")

  # Merge duplicate-node constraints by interval intersection.
  u <- sort(unique(nodes))
  keep <- logical(length(u))
  min_u <- numeric(length(u))
  max_u <- numeric(length(u))
  for (j in seq_along(u)) {
    idx <- which(nodes == u[j])
    mn <- max(age_min[idx], na.rm = TRUE)
    mx <- min(age_max[idx], na.rm = TRUE)
    if (is.finite(mn) && is.finite(mx) && mn <= mx) {
      keep[j] <- TRUE
      min_u[j] <- mn
      max_u[j] <- mx
    }
  }
  if (!any(keep)) stop("All calibration intervals became inconsistent after merge.")
  makeChronosCalib(phy, node = u[keep], age.min = min_u[keep], age.max = max_u[keep])
}

fit_score_chronos <- function(tr) {
  ph <- attr(tr, "PHIIC")
  if (is.list(ph) && is.finite(ph$PHIIC)) return(ph$PHIIC)
  pl <- attr(tr, "ploglik")
  if (is.finite(pl)) return(-pl)
  Inf
}

node_sig_internal <- function(tr) {
  n <- Ntip(tr)
  pp <- prop.part(tr)
  s <- character(tr$Nnode)
  for (i in seq_len(tr$Nnode)) {
    tips <- sort(tr$tip.label[pp[[i]]])
    s[i] <- paste(tips, collapse = "|")
  }
  names(s) <- as.character((n + 1):(n + tr$Nnode))
  s
}

node_heights_norm <- function(tr) {
  n <- Ntip(tr)
  d <- node.depth.edgelength(tr)
  maxd <- max(d[seq_len(n)])
  h <- 1 - (d / maxd) # root ~1, tips ~0
  names(h) <- as.character(seq_along(h))
  h
}

branching_tempo_error <- function(phy_ref, tr_est) {
  ref_sig <- node_sig_internal(phy_ref)
  est_sig <- node_sig_internal(tr_est)
  ref_h <- node_heights_norm(phy_ref)
  est_h <- node_heights_norm(tr_est)
  ref_df <- data.frame(sig = unname(ref_sig), h_ref = as.numeric(ref_h[names(ref_sig)]), stringsAsFactors = FALSE)
  est_df <- data.frame(sig = unname(est_sig), h_est = as.numeric(est_h[names(est_sig)]), stringsAsFactors = FALSE)
  j <- merge(ref_df, est_df, by = "sig", all = FALSE)
  e <- abs(j$h_est - j$h_ref)
  q75 <- quantile(j$h_ref, 0.75, na.rm = TRUE)
  early <- j$h_ref >= q75
  list(
    n_nodes = nrow(j),
    mae_all = mean(e),
    mae_early_q75 = mean(e[early]),
    median_early_q75 = median(e[early])
  )
}

metric_node_sig_internal <- function(tr) {
  n <- ape::Ntip(tr)
  pp <- ape::prop.part(tr)
  s <- character(tr$Nnode)
  for (i in seq_len(tr$Nnode)) {
    tips <- sort(tr$tip.label[pp[[i]]])
    s[i] <- paste(tips, collapse = "|")
  }
  names(s) <- as.character((n + 1):(n + tr$Nnode))
  s
}

metric_node_heights_norm <- function(tr) {
  n <- ape::Ntip(tr)
  d <- ape::node.depth.edgelength(tr)
  maxd <- max(d[seq_len(n)])
  h <- 1 - (d / maxd)
  names(h) <- as.character(seq_along(h))
  h
}

metric_node_ages <- function(tr) {
  n <- ape::Ntip(tr)
  d <- ape::node.depth.edgelength(tr)
  maxd <- max(d[seq_len(n)])
  age <- maxd - d
  names(age) <- as.character(seq_along(age))
  age
}

metric_wasserstein_1d <- function(x, y) {
  x <- sort(as.numeric(x))
  y <- sort(as.numeric(y))
  if (!length(x) || !length(y)) return(NA_real_)
  m <- max(length(x), length(y), 32L)
  p <- seq(0, 1, length.out = m)
  qx <- as.numeric(stats::quantile(x, p, type = 8, names = FALSE))
  qy <- as.numeric(stats::quantile(y, p, type = 8, names = FALSE))
  mean(abs(qx - qy))
}

metric_event_times_relative <- function(tr) {
  n <- ape::Ntip(tr)
  if (tr$Nnode < 2L) return(numeric(0))
  d <- ape::node.depth.edgelength(tr)
  max_tip <- max(d[seq_len(n)])
  if (!is.finite(max_tip) || max_tip <= 0) return(numeric(0))
  h <- 1 - (d / max_tip)
  hi <- h[(n + 1L):(n + tr$Nnode)]
  crown <- max(hi, na.rm = TRUE)
  if (!is.finite(crown) || crown <= 0) return(numeric(0))
  ev <- hi[hi < (crown - 1e-12)] / crown
  ev <- ev[is.finite(ev)]
  sort(pmax(0, pmin(1, ev)))
}

metric_burstiness_from_events <- function(ev) {
  if (!length(ev)) return(NA_real_)
  waits <- diff(c(0, ev, 1))
  if (length(waits) < 2L) return(NA_real_)
  m <- mean(waits)
  if (!is.finite(m) || m <= 0) return(NA_real_)
  stats::sd(waits) / m
}

metric_build_pulse_panel <- function(ref, min_tips = 8L, min_events = 4L) {
  ref_nodes <- (ape::Ntip(ref) + 1L):(ape::Ntip(ref) + ref$Nnode)
  panel <- list()
  for (nd in ref_nodes) {
    sub <- try(ape::extract.clade(ref, nd), silent = TRUE)
    if (inherits(sub, "try-error") || !inherits(sub, "phylo")) next
    nt <- ape::Ntip(sub)
    if (nt < min_tips || nt >= ape::Ntip(ref)) next
    ev <- metric_event_times_relative(sub)
    if (length(ev) < min_events) next
    b <- metric_burstiness_from_events(ev)
    if (!is.finite(b)) next
    panel[[length(panel) + 1L]] <- list(
      ref_node = nd,
      tips = sort(sub$tip.label),
      n_tips = nt,
      n_events = length(ev),
      ev_ref = ev,
      burst_ref = b,
      centroid_ref = mean(ev)
    )
  }
  panel
}

metric_score_pulse_panel <- function(ref, tr, panel,
                                     w_emd = 0.35,
                                     w_burst_loss = 0.55,
                                     w_centroid = 0.10,
                                     w_global = 0.20,
                                     coverage_penalty = 0.20) {
  global_ref <- metric_event_times_relative(ref)
  global_burst_ref <- metric_burstiness_from_events(global_ref)

  w_sum <- 0
  emd_sum <- 0
  burst_loss_sum <- 0
  centroid_sum <- 0
  matched <- 0L

  for (k in seq_along(panel)) {
    cl <- panel[[k]]
    node_c <- ape::getMRCA(tr, cl$tips)
    if (!is.finite(node_c)) next
    sub_c <- try(ape::extract.clade(tr, node_c), silent = TRUE)
    if (inherits(sub_c, "try-error") || !inherits(sub_c, "phylo")) next
    if (!setequal(sub_c$tip.label, cl$tips)) next
    ev_c <- metric_event_times_relative(sub_c)
    if (length(ev_c) < 2L) next
    emd <- metric_wasserstein_1d(cl$ev_ref, ev_c)
    burst_c <- metric_burstiness_from_events(ev_c)
    if (!is.finite(emd) || !is.finite(burst_c)) next
    burst_loss <- max(0, (cl$burst_ref - burst_c) / (cl$burst_ref + 1e-12))
    centroid_shift <- abs(cl$centroid_ref - mean(ev_c))
    w <- log1p(cl$n_tips) * sqrt(cl$n_events)
    w_sum <- w_sum + w
    emd_sum <- emd_sum + (w * emd)
    burst_loss_sum <- burst_loss_sum + (w * burst_loss)
    centroid_sum <- centroid_sum + (w * centroid_shift)
    matched <- matched + 1L
  }

  if (w_sum <= 0 || matched == 0L) return(NULL)

  mean_emd <- emd_sum / w_sum
  mean_burst_loss <- burst_loss_sum / w_sum
  mean_centroid <- centroid_sum / w_sum
  pulse_error <- (w_emd * mean_emd) + (w_burst_loss * mean_burst_loss) + (w_centroid * mean_centroid)

  global_est <- metric_event_times_relative(tr)
  global_emd <- metric_wasserstein_1d(global_ref, global_est)
  global_burst_est <- metric_burstiness_from_events(global_est)
  global_burst_loss <- max(0, (global_burst_ref - global_burst_est) / (global_burst_ref + 1e-12))
  global_error <- (0.35 * global_emd) + (0.65 * global_burst_loss)

  coverage <- matched / length(panel)
  selector_error <- ((1 - w_global) * pulse_error) + (w_global * global_error) + (coverage_penalty * (1 - coverage))
  selector_score <- 1 / (1 + selector_error)

  data.frame(
    matched_clades = matched,
    panel_clades = length(panel),
    coverage = coverage,
    mean_emd = mean_emd,
    mean_burst_loss = mean_burst_loss,
    mean_centroid_shift = mean_centroid,
    local_error = pulse_error,
    global_emd = global_emd,
    global_burst_loss = global_burst_loss,
    global_error = global_error,
    selector_error = selector_error,
    selector_score = selector_score,
    stringsAsFactors = FALSE
  )
}

metric_edge_signature <- function(tr) {
  n <- ape::Ntip(tr)
  all_tips <- ape::prop.part(tr)
  sig <- character(nrow(tr$edge))
  for (i in seq_len(nrow(tr$edge))) {
    ch <- tr$edge[i, 2]
    if (ch <= n) {
      sig[i] <- paste0("tip:", tr$tip.label[ch])
    } else {
      idx <- ch - n
      tips <- sort(tr$tip.label[all_tips[[idx]]])
      sig[i] <- paste0("clade:", paste(tips, collapse = "|"))
    }
  }
  sig
}

metric_rate_pair_metrics <- function(ref_tree, dated_tree) {
  ref_df <- data.frame(
    sig = metric_edge_signature(ref_tree),
    ref_branch = ref_tree$edge.length,
    stringsAsFactors = FALSE
  )
  est_sig <- metric_edge_signature(dated_tree)
  est_df <- data.frame(
    sig = est_sig,
    time_branch = dated_tree$edge.length,
    parent_node = dated_tree$edge[, 1],
    child_node = dated_tree$edge[, 2],
    stringsAsFactors = FALSE
  )

  merged <- merge(est_df, ref_df, by = "sig", all = FALSE)
  merged <- merged[is.finite(merged$time_branch) & merged$time_branch > 0 &
                     is.finite(merged$ref_branch) & merged$ref_branch > 0, , drop = FALSE]
  if (!nrow(merged)) {
    return(data.frame(
      rate_n_edges = 0L,
      rate_mean = NA_real_,
      rate_median = NA_real_,
      rate_cv = NA_real_,
      log_rate_sd = NA_real_,
      log_rate_iqr = NA_real_,
      extreme_rate_frac = NA_real_,
      parent_child_pairs = 0L,
      parent_child_jump_mean = NA_real_,
      parent_child_jump_q95 = NA_real_,
      rate_autocorr_spearman = NA_real_,
      rate_irregularity = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  merged$rate <- merged$ref_branch / merged$time_branch
  merged$log_rate <- log(merged$rate)

  node_to_sig <- setNames(est_sig, dated_tree$edge[, 2])
  parent_sig <- unname(node_to_sig[as.character(merged$parent_node)])
  rate_map <- setNames(merged$rate, merged$sig)
  parent_rate <- unname(rate_map[parent_sig])
  keep_pairs <- is.finite(parent_rate) & is.finite(merged$rate) & parent_rate > 0 & merged$rate > 0
  jump <- abs(log(merged$rate[keep_pairs]) - log(parent_rate[keep_pairs]))

  q <- stats::quantile(merged$log_rate, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  iqr_lr <- q[2] - q[1]
  if (is.finite(iqr_lr) && iqr_lr > 0) {
    lo <- q[1] - (1.5 * iqr_lr)
    hi <- q[2] + (1.5 * iqr_lr)
    extreme_frac <- mean(merged$log_rate < lo | merged$log_rate > hi)
  } else {
    extreme_frac <- 0
  }

  if (length(jump) >= 5) {
    autocorr <- suppressWarnings(stats::cor(log(parent_rate[keep_pairs]), log(merged$rate[keep_pairs]), method = "spearman"))
  } else {
    autocorr <- NA_real_
  }
  autocorr_penalty <- if (is.finite(autocorr)) {
    1 - max(autocorr, 0)
  } else if (stats::sd(merged$log_rate) < 1e-10) {
    0
  } else {
    1
  }

  data.frame(
    rate_n_edges = nrow(merged),
    rate_mean = mean(merged$rate),
    rate_median = stats::median(merged$rate),
    rate_cv = stats::sd(merged$rate) / mean(merged$rate),
    log_rate_sd = stats::sd(merged$log_rate),
    log_rate_iqr = stats::IQR(merged$log_rate),
    extreme_rate_frac = extreme_frac,
    parent_child_pairs = sum(keep_pairs),
    parent_child_jump_mean = if (length(jump)) mean(jump) else NA_real_,
    parent_child_jump_q95 = if (length(jump)) as.numeric(stats::quantile(jump, 0.95, names = FALSE)) else NA_real_,
    rate_autocorr_spearman = autocorr,
    rate_irregularity = stats::sd(merged$log_rate) +
      ifelse(length(jump), mean(jump), 0) +
      (2 * extreme_frac) +
      autocorr_penalty,
    stringsAsFactors = FALSE
  )
}

metric_fossil_gap_by_pairs <- function(dated_tree, calibrations, tol = 1e-4) {
  age_by_node <- metric_node_ages(dated_tree)
  detail <- calibrations
  detail$tree_mrca <- mapply(function(a, b) ape::getMRCA(dated_tree, c(a, b)), detail$taxonA, detail$taxonB)
  detail$node_age <- as.numeric(age_by_node[as.character(detail$tree_mrca)])
  detail$ghost_gap_ma_raw <- detail$node_age - detail$age_min
  detail$ghost_gap_ma_raw[abs(detail$ghost_gap_ma_raw) < tol] <- 0
  is_point <- is.finite(detail$age_max) & abs(detail$age_max - detail$age_min) < tol
  detail$ghost_gap_ma <- ifelse(is_point, abs(detail$ghost_gap_ma_raw), pmax(0, detail$ghost_gap_ma_raw))
  detail$ghost_relmin <- detail$ghost_gap_ma / detail$age_min
  detail$min_violation_ma <- pmax(0, detail$age_min - detail$node_age - tol)
  detail$max_violation_ma <- pmax(0, detail$node_age - detail$age_max - tol)

  total_min <- sum(detail$age_min, na.rm = TRUE)
  violation_sum <- sum((detail$min_violation_ma + detail$max_violation_ma)[!is_point], na.rm = TRUE)
  if (all(is_point)) {
    burden <- mean(detail$ghost_relmin, na.rm = TRUE)
    gap_mode <- "point_calibration_slack"
  } else if (all(!is_point)) {
    burden <- mean(detail$ghost_relmin, na.rm = TRUE) + (5 * violation_sum / total_min)
    gap_mode <- "minimum_window_gap"
  } else {
    burden <- mean(detail$ghost_relmin, na.rm = TRUE) + (5 * violation_sum / total_min)
    gap_mode <- "mixed_point_and_window"
  }
  data.frame(
    fossil_n_calibrations = nrow(detail),
    fossil_n_missing_node_age = sum(!is.finite(detail$node_age)),
    ghost_sum_ma = sum(detail$ghost_gap_ma, na.rm = TRUE),
    ghost_mean_ma = mean(detail$ghost_gap_ma, na.rm = TRUE),
    ghost_median_ma = stats::median(detail$ghost_gap_ma, na.rm = TRUE),
    ghost_rel_mean = mean(detail$ghost_relmin, na.rm = TRUE),
    ghost_rel_median = stats::median(detail$ghost_relmin, na.rm = TRUE),
    min_violation_count = sum(detail$min_violation_ma > 0, na.rm = TRUE),
    min_violation_sum_ma = sum(detail$min_violation_ma, na.rm = TRUE),
    max_violation_count = sum(detail$max_violation_ma > 0, na.rm = TRUE),
    max_violation_sum_ma = sum(detail$max_violation_ma, na.rm = TRUE),
    gap_mode = gap_mode,
    fossil_gap_burden = burden,
    stringsAsFactors = FALSE
  )
}

select_subset_tips <- function(phy, cal_pairs, subset_n = 400L, extreme_frac = 0.05, seed = 1L) {
  nt <- Ntip(phy)
  if (!is.finite(subset_n) || subset_n < 2L) stop("SUBSET_N must be >= 2")
  subset_n <- as.integer(subset_n)
  if (subset_n >= nt) return(phy$tip.label)
  set.seed(seed)

  all_tips <- phy$tip.label
  keep_fixed <- character(0)
  if (!is.null(cal_pairs) && nrow(cal_pairs)) {
    keep_fixed <- unique(c(as.character(cal_pairs$taxonA), as.character(cal_pairs$taxonB)))
    keep_fixed <- keep_fixed[keep_fixed %in% all_tips]
  }
  if (length(keep_fixed) > subset_n) {
    stop("Number of calibration taxa (", length(keep_fixed), ") exceeds SUBSET_N (", subset_n, "). Increase SUBSET_N.")
  }

  # Root-to-tip extremes from the phylogram.
  rt <- node.depth.edgelength(phy)[seq_len(nt)]
  names(rt) <- all_tips
  k_ext <- max(2L, floor(subset_n * extreme_frac))
  ord <- names(sort(rt))
  ext <- unique(c(head(ord, k_ext), tail(ord, k_ext)))
  ext <- ext[ext %in% all_tips]

  # Diversified spread via ladderized tip order.
  spread <- ladderize(phy)$tip.label
  keep <- unique(c(keep_fixed, ext))
  if (length(keep) < subset_n) {
    pool <- setdiff(spread, keep)
    if (length(pool)) {
      idx <- unique(round(seq(1, length(pool), length.out = min(subset_n - length(keep), length(pool)))))
      keep <- c(keep, pool[idx])
    }
  }
  if (length(keep) < subset_n) {
    pool2 <- setdiff(all_tips, keep)
    keep <- c(keep, sample(pool2, min(length(pool2), subset_n - length(keep)), replace = FALSE))
  }
  unique(keep)[seq_len(min(subset_n, length(unique(keep))))]
}

run_chronos_modelselect <- function(phy, calib, clock_switch_thresh = PLOG_CLOCK_SWITCH_THRESH) {
  best_phiic <- list(tree = NULL, score = Inf, ploglik = NA_real_, model = NA_character_, lambda = NA_real_, nb_rate_cat = NA_integer_, selector_note = "phiic_default")
  best_model_phiic <- setNames(vector("list", length(CHRONOS_MODELS)), CHRONOS_MODELS)
  best_model_ploglik <- setNames(vector("list", length(CHRONOS_MODELS)), CHRONOS_MODELS)

  eval_one <- function(model, lam, nb_rate_cat = NA_integer_) {
    out <- list(tree = NULL, score = Inf, ploglik = NA_real_)
    ctrl <- chronos.control()
    if (model == "discrete" && is.finite(nb_rate_cat)) ctrl$nb.rate.cat <- as.integer(nb_rate_cat)
    for (k in seq_len(N_RETRIES)) {
      tr <- try(chronos(phy, lambda = lam, model = model, calibration = calib, quiet = TRUE, control = ctrl), silent = TRUE)
      if (inherits(tr, "try-error")) next
      sc <- fit_score_chronos(tr)
      pl <- attr(tr, "ploglik")
      if (!is.finite(pl)) pl <- NA_real_
      if (is.finite(sc) && sc < out$score) out <- list(tree = tr, score = sc, ploglik = pl)
    }
    out
  }

  for (mdl in CHRONOS_MODELS) {
    k_grid <- if (mdl == "discrete") K_FIT_GRID else NA_integer_
    if (mdl != "discrete") k_grid <- NA_integer_
    for (kcat in k_grid) {
      for (lam in LAMBDA_GRID) {
        ev <- eval_one(mdl, lam, kcat)
        if (is.finite(ev$score) && ev$score < best_phiic$score) {
          best_phiic <- list(tree = ev$tree, score = ev$score, ploglik = ev$ploglik, model = mdl, lambda = lam, nb_rate_cat = kcat, selector_note = "phiic_default")
        }
        if (is.finite(ev$score)) {
          cur_phi <- best_model_phiic[[mdl]]
          if (is.null(cur_phi) || !is.finite(cur_phi$score) || ev$score < cur_phi$score) {
            best_model_phiic[[mdl]] <- list(tree = ev$tree, score = ev$score, ploglik = ev$ploglik, model = mdl, lambda = lam, nb_rate_cat = kcat)
          }
        }
        if (is.finite(ev$ploglik)) {
          cur <- best_model_ploglik[[mdl]]
          if (is.null(cur) || !is.finite(cur$score) || ev$ploglik > cur$score) {
            best_model_ploglik[[mdl]] <- list(tree = ev$tree, score = ev$ploglik, phiic = ev$score, model = mdl, lambda = lam, nb_rate_cat = kcat)
          }
        }
      }
    }
  }

  # Robust selector: start PHIIC, then switch using c1_n2_t2 rules.
  robust <- best_phiic
  clock_pl <- best_model_ploglik[["clock"]]
  corr_pl <- best_model_ploglik[["correlated"]]
  rel_pl <- best_model_ploglik[["relaxed"]]
  disc_pl <- best_model_ploglik[["discrete"]]
  candidates <- list(corr_pl, rel_pl, disc_pl)
  candidates <- candidates[sapply(candidates, function(x) !is.null(x) && is.finite(x$score))]
  best_nonclock <- if (length(candidates)) candidates[[which.max(sapply(candidates, function(x) x$score))]] else NULL

  if (identical(robust$model, "clock") &&
      !is.null(clock_pl) && is.finite(clock_pl$score) &&
      !is.null(best_nonclock) && is.finite(best_nonclock$score) &&
      (best_nonclock$score - clock_pl$score) >= clock_switch_thresh) {
    if (!is.null(corr_pl) && !is.null(rel_pl) &&
        is.finite(corr_pl$score) && is.finite(rel_pl$score) &&
        best_nonclock$model %in% c("correlated", "relaxed") &&
        abs(corr_pl$score - rel_pl$score) <= PLOG_TIE_EPS) {
      robust <- list(tree = corr_pl$tree, score = corr_pl$phiic, ploglik = corr_pl$score,
                     model = corr_pl$model, lambda = corr_pl$lambda, nb_rate_cat = corr_pl$nb_rate_cat,
                     selector_note = "switch_nonclock_tie_to_correlated")
    } else {
      robust <- list(tree = best_nonclock$tree, score = best_nonclock$phiic, ploglik = best_nonclock$score,
                     model = best_nonclock$model, lambda = best_nonclock$lambda, nb_rate_cat = best_nonclock$nb_rate_cat,
                     selector_note = "switch_nonclock_ploglik_gain")
    }
  }
  if (identical(robust$model, "correlated") || identical(robust$model, "relaxed")) {
    if (!is.null(corr_pl) && !is.null(rel_pl) && is.finite(corr_pl$score) && is.finite(rel_pl$score)) {
      if ((rel_pl$score - corr_pl$score) > PLOG_NONCLOCK_SWITCH_THRESH) {
        robust <- list(tree = rel_pl$tree, score = rel_pl$phiic, ploglik = rel_pl$score,
                       model = rel_pl$model, lambda = rel_pl$lambda, nb_rate_cat = rel_pl$nb_rate_cat,
                       selector_note = "switch_corr_to_relaxed_ploglik")
      } else if ((corr_pl$score - rel_pl$score) > PLOG_NONCLOCK_SWITCH_THRESH) {
        robust <- list(tree = corr_pl$tree, score = corr_pl$phiic, ploglik = corr_pl$score,
                       model = corr_pl$model, lambda = corr_pl$lambda, nb_rate_cat = corr_pl$nb_rate_cat,
                       selector_note = "switch_relaxed_to_corr_ploglik")
      }
    }
  }
  robust$model_fits <- best_model_phiic
  robust$model_ploglik_fits <- best_model_ploglik
  robust
}

run_chronos_fixed <- function(phy, calib, model, lam, nb_rate_cat = NA_integer_) {
  best <- list(tree = NULL, score = Inf, ploglik = NA_real_)
  ctrl <- chronos.control()
  if (identical(model, "clock")) {
    model <- "clock"
  }
  if (identical(model, "discrete") && is.finite(nb_rate_cat)) {
    ctrl$nb.rate.cat <- as.integer(nb_rate_cat)
  }
  for (k in seq_len(N_RETRIES)) {
    tr <- try(chronos(phy, lambda = lam, model = model, calibration = calib, quiet = TRUE, control = ctrl), silent = TRUE)
    if (inherits(tr, "try-error")) next
    sc <- fit_score_chronos(tr)
    pl <- attr(tr, "ploglik")
    if (!is.finite(pl)) pl <- NA_real_
    if (is.finite(sc) && sc < best$score) best <- list(tree = tr, score = sc, ploglik = pl)
  }
  best
}

# -------------------------
# Main
# -------------------------
OUT_DIR <- OUT_BASE_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
MAIN_DIR <- file.path(OUT_DIR, "main_files")
ALL_DIR <- file.path(OUT_DIR, "all_files")
TABLES_DIR <- file.path(ALL_DIR, "tables")
TREES_DIR <- file.path(ALL_DIR, "trees")
LOGS_DIR <- file.path(ALL_DIR, "logs")
CHECKPOINTS_DIR <- file.path(ALL_DIR, "checkpoints")
for (d in c(MAIN_DIR, ALL_DIR, TABLES_DIR, TREES_DIR, LOGS_DIR, CHECKPOINTS_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

if (isTRUE(CLEAN_PREVIOUS_PREFIX_OUTPUTS)) {
  old_files <- list.files(OUT_DIR, recursive = TRUE, full.names = TRUE, all.files = FALSE, no.. = TRUE)
  if (length(old_files)) {
    old_files <- old_files[file.exists(old_files) & !dir.exists(old_files)]
  }
  if (length(old_files)) {
    archive_dir <- file.path(OUT_DIR, "_archive", format(Sys.time(), "%Y%m%d_%H%M%S"))
    dir.create(archive_dir, showWarnings = FALSE, recursive = TRUE)
    for (f in old_files) {
      rel <- sub(paste0("^", OUT_DIR, "/?"), "", f)
      dst <- file.path(archive_dir, rel)
      dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
      file.rename(f, dst)
    }
  }
}

log_file <- file.path(LOGS_DIR, paste0("run_", OUT_PREFIX, ".log"))
msg <- function(...) {
  s <- paste0(...)
  cat(s, "\n")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s), file = log_file, append = TRUE)
}

msg("Start: ", OUT_PREFIX)
sel <- load_target_tree(TARGET_TREE_FILE, TARGET_TREE_NAME, TARGET_TREE_INDEX)
target_tree <- ladderize(sel$tree)
target_id <- sel$tree_id
msg("Target selected: ", target_id, " | Ntip=", Ntip(target_tree))

if (USE_CONGRUIF) {
  ref_tree <- read_any_tree(REFERENCE_TIME_TREE)
  if (inherits(ref_tree, "multiPhylo")) ref_tree <- ref_tree[[1]]
  cal_pairs <- build_calib_pairs_congruif(target_tree, ref_tree)
  if (!is.finite(ROOT_AGE)) {
    if (!is.ultrametric(ref_tree)) ref_tree <- force.ultrametric(ref_tree, method = "extend")
    ROOT_AGE <- max(node.depth.edgelength(ref_tree)[1:Ntip(ref_tree)])
  }
  msg("Calibration mode: congruify | pairs=", nrow(cal_pairs), " | root_age=", format(ROOT_AGE, digits = 12))
} else {
  cal_pairs <- normalize_calib_df(read.csv(MANUAL_CAL_CSV, stringsAsFactors = FALSE))
  msg("Calibration mode: manual_csv | pairs=", nrow(cal_pairs))
}

target_tree_full <- target_tree
cal_pairs_full <- cal_pairs

subset_tip_file <- NA_character_
subset_mode_applied <- FALSE
n_tips_original <- Ntip(target_tree)
if (isTRUE(USE_SUBSET) && Ntip(target_tree) > SUBSET_N) {
  keep_tips <- select_subset_tips(
    phy = target_tree,
    cal_pairs = cal_pairs,
    subset_n = SUBSET_N,
    extreme_frac = SUBSET_EXTREME_FRAC,
    seed = SUBSET_SEED
  )
  drop_tips <- setdiff(target_tree$tip.label, keep_tips)
  target_tree <- drop.tip(target_tree, drop_tips)
  cal_pairs <- cal_pairs[
    cal_pairs$taxonA %in% target_tree$tip.label & cal_pairs$taxonB %in% target_tree$tip.label,
    , drop = FALSE
  ]
  if (!nrow(cal_pairs)) stop("Subset pruning removed all calibration pairs. Increase SUBSET_N.")
  subset_mode_applied <- TRUE
  subset_tip_file <- file.path(TABLES_DIR, paste0(gsub("[^A-Za-z0-9_]+", "_", target_id), "_subset_tips_used.csv"))
  write.csv(data.frame(tip = target_tree$tip.label, stringsAsFactors = FALSE), subset_tip_file, row.names = FALSE)
  msg("Subset mode: applied | original tips=", n_tips_original,
      " | retained tips=", Ntip(target_tree),
      " | calibration pairs retained=", nrow(cal_pairs),
      " | subset tip file=", subset_tip_file,
      " | subset_tune_only=", SUBSET_TUNE_ON_SUBSET_ONLY)
} else {
  msg("Subset mode: not applied | USE_SUBSET=", USE_SUBSET,
      " | original tips=", n_tips_original,
      " | threshold SUBSET_N=", SUBSET_N)
}

calib <- build_chronos_calib(target_tree, cal_pairs, ROOT_AGE)
msg("Final calibration nodes on target tree: ", length(calib$node))
calib_full <- if (subset_mode_applied && isTRUE(SUBSET_TUNE_ON_SUBSET_ONLY)) build_chronos_calib(target_tree_full, cal_pairs_full, ROOT_AGE) else NULL
if (!is.null(calib_full)) {
  msg("Full-tree calibration nodes (for final run after subset tuning): ", length(calib_full$node))
}

safe_id <- gsub("[^A-Za-z0-9_]+", "_", target_id)
cal_csv_file <- file.path(TABLES_DIR, paste0(safe_id, "_calibrations_used.csv"))
cal_csv_file_full <- if (subset_mode_applied && isTRUE(SUBSET_TUNE_ON_SUBSET_ONLY)) {
  file.path(TABLES_DIR, paste0(safe_id, "_calibrations_used_full_tree.csv"))
} else NA_character_
model_fits_file <- file.path(TABLES_DIR, paste0("summary_", OUT_PREFIX, "_model_fits.csv"))
postfit_metrics_file <- file.path(TABLES_DIR, paste0("summary_", OUT_PREFIX, "_postfit_metrics.csv"))
interpretation_file <- file.path(TABLES_DIR, paste0("interpretation_", OUT_PREFIX, ".txt"))
rds_file <- file.path(TABLES_DIR, paste0("results_", OUT_PREFIX, ".rds"))
ckpt_file <- file.path(CHECKPOINTS_DIR, paste0("checkpoint_", OUT_PREFIX, ".rds"))
phylogram_used_file <- file.path(TREES_DIR, paste0(safe_id, "_phylogram_used.tree"))
phylogram_used <- if (subset_mode_applied && isTRUE(SUBSET_TUNE_ON_SUBSET_ONLY)) target_tree_full else target_tree
write.tree(phylogram_used, phylogram_used_file)

# Remove previous model-tree duplicates so each run keeps one tree per model.
old_model_trees <- Sys.glob(file.path(TREES_DIR, paste0(safe_id, "_chronos_dated_model*.tre")))
old_selected_trees <- Sys.glob(file.path(TREES_DIR, paste0(safe_id, "_chronos_dated_clockthresh*.tre")))
if (length(old_model_trees)) unlink(old_model_trees, force = TRUE)
if (length(old_selected_trees)) unlink(old_selected_trees, force = TRUE)

thresh_grid <- unique(as.numeric(CLOCK_SWITCH_SENSITIVITY))
thresh_grid <- thresh_grid[is.finite(thresh_grid)]
if (!length(thresh_grid)) thresh_grid <- PLOG_CLOCK_SWITCH_THRESH

fit_rows <- vector("list", length(thresh_grid))
fit_list <- vector("list", length(thresh_grid))
names(fit_list) <- paste0("clock_thresh_", thresh_grid)
model_candidates <- list()
model_candidate_i <- 1L

for (i in seq_along(thresh_grid)) {
  thr <- thresh_grid[i]
  fit_i <- run_chronos_modelselect(target_tree, calib, clock_switch_thresh = thr)
  if (is.null(fit_i$tree)) stop("chronos failed to return a dated tree at threshold=", thr)
  dated_tree_file_i <- file.path(TREES_DIR, paste0(safe_id, "_chronos_dated_clockthresh", thr, ".tre"))
  write.tree(fit_i$tree, dated_tree_file_i)

  for (mdl in CHRONOS_MODELS) {
    mfit <- fit_i$model_fits[[mdl]]
    mpl <- fit_i$model_ploglik_fits[[mdl]]
    model_candidates[[model_candidate_i]] <- list(
      threshold = thr,
      model = mdl,
      fit = mfit,
      fit_pl = mpl
    )
    model_candidate_i <- model_candidate_i + 1L
  }

  fit_list[[i]] <- fit_i
  fit_rows[[i]] <- data.frame(
    target_tree_file = TARGET_TREE_FILE,
    target_tree_id = target_id,
    n_tips = Ntip(target_tree),
    n_tips_original = n_tips_original,
    n_tips_tuning = Ntip(target_tree),
    n_tips_full = Ntip(target_tree_full),
    subset_mode = subset_mode_applied,
    subset_tune_on_subset_only = if (subset_mode_applied) SUBSET_TUNE_ON_SUBSET_ONLY else FALSE,
    subset_n = if (subset_mode_applied) SUBSET_N else NA_integer_,
    subset_extreme_frac = if (subset_mode_applied) SUBSET_EXTREME_FRAC else NA_real_,
    subset_seed = if (subset_mode_applied) SUBSET_SEED else NA_integer_,
    subset_tip_file = subset_tip_file,
    use_congruif = USE_CONGRUIF,
    reference_tree_file = if (USE_CONGRUIF) REFERENCE_TIME_TREE else "",
    root_age = ROOT_AGE,
    n_calib_pairs_input = nrow(cal_pairs),
    n_calib_pairs_full = nrow(cal_pairs_full),
    n_calib_nodes_final = length(calib$node),
    plog_clock_switch_thresh = thr,
    chronos_model = fit_i$model %||% "",
    chronos_lambda = fit_i$lambda %||% NA_real_,
    chronos_fit_score = fit_i$score %||% NA_real_,
    chronos_ploglik = fit_i$ploglik %||% NA_real_,
    chronos_nb_rate_cat = fit_i$nb_rate_cat %||% NA_integer_,
    selector_note = fit_i$selector_note %||% "",
    dated_tree_file = dated_tree_file_i,
    cal_csv_file = cal_csv_file,
    stringsAsFactors = FALSE
  )
}

write.csv(cal_pairs, cal_csv_file, row.names = FALSE)
if (!is.na(cal_csv_file_full)) {
  write.csv(cal_pairs_full, cal_csv_file_full, row.names = FALSE)
}
summary_sensitivity <- do.call(rbind, fit_rows)

# In subset-tuning mode, keep tuning tree outputs but run final selected settings on full tree.
summary_sensitivity$dated_tree_file_tuning <- summary_sensitivity$dated_tree_file
summary_sensitivity$dated_tree_file_full <- NA_character_
if (subset_mode_applied && isTRUE(SUBSET_TUNE_ON_SUBSET_ONLY)) {
  msg("Subset tuning mode: running final full-tree chronos with selected settings per threshold.")
  for (i in seq_len(nrow(summary_sensitivity))) {
    thr <- summary_sensitivity$plog_clock_switch_thresh[i]
    mdl <- as.character(summary_sensitivity$chronos_model[i])
    lam <- as.numeric(summary_sensitivity$chronos_lambda[i])
    kcat <- as.numeric(summary_sensitivity$chronos_nb_rate_cat[i])
    if (!is.finite(kcat)) kcat <- NA_integer_
    fit_full <- run_chronos_fixed(target_tree_full, calib_full, mdl, lam, kcat)
    if (is.null(fit_full$tree)) {
      warning("Full-tree chronos failed for threshold=", thr, " model=", mdl, " lambda=", lam, call. = FALSE)
      next
    }
    full_file <- file.path(TREES_DIR, paste0(safe_id, "_chronos_dated_fulltree_clockthresh", thr, ".tre"))
    write.tree(fit_full$tree, full_file)
    summary_sensitivity$dated_tree_file_full[i] <- full_file
    summary_sensitivity$dated_tree_file[i] <- full_file
  }
}

# Build one output tree per model (choose best PHIIC candidate across threshold runs).
model_rows <- list()
pulse_panel <- metric_build_pulse_panel(target_tree, min_tips = 8L, min_events = 4L)
if (!length(pulse_panel)) {
  msg("Post-fit pulse panel: no eligible clades found; pulse-preservation metrics will be NA.")
}
for (mdl in CHRONOS_MODELS) {
  cand <- Filter(function(x) identical(x$model, mdl) && !is.null(x$fit) && inherits(x$fit$tree, "phylo"), model_candidates)
  if (!length(cand)) next
  scores <- sapply(cand, function(x) x$fit$score %||% Inf)
  pick <- cand[[which.min(scores)]]
  out_tree <- file.path(TREES_DIR, paste0(safe_id, "_chronos_dated_model", mdl, ".tre"))
  write.tree(pick$fit$tree, out_tree)
  bt <- branching_tempo_error(target_tree, pick$fit$tree)
  pulse_default <- if (length(pulse_panel)) {
    metric_score_pulse_panel(
      target_tree, pick$fit$tree, pulse_panel,
      w_emd = 0.35, w_burst_loss = 0.55, w_centroid = 0.10,
      w_global = 0.20, coverage_penalty = 0.20
    )
  } else NULL
  pulse_burst <- if (length(pulse_panel)) {
    metric_score_pulse_panel(
      target_tree, pick$fit$tree, pulse_panel,
      w_emd = 0.20, w_burst_loss = 0.75, w_centroid = 0.05,
      w_global = 0.20, coverage_penalty = 0.20
    )
  } else NULL
  gap <- metric_fossil_gap_by_pairs(pick$fit$tree, cal_pairs)
  rate <- metric_rate_pair_metrics(target_tree, pick$fit$tree)
  model_rows[[length(model_rows) + 1L]] <- data.frame(
    target_tree_file = TARGET_TREE_FILE,
    target_tree_id = target_id,
    model = mdl,
    source_threshold = pick$threshold,
    best_phiic = pick$fit$score %||% NA_real_,
    best_phiic_lambda = pick$fit$lambda %||% NA_real_,
    best_phiic_nb_rate_cat = pick$fit$nb_rate_cat %||% NA_integer_,
    ploglik_at_best_phiic = pick$fit$ploglik %||% NA_real_,
    best_ploglik = if (is.null(pick$fit_pl)) NA_real_ else pick$fit_pl$score %||% NA_real_,
    best_ploglik_lambda = if (is.null(pick$fit_pl)) NA_real_ else pick$fit_pl$lambda %||% NA_real_,
    best_ploglik_nb_rate_cat = if (is.null(pick$fit_pl)) NA_integer_ else pick$fit_pl$nb_rate_cat %||% NA_integer_,
    tempo_n_nodes = bt$n_nodes,
    tempo_mae_all = bt$mae_all,
    tempo_mae_early_q75 = bt$mae_early_q75,
    tempo_median_early_q75 = bt$median_early_q75,
    pulse_default_selector_error = if (is.null(pulse_default)) NA_real_ else pulse_default$selector_error[1],
    pulse_default_selector_score = if (is.null(pulse_default)) NA_real_ else pulse_default$selector_score[1],
    pulse_default_coverage = if (is.null(pulse_default)) NA_real_ else pulse_default$coverage[1],
    pulse_default_mean_burst_loss = if (is.null(pulse_default)) NA_real_ else pulse_default$mean_burst_loss[1],
    pulse_burst_selector_error = if (is.null(pulse_burst)) NA_real_ else pulse_burst$selector_error[1],
    pulse_burst_selector_score = if (is.null(pulse_burst)) NA_real_ else pulse_burst$selector_score[1],
    pulse_burst_coverage = if (is.null(pulse_burst)) NA_real_ else pulse_burst$coverage[1],
    pulse_burst_mean_burst_loss = if (is.null(pulse_burst)) NA_real_ else pulse_burst$mean_burst_loss[1],
    gap_mode = gap$gap_mode[1],
    fossil_gap_burden = gap$fossil_gap_burden[1],
    ghost_mean_ma = gap$ghost_mean_ma[1],
    ghost_median_ma = gap$ghost_median_ma[1],
    log_rate_sd = rate$log_rate_sd[1],
    parent_child_jump_mean = rate$parent_child_jump_mean[1],
    extreme_rate_frac = rate$extreme_rate_frac[1],
    rate_autocorr_spearman = rate$rate_autocorr_spearman[1],
    rate_irregularity = rate$rate_irregularity[1],
    dated_tree_file = out_tree,
    stringsAsFactors = FALSE
  )
}
model_fits <- if (length(model_rows)) do.call(rbind, model_rows) else data.frame()
if (nrow(model_fits) > 0) {
  model_fits$tempo_delta_all_to_best <- model_fits$tempo_mae_all - min(model_fits$tempo_mae_all, na.rm = TRUE)
  model_fits$tempo_delta_early_to_best <- model_fits$tempo_mae_early_q75 - min(model_fits$tempo_mae_early_q75, na.rm = TRUE)
  model_fits$tempo_rank_all <- rank(model_fits$tempo_mae_all, ties.method = "min")
  model_fits$tempo_rank_early <- rank(model_fits$tempo_mae_early_q75, ties.method = "min")
  model_fits$tempo_composite <- model_fits$tempo_mae_all + model_fits$tempo_mae_early_q75
  model_fits$tempo_rank_composite <- rank(model_fits$tempo_composite, ties.method = "min")
  model_fits$rank_pulse_default <- rank(model_fits$pulse_default_selector_error, ties.method = "min", na.last = "keep")
  model_fits$rank_burst_loss <- rank(model_fits$pulse_default_mean_burst_loss, ties.method = "min", na.last = "keep")
  model_fits$rank_pulse_burst <- rank(model_fits$pulse_burst_selector_error, ties.method = "min", na.last = "keep")
  model_fits$rank_fossil_gap <- rank(model_fits$fossil_gap_burden, ties.method = "min", na.last = "keep")
  model_fits$rank_rate_irregularity <- rank(model_fits$rate_irregularity, ties.method = "min", na.last = "keep")
  model_fits$rank_pulse_family_mean <- rowMeans(cbind(
    model_fits$rank_pulse_default,
    model_fits$rank_burst_loss,
    model_fits$rank_pulse_burst
  ), na.rm = TRUE)
  model_fits$rank_pulse_family <- rank(model_fits$rank_pulse_family_mean, ties.method = "min", na.last = "keep")
  model_fits$rank_mean_3families <- rowMeans(cbind(
    model_fits$rank_pulse_family_mean,
    model_fits$rank_fossil_gap,
    model_fits$rank_rate_irregularity
  ), na.rm = TRUE)
  model_fits$rank_mean_3families_rank <- rank(model_fits$rank_mean_3families, ties.method = "min", na.last = "keep")
  model_fits$rank_mean_4families <- rowMeans(cbind(
    model_fits$rank_pulse_default,
    model_fits$rank_pulse_burst,
    model_fits$rank_fossil_gap,
    model_fits$rank_rate_irregularity
  ), na.rm = TRUE)
  model_fits$rank_mean_4families_rank <- rank(model_fits$rank_mean_4families, ties.method = "min", na.last = "keep")
}
write.csv(model_fits, model_fits_file, row.names = FALSE)

postfit_metrics <- if (nrow(model_fits) > 0) {
  data.frame(
    candidate = paste0("chronos_", model_fits$model),
    model = model_fits$model,
    dated_tree_file = model_fits$dated_tree_file,
    pulse_default_selector_error = model_fits$pulse_default_selector_error,
    burst_loss = model_fits$pulse_default_mean_burst_loss,
    pulse_burst_selector_error = model_fits$pulse_burst_selector_error,
    gap_mode = model_fits$gap_mode,
    fossil_gap_burden = model_fits$fossil_gap_burden,
    rate_irregularity = model_fits$rate_irregularity,
    rank_pulse_default = model_fits$rank_pulse_default,
    rank_burst_loss = model_fits$rank_burst_loss,
    rank_pulse_burst = model_fits$rank_pulse_burst,
    rank_pulse_family_mean = model_fits$rank_pulse_family_mean,
    rank_pulse_family = model_fits$rank_pulse_family,
    rank_fossil_gap = model_fits$rank_fossil_gap,
    rank_rate_irregularity = model_fits$rank_rate_irregularity,
    rank_mean_3families = model_fits$rank_mean_3families,
    rank_mean_3families_rank = model_fits$rank_mean_3families_rank,
    rank_mean_4families = model_fits$rank_mean_4families,
    rank_mean_4families_rank = model_fits$rank_mean_4families_rank,
    stringsAsFactors = FALSE
  )
} else {
  data.frame()
}
write.csv(postfit_metrics, postfit_metrics_file, row.names = FALSE)

pick_idx <- which(thresh_grid == PLOG_CLOCK_SWITCH_THRESH)[1]
if (!is.finite(pick_idx)) pick_idx <- 1L
summary_row <- summary_sensitivity[pick_idx, , drop = FALSE]

fav_by_thr <- paste(
  paste0("threshold ", summary_sensitivity$plog_clock_switch_thresh, ": ", summary_sensitivity$chronos_model),
  collapse = "; "
)
fav_default <- summary_row$chronos_model[1]
lambda_default <- summary_row$chronos_lambda[1]
if (nrow(model_fits) > 0) {
  tempo_pick_early <- model_fits[which.min(model_fits$tempo_mae_early_q75), , drop = FALSE]
  tempo_best_model_early <- tempo_pick_early$model[1]
  tempo_best_mae_early <- tempo_pick_early$tempo_mae_early_q75[1]
  tempo_pick_all <- model_fits[which.min(model_fits$tempo_mae_all), , drop = FALSE]
  tempo_best_model_all <- tempo_pick_all$model[1]
  tempo_best_mae_all <- tempo_pick_all$tempo_mae_all[1]
  tempo_pick_comp <- model_fits[which.min(model_fits$tempo_composite), , drop = FALSE]
  tempo_best_model_comp <- tempo_pick_comp$model[1]
  tempo_best_comp <- tempo_pick_comp$tempo_composite[1]
} else {
  tempo_best_model_early <- NA_character_
  tempo_best_mae_early <- NA_real_
  tempo_best_model_all <- NA_character_
  tempo_best_mae_all <- NA_real_
  tempo_best_model_comp <- NA_character_
  tempo_best_comp <- NA_real_
}
clock_near_tie <- NA
fit_selected_lambda_line <- if (is.finite(lambda_default)) {
  paste0("Within the fit-selected model (", fav_default, "), the selected lambda is ", format(lambda_default, digits = 6), ".")
} else {
  paste0("Within the fit-selected model (", fav_default, "), the selected lambda is not finite.")
}

lambda_lines <- c("Best lambda by PHIIC within each model:")
if (nrow(model_fits) > 0) {
  ord_lambda <- order(match(model_fits$model, CHRONOS_MODELS))
  mf_lambda <- model_fits[ord_lambda, , drop = FALSE]
  for (i in seq_len(nrow(mf_lambda))) {
    lambda_lines <- c(
      lambda_lines,
      paste0(
        " - ", mf_lambda$model[i],
        ": lambda=", format(mf_lambda$best_phiic_lambda[i], digits = 6),
        ifelse(is.finite(mf_lambda$best_phiic_nb_rate_cat[i]),
               paste0(" | nb_rate_cat=", mf_lambda$best_phiic_nb_rate_cat[i]), ""),
        " | PHIIC=", format(mf_lambda$best_phiic[i], digits = 8),
        " | ploglik=", format(mf_lambda$ploglik_at_best_phiic[i], digits = 8)
      )
    )
  }
}

postfit_best_default <- if (nrow(postfit_metrics)) postfit_metrics$candidate[which.min(postfit_metrics$pulse_default_selector_error)] else NA_character_
postfit_best_burst <- if (nrow(postfit_metrics)) postfit_metrics$candidate[which.min(postfit_metrics$pulse_burst_selector_error)] else NA_character_
postfit_best_burst_loss <- if (nrow(postfit_metrics)) postfit_metrics$candidate[which.min(postfit_metrics$burst_loss)] else NA_character_
postfit_best_gap <- if (nrow(postfit_metrics)) postfit_metrics$candidate[which.min(postfit_metrics$fossil_gap_burden)] else NA_character_
postfit_best_rate <- if (nrow(postfit_metrics)) postfit_metrics$candidate[which.min(postfit_metrics$rate_irregularity)] else NA_character_
postfit_best_overall <- if (nrow(postfit_metrics)) postfit_metrics$candidate[which.min(postfit_metrics$rank_mean_3families)] else NA_character_
postfit_gap_mode <- if (nrow(postfit_metrics)) unique(postfit_metrics$gap_mode) else NA_character_

recommended_model <- fav_default
recommendation_reason <- if (isTRUE(!is.na(postfit_best_overall) && identical(paste0("chronos_", fav_default), postfit_best_overall))) {
  "fit_and_postfit_agree"
} else {
  "fit_and_postfit_differ_report_both"
}

postfit_lines <- c(
  "Leaders by post-fit metric family (lower is better):",
  paste0(" - pulse preservation (default): ", postfit_best_default),
  paste0(" - burst loss: ", postfit_best_burst_loss),
  paste0(" - pulse preservation (burst): ", postfit_best_burst),
  paste0(" - gap burden (", postfit_gap_mode, "): ", postfit_best_gap),
  paste0(" - rate plausibility: ", postfit_best_rate),
  paste0(" - overall family-balanced mean rank: ", postfit_best_overall),
  "   pulse family contributes one-third of this overall rank by averaging the default pulse selector, burst loss, and the burst-priority pulse selector.",
  "",
  "Per-model post-fit comparison (lower is better):"
)
if (nrow(postfit_metrics) > 0) {
  ord_postfit <- order(postfit_metrics$rank_mean_3families_rank, postfit_metrics$rank_mean_3families)
  pf <- postfit_metrics[ord_postfit, , drop = FALSE]
  for (i in seq_len(nrow(pf))) {
    postfit_lines <- c(
      postfit_lines,
      paste0(
        " - ", pf$candidate[i],
        " | pulse_default=", format(pf$pulse_default_selector_error[i], digits = 6),
        " | burst_loss=", format(pf$burst_loss[i], digits = 6),
        " | pulse_burst=", format(pf$pulse_burst_selector_error[i], digits = 6),
        " | gap=", format(pf$fossil_gap_burden[i], digits = 6),
        " | rate=", format(pf$rate_irregularity[i], digits = 6),
        " | family_mean_rank=", format(pf$rank_mean_3families[i], digits = 4),
        " (rank ", pf$rank_mean_3families_rank[i], ")"
      )
    )
  }
}

interp <- c(
  if (subset_mode_applied && isTRUE(SUBSET_TUNE_ON_SUBSET_ONLY)) {
    paste0("Subset tuning mode was used: tuning on ", Ntip(target_tree), " tips; final dated tree(s) generated on full tree with ", n_tips_original, " tips.")
  } else {
    paste0("Full-tree mode was used: all tuning and final dating ran on ", Ntip(target_tree), " tips.")
  },
  "",
  "Clock fitting",
  paste0("The fit selector (default threshold = ", PLOG_CLOCK_SWITCH_THRESH, ") favors the ", fav_default, " model."),
  if (nrow(summary_sensitivity) >= 2 &&
      !is.na(summary_sensitivity$chronos_model[summary_sensitivity$plog_clock_switch_thresh == 1][1]) &&
      !is.na(summary_sensitivity$chronos_model[summary_sensitivity$plog_clock_switch_thresh == 2][1]) &&
      identical(summary_sensitivity$chronos_model[summary_sensitivity$plog_clock_switch_thresh == 1][1],
                summary_sensitivity$chronos_model[summary_sensitivity$plog_clock_switch_thresh == 2][1])) {
    paste0("Using threshold 1 and threshold 2, the favored model is still ",
           summary_sensitivity$chronos_model[summary_sensitivity$plog_clock_switch_thresh == 1][1], ".")
  } else {
    paste0("Using threshold 1 and threshold 2, favored models are: ", fav_by_thr, ".")
  },
  "",
  "Lambda tuning",
  fit_selected_lambda_line,
  lambda_lines,
  "",
  "Post-fit evaluation",
  "This layer is separate from clock fitting and lambda tuning. It compares the resulting chronograms under pulse preservation, gap burden, and rate plausibility.",
  postfit_lines,
  "",
  paste0(
    "In short: clock fitting favors ", fav_default,
    ". Post-fit evaluation favors ", postfit_best_overall,
    " overall, while ", postfit_best_rate,
    " is best on rate plausibility. Because fit and post-fit do not fully agree here, report both layers explicitly rather than collapsing them into a single claim."
  )
)
writeLines(interp, interpretation_file)

summary_row$favored_models_by_threshold <- fav_by_thr
summary_row$tempo_best_model_all <- tempo_best_model_all
summary_row$tempo_best_mae_all <- tempo_best_mae_all
summary_row$tempo_best_model_early_q75 <- tempo_best_model_early
summary_row$tempo_best_mae_early_q75 <- tempo_best_mae_early
summary_row$tempo_best_model_composite <- tempo_best_model_comp
summary_row$tempo_best_composite <- tempo_best_comp
summary_row$fit_vs_tempo_all_agree <- isTRUE(!is.na(tempo_best_model_all) && identical(fav_default, tempo_best_model_all))
summary_row$fit_vs_tempo_early_agree <- isTRUE(!is.na(tempo_best_model_early) && identical(fav_default, tempo_best_model_early))
summary_row$clock_near_tie_with_tempo_best <- clock_near_tie
summary_row$tempo_eq_abs_tol <- TEMPO_EQ_ABS_TOL
summary_row$tempo_eq_rel_tol <- TEMPO_EQ_REL_TOL
summary_row$recommended_model <- recommended_model
summary_row$recommendation_reason <- recommendation_reason
summary_row$postfit_metrics_file <- postfit_metrics_file
summary_row$postfit_best_pulse_default <- postfit_best_default
summary_row$postfit_best_burst_loss <- postfit_best_burst_loss
summary_row$postfit_best_pulse_burst <- postfit_best_burst
summary_row$postfit_best_gap_burden <- postfit_best_gap
summary_row$postfit_best_rate_plausibility <- postfit_best_rate
summary_row$postfit_best_overall <- postfit_best_overall
summary_row$postfit_gap_mode <- postfit_gap_mode
summary_row$interpretation_file <- interpretation_file
summary_row$n_tips_original <- n_tips_original
summary_row$subset_mode <- subset_mode_applied
summary_row$subset_n <- if (subset_mode_applied) SUBSET_N else NA_integer_
summary_row$subset_extreme_frac <- if (subset_mode_applied) SUBSET_EXTREME_FRAC else NA_real_
summary_row$subset_seed <- if (subset_mode_applied) SUBSET_SEED else NA_integer_
summary_row$subset_tune_on_subset_only <- if (subset_mode_applied) SUBSET_TUNE_ON_SUBSET_ONLY else FALSE
summary_row$subset_tip_file <- subset_tip_file
summary_row$cal_csv_file_tuning <- cal_csv_file
summary_row$cal_csv_file_full <- cal_csv_file_full
summary_row$phylogram_used_file <- phylogram_used_file
summary_row$n_tips_tuning <- Ntip(target_tree)
summary_row$n_tips_final <- if (subset_mode_applied && isTRUE(SUBSET_TUNE_ON_SUBSET_ONLY)) Ntip(target_tree_full) else Ntip(target_tree)
summary_row$dated_tree_file_tuning <- if ("dated_tree_file_tuning" %in% names(summary_sensitivity)) {
  summary_sensitivity$dated_tree_file_tuning[pick_idx]
} else summary_row$dated_tree_file
summary_row$dated_tree_file_full <- if ("dated_tree_file_full" %in% names(summary_sensitivity)) {
  summary_sensitivity$dated_tree_file_full[pick_idx]
} else NA_character_

saveRDS(list(summary = summary_row, summary_sensitivity = summary_sensitivity, model_fits = model_fits, fits = fit_list, calib = calib, cal_pairs = cal_pairs), rds_file)
saveRDS(list(summary = summary_row, summary_sensitivity = summary_sensitivity, model_fits = model_fits), ckpt_file)

# Copy key deliverables into a compact main_files folder.
main_copy <- c(
  cal_csv_file, cal_csv_file_full, model_fits_file, postfit_metrics_file, interpretation_file, subset_tip_file, phylogram_used_file,
  file.path(TREES_DIR, paste0(safe_id, "_chronos_dated_modelclock.tre")),
  file.path(TREES_DIR, paste0(safe_id, "_chronos_dated_modelcorrelated.tre")),
  file.path(TREES_DIR, paste0(safe_id, "_chronos_dated_modelrelaxed.tre")),
  file.path(TREES_DIR, paste0(safe_id, "_chronos_dated_modeldiscrete.tre"))
)
main_copy <- unique(main_copy[file.exists(main_copy)])
for (f in main_copy) file.copy(f, file.path(MAIN_DIR, basename(f)), overwrite = TRUE)

msg("Completed sensitivity grid: ", paste(thresh_grid, collapse = ","))
msg("Selected default threshold row: ", summary_row$plog_clock_switch_thresh,
    " | model=", summary_row$chronos_model,
    " | lambda=", summary_row$chronos_lambda,
    " | score=", format(summary_row$chronos_fit_score, digits = 6))
msg("Saved: ", model_fits_file)
msg("Saved: ", postfit_metrics_file)
msg("Saved: ", interpretation_file)
msg("Saved: ", rds_file)
