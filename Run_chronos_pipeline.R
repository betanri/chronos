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
TARGET_TREE_NAME <- ""   # e.g., "TERAP_ASTRAL_S1" (leave "" to use TARGET_TREE_INDEX)
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

# Output
OUT_DIR <- "chronos_empirical_out"
OUT_PREFIX <- "terap_empirical"
CLEAN_PREVIOUS_PREFIX_OUTPUTS <- TRUE

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

run_chronos_modelselect <- function(phy, calib, clock_switch_thresh = PLOG_CLOCK_SWITCH_THRESH) {
  best_phiic <- list(tree = NULL, score = Inf, ploglik = NA_real_, model = NA_character_, lambda = NA_real_, nb_rate_cat = NA_integer_, selector_note = "phiic_default")
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
  robust
}

# -------------------------
# Main
# -------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
TABLES_DIR <- file.path(OUT_DIR, "tables")
TREES_DIR <- file.path(OUT_DIR, "trees")
LOGS_DIR <- file.path(OUT_DIR, "logs")
CHECKPOINTS_DIR <- file.path(OUT_DIR, "checkpoints")
for (d in c(TABLES_DIR, TREES_DIR, LOGS_DIR, CHECKPOINTS_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

if (isTRUE(CLEAN_PREVIOUS_PREFIX_OUTPUTS)) {
  old_files <- c(
    Sys.glob(file.path(OUT_DIR, paste0("*", OUT_PREFIX, "*"))),
    Sys.glob(file.path(TABLES_DIR, paste0("*", OUT_PREFIX, "*"))),
    Sys.glob(file.path(TREES_DIR, paste0("*", OUT_PREFIX, "*"))),
    Sys.glob(file.path(LOGS_DIR, paste0("*", OUT_PREFIX, "*"))),
    Sys.glob(file.path(CHECKPOINTS_DIR, paste0("*", OUT_PREFIX, "*")))
  )
  old_files <- unique(old_files[file.exists(old_files)])
  if (length(old_files)) {
    archive_dir <- file.path(OUT_DIR, "_archive", format(Sys.time(), "%Y%m%d_%H%M%S"))
    dir.create(archive_dir, showWarnings = FALSE, recursive = TRUE)
    file.rename(old_files, file.path(archive_dir, basename(old_files)))
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

calib <- build_chronos_calib(target_tree, cal_pairs, ROOT_AGE)
msg("Final calibration nodes on target tree: ", length(calib$node))

safe_id <- gsub("[^A-Za-z0-9_]+", "_", target_id)
cal_csv_file <- file.path(TABLES_DIR, paste0(safe_id, "_calibrations_used.csv"))
summary_file <- file.path(TABLES_DIR, paste0("summary_", OUT_PREFIX, ".csv"))
summary_sensitivity_file <- file.path(TABLES_DIR, paste0("summary_", OUT_PREFIX, "_sensitivity.csv"))
rds_file <- file.path(TABLES_DIR, paste0("results_", OUT_PREFIX, ".rds"))
ckpt_file <- file.path(CHECKPOINTS_DIR, paste0("checkpoint_", OUT_PREFIX, ".rds"))

thresh_grid <- unique(as.numeric(CLOCK_SWITCH_SENSITIVITY))
thresh_grid <- thresh_grid[is.finite(thresh_grid)]
if (!length(thresh_grid)) thresh_grid <- PLOG_CLOCK_SWITCH_THRESH

fit_rows <- vector("list", length(thresh_grid))
fit_list <- vector("list", length(thresh_grid))
names(fit_list) <- paste0("clock_thresh_", thresh_grid)

for (i in seq_along(thresh_grid)) {
  thr <- thresh_grid[i]
  fit_i <- run_chronos_modelselect(target_tree, calib, clock_switch_thresh = thr)
  if (is.null(fit_i$tree)) stop("chronos failed to return a dated tree at threshold=", thr)
  dated_tree_file_i <- file.path(TREES_DIR, paste0(safe_id, "_chronos_dated_clockthresh", thr, ".tre"))
  write.tree(fit_i$tree, dated_tree_file_i)
  fit_list[[i]] <- fit_i
  fit_rows[[i]] <- data.frame(
    target_tree_file = TARGET_TREE_FILE,
    target_tree_id = target_id,
    n_tips = Ntip(target_tree),
    use_congruif = USE_CONGRUIF,
    reference_tree_file = if (USE_CONGRUIF) REFERENCE_TIME_TREE else "",
    root_age = ROOT_AGE,
    n_calib_pairs_input = nrow(cal_pairs),
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
summary_sensitivity <- do.call(rbind, fit_rows)
write.csv(summary_sensitivity, summary_sensitivity_file, row.names = FALSE)

pick_idx <- which(thresh_grid == PLOG_CLOCK_SWITCH_THRESH)[1]
if (!is.finite(pick_idx)) pick_idx <- 1L
summary_row <- summary_sensitivity[pick_idx, , drop = FALSE]
write.csv(summary_row, summary_file, row.names = FALSE)

saveRDS(list(summary = summary_row, summary_sensitivity = summary_sensitivity, fits = fit_list, calib = calib, cal_pairs = cal_pairs), rds_file)
saveRDS(list(summary = summary_row, summary_sensitivity = summary_sensitivity), ckpt_file)

msg("Completed sensitivity grid: ", paste(thresh_grid, collapse = ","))
msg("Selected default threshold row: ", summary_row$plog_clock_switch_thresh,
    " | model=", summary_row$chronos_model,
    " | lambda=", summary_row$chronos_lambda,
    " | score=", format(summary_row$chronos_fit_score, digits = 6))
msg("Saved: ", summary_file)
msg("Saved: ", summary_sensitivity_file)
msg("Saved: ", rds_file)
