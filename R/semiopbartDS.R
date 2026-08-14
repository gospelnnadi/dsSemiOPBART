# dsSemiOPBARTBase.R
# ---------------------------------------------------------------------------
# ARCHITECTURE D, REDESIGNED FROM PARTITION+ROTATION TO SWAP:
#
# Previously: ONE global ensemble of num_tree trees, PARTITIONED across
# sites (site1 owns 4, site2 owns 7, ...), each site's "shadow forest"
# reconstructed the other sites' trees via broadcast so f(x) at any site
# = owned_forest$do_predict + shadow_forest$do_predict.
#
# Now: EVERY site keeps its OWN COMPLETE, permanent num_tree-tree forest,
# always (no partition, no shadow forest, no reconstruction needed for
# prediction -- f(x) at a site is just owned_forest$do_predict(X)).
# Periodically (every swap_every sweeps), each site exchanges a FIXED
# n_swap of its own current trees with a neighbor along a ring topology
# (ds_semiopbart.R's swap_trees()) -- a strict swap, so every site's tree
# count stays exactly num_tree at all times. The federation's aggregate
# tree count is therefore num_tree * n_sites, but each site's own forest
# is independently self-sufficient for prediction there -- nothing about
# another site's data or trees needs to be present anywhere to predict.
#
# "Except the received trees": SoftBart's do_gibbs()/do_gibbs_weighted()
# always backfits every tree in the forest -- there is no subset-grow
# primitive in the C++ layer. Excluding just-received trees from a grow
# sweep is done in pure R instead, using the already-patched get_trees()/
# set_trees(): snapshot the received slots before growing, grow the whole
# forest as normal, then restore exactly those slots to their pre-grow
# state. The freeze lasts exactly ONE grow call (the sweep right after
# receipt) -- from the sweep after that, a received tree resumes normal
# local growth like any other, until it's next selected for a swap-out
# (which can happen from ANY current slot, not just "originally grown
# here" ones -- so trees can migrate through more than one site over
# time, not just bounce between an original pair).
# ---------------------------------------------------------------------------

library(usethis)
library(devtools)
#devtools::load_all("SoftBart")
#devtools::load_all("~/Datashield_package_dev/FedBART/dsSemiOPBART/SoftBart")
library(truncnorm)
#library(SoftBart)

#' @name semiOPBARTLocalInitDS
#' One-time local setup: allocate this site's OWN, COMPLETE num_tree-tree
#' forest from the already-prepared TRAIN object. No ownership partition
#' argument anymore -- every site always gets the full num_tree trees.
#'
#' @param data.name  the TRAIN object (from ds.semiOPBARTSplit())
#' @param num_tree   this (and every) site's own local forest size
#' @param k          SoftBart leaf-shrinkage hyperparameter
#' @param allow_local_ecdf  FALSE (default): refuses to run under
#'   normalize_method = "local_ecdf" -- see the note below for why. TRUE:
#'   proceeds anyway, with a loud warning instead of a stop -- an
#'   explicit, informed override, not a silent relaxation of the check.
#' @export
semiOPBARTLocalInitDS <- function(data.name, num_tree = 20, k = 1, nfilter = 5,
                                  state_name = ".semiOPBART_state",
                                  allow_local_ecdf = FALSE,
                                  seed = 35) {
  set.seed(seed)
  s <- get(data.name, envir = parent.frame())
  if (length(s$Y) < nfilter) stop("site n below disclosure threshold")

  # NOT a hard requirement anymore -- an OPT-IN override, off by default.
  # Trees are exchanged between sites via set_trees() (see
  # semiOPBARTLocalSwapDS()). A tree's cutpoint (`val`, in [0,1]) only
  # refers to the same raw covariate value at every site if every site
  # normalized onto the SAME scale. Under "local_ecdf", each site fits its
  # own ECDF from its own training data -- a cutpoint of 0.7 grown at one
  # site means a different raw value at another, so a swapped-in tree
  # evaluated there is evaluating a DIFFERENT covariate distribution than
  # the one it was grown against, not an approximation of the same one.
  # This can quietly degrade accuracy (a swapped tree's splits land in the
  # wrong place relative to the receiving site's own data) without ever
  # producing an error -- allow_local_ecdf = TRUE accepts that risk
  # explicitly, it does not fix it. "federated_minmax" and
  # "federated_ecdf" remain the only two methods that keep tree-sharing
  # actually coherent (both are a single shared, public function applied
  # identically everywhere).
  norm_method <- s$norm_info$method
  local_ecdf_used <- is.null(norm_method) || !(norm_method %in% c("federated_minmax", "federated_ecdf"))
  if (local_ecdf_used && !allow_local_ecdf) {
    stop("Architecture D defaults to requiring normalize_method = ",
         "'federated_minmax' or 'federated_ecdf' -- this site's train ",
         "object was normalized with '",
         if (is.null(norm_method)) "unknown" else norm_method, "'. Trees are ",
         "swapped between sites in this architecture, and a swapped tree's ",
         "cutpoints are only meaningful if every site normalizes onto the ",
         "same scale; local_ecdf does not guarantee that. Either re-run ",
         "ds.semiOPBARTPrepare() with normalize_method = 'federated_minmax' ",
         "or 'federated_ecdf', or pass allow_local_ecdf = TRUE to proceed ",
         "anyway with that risk accepted explicitly.")
  } else if (local_ecdf_used) {
    message("semiOPBARTLocalInitDS(): proceeding under normalize_method = '",
            if (is.null(norm_method)) "unknown" else norm_method, "' with ",
            "allow_local_ecdf = TRUE -- swapped trees' cutpoints will NOT ",
            "mean the same covariate value at every site. This can silently ",
            "reduce accuracy; it will not produce an error anywhere.")
  }

  hypers_owned <- SoftBart::Hypers(X = s$X, Y = as.numeric(s$Y))
  hypers_owned$sigma_mu <- 3 / k / sqrt(num_tree)
  hypers_owned$sigma     <- 1     # fixed, probit identifiability
  hypers_owned$sigma_hat <- 1
  hypers_owned$num_tree  <- num_tree
  # CONFIRMED against the real smopbart() source: without this, any
  # multi-level categorical x_feature's dummy columns would compete
  # independently in the Dirichlet variable-selection prior instead of
  # being treated as one logical variable, silently diverging from the
  # real algorithm.
  hypers_owned$group <- dummy_assign(s$dv)

  opts <- SoftBart::Opts()
  opts$update_sigma <- FALSE

  owned_forest <- SoftBart::MakeForest(hypers_owned, opts, FALSE)

  assign(state_name, list(
    X = s$X, Y = as.numeric(s$Y), W = s$W, dv = s$dv, norm_info = s$norm_info,
    num_tree = num_tree, k = k,
    hypers_owned = hypers_owned, opts = opts,
    owned_forest = owned_forest,
    frozen_idx = integer(0),   # LOCAL 0-based slots to skip on the NEXT
                                # grow call only -- see file header
    fx_local = rep(0, length(s$Y)), n = length(s$Y)
  ), envir = parent.frame())

  list(n = length(s$Y), p = ncol(s$W), J = nlevels(s$Y))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' @name semiOPBARTLocalInitEDS
#' Architecture E's init: every site owns ALL num_tree trees, permanently,
#' and NEVER swaps or exchanges any of them with any other site -- f(x) is
#' 100% local, forever. Drops the federated_minmax/federated_ecdf hard
#' requirement entirely: that requirement exists only because D exchanges
#' grown trees between sites, so a cutpoint has to mean the same raw value
#' everywhere. E never shares a tree with anyone, so there is no
#' shared-scale requirement to enforce here. Any normalize_method
#' (including local_ecdf) is fine.
#'
#' Deliberately assigns to the SAME ".semiOPBART_state" slot
#' semiOPBARTLocalInitDS() uses, so semiOPBARTLocalSyncDS()/
#' semiOPBARTLocalGrowDS()/semiOPBARTLocalThresholdStatsDS() below are
#' reused UNCHANGED for Architecture E's per-sweep loop. CAVEAT: this
#' means Architecture D and E must not be trained in the same Opal
#' session at the same site without re-initializing in between --
#' whichever ran last owns the slot.
#' @export
semiOPBARTLocalInitEDS <- function(data.name, num_tree = 20, k = 1, nfilter = 5,
                                   state_name = ".semiOPBART_state", seed = 35) {
  set.seed(seed)
  s <- get(data.name, envir = parent.frame())
  if (length(s$Y) < nfilter) stop("site n below disclosure threshold")

  hypers_owned <- SoftBart::Hypers(X = s$X, Y = as.numeric(s$Y))
  hypers_owned$sigma_mu <- 3 / k / sqrt(num_tree)
  hypers_owned$sigma     <- 1
  hypers_owned$sigma_hat <- 1
  hypers_owned$num_tree  <- num_tree
  hypers_owned$group     <- dummy_assign(s$dv)

  opts <- SoftBart::Opts()
  opts$update_sigma <- FALSE

  owned_forest <- SoftBart::MakeForest(hypers_owned, opts, FALSE)

  assign(state_name, list(
    X = s$X, Y = as.numeric(s$Y), W = s$W, dv = s$dv, norm_info = s$norm_info,
    num_tree = num_tree, k = k,
    hypers_owned = hypers_owned, opts = opts,
    owned_forest = owned_forest, frozen_idx = integer(0),
    fx_local = rep(0, length(s$Y)), n = length(s$Y)
  ), envir = parent.frame())

  list(n = length(s$Y), p = ncol(s$W), J = nlevels(s$Y))
}

#' @name semiOPBARTLocalSyncDS
#' Resample this site's latent Z against its OWN current forest, and
#' report this site's contribution to the pooled theta update.
#'
#' CONFIRMED against the real smopbart() source (theta_hat/theta_sigma
#' computed from Z_star = Z - fx_train[i-1,], i.e. the PREVIOUS sweep's
#' f(x), strictly BEFORE this sweep's tree growth): W_gram/W_resid here
#' are computed against s$fx_local as of right after this Z draw -- before
#' semiOPBARTLocalGrowDS() runs and changes it. Every site contributes to
#' theta here, every sweep -- theta is a pooled linear parameter with
#' nothing to do with which trees a site currently holds.
#'
#' @param theta_Serialize Serialize-encoded current global linear coefficients
#' @param us_Serialize    Serialize-encoded current global ordinal thresholds
#' @export
semiOPBARTLocalSyncDS <- function(theta_Serialize, us_Serialize,
                                  state_name = ".semiOPBART_state",
                                  seed = 35) {
  
  set.seed(seed)
  s <- get(state_name, envir = parent.frame())
  theta <- semiOPBART_fromSerialize(theta_Serialize)
  us    <- semiOPBART_fromSerialize(us_Serialize)

  lower_us <- c(-Inf, us); upper_us <- c(us, Inf)
  y_idx <- as.integer(s$Y)
  lower <- lower_us[y_idx]; upper <- upper_us[y_idx]

  s$fx_local <- as.numeric(s$owned_forest$do_predict(s$X))   # f(x) as of
                                     # BEFORE this sweep's growth -- matches
                                     # real smopbart()'s fx_train[i-1,] exactly

  if (any(lower>= upper)) {
    stop(
        "Invalid truncation interval: ",
        sum(lower >= upper),
        " observations."
    )
  }
  s$Z <- truncnorm::rtruncnorm(n = s$n, a = lower, b = upper,
                                mean = as.numeric(s$fx_local + s$W %*% theta), sd = 1)
  assign(state_name, s, envir = parent.frame())


  # resid <- s$Z - s$fx_local
  # cat("Site:", s$site_id, "\n")
  # cat("anyNA(theta):", anyNA(theta), "\n")
  # cat("anyNA(fx_local):", anyNA(s$fx_local), "\n")
  # cat("anyNA(Z):", anyNA(s$Z), "\n")
  # cat("anyNA(resid):", anyNA(resid), "\n")

  list(W_gram  = crossprod(s$W),
       W_resid = crossprod(s$W, s$Z - s$fx_local))
}

#' @name semiOPBARTLocalGrowDS
#' Grow/update this site's OWN forest -- every tree EXCEPT whichever local
#' slots were just overwritten by the previous swap event (s$frozen_idx),
#' which sit out exactly this one grow call -- then export the FULL
#' current forest (not a subset) via the patched get_trees() method, for
#' the client to use as this sweep's swap-export pool if a swap is due.
#'
#' @param theta_Serialize Serialize-encoded theta -- the JUST-UPDATED value from
#'   this same sweep's theta update (see semiOPBARTLocalSyncDS()'s
#'   docstring), matching real smopbart()'s
#'   Z_tilde = Z - W %*% theta_p[i,] (the new theta) exactly.
#' @export
semiOPBARTLocalGrowDS <- function(theta_Serialize, state_name = ".semiOPBART_state", seed = 35) {
  set.seed(seed)
  s <- get(state_name, envir = parent.frame())
  theta <- semiOPBART_fromSerialize(theta_Serialize)

  Z_tilde <- as.numeric(s$Z - s$W %*% theta)

  # snapshot the just-received (frozen) slots BEFORE growing -- do_gibbs()
  # has no subset-grow primitive, so it always touches every tree; "don't
  # grow the received trees this sweep" is implemented by growing
  # everything, then restoring exactly these slots to their pre-grow
  # state via the already-patched set_trees().
  frozen_snapshot <- if (length(s$frozen_idx))
    s$owned_forest$get_trees(s$frozen_idx) else NULL

  fx_new <- as.numeric(s$owned_forest$do_gibbs(s$X, Z_tilde, s$X, 1))

  if (!is.null(frozen_snapshot)) {
    s$owned_forest$set_trees(frozen_snapshot, s$frozen_idx)
    fx_new <- as.numeric(s$owned_forest$do_predict(s$X))   # reflects the
                                     # restored (frozen-excluded) forest,
                                     # not the pre-restore grow output
  }

  s$fx_local  <- fx_new
  s$frozen_idx <- integer(0)   # freeze lasts exactly one grow call -- see
                                # file header
  assign(state_name, s, envir = parent.frame())

  list(n = s$n,
       trees_Serialize = semiOPBART_treesToSerialize(s$owned_forest, 0:(s$num_tree - 1)))
}

#' @name semiOPBARTLocalSwapDS
#' Overwrite this site's OWN chosen local slots with trees received from a
#' swap partner (ds_semiopbart.R's swap_trees()), and freeze those exact
#' slots for the NEXT grow call only. `target_idx` is THIS site's own
#' previously-chosen export indices (the ones it just gave away), so the
#' incoming trees land in the SAME positions that were vacated -- a
#' strict in-place swap, never changing this site's total tree count.
#'
#' @param target_idx_Serialize     Serialize-encoded LOCAL 0-based indices (this
#'   site's own chosen swap-out slots) to overwrite
#' @param incoming_trees_Serialize Serialize-encoded tree list, same length as
#'   target_idx, from the swap partner's most recent grow export
#' @export
semiOPBARTLocalSwapDS <- function(
    target_idx_Serialize,
    incoming_trees_Serialize,
    state_name = ".semiOPBART_state",
    seed = 35
) {
    set.seed(seed)

    s <- get(state_name, envir = parent.frame())

    target_idx <- semiOPBART_fromSerialize(
        target_idx_Serialize
    )

    if (!exists(
        "semiOPBART_treesFromSerialize",
        mode = "function",
        inherits = TRUE
    )) {
        stop(
            "semiOPBART_treesFromSerialize is not available ",
            "as a function on this DataSHIELD server."
        )
    }

    f <- get(
        "semiOPBART_treesFromSerialize",
        mode = "function",
        inherits = TRUE
    )

    f(
        s$owned_forest,
        incoming_trees_Serialize,
        target_idx
    )

    s$frozen_idx <- target_idx

    assign(
        state_name,
        s,
        envir = parent.frame()
    )

    list(
        n = s$n,
        n_swapped = length(target_idx)
    )
}


#' @name semiOPBARTLocalThresholdStatsDS
#' Local contribution to the ordered-probit threshold update.
#'
#' @param J  total number of ordinal categories (analyst-specified, from
#'   ds.semiOPBARTPrepare()'s `levels`) -- REQUIRED, not derived from
#'   `sort(unique(s$Y))`. Using locally-observed categories was a real bug:
#'   if one site never happened to observe some middle category, its
#'   result would have a different, shorter set of names than other
#'   sites', silently misaligning update_thresholds_from_site_stats()'s
#'   per-category loop across sites (it indexes by position, not by
#'   category value). Every site now always returns exactly J entries,
#'   NA for any category it didn't observe enough of (or at all).
#' @param threshold_method  "percentile" (default): disclosure-safe
#'   5th/95th percentile surrogate for the real algorithm's exact
#'   `max(Z[Y==j])`/`min(Z[Y==j+1])`. "exact": use the REAL algorithm's
#'   exact min/max, matching smopbart() bit-for-bit in what value it
#'   would compute -- but this releases a literal minimum or maximum
#'   latent value from this site's data, which is a much sharper
#'   disclosure risk than a percentile. Opt-in only.
#' @export
semiOPBARTLocalThresholdStatsDS <- function(nfilter = 5, J,
                                             threshold_method = "percentile",
                                             state_name = ".semiOPBART_state") {
  s <- get(state_name, envir = parent.frame())
  out <- lapply(seq_len(J), function(j) {
    zj <- s$Z[s$Y == j]
    if (length(zj) < nfilter) return(c(lo = NA, hi = NA, n = length(zj)))
    if (identical(threshold_method, "exact")) {
      c(lo = min(zj), hi = max(zj), n = length(zj))
    } else {
      c(lo = unname(quantile(zj, 0.05)), hi = unname(quantile(zj, 0.95)), n = length(zj))
    }
  })
  names(out) <- as.character(seq_len(J))
  out
}

# semiopbartDS.R - Add this new server function

#' @name semiOPBARTLocalMCMC
#' Run a FULL LOCAL MCMC loop on the server side for `num_local_sweeps`
#' iterations, with NO client communication in between. The client only
#' intervenes periodically to:
#'   1. Collect sufficient statistics for theta aggregation
#'   2. Exchange trees (Architecture D)
#'   3. Send back updated global theta/us
#'
#' This is the key function that reduces communication overhead -
#' the ENTIRE local Gibbs loop runs on the server, not iteration-by-iteration.
#'
#' @param num_local_sweeps  number of local Gibbs iterations to run
#'   before the next client intervention
#' @param theta_Serialize   Serialize-encoded global theta from last sync
#' @param us_Serialize      Serialize-encoded global us from last sync
#' @param perform_swap      logical - should tree swap happen after this block?
#' @param n_swap            number of trees to swap if perform_swap = TRUE
#' @param swap_targets_Serialize  Serialize-encoded local indices to swap out
#'   (if perform_swap = TRUE)
#' @param state_name        where to find/save local state
#' @export
semiOPBARTLocalMCMC <- function(num_local_sweeps, 
                                 theta_Serialize,
                                 us_Serialize,
                                 perform_swap = FALSE,
                                 n_swap = 0,
                                 swap_targets_Serialize = NULL,
                                 swap_incoming_Serialize = NULL,
                                 nfilter_threshold = 5,
                                 J,
                                 threshold_method = "exact",
                                 state_name = ".semiOPBART_state",
                                 seed = 35) {
  set.seed(seed)
  s <- get(state_name, envir = parent.frame())
  
  # Decode global parameters from last synchronization
  global_theta <- semiOPBART_fromSerialize(theta_Serialize)
  global_us <- semiOPBART_fromSerialize(us_Serialize)
  
  # Use local copies - these will evolve independently during local sweeps
  local_theta <- global_theta
  local_us <- global_us
  
  # Store sufficient statistics for aggregation at next sync
  WtW_accum <- matrix(0, ncol(s$W), ncol(s$W))
  WtZr_accum <- numeric(ncol(s$W))
  
  # Track Z and threshold stats for aggregation
  th_stats_list <- list()
  
  # ======================================================================
  # LOCAL GIBBS LOOP - runs ENTIRELY on the server, NO client communication
  # ======================================================================
  
  for (t in seq_len(num_local_sweeps)) {
    
    # ---- 1. Update latent Z ----
    lower_us <- c(-Inf, local_us)
    upper_us <- c(local_us, Inf)
    y_idx <- as.integer(s$Y)
    lower <- lower_us[y_idx]
    upper <- upper_us[y_idx]
    
    # f(x) from current forest
    fx_local <- as.numeric(s$owned_forest$do_predict(s$X))
    
    if (any(lower >= upper)) {
      stop("Invalid truncation interval for Z sampling")
    }
    
    Z <- truncnorm::rtruncnorm(n = s$n, a = lower, b = upper,
                                mean = as.numeric(fx_local + s$W %*% local_theta), 
                                sd = 1)
    
    # ---- 2. Accumulate sufficient statistics for theta ----
    # These will be sent to client for aggregation at next sync
    WtW_accum <- WtW_accum + crossprod(s$W)
    WtZr_accum <- WtZr_accum + crossprod(s$W, Z - fx_local)
    
    # ---- 3. Update theta locally (for the local Gibbs sweep) ----
    Z_star <- Z - fx_local
    theta_hat <- solve(crossprod(s$W)) %*% (crossprod(s$W, Z_star))
    theta_sigma <- solve(crossprod(s$W))
    local_theta <- as.numeric(mvtnorm::rmvnorm(1, theta_hat, theta_sigma))
    
    # ---- 4. Grow trees locally ----
    Z_tilde <- Z - s$W %*% local_theta
    s$owned_forest$do_gibbs(s$X, Z_tilde, s$X, 1)
    
    # ---- 5. Update thresholds locally ----
    # Accumulate stats for global threshold update
    th_stats <- lapply(seq_len(J), function(j) {
      zj <- Z[s$Y == j]
      if (length(zj) < nfilter_threshold) {
        return(c(lo = NA, hi = NA, n = length(zj)))
      }
      if (identical(threshold_method, "exact")) {
        c(lo = min(zj), hi = max(zj), n = length(zj))
      } else {
        c(lo = unname(quantile(zj, 0.05)), 
          hi = unname(quantile(zj, 0.95)), 
          n = length(zj))
      }
    })
    names(th_stats) <- as.character(seq_len(J))
    th_stats_list[[t]] <- th_stats
    
    # Local threshold update (using this site's own Z)
    # FIX: guard against length(local_us) == 1 (binary outcome, J = 2) --
    # 2:length(local_us) would otherwise be 2:1, a reversed/out-of-range
    # sequence. Previously harmless only because local_us[1] gets forced
    # back to 0 below regardless, and the nfilter_threshold guard happened
    # to skip the resulting out-of-range comparisons -- made explicit here.
    if (length(local_us) >= 2) {
      for (j in 2:length(local_us)) {
        z_j <- Z[s$Y == j]
        z_j1 <- Z[s$Y == j + 1]
        
        if (length(z_j) < nfilter_threshold || length(z_j1) < nfilter_threshold) {
          next
        }
        
        if (j == length(local_us)) {
          lower_bound <- max(max(z_j), local_us[j - 1])
          upper_bound <- min(min(z_j1), Inf)
        } else {
          lower_bound <- max(max(z_j), local_us[j - 1])
          upper_bound <- min(min(z_j1), local_us[j + 1])
        }
        
        if (is.finite(lower_bound) && is.finite(upper_bound) && 
            upper_bound > lower_bound) {
          local_us[j] <- stats::runif(1, lower_bound, upper_bound)
        }
      }
    }
    local_us[1] <- 0
    
    # ---- 6. Update the stored state ----
    s$Z <- Z
    s$fx_local <- fx_local
    assign(state_name, s, envir = parent.frame())
  }
  
  # ======================================================================
  # END OF LOCAL LOOP - now aggregate stats for client
  # ======================================================================
  
  # Aggregate threshold stats across all local sweeps
  # Use the MOST RECENT stats (or average - design choice)
  # Using most recent gives more weight to current state
  final_th_stats <- th_stats_list[[length(th_stats_list)]]
  
  # Export current forest for potential swapping
  trees_Serialize <- semiOPBART_treesToSerialize(s$owned_forest, 
                                                  0:(s$num_tree - 1))
  
  # Return sufficient statistics for client aggregation
  list(
    # For theta aggregation
    WtW = WtW_accum,
    WtZr = WtZr_accum,
    # For threshold aggregation
    th_stats = final_th_stats,
    # For tree swapping
    trees_Serialize = trees_Serialize,
    # Final local theta/us (after the last local sweep)
    local_theta = local_theta,
    local_us = local_us,
    n = s$n,
    num_sweeps = num_local_sweeps
  )
}




# semiopbartDS.R - Fixed server function

#' @name semiOPBARTLocalMCMCImproved
#' Enhanced local MCMC with support for adaptive communication patterns
#' and parameter-specific updates.
#' 
#' @param num_local_sweeps  number of local Gibbs iterations (must be integer)
#' @param theta_Serialize   Serialize-encoded global theta
#' @param us_Serialize      Serialize-encoded global us
#' @param update_theta      logical - should theta be updated locally?
#' @param update_us         logical - should us be updated locally?
#' @param perform_swap      logical - should tree swap happen?
#' @param n_swap            number of trees to swap
#' @param swap_targets_Serialize  local indices to swap out
#' @param swap_incoming_Serialize  incoming trees
#' @param state_name        where to find/save local state
#' @param adapt_rate        how much to adapt theta toward local MLE (0-1)
#' @export
semiOPBARTLocalMCMCImproved <- function(num_local_sweeps, 
                                         theta_Serialize,
                                         us_Serialize,
                                         update_theta = TRUE,
                                         update_us = TRUE,
                                         perform_swap = FALSE,
                                         n_swap = 0,
                                         swap_targets_Serialize = NULL,
                                         swap_incoming_Serialize = NULL,
                                         nfilter_threshold = 5,
                                         J,
                                         threshold_method = "percentile",
                                         state_name = ".semiOPBART_state",
                                         adapt_rate = 0.1,
                                         seed = 35) {
  set.seed(seed)
  
  # ---- FIX: Ensure num_local_sweeps is an integer ----
  num_local_sweeps <- as.integer(round(num_local_sweeps))
  
  # If num_local_sweeps is 0 or negative, do nothing
  if (num_local_sweeps <= 0) {
    s <- get(state_name, envir = parent.frame())
    trees_Serialize <- semiOPBART_treesToSerialize(s$owned_forest, 0:(s$num_tree - 1))
    return(list(
      WtW = matrix(0, ncol(s$W), ncol(s$W)),
      WtZr = numeric(ncol(s$W)),
      th_stats = NULL,
      trees_Serialize = trees_Serialize,
      local_theta = numeric(ncol(s$W)),
      local_us = numeric(length(us_Serialize) - 1),
      theta_history = matrix(0, 1, ncol(s$W)),
      theta_change = 0,
      n = s$n,
      num_sweeps = 0
    ))
  }
  
  s <- get(state_name, envir = parent.frame())
  
  # Decode global parameters
  global_theta <- semiOPBART_fromSerialize(theta_Serialize)
  global_us <- semiOPBART_fromSerialize(us_Serialize)
  
  # Use local copies
  local_theta <- global_theta
  local_us <- global_us
  
  # Store sufficient statistics
  WtW_accum <- matrix(0, ncol(s$W), ncol(s$W))
  WtZr_accum <- numeric(ncol(s$W))
  th_stats_list <- list()
  
  # Track local theta evolution for adaptive communication
  theta_history <- matrix(NA_real_, num_local_sweeps, length(local_theta))
  theta_history[1, ] <- local_theta
  
  for (t in seq_len(num_local_sweeps)) {
    
    # ---- 1. Update latent Z ----
    lower_us <- c(-Inf, local_us)
    upper_us <- c(local_us, Inf)
    y_idx <- as.integer(s$Y)
    lower <- lower_us[y_idx]
    upper <- upper_us[y_idx]
    
    fx_local <- as.numeric(s$owned_forest$do_predict(s$X))
    
    if (any(lower >= upper)) {
      stop("Invalid truncation interval for Z sampling")
    }
    
    Z <- truncnorm::rtruncnorm(n = s$n, a = lower, b = upper,
                                mean = as.numeric(fx_local + s$W %*% local_theta), 
                                sd = 1)
    
    # ---- 2. Accumulate sufficient statistics ----
    WtW_accum <- WtW_accum + crossprod(s$W)
    WtZr_accum <- WtZr_accum + crossprod(s$W, Z - fx_local)
    
    # ---- 3. Update theta (if allowed) ----
    if (update_theta) {
      Z_star <- Z - fx_local
      # Add small regularization to avoid singular matrix
      WtW <- crossprod(s$W)
      diag(WtW) <- diag(WtW) + 1e-10
      theta_hat <- solve(WtW) %*% (crossprod(s$W, Z_star))
      theta_sigma <- solve(WtW)
      local_theta <- as.numeric(mvtnorm::rmvnorm(1, theta_hat, theta_sigma))
      
      # Adaptive mixing: blend local update with global to prevent drift
      # during long communication gaps
      if (t > num_local_sweeps * 0.5) {
        local_theta <- (1 - adapt_rate) * local_theta + 
                       adapt_rate * as.numeric(theta_hat)
      }
    }
    theta_history[t, ] <- local_theta
    
    # ---- 4. Grow trees ----
    # FIX: honor s$frozen_idx, the same protection semiOPBARTLocalGrowDS()
    # already implements for the single-sweep-per-communication path.
    # Trees just received via semiOPBARTLocalSwapDS() must sit out exactly
    # one grow call, or do_gibbs() (which always touches every tree - there
    # is no subset-grow primitive) silently overwrites them before they
    # ever contribute to a single do_predict() call, nullifying the swap.
    # Only sweep t==1 of a block can immediately follow a swap event (swap
    # happens client-side between blocks - see ds_semiopbart.R), so the
    # snapshot/restore only needs to run there; frozen_idx is cleared once
    # consumed so later sweeps in this same block grow normally.
    Z_tilde <- Z - s$W %*% local_theta
    frozen_snapshot <- if (t == 1 && length(s$frozen_idx))
      s$owned_forest$get_trees(s$frozen_idx) else NULL
    s$owned_forest$do_gibbs(s$X, Z_tilde, s$X, 1)
    if (!is.null(frozen_snapshot)) {
      s$owned_forest$set_trees(frozen_snapshot, s$frozen_idx)
      s$frozen_idx <- integer(0)
      fx_local <- as.numeric(s$owned_forest$do_predict(s$X))  # reflects the
                                       # restored (frozen-excluded) forest,
                                       # not the pre-restore grow output --
                                       # matches semiOPBARTLocalGrowDS()
    }
    
    # ---- 5. Update thresholds (if allowed) ----
    if (update_us) {
      # Accumulate stats for global threshold update
      th_stats <- lapply(seq_len(J), function(j) {
        zj <- Z[s$Y == j]
        if (length(zj) < nfilter_threshold) {
          return(c(lo = NA, hi = NA, n = length(zj)))
        }
        if (identical(threshold_method, "exact")) {
          c(lo = min(zj), hi = max(zj), n = length(zj))
        } else {
          c(lo = unname(quantile(zj, 0.05)), 
            hi = unname(quantile(zj, 0.95)), 
            n = length(zj))
        }
      })
      names(th_stats) <- as.character(seq_len(J))
      th_stats_list[[t]] <- th_stats
      
      # Local threshold update
      # FIX: guard against length(local_us) == 1 (binary outcome, J = 2) --
      # see the identical fix/comment in semiOPBARTLocalMCMC() above.
      if (length(local_us) >= 2) {
        for (j in 2:length(local_us)) {
          z_j <- Z[s$Y == j]
          z_j1 <- Z[s$Y == j + 1]
          
          if (length(z_j) < nfilter_threshold || length(z_j1) < nfilter_threshold) {
            next
          }
          
          if (j == length(local_us)) {
            lower_bound <- max(max(z_j), local_us[j - 1])
            upper_bound <- min(min(z_j1), Inf)
          } else {
            lower_bound <- max(max(z_j), local_us[j - 1])
            upper_bound <- min(min(z_j1), local_us[j + 1])
          }
          
          if (is.finite(lower_bound) && is.finite(upper_bound) && 
              upper_bound > lower_bound) {
            local_us[j] <- stats::runif(1, lower_bound, upper_bound)
          }
        }
      }
      local_us[1] <- 0
    }
    
    # ---- 6. Update state ----
    s$Z <- Z
    s$fx_local <- fx_local
    assign(state_name, s, envir = parent.frame())
  }
  
  # ---- Return results ----
  final_th_stats <- if (update_us && length(th_stats_list) > 0) {
    th_stats_list[[length(th_stats_list)]]
  } else {
    NULL
  }
  
  trees_Serialize <- semiOPBART_treesToSerialize(s$owned_forest, 
                                                  0:(s$num_tree - 1))
  
  # Compute theta change for adaptive scheduling
  if (num_local_sweeps > 1 && !any(is.na(theta_history[1, ])) && !any(is.na(theta_history[num_local_sweeps, ]))) {
    theta_change <- mean(abs(theta_history[num_local_sweeps, ] - 
                             theta_history[1, ]) / 
                          (abs(theta_history[1, ]) + 1e-10))
  } else {
    theta_change <- 0
  }
  
  list(
    WtW = WtW_accum,
    WtZr = WtZr_accum,
    th_stats = final_th_stats,
    trees_Serialize = trees_Serialize,
    local_theta = local_theta,
    local_us = local_us,
    theta_history = theta_history,
    theta_change = theta_change,
    n = s$n,
    num_sweeps = num_local_sweeps
  )
}

#' @name semiOPBARTLocalThresholdUpdateDS
#' Local-only threshold update for Architecture D's non-communication
#' iterations. Each site updates its thresholds using its OWN Z and
#' the current local us value, with NO cross-site aggregation.
#' 
#' This mirrors the real smopbart()'s threshold update but uses
#' ONLY local data, and does NOT return any statistics to the client.
#' 
#' @param J  total number of ordinal categories (analyst-specified)
#' @param us_Serialize  Serialize-encoded current threshold vector
#'   (length J-1) - this is the site's LAST communicated value
#' @export
semiOPBARTLocalThresholdUpdateDS <- function(nfilter = 5, J,
                                              us_Serialize,
                                              state_name = ".semiOPBART_state") {
  s <- get(state_name, envir = parent.frame())
  us_prev <- semiOPBART_fromSerialize(us_Serialize)
  
  # Local-only threshold update - identical math to the real smopbart()
  # but using only this site's own Z, no aggregation across sites.
  new_us <- us_prev
  
  # FIX: guard against length(us_prev) == 1 (binary outcome, J = 2) --
  # see the identical fix/comment on the block-MCMC threshold loops above.
  if (length(us_prev) >= 2) {
    for (j in 2:length(us_prev)) {
      z_j  <- s$Z[s$Y == j]
      z_j1 <- s$Z[s$Y == j + 1]
      
      if (length(z_j) < nfilter || length(z_j1) < nfilter) {
        # Not enough data at this site for this category - keep previous value
        next
      }
      
      if (j == length(us_prev)) {
        lower_bound <- max(max(z_j), new_us[j - 1])
        upper_bound <- min(min(z_j1), Inf)
      } else {
        lower_bound <- max(max(z_j), new_us[j - 1])
        upper_bound <- min(min(z_j1), us_prev[j + 1])
      }
      
      if (is.finite(lower_bound) && is.finite(upper_bound) && 
          upper_bound > lower_bound) {
        new_us[j] <- stats::runif(1, lower_bound, upper_bound)
      }
    }
  }
  new_us[1] <- 0
  
  # Store the updated thresholds locally for future iterations
  s$us_local <- new_us
  assign(state_name, s, envir = parent.frame())
  
  # Return nothing - this is a local update only
  invisible(NULL)
}