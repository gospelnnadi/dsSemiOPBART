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
library(SoftBart)
library(truncnorm)


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

  message(sprintf(
  "[semiOPBARTLocalInitEDS] Initialized with p=%d, J=%d, num_tree=%d, k=%.3g",
  ncol(s$W), nlevels(s$Y), num_tree, k
))

  assign(state_name, list(
    X = s$X, Y = as.numeric(s$Y), W = s$W, dv = s$dv, norm_info = s$norm_info,
    num_tree = num_tree, k = k,
    hypers_owned = hypers_owned, opts = opts,
    owned_forest = owned_forest, frozen_idx = integer(0),
    fx_local = rep(0, length(s$Y)), n = length(s$Y)
  ), envir = parent.frame())

  list(n = length(s$Y), p = ncol(s$W), J = nlevels(s$Y))
}




#' @name semiOPBARTLocalVarCountsDS
#' Variable importance for Architectures D and E: a snapshot of how many
#' times each X feature is currently used as a splitting variable across
#' `s$owned_forest`, via the same `$get_counts()` SoftBart forest method
#' the non-federated smopbart() already uses for this (see semiOPBART.R,
#' `varcounts[i,] = smopbart_forest$get_counts()`) -- just called once,
#' against the FINAL trained forest, rather than accumulated every sweep,
#' to avoid changing the already-tested semiOPBARTLocalMCMC(Improved)()
#' hot loop for this. A per-feature split-count summary is model
#' structure, not row-level data -- the same disclosure footprint as
#' n_swapped above or get_sigma_mu() in the non-federated version.
#'
#' Architecture F needs no equivalent of this: its local smopbart() call
#' already computes var_counts per sweep and it's now returned directly by
#' semiOPBARTLocalFitFDS() (localFitDS.R) as var_counts_mean.
#'
#' @param state_name  which run's state to read counts from (D or E; each
#'   run_tag/architecture combination has its own state_name -- see
#'   semiOPBART_run_names())
#' @export
semiOPBARTLocalVarCountsDS <- function(state_name = ".semiOPBART_state") {
  s <- get(state_name, envir = parent.frame())

  if (!is.function(s$owned_forest$get_counts))
    stop("semiOPBARTLocalVarCountsDS: this site's forest object has no ",
         "get_counts() method -- check the installed SoftBart version.")

  counts <- as.numeric(s$owned_forest$get_counts())
  feat_names <- if (!is.null(colnames(s$X))) colnames(s$X) else NULL
  if (!is.null(feat_names) && length(feat_names) == length(counts))
    names(counts) <- feat_names

  list(n = s$n, var_counts = counts)
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
#' KNOWN DEVIATION FROM THE REFERENCE ALGORITHM (confirmed while building
#' ds.semiOPBARTTrainRef() -- see that function's docstring for the full
#' comparison): within each sweep of the loop below, Z is sampled BEFORE
#' thresholds are updated, using stale, previous-block thresholds instead
#' of freshly-updated ones -- the real smopbart() does this in the
#' opposite order. do_gibbs()'s return value is also discarded here in
#' favor of a separate do_predict() call made BEFORE growth (harmless
#' today only because the next sweep re-derives f(x) fresh rather than
#' trusting the stored value, not because it's a faithful reproduction of
#' `fx_train[i,] = smopbart_forest$do_gibbs(...)`). Neither has been
#' patched here -- this function still behaves exactly as documented
#' above/below -- but if you need behavior verified against the reference
#' algorithm, use ds.semiOPBARTTrainRef() (sync_every = 1, correct order,
#' assembled from semiOPBARTLocalSyncDS()/GrowDS()/ThresholdStatsDS(),
#' which DO already get this right) rather than this function.
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
semiOPBARTLocalMCMC <- function(
    num_local_sweeps,
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
    seed = 35,
    sd = 1
) {

    set.seed(seed)

    caller_env <- parent.frame()

    s <- get(
        state_name,
        envir = caller_env,
        inherits = FALSE
    )

    # Parameters from previous synchronization
    global_theta <- semiOPBART_fromSerialize(theta_Serialize)
    global_us <- semiOPBART_fromSerialize(us_Serialize)

    # Current local state
    local_theta <- global_theta
    local_us <- global_us



    y_idx <- s$Y_idx


    Y_raw <- s$Y

if (is.factor(Y_raw)) {

    # Preserve the existing ordered factor levels
    Y_levels <- levels(Y_raw)
    y_idx <- as.integer(Y_raw)

} else {

    # Explicitly map observed categories to 1,...,J
    Y_levels <- sort(unique(Y_raw))

    if (anyNA(Y_raw)) {
        stop("Y contains NA values")
    }

    y_idx <- match(Y_raw, Y_levels)
}

J <- length(Y_levels)

if (anyNA(y_idx)) {
    stop("Y_idx contains NA after recoding")
}

if (!all(y_idx %in% seq_len(J))) {
    stop(sprintf(
        "Invalid Y_idx: values are %s but expected 1:%d",
        paste(sort(unique(y_idx)), collapse = ", "),
        J
    ))
}

s$Y_idx <- y_idx
s$Y_levels <- Y_levels
s$J <- J

    WtW <- crossprod(s$W)

    #WtZr_accum <- numeric(ncol(s$W))

    th_stats_list <- vector(
        "list",
        num_local_sweeps
    )

    for (t in seq_len(num_local_sweeps)) {

        # ==========================================================
        # 1. UPDATE THRESHOLDS USING PREVIOUS Z
        # ==========================================================

        Z_previous <- s$Z
        previous_us <- local_us

        if (length(local_us) >= 2) {

            for (j in 2:length(local_us)) {

                z_j <- Z_previous[y_idx == j]
                z_j1 <- Z_previous[y_idx == j + 1]

                if (length(z_j) < nfilter_threshold ||
                    length(z_j1) < nfilter_threshold) {
                    next
                }

                lower_bound <- max(
                    z_j,
                    local_us[j - 1]
                )

                if (j == length(local_us)) {

                    upper_bound <- min(z_j1)

                } else {

                    upper_bound <- min(
                        z_j1,
                        previous_us[j + 1]
                    )
                }

                if (is.finite(lower_bound) &&
                    is.finite(upper_bound) &&
                    upper_bound > lower_bound) {

                    local_us[j] <- runif(
                        1,
                        lower_bound,
                        upper_bound
                    )
                }
            }
        }

        local_us[1] <- 0

         # ==========================================================
        # 2. fx_(i-1)
        # ==========================================================

        fx_local <- as.numeric(
            s$owned_forest$do_predict(s$X)
        )

        
        # ==========================================================
        # 3. SAMPLE Z_i
        # ==========================================================

        lower_us <- c(-Inf, local_us)
        upper_us <- c(local_us, Inf)
        lower <- lower_us[y_idx]
        upper <- upper_us[y_idx]


        if (any(!is.finite(lower) &
               lower != -Inf)) {
            stop("Invalid lower truncation bounds")
        }

        if (any(lower >= upper)) {
            stop("Invalid truncation interval for Z sampling")
        }
        
        Z <- truncnorm::rtruncnorm(
            n = s$n,
            a = lower,
            b = upper,
            mean = as.numeric(
                fx_local + s$W %*% local_theta
            ),
            sd = sd #1 
        )

        # ==========================================================
        # 4. THETA_i
        # ==========================================================

        Z_star <- Z - fx_local

        WtZr <- crossprod(
            s$W,
            Z_star
        )

        theta_hat <- solve(
            WtW,
            WtZr
        )

        theta_sigma <- solve(WtW)

        local_theta <- as.numeric(
            mvtnorm::rmvnorm(
                n = 1,
                mean = theta_hat,
                sigma = theta_sigma
            )
        )

        # ==========================================================
        # 5. FOREST UPDATE -> f_i
        # ==========================================================

        Z_tilde <- Z -
            s$W %*% local_theta

        fx_new <- as.numeric(
            s$owned_forest$do_gibbs(
                s$X,
                Z_tilde,
                s$X,
                1
            )
        )

        
        # ==========================================================
        # 6. THRESHOLD STATISTICS
        # ==========================================================

        th_stats <- lapply(seq_len(J), function(j) {

            zj <- Z[y_idx == j]

            if (length(zj) < nfilter_threshold) {

                return(
                    c(
                        lo = NA_real_,
                        hi = NA_real_,
                        n = length(zj)
                    )
                )
            }

            if (identical(threshold_method, "exact")) {

                c(
                    lo = min(zj),
                    hi = max(zj),
                    n = length(zj)
                )

            } else {

                c(
                    lo = unname(
                        quantile(zj, 0.05)
                    ),
                    hi = unname(
                        quantile(zj, 0.95)
                    ),
                    n = length(zj)
                )
            }
        })

        names(th_stats) <-
            as.character(seq_len(J))

        th_stats_list[[t]] <- th_stats

        # ==========================================================
        # 7. SAVE STATE
        # ==========================================================

        s$Z <- Z
        s$fx_local <- fx_new
        s$theta <- local_theta
        s$us <- local_us

        assign(
            state_name,
            s,
            envir = caller_env
        )


        # sufficient statistic
        #WtZr_accum <- WtZr_accum + WtZr
    }
 

    final_th_stats <-
        th_stats_list[[num_local_sweeps]]

    list(
        WtW = WtW,
        WtZr =  WtZr, 
        th_stats = final_th_stats,
        local_theta = local_theta,
        local_us = local_us,
        n = s$n,
        num_sweeps = num_local_sweeps
    )
}


# semiopbartDS.R - Add this function

#' @name semiOPBARTLocalCheckExistsDS
#' Lightweight, disclosure-safe check: does a given object name exist
#' in the server environment? Returns TRUE/FALSE only, no data.
#' 
#' This is the single shared entry point every trainer should call
#' before attempting semiOPBARTLocalInitDS/EDS/etc., so auto-demoted
#' sites (which have D_holdout but no D_train) are silently excluded
#' rather than crashing the whole training call.
#' 
#' @param object.name  name of the object to check (e.g. "D_train")
#' @export
semiOPBARTLocalCheckExistsDS <- function(object.name) {
  exists(object.name, envir = parent.frame(), inherits = FALSE)
}