# dsSemiOPBARTExternalPredict.R
# ---------------------------------------------------------------------------
# OPT-IN prediction at sites that never trained, Architecture D.
#
# UPDATED for the swap-based redesign: no pooled "global" tree ensemble
# exists anymore -- every site keeps its own complete, permanent forest.
# semiOPBARTLocalExportForestDS() below exports exactly ONE site's own
# current forest, at the orchestrator/user's explicit choice of which
# site (ds.semiOPBARTExportReferenceForest(), ds_externalPredict.R).
# semiOPBARTLocalPredictExternalDS() itself needed NO change: it already
# builds a fresh eval_forest from whatever final_trees_Serialize it's handed,
# regardless of whether that Serialize came from a pooled reconstruction (the
# old design) or one site's own forest (now) -- it was already source-
# agnostic.
#
# ds.semiOPBARTPredict() (dsSemiOPBARTTrainTest.R) is UNCHANGED: it still
# requires .semiOPBART_state (i.e. the site ran ds.semiOPBARTTrain()
# itself) and fails otherwise. That remains the default, lower-disclosure
# path for sites that trained.
#
# PRECONDITION: this only works if normalize_method was "federated_minmax"
# or "federated_ecdf" for the original ds.semiOPBARTPrepare() call.
# "local_ecdf" normalization is inherently site-specific -- there is no
# single global reference to ship, so training itself already refuses to
# run under "local_ecdf" (semiOPBARTLocalInitDS()'s own hard requirement).
# ---------------------------------------------------------------------------

#source("serializeDS.R")

# ---- server side (Opal) ----------------------------------------------------

#' Export THIS site's own current forest (all num_tree trees), for the
#' client to ship elsewhere as a reference-site model. Requires
#' ds.semiOPBARTTrain() to have already run here; errors otherwise.
#' @export
semiOPBARTLocalExportForestDS <- function(
    state_name = ".semiOPBART_state"
) {

  caller_env <- parent.frame()

  if (!exists(
      state_name,
      envir = caller_env,
      inherits = FALSE
  )) {
    stop(
      "semiOPBARTLocalExportForestDS(): no cached local state ('",
      state_name,
      "') at this site -- run training here first."
    )
  }

  s <- get(
    state_name,
    envir = caller_env,
    inherits = FALSE
  )

  if (is.null(s$owned_forest)) {
    stop(
      "this site's local forest is missing -- did training actually run here?"
    )
  }

  # ---------------------------------------------------------------
  # Normalization reference
  # ---------------------------------------------------------------
  #
  # Expected location:
  #   s$norm_info$local_ecdf
  #   s$norm_info$global_ecdf
  #
  # We export whichever reference actually exists.
  #
  norm_info <- s$norm_info

local_ecdf <- NULL
global_ecdf <- NULL

if (!is.null(norm_info)) {
  if (identical(norm_info$method, "local_ecdf") && !is.null(norm_info$ecdfs)) {
    local_ecdf <- norm_info$ecdfs
  }
  if (identical(norm_info$method, "federated_ecdf") && !is.null(norm_info$global_ecdf)) {
    global_ecdf <- norm_info$global_ecdf
  }
}

  if (is.null(local_ecdf) && is.null(global_ecdf)) {
    stop(
      "semiOPBARTLocalExportForestDS(): neither local_ecdf nor ",
      "global_ecdf is available in this site's normalization state. ",
      "Cannot export an external prediction reference."
    )
  }

  # ---------------------------------------------------------------
  # Forest
  # ---------------------------------------------------------------

  trees_Serialize <- semiOPBART_treesToSerialize(s$owned_forest)

  list(
    trees_Serialize = trees_Serialize,
    num_tree = s$num_tree,

    # Ship normalization alongside the forest.
    local_ecdf_Serialize =
      if (!is.null(local_ecdf))
        semiOPBART_toSerialize(local_ecdf)
      else
        NULL,

    global_ecdf_Serialize =
      if (!is.null(global_ecdf))
        semiOPBART_toSerialize(global_ecdf)
      else
        NULL,

    has_local_ecdf = !is.null(local_ecdf),
    has_global_ecdf = !is.null(global_ecdf)
  )
}

#' Predict at a site with NO .semiOPBART_state (never ran Init/training).
#' Builds a fresh hypers/opts template locally from whatever X this site
#' already has (either already prepared, or just re-prepared by
#' ds.semiOPBARTPredictExternal() below using the shipped global
#' normalization reference) -- unlike semiOPBARTLocalPredictDS(), this
#' does NOT read .semiOPBART_state at all, which is precisely what makes
#' it usable at a site that never trained.
#'
#' @param theta_draws_Serialize, us_draws_Serialize, final_trees_Serialize  same as
#'   semiOPBARTLocalPredictDS() -- final_trees_Serialize here is ONE reference
#'   site's own forest (from semiOPBARTLocalExportForestDS() above), not
#'   a pooled reconstruction; this function doesn't need to know or care
#'   which, it just loads whatever tree Serialize it's handed
#' @param data.name_test  a PREPARED (list(X, W, Y, ...)) object at this
#'   site, on the SAME normalized scale as the trees were grown on --
#'   i.e. built using the same normalization reference, whether that
#'   happened via the original shared ds.semiOPBARTPrepare() call or via
#'   ds.semiOPBARTPredictExternal()'s own re-prepare step
#' @param num_tree, k  must match what training used, for sigma_mu scaling
#'   to be consistent with the shipped tree ensemble
#' @export
semiOPBARTLocalPredictExternalDS <- function(
    theta_draws_Serialize,
    us_draws_Serialize,
    final_trees_Serialize,
    ecdf_reference_Serialize,
    data.name_test,
    num_tree,
    k = 1,
    newobj_pred = "semiOPBART_pred",
    levels_used = 0:4,
    nfilter = 5,
    seed = 35
) {

  set.seed(seed)

  caller_env <- parent.frame()

  s_test <- get(
    data.name_test,
    envir = caller_env,
    inherits = FALSE
  )

  if (length(s_test$Y) < nfilter) {
    stop("test n below disclosure threshold")
  }

  theta_draws <-
    semiOPBART_fromSerialize(theta_draws_Serialize)

  us_draws <-
    semiOPBART_fromSerialize(us_draws_Serialize)

  ecdf_reference <-
    semiOPBART_fromSerialize(ecdf_reference_Serialize)

  n_final_trees <-
    length(
      Serializelite::fromSerialize(
        final_trees_Serialize,
        simplifyVector = FALSE
      )
    )

  # ---------------------------------------------------------------
  # Apply exported ECDF to X
  # ---------------------------------------------------------------

  if (length(ecdf_reference) != ncol(s_test$X)) {
    stop(
      "ECDF reference has ",
      length(ecdf_reference),
      " columns, but test X has ",
      ncol(s_test$X)
    )
  }

  for (j in seq_len(ncol(s_test$X))) {
    s_test$X[, j] <-
      ecdf_reference[[j]](
        s_test$X[, j]
      )
  }

  # ---------------------------------------------------------------
  # Reconstruct reference forest
  # ---------------------------------------------------------------

  eval_forest <-
    semiOPBART_treesFromSerialize(
      final_trees_Serialize
    )

  fx_new <-
    as.numeric(
      eval_forest$do_predict(s_test$X)
    )

  # ---------------------------------------------------------------
  # Posterior prediction
  # ---------------------------------------------------------------

  hw_draws <-
    s_test$W %*% t(theta_draws)

  n_draws <-
    nrow(theta_draws)

  combined_draws <-
    matrix(
      fx_new,
      nrow = nrow(s_test$X),
      ncol = n_draws
    ) +
    hw_draws

  J <- length(levels_used)

  if (length(us_draws) == 0L) {
    stop("us posterior draws are empty")
  }

  prob_list <- vector("list", J)

  for (j in seq_len(J)) {

    if (j == 1L) {

      prob_list[[j]] <-
        pnorm(
          us_draws[, j] -
            t(combined_draws)
        )

    } else if (j == J) {

      prob_list[[j]] <-
        1 -
        pnorm(
          us_draws[, j - 1L] -
            t(combined_draws)
        )

    } else {

      prob_list[[j]] <-
        pnorm(
          us_draws[, j] -
            t(combined_draws)
        ) -
        pnorm(
          us_draws[, j - 1L] -
            t(combined_draws)
        )
    }
  }

  n_test <-
    nrow(s_test$X)

  prob_matrix <-
    matrix(
      vapply(
        prob_list,
        colMeans,
        numeric(n_test)
      ),
      nrow = n_test,
      ncol = J
    )

  colnames(prob_matrix) <-
    paste0(
      "P_",
      levels_used
    )

  map_score <-
    levels_used[
      apply(
        prob_matrix,
        1,
        which.max
      )
    ]

  assign(
    newobj_pred,
    list(
      prob = prob_matrix,
      map = map_score,
      truth = s_test$Y,
      used_external_shipping = TRUE
    ),
    envir = caller_env
  )

  list(
    n_test = n_test
  )
}


# semiOPBARTLocalExportForestDS <- function(state_name = ".semiOPBART_state") {
#   if (!exists(state_name, envir = parent.frame()))
#     stop("semiOPBARTLocalExportForestDS(): no cached local state ('", state_name,
#          "') at this site -- run ds.semiOPBARTTrain() here first (this site ",
#          "must be one that actually trained).")
#   s <- get(state_name, envir = parent.frame())
#   if (is.null(s$owned_forest))
#     stop("this site's local forest is missing -- did ds.semiOPBARTTrain() ",
#          "actually run here?")

#   list(trees_Serialize = semiOPBART_treesToSerialize(s$owned_forest, 0:(s$num_tree - 1)),
#        num_tree = s$num_tree)
# }

# #' Predict at a site with NO .semiOPBART_state (never ran Init/training).
# #' Builds a fresh hypers/opts template locally from whatever X this site
# #' already has (either already prepared, or just re-prepared by
# #' ds.semiOPBARTPredictExternal() below using the shipped global
# #' normalization reference) -- unlike semiOPBARTLocalPredictDS(), this
# #' does NOT read .semiOPBART_state at all, which is precisely what makes
# #' it usable at a site that never trained.
# #'
# #' @param theta_draws_Serialize, us_draws_Serialize, final_trees_Serialize  same as
# #'   semiOPBARTLocalPredictDS() -- final_trees_Serialize here is ONE reference
# #'   site's own forest (from semiOPBARTLocalExportForestDS() above), not
# #'   a pooled reconstruction; this function doesn't need to know or care
# #'   which, it just loads whatever tree Serialize it's handed
# #' @param data.name_test  a PREPARED (list(X, W, Y, ...)) object at this
# #'   site, on the SAME normalized scale as the trees were grown on --
# #'   i.e. built using the same normalization reference, whether that
# #'   happened via the original shared ds.semiOPBARTPrepare() call or via
# #'   ds.semiOPBARTPredictExternal()'s own re-prepare step
# #' @param num_tree, k  must match what training used, for sigma_mu scaling
# #'   to be consistent with the shipped tree ensemble
# #' @export
# semiOPBARTLocalPredictExternalDS <- function(theta_draws_Serialize, us_draws_Serialize,
#                                               final_trees_Serialize,ecdf_reference_Serialize, data.name_test,
#                                               num_tree, k = 1, 
#                                               newobj_pred = "semiOPBART_pred",
#                                               levels_used = 0:4,
#                                               nfilter = 5, seed = 35) {
  
#   set.seed(seed)
#   s_test <- get(data.name_test, envir = parent.frame())

#   if (length(s_test$Y) < nfilter) stop("test n below disclosure threshold")
  

#   theta_draws <- semiOPBART_fromSerialize(theta_draws_Serialize)
#   us_draws    <- semiOPBART_fromSerialize(us_draws_Serialize)
#   n_final_trees <- length(Serializelite::fromSerialize(final_trees_Serialize, simplifyVector = FALSE))
#   ecdf_reference <- semiOPBART_fromSerialize(ecdf_reference_Serialize)
#   apply_norm <- function(X) { for (j in seq_len(ncol(X))) X[, j] <- ecdf_reference[[j]](X[, j]); X }
#   s_test$X <- apply_norm(s_test$X)

#   # # built FRESH from local X/Y -- no .semiOPBART_state dependency, which is
#   # # the entire point of this function existing separately from
#   # # semiOPBARTLocalPredictDS()
#   # hypers_template <- SoftBart::Hypers(X = s_test$X, Y = rep(1, nrow(s_test$X)))
#   # hypers_template$sigma_mu <- 3 / k / sqrt(num_tree)
#   # hypers_template$sigma     <- 1
#   # hypers_template$sigma_hat <- 1
#   # hypers_template$num_tree  <- n_final_trees
#   # # same confirmed fix as semiOPBARTLocalInitDS() (semiopbartDS.R):
#   # # without this, multi-level categorical x_features would be mishandled
#   # # by the Dirichlet variable-selection prior
#   # hypers_template$group <- dummy_assign(s_test$dv)
#   # opts <- SoftBart::Opts(); opts$update_sigma <- FALSE

#   # eval_forest <- SoftBart::MakeForest(hypers_template, opts, FALSE)
#   # semiOPBART_treesFromSerialize(eval_forest, final_trees_Serialize, seq_len(n_final_trees) - 1)
  
#   eval_forest <- semiOPBART_treesFromSerialize(final_trees_Serialize)

#   fx_new <- as.numeric(eval_forest$do_predict(s_test$X))
#   hw_draws <- s_test$W %*% t(theta_draws)
#   n_draws  <- nrow(theta_draws)
#   combined_draws <- matrix(fx_new, nrow(s_test$X), n_draws) + hw_draws

  
#   J <- length(levels_used)
#   prob_list <- vector("list", J)
#   for (j in seq_len(J)) {
#     if (j == 1) prob_list[[j]] <- pnorm(us_draws[, j] - t(combined_draws))
#     else if (j == J) prob_list[[j]] <- 1 - pnorm(us_draws[, j - 1] - t(combined_draws))
#     else prob_list[[j]] <- pnorm(us_draws[, j] - t(combined_draws)) -
#                             pnorm(us_draws[, j - 1] - t(combined_draws))
#   }
#   n_test <- nrow(s_test$X)
#   prob_matrix <- matrix(vapply(prob_list, colMeans, numeric(n_test)), n_test, J)
#   colnames(prob_matrix) <- paste0("P_", levels_used)
#   map_score <- levels_used[apply(prob_matrix, 1, which.max)]

#   assign(newobj_pred, list(prob = prob_matrix, map = map_score, truth = s_test$Y,
#                             used_external_shipping = TRUE),
#          envir = parent.frame())
#   list(n_test = n_test)
# }



