# dsSemiOPBARTTrainTest.R
# ---------------------------------------------------------------------------
# TRAIN (ds.semiOPBARTTrain(), an alias defined in dsSemiOPBARTClient.R) and
# TEST/PREDICT (this file) are separate functions, server- and client-side.
#
# UPDATED for the swap-based redesign (semiopbartDS.R/ds_semiopbart.R):
# there is no longer a single pooled "global" tree ensemble to reconstruct
# and ship to a predicting site -- every site keeps its OWN complete,
# permanent num_tree-tree forest, continuously enriched by tree swaps
# during training. Prediction at an ALREADY-TRAINED site is therefore now
# just "read this site's own persisted .semiOPBART_state$owned_forest and
# call do_predict()" -- no final_trees_Serialize/final_ownership reassembly at
# all, which also removes an entire class of payload (previously: every
# tree in the whole federation, shipped to every predicting site).
#
# CONSEQUENCE, stated plainly: there is now no never-trained-site
# prediction path here (there wasn't a clean one before either, strictly
# speaking, but the pooled tree set at least made it POSSIBLE to points a
# site at trees it never grew). If you need to score a genuinely
# never-trained site, that needs its own design decision (e.g. picking
# one training site's forest as a reference, á la the reference-site
# mechanism explored earlier for Option C) -- not attempted here.
# ---------------------------------------------------------------------------

#source("serializeDS.R")

# =====================  TEST / PREDICT  =====================================

# ---- server side (Opal) ----------------------------------------------------



#' Predict AT THIS SITE ONLY, using this site's OWN cached local forest
#' (.semiOPBART_state$owned_forest -- never any other site's, and never a
#' shadow forest, which is always NULL under Architecture E) combined
#' with the FEDERATED theta/us shipped in from ds.semiOPBARTTrainE()'s
#' output -- not any local re-derivation. Requires
#' semiOPBARTLocalInitEDS()/ds.semiOPBARTTrainE() to have already run
#' here; errors otherwise (Architecture E has no never-trained-site path).
#'
#' Same convention as semiOPBARTLocalPredictDS() above: data.name_test is
#' a PREPARED list(X, W, Y, ...) object (from ds.semiOPBARTPrepare()), not
#' a raw data.frame, so no outcome_col/dv/ecdf re-derivation is needed
#' here -- s_test$Y is already a properly-leveled factor and s_test$X/$W
#' are already normalized/dummy-encoded exactly the same way training's
#' s$X/$W were.
#'
#' @param theta_Serialize, us_Serialize  Serialize-encoded EITHER a plain vector (point
#'   estimate: ds.semiOPBARTPredictE(prediction_method = "point"), the
#'   default) OR a matrix, one row per retained draw (ds.semiOPBARTPredictE
#'   (prediction_method = "draws")) -- is.matrix() after decoding tells
#'   these apart; no separate argument or server function needed for
#'   either choice.
#' @export
semiOPBARTLocalPredictEDS <- function(theta_Serialize, us_Serialize,
                                      data.name_test,
                                      newobj_pred = "semiOPBART_pred_E",
                                      nfilter = 5, state_name = ".semiOPBART_state", seed = 35) {
  set.seed(seed)
  if (!exists(state_name, envir = parent.frame()))
    stop("semiOPBARTLocalPredictEDS(): no cached local state ('", state_name,
         "') at this site -- run ds.semiOPBARTTrainE() here first, with the ",
         "SAME state_name this predict call is using. Architecture E has no ",
         "never-trained-site path -- every site that wants a prediction ",
         "must have trained locally.")
  s_train <- get(state_name, envir = parent.frame())
  if (is.null(s_train$owned_forest))
    stop("this site's local forest is missing -- did ds.semiOPBARTTrainE() ",
         "actually run here (not just semiOPBARTLocalInitEDS())?")

  s_test <- get(data.name_test, envir = parent.frame())
  if (length(s_test$Y) < nfilter) stop("test n below disclosure threshold at this site")

  theta <- semiOPBART_fromSerialize(theta_Serialize)
  us    <- semiOPBART_fromSerialize(us_Serialize)

  fx_new <- as.numeric(s_train$owned_forest$do_predict(s_test$X))
  levels_used <- levels(s_test$Y)
  J <- length(levels_used)
  n_test <- nrow(s_test$X)

  if (is.matrix(theta)) {
    # full posterior draws -- average class probabilities over every
    # retained draw, same as semiOPBARTLocalPredictDS() (Architecture D)
    # already does. fx_new stays a single point (this site's forest is
    # ONE object, not one-per-draw -- see ds.semiOPBARTTrainE()'s
    # docstring), only theta/us vary draw to draw.
    hw_draws <- s_test$W %*% t(theta)
    n_draws  <- nrow(theta)
    combined_draws <- matrix(fx_new, n_test, n_draws) + hw_draws

    prob_list <- vector("list", J)
    for (j in seq_len(J)) {
      prob_list[[j]] <- if (j == 1) pnorm(us[, j] - t(combined_draws))
        else if (j == J) 1 - pnorm(us[, j - 1] - t(combined_draws))
        else pnorm(us[, j] - t(combined_draws)) - pnorm(us[, j - 1] - t(combined_draws))
    }
    prob_matrix <- matrix(vapply(prob_list, colMeans, numeric(n_test)), n_test, J)
  } else {
    hw_new <- as.numeric(s_test$W %*% theta)
    combined <- fx_new + hw_new
    prob_matrix <- matrix(NA_real_, n_test, J)
    for (j in seq_len(J)) {
      prob_matrix[, j] <- if (j == 1) pnorm(us[j] - combined)
        else if (j == J) 1 - pnorm(us[j - 1] - combined)
        else pnorm(us[j] - combined) - pnorm(us[j - 1] - combined)
    }
  }
  colnames(prob_matrix) <- paste0("P_", levels_used)
  map_score <- levels_used[apply(prob_matrix, 1, which.max)]

  base::assign(newobj_pred, list(prob = prob_matrix, map = map_score, truth = s_test$Y),
         envir = parent.frame())
  list(n_test = n_test)
}



#' @name semiOPBARTLocalPredictFDS
#' Predict AT THIS SITE ONLY, using this site's OWN cached local forest
#' (.semiOPBART_local_fit_F$owned_forest -- never any other site's) combined
#' with Architecture F's FEDERATED theta/us (from
#' ds.semiOPBARTCombineF()) -- point estimate only, see
#' ds.semiOPBARTPredictF()'s docstring for why. Requires
#' semiOPBARTLocalFitFDS() to have already run here; errors otherwise
#' (Architecture F has no never-trained-site path either).
#'
#' Unlike semiOPBARTLocalPredictEDS(), data.name_test here is a plain,
#' still-raw data.frame -- Architecture F's Fit step never routed through
#' dsSemiOPBARTDataPrep.R's normalization pipeline (it calls smopbart()
#' directly on a raw data.frame, same as the original design), so its
#' own dv/ecdfs (not any shared normalization) are what a new row must
#' go through here.
#' @export
#' @name semiOPBARTLocalPredictFDS
#' Predict AT THIS SITE ONLY, using this site's OWN cached local forest
#' (.semiOPBART_local_fit_F$owned_forest -- never any other site's) combined
#' with Architecture F's FEDERATED theta/us (from
#' ds.semiOPBARTCombineF()) -- point estimate only, see
#' ds.semiOPBARTPredictF()'s docstring for why. Requires
#' semiOPBARTLocalFitFDS() to have already run here; errors otherwise
#' (Architecture F has no never-trained-site path either).
#'
#' data.name_test is a PREPARED (list(X, W, Y, ...)) object here, same
#' convention as semiOPBARTLocalPredictDS()/PredictEDS() -- Architecture
#' F's Fit step now goes through smopbart()'s "accept a prepared object"
#' branch (semiOPBART.R), so df_test$X/$W are ALREADY dummy-encoded and
#' normalized under whichever method this run used; re-running them
#' through fit$dv/fit$ecdfs here would double-process them (and, for
#' federated_minmax/federated_ecdf, fit$ecdfs doesn't even exist -- see
#' semiOPBART.R's matching fix). Use df_test$X/$W directly, exactly like
#' semiOPBARTLocalPredictDS()/PredictEDS() already do.
#' @export
semiOPBARTLocalPredictFDS <- function(theta_Serialize, us_Serialize,
                                      data.name_test, outcome_col,
                                      newobj_pred = "semiOPBART_pred_F",
                                      nfilter = 5, state_name = ".semiOPBART_local_fit_F", seed = 35) {
  set.seed(seed)
  if (!exists(state_name, envir = parent.frame()))
    stop("semiOPBARTLocalPredictFDS(): no cached local fit ('", state_name,
         "') at this site -- run ds.semiOPBARTFitF() here first, with the ",
         "SAME state_name this predict call is using. Architecture F has no ",
         "never-trained-site path -- every site that wants a prediction ",
         "must fit locally.")
  fit <- get(state_name, envir = parent.frame())

  df_test <- eval(parse(text = data.name_test), envir = parent.frame())
  if (nrow(df_test$X) < nfilter)
    stop("test n below disclosure threshold at this site")

  X_new <- as.matrix(df_test$X)
  W_new <- as.matrix(df_test$W)
  levels_used <- levels(df_test$Y)

  theta <- semiOPBART_fromSerialize(theta_Serialize)
  us    <- semiOPBART_fromSerialize(us_Serialize)

  fx_new <- as.numeric(fit$owned_forest$do_predict(X_new))
  hw_new <- as.numeric(W_new %*% theta)     # FEDERATED theta, not fit$theta_p
  combined <- fx_new + hw_new

  J <- length(levels_used)
  prob_matrix <- matrix(NA_real_, length(combined), J)
  for (j in seq_len(J)) {
    prob_matrix[, j] <- if (j == 1) pnorm(us[j] - combined)
    else if (j == J) 1 - pnorm(us[j - 1] - combined)
    else pnorm(us[j] - combined) - pnorm(us[j - 1] - combined)
  }
  colnames(prob_matrix) <- paste0("P_", levels_used)
  map_score <- levels_used[apply(prob_matrix, 1, which.max)]

  base::assign(newobj_pred, list(prob = prob_matrix, map = map_score,
                           truth = df_test$Y),
         envir = parent.frame())
  list(n_test = length(combined))
}


