# # dsSemiOPBARTLocalFitF.R (server side)
# # ---------------------------------------------------------------------------
# # Architecture F's fit -- see ds_localFit.R's header. formula_str/
# # linear_formula_str arrive Serialize-encoded (they contain `~`/`+`/spaces,
# # which break Opal's call() reconstruction unwrapped); data.name stays a
# # plain, unwrapped string, same convention as every other server function
# # in this codebase.
# # ---------------------------------------------------------------------------

# #' @name semiOPBARTLocalFitFDS
# #' Fit the unmodified single-site smopbart() locally and return ONLY
# #' posterior summaries -- never raw MCMC draws, never row-level predictions.
# #'
# #' @param data.name  site-side raw data.frame name
# #' @export
# semiOPBARTLocalFitFDS <- function(formula_str, linear_formula_str, data.name,
#                                   num_tree = 20, k = 1,
#                                   num_burn = 1000, num_save = 1000,
#                                   nfilter = 5, state_name = ".semiOPBART_local_fit_F", seed = 35) {
#   set.seed(seed)
#   formula_str        <- semiOPBART_fromSerialize(formula_str)
#   linear_formula_str <- semiOPBART_fromSerialize(linear_formula_str)

#   df <- eval(parse(text = data.name), envir = parent.frame())
#   if (nrow(df$X) < nfilter) stop("site n below disclosure threshold")

#   formula <- stats::as.formula(formula_str)
#   linear_formula <- stats::as.formula(linear_formula_str)
#   print(paste0("semiOPBARTLocalFitFDS: colnames(df$X) = ", paste(colnames(df$X), collapse = ", "), ", colnames(df$W) = ", paste(colnames(df$W), collapse = ", "), ", length(df$Y) = ", length(df$Y)))

#   fit <- smopbart(formula = formula, linear_formula = linear_formula,
#                   train_df = df, test_df = df,
#                   num_tree = num_tree, k = k,
#                   opts = { o <- SoftBart::Opts(); o$num_burn <- num_burn;
#                   o$num_save <- num_save; o },
#                   verbose = FALSE, seed=seed)
#   print("\nsemiOPBARTLocalFitFDS: fit complete, returning posterior summaries only\n")
#   fit$feat_names_x <- attr(fit$dv$terms, "term.labels")
#   fit$feat_names_w <- all.vars(fit$linear_formula)
#   fit$num_tree_used <- num_tree
#   fit$k_used <- k
#   # state_name defaults to the old fixed slot, but the client now always
#   # passes a run-tagged name -- keeps F's cached local fit from colliding
#   # with ANOTHER run (different pathology/model) at the same site, not
#   # just with D/E's own state (which lives under a completely different
#   # slot regardless)
#   assign(state_name, fit, envir = parent.frame())

# #   list(n          = nrow(df$X),
# #        theta_mean = colMeans(fit$theta_p),
# #        theta_cov  = stats::cov(fit$theta_p),   # for inverse-variance combine
# #        us_mean    = colMeans(fit$us_p))
# # var_counts (from smopbart_forest$get_counts() every post-burn-in sweep,
#   # column-named by X feature -- see semiOPBART.R) was already being
#   # computed here and simply never returned to the client. It's a
#   # per-feature split-count SUMMARY (mean across the kept posterior
#   # sweeps), not row-level data, so it's safe to return like theta_mean/
#   # us_mean already are -- this is what powers the explainability module's
#   # Architecture F variable-importance plot (dsSemiOPBARTExplain.R).
#   var_counts_mean <- colMeans(fit$var_counts)

#   list(n          = nrow(df$X),
#        theta_mean = colMeans(fit$theta_p),
#        theta_cov  = stats::cov(fit$theta_p),   # for inverse-variance combine
#        us_mean    = colMeans(fit$us_p),
#        var_counts_mean = var_counts_mean)
# }


# dsSemiOPBARTLocalFitF.R (server side)
# ---------------------------------------------------------------------------
# Architecture F's fit -- see ds_localFit.R's header. formula_str/
# linear_formula_str arrive Serialize-encoded (they contain `~`/`+`/spaces,
# which break Opal's call() reconstruction unwrapped); data.name stays a
# plain, unwrapped string, same convention as every other server function
# in this codebase.
# ---------------------------------------------------------------------------

#' @name semiOPBARTLocalFitFDS
#' Fit the unmodified single-site smopbart() locally and return ONLY
#' posterior summaries -- never raw MCMC draws, never row-level predictions.
#'
#' @param data.name  site-side raw data.frame name
#' @export
semiOPBARTLocalFitFDS <- function(formula_str, linear_formula_str, data.name,
                                  num_tree = 20, k = 1,
                                  num_burn = 1000, num_save = 1000,
                                  nfilter = 5, state_name = ".semiOPBART_local_fit_F", seed = 35, sd = 1) {
  set.seed(seed)
  formula_str        <- semiOPBART_fromSerialize(formula_str)
  linear_formula_str <- semiOPBART_fromSerialize(linear_formula_str)

  df <- eval(parse(text = data.name), envir = parent.frame())
  if (nrow(df$X) < nfilter) stop("site n below disclosure threshold")

  formula <- stats::as.formula(formula_str)
  linear_formula <- stats::as.formula(linear_formula_str)

  fit <- smopbart(formula = formula, linear_formula = linear_formula,
                  train_df = df, test_df = df,
                  num_tree = num_tree, k = k,
                  opts = { o <- SoftBart::Opts(); o$num_burn <- num_burn;
                  o$num_save <- num_save; o },
                  verbose = FALSE, seed=seed, sd = sd)
  fit$feat_names_x <- attr(fit$dv$terms, "term.labels")
  fit$feat_names_w <- all.vars(fit$linear_formula)
  fit$num_tree_used <- num_tree
  fit$k_used <- k
  # state_name defaults to the old fixed slot, but the client now always
  # passes a run-tagged name -- keeps F's cached local fit from colliding
  # with ANOTHER run (different pathology/model) at the same site, not
  # just with D/E's own state (which lives under a completely different
  # slot regardless)
  assign(state_name, fit, envir = parent.frame())

  # var_counts (from smopbart_forest$get_counts() every post-burn-in sweep,
  # column-named by X feature -- see semiOPBART.R) was already being
  # computed here and simply never returned to the client. It's a
  # per-feature split-count SUMMARY (mean across the kept posterior
  # sweeps), not row-level data, so it's safe to return like theta_mean/
  # us_mean already are -- this is what powers the explainability module's
  # Architecture F variable-importance plot (dsSemiOPBARTExplain.R).
  var_counts_mean <- colMeans(fit$var_counts)

  # us_cov: theta_cov's counterpart for the thresholds, computed the same
  # way (stats::cov() of the posterior draws) -- previously never
  # computed at all, meaning us_mean has only ever been combined by
  # simple n-weighted averaging across sites (see ds.semiOPBARTCombineF()),
  # never inverse-variance or random-effects, purely because there was no
  # variance to weight by. fit$us_p[,1] is EXCLUDED here: it's always
  # fixed at exactly 0 by construction (the ordinal probit's
  # identifiability constraint, see semiOPBART.R's `us_p[,1] = 0`), so
  # its variance is exactly 0 -- including it would make us_cov singular.
  us_cov <- if (ncol(fit$us_p) >= 2) stats::cov(fit$us_p[, -1, drop = FALSE]) else NULL

  list(n          = nrow(df$X),
       theta_mean = colMeans(fit$theta_p),
       theta_cov  = stats::cov(fit$theta_p),   # for inverse-variance combine
       us_mean    = colMeans(fit$us_p),
       us_cov     = us_cov,                    # for inverse-variance/random-effects
                                                 # combine of the FREE thresholds only
       var_counts_mean = var_counts_mean,
       site_n     = nrow(df$X),
        # sufficient statistics for federated combination
       WtW = fit$WtW,
       WtZr = fit$WtZr,

       # threshold statistics for federated threshold update
       th_stats = fit$th_stats)
}