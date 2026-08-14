# dsSemiOPBARTLocalCombineF.R (server side)
# ---------------------------------------------------------------------------
# Architecture F needs no server-side combine step (see ds_localCombine.R
# -- pure client arithmetic). Only a shared helper lives here, used by
# semiOPBARTLocalPredictFDS() (trainTestDS.R) to rebuild a design matrix
# from F's still-raw test data.frame the SAME way dsSemiOPBARTDataPrep.R's
# Transform step does (model.matrix + drop intercept), not by raw column
# selection -- raw selection is wrong whenever a w_feature is categorical.
# ---------------------------------------------------------------------------

.semiOPBART_buildW <- function(df, linear_formula) {
  W <- model.matrix(linear_formula, df)
  keep <- which(colnames(W) != "(Intercept)")
  as.matrix(W[, keep, drop = FALSE])
}
