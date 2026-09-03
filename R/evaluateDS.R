# dsSemiOPBARTEvaluate.R
# ---------------------------------------------------------------------------
# Federated evaluation, mirroring evaluate_semiopbart()'s full metric set
# (accuracy, balanced accuracy, macro F1, kappa, precision, recall, plus the
# binary 0-vs-non-0 versions and both confusion matrices).
#
# KEY DESIGN POINT: every one of those metrics is a function of the
# confusion matrix alone -- and confusion-matrix cells are just counts,
# which sum safely across sites (subject to the usual nfilter floor). So
# unlike training, evaluation does NOT need an iterative round-trip
# protocol at all: one local confusion-matrix computation per site, one
# aggregation, done. This file is intentionally generic -- it works
# identically whether `newobj_pred` was produced by
# dsSemiOPBARTTrainTest.R (Architecture D) or by the predict addendum in
# dsSemiOPBARTLocalCombine.R (Option C), which is what makes the two
# directly comparable in the orchestrator.
# ---------------------------------------------------------------------------

# ---- server side (Opal) ----------------------------------------------------

#' Build this site's local multiclass and binary (0 vs non-0) confusion
#' matrices from a locally stored prediction object, with nfilter cell
#' suppression.
#'
#' @param pred_obj  name of the local object created by a predict function
#'   (must contain $map, and, absent a separate validation_col, $truth --
#'   both factors with the SAME levels)
#' @param levels_Serialize  Serialize-encoded full ordinal level set, analyst-
#'   specified (as in semiOPBARTLocalPrepareDS()), so an empty local cell
#'   is still a real zero and not a missing level. Previously passed as a
#'   raw vector call() argument -- fixed for consistency with the
#'   Serialize-only wire format used everywhere else in this codebase.
#' @param nfilter   cells with count < nfilter are reported as NA rather
#'   than the true count (still summed downstream as 0 contribution -- see
#'   note in ds.semiOPBARTEvaluate())
#' @param validation_col   "null" (default, as a literal string sentinel --
#'   see below): validate against `pred_obj$truth`, i.e. the same label
#'   the model was trained on -- unchanged behaviour. A column name (e.g.
#'   "USPDST"): validate against THAT column instead, looked up in
#'   `validation_data.name`. If the column doesn't exist at this site, or
#'   exists but has fewer than `nfilter` non-missing values overlapping
#'   with predicted rows, this function returns NULL rather than erroring
#'   -- signalling "no validation possible at this site", which
#'   ds.semiOPBARTEvaluate() aggregates around gracefully rather than
#'   failing the whole run over one site's missing column.
#' @param validation_data.name  "null" (default) or the name of the
#'   site-side object to look `validation_col` up in -- ignored if
#'   validation_col is "null". Both this and validation_col use the
#'   string "null" as a sentinel (matching the pattern already used for
#'   owned_tree_ids_Serialize elsewhere) rather than a bare R NULL, since NULL
#'   is not itself a string and this codebase's wire format requires one.
#' @param dichotomize_threshold  "null" (default): binary evaluation uses
#'   the ORIGINAL baseline-category-vs-rest split (`levels[1]` vs
#'   everything else). A numeric value, as a STRING (e.g. "2"): binary
#'   evaluation instead dichotomizes at THIS threshold -- truth/prediction
#'   threshold becomes "1", below becomes "0". Applies to whichever >=
#'   truth source this call is using (the training label, or a
#'   validation_col) -- pass a DIFFERENT dichotomize_threshold in separate
#'   calls if the label and the validation column need different cut
#'   points (this is how ds.semiOPBARTEvaluate() exposes it: as two
#'   independent arguments, since label-based and validation-based
#'   evaluation are already two separate calls in this codebase's design).
#' @export
semiOPBARTLocalConfusionDS <- function(pred_obj, levels_Serialize, nfilter = 5,
                                        validation_col = "null",
                                        validation_data.name = "null",
                                        dichotomize_threshold = "null") {
  # a site can legitimately have no prediction object at all -- e.g. it
  # never got a Group A/Group B split object big enough to clear nfilter
  # in the first place (semiOPBARTLocalSplitDS()'s auto-demotion), so
  # predict was never even attempted there. That is exactly the same
  # "this site can't participate" case validation_data.name/validation_col
  # below already handle by returning NULL gracefully -- pred_obj needs
  # the identical guard, not a hard stop, or one such site crashes
  # the WHOLE evaluation run instead of just being excluded from it.
  if (!exists(pred_obj, envir = parent.frame(), inherits = FALSE))
    return(NULL)   # no prediction object at this site -- nothing to evaluate
  p <- get(pred_obj, envir = parent.frame())
  levels <- semiOPBART_fromSerialize(levels_Serialize)
  if (identical(validation_col, "null")) validation_col <- NULL
  if (identical(validation_data.name, "null")) validation_data.name <- NULL

  if (is.null(validation_col)) {
    truth_raw <- p$truth
  } else {
    if (is.null(validation_data.name) ||
        !exists(validation_data.name, envir = parent.frame(), inherits = FALSE))
      return(NULL)   # no such object at this site -- no validation possible
    vdf <- get(validation_data.name, envir = parent.frame())
    if (!validation_col %in% names(vdf)) return(NULL)   # column absent here
    truth_raw <- vdf[[validation_col]]
    if (length(truth_raw) != length(p$map))
      return(NULL)   # row-count mismatch between predictions and this
                      # validation source -- cannot align them safely
  }

  truth_f <- factor(as.character(truth_raw), levels = levels)
  pred_f  <- factor(as.character(p$map),     levels = levels)
  keep <- !is.na(truth_f)   # rows with no validation value simply don't
                             # contribute to the confusion matrix, rather
                             # than being coerced into a spurious category
  if (sum(keep) < nfilter) return(NULL)
  truth_f <- truth_f[keep]; pred_f <- pred_f[keep]

  cm <- table(truth_f, pred_f)
  cm_reported <- cm
  cm_reported[cm_reported > 0 & cm_reported < nfilter] <- NA

  print(paste0("truth_f : ", truth_f))
  print(paste0("pred_f : ", pred_f))

  if (identical(dichotomize_threshold, "null")) {
    # ORIGINAL behaviour, unchanged: baseline category (levels[1]) vs rest
    truth_bin <- factor(ifelse(as.character(truth_f) == levels[1], "0", "1"),
                         levels = c("0", "1"))
    pred_bin  <- factor(ifelse(as.character(pred_f)  == levels[1], "0", "1"),
                         levels = c("0", "1"))
  } else {
    thr <- as.numeric(dichotomize_threshold)
    truth_num <- suppressWarnings(as.numeric(as.character(truth_f)))
    pred_num  <- suppressWarnings(as.numeric(as.character(pred_f)))
    if (anyNA(truth_num) || anyNA(pred_num))
      stop("dichotomize_threshold requires numeric-coercible levels -- ",
           "got non-numeric values in truth or predicted labels")
    truth_bin <- factor(ifelse(truth_num >= thr, "1", "0"), levels = c("0", "1"))
    pred_bin  <- factor(ifelse(pred_num  >= thr, "1", "0"), levels = c("0", "1"))
  }
  cm_bin <- table(truth_bin, pred_bin)
  cm_bin_reported <- cm_bin
  cm_bin_reported[cm_bin_reported > 0 & cm_bin_reported < nfilter] <- NA

  list(confusion = cm, confusion_reported = cm_reported,
       confusion_bin = cm_bin, confusion_bin_reported = cm_bin_reported,
       n = sum(cm))
}






