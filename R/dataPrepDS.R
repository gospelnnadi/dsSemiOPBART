# dsSemiOPBARTDataPrep.R
# ---------------------------------------------------------------------------
# Data preparation is now three explicit steps rather than one, so that
# normalization can be swapped independently of feature engineering:
#
#   1. TRANSFORM  -- vetted feature engineering + dummy encoding + outcome
#                    factor coercion. Produces UNNORMALISED X, plus W and Y.
#   2. RANGE       -- (only if normalize_method = "federated_minmax") each
#                    site reports its local per-column min/max on the
#                    transformed (not yet normalised) X; client aggregates
#                    to a single GLOBAL min/max per column and broadcasts it
#                    back. Skipped entirely for "local_ecdf".
#   3. NORMALIZE   -- applies either the ORIGINAL per-site local ECDF
#                    (matches the non-federated smopbart() exactly, but is
#                    only a same-site approximation of the pooled
#                    distribution -- the bias flagged early in this
#                    conversation) or the newly-broadcast federated
#                    min/max (every site normalises against the SAME global
#                    min/max, removing that bias, at a much smaller
#                    disclosure cost than the quantile-grid alternative
#                    floated earlier: two numbers per column instead of a
#                    whole grid).
#
# ds.semiOPBARTPrepare() below still offers the single-call convenience
# interface from before; it just now orchestrates three round trips instead
# of one when normalize_method = "federated_minmax".
#
# DISCLOSURE NOTE ON federated_minmax, READ BEFORE USING IT:
# the global min and global max ARE, by construction, an exact covariate
# value belonging to one real patient at one real site (whichever site
# happens to hold the true extreme). That is a smaller disclosure surface
# than publishing a full quantile grid (two numbers per column instead of
# 20-50), but each of those two numbers is maximally identifying FOR THAT
# ONE RECORD -- it is by definition the most unusual value for that
# covariate in the whole federation. If any covariate is itself sensitive
# or rare (e.g. an unusually high/low lab value that could re-identify a
# patient in combination with other public knowledge), don't default to
# this option for that covariate. It is offered as opt-in, per this
# request, not as the new default.
# ---------------------------------------------------------------------------

#source("serializeDS.R")

# ---- vetted transform registry (extend by adding named entries here) ------

semiOPBARTTransformRegistry <- list(
  rheumatoid_ip_ep = function(d) {
    within(d, {
      bclass_IL6i  <- as.integer(bDMARD == 2)
      CRP_collapsed <- as.integer(10 * CRP < ESR)
      CSR_ESR_GAP  <- abs(ESR - SJC28)
      CRP_noscale  <- CRP  +
        (0.2 * CSR_ESR_GAP * CRP_collapsed) +
        (0.6 * CSR_ESR_GAP * bclass_IL6i * as.integer(CRP == 0))
      CRP <- CRP_noscale
    })
  },
  none = function(d) d
)


#
# ---- server side (Opal) ----------------------------------------------------
#' semiOPBARTLocalClearStateDS
#' Remove semiOPBART state from the site.
#'
#' PREVIOUSLY this hardcoded a fixed list of DEFAULT parameter names
#' ("semiOPBART_train" etc.) -- which are NOT the names
#' run_semiOPBART_orchestrator() actually uses ("D_train" etc., let alone
#' any run-tagged variant of them), so calling ds.semiOPBARTClearState()
#' silently cleared nothing real and stale state from a previous run (a
#' different pathology/model, or a different architecture) could leak
#' into the next one -- exactly the "federated_ecdf columns don't match"
#' failure this was built to fix. Now it removes EXACTLY the names it's
#' given, nothing guessed.
#'
#' @param names_Serialize  Serialize-encoded character vector of object
#'   names to remove, built client-side by ds.semiOPBARTClearState() from
#'   whatever run_tag/architecture-specific names that run actually used
#' @export
semiOPBARTLocalClearStateDS <- function(names_Serialize = "null") {
  rm_if_exists <- function(nm, envir) {
    if (exists(nm, envir = envir, inherits = FALSE)) {
      rm(list = nm, envir = envir)
      TRUE
    } else {
      FALSE
    }
  }

  names_to_clear <- if (identical(names_Serialize, "null")) {
    # no explicit list given -- fall back to the legacy default-parameter
    # names, purely defensive (kept for anyone still relying on the old
    # unnamed/untagged single-run convention)
    c("semiOPBART_transformed", "semiOPBART_train", "semiOPBART_test",
      "semiOPBART_holdout", "federated_ecdf", "global_range", ".semiOPBART_state")
  } else {
    semiOPBART_fromSerialize(names_Serialize)
  }

  cleared_frame <- vapply(names_to_clear, rm_if_exists, logical(1), envir = parent.frame())

  list(cleared = sum(cleared_frame), names_cleared = names_to_clear[cleared_frame])
}

# =====================  1. TRANSFORM  ========================================

# ---- server side (Opal) ----------------------------------------------------
#' @name semiOPBARTLocalTransformDS
#' @param data.name         raw site-side data frame name
#' @param outcome_col       ordinal outcome column name
#' @param levels_Serialize       Serialize-encoded full outcome level set, e.g. 0:4
#'   (analyst-specified, not data-derived -- see prior discussion)
#' @param x_features_Serialize   Serialize-encoded character vector, nonparametric
#'   (tree) feature names
#' @param w_features_Serialize   Serialize-encoded character vector, parametric
#'   (linear) feature names
#' @param transform_recipe  name of a registered, vetted transform, or "none"
#' @param newobj            name under which the UNNORMALISED transformed
#'                          object (a list with $X, $W, $Y, $dv) is stored
#' @param nfilter           minimum site n to proceed at all
#' @export
semiOPBARTLocalTransformDS <- function(data.name, outcome_col, levels_Serialize,
                                        x_features_Serialize, w_features_Serialize,
                                        transform_recipe = "none",
                                        newobj = "semiOPBART_transformed",
                                        nfilter = 5) {
  df <- eval(parse(text = data.name), envir = parent.frame())

  
  if (nrow(df) < nfilter) stop("site n below disclosure threshold")

  levels     <- semiOPBART_fromSerialize(levels_Serialize)
  x_features <- semiOPBART_fromSerialize(x_features_Serialize)
  w_features <- semiOPBART_fromSerialize(w_features_Serialize)

  if (!transform_recipe %in% names(semiOPBARTTransformRegistry))
    stop("transform_recipe '", transform_recipe, "' is not registered/vetted. ",
         "Options: ", paste(names(semiOPBARTTransformRegistry), collapse = ", "))
  df <- semiOPBARTTransformRegistry[[transform_recipe]](df)


  # ---- NEW: use only local available features ----
  x_features <- intersect(x_features, names(df))
  w_features <- intersect(w_features, names(df))

  if (length(x_features) == 0)
    stop("No x_features available in local dataframe.")

  if (length(w_features) == 0)
    stop("No w_features available in local dataframe.")

  # ---- PROTECTIVE OUTCOME HANDLING ---- 
  # The outcome is optional at a site. If it is absent, Y is explicitly NULL.
  outcome_available <- !is.null(outcome_col) && length(outcome_col) == 1L && nzchar(outcome_col) && outcome_col %in% names(df) 
  if (outcome_available) { 
    df[[outcome_col]] <- factor(df[[outcome_col]], levels = levels) 
    Y <- df[[outcome_col]] 
    cat( "TRANSFORM: outcome available |", "column =", outcome_col, "| length(Y) =", length(Y), "| NA(Y) =", sum(is.na(Y)), "\n" ) 
    } else { 
      Y <- NULL 
      cat( "TRANSFORM: outcome unavailable locally |", "requested outcome_col =", outcome_col, "| Y = NULL\n" ) 
      }
  #df[[outcome_col]] <- factor(df[[outcome_col]], levels = levels)
  # UNLIKE THE ORIGINAL FILTER HERE: unlabeled rows (outcome_col == NA) are
  # KEPT, not dropped. A site that turns out to have too few labeled rows
  # to train can still be used for PREDICTION -- see semiOPBARTLocalSplitDS()
  # below, which is where labeled-only filtering now happens, specifically
  # for the TRAIN object. The TEST object may legitimately contain
  # unlabeled rows (no training label at all, only prediction + optional
  # validation against a separate variable -- see dsSemiOPBARTEvaluate.R).

  dv <- caret::dummyVars(as.formula(paste("~", paste(x_features, collapse = " + "))), df)
  X  <- suppressWarnings(predict(dv, df))    # UNNORMALISED -- normalization
                                              # happens in step 3, not here
  W  <- as.matrix(df[, w_features, drop = FALSE])
  #Y  <- df[[outcome_col]]

  print(paste0("semiOPBARTLocalTransformDS: colnames(X) = ", paste(colnames(X), collapse = ", "), ", colnames(W) = ", paste(colnames(W), collapse = ", "), ", length(Y) = ", length(Y)))

  base::assign(newobj, list(X = X, W = W, Y = Y, dv = dv), envir = parent.frame())
  print(paste0("semiOPBARTLocalTransformDS: base::assigned newobj '", newobj, "' in parent.frame()"))

  # column names ARE returned (server -> client direction, not subject to
  # the inbound Serialize restriction) so the client can verify every site
  # produced the SAME dummy-encoded column set before proceeding -- a
  # cheap, useful consistency check given that dummyVars' output columns
  # depend on which factor levels are actually present locally.
  list(n = nrow(df), x_cols = colnames(X), w_cols = colnames(W))
}




# =====================  2. FEDERATED RANGE  ==================================
# Only needed for normalize_method = "federated_minmax". Skip entirely for
# "local_ecdf".

# ---- server side (Opal) ----------------------------------------------------
#' @name semiOPBARTLocalFeatureRangeDS
#' Report this site's local per-column min/max on X. As of the
#' normalization-ordering fƒix, this is called on the TRAIN object (across
#' train_sites only), not the full pre-split transformed object -- a
#' site's own held-out test/holdout rows no longer contribute to the
#' global range it will later be normalized against, mirroring the same
#' train-only principle applied to local_ecdf.
#' @export
semiOPBARTLocalFeatureRangeDS <- function(data.name = "semiOPBART_train") {
  s <- eval(parse(text = data.name), envir = parent.frame())
  

  list(cols = colnames(s$X),
       min  = apply(s$X, 2, min, na.rm = TRUE),
       max  = apply(s$X, 2, max, na.rm = TRUE))
}

# =====================  2b. FEDERATED HISTOGRAM (for federated_ecdf) =========
#
# federated_minmax anchors only the two ENDPOINTS of the normalized scale
# -- it says nothing about the SHAPE of the pooled distribution in
# between. If sites have very different marginal distributions of a
# covariate, min-max normalization keeps tree-sharing coherent (any
# consistent global rescaling does) but doesn't give the roughly-uniform
# spread that the real algorithm's local ecdf() normalization is designed
# to produce. This section builds a genuinely federated approximation to
# the POOLED ecdf, without ever moving a raw value: sites agree on shared
# bin edges (from the already-computed global range), report LOCAL COUNTS
# per bin, and those counts are summed across sites -- an aggregate
# statistic, structurally identical to any ds.glm-style sufficient
# statistic. A bin count is a much smaller disclosure surface than an
# exact min/max: it summarizes a whole interval of values, not one
# identifiable extreme.

# ---- server side (Opal) ----------------------------------------------------

#' Report this site's local counts within a SHARED set of bin edges (the
#' same edges at every site, so summing is meaningful).
#'
#' @param data.name        the TRAIN object
#' @param bin_edges_Serialize   Serialize-encoded named list, one numeric vector of
#'   bin edges (length num_bins + 1) per column -- built once, client-side,
#'   from the global range, and broadcast to every site so counts line up
#' @export
semiOPBARTLocalHistogramDS <- function(data.name = "semiOPBART_train",
                                        bin_edges_Serialize) {
  s <- eval(parse(text = data.name), envir = parent.frame())
    
  
  bin_edges <- semiOPBART_fromSerialize(bin_edges_Serialize)
  cols <- colnames(s$X)
  counts <- lapply(cols, function(cn) {
    edges <- bin_edges[[cn]]
    as.integer(table(cut(s$X[, cn], breaks = edges, include.lowest = TRUE)))
  })
  names(counts) <- cols
  list(cols = cols, counts = counts)
}

# =====================  3. NORMALIZE (runs AFTER split, not before) =========
#
# CONFIRMED against the real smopbart() source: ECDFs are fit on X_train
# ONLY, then applied out-of-sample to X_test. This section replaces an
# earlier version that normalised BEFORE splitting -- which let a site's
# own test/holdout rows influence that SAME site's local_ecdf scale, a
# real (if mild) fidelity deviation, not a privacy issue. Normalize now
# takes the ALREADY-SPLIT train/test/holdout object names together, fits
# on train only (local_ecdf) or applies the shipped global range
# (federated_minmax) to whichever of the three objects exist at this site.

# ---- server side (Opal) ----------------------------------------------------
#' @name semiOPBARTLocalNormalizeSplitDS
#' @param train.name, test.name, holdout.name  names of the (possibly
#'   partially absent -- e.g. a role="test" site has no train object)
#'   SPLIT objects at this site, from semiOPBARTLocalSplitDS()/
#'   semiOPBARTLocalAsTestDS()
#' @param normalize_method  "local_ecdf": requires a train object at this
#'   site (fits ECDFs there, applies to test/holdout too, exactly matching
#'   the real algorithm's train-only fit). "federated_minmax": requires
#'   global_min_Serialize/global_max_Serialize, applies to whichever objects exist,
#'   train not required (this is what makes never-trained-site prediction
#'   possible at all -- see dsSemiOPBARTExternalPredict.R). "federated_ecdf":
#'   requires federated_ecdf_Serialize from ds.semiOPBARTComputeGlobalECDF() --
#'   a genuinely pooled (all-sites) empirical CDF approximation, built
#'   from aggregated histogram bin counts, not just two endpoints.
#' @param federated_ecdf_Serialize  Serialize-encoded list(bin_edges=, cum_frac=) per
#'   column, from ds.semiOPBARTComputeGlobalECDF() -- required (and only
#'   used) when normalize_method = "federated_ecdf"
#' @export
semiOPBARTLocalNormalizeSplitDS <- function(train.name = "semiOPBART_train",
                                             test.name = "semiOPBART_test",
                                             holdout.name = "semiOPBART_holdout",
                                             normalize_method = "local_ecdf",
                                             federated_ecdf_Serialize = "null",
                                             nfilter = 5) {

  if (identical(federated_ecdf_Serialize, "null")) federated_ecdf_Serialize <- NULL
  caller_env <- parent.frame()
  has <- function(nm)!identical(nm, "null") && exists(nm, envir = caller_env) #exists(parse(text = nm), envir = parent.frame()) #
  has_train <- has(train.name); has_test <- has(test.name); has_holdout <- has(holdout.name)

  if (normalize_method == "local_ecdf") {
    # if (!has_train)
    #   stop("local_ecdf normalization needs a TRAIN object at this site to ",
    #        "fit ECDFs on (this site has none -- fully demoted or role='test'); ",
    #        "use normalize_method = 'federated_minmax' for sites that never train")
    if (has_train){
    s_train <- eval(parse(text = train.name), envir = parent.frame())
    print(paste0("semiOPBARTLocalNormalizeSplitDS: dim(s_train$X) = ", paste(dim(s_train$X), collapse = ", ")))

    X_train <- s_train$X
    ecdfs <- setNames( lapply(seq_len(ncol(X_train)), function(j) {
      u <- unique(X_train[, j])
      if (length(u) == 1) return(identity)
      if (length(u) == 2) {
        a <- min(X_train[, j]); b <- max(X_train[, j])
        return(function(y) (y - a) / (b - a))
      }
      ecdf(X_train[, j])

    }),
    colnames(X_train)
    )
    print("names(ecdfs):")
    print(names(ecdfs))
    apply_norm <- function(X) { for (j in seq_len(ncol(X))) X[, j] <- ecdfs[[j]](X[, j]); X }
    norm_info <- list(method = "local_ecdf", ecdfs = ecdfs)
    }
  } else if (normalize_method == "federated_ecdf") {
    if (is.null(federated_ecdf_Serialize))
      stop("federated_ecdf requires federated_ecdf_Serialize from ",
           "ds.semiOPBARTComputeGlobalECDF()")
    federated_ecdf <- semiOPBART_fromSerialize(federated_ecdf_Serialize)
    ecdf_fns <- setNames(lapply(federated_ecdf$cols, function(cn) {
      edges <- federated_ecdf$bin_edges[[cn]]; cf <- federated_ecdf$cum_frac[[cn]]
      function(y) pmin(pmax(stats::approx(edges, cf, xout = y, method = "linear",
                                           rule = 2)$y, 0), 1)
    }), federated_ecdf$cols)


    apply_norm <- function(X) {

        print("setdiff(federated_ecdf$cols, colnames(X))")
        print(setdiff(federated_ecdf$cols, colnames(X)))

        common_cols <- intersect(colnames(X), names(ecdf_fns))

        print("intersect(colnames(X), federated_ecdf$cols)")
        print(common_cols)

        if (length(common_cols) == 0) {
          warning(
            "No intersecting columns between federated_ecdf and this site's X; ",
            "no normalization applied"
          )
         return(X)
        }

        for (j in common_cols) {
        X[, j] <- ecdf_fns[[j]](X[, j])
        }

      X
      }

    norm_info <- list(method = "federated_ecdf", federated_ecdf = federated_ecdf)

  } else stop("normalize_method must be 'local_ecdf', ",
              "or 'federated_ecdf'")
 

    update_obj <- function(nm) {

        s <- get(
            nm,
            envir = caller_env,
            inherits = FALSE
        )

        s$X <- apply_norm(s$X)
        s$norm_info <- norm_info
      cat(
    "SPLIT:",
    "nm =", nm,
    "| dim(X) =", paste(dim(s$X), collapse = "x"),
    "\n"
    )
        base::assign(
            nm,
            s,
            envir = caller_env
        )
      s <- eval(parse(text = nm), envir = caller_env)
    print(paste0("semiOPBARTLocalNormalizeSplitDS: updated ", nm, ": dim(s$X) = ", paste(dim(s$X), collapse = ", "), ", length(s$Y) = ", length(s$Y)))
     s
    }

  if (has_train)   update_obj(train.name)
  if (has_test) { 
    s_train <-  update_obj(test.name)
    print(paste0( "semiOPBARTLocalSplitDS:  s$Y labels = ", s_train$Y ))
  }
  if (has_holdout && has_train) update_obj(holdout.name)
  


  list(normalize_method = normalize_method,
       has_train = has_train, has_test = has_test, has_holdout = has_holdout)
}



# ---- server side (Opal) ----------------------------------------------------
#' @name semiOPBARTLocalSplitDS
#' @param data.name    the prepared object
#' @param outcome_col  outcome column name (documentary; labels are in s$Y)
#' @param train_ratio  NULL -> ALL labeled rows become training data;
#'                     numeric in (0,1) -> stratified split
#' @param seed         random seed
#' @param balance_classes
#'                     if TRUE, balance classes 0-4 in TRAIN ONLY
#' @param newobj_train
#' @param newobj_test
#' @param newobj_holdout
#' @param nfilter      minimum labeled row count required to train
#'
#' Group A:
#'   newobj_test = held-out labeled rows only
#'
#' Group B:
#'   newobj_holdout = Group A + all unlabeled rows
#'
#' Balancing is performed ONLY on the training population.
#' Test and holdout populations are never replicated or modified.
#'
#' @export
semiOPBARTLocalSplitDS <- function(data.name, outcome_col, train_ratio = NULL,
                                    seed = NULL,
                                    newobj_train = "semiOPBART_train",
                                    newobj_test = "semiOPBART_test",
                                    newobj_holdout = "semiOPBART_holdout",
                                    nfilter = 5,levels=0:4,
                                    balance_classes = TRUE, #FALSE,
                                    balance_target = 5 , 
                                    balance_seed = NULL) {

  s <- eval(parse(text = data.name), envir = parent.frame())

  labeled_idx   <- which(!is.na(s$Y))
  unlabeled_idx <- which(is.na(s$Y))
  n_labeled <- length(labeled_idx)
  s$level <- levels

  cat(
    "SOURCE:",
    data.name,
    "| exists =", exists(data.name, parent.frame()),
    "\n"
  )

  if (exists(data.name, parent.frame())) {
    tmp <- get(data.name, parent.frame())

    cat(
      "SOURCE DIM:",
      paste(dim(tmp$X), collapse = "x"),
      "| W:",
      paste(dim(tmp$W), collapse = "x"),
      "| Y:",
      length(tmp$Y),
      "\n"
    )
  }

  print(
    paste0(
      "semiOPBARTLocalSplitDS: n_labeled = ", n_labeled,
      ", n_unlabeled = ", length(unlabeled_idx),
      ", dim(s$X) = ", paste(dim(s$X), collapse = "x"),
      ", dim(s$W) = ", paste(dim(s$W), collapse = "x"),
      ", length(s$Y) = ", length(s$Y)
    )
  )
 
  slice <- function(idx) {
    list(
      X = s$X[idx, , drop = FALSE],
      W = s$W[idx, , drop = FALSE],
      Y = s$Y[idx],
      dv = s$dv,
      norm_info = s$norm_info
    )
  }

  # ------------------------------------------------------------
  # Too few labeled observations -> auto-demote to holdout
  # ------------------------------------------------------------

  if (n_labeled < nfilter) {

    base::assign(
      newobj_holdout,
      s,
      envir = parent.frame()
    )

    print(
      paste0(
        "semiOPBARTLocalSplitDS: auto-demoted, n_labeled = ",
        n_labeled,
        " < nfilter = ",
        nfilter
      )
    )

    return(
      list(
        n_train = 0L,
        n_test = 0L,
        n_holdout = length(s$Y),
        split_applied = FALSE,
        auto_demoted = TRUE,
        reason = paste0(
          "only ", n_labeled,
          " labeled rows, need >= ", nfilter
        )
      )
    )
  }

  # ------------------------------------------------------------
  # No train/test split requested
  # ------------------------------------------------------------

  if (is.null(train_ratio)) {

    print(
      "semiOPBARTLocalSplitDS: train_ratio is NULL, using ALL labeled rows for training"
    )

    s_train <- slice(labeled_idx)

    # ----------------------------------------------------------
    # BALANCE TRAINING SET ONLY
    # ----------------------------------------------------------

    if (balance_classes) {

      if (is.null(balance_seed))
        balance_seed <- seed

      bal <- semiOPBART_balance_train(
        s = s_train,
        levels =levels,
        target = balance_target,
        seed = balance_seed
      )

      cat(
        "TRAIN BALANCE BEFORE:",
        paste(
          names(bal$counts_before),
          as.integer(bal$counts_before),
          sep = "=",
          collapse = ", "
        ),
        "\n"
      )

      cat(
        "TRAIN BALANCE AFTER:",
        paste(
          names(bal$counts_after),
          as.integer(bal$counts_after),
          sep = "=",
          collapse = ", "
        ),
        "\n"
      )

      s_train <- bal$data
    }

    base::assign(
      newobj_train,
      s_train,
      envir = parent.frame()
    )

    n_holdout <- 0L

    if (length(unlabeled_idx) >= nfilter) {

      base::assign(
        newobj_holdout,
        slice(unlabeled_idx),
        envir = parent.frame()
      )

      n_holdout <- length(unlabeled_idx)
    }

    return(
      list(
        n_train = length(s_train$Y),
        n_test = 0L,
        n_holdout = n_holdout,
        split_applied = FALSE,
        auto_demoted = FALSE
      )
    )
  }

  # ------------------------------------------------------------
  # Standard train/test split
  # ------------------------------------------------------------

  stopifnot(
    "train_ratio must be in (0,1)" =
      train_ratio > 0 && train_ratio < 1
  )

  if (!is.null(seed))
    set.seed(seed)

  train_sub <- unlist(
    caret::createDataPartition(
      s$Y[labeled_idx],
      p = train_ratio,
      list = TRUE
    )
  )

  train_idx   <- labeled_idx[train_sub]
  test_idx    <- setdiff(labeled_idx, train_idx)
  holdout_idx <- c(test_idx, unlabeled_idx)

  print(
    paste0(
      "semiOPBARTLocalSplitDS: train_idx = ",
      length(train_idx),
      ", test_idx = ",
      length(test_idx),
      ", holdout_idx = ",
      length(holdout_idx)
    )
  )

  print("semiOPBARTLocalSplitDS: test_Y index =")
  print(test_idx)

  # ------------------------------------------------------------
  # TRAIN
  # ------------------------------------------------------------

  s_train <- slice(train_idx)

  if (length(s_train$Y) < nfilter) {

    stop(
      "train subset below disclosure threshold at this site even though ",
      "labeled count passed the initial check -- train_ratio too low"
    )
  }

  # ------------------------------------------------------------
  # BALANCE TRAINING SET ONLY
  # ------------------------------------------------------------

  if (balance_classes) {

    if (is.null(balance_seed))
      balance_seed <- seed

    bal <- semiOPBART_balance_train(
      s = s_train,
      levels = levels,
      target = balance_target,
      seed = balance_seed
    )

    cat(
      "TRAIN BALANCE BEFORE:",
      paste(
        names(bal$counts_before),
        as.integer(bal$counts_before),
        sep = "=",
        collapse = ", "
      ),
      "\n"
    )

    cat(
      "TRAIN BALANCE AFTER:",
      paste(
        names(bal$counts_after),
        as.integer(bal$counts_after),
        sep = "=",
        collapse = ", "
      ),
      "\n"
    )

    s_train <- bal$data
  }

  base::assign(
    newobj_train,
    s_train,
    envir = parent.frame()
  )

  # ------------------------------------------------------------
  # TEST -- NEVER BALANCE
  # ------------------------------------------------------------

  n_test <- 0L

  if (length(test_idx) >= nfilter) {

    base::assign(
      newobj_test,
      slice(test_idx),
      envir = parent.frame()
    )

    

    n_test <- length(test_idx)
  }

  # ------------------------------------------------------------
  # HOLDOUT -- NEVER BALANCE
  # ------------------------------------------------------------

  n_holdout <- 0L

  if (length(holdout_idx) >= nfilter) {

    base::assign(
      newobj_holdout,
      slice(holdout_idx),
      envir = parent.frame()
    )

    n_holdout <- length(holdout_idx)
  }

  # ------------------------------------------------------------
  # Diagnostics
  # ------------------------------------------------------------
  cat(
    "LocalSPLIT:",
    "data.name =", data.name,
    "| newobj_train =", newobj_train,
    "| newobj_test =", newobj_test,
    "| newobj_holdout =", newobj_holdout,
    "| dim(X) =", paste(dim(s$X), collapse = "x"),
    "| dim(W) =", paste(dim(s$W), collapse = "x"),
    "| length(Y) =", length(s$Y),
    "\n"
  )

  cat(
    "LOADED:",
    data.name,
    "| dim(X) =", paste(dim(s$X), collapse = "x"),
    "| dim(W) =", paste(dim(s$W), collapse = "x"),
    "| length(Y) =", length(s$Y),
    "| pid =", Sys.getpid(),
    "| host =", Sys.info()[["nodename"]],
    "\n"
  )

  list(
    n_train = length(s_train$Y),
    n_test = n_test,
    n_holdout = n_holdout,
    split_applied = TRUE,
    auto_demoted = FALSE
  )
}

# ---- server side (Opal) ----------------------------------------------------
#' Balance training indices only
#'
#' Replicates observations WITH replacement so that every class in
#' `levels` has the same number of observations as the largest class.
#'
#' IMPORTANT:
#'   This function returns row indices only.
#'   The caller applies these indices to X, W and Y simultaneously.
#'
#' @export
semiOPBART_balance_train <- function(s,
                                     levels = 0:4,
                                     target = 8,
                                     seed = NULL) {

  if (!is.null(seed))
    set.seed(seed)

  y <- as.character(s$Y)
  lev <- as.character(levels)

  # Only labeled observations should reach this function
  keep <- which(!is.na(y))

  if (!length(keep))
    stop("No labeled observations available for balancing.")

  counts_before <- table(
    factor(y[keep], levels = lev)
  )

  # Make sure all five classes exist
  if (any(counts_before == 0)) {
    stop(
      "Cannot balance training data: missing class(es): ",
      paste(
        lev[counts_before == 0],
        collapse = ", "
      )
    )
  }

  # Never downsample classes already above target
  sampled_idx <- lapply(lev, function(cl) {

    idx <- keep[y[keep] == cl]

    if (length(idx) >= target) {
      idx
    } else {
      c(
        idx,
        sample(
          idx,
          size = target - length(idx),
          replace = TRUE
        )
      )
    }
  })

  balanced_idx <- unlist(sampled_idx, use.names = FALSE)

  # Shuffle so replicated observations are not grouped by class
  balanced_idx <- sample(balanced_idx)

  s_balanced <- list(
    X = s$X[balanced_idx, , drop = FALSE],
    W = s$W[balanced_idx, , drop = FALSE],
    Y = s$Y[balanced_idx],
    dv = s$dv,
    norm_info = s$norm_info
  )

  counts_after <- table(
    factor(
      as.character(s_balanced$Y),
      levels = lev
    )
  )

  list(
    data = s_balanced,
    idx = balanced_idx,
    counts_before = counts_before,
    counts_after = counts_after
  )
}










# ============================================================
# POOLED RARE-CLASS FILTERING
# ============================================================
#' @export
semiOPBARTLocalClassCountsDS <- function(data.name, outcome_col, levels,
                                         nfilter = 5) {
  df <- eval(parse(text = data.name), envir = parent.frame())

  if (nrow(df) < nfilter) {
    stop("site n below disclosure threshold")
  }

  if (!outcome_col %in% names(df)) {
    stop("Outcome column not found: ", outcome_col)
  }

  lev <- as.character(levels)
  y   <- as.character(df[[outcome_col]])

  # IMPORTANT:
  # Always count against the analyst-supplied levels.
  # This guarantees identical output structure across sites.
  counts <- table(
    factor(y, levels = lev),
    useNA = "no"
  )

  list(
    levels = lev,
    counts = as.numeric(counts)
  )
}

#' @export
semiOPBARTLocalRemoveClassesDS <- function(data.name,
                                           outcome_col,
                                           remove_levels,
                                           newobj = "semiOPBART_filtered",
                                           nfilter = 5) {

  s <- eval(parse(text = data.name), envir = parent.frame())

  if (!is.data.frame(s))
    stop("data.name must refer to a data.frame")

  if (!outcome_col %in% names(s))
    stop("outcome_col '", outcome_col, "' not found in data")

  y_chr <- as.character(s[[outcome_col]])
  remove_levels <- as.character(remove_levels)

  keep <- !is.na(y_chr) & !(y_chr %in% remove_levels)

  # Keep unlabeled rows. They are not part of the class-count
  # decision and may be needed later for the holdout population.
  keep[is.na(y_chr)] <- TRUE

  s_filtered <- s[keep, , drop = FALSE]

  base::assign(
    newobj,
    s_filtered,
    envir = parent.frame()
  )

  list(
    n_before = nrow(s),
    n_after = nrow(s_filtered),
    n_removed = sum(!keep),
    removed_classes = remove_levels
  )
}












#' @name semiOPBARTLocalAsTestDS
#' @param newobj_test     Group A (labeled rows only) -- created ONLY if
#'   this site has >= nfilter labeled rows, since a role="test" site
#'   (typically an external validation cohort) may well have real outcomes
#'   available even though it never trained
#' @param newobj_holdout  Group B -- ALL of this site's rows, unconditionally
#' @export
semiOPBARTLocalAsTestDS <- function(data.name, newobj_test = "semiOPBART_test",
                                     newobj_holdout = "semiOPBART_holdout",
                                     nfilter = 5) {

  s <- eval(parse(text = data.name), envir = parent.frame())
  if (length(s$Y) < nfilter) stop("site n below disclosure threshold")

  labeled_idx <- which(!is.na(s$Y))
  n_test <- 0L
  if (length(labeled_idx) >= nfilter) {
    base::assign(newobj_test, list(X = s$X[labeled_idx, , drop = FALSE],
                              W = s$W[labeled_idx, , drop = FALSE],
                              Y = s$Y[labeled_idx], dv = s$dv, norm_info = s$norm_info),
           envir = parent.frame())

    n_test <- length(labeled_idx)
  }

  base::assign(newobj_holdout, s, envir = parent.frame())   # ALL rows, unconditionally
  list(n_test = n_test, n_holdout = length(s$Y))
}






