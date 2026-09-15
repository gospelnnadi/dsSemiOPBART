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
library(SoftBart)

#source("serializeDS.R")

#' Store one chunk of an external reference on the server.
#'
#' @export
semiOPBARTLocalStoreExternalReferenceChunkDS <- function(
    reference_name,
    kind,
    chunk,
    chunk_index,
    n_chunks
) {

  caller_env <- parent.frame()

  if (!kind %in% c("trees", "ecdf")) {
    stop(
      "kind must be either 'trees' or 'ecdf'"
    )
  }

  if (
    length(chunk) != 1L ||
    is.na(chunk) ||
    !nzchar(chunk)
  ) {
    stop(
      "Received invalid chunk: it must be one non-empty non-NA string"
    )
  }

  if (
    chunk_index < 1L ||
    chunk_index > n_chunks
  ) {
    stop("Invalid chunk index")
  }

  staged_name <- paste0(
    ".semiOPBART_external_reference_",
    reference_name
  )

  if (!exists(
      staged_name,
      envir = caller_env,
      inherits = FALSE
  )) {

    ref <- list(
      tree_chunks = vector("list", n_chunks),
      ecdf_chunks = vector("list", n_chunks),
      n_tree_chunks = NA_integer_,
      n_ecdf_chunks = NA_integer_
    )

  } else {

    ref <- get(
      staged_name,
      envir = caller_env,
      inherits = FALSE
    )

    # Expand if the first stream had fewer chunks.
    if (
      kind == "trees" &&
      length(ref$tree_chunks) < n_chunks
    ) {
      length(ref$tree_chunks) <- n_chunks
    }

    if (
      kind == "ecdf" &&
      length(ref$ecdf_chunks) < n_chunks
    ) {
      length(ref$ecdf_chunks) <- n_chunks
    }
  }

  if (kind == "trees") {

    ref$tree_chunks[[chunk_index]] <- chunk
    ref$n_tree_chunks <- as.integer(n_chunks)

  } else {

    ref$ecdf_chunks[[chunk_index]] <- chunk
    ref$n_ecdf_chunks <- as.integer(n_chunks)
  }

  base::assign(
    staged_name,
    ref,
    envir = caller_env
  )

  list(
    ok = TRUE,
    reference_name = reference_name,
    kind = kind,
    chunk_index = chunk_index,
    n_chunks = n_chunks
  )
}

# ---- server side (Opal) ----------------------------------------------------
#' Store an externally shipped forest + normalization reference locally.
#'
#' The payload is transferred in chunks so a large Serialize string does not
#' have to be embedded as one huge literal in a DataSHIELD call.
#'
#' @export
semiOPBARTLocalStoreExternalReferenceDS <- function(
    reference_name,
    trees_chunk,
    ecdf_chunk,
    chunk_index,
    n_chunks
) {

  caller_env <- parent.frame()

  state_name <- paste0(
    ".semiOPBART_external_reference_",
    reference_name
  )

  if (
    !exists(
      state_name,
      envir = caller_env,
      inherits = FALSE
    )
  ) {
    ref <- list(
      trees_chunks = vector("list", n_chunks),
      ecdf_chunks = vector("list", n_chunks),
      n_chunks = n_chunks
    )
  } else {
    ref <- get(
      state_name,
      envir = caller_env,
      inherits = FALSE
    )
  }

  ref$trees_chunks[[chunk_index]] <- trees_chunk
  ref$ecdf_chunks[[chunk_index]] <- ecdf_chunk

  base::assign(
    state_name,
    ref,
    envir = caller_env
  )

  list(
    ok = TRUE,
    reference_name = reference_name,
    chunk_index = chunk_index
  )
}

#' Assemble a complete external reference from staged chunks.
#'
#' @export
semiOPBARTLocalAssembleExternalReferenceDS <- function(
    reference_name
) {

  caller_env <- parent.frame()

  staged_name <- paste0(
    ".semiOPBART_external_reference_",
    reference_name
  )

  if (!exists(
      staged_name,
      envir = caller_env,
      inherits = FALSE
  )) {
    stop(
      "No staged external reference found for '",
      reference_name,
      "'"
    )
  }

  ref <- get(
    staged_name,
    envir = caller_env,
    inherits = FALSE
  )

  # ---------------------------------------------------------------
  # Validate tree chunks
  # ---------------------------------------------------------------

  if (
    is.na(ref$n_tree_chunks) ||
    ref$n_tree_chunks < 1L
  ) {
    stop("No tree chunks were staged")
  }

  tree_missing <- vapply(
    ref$tree_chunks[seq_len(ref$n_tree_chunks)],
    function(x) {
      is.null(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(x)
    },
    logical(1)
  )

  if (any(tree_missing)) {
    stop(
      "Missing/invalid tree chunk(s): ",
      paste(
        which(tree_missing),
        collapse = ", "
      )
    )
  }

  # ---------------------------------------------------------------
  # Validate ECDF chunks
  # ---------------------------------------------------------------

  if (
    is.na(ref$n_ecdf_chunks) ||
    ref$n_ecdf_chunks < 1L
  ) {
    stop("No ECDF chunks were staged")
  }

  ecdf_missing <- vapply(
    ref$ecdf_chunks[seq_len(ref$n_ecdf_chunks)],
    function(x) {
      is.null(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(x)
    },
    logical(1)
  )

  if (any(ecdf_missing)) {
    stop(
      "Missing/invalid ECDF chunk(s): ",
      paste(
        which(ecdf_missing),
        collapse = ", "
      )
    )
  }

  # ---------------------------------------------------------------
  # Reconstruct EXACT strings
  # ---------------------------------------------------------------

  trees_Serialize <- paste0(
    unlist(
      ref$tree_chunks[seq_len(ref$n_tree_chunks)],
      use.names = FALSE
    ),
    collapse = ""
  )

  ecdf_reference_Serialize <- paste0(
    unlist(
      ref$ecdf_chunks[seq_len(ref$n_ecdf_chunks)],
      use.names = FALSE
    ),
    collapse = ""
  )

  # ---------------------------------------------------------------
  # Validate reconstructed serialization
  # ---------------------------------------------------------------

  if (!nzchar(trees_Serialize)) {
    stop("Reconstructed tree serialization is empty")
  }

  if (!nzchar(ecdf_reference_Serialize)) {
    stop("Reconstructed ECDF serialization is empty")
  }

  if (!startsWith(trees_Serialize, "580a")) {
    stop(
      "Reconstructed tree serialization has invalid header: ",
      substr(trees_Serialize, 1L, 32L)
    )
  }

  if (!startsWith(ecdf_reference_Serialize, "580a")) {
    stop(
      "Reconstructed ECDF serialization has invalid header: ",
      substr(ecdf_reference_Serialize, 1L, 32L)
    )
  }

  # ---------------------------------------------------------------
  # Save final server-side reference
  # ---------------------------------------------------------------

  final_state_name <- paste0(
    ".semiOPBART_external_reference_final_",
    reference_name
  )

  base::assign(
    final_state_name,
    list(
      trees_Serialize = trees_Serialize,
      ecdf_reference_Serialize =
        ecdf_reference_Serialize
    ),
    envir = caller_env
  )

  # Remove temporary staging object.
  rm(
    list = staged_name,
    envir = caller_env
  )

  message(
    "[semiOPBARTLocalAssembleExternalReferenceDS] ",
    "assembled forest length=",
    nchar(trees_Serialize),
    ", ECDF length=",
    nchar(ecdf_reference_Serialize)
  )

  list(
    ok = TRUE,
    state_name = final_state_name,
    trees_length = nchar(trees_Serialize),
    ecdf_length = nchar(ecdf_reference_Serialize)
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
#'   semiOPBARTLocalPredictDS() -- final_trees_Serialize here is ONE reference site's own forest
#'   (from semiOPBARTLocalExportForestDS() above), not a pooled reconstruction; this function doesn't need to know or care which, it just loads whatever tree Serialize it's handed
#' @param data.name_test  a PREPARED (list(X, W, Y, ...)) object at this site, on the SAME normalized scale as the trees were grown on --
#'   i.e. built using the same normalization reference, whether that happened via the original shared ds.semiOPBARTPrepare() call or via ds.semiOPBARTPredictExternal()'s own re-prepare step
#' @param num_tree, k  must match what training used, for sigma_mu scaling to be consistent with the shipped tree ensemble
#' @export
semiOPBARTLocalPredictExternalStoredDS <- function(
    theta_draws_Serialize,
    us_draws_Serialize,
    reference_state_name,
    data.name_test,
    normalization_used = c("local_ecdf", "federated_ecdf"),
    num_tree,
    k = 1,
    newobj_pred = "semiOPBART_pred",
    levels_used_serialized,
    nfilter = 5,
    seed = 35
) {
  set.seed(seed)
  normalization_used <- match.arg(normalization_used)
  
  caller_env <- parent.frame()

  if (!exists(
      reference_state_name,
      envir = caller_env,
      inherits = FALSE
  )) {
    stop(
      "External reference state '",
      reference_state_name,
      "' not found"
    )
  }

  ref <- get(
    reference_state_name,
    envir = caller_env,
    inherits = FALSE
  )

  if (
    is.null(ref$trees_Serialize) ||
    is.null(ref$ecdf_reference_Serialize)
  ) {
    stop("External reference is incomplete")
  }

  semiOPBARTLocalPredictExternalDS(
    theta_draws_Serialize =
      theta_draws_Serialize,

    us_draws_Serialize =
      us_draws_Serialize,

    final_trees_Serialize =
      ref$trees_Serialize,

    ecdf_reference_Serialize =
      ref$ecdf_reference_Serialize,
    normalization_used=normalization_used,

    data.name_test =
      data.name_test,

    num_tree =
      num_tree,

    k =
      k,

    newobj_pred =
      newobj_pred,

    levels_used_serialized =
      levels_used_serialized,

    nfilter =
      nfilter,

    seed =
      seed,
    caller_env =
      caller_env
  )
}

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
  #   s$norm_info$federated_ecdf
  #
  # We export whichever reference actually exists.
  #
  norm_info <- s$norm_info

local_ecdf <- NULL
federated_ecdf <- NULL

if (!is.null(norm_info)) {
  if (identical(norm_info$method, "local_ecdf") && !is.null(norm_info$ecdfs)) {
    local_ecdf <- norm_info$ecdfs
  }
  if (identical(norm_info$method, "federated_ecdf") && !is.null(norm_info$federated_ecdf)) {
    federated_ecdf <- norm_info$federated_ecdf
  }
}

  if (is.null(local_ecdf) && is.null(federated_ecdf)) {
    stop(
      "semiOPBARTLocalExportForestDS(): neither local_ecdf nor ",
      "federated_ecdf is available in this site's normalization state. ",
      "Cannot export an external prediction reference."
    )
  }
  
  # ---------------------------------------------------------------
  # Forest
  # ---------------------------------------------------------------

  trees_Serialize <- semiOPBART_treesToSerialize(s$owned_forest$get_trees(seq_len(s$num_tree) - 1))

  list(
    trees_Serialize = trees_Serialize,
    num_tree = s$num_tree,

    # Ship normalization alongside the forest.
    local_ecdf_Serialize =
      if (!is.null(local_ecdf))
        semiOPBART_toSerialize(local_ecdf)
      else
        NULL,

    federated_ecdf_Serialize =
      if (!is.null(federated_ecdf))
        semiOPBART_toSerialize(federated_ecdf)
      else
        NULL,

    has_local_ecdf = !is.null(local_ecdf),
    has_federated_ecdf = !is.null(federated_ecdf)
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
    normalization_used = c("local_ecdf", "federated_ecdf"),
    data.name_test,
    num_tree,
    k = 1,
    newobj_pred = "semiOPBART_pred",
    levels_used_serialized = NULL,
    nfilter = 5,
    seed = 35,
    caller_env = parent.frame()
) {
    print("semiOPBARTLocalPredictExternalDS(): starting external prediction")

  set.seed(seed)
  normalization_used <- match.arg(normalization_used)
  s_test <- get(
    data.name_test,
    envir = caller_env,
    inherits = FALSE
  )
 
  if (is.null(levels_used_serialized)) {
    stop("semiOPBARTLocalPredictExternalDS(): levels_used_serialized is NULL")
  }
  levels_used <- semiOPBART_fromSerialize(levels_used_serialized)  
  theta_draws <-
    semiOPBART_fromSerialize(theta_draws_Serialize)

  us_draws <-
    semiOPBART_fromSerialize(us_draws_Serialize)

  ecdf_reference <-
    semiOPBART_fromSerialize(ecdf_reference_Serialize)

# ---------------------------------------------------------------
# Normalize posterior parameter shapes
#
# A single posterior draw may arrive as a vector.
# Internally we always use:
#
#   theta_draws: n_draws x n_theta
#   us_draws:    n_draws x n_thresholds
# ---------------------------------------------------------------

if (is.null(dim(theta_draws))) {

  if (!is.numeric(theta_draws) || length(theta_draws) == 0L) {
    stop("theta_draws must be a non-empty numeric vector or matrix")
  }

  theta_draws <- matrix(
    theta_draws,
    nrow = 1L,
    ncol = length(theta_draws)
  )
}

if (is.null(dim(us_draws))) {

  if (!is.numeric(us_draws) || length(us_draws) == 0L) {
    stop("us_draws must be a non-empty numeric vector or matrix")
  }

  us_draws <- matrix(
    us_draws,
    nrow = 1L,
    ncol = length(us_draws)
  )
}

cat(
  "[PREDICT] theta_draws dimensions = ",
  paste(dim(theta_draws), collapse = " x "),
  "\n",
  sep = ""
)

cat(
  "[PREDICT] us_draws dimensions = ",
  paste(dim(us_draws), collapse = " x "),
  "\n",
  sep = ""
)
  # ---------------------------------------------------------------
  # Apply exported ECDF to X
  # ---------------------------------------------------------------

expected_cols <- if ((normalization_used == "local_ecdf")) {
  names(ecdf_reference)
} else {
  ecdf_reference$cols
}
actual_cols <- colnames(s_test$X)

cat("[PREDICT] ECDF columns =", length(expected_cols), "\n")
cat("[PREDICT] test X columns =", length(actual_cols), "\n")

if (!identical(actual_cols, expected_cols)) {

    missing_cols <- setdiff(expected_cols, actual_cols)
    extra_cols <- setdiff(actual_cols, expected_cols)

    stop(
        "Federated ECDF/X column mismatch. ",
        "Expected: ", paste(expected_cols, collapse = ", "),
        "; actual: ", paste(actual_cols, collapse = ", "),
        "; missing: ",
        if (length(missing_cols))
            paste(missing_cols, collapse = ", ")
        else
            "NONE",
        "; extra: ",
        if (length(extra_cols))
            paste(extra_cols, collapse = ", ")
        else
            "NONE"
    )
}
cat(
    "[PREDICT] ECDF columns match test X: ",
    length(expected_cols),
    " columns\n",
    sep = ""
)

if (normalization_used == "federated_ecdf") {
  print("semiOPBARTLocalPredictExternalDS(): using federated ECDF normalization")
if (!is.list(ecdf_reference) ||
    is.null(ecdf_reference$cols) ||
    is.null(ecdf_reference$bin_edges) ||
    is.null(ecdf_reference$cum_frac)) {
    stop(
        "Invalid serialized global ECDF reference: ",
        "expected cols, bin_edges and cum_frac"
    )
}
    ecdf_fns <- setNames(lapply(ecdf_reference$cols, function(cn) {
      edges <- ecdf_reference$bin_edges[[cn]]; cf <- ecdf_reference$cum_frac[[cn]]
      function(y) pmin(pmax(stats::approx(edges, cf, xout = y, method = "linear",
                                           rule = 2)$y, 0), 1)
    }), ecdf_reference$cols)

  for (j in seq_len(ncol(s_test$X))) {
        s_test$X[, j] <- ecdf_fns[[j]](s_test$X[, j])
        }

}else if (normalization_used == "local_ecdf") {
    print("semiOPBARTLocalPredictExternalDS(): using local ECDF normalization")
  for (j in seq_len(ncol(s_test$X))) {
    s_test$X[, j] <-
      ecdf_reference[[j]](
        s_test$X[, j]
      )
  }
}else {
    stop("semiOPBARTLocalPredictExternalDS(): unknown normalization_used: ", normalization_used)
}
  print("semiOPBARTLocalPredictExternalDS(): applied ECDF reference to test X")


  # ---------------------------------------------------------------
  # Reconstruct reference forest
  # ---------------------------------------------------------------

  ref_forest <- tryCatch(

    semiOPBART_treesFromSerialize_check(
      final_trees_Serialize
    ),

    error = function(e) {

      stop(
        "semiOPBARTLocalPredictExternalDS(): ",
        "semiOPBART_treesFromSerialize() failed: ",
        conditionMessage(e)
      )

    }

  )

  if (is.null(ref_forest)) {

    stop(
      "semiOPBARTLocalPredictExternalDS(): ",
      "semiOPBART_treesFromSerialize() returned NULL. ",
      "The reconstructed serialized forest is syntactically present ",
      "but could not be converted back into a forest."
    )

  }

  print("semiOPBARTLocalPredictExternalDS(): reconstructed reference forest successfully")

hypers_template <- tryCatch(
   SoftBart::Hypers(
    X = s_test$X,
    Y = rep(1, nrow(s_test$X)),
    sigma_hat = 1,
    normalize_Y = FALSE
),
    error = function(e) {
        stop(
            "SoftBart::Hypers() failed: ",
            conditionMessage(e)
        )
    }
)

  hypers_template$sigma_mu <- 3 / k / sqrt(num_tree)
  hypers_template$sigma     <- 1
  hypers_template$sigma_hat <- 1
  hypers_template$num_tree  <- num_tree
  # same confirmed fix as semiOPBARTLocalInitDS() (semiopbartDS.R):
  # without this, multi-level categorical x_features would be mishandled
  # by the Dirichlet variable-selection prior
  hypers_template$group <- dummy_assign(s_test$dv)
  opts <- SoftBart::Opts(); opts$update_sigma <- FALSE
  eval_forest <- SoftBart::MakeForest(hypers_template, opts, FALSE)



  semiOPBART_treesFromSerialize(eval_forest, final_trees_Serialize, seq_len(num_tree) - 1)




print("semiOPBARTLocalPredictExternalDS(): reconstructed reference forest and set trees")
if (!is.function(eval_forest$do_predict) &&
    !methods::hasMethod("do_predict", class(eval_forest))) {
  stop("Reconstructed external forest has no usable do_predict()")
}
 
  fx_new <- tryCatch(
  {
    print("semiOPBARTLocalPredictExternalDS(): calling do_predict()")
    out <- eval_forest$do_predict(s_test$X)
    print("semiOPBARTLocalPredictExternalDS(): do_predict() returned")
    as.numeric(out)
  },
  error = function(e) {
    stop(
      "semiOPBARTLocalPredictExternalDS(): do_predict() failed: ",
      conditionMessage(e)
    )
  }
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
  print("semiOPBARTLocalPredictExternalDS(): computed posterior probabilities and MAP scores")

  base::assign(
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
    n_test = n_test, prob = prob_matrix,
      map = map_score,
      truth = s_test$Y, used_external_shipping = TRUE
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

#   base::assign(newobj_pred, list(prob = prob_matrix, map = map_score, truth = s_test$Y,
#                             used_external_shipping = TRUE),
#          envir = parent.frame())
#   list(n_test = n_test)
# }



