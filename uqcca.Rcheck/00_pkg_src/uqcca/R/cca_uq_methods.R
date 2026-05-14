`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.uqcca_cache <- new.env(parent = emptyenv())

option_or_env <- function(option_name, env_name, default = NULL) {
  value <- getOption(option_name, default = NULL)
  if (is.null(value) || !length(value) || !nzchar(as.character(value)[1])) {
    value <- Sys.getenv(env_name, unset = "")
  }

  if (!length(value) || !nzchar(as.character(value)[1])) {
    return(default)
  }

  as.character(value)[1]
}

default_ccar3_path <- function() {
  option_or_env("uqcca.ccar3_path", "UQCCA_CCAR3_PATH", default = NULL)
}

default_ccar3_code_path <- function() {
  option_or_env("uqcca.ccar3_code_path", "UQCCA_CCAR3_CODE_PATH", default = NULL)
}

ensure_feature_names <- function(x, prefix) {
  if (is.null(colnames(x))) {
    colnames(x) <- paste0(prefix, seq_len(ncol(x)))
  }
  x
}

as_numeric_matrix <- function(x, name) {
  if (is.null(dim(x))) {
    stop(sprintf("`%s` must be a matrix-like object.", name), call. = FALSE)
  }

  mat <- data.matrix(x)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  if (anyNA(mat)) {
    stop(sprintf("`%s` contains missing values after coercion to a numeric matrix.", name),
      call. = FALSE
    )
  }

  mat
}

compute_preprocess_params <- function(x, mode = c("scale", "center", "none")) {
  mode <- match.arg(mode)
  x <- as_numeric_matrix(x, "x")

  center <- if (mode %in% c("scale", "center")) colMeans(x) else rep(0, ncol(x))
  scale_values <- if (mode == "scale") apply(x, 2, stats::sd) else rep(1, ncol(x))
  scale_values[!is.finite(scale_values) | scale_values == 0] <- 1

  names(center) <- colnames(x)
  names(scale_values) <- colnames(x)

  list(mode = mode, center = as.numeric(center), scale = as.numeric(scale_values))
}

apply_preprocess <- function(x, params, name) {
  x <- as_numeric_matrix(x, name)

  if (params$mode %in% c("scale", "center")) {
    x <- sweep(x, 2, params$center, "-")
  }
  if (identical(params$mode, "scale")) {
    x <- sweep(x, 2, params$scale, "/")
  }

  x[!is.finite(x)] <- 0
  x
}

prepare_xy <- function(X, Y, mode = c("scale", "center", "none")) {
  mode <- match.arg(mode)
  x_params <- compute_preprocess_params(X, mode = mode)
  y_params <- compute_preprocess_params(Y, mode = mode)

  list(
    X = apply_preprocess(X, x_params, "X"),
    Y = apply_preprocess(Y, y_params, "Y"),
    x_params = x_params,
    y_params = y_params,
    mode = mode
  )
}

apply_xy_preprocess <- function(X, Y, preprocess) {
  list(
    X = apply_preprocess(X, preprocess$x_params, "X"),
    Y = apply_preprocess(Y, preprocess$y_params, "Y")
  )
}

extract_alcohol_blocks <- function(alcohol) {
  if (!is.null(alcohol$gene) && !is.null(alcohol$meth)) {
    return(list(X = alcohol$gene, Y = alcohol$meth, disorder = alcohol$disorder))
  }

  if (!is.null(alcohol$Xlist) && length(alcohol$Xlist) >= 2) {
    return(list(X = alcohol$Xlist[[1]], Y = alcohol$Xlist[[2]], disorder = alcohol$Y))
  }

  if (!is.null(alcohol$X1) && !is.null(alcohol$X2)) {
    return(list(X = alcohol$X1, Y = alcohol$X2, disorder = alcohol$Y))
  }

  if (length(alcohol) >= 2) {
    disorder <- if (length(alcohol) >= 3) alcohol[[3]] else NULL
    return(list(X = alcohol[[1]], Y = alcohol[[2]], disorder = disorder))
  }

  stop("Could not determine how to extract the alcohol example matrices.", call. = FALSE)
}

#' Load the alcohol multi-omics example.
#'
#' @param data_path Optional path to an `.rda` file containing `alcohol`.
#' @return A list with raw `X`, raw `Y`, the `disorder` labels, and metadata.
#' @export
load_alcohol_example <- function(data_path = NULL) {
  if (!is.null(data_path)) {
    alcohol_env <- new.env(parent = emptyenv())
    load(data_path, envir = alcohol_env)
    if (!exists("alcohol", envir = alcohol_env, inherits = FALSE)) {
      stop("`data_path` did not create an object called `alcohol`.", call. = FALSE)
    }
    alcohol <- get("alcohol", envir = alcohol_env, inherits = FALSE)
    source_name <- normalizePath(data_path)
  } else {
    if (!requireNamespace("CVR", quietly = TRUE)) {
      stop(
        "Package 'CVR' is required to load the alcohol example, or provide `data_path`.",
        call. = FALSE
      )
    }
    data("alcohol", package = "CVR", envir = environment())
    source_name <- "CVR::alcohol"
  }

  blocks <- extract_alcohol_blocks(alcohol)
  X <- ensure_feature_names(as_numeric_matrix(blocks$X, "alcohol X"), "gene_")
  Y <- ensure_feature_names(as_numeric_matrix(blocks$Y, "alcohol Y"), "meth_")

  list(
    X = X,
    Y = Y,
    disorder = blocks$disorder,
    n = nrow(X),
    p = ncol(X),
    q = ncol(Y),
    source = source_name
  )
}

get_ccar3_api <- function(ccar3_path = NULL,
                          prefer_source = TRUE,
                          quiet = TRUE) {
  ccar3_path <- ccar3_path %||% default_ccar3_path()
  normalized_path <- if (!is.null(ccar3_path) && nzchar(ccar3_path)) {
    normalizePath(ccar3_path, mustWork = FALSE)
  } else {
    ""
  }
  cache_key <- paste(normalized_path, prefer_source, sep = "::")

  if (exists(cache_key, envir = .uqcca_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .uqcca_cache, inherits = FALSE))
  }

  api <- NULL

  if (prefer_source && nzchar(normalized_path) && dir.exists(normalized_path)) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      warning(
        "Package 'pkgload' is not available; falling back to the installed 'ccar3' package.",
        call. = FALSE
      )
    } else {
      api <- tryCatch(
        pkgload::load_all(
          normalized_path,
          quiet = quiet,
          export_all = FALSE,
          helpers = FALSE,
          attach_testthat = FALSE
        ),
        error = function(e) NULL
      )
      if (is.list(api) && !is.null(api$env)) {
        api <- api$env
      }
    }
  }

  if (is.null(api)) {
    if (!requireNamespace("ccar3", quietly = TRUE)) {
      stop(
        paste(
          "Could not load 'ccar3'.",
          "Install the package, set `options(uqcca.ccar3_path = \"/path/to/ccar3\")`,",
          "or define the `UQCCA_CCAR3_PATH` environment variable."
        ),
        call. = FALSE
      )
    }
    api <- asNamespace("ccar3")
  }

  assign(cache_key, api, envir = .uqcca_cache)
  api
}

call_ccar3_fn <- function(fun, args) {
  args <- args[!vapply(args, is.null, logical(1))]
  allowed <- intersect(names(args), names(formals(fun)))
  do.call(fun, args[allowed])
}

get_ccar3_fn <- function(api, name) {
  get(name, envir = api, inherits = TRUE)
}

new_cca_fit <- function(fit, method, preprocess, settings, x_names, y_names) {
  cor_values <- fit$cor %||% fit$Lambda %||% numeric()
  lambda_value <- fit$lambda %||% fit$lambda_x %||% fit$lambda.opt %||% settings$lambda
  B_value <- fit$B %||% fit$B_opt %||% fit$Bhat
  cv_summary <- fit$cv_summary %||% fit$resultsx %||% fit$cv.scores

  structure(
    list(
      method = method,
      U = fit$U,
      V = fit$V,
      cor = as.numeric(cor_values),
      lambda = lambda_value,
      B = B_value,
      cv_summary = cv_summary,
      preprocess = preprocess,
      settings = settings,
      x_names = x_names,
      y_names = y_names,
      fit = fit
    ),
    class = "cca_uq_fit"
  )
}

effective_ecca_preprocess_mode <- function(mode) {
  # ecca always centers internally, so "none" can only mean "no extra scaling".
  if (identical(mode, "scale")) "scale" else "center"
}

build_ecca_preprocess <- function(X, Y, preprocess_mode) {
  effective_mode <- effective_ecca_preprocess_mode(preprocess_mode)

  list(
    mode = effective_mode,
    x_params = compute_preprocess_params(X, mode = effective_mode),
    y_params = compute_preprocess_params(Y, mode = effective_mode)
  )
}

#' Fit sparse CCA with cross-validated `ccar3::ecca`
#'
#' Fits sparse canonical directions on two paired data blocks and stores the
#' selected penalty, canonical correlations, preprocessing recipe, and feature
#' names needed for downstream uncertainty quantification.
#'
#' @param X,Y Numeric matrices with rows representing the shared samples and
#'   columns representing features from the two views.
#' @param r Target rank, i.e. the number of canonical components to estimate.
#' @param lambdas Candidate penalty values passed to `ccar3::ecca.cv`.
#' @param kfolds Number of folds used by cross-validation when tuning
#'   `lambda`.
#' @param preprocess_mode Preprocessing applied before fitting. `"scale"`
#'   centers and scales each column, `"center"` only centers, and `"none"`
#'   skips any extra scaling. `ecca` still centers internally.
#' @param ccar3_api Optional namespace or environment containing the `ccar3`
#'   fitting functions. When `NULL`, the package looks for an installed
#'   `ccar3` package or a local source checkout.
#' @param ccar3_path Optional path to a local `ccar3` source checkout.
#' @param prefer_source Whether to prefer loading `ccar3_path` with
#'   `pkgload::load_all()` before falling back to an installed `ccar3`
#'   package.
#' @param solver,LW_Sy,thresh_0,matrix_free_threshold,cg_tol,cg_maxiter
#'   Compatibility arguments retained in the stored settings list.
#' @param parallelize Whether to request parallel cross-validation from the
#'   underlying `ccar3` routine.
#' @param rho ADMM penalty parameter forwarded to `ccar3`.
#' @param niter Maximum number of iterations for the sparse CCA solver.
#' @param thresh Convergence tolerance passed to `ccar3`.
#' @param verbose Whether to print progress from the underlying fit.
#' @param nb_cores Optional number of worker cores for cross-validation.
#' @param select Cross-validation selection rule such as `"lambda.min"` or
#'   `"lambda.1se"` when supported by the local `ccar3`.
#' @param maxiter_cv Optional iteration cap for the inner CV fits.
#' @param set_seed_cv Optional seed forwarded to the `ccar3` CV routine.
#' @param scoring_method Cross-validation loss used by `ecca.cv`.
#' @param cv_use_median Whether to aggregate fold losses with the median
#'   instead of the mean when supported by `ccar3`.
#' @param epsilon_sv Small singular-value guard passed to `ccar3`.
#' @param ridge_whiten Ridge regularization used during whitening for
#'   numerical stability.
#'
#' @return A `cca_uq_fit` list containing the fitted loadings, canonical
#'   correlations, chosen penalty, preprocessing recipe, and bookkeeping for
#'   downstream bootstrap and cross-fitted UQ.
#' @export
fit_ecca_cv <- function(X, Y,
                        r = 2,
                        lambdas = 10^seq(-2, 0, length.out = 8),
                        kfolds = 5,
                        preprocess_mode = c("scale", "center", "none"),
                        ccar3_api = NULL,
                        ccar3_path = NULL,
                        prefer_source = TRUE,
                        solver = "ADMM",
                        parallelize = FALSE,
                        LW_Sy = TRUE,
                        rho = 1,
                        niter = 5000,
                        thresh = 1e-4,
                        thresh_0 = 1e-6,
                        matrix_free_threshold = 4000L,
                        cg_tol = 1e-6,
                        cg_maxiter = NULL,
                        verbose = FALSE,
                        nb_cores = NULL,
                        select = "lambda.min",
                        maxiter_cv = NULL,
                        set_seed_cv = NULL,
                        scoring_method = c("mse", "trace"),
                        cv_use_median = FALSE,
                        epsilon_sv = 1e-8,
                        ridge_whiten = 1e-8) {
  preprocess_mode <- match.arg(preprocess_mode)
  scoring_method <- match.arg(scoring_method)
  X <- ensure_feature_names(as_numeric_matrix(X, "X"), "X_")
  Y <- ensure_feature_names(as_numeric_matrix(Y, "Y"), "Y_")
  api <- ccar3_api %||% get_ccar3_api(
    ccar3_path = ccar3_path,
    prefer_source = prefer_source,
    quiet = !isTRUE(verbose)
  )
  effective_mode <- effective_ecca_preprocess_mode(preprocess_mode)
  preprocess <- build_ecca_preprocess(X, Y, preprocess_mode)
  maxiter_cv <- maxiter_cv %||% min(as.integer(niter), 300L)

  fit <- call_ccar3_fn(
    get_ccar3_fn(api, "ecca.cv"),
    list(
      X = X,
      Y = Y,
      r = r,
      lambdas = lambdas,
      nfold = kfolds,
      select = select,
      standardize = identical(effective_mode, "scale"),
      rho = rho,
      eps = thresh,
      maxiter = niter,
      maxiter_cv = maxiter_cv,
      parallel = parallelize,
      set_seed_cv = set_seed_cv,
      scoring_method = scoring_method,
      cv_use_median = cv_use_median,
      epsilon_sv = epsilon_sv,
      ridge_whiten = ridge_whiten,
      verbose = verbose,
      nb_cores = nb_cores
    )
  )

  settings <- list(
    fitter = "ecca",
    r = r,
    lambdas = lambdas,
    kfolds = kfolds,
    parallelize = parallelize,
    rho = rho,
    niter = niter,
    thresh = thresh,
    nb_cores = nb_cores,
    preprocess_mode = effective_mode,
    requested_preprocess_mode = preprocess_mode,
    select = select,
    maxiter_cv = maxiter_cv,
    set_seed_cv = set_seed_cv,
    scoring_method = scoring_method,
    cv_use_median = cv_use_median,
    epsilon_sv = epsilon_sv,
    ridge_whiten = ridge_whiten,
    solver = solver,
    LW_Sy = LW_Sy,
    thresh_0 = thresh_0,
    matrix_free_threshold = matrix_free_threshold,
    cg_tol = cg_tol,
    cg_maxiter = cg_maxiter,
    ccar3_path = ccar3_path,
    prefer_source = prefer_source
  )

  new_cca_fit(
    fit = fit,
    method = "ecca.cv",
    preprocess = preprocess,
    settings = settings,
    x_names = colnames(X),
    y_names = colnames(Y)
  )
}

#' Fit sparse CCA with a fixed `ccar3::ecca` penalty
#'
#' @inheritParams fit_ecca_cv
#' @param lambda Single penalty value passed to `ccar3::ecca`.
#'
#' @return A `cca_uq_fit` object using the supplied `lambda`.
#' @export
fit_ecca_fixed_lambda <- function(X, Y,
                                  lambda,
                                  r = 2,
                                  preprocess_mode = c("scale", "center", "none"),
                                  ccar3_api = NULL,
                                  ccar3_path = NULL,
                                  prefer_source = TRUE,
                                  solver = "ADMM",
                                  LW_Sy = TRUE,
                                  rho = 1,
                                  niter = 10000,
                                  thresh = 1e-4,
                                  thresh_0 = 1e-6,
                                  matrix_free_threshold = 4000L,
                                  cg_tol = 1e-6,
                                  cg_maxiter = NULL,
                                  verbose = FALSE,
                                  epsilon_sv = 1e-8,
                                  ridge_whiten = 1e-8) {
  preprocess_mode <- match.arg(preprocess_mode)
  X <- ensure_feature_names(as_numeric_matrix(X, "X"), "X_")
  Y <- ensure_feature_names(as_numeric_matrix(Y, "Y"), "Y_")
  api <- ccar3_api %||% get_ccar3_api(
    ccar3_path = ccar3_path,
    prefer_source = prefer_source,
    quiet = !isTRUE(verbose)
  )
  effective_mode <- effective_ecca_preprocess_mode(preprocess_mode)
  preprocess <- build_ecca_preprocess(X, Y, preprocess_mode)

  fit <- call_ccar3_fn(
    get_ccar3_fn(api, "ecca"),
    list(
      X = X,
      Y = Y,
      lambda = lambda,
      r = r,
      standardize = identical(effective_mode, "scale"),
      rho = rho,
      maxiter = niter,
      eps = thresh,
      epsilon_sv = epsilon_sv,
      ridge_whiten = ridge_whiten,
      verbose = verbose
    )
  )

  settings <- list(
    fitter = "ecca",
    lambda = lambda,
    r = r,
    rho = rho,
    niter = niter,
    thresh = thresh,
    preprocess_mode = effective_mode,
    requested_preprocess_mode = preprocess_mode,
    epsilon_sv = epsilon_sv,
    ridge_whiten = ridge_whiten,
    solver = solver,
    LW_Sy = LW_Sy,
    thresh_0 = thresh_0,
    matrix_free_threshold = matrix_free_threshold,
    cg_tol = cg_tol,
    cg_maxiter = cg_maxiter,
    ccar3_path = ccar3_path,
    prefer_source = prefer_source
  )

  new_cca_fit(
    fit = fit,
    method = "ecca",
    preprocess = preprocess,
    settings = settings,
    x_names = colnames(X),
    y_names = colnames(Y)
  )
}

#' Backward-compatible alias for [fit_ecca_cv()].
#'
#' @param ... Arguments passed to [fit_ecca_cv()].
#' @return A `cca_uq_fit` object returned by [fit_ecca_cv()].
#' @export
fit_cca_rrr_cv <- function(...) {
  fit_ecca_cv(...)
}

#' Backward-compatible alias for [fit_ecca_fixed_lambda()].
#'
#' @param ... Arguments passed to [fit_ecca_fixed_lambda()].
#' @return A `cca_uq_fit` object returned by [fit_ecca_fixed_lambda()].
#' @export
fit_cca_rrr_fixed_lambda <- function(...) {
  fit_ecca_fixed_lambda(...)
}

canonical_scores <- function(X, Y, fit) {
  transformed <- apply_xy_preprocess(X, Y, fit$preprocess)
  list(
    X = transformed$X %*% fit$U,
    Y = transformed$Y %*% fit$V,
    X_preprocessed = transformed$X,
    Y_preprocessed = transformed$Y
  )
}

make_identity_preprocess <- function(X, Y) {
  make_params <- function(feature_names) {
    n_features <- length(feature_names)
    center <- rep(0, n_features)
    scale <- rep(1, n_features)
    names(center) <- feature_names
    names(scale) <- feature_names
    list(mode = "none", center = center, scale = scale)
  }

  list(
    mode = "none",
    x_params = make_params(colnames(X)),
    y_params = make_params(colnames(Y))
  )
}

preprocess_once_for_crossfit <- function(X, Y, reference_fit) {
  transformed <- apply_xy_preprocess(X, Y, reference_fit$preprocess)

  diagnostics <- data.frame(
    block = c("X", "Y"),
    preprocess_mode = reference_fit$preprocess$mode,
    mean_abs_colmean = c(
      mean(abs(colMeans(transformed$X))),
      mean(abs(colMeans(transformed$Y)))
    ),
    max_abs_colmean = c(
      max(abs(colMeans(transformed$X))),
      max(abs(colMeans(transformed$Y)))
    ),
    stringsAsFactors = FALSE
  )

  list(
    X = transformed$X,
    Y = transformed$Y,
    preprocess = make_identity_preprocess(transformed$X, transformed$Y),
    diagnostics = diagnostics
  )
}

compute_joint_rotation <- function(U, V, U_ref, V_ref) {
  target <- crossprod(U, U_ref) + crossprod(V, V_ref)
  if (all(abs(target) < 1e-12)) {
    return(diag(ncol(U)))
  }

  svd_target <- svd(target)
  svd_target$u %*% t(svd_target$v)
}

component_similarity_matrix <- function(U, V, U_ref, V_ref, absolute = FALSE) {
  x_similarity <- crossprod(U, U_ref)
  y_similarity <- crossprod(V, V_ref)
  if (absolute) {
    return(abs(x_similarity) + abs(y_similarity))
  }
  x_similarity + y_similarity
}

enumerate_permutations <- function(index) {
  if (length(index) <= 1L) {
    return(list(index))
  }

  output <- vector("list", 0L)
  for (i in seq_along(index)) {
    remaining <- index[-i]
    tails <- enumerate_permutations(remaining)
    output <- c(
      output,
      lapply(tails, function(tail) c(index[i], tail))
    )
  }
  output
}

best_component_permutation <- function(similarity) {
  n_components <- nrow(similarity)
  index <- seq_len(n_components)

  if (n_components <= 1L) {
    return(index)
  }

  if (n_components <= 7L) {
    candidates <- enumerate_permutations(index)
    scores <- vapply(
      candidates,
      function(candidate) sum(similarity[cbind(candidate, index)]),
      numeric(1)
    )
    return(candidates[[which.max(scores)]])
  }

  chosen <- integer(n_components)
  available <- rep(TRUE, n_components)
  for (component in index) {
    scores <- similarity[, component]
    scores[!available] <- -Inf
    match_idx <- which.max(scores)
    chosen[component] <- match_idx
    available[match_idx] <- FALSE
  }

  chosen
}

align_fit_to_reference <- function(fit, reference_fit) {
  pre_similarity <- component_similarity_matrix(
    U = fit$U,
    V = fit$V,
    U_ref = reference_fit$U,
    V_ref = reference_fit$V
  )

  rotation <- compute_joint_rotation(
    U = fit$U,
    V = fit$V,
    U_ref = reference_fit$U,
    V_ref = reference_fit$V
  )

  fit$U <- fit$U %*% rotation
  fit$V <- fit$V %*% rotation
  rotated_similarity <- component_similarity_matrix(
    U = fit$U,
    V = fit$V,
    U_ref = reference_fit$U,
    V_ref = reference_fit$V,
    absolute = TRUE
  )

  permutation <- best_component_permutation(rotated_similarity)
  fit$U <- fit$U[, permutation, drop = FALSE]
  fit$V <- fit$V[, permutation, drop = FALSE]
  if (!is.null(fit$cor) && length(fit$cor) >= length(permutation)) {
    fit$cor <- as.numeric(fit$cor[permutation])
  }

  sign_correction <- sign(
    colSums(fit$U * reference_fit$U) + colSums(fit$V * reference_fit$V)
  )
  sign_correction[!is.finite(sign_correction) | sign_correction == 0] <- 1

  fit$U <- sweep(fit$U, 2, sign_correction, "*")
  fit$V <- sweep(fit$V, 2, sign_correction, "*")
  if (!is.null(fit$cor)) {
    fit$cor <- abs(fit$cor)
  }

  if (is.list(fit$fit)) {
    fit$fit$U <- fit$U
    fit$fit$V <- fit$V
    if (!is.null(fit$fit$cor)) {
      fit$fit$cor <- fit$cor
    }
    if (!is.null(fit$fit$Lambda) && !is.null(fit$cor)) {
      fit$fit$Lambda <- fit$cor
    }
  }

  aligned_similarity <- component_similarity_matrix(
    U = fit$U,
    V = fit$V,
    U_ref = reference_fit$U,
    V_ref = reference_fit$V
  )

  fit$alignment <- list(
    rotation = rotation,
    permutation = permutation,
    sign = sign_correction,
    pre_similarity = pre_similarity,
    rotated_abs_similarity = rotated_similarity,
    aligned_similarity = aligned_similarity,
    aligned_abs_similarity = abs(aligned_similarity)
  )

  fit
}

refit_from_reference <- function(X, Y,
                                 reference_fit,
                                 mode = c("fixed_lambda", "retune_cv"),
                                 ccar3_api = NULL,
                                 overrides = list()) {
  mode <- match.arg(mode)
  base <- reference_fit$settings
  fitter <- base$fitter %||% if (reference_fit$method %in% c("ecca", "ecca.cv")) "ecca" else "cca_rrr"

  if (mode == "fixed_lambda") {
    args <- modifyList(
      list(
        X = X,
        Y = Y,
        lambda = reference_fit$lambda,
        r = base$r,
        preprocess_mode = base$preprocess_mode,
        ccar3_api = ccar3_api,
        ccar3_path = base$ccar3_path,
        prefer_source = base$prefer_source,
        solver = base$solver,
        LW_Sy = base$LW_Sy,
        rho = base$rho,
        niter = base$niter,
        thresh = base$thresh,
        thresh_0 = base$thresh_0,
        matrix_free_threshold = base$matrix_free_threshold,
        cg_tol = base$cg_tol,
        cg_maxiter = base$cg_maxiter,
        epsilon_sv = base$epsilon_sv,
        ridge_whiten = base$ridge_whiten,
        verbose = FALSE
      ),
      overrides,
      keep.null = TRUE
    )
    if (identical(fitter, "ecca")) {
      return(do.call(fit_ecca_fixed_lambda, args))
    }
    return(do.call(fit_cca_rrr_fixed_lambda, args))
  }

  args <- modifyList(
    list(
      X = X,
      Y = Y,
      r = base$r,
      lambdas = base$lambdas,
      kfolds = base$kfolds,
      preprocess_mode = base$preprocess_mode,
      ccar3_api = ccar3_api,
      ccar3_path = base$ccar3_path,
      prefer_source = base$prefer_source,
      solver = base$solver,
      parallelize = FALSE,
      LW_Sy = base$LW_Sy,
      rho = base$rho,
      niter = base$niter,
      thresh = base$thresh,
      thresh_0 = base$thresh_0,
      matrix_free_threshold = base$matrix_free_threshold,
      cg_tol = base$cg_tol,
      cg_maxiter = base$cg_maxiter,
      select = base$select,
      maxiter_cv = base$maxiter_cv,
      set_seed_cv = base$set_seed_cv,
      scoring_method = base$scoring_method,
      cv_use_median = base$cv_use_median,
      epsilon_sv = base$epsilon_sv,
      ridge_whiten = base$ridge_whiten,
      verbose = FALSE,
      nb_cores = base$nb_cores
    ),
    overrides,
    keep.null = TRUE
  )

  if (identical(fitter, "ecca")) {
    return(do.call(fit_ecca_cv, args))
  }
  do.call(fit_cca_rrr_cv, args)
}

summarize_bootstrap_matrix <- function(boot_array, alpha = 0.05) {
  list(
    mean = apply(boot_array, c(1, 2), mean, na.rm = TRUE),
    sd = apply(boot_array, c(1, 2), stats::sd, na.rm = TRUE),
    lower = apply(boot_array, c(1, 2), stats::quantile, probs = alpha / 2, na.rm = TRUE),
    upper = apply(boot_array, c(1, 2), stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  )
}

summarize_bootstrap_vector <- function(boot_matrix, alpha = 0.05) {
  list(
    mean = colMeans(boot_matrix, na.rm = TRUE),
    sd = apply(boot_matrix, 2, stats::sd, na.rm = TRUE),
    lower = apply(boot_matrix, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE),
    upper = apply(boot_matrix, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  )
}

selection_frequency <- function(boot_array, threshold = 1e-6) {
  apply(abs(boot_array) > threshold, c(1, 2), mean, na.rm = TRUE)
}

sign_stability <- function(boot_array, reference, threshold = 1e-6) {
  aligned_sign <- sign(boot_array)
  aligned_sign[abs(boot_array) <= threshold | !is.finite(boot_array)] <- 0

  reference_sign <- sign(reference)
  reference_sign[abs(reference) <= threshold | !is.finite(reference)] <- 0
  reference_array <- array(
    reference_sign,
    dim = dim(boot_array)
  )

  apply(aligned_sign == reference_array, c(1, 2), mean, na.rm = TRUE)
}

orthonormal_basis <- function(x) {
  x <- as_numeric_matrix(x, "x")
  qr_x <- qr(x)
  rank_x <- qr_x$rank
  if (rank_x < 1L) {
    stop("`x` must span a non-empty subspace.", call. = FALSE)
  }
  qr.Q(qr_x)[, seq_len(rank_x), drop = FALSE]
}

compute_principal_angles <- function(estimate, reference) {
  q_estimate <- orthonormal_basis(estimate)
  q_reference <- orthonormal_basis(reference)
  overlap <- svd(crossprod(q_reference, q_estimate), nu = 0, nv = 0)$d
  overlap <- pmin(pmax(overlap, -1), 1)
  acos(overlap)
}

procrustes_align_matrix <- function(estimate, reference) {
  estimate <- as_numeric_matrix(estimate, "estimate")
  reference <- as_numeric_matrix(reference, "reference")

  target <- crossprod(estimate, reference)
  svd_target <- svd(target)
  rotation <- svd_target$u %*% t(svd_target$v)

  list(
    aligned = estimate %*% rotation,
    rotation = rotation
  )
}

compute_subspace_distance <- function(estimate, reference) {
  estimate <- as_numeric_matrix(estimate, "estimate")
  reference <- as_numeric_matrix(reference, "reference")

  if (nrow(estimate) != nrow(reference)) {
    stop("`estimate` and `reference` must have the same number of rows.", call. = FALSE)
  }

  estimate_qr <- qr(estimate)
  reference_qr <- qr(reference)
  estimate_rank <- estimate_qr$rank
  reference_rank <- reference_qr$rank

  if (estimate_rank < 1L || reference_rank < 1L) {
    stop("Both `estimate` and `reference` must span a non-empty subspace.", call. = FALSE)
  }

  estimate_basis <- qr.Q(estimate_qr)[, seq_len(estimate_rank), drop = FALSE]
  reference_basis <- qr.Q(reference_qr)[, seq_len(reference_rank), drop = FALSE]

  norm(
    tcrossprod(estimate_basis) - tcrossprod(reference_basis),
    type = "F"
  )
}

compute_subspace_distances <- function(boot_array, reference) {
  vapply(seq_len(dim(boot_array)[3]), function(index) {
    compute_subspace_distance(boot_array[, , index], reference)
  }, numeric(1))
}

#' Bootstrap uncertainty quantification for sparse CCA loadings
#'
#' Refits sparse CCA on bootstrap resamples, aligns each bootstrap solution to
#' the full-data fit, and stores aligned loading arrays together with canonical
#' correlations for later summarization.
#'
#' @param X,Y Numeric matrices with shared rows representing the paired
#'   samples.
#' @param reference_fit A fitted `cca_uq_fit` object, typically produced by
#'   [fit_ecca_cv()] or [fit_ecca_fixed_lambda()].
#' @param B Number of bootstrap replicates to request.
#' @param seed Random seed for bootstrap resampling.
#' @param refit_mode Whether each bootstrap replicate reuses the selected
#'   penalty (`"fixed_lambda"`) or re-runs cross-validation (`"retune_cv"`).
#' @param refit_args Optional named overrides passed to each bootstrap refit.
#' @param ccar3_api Optional loaded `ccar3` namespace or environment.
#' @param progress Whether to emit bootstrap progress messages.
#'
#' @return A list with the reference fit, aligned bootstrap arrays for `U` and
#'   `V`, aligned bootstrap canonical correlations, and replicate counts.
#' @export
bootstrap_cca_uq <- function(X, Y,
                             reference_fit,
                             B = 100,
                             seed = 1,
                             refit_mode = c("fixed_lambda", "retune_cv"),
                             refit_args = list(),
                             ccar3_api = NULL,
                             progress = interactive()) {
  refit_mode <- match.arg(refit_mode)
  X <- as_numeric_matrix(X, "X")
  Y <- as_numeric_matrix(Y, "Y")
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  r <- ncol(reference_fit$U)
  api <- ccar3_api %||% get_ccar3_api(
    ccar3_path = reference_fit$settings$ccar3_path,
    prefer_source = reference_fit$settings$prefer_source
  )

  set.seed(seed)
  U_array <- array(NA_real_, dim = c(p, r, B))
  V_array <- array(NA_real_, dim = c(q, r, B))
  cor_matrix <- matrix(NA_real_, nrow = B, ncol = r)

  for (b in seq_len(B)) {
    index <- sample.int(n, size = n, replace = TRUE)
    boot_fit <- tryCatch(
      refit_from_reference(
        X = X[index, , drop = FALSE],
        Y = Y[index, , drop = FALSE],
        reference_fit = reference_fit,
        mode = refit_mode,
        ccar3_api = api,
        overrides = refit_args
      ),
      error = function(e) NULL
    )

    if (!is.null(boot_fit)) {
      boot_fit <- align_fit_to_reference(boot_fit, reference_fit)
      U_array[, , b] <- boot_fit$U
      V_array[, , b] <- boot_fit$V
      cor_matrix[b, ] <- boot_fit$cor
    }

    if (progress && (b %% 10 == 0 || b == B)) {
      message(sprintf("Bootstrap replicate %d/%d complete.", b, B))
    }
  }

  valid <- stats::complete.cases(cor_matrix[, 1])
  if (!any(valid)) {
    stop("Bootstrap failed on every replicate.", call. = FALSE)
  }

  list(
    reference_fit = reference_fit,
    U_array = U_array[, , valid, drop = FALSE],
    V_array = V_array[, , valid, drop = FALSE],
    cor_matrix = cor_matrix[valid, , drop = FALSE],
    refit_mode = refit_mode,
    n_success = sum(valid),
    n_requested = B,
    seed = seed
  )
}

build_bootstrap_loading_table <- function(bootstrap_result,
                                          side = c("x", "y"),
                                          component = 1,
                                          top_n = 10,
                                          alpha = 0.05,
                                          threshold = 1e-6) {
  full_table <- build_bootstrap_loading_table_all(
    bootstrap_result = bootstrap_result,
    side = side,
    alpha = alpha,
    threshold = threshold
  )
  table <- full_table[full_table$component == component, , drop = FALSE]
  table <- table[order(abs(table$estimate), decreasing = TRUE), , drop = FALSE]
  utils::head(table, top_n)
}

#' Summarize bootstrap loading uncertainty for every feature and component
#'
#' @param bootstrap_result Output from [bootstrap_cca_uq()].
#' @param side Which data block to summarize, `"x"` or `"y"`.
#' @param alpha Two-sided error level used for empirical bootstrap intervals.
#' @param threshold Absolute value threshold used when computing selection and
#'   sign stability.
#'
#' @return A data frame with one row per feature/component pair containing the
#'   point estimate, bootstrap mean and standard deviation, empirical interval,
#'   selection frequency, and sign stability.
#' @export
build_bootstrap_loading_table_all <- function(bootstrap_result,
                                              side = c("x", "y"),
                                              alpha = 0.05,
                                              threshold = 1e-6) {
  side <- match.arg(side)
  fit <- bootstrap_result$reference_fit
  arr <- if (side == "x") bootstrap_result$U_array else bootstrap_result$V_array
  names <- if (side == "x") fit$x_names else fit$y_names
  summary <- summarize_bootstrap_matrix(arr, alpha = alpha)
  stability <- selection_frequency(arr, threshold = threshold)
  sign_match <- sign_stability(arr, if (side == "x") fit$U else fit$V, threshold = threshold)
  estimates <- if (side == "x") fit$U else fit$V
  direction_prefix <- if (side == "x") "U" else "V"

  tables <- lapply(seq_len(ncol(estimates)), function(component) {
    data.frame(
      side = side,
      feature = names,
      component = component,
      direction = paste0(direction_prefix, component),
      estimate = estimates[, component],
      boot_mean = summary$mean[, component],
      boot_sd = summary$sd[, component],
      ci_lower = summary$lower[, component],
      ci_upper = summary$upper[, component],
      selection_frequency = stability[, component],
      sign_stability = sign_match[, component],
      stringsAsFactors = FALSE
    )
  })

  table <- do.call(rbind, tables)
  rownames(table) <- NULL
  table[order(table$component, -abs(table$estimate)), , drop = FALSE]
}

#' Summarize bootstrap uncertainty for canonical correlations
#'
#' @inheritParams build_bootstrap_loading_table_all
#'
#' @return A data frame with bootstrap means, standard deviations, and
#'   empirical confidence intervals for each canonical correlation.
#' @export
build_bootstrap_correlation_table <- function(bootstrap_result, alpha = 0.05) {
  fit <- bootstrap_result$reference_fit
  summary <- summarize_bootstrap_vector(bootstrap_result$cor_matrix, alpha = alpha)

  data.frame(
    component = seq_along(fit$cor),
    estimate = fit$cor,
    boot_mean = summary$mean,
    boot_sd = summary$sd,
    ci_lower = summary$lower,
    ci_upper = summary$upper,
    stringsAsFactors = FALSE
  )
}

#' Summarize bootstrap uncertainty for a block of canonical components
#'
#' Aligns a possibly non-identifiable block of bootstrap components to the
#' reference fit with a Procrustes rotation and reports subspace distances and
#' principal-angle summaries.
#'
#' @inheritParams build_bootstrap_loading_table_all
#' @param components Integer vector of component indices forming the block to
#'   align and summarize.
#'
#' @return A list containing the reference stacked block, Procrustes-aligned
#'   bootstrap blocks, per-bootstrap subspace metrics, and a one-row summary.
#' @export
summarize_bootstrap_subspace <- function(bootstrap_result,
                                         components = seq_len(ncol(bootstrap_result$reference_fit$U))) {
  components <- sort(unique(as.integer(components)))
  fit <- bootstrap_result$reference_fit
  ref_block <- rbind(
    fit$U[, components, drop = FALSE],
    fit$V[, components, drop = FALSE]
  )
  n_boot <- dim(bootstrap_result$U_array)[3]
  aligned_blocks <- array(
    NA_real_,
    dim = c(nrow(ref_block), length(components), n_boot)
  )
  feature_labels <- c(
    paste0("X::", fit$x_names),
    paste0("Y::", fit$y_names)
  )
  rownames(ref_block) <- feature_labels

  metrics <- lapply(seq_len(n_boot), function(index) {
    u_block <- matrix(
      bootstrap_result$U_array[, components, index],
      nrow = nrow(bootstrap_result$U_array),
      ncol = length(components)
    )
    v_block <- matrix(
      bootstrap_result$V_array[, components, index],
      nrow = nrow(bootstrap_result$V_array),
      ncol = length(components)
    )
    estimate_block <- rbind(
      u_block,
      v_block
    )
    aligned <- procrustes_align_matrix(estimate_block, ref_block)
    principal_angles <- compute_principal_angles(estimate_block, ref_block)
    aligned_blocks[, , index] <- aligned$aligned

    data.frame(
      bootstrap = index,
      component_block = paste(components, collapse = ","),
      subspace_distance = compute_subspace_distance(estimate_block, ref_block),
      mean_principal_angle = mean(principal_angles),
      max_principal_angle = max(principal_angles),
      stringsAsFactors = FALSE
    )
  })

  metrics <- do.call(rbind, metrics)
  summary <- data.frame(
    component_block = paste(components, collapse = ","),
    mean_subspace_distance = mean(metrics$subspace_distance, na.rm = TRUE),
    sd_subspace_distance = stats::sd(metrics$subspace_distance, na.rm = TRUE),
    mean_principal_angle = mean(metrics$mean_principal_angle, na.rm = TRUE),
    mean_max_principal_angle = mean(metrics$max_principal_angle, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(
    reference_block = ref_block,
    aligned_blocks = aligned_blocks,
    metrics = metrics,
    summary = summary,
    feature_labels = feature_labels
  )
}

plot_coefficient_ci_bars <- function(coefficient_table,
                                     output_path,
                                     title_prefix,
                                     estimate_col = "estimate",
                                     lower_col = "ci_lower",
                                     upper_col = "ci_upper",
                                     fill = "#4C78A8",
                                     features_per_page = 40) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to draw loading CI plots.", call. = FALSE)
  }

  coefficient_table <- as.data.frame(coefficient_table)
  required_columns <- c("feature", "component", estimate_col, lower_col, upper_col)
  missing_columns <- setdiff(required_columns, names(coefficient_table))
  if (length(missing_columns)) {
    stop(
      sprintf("Missing required columns for CI plotting: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!nrow(coefficient_table)) {
    stop("`coefficient_table` is empty.", call. = FALSE)
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(output_path, width = 12, height = 10, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  features_per_page <- max(1L, as.integer(features_per_page))
  components <- sort(unique(coefficient_table$component))
  for (component in components) {
    component_table <- coefficient_table[coefficient_table$component == component, , drop = FALSE]
    component_table <- component_table[
      order(abs(component_table[[estimate_col]]), decreasing = TRUE),
      ,
      drop = FALSE
    ]
    n_pages <- ceiling(nrow(component_table) / features_per_page)

    for (page in seq_len(n_pages)) {
      start <- (page - 1) * features_per_page + 1
      stop <- min(page * features_per_page, nrow(component_table))
      page_table <- component_table[start:stop, , drop = FALSE]
      page_table <- page_table[order(page_table[[estimate_col]]), , drop = FALSE]
      page_table$plot_estimate <- page_table[[estimate_col]]
      page_table$plot_lower <- page_table[[lower_col]]
      page_table$plot_upper <- page_table[[upper_col]]
      page_table$feature <- factor(page_table$feature, levels = page_table$feature)
      page_label <- if ("direction" %in% names(page_table)) {
        unique(page_table$direction)
      } else {
        paste("component", component)
      }

      plot_object <- ggplot2::ggplot(
        page_table,
        ggplot2::aes(x = feature, y = plot_estimate)
      ) +
        ggplot2::geom_col(fill = fill, alpha = 0.85) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = plot_lower, ymax = plot_upper),
          width = 0.2,
          linewidth = 0.3,
          color = "black"
        ) +
        ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "#B22222") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = sprintf(
            "%s: %s (features %d-%d of %d)",
            title_prefix,
            page_label,
            start,
            stop,
            nrow(component_table)
          ),
          x = NULL,
          y = "Coefficient"
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(face = "bold")
        )

      print(plot_object)
    }
  }

  invisible(output_path)
}

plot_bootstrap_loading_cis <- function(loading_table,
                                       output_path,
                                       title_prefix,
                                       fill = "#4C78A8",
                                       features_per_page = 40) {
  plot_coefficient_ci_bars(
    coefficient_table = loading_table,
    output_path = output_path,
    title_prefix = title_prefix,
    estimate_col = "estimate",
    lower_col = "ci_lower",
    upper_col = "ci_upper",
    fill = fill,
    features_per_page = features_per_page
  )
}

plot_crossfit_loading_cis <- function(loading_table,
                                      output_path,
                                      title_prefix,
                                      fill = "#54A24B",
                                      features_per_page = 40) {
  plot_coefficient_ci_bars(
    coefficient_table = loading_table,
    output_path = output_path,
    title_prefix = title_prefix,
    estimate_col = "beta",
    lower_col = "ci_lower",
    upper_col = "ci_upper",
    fill = fill,
    features_per_page = features_per_page
  )
}

plot_crossfit_variance_diagnostics <- function(variance_table,
                                               output_path,
                                               title_prefix = "Cross-fit variance shrinkage diagnostic") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to draw variance diagnostics.", call. = FALSE)
  }

  variance_table <- as.data.frame(variance_table)
  if (!nrow(variance_table)) {
    stop("`variance_table` is empty.", call. = FALSE)
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  variance_table$side <- factor(variance_table$side, levels = c("x", "y"), labels = c("X block", "Y block"))
  variance_table$component <- factor(variance_table$component, levels = sort(unique(variance_table$component)))
  variance_table$raw_plot <- pmax(variance_table$raw_s2, 1e-12)
  variance_table$shrunken_plot <- pmax(variance_table$shrunken_s2, 1e-12)

  plot_object <- ggplot2::ggplot(
    variance_table,
    ggplot2::aes(x = raw_plot, y = shrunken_plot)
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2,
      color = "#B22222"
    ) +
    ggplot2::geom_point(
      alpha = 0.2,
      size = 0.7,
      color = "#2E6F95"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_grid(side ~ component) +
    ggplot2::labs(
      title = title_prefix,
      subtitle = "Diagonal = unchanged variances; points above the line indicate larger reported variances.",
      x = "Raw residual variance",
      y = "Shrunken residual variance"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  grDevices::pdf(output_path, width = 12, height = 8, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_object)

  invisible(output_path)
}

#' Write CSV and PDF bootstrap reports
#'
#' @param results A list containing at least a `bootstrap` element produced by
#'   [bootstrap_cca_uq()].
#' @param output_dir Directory where the report artifacts should be written.
#' @param alpha,threshold Passed to [build_bootstrap_loading_table_all()].
#' @param features_per_page Number of coefficients to show on each PDF page.
#' @param verbose Whether to print the output locations.
#'
#' @return A list of generated tables and artifact paths.
#' @export
write_bootstrap_loading_reports <- function(results,
                                            output_dir,
                                            alpha = 0.05,
                                            threshold = 1e-6,
                                            features_per_page = 40,
                                            verbose = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  gene_table <- build_bootstrap_loading_table_all(
    bootstrap_result = results$bootstrap,
    side = "x",
    alpha = alpha,
    threshold = threshold
  )
  methylation_table <- build_bootstrap_loading_table_all(
    bootstrap_result = results$bootstrap,
    side = "y",
    alpha = alpha,
    threshold = threshold
  )

  gene_csv <- file.path(output_dir, "bootstrap_gene_loading_cis.csv")
  methylation_csv <- file.path(output_dir, "bootstrap_methylation_loading_cis.csv")
  utils::write.csv(gene_table, gene_csv, row.names = FALSE)
  utils::write.csv(methylation_table, methylation_csv, row.names = FALSE)

  gene_pdf <- file.path(output_dir, "bootstrap_gene_loading_cis.pdf")
  methylation_pdf <- file.path(output_dir, "bootstrap_methylation_loading_cis.pdf")

  plot_bootstrap_loading_cis(
    loading_table = gene_table,
    output_path = gene_pdf,
    title_prefix = "Alcohol example bootstrap CI plot for gene directions",
    fill = "#4C78A8",
    features_per_page = features_per_page
  )
  plot_bootstrap_loading_cis(
    loading_table = methylation_table,
    output_path = methylation_pdf,
    title_prefix = "Alcohol example bootstrap CI plot for methylation directions",
    fill = "#F58518",
    features_per_page = features_per_page
  )

  if (verbose) {
    message(sprintf("Saved full bootstrap CI tables to %s and %s.", gene_csv, methylation_csv))
    message(sprintf("Saved bootstrap loading CI plots to %s and %s.", gene_pdf, methylation_pdf))
  }

  list(
    gene_table = gene_table,
    methylation_table = methylation_table,
    gene_csv = gene_csv,
    methylation_csv = methylation_csv,
    gene_pdf = gene_pdf,
    methylation_pdf = methylation_pdf
  )
}

#' Write CSV and PDF cross-fitted inference reports
#'
#' @param results A list containing at least a `crossfit` element produced by
#'   one of the cross-fitted inference routines.
#' @param output_dir Directory where the report artifacts should be written.
#' @param features_per_page Number of coefficients to show on each PDF page.
#' @param verbose Whether to print the output locations.
#'
#' @return A list of generated tables and artifact paths.
#' @export
write_crossfit_loading_reports <- function(results,
                                           output_dir,
                                           features_per_page = 40,
                                           verbose = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  gene_table <- results$crossfit$x_results
  methylation_table <- results$crossfit$y_results
  variance_table <- results$crossfit$variance_diagnostics
  alignment_table <- results$crossfit$alignment_diagnostics
  preprocess_table <- results$crossfit$preprocess_diagnostics
  shrinkage_label <- variance_shrinkage_label(results$crossfit$variance_shrinkage %||% "james-stein")

  gene_csv <- file.path(output_dir, "crossfit_gene_regression_cis.csv")
  methylation_csv <- file.path(output_dir, "crossfit_methylation_regression_cis.csv")
  alignment_csv <- file.path(output_dir, "crossfit_alignment_diagnostics.csv")
  preprocess_csv <- file.path(output_dir, "crossfit_preprocess_diagnostics.csv")
  utils::write.csv(gene_table, gene_csv, row.names = FALSE)
  utils::write.csv(methylation_table, methylation_csv, row.names = FALSE)
  utils::write.csv(alignment_table, alignment_csv, row.names = FALSE)
  utils::write.csv(preprocess_table, preprocess_csv, row.names = FALSE)
  variance_csv <- NULL
  if (!is.null(variance_table) && nrow(variance_table)) {
    variance_csv <- file.path(output_dir, "crossfit_variance_shrinkage.csv")
    utils::write.csv(variance_table, variance_csv, row.names = FALSE)
  }

  gene_pdf <- file.path(output_dir, "crossfit_gene_regression_cis.pdf")
  methylation_pdf <- file.path(output_dir, "crossfit_methylation_regression_cis.pdf")
  variance_pdf <- NULL

  plot_crossfit_loading_cis(
    loading_table = gene_table,
    output_path = gene_pdf,
    title_prefix = sprintf(
      "Alcohol example cross-fit CCT-inverted CI plot for gene coefficients (%s)",
      shrinkage_label
    ),
    fill = "#54A24B",
    features_per_page = features_per_page
  )
  plot_crossfit_loading_cis(
    loading_table = methylation_table,
    output_path = methylation_pdf,
    title_prefix = sprintf(
      "Alcohol example cross-fit CCT-inverted CI plot for methylation coefficients (%s)",
      shrinkage_label
    ),
    fill = "#E45756",
    features_per_page = features_per_page
  )
  if (!is.null(variance_table) && nrow(variance_table)) {
    variance_pdf <- file.path(output_dir, "crossfit_variance_shrinkage.pdf")
    plot_crossfit_variance_diagnostics(
      variance_table = variance_table,
      output_path = variance_pdf,
      title_prefix = sprintf(
        "Alcohol example cross-fit variance shrinkage diagnostic (%s)",
        shrinkage_label
      )
    )
  }

  if (verbose) {
    message(sprintf("Saved full cross-fit CI tables to %s and %s.", gene_csv, methylation_csv))
    if (!is.null(variance_csv)) {
      message(sprintf("Saved cross-fit variance diagnostics to %s.", variance_csv))
    }
    message(sprintf("Saved cross-fit alignment diagnostics to %s.", alignment_csv))
    message(sprintf("Saved cross-fit preprocess diagnostics to %s.", preprocess_csv))
    message(sprintf("Saved cross-fit CI plots to %s and %s.", gene_pdf, methylation_pdf))
    if (!is.null(variance_pdf)) {
      message(sprintf("Saved cross-fit variance shrinkage plot to %s.", variance_pdf))
    }
  }

  list(
    gene_table = gene_table,
    methylation_table = methylation_table,
    variance_table = variance_table,
    alignment_table = alignment_table,
    preprocess_table = preprocess_table,
    gene_csv = gene_csv,
    methylation_csv = methylation_csv,
    variance_csv = variance_csv,
    alignment_csv = alignment_csv,
    preprocess_csv = preprocess_csv,
    gene_pdf = gene_pdf,
    methylation_pdf = methylation_pdf,
    variance_pdf = variance_pdf
  )
}

batch_simple_reg <- function(Ymat, x) {
  n <- nrow(Ymat)
  x <- as.numeric(x)
  x_centered <- x - mean(x)
  sxx <- sum(x_centered^2)
  df <- n - 2

  if (!is.finite(sxx) || sxx <= .Machine$double.eps || df <= 0) {
    return(list(
      beta1 = rep(NA_real_, ncol(Ymat)),
      s2 = rep(NA_real_, ncol(Ymat)),
      df = max(df, 1),
      sxx = NA_real_
    ))
  }

  ybar <- colMeans(Ymat)
  slopes <- colSums((Ymat - matrix(ybar, n, ncol(Ymat), byrow = TRUE)) * x_centered) / sxx
  intercepts <- ybar - slopes * mean(x)

  fitted <- outer(x, slopes) + matrix(1, n, 1) %*% t(intercepts)
  residuals <- Ymat - fitted
  s2 <- colSums(residuals^2) / df

  list(beta1 = slopes, s2 = s2, df = df, sxx = sxx)
}

match_variance_shrinkage <- function(method) {
  choices <- c("james-stein", "mr.mashr", "none")

  if (is.null(method) || identical(method, choices)) {
    return(choices[1L])
  }
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("`variance_shrinkage` must be a single non-missing string.", call. = FALSE)
  }

  normalized <- tolower(trimws(method))
  alias_map <- c(
    "james-stein" = "james-stein",
    "james_stein" = "james-stein",
    "jamesstein" = "james-stein",
    "mr.mashr" = "mr.mashr",
    "mr_mashr" = "mr.mashr",
    "mrmashr" = "mr.mashr",
    "mr-mashr" = "mr.mashr",
    "mr.ash" = "mr.mashr",
    "mr.ashr" = "mr.mashr",
    "mr_ash" = "mr.mashr",
    "mrash" = "mr.mashr",
    "mr-ash" = "mr.mashr",
    "none" = "none",
    "no-shrinkage" = "none",
    "no_shrinkage" = "none"
  )

  if (normalized %in% names(alias_map)) {
    return(unname(alias_map[[normalized]]))
  }

  stop(
    sprintf(
      "Unknown variance shrinkage method `%s`. Use one of `james-stein`, `mr.mashr`, or `none`.",
      method
    ),
    call. = FALSE
  )
}

variance_shrinkage_label <- function(method) {
  method <- match_variance_shrinkage(method)
  switch(
    method,
    "james-stein" = "James-Stein variance shrinkage",
    "mr.mashr" = "mr.mashr variance shrinkage",
    "none" = "No variance shrinkage"
  )
}

sanitize_path_label <- function(label) {
  sanitized <- gsub("[^A-Za-z0-9]+", "_", label)
  sanitized <- gsub("^_+|_+$", "", sanitized)
  if (!nzchar(sanitized)) "output" else sanitized
}

resolve_mr_ash_backend <- function(preferred = NULL) {
  candidates <- unique(c(preferred, "mr.mashr", "mr.ash.alpha", "mr.ash", "mrash.alpha"))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  available <- candidates[
    vapply(candidates, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]

  if (!length(available)) {
    return(NULL)
  }

  available[[1]]
}

extract_object_component <- function(object, path) {
  parts <- strsplit(path, "\\$", perl = TRUE)[[1]]
  value <- object

  for (part in parts) {
    if (is.null(value)) {
      return(NULL)
    }

    if (is.list(value) && !is.null(value[[part]])) {
      value <- value[[part]]
      next
    }

    if (is.environment(value) && exists(part, envir = value, inherits = FALSE)) {
      value <- get(part, envir = value, inherits = FALSE)
      next
    }

    if (methods::is(value, "S4") && part %in% methods::slotNames(value)) {
      value <- methods::slot(value, part)
      next
    }

    return(NULL)
  }

  value
}

extract_numeric_component <- function(object, candidates) {
  for (candidate in candidates) {
    value <- extract_object_component(object, candidate)
    if (is.null(value)) {
      next
    }

    if (is.matrix(value) || is.array(value)) {
      value <- as.numeric(value)
    }

    if (is.numeric(value)) {
      return(as.numeric(value))
    }
  }

  NULL
}

fit_mr_ash_variances <- function(Ymat, x, mr_ash_args = list()) {
  backend <- resolve_mr_ash_backend(
    preferred = mr_ash_args$backend %||% mr_ash_args$package %||% NULL
  )
  if (is.null(backend)) {
    stop(
      paste(
        "Variance shrinkage method `mr.mashr` requires the `mr.mashr` package",
        "(older backends `mr.ash.alpha`, `mr.ash`, and `mrash.alpha` are also accepted),",
        "but none is installed."
      ),
      call. = FALSE
    )
  }

  fit_args <- mr_ash_args[setdiff(names(mr_ash_args), c("backend", "package"))]
  if (identical(backend, "mr.mashr")) {
    backend_ns <- asNamespace(backend)
    sumstats <- backend_ns$compute_univariate_sumstats(
      X = matrix(as.numeric(x), ncol = 1L),
      Y = as.matrix(Ymat),
      standardize = FALSE,
      standardize.response = FALSE
    )
    grid <- tryCatch(
      backend_ns$autoselect.mixsd(sumstats, mult = sqrt(2))^2,
      error = function(e) numeric()
    )
    grid <- unique(as.numeric(grid[is.finite(grid) & grid > 0]))
    if (!length(grid)) {
      grid <- 1
    }

    default_s0 <- backend_ns$compute_canonical_covs(
      ncol(Ymat),
      singletons = TRUE,
      hetgrid = c(0, 0.25, 0.5, 0.75, 1)
    )
    default_s0 <- backend_ns$expand_covs(default_s0, grid, zeromat = TRUE)
    fit_fun <- get("mr.mash", envir = backend_ns, inherits = TRUE)
    candidate_args <- modifyList(
      list(
        X = matrix(as.numeric(x), ncol = 1L),
        Y = as.matrix(Ymat),
        S0 = default_s0,
        standardize = FALSE,
        verbose = FALSE,
        update_V = TRUE,
        nthreads = 1L
      ),
      fit_args,
      keep.null = TRUE
    )
  } else {
    fit_fun <- get("mr.ash", envir = asNamespace(backend), inherits = TRUE)
    candidate_args <- modifyList(
      list(
        X = matrix(as.numeric(x), ncol = 1L),
        Y = as.matrix(Ymat),
        x = matrix(as.numeric(x), ncol = 1L),
        y = as.matrix(Ymat),
        intercept = TRUE,
        standardize = FALSE,
        standardizeX = FALSE,
        standardizeY = FALSE,
        center = FALSE,
        verbose = FALSE
      ),
      fit_args,
      keep.null = TRUE
    )
  }
  fit_formals <- names(formals(fit_fun))
  fit_args <- candidate_args[intersect(names(candidate_args), fit_formals)]

  fit <- tryCatch(
    do.call(fit_fun, fit_args),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    stop(
      sprintf("mr.mashr variance shrinkage failed: %s", conditionMessage(fit)),
      call. = FALSE
    )
  }

  sigma2 <- NULL
  if (identical(backend, "mr.mashr")) {
    residual_cov <- extract_object_component(fit, "V")
    if (is.matrix(residual_cov) &&
        nrow(residual_cov) == ncol(residual_cov) &&
        ncol(residual_cov) == ncol(Ymat)) {
      sigma2 <- diag(residual_cov)
    }
  }

  if (is.null(sigma2)) {
    sigma2 <- extract_numeric_component(
      fit,
      c(
        "sigma2",
        "sigma2hat",
        "sigma2_hat",
        "sigma.sq",
        "residual_variance",
        "residual_variances",
        "fit$sigma2",
        "fitted$sigma2",
        "model$sigma2"
      )
    )
  }
  if (is.null(sigma2)) {
    stop("Could not extract residual variances from the mr.mashr fit.", call. = FALSE)
  }

  if (length(sigma2) == 1L) {
    sigma2 <- rep(sigma2, ncol(Ymat))
  }
  if (length(sigma2) != ncol(Ymat)) {
    stop(
      sprintf(
        "mr.mashr returned %d residual variances for %d responses.",
        length(sigma2),
        ncol(Ymat)
      ),
      call. = FALSE
    )
  }

  positive <- sigma2[is.finite(sigma2) & sigma2 > 0]
  if (!length(positive)) {
    stop("mr.mashr did not return any positive residual variance estimates.", call. = FALSE)
  }

  sigma2[!is.finite(sigma2) | sigma2 <= 0] <- mean(positive)
  list(
    s2_shrunk = as.numeric(sigma2),
    lambda = NA_real_,
    mu = mean(sigma2),
    backend = backend,
    method = "mr.mashr"
  )
}

shrink_variances_james_stein <- function(s2) {
  mu <- mean(s2, na.rm = TRUE)
  variance <- stats::var(s2, na.rm = TRUE)
  lambda <- if (is.na(variance) || variance <= 0) {
    1
  } else {
    min(1, max(0, variance / (variance + mu^2)))
  }

  list(
    s2_shrunk = (1 - lambda) * s2 + lambda * mu,
    lambda = lambda,
    mu = mu,
    backend = NA_character_,
    method = "james-stein"
  )
}

shrink_variances_none <- function(s2) {
  list(
    s2_shrunk = s2,
    lambda = 0,
    mu = mean(s2, na.rm = TRUE),
    backend = NA_character_,
    method = "none"
  )
}

shrink_variances <- function(s2,
                             method = c("james-stein", "mr.mashr", "none"),
                             Ymat = NULL,
                             x = NULL,
                             mr_mashr_args = list(),
                             mr_ash_args = NULL) {
  method <- match_variance_shrinkage(method)
  if (is.null(mr_ash_args)) {
    mr_ash_args <- mr_mashr_args
  } else if (!length(mr_mashr_args)) {
    mr_mashr_args <- mr_ash_args
  } else {
    mr_mashr_args <- modifyList(mr_ash_args, mr_mashr_args, keep.null = TRUE)
  }

  if (identical(method, "james-stein")) {
    return(shrink_variances_james_stein(s2))
  }
  if (identical(method, "none")) {
    return(shrink_variances_none(s2))
  }
  if (is.null(Ymat) || is.null(x)) {
    stop("`Ymat` and `x` are required for `mr.mashr` variance shrinkage.", call. = FALSE)
  }

  fit_mr_ash_variances(Ymat = Ymat, x = x, mr_ash_args = mr_mashr_args)
}

build_ci <- function(beta1, s2_shrunk, sxx, df, alpha = 0.05) {
  if (!is.finite(sxx) || sxx <= .Machine$double.eps) {
    return(list(
      se = rep(NA_real_, length(beta1)),
      lower = rep(NA_real_, length(beta1)),
      upper = rep(NA_real_, length(beta1)),
      p_value = rep(NA_real_, length(beta1))
    ))
  }

  se <- sqrt(s2_shrunk / sxx)
  tcrit <- stats::qt(1 - alpha / 2, df = df)
  tval <- beta1 / se

  list(
    se = se,
    lower = beta1 - tcrit * se,
    upper = beta1 + tcrit * se,
    p_value = 2 * stats::pt(-abs(tval), df = df)
  )
}

make_balanced_folds <- function(n, K, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  sample(rep(seq_len(K), length.out = n))
}

aggregate_cauchy_pvalues <- function(pvals) {
  pvals <- as.numeric(pvals)
  pvals <- pvals[is.finite(pvals)]
  if (!length(pvals)) {
    return(NA_real_)
  }

  # Stabilize the tangent transform near 0 and 1.
  pvals <- pmin(pmax(pvals, 1e-15), 1 - 1e-15)
  statistic <- mean(tan((0.5 - pvals) * pi))
  0.5 - atan(statistic) / pi
}

#' Aggregate p-values with the Cauchy combination rule.
#'
#' @param pvalue_list List of numeric p-value matrices, one per fold or
#'   resampling split.
#' @return A numeric matrix of aggregated p-values with the same dimensions as
#'   the input matrices.
aggregate_pvalues_cauchy <- function(pvalue_list) {
  n_rows <- nrow(pvalue_list[[1]])
  n_cols <- ncol(pvalue_list[[1]])
  output <- matrix(NA_real_, nrow = n_rows, ncol = n_cols)

  for (i in seq_len(n_rows)) {
    for (j in seq_len(n_cols)) {
      pvals <- vapply(pvalue_list, function(mat) mat[i, j], numeric(1))
      output[i, j] <- aggregate_cauchy_pvalues(pvals)
    }
  }

  output
}

cauchy_combined_pvalue <- function(b0, beta, se, df) {
  beta <- as.numeric(beta)
  se <- as.numeric(se)
  df <- as.numeric(df)

  if (length(df) == 1L && length(beta) > 1L) {
    df <- rep(df, length(beta))
  }

  keep <- is.finite(beta) & is.finite(se) & se > 0 & is.finite(df) & df > 0
  if (!any(keep)) {
    return(NA_real_)
  }

  tval <- (beta[keep] - b0) / se[keep]
  pvals <- 2 * stats::pt(-abs(tval), df = df[keep])
  aggregate_cauchy_pvalues(pvals)
}

weighted_mean_safe <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) {
    return(mean(x, na.rm = TRUE))
  }

  sum(x[keep] * w[keep]) / sum(w[keep])
}

select_cauchy_ci_center <- function(beta, se, df) {
  weights <- 1 / pmax(se^2, 1e-12)
  center <- weighted_mean_safe(beta, weights)
  p_center <- cauchy_combined_pvalue(center, beta = beta, se = se, df = df)

  if (is.finite(p_center) && p_center > 0) {
    return(list(center = center, p_center = p_center))
  }

  search_interval <- c(
    min(beta - 8 * se, na.rm = TRUE),
    max(beta + 8 * se, na.rm = TRUE)
  )
  if (!all(is.finite(search_interval)) || diff(search_interval) <= 0) {
    scale <- max(stats::median(se, na.rm = TRUE), 1)
    search_interval <- c(center - scale, center + scale)
  }

  objective <- function(b0) {
    pval <- cauchy_combined_pvalue(b0, beta = beta, se = se, df = df)
    if (!is.finite(pval)) {
      return(.Machine$double.xmax)
    }
    -pval
  }

  optimum <- tryCatch(
    stats::optimize(objective, interval = search_interval),
    error = function(e) NULL
  )
  if (is.null(optimum)) {
    return(list(center = center, p_center = p_center))
  }

  list(center = optimum$minimum, p_center = -optimum$objective)
}

invert_cauchy_ci <- function(beta,
                             se,
                             df,
                             alpha = 0.05,
                             tol = 1e-6,
                             max_expansions = 60L) {
  beta <- as.numeric(beta)
  se <- as.numeric(se)
  df <- as.numeric(df)

  if (length(df) == 1L && length(beta) > 1L) {
    df <- rep(df, length(beta))
  }

  keep <- is.finite(beta) & is.finite(se) & se > 0 & is.finite(df) & df > 0
  beta <- beta[keep]
  se <- se[keep]
  df <- df[keep]

  if (!length(beta)) {
    return(list(lower = NA_real_, upper = NA_real_, center = NA_real_, p_center = NA_real_))
  }

  center_info <- select_cauchy_ci_center(beta = beta, se = se, df = df)
  center <- center_info$center
  p_center <- center_info$p_center

  if (!is.finite(center)) {
    center <- mean(beta, na.rm = TRUE)
    p_center <- cauchy_combined_pvalue(center, beta = beta, se = se, df = df)
  }

  if (!is.finite(p_center) || p_center <= alpha) {
    return(list(lower = NA_real_, upper = NA_real_, center = center, p_center = p_center))
  }

  target <- function(b0) {
    cauchy_combined_pvalue(b0, beta = beta, se = se, df = df) - alpha
  }

  initial_span <- max(abs(beta - center), na.rm = TRUE)
  initial_span <- max(initial_span, stats::median(se, na.rm = TRUE), 1e-3)

  locate_endpoint <- function(direction = c(-1, 1)) {
    direction <- match.arg(as.character(direction), choices = c("-1", "1"))
    direction <- as.numeric(direction)
    step <- initial_span
    outer <- center + direction * step
    outer_value <- target(outer)
    iter <- 0L

    while (is.finite(outer_value) && outer_value > 0 && iter < max_expansions) {
      step <- step * 2
      outer <- center + direction * step
      outer_value <- target(outer)
      iter <- iter + 1L
    }

    if (!is.finite(outer_value)) {
      return(NA_real_)
    }
    if (outer_value > 0) {
      return(if (direction < 0) -Inf else Inf)
    }

    root <- tryCatch(
      stats::uniroot(target, interval = sort(c(center, outer)), tol = tol)$root,
      error = function(e) NA_real_
    )
    root
  }

  list(
    lower = locate_endpoint(-1),
    upper = locate_endpoint(1),
    center = center,
    p_center = p_center
  )
}

invert_cauchy_ci_matrix <- function(beta_list, se_list, df_list, alpha = 0.05) {
  beta_array <- simplify2array(beta_list)
  se_array <- simplify2array(se_list)
  n_features <- dim(beta_array)[1]
  n_components <- dim(beta_array)[2]

  lower <- matrix(NA_real_, nrow = n_features, ncol = n_components)
  upper <- matrix(NA_real_, nrow = n_features, ncol = n_components)
  center <- matrix(NA_real_, nrow = n_features, ncol = n_components)

  for (component in seq_len(n_components)) {
    df_component <- vapply(df_list, function(x) x[component], numeric(1))
    for (feature in seq_len(n_features)) {
      ci <- invert_cauchy_ci(
        beta = beta_array[feature, component, ],
        se = se_array[feature, component, ],
        df = df_component,
        alpha = alpha
      )
      lower[feature, component] <- ci$lower
      upper[feature, component] <- ci$upper
      center[feature, component] <- ci$center
    }
  }

  list(lower = lower, upper = upper, center = center)
}

average_matrices <- function(matrices) {
  matrix_array <- simplify2array(matrices)
  apply(matrix_array, c(1, 2), mean, na.rm = TRUE)
}

safe_correlation <- function(x, y) {
  if (all(!is.finite(x)) || all(!is.finite(y))) {
    return(NA_real_)
  }
  sx <- stats::sd(x, na.rm = TRUE)
  sy <- stats::sd(y, na.rm = TRUE)
  if (!is.finite(sx) || !is.finite(sy) ||
      sx <= .Machine$double.eps ||
      sy <= .Machine$double.eps) {
    return(NA_real_)
  }
  stats::cor(x, y)
}

build_inference_table <- function(feature_names,
                                  beta,
                                  lower,
                                  upper,
                                  pval,
                                  padj,
                                  alpha,
                                  side,
                                  variance_shrinkage = NULL,
                                  estimate_name = "beta",
                                  inference_target = NULL,
                                  selected = NULL,
                                  backend = NULL) {
  tables <- lapply(seq_len(ncol(beta)), function(component) {
    table <- data.frame(
      side = side,
      feature = feature_names,
      component = component,
      ci_lower = lower[, component],
      ci_upper = upper[, component],
      p_value = pval[, component],
      p_adjusted = padj[, component],
      significant = padj[, component] < alpha,
      stringsAsFactors = FALSE
    )
    table[[estimate_name]] <- beta[, component]
    if (!identical(estimate_name, "beta")) {
      table$beta <- beta[, component]
    }

    if (!is.null(variance_shrinkage)) {
      table$variance_shrinkage <- variance_shrinkage
    }
    if (!is.null(inference_target)) {
      table$inference_target <- inference_target
    }
    if (!is.null(selected)) {
      table$selected <- selected[, component]
    }
    if (!is.null(backend)) {
      table$backend <- backend[, component]
    }

    table
  })

  do.call(rbind, tables)
}

adjust_pvalue_matrix <- function(p_matrix) {
  adjusted <- p_matrix
  for (component in seq_len(ncol(p_matrix))) {
    valid <- is.finite(p_matrix[, component])
    adjusted[, component] <- NA_real_
    if (any(valid)) {
      adjusted[valid, component] <- stats::p.adjust(p_matrix[valid, component], method = "BH")
    }
  }
  adjusted
}

validate_optional_covariates <- function(covariates, n, name = "covariates") {
  if (is.null(covariates)) {
    return(NULL)
  }

  covariates <- as_numeric_matrix(covariates, name)
  if (nrow(covariates) != n) {
    stop(sprintf("`%s` must have %d rows.", name, n), call. = FALSE)
  }

  covariates
}

build_design_matrix <- function(score, covariates = NULL) {
  score <- as.numeric(score)
  if (!length(score)) {
    stop("`score` must be non-empty.", call. = FALSE)
  }

  design <- cbind(`(Intercept)` = 1, score = score)
  if (!is.null(covariates)) {
    covariates <- as.matrix(covariates)
    if (is.null(colnames(covariates))) {
      colnames(covariates) <- paste0("covariate_", seq_len(ncol(covariates)))
    }
    design <- cbind(design, covariates)
  }

  design
}

fit_feature_score_regression <- function(feature_matrix,
                                         score,
                                         covariates = NULL,
                                         alpha = 0.05,
                                         use_limma = TRUE) {
  feature_matrix <- as_numeric_matrix(feature_matrix, "feature_matrix")
  design <- build_design_matrix(score, covariates = covariates)
  n <- nrow(design)
  coef_index <- 2L
  backend <- "ols"

  if (use_limma && requireNamespace("limma", quietly = TRUE)) {
    fit <- limma::lmFit(t(feature_matrix), design)
    fit <- limma::eBayes(fit)
    coef <- fit$coefficients[, coef_index]
    se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, coef_index]
    df <- fit$df.total
    crit <- stats::qt(1 - alpha / 2, df = pmax(df, 1))
    lower <- coef - crit * se
    upper <- coef + crit * se
    p_value <- fit$p.value[, coef_index]
    backend <- "limma"

    return(list(
      coef = coef,
      se = se,
      lower = lower,
      upper = upper,
      p_value = p_value,
      df = df,
      backend = rep(backend, ncol(feature_matrix))
    ))
  }

  if (qr(design)$rank < ncol(design)) {
    stop("Feature-score regression design matrix is rank-deficient.", call. = FALSE)
  }

  coef_matrix <- qr.coef(qr(design), feature_matrix)
  fitted_values <- qr.fitted(qr(design), feature_matrix)
  residuals <- feature_matrix - fitted_values
  df <- n - ncol(design)
  if (df <= 0) {
    stop("Not enough degrees of freedom for feature-score regression.", call. = FALSE)
  }

  s2 <- colSums(residuals^2) / df
  xtx_inv <- solve(crossprod(design))
  coef <- coef_matrix[coef_index, ]
  se <- sqrt(pmax(s2, 0) * xtx_inv[coef_index, coef_index])
  statistic <- coef / se
  p_value <- 2 * stats::pt(-abs(statistic), df = df)
  crit <- stats::qt(1 - alpha / 2, df = df)
  lower <- coef - crit * se
  upper <- coef + crit * se

  list(
    coef = coef,
    se = se,
    lower = lower,
    upper = upper,
    p_value = p_value,
    df = rep(df, ncol(feature_matrix)),
    backend = rep(backend, ncol(feature_matrix))
  )
}

fit_conditional_score_regression <- function(response,
                                             predictors,
                                             alpha = 0.05,
                                             method = c("auto", "ols", "split_lasso"),
                                             split_seed = 1,
                                             selection_fraction = 0.5,
                                             lasso_select = c("lambda.1se", "lambda.min")) {
  method <- match.arg(method)
  lasso_select <- match.arg(lasso_select)
  predictors <- as_numeric_matrix(predictors, "predictors")
  response <- as.numeric(response)
  n <- nrow(predictors)
  p <- ncol(predictors)

  if (length(response) != n) {
    stop("`response` must have the same number of rows as `predictors`.", call. = FALSE)
  }

  predictor_centered <- scale(predictors, center = TRUE, scale = FALSE)
  response_centered <- response - mean(response)
  ols_feasible <- n > (p + 1L) && qr(predictor_centered)$rank == p

  if (identical(method, "auto")) {
    method <- if (ols_feasible) "ols" else "split_lasso"
  }

  if (identical(method, "ols")) {
    if (!ols_feasible) {
      stop("OLS conditional regression is not feasible in the current dimension.", call. = FALSE)
    }

    design <- cbind(`(Intercept)` = 1, predictor_centered)
    coef_matrix <- qr.coef(qr(design), response_centered)
    fitted_values <- qr.fitted(qr(design), response_centered)
    residuals <- response_centered - fitted_values
    df <- n - ncol(design)
    s2 <- sum(residuals^2) / df
    xtx_inv <- solve(crossprod(design))
    beta <- coef_matrix[-1]
    se <- sqrt(pmax(s2, 0) * diag(xtx_inv)[-1])
    statistic <- beta / se
    p_value <- 2 * stats::pt(-abs(statistic), df = df)
    crit <- stats::qt(1 - alpha / 2, df = df)

    return(list(
      beta = beta,
      se = se,
      lower = beta - crit * se,
      upper = beta + crit * se,
      p_value = p_value,
      selected = rep(TRUE, p),
      backend = rep("ols", p),
      diagnostics = data.frame(
        backend = "ols",
        n_selected = p,
        selection_n = n,
        inference_n = n,
        stringsAsFactors = FALSE
      )
    ))
  }

  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for split-lasso conditional inference.", call. = FALSE)
  }

  selection_n <- max(2L, floor(selection_fraction * n))
  selection_n <- min(selection_n, n - 2L)
  if (selection_n < 2L || (n - selection_n) < 3L) {
    stop("Not enough samples for split-lasso conditional inference.", call. = FALSE)
  }

  set.seed(split_seed)
  selection_index <- sample.int(n, size = selection_n, replace = FALSE)
  inference_index <- setdiff(seq_len(n), selection_index)

  x_selection <- predictor_centered[selection_index, , drop = FALSE]
  y_selection <- response_centered[selection_index]
  x_inference <- predictor_centered[inference_index, , drop = FALSE]
  y_inference <- response_centered[inference_index]
  n_inference <- length(inference_index)

  cv_fit <- glmnet::cv.glmnet(
    x = x_selection,
    y = y_selection,
    alpha = 1,
    intercept = TRUE,
    standardize = FALSE
  )
  lambda_name <- if (identical(lasso_select, "lambda.min")) "lambda.min" else "lambda.1se"
  coef_selection <- as.matrix(stats::coef(cv_fit, s = lambda_name))
  selected <- which(abs(coef_selection[-1, 1]) > 0)

  beta <- rep(0, p)
  se <- rep(NA_real_, p)
  lower <- rep(NA_real_, p)
  upper <- rep(NA_real_, p)
  p_value <- rep(1, p)
  selected_mask <- rep(FALSE, p)
  backend <- rep("split_lasso", p)

  if (!length(selected)) {
    return(list(
      beta = beta,
      se = se,
      lower = lower,
      upper = upper,
      p_value = p_value,
      selected = selected_mask,
      backend = backend,
      diagnostics = data.frame(
        backend = "split_lasso",
        n_selected = 0L,
        selection_n = selection_n,
        inference_n = n_inference,
        stringsAsFactors = FALSE
      )
    ))
  }

  max_selected <- max(1L, n_inference - 2L)
  if (length(selected) > max_selected) {
    ranked <- order(abs(coef_selection[selected + 1L, 1]), decreasing = TRUE)
    selected <- selected[ranked[seq_len(max_selected)]]
  }

  selected_mask[selected] <- TRUE
  design_inference <- cbind(`(Intercept)` = 1, x_inference[, selected, drop = FALSE])
  if (qr(design_inference)$rank < ncol(design_inference)) {
    selected <- selected[seq_len(max(1L, min(length(selected), n_inference - 2L)))]
    selected_mask <- rep(FALSE, p)
    selected_mask[selected] <- TRUE
    design_inference <- cbind(`(Intercept)` = 1, x_inference[, selected, drop = FALSE])
  }

  if (qr(design_inference)$rank < ncol(design_inference)) {
    return(list(
      beta = beta,
      se = se,
      lower = lower,
      upper = upper,
      p_value = p_value,
      selected = selected_mask,
      backend = rep("split_lasso_rank_deficient", p),
      diagnostics = data.frame(
        backend = "split_lasso_rank_deficient",
        n_selected = sum(selected_mask),
        selection_n = selection_n,
        inference_n = n_inference,
        stringsAsFactors = FALSE
      )
    ))
  }

  coef_inference <- qr.coef(qr(design_inference), y_inference)
  fitted_values <- qr.fitted(qr(design_inference), y_inference)
  residuals <- y_inference - fitted_values
  df <- n_inference - ncol(design_inference)
  s2 <- sum(residuals^2) / df
  xtx_inv <- solve(crossprod(design_inference))
  beta_selected <- coef_inference[-1]
  se_selected <- sqrt(pmax(s2, 0) * diag(xtx_inv)[-1])
  statistic <- beta_selected / se_selected
  p_selected <- 2 * stats::pt(-abs(statistic), df = df)
  crit <- stats::qt(1 - alpha / 2, df = df)

  beta[selected] <- beta_selected
  se[selected] <- se_selected
  lower[selected] <- beta_selected - crit * se_selected
  upper[selected] <- beta_selected + crit * se_selected
  p_value[selected] <- p_selected

  list(
    beta = beta,
    se = se,
    lower = lower,
    upper = upper,
    p_value = p_value,
    selected = selected_mask,
    backend = backend,
    diagnostics = data.frame(
      backend = "split_lasso",
      n_selected = sum(selected_mask),
      selection_n = selection_n,
      inference_n = n_inference,
      stringsAsFactors = FALSE
    )
  )
}

compute_crossfit_scores <- function(X, Y,
                                    reference_fit,
                                    K = 5,
                                    seed = 1,
                                    fit_mode = c("fixed_lambda", "retune_cv"),
                                    fit_args = list(),
                                    ccar3_api = NULL,
                                    progress = interactive()) {
  fit_mode <- match.arg(fit_mode)
  X <- ensure_feature_names(as_numeric_matrix(X, "X"), "X_")
  Y <- ensure_feature_names(as_numeric_matrix(Y, "Y"), "Y_")
  n <- nrow(X)
  r <- ncol(reference_fit$U)
  api <- ccar3_api %||% get_ccar3_api(
    ccar3_path = reference_fit$settings$ccar3_path,
    prefer_source = reference_fit$settings$prefer_source
  )
  preprocessed <- preprocess_once_for_crossfit(X, Y, reference_fit)
  X_cf <- preprocessed$X
  Y_cf <- preprocessed$Y
  refit_args <- modifyList(fit_args, list(preprocess_mode = "none"), keep.null = TRUE)
  folds <- make_balanced_folds(n, K, seed = seed)

  x_scores <- matrix(NA_real_, nrow = n, ncol = r)
  y_scores <- matrix(NA_real_, nrow = n, ncol = r)
  alignment_diagnostics <- vector("list", K)
  fold_correlations <- matrix(NA_real_, nrow = K, ncol = r)
  successful_folds <- logical(K)

  for (fold in seq_len(K)) {
    train <- which(folds != fold)
    test <- which(folds == fold)

    fit_fold <- tryCatch(
      refit_from_reference(
        X = X_cf[train, , drop = FALSE],
        Y = Y_cf[train, , drop = FALSE],
        reference_fit = reference_fit,
        mode = fit_mode,
        ccar3_api = api,
        overrides = refit_args
      ),
      error = function(e) NULL
    )

    if (is.null(fit_fold)) {
      next
    }

    fit_fold <- align_fit_to_reference(fit_fold, reference_fit)
    scores <- canonical_scores(X_cf[test, , drop = FALSE], Y_cf[test, , drop = FALSE], fit_fold)
    x_scores[test, ] <- scores$X
    y_scores[test, ] <- scores$Y
    successful_folds[fold] <- TRUE

    alignment_diagnostics[[fold]] <- data.frame(
      fold = fold,
      component = seq_len(r),
      matched_fold_component = fit_fold$alignment$permutation,
      sign = fit_fold$alignment$sign,
      signed_similarity = diag(fit_fold$alignment$aligned_similarity),
      abs_similarity = diag(fit_fold$alignment$aligned_abs_similarity),
      stringsAsFactors = FALSE
    )

    fold_correlations[fold, ] <- vapply(seq_len(r), function(component) {
      safe_correlation(scores$X[, component], scores$Y[, component])
    }, numeric(1))

    if (progress) {
      message(sprintf("Cross-fit score fold %d/%d complete.", fold, K))
    }
  }

  if (!any(successful_folds)) {
    stop("Cross-fitted score construction failed on every fold.", call. = FALSE)
  }

  list(
    x_scores = x_scores,
    y_scores = y_scores,
    X_preprocessed = X_cf,
    Y_preprocessed = Y_cf,
    fold_correlations = fold_correlations[successful_folds, , drop = FALSE],
    alignment_diagnostics = do.call(rbind, alignment_diagnostics[successful_folds]),
    preprocess_diagnostics = preprocessed$diagnostics,
    folds = folds,
    fit_mode = fit_mode,
    K = sum(successful_folds)
  )
}

#' Cross-fitted marginal feature-score regression for sparse CCA
#'
#' Uses held-out canonical scores to estimate marginal feature associations with
#' the opposite-view canonical variates, optionally with a common covariate
#' adjustment and limma-style variance moderation.
#'
#' @inheritParams bootstrap_cca_uq
#' @param K Number of folds used for cross-fitting.
#' @param alpha Two-sided error level used for confidence intervals and
#'   significance flags.
#' @param fit_mode Whether each fold refit reuses the reference penalty
#'   (`"fixed_lambda"`) or re-runs cross-validation (`"retune_cv"`).
#' @param fit_args Optional named overrides passed to each fold refit.
#' @param covariates Optional matrix or data frame of additional covariates to
#'   include in every univariate regression.
#' @param use_limma Whether to use `limma` empirical-Bayes moderation when it
#'   is available.
#'
#' @return A list containing marginal regression summaries for the X and Y
#'   blocks, cross-fitted canonical scores, fold diagnostics, and metadata about
#'   the inference backend.
#' @export
crossfit_cca_marginal_inference <- function(X, Y,
                                            reference_fit,
                                            K = 5,
                                            seed = 1,
                                            alpha = 0.05,
                                            fit_mode = c("fixed_lambda", "retune_cv"),
                                            fit_args = list(),
                                            ccar3_api = NULL,
                                            progress = interactive(),
                                            covariates = NULL,
                                            use_limma = TRUE) {
  fit_mode <- match.arg(fit_mode)
  X <- ensure_feature_names(as_numeric_matrix(X, "X"), "X_")
  Y <- ensure_feature_names(as_numeric_matrix(Y, "Y"), "Y_")
  covariates <- validate_optional_covariates(covariates, nrow(X))

  scores <- compute_crossfit_scores(
    X = X,
    Y = Y,
    reference_fit = reference_fit,
    K = K,
    seed = seed,
    fit_mode = fit_mode,
    fit_args = fit_args,
    ccar3_api = ccar3_api,
    progress = progress
  )

  p <- ncol(X)
  q <- ncol(Y)
  r <- ncol(reference_fit$U)
  alpha_x <- se_x <- lower_x <- upper_x <- pval_x <- matrix(NA_real_, nrow = p, ncol = r)
  alpha_y <- se_y <- lower_y <- upper_y <- pval_y <- matrix(NA_real_, nrow = q, ncol = r)
  backend_x <- matrix(NA_character_, nrow = p, ncol = r)
  backend_y <- matrix(NA_character_, nrow = q, ncol = r)

  for (component in seq_len(r)) {
    valid_y <- is.finite(scores$y_scores[, component])
    if (sum(valid_y) > 2L) {
      fit_x <- fit_feature_score_regression(
        feature_matrix = scores$X_preprocessed[valid_y, , drop = FALSE],
        score = scores$y_scores[valid_y, component],
        covariates = if (is.null(covariates)) NULL else covariates[valid_y, , drop = FALSE],
        alpha = alpha,
        use_limma = use_limma
      )
      alpha_x[, component] <- fit_x$coef
      se_x[, component] <- fit_x$se
      lower_x[, component] <- fit_x$lower
      upper_x[, component] <- fit_x$upper
      pval_x[, component] <- fit_x$p_value
      backend_x[, component] <- fit_x$backend
    }

    valid_x <- is.finite(scores$x_scores[, component])
    if (sum(valid_x) > 2L) {
      fit_y <- fit_feature_score_regression(
        feature_matrix = scores$Y_preprocessed[valid_x, , drop = FALSE],
        score = scores$x_scores[valid_x, component],
        covariates = if (is.null(covariates)) NULL else covariates[valid_x, , drop = FALSE],
        alpha = alpha,
        use_limma = use_limma
      )
      alpha_y[, component] <- fit_y$coef
      se_y[, component] <- fit_y$se
      lower_y[, component] <- fit_y$lower
      upper_y[, component] <- fit_y$upper
      pval_y[, component] <- fit_y$p_value
      backend_y[, component] <- fit_y$backend
    }
  }

  adjusted_p_x <- adjust_pvalue_matrix(pval_x)
  adjusted_p_y <- adjust_pvalue_matrix(pval_y)

  list(
    x_results = build_inference_table(
      feature_names = colnames(X),
      beta = alpha_x,
      lower = lower_x,
      upper = upper_x,
      pval = pval_x,
      padj = adjusted_p_x,
      alpha = alpha,
      side = "x",
      estimate_name = "alpha",
      inference_target = "marginal",
      backend = backend_x
    ),
    y_results = build_inference_table(
      feature_names = colnames(Y),
      beta = alpha_y,
      lower = lower_y,
      upper = upper_y,
      pval = pval_y,
      padj = adjusted_p_y,
      alpha = alpha,
      side = "y",
      estimate_name = "alpha",
      inference_target = "marginal",
      backend = backend_y
    ),
    x_scores = scores$x_scores,
    y_scores = scores$y_scores,
    X_preprocessed = scores$X_preprocessed,
    Y_preprocessed = scores$Y_preprocessed,
    fold_correlations = scores$fold_correlations,
    alignment_diagnostics = scores$alignment_diagnostics,
    preprocess_diagnostics = scores$preprocess_diagnostics,
    folds = scores$folds,
    fit_mode = fit_mode,
    K = scores$K,
    alpha = alpha,
    variance_diagnostics = data.frame(),
    variance_shrinkage = NULL,
    covariates_included = !is.null(covariates),
    backend = if (use_limma && requireNamespace("limma", quietly = TRUE)) "limma" else "ols"
  )
}

#' Cross-fitted conditional regression for sparse CCA
#'
#' Estimates which variables in one data block conditionally predict the
#' opposite-view canonical score using either ordinary least squares or a
#' sample-split lasso plus post-selection OLS fallback in high dimensions.
#'
#' @inheritParams crossfit_cca_marginal_inference
#' @param method Conditional regression backend: automatic selection, direct
#'   OLS, or sample-split lasso followed by OLS on the selected features.
#' @param selection_fraction Fraction of samples allocated to the lasso
#'   selection split when `method = "split_lasso"`.
#' @param lasso_select Cross-validation rule used to choose the lasso penalty.
#'
#' @return A list containing conditional regression summaries for the X and Y
#'   blocks, cross-fitted canonical scores, fold diagnostics, and backend
#'   details.
#' @export
crossfit_cca_conditional_inference <- function(X, Y,
                                               reference_fit,
                                               K = 5,
                                               seed = 1,
                                               alpha = 0.05,
                                               fit_mode = c("fixed_lambda", "retune_cv"),
                                               fit_args = list(),
                                               ccar3_api = NULL,
                                               progress = interactive(),
                                               method = c("auto", "ols", "split_lasso"),
                                               selection_fraction = 0.5,
                                               lasso_select = c("lambda.1se", "lambda.min")) {
  fit_mode <- match.arg(fit_mode)
  method <- match.arg(method)
  lasso_select <- match.arg(lasso_select)
  X <- ensure_feature_names(as_numeric_matrix(X, "X"), "X_")
  Y <- ensure_feature_names(as_numeric_matrix(Y, "Y"), "Y_")

  scores <- compute_crossfit_scores(
    X = X,
    Y = Y,
    reference_fit = reference_fit,
    K = K,
    seed = seed,
    fit_mode = fit_mode,
    fit_args = fit_args,
    ccar3_api = ccar3_api,
    progress = progress
  )

  p <- ncol(X)
  q <- ncol(Y)
  r <- ncol(reference_fit$U)
  beta_x <- se_x <- lower_x <- upper_x <- pval_x <- matrix(NA_real_, nrow = p, ncol = r)
  beta_y <- se_y <- lower_y <- upper_y <- pval_y <- matrix(NA_real_, nrow = q, ncol = r)
  selected_x <- matrix(FALSE, nrow = p, ncol = r)
  selected_y <- matrix(FALSE, nrow = q, ncol = r)
  backend_x <- matrix(NA_character_, nrow = p, ncol = r)
  backend_y <- matrix(NA_character_, nrow = q, ncol = r)
  diagnostics <- vector("list", 2L * r)

  for (component in seq_len(r)) {
    valid_y <- is.finite(scores$y_scores[, component])
    if (sum(valid_y) > 3L) {
      fit_x <- fit_conditional_score_regression(
        response = scores$y_scores[valid_y, component],
        predictors = scores$X_preprocessed[valid_y, , drop = FALSE],
        alpha = alpha,
        method = method,
        split_seed = seed + component,
        selection_fraction = selection_fraction,
        lasso_select = lasso_select
      )
      beta_x[, component] <- fit_x$beta
      se_x[, component] <- fit_x$se
      lower_x[, component] <- fit_x$lower
      upper_x[, component] <- fit_x$upper
      pval_x[, component] <- fit_x$p_value
      selected_x[, component] <- fit_x$selected
      backend_x[, component] <- fit_x$backend
      diagnostics[[2L * component - 1L]] <- cbind(
        side = "x",
        component = component,
        fit_x$diagnostics
      )
    }

    valid_x <- is.finite(scores$x_scores[, component])
    if (sum(valid_x) > 3L) {
      fit_y <- fit_conditional_score_regression(
        response = scores$x_scores[valid_x, component],
        predictors = scores$Y_preprocessed[valid_x, , drop = FALSE],
        alpha = alpha,
        method = method,
        split_seed = seed + r + component,
        selection_fraction = selection_fraction,
        lasso_select = lasso_select
      )
      beta_y[, component] <- fit_y$beta
      se_y[, component] <- fit_y$se
      lower_y[, component] <- fit_y$lower
      upper_y[, component] <- fit_y$upper
      pval_y[, component] <- fit_y$p_value
      selected_y[, component] <- fit_y$selected
      backend_y[, component] <- fit_y$backend
      diagnostics[[2L * component]] <- cbind(
        side = "y",
        component = component,
        fit_y$diagnostics
      )
    }
  }

  adjusted_p_x <- adjust_pvalue_matrix(pval_x)
  adjusted_p_y <- adjust_pvalue_matrix(pval_y)
  diagnostics <- diagnostics[!vapply(diagnostics, is.null, logical(1))]

  list(
    x_results = build_inference_table(
      feature_names = colnames(X),
      beta = beta_x,
      lower = lower_x,
      upper = upper_x,
      pval = pval_x,
      padj = adjusted_p_x,
      alpha = alpha,
      side = "x",
      estimate_name = "beta",
      inference_target = "conditional",
      selected = selected_x,
      backend = backend_x
    ),
    y_results = build_inference_table(
      feature_names = colnames(Y),
      beta = beta_y,
      lower = lower_y,
      upper = upper_y,
      pval = pval_y,
      padj = adjusted_p_y,
      alpha = alpha,
      side = "y",
      estimate_name = "beta",
      inference_target = "conditional",
      selected = selected_y,
      backend = backend_y
    ),
    x_scores = scores$x_scores,
    y_scores = scores$y_scores,
    X_preprocessed = scores$X_preprocessed,
    Y_preprocessed = scores$Y_preprocessed,
    fold_correlations = scores$fold_correlations,
    alignment_diagnostics = scores$alignment_diagnostics,
    preprocess_diagnostics = scores$preprocess_diagnostics,
    folds = scores$folds,
    fit_mode = fit_mode,
    K = scores$K,
    alpha = alpha,
    variance_diagnostics = data.frame(),
    variance_shrinkage = NULL,
    conditional_diagnostics = if (length(diagnostics)) do.call(rbind, diagnostics) else data.frame(),
    backend = method
  )
}

#' Cross-fitted score-loading regression with Cauchy-combined p-values
#'
#' This legacy cross-fitted routine regresses held-out canonical scores onto the
#' opposite view feature matrix one fold at a time, shrinks foldwise residual
#' variances, and combines fold-specific p-values with the Cauchy rule.
#'
#' @inheritParams crossfit_cca_marginal_inference
#' @param variance_shrinkage Residual-variance shrinkage rule applied within
#'   each fold.
#' @param mr_mashr_args,mr_ash_args Optional tuning lists forwarded to the
#'   `mr.ash` / `mr.mashr` backend when it is requested.
#'
#' @return A list containing X-side and Y-side inference tables, fold
#'   correlations, preprocessing diagnostics, and variance-shrinkage diagnostics.
#' @export
crossfit_cca_inference <- function(X, Y,
                                   reference_fit,
                                   K = 5,
                                   seed = 1,
                                   alpha = 0.05,
                                   fit_mode = c("fixed_lambda", "retune_cv"),
                                   fit_args = list(),
                                   ccar3_api = NULL,
                                   progress = interactive(),
                                   variance_shrinkage = c("james-stein", "mr.mashr", "none"),
                                   mr_mashr_args = list(),
                                   mr_ash_args = NULL) {
  fit_mode <- match.arg(fit_mode)
  variance_shrinkage <- match_variance_shrinkage(variance_shrinkage)
  if (is.null(mr_ash_args)) {
    mr_ash_args <- mr_mashr_args
  } else if (!length(mr_mashr_args)) {
    mr_mashr_args <- mr_ash_args
  } else {
    mr_mashr_args <- modifyList(mr_ash_args, mr_mashr_args, keep.null = TRUE)
  }
  X <- ensure_feature_names(as_numeric_matrix(X, "X"), "X_")
  Y <- ensure_feature_names(as_numeric_matrix(Y, "Y"), "Y_")
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  r <- ncol(reference_fit$U)
  api <- ccar3_api %||% get_ccar3_api(
    ccar3_path = reference_fit$settings$ccar3_path,
    prefer_source = reference_fit$settings$prefer_source
  )
  preprocessed <- preprocess_once_for_crossfit(X, Y, reference_fit)
  X_cf <- preprocessed$X
  Y_cf <- preprocessed$Y
  refit_args <- modifyList(fit_args, list(preprocess_mode = "none"), keep.null = TRUE)
  folds <- make_balanced_folds(n, K, seed = seed)

  pvalues_x <- betas_x <- ses_x <- vector("list", K)
  pvalues_y <- betas_y <- ses_y <- vector("list", K)
  df_x <- df_y <- vector("list", K)
  alignment_diagnostics <- variance_diagnostics <- vector("list", K)
  fold_correlations <- matrix(NA_real_, nrow = K, ncol = r)

  for (fold in seq_len(K)) {
    train <- which(folds != fold)
    test <- which(folds == fold)

    fit_fold <- tryCatch(
      refit_from_reference(
        X = X_cf[train, , drop = FALSE],
        Y = Y_cf[train, , drop = FALSE],
        reference_fit = reference_fit,
        mode = fit_mode,
        ccar3_api = api,
        overrides = refit_args
      ),
      error = function(e) NULL
    )

    if (is.null(fit_fold)) {
      next
    }

    fit_fold <- align_fit_to_reference(fit_fold, reference_fit)
    scores <- canonical_scores(X_cf[test, , drop = FALSE], Y_cf[test, , drop = FALSE], fit_fold)

    alignment_diagnostics[[fold]] <- data.frame(
      fold = fold,
      component = seq_len(r),
      matched_fold_component = fit_fold$alignment$permutation,
      sign = fit_fold$alignment$sign,
      signed_similarity = diag(fit_fold$alignment$aligned_similarity),
      abs_similarity = diag(fit_fold$alignment$aligned_abs_similarity),
      stringsAsFactors = FALSE
    )

    fold_correlations[fold, ] <- vapply(seq_len(r), function(component) {
      safe_correlation(scores$X[, component], scores$Y[, component])
    }, numeric(1))

    fold_p_x <- fold_b_x <- matrix(NA_real_, nrow = p, ncol = r)
    fold_p_y <- fold_b_y <- matrix(NA_real_, nrow = q, ncol = r)
    fold_se_x <- matrix(NA_real_, nrow = p, ncol = r)
    fold_se_y <- matrix(NA_real_, nrow = q, ncol = r)
    fold_df_x <- rep(NA_real_, r)
    fold_df_y <- rep(NA_real_, r)
    fold_variance_tables <- vector("list", 2L * r)

    for (component in seq_len(r)) {
      reg_x <- batch_simple_reg(scores$X_preprocessed, scores$Y[, component])
      shr_x <- shrink_variances(
        reg_x$s2,
        method = variance_shrinkage,
        Ymat = scores$X_preprocessed,
        x = scores$Y[, component],
        mr_mashr_args = mr_mashr_args
      )
      ci_x <- build_ci(reg_x$beta1, shr_x$s2_shrunk, reg_x$sxx, reg_x$df, alpha = alpha)
      fold_b_x[, component] <- reg_x$beta1
      fold_se_x[, component] <- ci_x$se
      fold_p_x[, component] <- ci_x$p_value
      fold_df_x[component] <- reg_x$df
      fold_variance_tables[[2L * component - 1L]] <- data.frame(
        fold = fold,
        side = "x",
        component = component,
        feature = colnames(X),
        raw_s2 = reg_x$s2,
        shrunken_s2 = shr_x$s2_shrunk,
        shrinkage_ratio = shr_x$s2_shrunk / pmax(reg_x$s2, .Machine$double.eps),
        variance_shrinkage = variance_shrinkage,
        shrink_backend = shr_x$backend %||% NA_character_,
        shrink_lambda = shr_x$lambda,
        shrink_target = shr_x$mu,
        stringsAsFactors = FALSE
      )

      reg_y <- batch_simple_reg(scores$Y_preprocessed, scores$X[, component])
      shr_y <- shrink_variances(
        reg_y$s2,
        method = variance_shrinkage,
        Ymat = scores$Y_preprocessed,
        x = scores$X[, component],
        mr_mashr_args = mr_mashr_args
      )
      ci_y <- build_ci(reg_y$beta1, shr_y$s2_shrunk, reg_y$sxx, reg_y$df, alpha = alpha)
      fold_b_y[, component] <- reg_y$beta1
      fold_se_y[, component] <- ci_y$se
      fold_p_y[, component] <- ci_y$p_value
      fold_df_y[component] <- reg_y$df
      fold_variance_tables[[2L * component]] <- data.frame(
        fold = fold,
        side = "y",
        component = component,
        feature = colnames(Y),
        raw_s2 = reg_y$s2,
        shrunken_s2 = shr_y$s2_shrunk,
        shrinkage_ratio = shr_y$s2_shrunk / pmax(reg_y$s2, .Machine$double.eps),
        variance_shrinkage = variance_shrinkage,
        shrink_backend = shr_y$backend %||% NA_character_,
        shrink_lambda = shr_y$lambda,
        shrink_target = shr_y$mu,
        stringsAsFactors = FALSE
      )
    }

    pvalues_x[[fold]] <- fold_p_x
    betas_x[[fold]] <- fold_b_x
    ses_x[[fold]] <- fold_se_x
    df_x[[fold]] <- fold_df_x
    pvalues_y[[fold]] <- fold_p_y
    betas_y[[fold]] <- fold_b_y
    ses_y[[fold]] <- fold_se_y
    df_y[[fold]] <- fold_df_y
    variance_diagnostics[[fold]] <- do.call(rbind, fold_variance_tables)

    if (progress) {
      message(sprintf("Cross-fit fold %d/%d complete.", fold, K))
    }
  }

  valid_folds <- which(vapply(pvalues_x, Negate(is.null), logical(1)))
  if (!length(valid_folds)) {
    stop("Cross-fitted inference failed on every fold.", call. = FALSE)
  }

  pvalues_x <- pvalues_x[valid_folds]
  betas_x <- betas_x[valid_folds]
  ses_x <- ses_x[valid_folds]
  df_x <- df_x[valid_folds]
  pvalues_y <- pvalues_y[valid_folds]
  betas_y <- betas_y[valid_folds]
  ses_y <- ses_y[valid_folds]
  df_y <- df_y[valid_folds]
  alignment_diagnostics <- alignment_diagnostics[valid_folds]
  variance_diagnostics <- variance_diagnostics[valid_folds]

  aggregated_p_x <- aggregate_pvalues_cauchy(pvalues_x)
  aggregated_p_y <- aggregate_pvalues_cauchy(pvalues_y)
  avg_beta_x <- average_matrices(betas_x)
  avg_beta_y <- average_matrices(betas_y)
  cct_ci_x <- invert_cauchy_ci_matrix(betas_x, ses_x, df_x, alpha = alpha)
  cct_ci_y <- invert_cauchy_ci_matrix(betas_y, ses_y, df_y, alpha = alpha)

  adjusted_p_x <- aggregated_p_x
  adjusted_p_y <- aggregated_p_y
  for (component in seq_len(r)) {
    adjusted_p_x[, component] <- stats::p.adjust(aggregated_p_x[, component], method = "BH")
    adjusted_p_y[, component] <- stats::p.adjust(aggregated_p_y[, component], method = "BH")
  }

  list(
    x_results = build_inference_table(
      feature_names = colnames(X),
      beta = avg_beta_x,
      lower = cct_ci_x$lower,
      upper = cct_ci_x$upper,
      pval = aggregated_p_x,
      padj = adjusted_p_x,
      alpha = alpha,
      side = "x",
      variance_shrinkage = variance_shrinkage
    ),
    y_results = build_inference_table(
      feature_names = colnames(Y),
      beta = avg_beta_y,
      lower = cct_ci_y$lower,
      upper = cct_ci_y$upper,
      pval = aggregated_p_y,
      padj = adjusted_p_y,
      alpha = alpha,
      side = "y",
      variance_shrinkage = variance_shrinkage
    ),
    fold_correlations = fold_correlations[valid_folds, , drop = FALSE],
    alignment_diagnostics = do.call(rbind, alignment_diagnostics),
    variance_diagnostics = do.call(rbind, variance_diagnostics),
    preprocess_diagnostics = preprocessed$diagnostics,
    folds = folds,
    fit_mode = fit_mode,
    variance_shrinkage = variance_shrinkage,
    ci_method = "cauchy_inversion",
    K = length(valid_folds),
    alpha = alpha
  )
}

extract_crossfit_features <- function(crossfit_result,
                                      side = c("x", "y"),
                                      component = 1,
                                      top_n = 10,
                                      significant_only = TRUE) {
  side <- match.arg(side)
  table <- if (side == "x") crossfit_result$x_results else crossfit_result$y_results
  table <- table[table$component == component, , drop = FALSE]
  table <- table[order(table$p_adjusted, -abs(table$beta), table$p_value), , drop = FALSE]

  if (significant_only) {
    significant <- table[table$significant, , drop = FALSE]
    if (nrow(significant) > 0) {
      table <- significant
    }
  }

  utils::head(table, top_n)
}

#' Compare bootstrap and cross-fitted uncertainty summaries
#'
#' @param reference_fit A fitted `cca_uq_fit` object.
#' @param bootstrap_result Output from [bootstrap_cca_uq()].
#' @param crossfit_result Output from one of the cross-fitted inference
#'   routines.
#' @param component Component to highlight in the returned top-feature tables.
#' @param top_n Number of top rows to retain in the component summaries.
#' @param alpha Two-sided error level used for the bootstrap correlation
#'   summary.
#'
#' @return A named list combining the main bootstrap and cross-fitted summaries
#'   for quick inspection.
#' @export
compare_uq_methods <- function(reference_fit,
                               bootstrap_result,
                               crossfit_result,
                               component = 1,
                               top_n = 10,
                               alpha = 0.05) {
  bootstrap_cor <- summarize_bootstrap_vector(bootstrap_result$cor_matrix, alpha = alpha)

  list(
    fit_summary = data.frame(
      component = seq_along(reference_fit$cor),
      canonical_correlation = reference_fit$cor,
      lambda = reference_fit$lambda,
      stringsAsFactors = FALSE
    ),
    bootstrap_correlations = data.frame(
      component = seq_along(reference_fit$cor),
      estimate = reference_fit$cor,
      boot_mean = bootstrap_cor$mean,
      boot_sd = bootstrap_cor$sd,
      ci_lower = bootstrap_cor$lower,
      ci_upper = bootstrap_cor$upper,
      stringsAsFactors = FALSE
    ),
    crossfit_correlations = data.frame(
      component = seq_len(ncol(crossfit_result$fold_correlations)),
      fold_mean = colMeans(crossfit_result$fold_correlations, na.rm = TRUE),
      fold_sd = apply(crossfit_result$fold_correlations, 2, stats::sd, na.rm = TRUE),
      stringsAsFactors = FALSE
    ),
    bootstrap_x = build_bootstrap_loading_table(
      bootstrap_result,
      side = "x",
      component = component,
      top_n = top_n,
      alpha = alpha
    ),
    bootstrap_y = build_bootstrap_loading_table(
      bootstrap_result,
      side = "y",
      component = component,
      top_n = top_n,
      alpha = alpha
    ),
    crossfit_x = extract_crossfit_features(
      crossfit_result,
      side = "x",
      component = component,
      top_n = top_n,
      significant_only = TRUE
    ),
    crossfit_y = extract_crossfit_features(
      crossfit_result,
      side = "y",
      component = component,
      top_n = top_n,
      significant_only = TRUE
    )
  )
}

#' Run the alcohol example end-to-end
#'
#' Convenience wrapper that loads the alcohol data, fits sparse CCA, runs the
#' bootstrap and cross-fitted UQ procedures, and assembles a side-by-side
#' summary ready for reporting.
#'
#' @inheritParams fit_ecca_cv
#' @inheritParams bootstrap_cca_uq
#' @inheritParams crossfit_cca_inference
#' @param n_boot Number of bootstrap replicates.
#' @param bootstrap_refit_mode Bootstrap refitting mode passed to
#'   [bootstrap_cca_uq()].
#' @param crossfit_fit_mode Foldwise refitting mode passed to the cross-fitted
#'   inference routine.
#' @param crossfit_variance_shrinkage Residual-variance shrinkage rule used by
#'   [crossfit_cca_inference()].
#' @param data_path Optional path to an `.rda` file containing `alcohol`.
#'
#' @return A list containing the loaded data, fitted sparse CCA model,
#'   bootstrap output, cross-fitted output, and a comparison summary.
#' @export
run_alcohol_uq_comparison <- function(r = 2,
                                      lambdas = 10^seq(-2, 0, length.out = 6),
                                      kfolds = 5,
                                      n_boot = 25,
                                      K = 5,
                                      seed = 2025,
                                      preprocess_mode = c("scale", "center", "none"),
                                      bootstrap_refit_mode = c("fixed_lambda", "retune_cv"),
                                      crossfit_fit_mode = c("fixed_lambda", "retune_cv"),
                                      crossfit_variance_shrinkage = c("james-stein", "mr.mashr", "none"),
                                      mr_mashr_args = list(),
                                      mr_ash_args = NULL,
                                      data_path = NULL,
                                      ccar3_path = NULL,
                                      prefer_source = TRUE,
                                      parallelize = FALSE,
                                      nb_cores = NULL,
                                      alpha = 0.05,
                                      verbose = TRUE) {
  preprocess_mode <- match.arg(preprocess_mode)
  bootstrap_refit_mode <- match.arg(bootstrap_refit_mode)
  crossfit_fit_mode <- match.arg(crossfit_fit_mode)
  crossfit_variance_shrinkage <- match_variance_shrinkage(crossfit_variance_shrinkage)
  if (is.null(mr_ash_args)) {
    mr_ash_args <- mr_mashr_args
  } else if (!length(mr_mashr_args)) {
    mr_mashr_args <- mr_ash_args
  } else {
    mr_mashr_args <- modifyList(mr_ash_args, mr_mashr_args, keep.null = TRUE)
  }

  data <- load_alcohol_example(data_path = data_path)
  api <- get_ccar3_api(
    ccar3_path = ccar3_path,
    prefer_source = prefer_source,
    quiet = !isTRUE(verbose)
  )

  if (verbose) {
    message(sprintf(
      "Running alcohol example with n=%d, p=%d, q=%d.",
      data$n, data$p, data$q
    ))
  }

  fit <- fit_ecca_cv(
    X = data$X,
    Y = data$Y,
    r = r,
    lambdas = lambdas,
    kfolds = kfolds,
    preprocess_mode = preprocess_mode,
    ccar3_api = api,
    ccar3_path = ccar3_path,
    prefer_source = prefer_source,
    parallelize = parallelize,
    nb_cores = nb_cores,
    verbose = FALSE
  )

  bootstrap <- bootstrap_cca_uq(
    X = data$X,
    Y = data$Y,
    reference_fit = fit,
    B = n_boot,
    seed = seed,
    refit_mode = bootstrap_refit_mode,
    ccar3_api = api,
    progress = verbose
  )

  crossfit <- crossfit_cca_inference(
    X = data$X,
    Y = data$Y,
    reference_fit = fit,
    K = K,
    seed = seed + 1,
    alpha = alpha,
    fit_mode = crossfit_fit_mode,
    ccar3_api = api,
    progress = verbose,
    variance_shrinkage = crossfit_variance_shrinkage,
    mr_mashr_args = mr_mashr_args
  )

  comparison <- compare_uq_methods(
    reference_fit = fit,
    bootstrap_result = bootstrap,
    crossfit_result = crossfit,
    component = 1,
    top_n = 10,
    alpha = alpha
  )

  list(
    data = data,
    fit = fit,
    bootstrap = bootstrap,
    crossfit = crossfit,
    comparison = comparison
  )
}

#' Print a concise summary for [run_alcohol_uq_comparison()]
#'
#' @param results Output from [run_alcohol_uq_comparison()].
#' @param top_n Number of top features to print for each component summary.
#' @param digits Number of digits used when rounding numeric output.
#'
#' @return Invisibly returns `results`.
#' @export
print_alcohol_uq_summary <- function(results, top_n = 10, digits = 3) {
  comparison <- results$comparison
  format_table <- function(x) {
    x <- as.data.frame(x)
    numeric_cols <- vapply(x, is.numeric, logical(1))
    x[numeric_cols] <- lapply(x[numeric_cols], round, digits = digits)
    x
  }

  cat("Sparse CCA on the alcohol example via ccar3::ecca.cv\n")
  cat(sprintf(
    "Data source: %s | n = %d, p = %d, q = %d\n",
    results$data$source,
    results$data$n,
    results$data$p,
    results$data$q
  ))
  cat(sprintf("Selected lambda: %s\n\n", signif(results$fit$lambda, digits)))
  cat(sprintf(
    "Cross-fit variance shrinkage: %s\n\n",
    variance_shrinkage_label(results$crossfit$variance_shrinkage %||% "james-stein")
  ))

  cat("Canonical correlations\n")
  print(format_table(comparison$bootstrap_correlations))
  cat("\n")

  cat("Cross-fit fold correlations\n")
  print(format_table(comparison$crossfit_correlations))
  cat("\n")

  cat("Cross-fit preprocessing diagnostics\n")
  print(format_table(results$crossfit$preprocess_diagnostics))
  cat("\n")

  cat("Cross-fit alignment diagnostics\n")
  print(format_table(results$crossfit$alignment_diagnostics))
  cat("\n")

  cat(sprintf("Top bootstrap-stable gene loadings for component 1 (top %d)\n", top_n))
  print(format_table(comparison$bootstrap_x[seq_len(min(top_n, nrow(comparison$bootstrap_x))), ]))
  cat("\n")

  cat(sprintf("Top bootstrap-stable methylation loadings for component 1 (top %d)\n", top_n))
  print(format_table(comparison$bootstrap_y[seq_len(min(top_n, nrow(comparison$bootstrap_y))), ]))
  cat("\n")

  cat(sprintf("Top cross-fit gene findings for component 1 (top %d)\n", top_n))
  print(format_table(comparison$crossfit_x[seq_len(min(top_n, nrow(comparison$crossfit_x))), ]))
  cat("\n")

  cat(sprintf("Top cross-fit methylation findings for component 1 (top %d)\n", top_n))
  print(format_table(comparison$crossfit_y[seq_len(min(top_n, nrow(comparison$crossfit_y))), ]))

  invisible(results)
}
