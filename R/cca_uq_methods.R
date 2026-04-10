`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.uqcca_cache <- new.env(parent = emptyenv())

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

get_ccar3_api <- function(ccar3_path = "/Users/clairedonnat/Documents/ccar3",
                          prefer_source = TRUE,
                          quiet = TRUE) {
  normalized_path <- normalizePath(ccar3_path, mustWork = FALSE)
  cache_key <- paste(normalized_path, prefer_source, sep = "::")

  if (exists(cache_key, envir = .uqcca_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .uqcca_cache, inherits = FALSE))
  }

  api <- NULL

  if (prefer_source && dir.exists(normalized_path)) {
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
        "Could not load 'ccar3'. Install the package or provide a readable `ccar3_path`.",
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

#' Fit sparse CCA with `ccar3::ecca.cv`.
fit_ecca_cv <- function(X, Y,
                        r = 2,
                        lambdas = 10^seq(-2, 0, length.out = 8),
                        kfolds = 5,
                        preprocess_mode = c("scale", "center", "none"),
                        ccar3_api = NULL,
                        ccar3_path = "/Users/clairedonnat/Documents/ccar3",
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

#' Refit sparse CCA with a fixed penalty using `ccar3::ecca`.
fit_ecca_fixed_lambda <- function(X, Y,
                                  lambda,
                                  r = 2,
                                  preprocess_mode = c("scale", "center", "none"),
                                  ccar3_api = NULL,
                                  ccar3_path = "/Users/clairedonnat/Documents/ccar3",
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

#' Backward-compatible alias for the previous CV fitter name.
fit_cca_rrr_cv <- function(...) {
  fit_ecca_cv(...)
}

#' Backward-compatible alias for the previous fixed-lambda fitter name.
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

#' Bootstrap uncertainty quantification for sparse CCA loadings.
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
      stringsAsFactors = FALSE
    )
  })

  table <- do.call(rbind, tables)
  rownames(table) <- NULL
  table[order(table$component, -abs(table$estimate)), , drop = FALSE]
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
  variance_csv <- file.path(output_dir, "crossfit_variance_shrinkage.csv")
  alignment_csv <- file.path(output_dir, "crossfit_alignment_diagnostics.csv")
  preprocess_csv <- file.path(output_dir, "crossfit_preprocess_diagnostics.csv")
  utils::write.csv(gene_table, gene_csv, row.names = FALSE)
  utils::write.csv(methylation_table, methylation_csv, row.names = FALSE)
  utils::write.csv(variance_table, variance_csv, row.names = FALSE)
  utils::write.csv(alignment_table, alignment_csv, row.names = FALSE)
  utils::write.csv(preprocess_table, preprocess_csv, row.names = FALSE)

  gene_pdf <- file.path(output_dir, "crossfit_gene_regression_cis.pdf")
  methylation_pdf <- file.path(output_dir, "crossfit_methylation_regression_cis.pdf")
  variance_pdf <- file.path(output_dir, "crossfit_variance_shrinkage.pdf")

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
  plot_crossfit_variance_diagnostics(
    variance_table = variance_table,
    output_path = variance_pdf,
    title_prefix = sprintf(
      "Alcohol example cross-fit variance shrinkage diagnostic (%s)",
      shrinkage_label
    )
  )

  if (verbose) {
    message(sprintf("Saved full cross-fit CI tables to %s and %s.", gene_csv, methylation_csv))
    message(sprintf("Saved cross-fit variance diagnostics to %s.", variance_csv))
    message(sprintf("Saved cross-fit alignment diagnostics to %s.", alignment_csv))
    message(sprintf("Saved cross-fit preprocess diagnostics to %s.", preprocess_csv))
    message(sprintf("Saved cross-fit CI plots to %s and %s.", gene_pdf, methylation_pdf))
    message(sprintf("Saved cross-fit variance shrinkage plot to %s.", variance_pdf))
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
                                  variance_shrinkage = NULL) {
  tables <- lapply(seq_len(ncol(beta)), function(component) {
    table <- data.frame(
      side = side,
      feature = feature_names,
      component = component,
      beta = beta[, component],
      ci_lower = lower[, component],
      ci_upper = upper[, component],
      p_value = pval[, component],
      p_adjusted = padj[, component],
      significant = padj[, component] < alpha,
      stringsAsFactors = FALSE
    )

    if (!is.null(variance_shrinkage)) {
      table$variance_shrinkage <- variance_shrinkage
    }

    table
  })

  do.call(rbind, tables)
}

#' Cross-fitted inference for sparse CCA loadings.
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

#' Compare bootstrap and cross-fitted uncertainty quantification outputs.
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

#' Run the alcohol example with bootstrap and cross-fitted UQ.
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
                                      ccar3_path = "/Users/clairedonnat/Documents/ccar3",
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
