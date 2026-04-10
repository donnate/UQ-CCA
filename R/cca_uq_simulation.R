.uqcca_sim_cache <- new.env(parent = emptyenv())

rbind_fill <- function(data_frames) {
  data_frames <- data_frames[!vapply(data_frames, is.null, logical(1))]
  if (!length(data_frames)) {
    return(data.frame())
  }

  all_names <- unique(unlist(lapply(data_frames, names), use.names = FALSE))
  aligned <- lapply(data_frames, function(df) {
    missing <- setdiff(all_names, names(df))
    for (name in missing) {
      df[[name]] <- NA
    }
    df[all_names]
  })

  do.call(rbind, aligned)
}

get_rrr_simulation_api <- function(ccar3_code_path = "/Users/clairedonnat/Documents/CCAR3_code") {
  ccar3_code_path <- normalizePath(ccar3_code_path, mustWork = TRUE)
  sim_file <- normalizePath(
    file.path(ccar3_code_path, "experiments", "simulations", "generate_example_rrr.R"),
    mustWork = TRUE
  )

  if (exists(sim_file, envir = .uqcca_sim_cache, inherits = FALSE)) {
    return(get(sim_file, envir = .uqcca_sim_cache, inherits = FALSE))
  }

  env <- new.env(parent = globalenv())
  previous_wd <- getwd()
  on.exit(setwd(previous_wd), add = TRUE)
  setwd(ccar3_code_path)
  sys.source(sim_file, envir = env, chdir = FALSE)
  assign(sim_file, env, envir = .uqcca_sim_cache)
  env
}

resolve_theta_matrix <- function(theta_strength = "high", r = 2) {
  if (is.matrix(theta_strength)) {
    if (!all(dim(theta_strength) == c(r, r))) {
      stop("`theta_strength` matrix must be r x r.", call. = FALSE)
    }
    return(theta_strength)
  }

  if (is.numeric(theta_strength)) {
    if (length(theta_strength) == 1L) {
      return(diag(rep(theta_strength, r), nrow = r))
    }
    if (length(theta_strength) == r) {
      return(diag(theta_strength, nrow = r))
    }
    stop("Numeric `theta_strength` must have length 1 or r.", call. = FALSE)
  }

  theta_strength <- match.arg(theta_strength, c("high", "medium", "low"))
  if (theta_strength == "high") {
    return(diag(seq(0.9, 0.75, length.out = r), nrow = r))
  }
  if (theta_strength == "medium") {
    return(diag(seq(0.7, 0.55, length.out = r), nrow = r))
  }

  diag(seq(0.5, 0.35, length.out = r), nrow = r)
}

generate_rrr_simulation_data <- function(n = 100,
                                         p = 100,
                                         q = 20,
                                         r = 2,
                                         nnzeros = 10,
                                         theta_strength = "high",
                                         r_pca = 5,
                                         lambda_pca = 1,
                                         overlapping_amount = 1,
                                         normalize_diagonal = TRUE,
                                         n_new = 5000,
                                         ccar3_code_path = "/Users/clairedonnat/Documents/CCAR3_code",
                                         generator_api = NULL) {
  api <- generator_api %||% get_rrr_simulation_api(ccar3_code_path = ccar3_code_path)
  theta <- resolve_theta_matrix(theta_strength = theta_strength, r = r)

  sim <- api$generate_example_sparse_U(
    n = n,
    p1 = p,
    p2 = q,
    r_pca = r_pca,
    nnzeros = nnzeros,
    theta = theta,
    lambda_pca = lambda_pca,
    r = r,
    overlapping_amount = overlapping_amount,
    normalize_diagonal = normalize_diagonal,
    n_new = n_new
  )

  colnames(sim$X) <- paste0("X_", seq_len(ncol(sim$X)))
  colnames(sim$Y) <- paste0("Y_", seq_len(ncol(sim$Y)))
  colnames(sim$Xnew) <- colnames(sim$X)
  colnames(sim$Ynew) <- colnames(sim$Y)
  rownames(sim$u) <- colnames(sim$X)
  rownames(sim$v) <- colnames(sim$Y)

  truth_beta_x <- sim$Sigmax %*% sim$u %*% theta
  truth_beta_y <- sim$Sigmay %*% sim$v %*% theta
  rownames(truth_beta_x) <- colnames(sim$X)
  rownames(truth_beta_y) <- colnames(sim$Y)

  list(
    X = sim$X,
    Y = sim$Y,
    Xnew = sim$Xnew,
    Ynew = sim$Ynew,
    truth = list(
      U = sim$u,
      V = sim$v,
      beta_x = truth_beta_x,
      beta_y = truth_beta_y,
      theta = theta
    ),
    generator = sim,
    settings = list(
      n = n,
      p = p,
      q = q,
      r = r,
      nnzeros = nnzeros,
      theta_strength = theta_strength,
      theta = theta,
      r_pca = r_pca,
      lambda_pca = lambda_pca,
      overlapping_amount = overlapping_amount,
      normalize_diagonal = normalize_diagonal,
      n_new = n_new,
      ccar3_code_path = ccar3_code_path
    )
  )
}

build_truth_coefficient_table <- function(truth_matrix,
                                          side = c("x", "y"),
                                          parameter_type = c("canonical_direction", "regression_coefficient"),
                                          truth_tol = 1e-8) {
  side <- match.arg(side)
  parameter_type <- match.arg(parameter_type)
  feature_names <- rownames(truth_matrix) %||%
    paste0(if (side == "x") "X_" else "Y_", seq_len(nrow(truth_matrix)))

  direction_prefix <- if (parameter_type == "canonical_direction") {
    if (side == "x") "U" else "V"
  } else if (side == "x") {
    "RX"
  } else {
    "RY"
  }

  tables <- lapply(seq_len(ncol(truth_matrix)), function(component) {
    values <- truth_matrix[, component]
    data.frame(
      side = side,
      parameter_type = parameter_type,
      feature = feature_names,
      component = component,
      direction = paste0(direction_prefix, component),
      truth_value = values,
      abs_truth = abs(values),
      true_nonzero = abs(values) > truth_tol,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, tables)
}

ci_detects_signal <- function(lower, upper) {
  is.finite(lower) & is.finite(upper) & (lower > 0 | upper < 0)
}

combine_bootstrap_ci_tables <- function(bootstrap_result,
                                        alpha = 0.05,
                                        threshold = 1e-6) {
  x_table <- build_bootstrap_loading_table_all(
    bootstrap_result = bootstrap_result,
    side = "x",
    alpha = alpha,
    threshold = threshold
  )
  y_table <- build_bootstrap_loading_table_all(
    bootstrap_result = bootstrap_result,
    side = "y",
    alpha = alpha,
    threshold = threshold
  )

  combined <- rbind(x_table, y_table)
  combined$parameter_type <- "canonical_direction"
  combined$method <- "bootstrap"
  combined$estimate_value <- combined$estimate
  combined$discovered <- ci_detects_signal(combined$ci_lower, combined$ci_upper)
  combined
}

combine_crossfit_ci_tables <- function(crossfit_result) {
  combined <- rbind(crossfit_result$x_results, crossfit_result$y_results)
  combined$parameter_type <- "regression_coefficient"
  combined$method <- "crossfit"
  combined$estimate_value <- combined$beta
  combined$discovered <- ci_detects_signal(combined$ci_lower, combined$ci_upper)
  combined
}

evaluate_ci_outcomes <- function(estimated_table,
                                 truth_table,
                                 replication,
                                 settings = list()) {
  merged <- merge(
    estimated_table,
    truth_table[, c("side", "feature", "component", "truth_value", "abs_truth", "true_nonzero")],
    by = c("side", "feature", "component"),
    all.x = TRUE,
    sort = FALSE
  )

  merged$replication <- replication
  merged$true_positive <- merged$discovered & merged$true_nonzero
  merged$false_positive <- merged$discovered & !merged$true_nonzero
  merged$false_negative <- (!merged$discovered) & merged$true_nonzero

  if (length(settings)) {
    for (name in names(settings)) {
      merged[[name]] <- settings[[name]]
    }
  }

  merged
}

summarize_ci_outcomes_by <- function(outcomes, group_cols) {
  split_index <- interaction(outcomes[group_cols], drop = TRUE, lex.order = TRUE)
  groups <- split(outcomes, split_index)

  summaries <- lapply(groups, function(df) {
    n_nonzero <- sum(df$true_nonzero, na.rm = TRUE)
    n_zero <- sum(!df$true_nonzero, na.rm = TRUE)
    discoveries <- sum(df$discovered, na.rm = TRUE)
    true_positives <- sum(df$true_positive, na.rm = TRUE)
    false_positives <- sum(df$false_positive, na.rm = TRUE)

    summary_row <- df[1, group_cols, drop = FALSE]
    summary_row$n_total <- nrow(df)
    summary_row$n_nonzero <- n_nonzero
    summary_row$n_zero <- n_zero
    summary_row$discoveries <- discoveries
    summary_row$true_positives <- true_positives
    summary_row$false_positives <- false_positives
    summary_row$power <- if (n_nonzero > 0) true_positives / n_nonzero else NA_real_
    summary_row$false_discovery_proportion <- if (discoveries > 0) {
      false_positives / discoveries
    } else {
      NA_real_
    }
    summary_row$false_positive_rate <- if (n_zero > 0) false_positives / n_zero else NA_real_
    summary_row
  })

  out <- do.call(rbind, summaries)
  rownames(out) <- NULL
  out
}

assign_strength_bins <- function(abs_truth, n_bins = 6) {
  abs_truth <- abs_truth[is.finite(abs_truth)]
  if (!length(abs_truth)) {
    return(list(labels = character(), left = numeric(), right = numeric()))
  }

  value_range <- range(abs_truth)
  if (diff(value_range) <= .Machine$double.eps) {
    label <- sprintf("[%.3f, %.3f]", value_range[1], value_range[2])
    return(list(
      labels = rep(label, length(abs_truth)),
      left = rep(value_range[1], length(abs_truth)),
      right = rep(value_range[2], length(abs_truth))
    ))
  }

  breaks <- seq(value_range[1], value_range[2], length.out = n_bins + 1L)
  breaks[1] <- breaks[1] - 1e-12
  breaks[length(breaks)] <- breaks[length(breaks)] + 1e-12

  bins <- cut(abs_truth, breaks = breaks, include.lowest = TRUE, right = TRUE)
  left <- breaks[as.integer(bins)]
  right <- breaks[as.integer(bins) + 1L]
  labels <- sprintf("[%.3f, %.3f]", left, right)

  list(labels = labels, left = left, right = right)
}

summarize_power_by_strength <- function(outcomes, n_bins = 6) {
  signal_only <- outcomes[outcomes$true_nonzero, , drop = FALSE]
  if (!nrow(signal_only)) {
    return(data.frame())
  }

  group_cols <- c("method", "parameter_type", "side", "component")
  split_index <- interaction(signal_only[group_cols], drop = TRUE, lex.order = TRUE)
  groups <- split(signal_only, split_index)

  summaries <- lapply(groups, function(df) {
    bins <- assign_strength_bins(df$abs_truth, n_bins = n_bins)
    df$strength_bin <- bins$labels
    df$bin_left <- bins$left
    df$bin_right <- bins$right

    bin_split <- split(df, df$strength_bin)
    do.call(rbind, lapply(bin_split, function(bin_df) {
      out <- bin_df[1, group_cols, drop = FALSE]
      out$strength_bin <- bin_df$strength_bin[1]
      out$bin_left <- bin_df$bin_left[1]
      out$bin_right <- bin_df$bin_right[1]
      out$bin_midpoint <- mean(c(out$bin_left, out$bin_right))
      out$n_in_bin <- nrow(bin_df)
      out$power <- mean(bin_df$discovered, na.rm = TRUE)
      out
    }))
  })

  power_df <- do.call(rbind, summaries)
  rownames(power_df) <- NULL
  power_df[order(power_df$parameter_type, power_df$side, power_df$component, power_df$bin_left), , drop = FALSE]
}

build_subspace_recovery_table <- function(fit,
                                          truth,
                                          replication,
                                          settings = list()) {
  tables <- lapply(
    c("x", "y"),
    function(side) {
      estimate <- if (side == "x") fit$U else fit$V
      reference <- if (side == "x") truth$U else truth$V
      rank <- ncol(reference)
      max_distance <- sqrt(2 * rank)
      distance <- compute_subspace_distance(estimate, reference)

      data.frame(
        side = side,
        rank = rank,
        subspace_distance = distance,
        normalized_subspace_distance = distance / max_distance,
        subspace_similarity = 1 - distance / max_distance,
        replication = replication,
        stringsAsFactors = FALSE
      )
    }
  )

  out <- do.call(rbind, tables)
  if (length(settings)) {
    for (name in names(settings)) {
      out[[name]] <- settings[[name]]
    }
  }

  rownames(out) <- NULL
  out
}

summarize_subspace_recovery_by <- function(subspace_table, group_cols) {
  if (!nrow(subspace_table)) {
    return(data.frame())
  }

  split_index <- interaction(subspace_table[group_cols], drop = TRUE, lex.order = TRUE)
  groups <- split(subspace_table, split_index)

  summaries <- lapply(groups, function(df) {
    summary_row <- df[1, group_cols, drop = FALSE]
    summary_row$n_replications <- nrow(df)
    summary_row$mean_subspace_distance <- mean(df$subspace_distance, na.rm = TRUE)
    summary_row$sd_subspace_distance <- stats::sd(df$subspace_distance, na.rm = TRUE)
    summary_row$median_subspace_distance <- stats::median(df$subspace_distance, na.rm = TRUE)
    summary_row$min_subspace_distance <- min(df$subspace_distance, na.rm = TRUE)
    summary_row$max_subspace_distance <- max(df$subspace_distance, na.rm = TRUE)
    summary_row$mean_normalized_subspace_distance <- mean(df$normalized_subspace_distance, na.rm = TRUE)
    summary_row$sd_normalized_subspace_distance <- stats::sd(df$normalized_subspace_distance, na.rm = TRUE)
    summary_row$mean_subspace_similarity <- mean(df$subspace_similarity, na.rm = TRUE)
    summary_row
  })

  out <- do.call(rbind, summaries)
  rownames(out) <- NULL
  out
}

plot_simulation_power_fdr <- function(summary_table, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to draw simulation plots.", call. = FALSE)
  }

  plot_df <- rbind(
    transform(summary_table, metric = "Power", value = power),
    transform(summary_table, metric = "FDP", value = false_discovery_proportion),
    transform(summary_table, metric = "FPR", value = false_positive_rate)
  )
  plot_df <- plot_df[is.finite(plot_df$value), , drop = FALSE]

  plot_df$label <- paste(plot_df$method, plot_df$side, sep = " / ")

  plot_object <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = label, y = value, fill = method)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_grid(metric ~ parameter_type, scales = "free_y") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "CI-based power and false discovery summary",
      x = NULL,
      y = "Rate"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(output_path, plot_object, width = 10, height = 7)
  invisible(output_path)
}

plot_power_by_strength <- function(power_by_strength, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to draw simulation plots.", call. = FALSE)
  }

  if (!nrow(power_by_strength)) {
    stop("`power_by_strength` is empty.", call. = FALSE)
  }
  power_by_strength <- power_by_strength[
    is.finite(power_by_strength$bin_midpoint) & is.finite(power_by_strength$power),
    ,
    drop = FALSE
  ]

  power_by_strength$series <- paste(power_by_strength$method, power_by_strength$side, sep = " / ")

  plot_object <- ggplot2::ggplot(
    power_by_strength,
    ggplot2::aes(x = bin_midpoint, y = power, color = series)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_grid(parameter_type + side ~ component, scales = "free_x") +
    ggplot2::labs(
      title = "Power as a function of true coefficient magnitude",
      x = "True coefficient magnitude bin midpoint",
      y = "Power",
      color = "Method / side"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  ggplot2::ggsave(output_path, plot_object, width = 11, height = 7)
  invisible(output_path)
}

plot_subspace_distance_summary <- function(subspace_table, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to draw simulation plots.", call. = FALSE)
  }

  if (!nrow(subspace_table)) {
    stop("`subspace_table` is empty.", call. = FALSE)
  }

  plot_df <- subspace_table
  plot_df$side <- factor(plot_df$side, levels = c("x", "y"), labels = c("X subspace", "Y subspace"))

  plot_object <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = side, y = subspace_distance, fill = side)
  ) +
    ggplot2::geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = NA) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.08, height = 0),
      size = 2,
      alpha = 0.8
    ) +
    ggplot2::labs(
      title = "Estimated-versus-true subspace distance",
      subtitle = "Projection Frobenius distance across simulation replications; smaller is better",
      x = NULL,
      y = "Subspace distance"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(output_path, plot_object, width = 8, height = 5)
  invisible(output_path)
}

run_rrr_uq_simulation <- function(n_sims = 20,
                                  n = 500,
                                  p = 100,
                                  q = 20,
                                  r = 2,
                                  nnzeros = 10,
                                  theta_strength = "high",
                                  r_pca = 5,
                                  lambda_pca = 1,
                                  overlapping_amount = 1,
                                  normalize_diagonal = TRUE,
                                  n_new = 5000,
                                  lambdas = 10^seq(-3, 1, length.out = 30),
                                  kfolds = 5,
                                  n_boot = 25,
                                  K = 2,
                                  alpha = 0.05,
                                  truth_tol = 1e-8,
                                  bootstrap_threshold = 1e-6,
                                  strength_bins = 6,
                                  seed = 1,
                                  preprocess_mode = "center",
                                  crossfit_variance_shrinkage = c("james-stein", "mr.mashr", "none"),
                                  mr_mashr_args = list(),
                                  mr_ash_args = NULL,
                                  ccar3_path = "/Users/clairedonnat/Documents/ccar3",
                                  ccar3_code_path = "/Users/clairedonnat/Documents/CCAR3_code",
                                  prefer_source = TRUE,
                                  parallelize = FALSE,
                                  verbose = TRUE) {
  crossfit_variance_shrinkage <- match_variance_shrinkage(crossfit_variance_shrinkage)
  if (is.null(mr_ash_args)) {
    mr_ash_args <- mr_mashr_args
  } else if (!length(mr_mashr_args)) {
    mr_mashr_args <- mr_ash_args
  } else {
    mr_mashr_args <- modifyList(mr_ash_args, mr_mashr_args, keep.null = TRUE)
  }
  ccar3_api <- get_ccar3_api(
    ccar3_path = ccar3_path,
    prefer_source = prefer_source,
    quiet = !isTRUE(verbose)
  )
  generator_api <- get_rrr_simulation_api(ccar3_code_path = ccar3_code_path)

  outcomes <- vector("list", n_sims)
  subspace_recovery <- vector("list", n_sims)
  simulations <- vector("list", n_sims)

  for (replication in seq_len(n_sims)) {
    if (verbose) {
      message(sprintf("Simulation replicate %d/%d", replication, n_sims))
    }

    set.seed(seed + replication - 1L)
    sim_data <- generate_rrr_simulation_data(
      n = n,
      p = p,
      q = q,
      r = r,
      nnzeros = nnzeros,
      theta_strength = theta_strength,
      r_pca = r_pca,
      lambda_pca = lambda_pca,
      overlapping_amount = overlapping_amount,
      normalize_diagonal = normalize_diagonal,
      n_new = n_new,
      ccar3_code_path = ccar3_code_path,
      generator_api = generator_api
    )

    settings <- list(
      n = n,
      p = p,
      q = q,
      r = r,
      nnzeros = nnzeros,
      theta_strength = if (is.character(theta_strength)) theta_strength else "custom",
      crossfit_variance_shrinkage = crossfit_variance_shrinkage
    )

    fit <- fit_ecca_cv(
      X = sim_data$X,
      Y = sim_data$Y,
      r = r,
      lambdas = lambdas,
      kfolds = kfolds,
      preprocess_mode = preprocess_mode,
      ccar3_api = ccar3_api,
      ccar3_path = ccar3_path,
      prefer_source = prefer_source,
      parallelize = parallelize,
      verbose = FALSE
    )
    fit <- align_fit_to_reference(fit, list(U = sim_data$truth$U, V = sim_data$truth$V))
    subspace_recovery[[replication]] <- build_subspace_recovery_table(
      fit = fit,
      truth = sim_data$truth,
      replication = replication,
      settings = settings
    )

    bootstrap <- bootstrap_cca_uq(
      X = sim_data$X,
      Y = sim_data$Y,
      reference_fit = fit,
      B = n_boot,
      seed = seed * 1000L + replication,
      refit_mode = "fixed_lambda",
      ccar3_api = ccar3_api,
      progress = FALSE
    )
    crossfit <- crossfit_cca_inference(
      X = sim_data$X,
      Y = sim_data$Y,
      reference_fit = fit,
      K = K,
      seed = seed * 2000L + replication,
      alpha = alpha,
      fit_mode = "fixed_lambda",
      ccar3_api = ccar3_api,
      progress = FALSE,
      variance_shrinkage = crossfit_variance_shrinkage,
      mr_mashr_args = mr_mashr_args
    )

    bootstrap_truth <- rbind(
      build_truth_coefficient_table(
        truth_matrix = sim_data$truth$U,
        side = "x",
        parameter_type = "canonical_direction",
        truth_tol = truth_tol
      ),
      build_truth_coefficient_table(
        truth_matrix = sim_data$truth$V,
        side = "y",
        parameter_type = "canonical_direction",
        truth_tol = truth_tol
      )
    )
    crossfit_truth <- rbind(
      build_truth_coefficient_table(
        truth_matrix = sim_data$truth$beta_x,
        side = "x",
        parameter_type = "regression_coefficient",
        truth_tol = truth_tol
      ),
      build_truth_coefficient_table(
        truth_matrix = sim_data$truth$beta_y,
        side = "y",
        parameter_type = "regression_coefficient",
        truth_tol = truth_tol
      )
    )

    bootstrap_outcomes <- evaluate_ci_outcomes(
      estimated_table = combine_bootstrap_ci_tables(
        bootstrap_result = bootstrap,
        alpha = alpha,
        threshold = bootstrap_threshold
      ),
      truth_table = bootstrap_truth,
      replication = replication,
      settings = settings
    )
    crossfit_outcomes <- evaluate_ci_outcomes(
      estimated_table = combine_crossfit_ci_tables(crossfit_result = crossfit),
      truth_table = crossfit_truth,
      replication = replication,
      settings = settings
    )

    outcomes[[replication]] <- rbind_fill(list(bootstrap_outcomes, crossfit_outcomes))
    simulations[[replication]] <- list(
      data = sim_data,
      fit = fit,
      bootstrap = bootstrap,
      crossfit = crossfit
    )
  }

  outcomes_df <- do.call(rbind, outcomes)
  rownames(outcomes_df) <- NULL
  subspace_df <- do.call(rbind, subspace_recovery)
  rownames(subspace_df) <- NULL

  list(
    outcomes = outcomes_df,
    subspace = subspace_df,
    summary = list(
      overall = summarize_ci_outcomes_by(
        outcomes_df,
        group_cols = c("method", "parameter_type", "side")
      ),
      by_component = summarize_ci_outcomes_by(
        outcomes_df,
        group_cols = c("method", "parameter_type", "side", "component")
      ),
      by_replication = summarize_ci_outcomes_by(
        outcomes_df,
        group_cols = c("method", "parameter_type", "side", "replication")
      ),
      subspace = summarize_subspace_recovery_by(
        subspace_df,
        group_cols = c("side")
      )
    ),
    power_by_strength = summarize_power_by_strength(outcomes_df, n_bins = strength_bins),
    simulations = simulations,
    settings = list(
      n_sims = n_sims,
      n = n,
      p = p,
      q = q,
      r = r,
      nnzeros = nnzeros,
      theta_strength = theta_strength,
      r_pca = r_pca,
      lambda_pca = lambda_pca,
      overlapping_amount = overlapping_amount,
      normalize_diagonal = normalize_diagonal,
      n_new = n_new,
      lambdas = lambdas,
      kfolds = kfolds,
      n_boot = n_boot,
      K = K,
      alpha = alpha,
      truth_tol = truth_tol,
      bootstrap_threshold = bootstrap_threshold,
      strength_bins = strength_bins,
      seed = seed,
      preprocess_mode = preprocess_mode,
      ccar3_path = ccar3_path,
      ccar3_code_path = ccar3_code_path,
      prefer_source = prefer_source
    )
  )
}

write_rrr_uq_simulation_reports <- function(simulation_results,
                                            output_dir,
                                            verbose = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  outcomes_csv <- file.path(output_dir, "simulation_ci_outcomes.csv")
  summary_csv <- file.path(output_dir, "simulation_ci_summary.csv")
  component_csv <- file.path(output_dir, "simulation_ci_summary_by_component.csv")
  replication_csv <- file.path(output_dir, "simulation_ci_summary_by_replication.csv")
  strength_csv <- file.path(output_dir, "simulation_power_by_strength.csv")
  subspace_csv <- file.path(output_dir, "simulation_subspace_distance.csv")
  subspace_summary_csv <- file.path(output_dir, "simulation_subspace_distance_summary.csv")

  utils::write.csv(simulation_results$outcomes, outcomes_csv, row.names = FALSE)
  utils::write.csv(simulation_results$summary$overall, summary_csv, row.names = FALSE)
  utils::write.csv(simulation_results$summary$by_component, component_csv, row.names = FALSE)
  utils::write.csv(simulation_results$summary$by_replication, replication_csv, row.names = FALSE)
  utils::write.csv(simulation_results$power_by_strength, strength_csv, row.names = FALSE)
  utils::write.csv(simulation_results$subspace, subspace_csv, row.names = FALSE)
  utils::write.csv(simulation_results$summary$subspace, subspace_summary_csv, row.names = FALSE)

  summary_pdf <- file.path(output_dir, "simulation_power_false_discovery_summary.pdf")
  strength_pdf <- file.path(output_dir, "simulation_power_by_strength.pdf")
  subspace_pdf <- file.path(output_dir, "simulation_subspace_distance.pdf")
  plot_simulation_power_fdr(simulation_results$summary$overall, summary_pdf)
  plot_power_by_strength(simulation_results$power_by_strength, strength_pdf)
  plot_subspace_distance_summary(simulation_results$subspace, subspace_pdf)

  if (verbose) {
    message(sprintf("Saved simulation CI outcomes to %s.", outcomes_csv))
    message(sprintf("Saved simulation summaries to %s, %s, and %s.", summary_csv, component_csv, replication_csv))
    message(sprintf("Saved power-by-strength data to %s.", strength_csv))
    message(sprintf("Saved subspace-distance data to %s and %s.", subspace_csv, subspace_summary_csv))
    message(sprintf("Saved simulation plots to %s, %s, and %s.", summary_pdf, strength_pdf, subspace_pdf))
  }

  list(
    outcomes_csv = outcomes_csv,
    summary_csv = summary_csv,
    component_csv = component_csv,
    replication_csv = replication_csv,
    strength_csv = strength_csv,
    subspace_csv = subspace_csv,
    subspace_summary_csv = subspace_summary_csv,
    summary_pdf = summary_pdf,
    strength_pdf = strength_pdf,
    subspace_pdf = subspace_pdf
  )
}
