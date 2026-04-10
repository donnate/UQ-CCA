setwd("~/Documents/UQ-CCA/")
find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)

  repeat {
    helper_path <- file.path(path, "R", "cca_uq_methods.R")
    if (file.exists(helper_path)) {
      return(path)
    }

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find the project root containing R/cca_uq_methods.R.", call. = FALSE)
    }
    path <- parent
  }
}

parse_integer_values <- function(value, default) {
  if (missing(value) || !nzchar(value)) {
    return(as.integer(default))
  }

  pieces <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  pieces <- pieces[nzchar(pieces)]
  values <- suppressWarnings(as.integer(pieces))

  if (!length(values) || anyNA(values) || any(values <= 0L)) {
    stop("`UQCCA_SIM_Q_VALUES` must be a comma-separated list of positive integers.", call. = FALSE)
  }

  unique(values)
}

plot_q_sweep_subspace_summary <- function(subspace_summary, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to draw q-sweep plots.", call. = FALSE)
  }

  if (!nrow(subspace_summary)) {
    stop("`subspace_summary` is empty.", call. = FALSE)
  }

  plot_df <- subspace_summary
  plot_df$side <- factor(plot_df$side, levels = c("x", "y"), labels = c("X subspace", "Y subspace"))
  side_has_path <- vapply(
    split(plot_df$q, plot_df$side),
    function(values) length(unique(values)) > 1L,
    logical(1)
  )

  plot_object <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = q, y = mean_subspace_distance, color = side)
  )
  if (any(side_has_path)) {
    plot_object <- plot_object + ggplot2::geom_line(linewidth = 0.8)
  }

  plot_object <- plot_object +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = pmax(mean_subspace_distance - sd_subspace_distance, 0),
        ymax = mean_subspace_distance + sd_subspace_distance
      ),
      width = 0.8
    ) +
    ggplot2::scale_x_continuous(breaks = sort(unique(plot_df$q))) +
    ggplot2::labs(
      title = "Subspace distance across q values",
      subtitle = "Error bars show +/-1 SD across replications; smaller is better",
      x = "q",
      y = "Mean subspace distance",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  ggplot2::ggsave(output_path, plot_object, width = 8, height = 5)
  invisible(output_path)
}

project_root <- find_project_root()
source(file.path(project_root, "R", "cca_uq_methods.R"))
source(file.path(project_root, "R", "cca_uq_simulation.R"))

ccar3_path <- Sys.getenv("UQCCA_CCAR3_PATH", unset = "/Users/clairedonnat/Documents/ccar3")
ccar3_code_path <- Sys.getenv("UQCCA_CCAR3_CODE_PATH", unset = "/Users/clairedonnat/Documents/CCAR3_code")
output_dir <- Sys.getenv(
  "UQCCA_SIM_Q_SWEEP_OUTPUT_DIR",
  unset = file.path(project_root, "results", "simulation_rrr_uq_by_q")
)

n_sims <- as.integer(Sys.getenv("UQCCA_SIM_N_REPS", unset = "10"))
n <- as.integer(Sys.getenv("UQCCA_SIM_N", unset = "500"))
p <- as.integer(Sys.getenv("UQCCA_SIM_P", unset = "100"))
q_values <- parse_integer_values(Sys.getenv("UQCCA_SIM_Q_VALUES", unset = "10,20,40,80"), c(10L, 20L, 40L, 80L))
r <- as.integer(Sys.getenv("UQCCA_SIM_R", unset = "2"))
nnzeros <- as.integer(Sys.getenv("UQCCA_SIM_NNZEROS", unset = "10"))
r_pca <- as.integer(Sys.getenv("UQCCA_SIM_R_PCA", unset = "5"))
n_boot <- as.integer(Sys.getenv("UQCCA_SIM_BOOT", unset = "10"))
k_outer <- as.integer(Sys.getenv("UQCCA_SIM_K", unset = "2"))
kfolds <- as.integer(Sys.getenv("UQCCA_SIM_KFOLDS", unset = "3"))
strength_bins <- as.integer(Sys.getenv("UQCCA_SIM_STRENGTH_BINS", unset = "6"))
theta_strength <- Sys.getenv("UQCCA_SIM_THETA_STRENGTH", unset = "high")
seed <- as.integer(Sys.getenv("UQCCA_SIM_SEED", unset = "2025"))
lambdas <- 10^seq(-3, 0, length.out = 10)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

simulation_results_by_q <- vector("list", length(q_values))
artifacts_by_q <- vector("list", length(q_values))
combined_summary <- vector("list", length(q_values))
combined_subspace <- vector("list", length(q_values))
combined_subspace_summary <- vector("list", length(q_values))
artifact_rows <- vector("list", length(q_values))

for (index in seq_along(q_values)) {
  q_value <- q_values[[index]]
  cat(sprintf("\nRunning q=%d (%d/%d)\n", q_value, index, length(q_values)))

  simulation_results_by_q[[index]] <- run_rrr_uq_simulation(
    n_sims = n_sims,
    n = n,
    p = p,
    q = q_value,
    r = r,
    nnzeros = nnzeros,
    theta_strength = theta_strength,
    r_pca = r_pca,
    lambdas = lambdas,
    kfolds = kfolds,
    n_boot = n_boot,
    K = k_outer,
    strength_bins = strength_bins,
    seed = seed + index - 1L,
    preprocess_mode = "center",
    ccar3_path = ccar3_path,
    ccar3_code_path = ccar3_code_path,
    prefer_source = TRUE,
    parallelize = FALSE,
    verbose = TRUE
  )

  q_output_dir <- file.path(output_dir, sprintf("q_%03d", q_value))
  artifacts_by_q[[index]] <- write_rrr_uq_simulation_reports(
    simulation_results = simulation_results_by_q[[index]],
    output_dir = q_output_dir,
    verbose = TRUE
  )

  combined_summary[[index]] <- cbind(q = q_value, simulation_results_by_q[[index]]$summary$overall)
  combined_subspace[[index]] <- cbind(q = q_value, simulation_results_by_q[[index]]$subspace)
  combined_subspace_summary[[index]] <- cbind(q = q_value, simulation_results_by_q[[index]]$summary$subspace)
  artifact_rows[[index]] <- data.frame(
    q = q_value,
    output_dir = q_output_dir,
    summary_csv = artifacts_by_q[[index]]$summary_csv,
    subspace_csv = artifacts_by_q[[index]]$subspace_csv,
    subspace_summary_csv = artifacts_by_q[[index]]$subspace_summary_csv,
    stringsAsFactors = FALSE
  )
}

combined_summary_df <- do.call(rbind, combined_summary)
combined_subspace_df <- do.call(rbind, combined_subspace)
combined_subspace_summary_df <- do.call(rbind, combined_subspace_summary)
artifact_table <- do.call(rbind, artifact_rows)

summary_csv <- file.path(output_dir, "simulation_q_sweep_ci_summary.csv")
subspace_csv <- file.path(output_dir, "simulation_q_sweep_subspace_distance.csv")
subspace_summary_csv <- file.path(output_dir, "simulation_q_sweep_subspace_distance_summary.csv")
artifact_csv <- file.path(output_dir, "simulation_q_sweep_artifacts.csv")
subspace_pdf <- file.path(output_dir, "simulation_q_sweep_subspace_distance.pdf")

utils::write.csv(combined_summary_df, summary_csv, row.names = FALSE)
utils::write.csv(combined_subspace_df, subspace_csv, row.names = FALSE)
utils::write.csv(combined_subspace_summary_df, subspace_summary_csv, row.names = FALSE)
utils::write.csv(artifact_table, artifact_csv, row.names = FALSE)
plot_q_sweep_subspace_summary(combined_subspace_summary_df, subspace_pdf)

cat("\nQ sweep CI summary\n")
print(combined_summary_df)
cat("\nQ sweep subspace summary\n")
print(combined_subspace_summary_df)
cat("\nArtifacts\n")
cat(sprintf("Combined summary: %s\n", summary_csv))
cat(sprintf("Combined subspace distances: %s\n", subspace_csv))
cat(sprintf("Combined subspace summary: %s\n", subspace_summary_csv))
cat(sprintf("Per-q artifact index: %s\n", artifact_csv))
cat(sprintf("Subspace-distance comparison plot: %s\n", subspace_pdf))

if (interactive()) {
  assign(
    "rrr_uq_q_sweep_results",
    list(
      q_values = q_values,
      runs = simulation_results_by_q,
      artifacts = artifacts_by_q,
      combined_summary = combined_summary_df,
      combined_subspace = combined_subspace_df,
      combined_subspace_summary = combined_subspace_summary_df
    ),
    envir = .GlobalEnv
  )
}

invisible(
  list(
    q_values = q_values,
    runs = simulation_results_by_q,
    artifacts = artifacts_by_q,
    combined_summary = combined_summary_df,
    combined_subspace = combined_subspace_df,
    combined_subspace_summary = combined_subspace_summary_df,
    artifact_table = artifact_table
  )
)
