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

project_root <- find_project_root()
source(file.path(project_root, "R", "cca_uq_methods.R"))
source(file.path(project_root, "R", "cca_uq_simulation.R"))

ccar3_path <- Sys.getenv("UQCCA_CCAR3_PATH", unset = "/Users/clairedonnat/Documents/ccar3")
ccar3_code_path <- Sys.getenv("UQCCA_CCAR3_CODE_PATH", unset = "/Users/clairedonnat/Documents/CCAR3_code")
output_dir <- Sys.getenv(
  "UQCCA_SIM_OUTPUT_DIR",
  unset = file.path(project_root, "results", "simulation_rrr_uq")
)

n_sims <- as.integer(Sys.getenv("UQCCA_SIM_N_REPS", unset = "10"))
n <- as.integer(Sys.getenv("UQCCA_SIM_N", unset = "100"))
p <- as.integer(Sys.getenv("UQCCA_SIM_P", unset = "100"))
q <- as.integer(Sys.getenv("UQCCA_SIM_Q", unset = "20"))
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

simulation_results <- run_rrr_uq_simulation(
  n_sims = n_sims,
  n = n,
  p = p,
  q = q,
  r = r,
  nnzeros = nnzeros,
  theta_strength = theta_strength,
  r_pca = r_pca,
  lambdas = lambdas,
  kfolds = kfolds,
  n_boot = n_boot,
  K = k_outer,
  strength_bins = strength_bins,
  seed = seed,
  preprocess_mode = "center",
  ccar3_path = ccar3_path,
  ccar3_code_path = ccar3_code_path,
  prefer_source = TRUE,
  parallelize = FALSE,
  verbose = TRUE
)

artifacts <- write_rrr_uq_simulation_reports(
  simulation_results = simulation_results,
  output_dir = output_dir,
  verbose = TRUE
)

cat("Simulation summary\n")
print(simulation_results$summary$overall)
cat("\nSubspace recovery summary\n")
print(simulation_results$summary$subspace)
cat("\nArtifacts\n")
cat(sprintf("Outcomes: %s\n", artifacts$outcomes_csv))
cat(sprintf("Summary: %s\n", artifacts$summary_csv))
cat(sprintf("By component: %s\n", artifacts$component_csv))
cat(sprintf("By replication: %s\n", artifacts$replication_csv))
cat(sprintf("Power by strength: %s\n", artifacts$strength_csv))
cat(sprintf("Subspace distances: %s\n", artifacts$subspace_csv))
cat(sprintf("Subspace summary: %s\n", artifacts$subspace_summary_csv))
cat(sprintf("Summary plot: %s\n", artifacts$summary_pdf))
cat(sprintf("Power-by-strength plot: %s\n", artifacts$strength_pdf))
cat(sprintf("Subspace-distance plot: %s\n", artifacts$subspace_pdf))

if (interactive()) {
  assign("rrr_uq_simulation_results", simulation_results, envir = .GlobalEnv)
}

invisible(simulation_results)
