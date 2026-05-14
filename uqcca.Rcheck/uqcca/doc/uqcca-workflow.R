## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = FALSE)

## -----------------------------------------------------------------------------
# devtools::install("/path/to/UQ-CCA")

## -----------------------------------------------------------------------------
# options(uqcca.ccar3_path = "/path/to/ccar3")

## -----------------------------------------------------------------------------
# Sys.setenv(UQCCA_CCAR3_PATH = "/path/to/ccar3")

## -----------------------------------------------------------------------------
# options(uqcca.ccar3_code_path = "/path/to/CCAR3_code")

## -----------------------------------------------------------------------------
# library(uqcca)
# 
# alcohol <- load_alcohol_example()
# str(alcohol, max.level = 1)

## -----------------------------------------------------------------------------
# fit <- fit_ecca_cv(
#   X = alcohol$X,
#   Y = alcohol$Y,
#   r = 2,
#   lambdas = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2),
#   kfolds = 5,
#   preprocess_mode = "center"
# )

## -----------------------------------------------------------------------------
# boot <- bootstrap_cca_uq(
#   X = alcohol$X,
#   Y = alcohol$Y,
#   reference_fit = fit,
#   B = 50,
#   seed = 2025,
#   refit_mode = "fixed_lambda"
# )
# 
# gene_boot <- build_bootstrap_loading_table_all(boot, side = "x")
# meth_boot <- build_bootstrap_loading_table_all(boot, side = "y")
# rho_boot <- build_bootstrap_correlation_table(boot)
# subspace_boot <- summarize_bootstrap_subspace(boot, components = 1:2)

## -----------------------------------------------------------------------------
# write_bootstrap_loading_reports(
#   results = list(bootstrap = boot),
#   output_dir = "results/alcohol_bootstrap"
# )

## -----------------------------------------------------------------------------
# crossfit <- crossfit_cca_inference(
#   X = alcohol$X,
#   Y = alcohol$Y,
#   reference_fit = fit,
#   K = 5,
#   seed = 2026,
#   variance_shrinkage = "james-stein"
# )

## -----------------------------------------------------------------------------
# marginal <- crossfit_cca_marginal_inference(
#   X = alcohol$X,
#   Y = alcohol$Y,
#   reference_fit = fit,
#   K = 5,
#   seed = 2026
# )
# 
# conditional <- crossfit_cca_conditional_inference(
#   X = alcohol$X,
#   Y = alcohol$Y,
#   reference_fit = fit,
#   K = 5,
#   seed = 2026,
#   method = "auto"
# )

## -----------------------------------------------------------------------------
# results <- run_alcohol_uq_comparison(
#   n_boot = 25,
#   K = 5,
#   preprocess_mode = "center"
# )
# 
# print_alcohol_uq_summary(results)

## -----------------------------------------------------------------------------
# sim <- run_rrr_uq_simulation(
#   n_sims = 20,
#   n = 500,
#   p = 100,
#   q = 20,
#   r = 2,
#   nnzeros = 10,
#   n_boot = 25,
#   K = 2,
#   preprocess_mode = "center"
# )
# 
# sim$summary$overall
# sim$summary$subspace

## -----------------------------------------------------------------------------
# write_rrr_uq_simulation_reports(
#   simulation_results = sim,
#   output_dir = "results/simulation_rrr_uq"
# )

## -----------------------------------------------------------------------------
# system.file("scripts", package = "uqcca")

