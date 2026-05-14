suppressPackageStartupMessages(library(uqcca))

env_or_option <- function(env_name, option_name, default = NULL) {
  value <- Sys.getenv(env_name, unset = "")
  if (!nzchar(value)) {
    value <- getOption(option_name, default = default)
  }

  if (is.null(value) || !length(value) || !nzchar(as.character(value)[1])) {
    return(default)
  }

  as.character(value)[1]
}

parse_variance_shrinkage_methods <- function(value) {
  methods <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  methods <- methods[nzchar(methods)]
  if (!length(methods)) {
    methods <- "james-stein"
  }

  unique(methods)
}

ccar3_path <- env_or_option("UQCCA_CCAR3_PATH", "uqcca.ccar3_path")
n_boot <- as.integer(Sys.getenv("UQCCA_N_BOOT", unset = "25"))
k_outer <- as.integer(Sys.getenv("UQCCA_K", unset = "5"))
features_per_page <- as.integer(Sys.getenv("UQCCA_FEATURES_PER_PAGE", unset = "40"))
variance_shrinkage_methods <- parse_variance_shrinkage_methods(
  Sys.getenv("UQCCA_VARIANCE_SHRINKAGE", unset = "james-stein,none,mr.mashr")
)
output_dir <- Sys.getenv(
  "UQCCA_OUTPUT_DIR",
  unset = file.path(getwd(), "results", "alcohol_uq")
)

data <- load_alcohol_example()
fit <- fit_ecca_cv(
  X = data$X,
  Y = data$Y,
  r = 2,
  lambdas = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2),
  kfolds = 5,
  preprocess_mode = "center",
  ccar3_path = ccar3_path,
  prefer_source = TRUE,
  parallelize = FALSE,
  verbose = TRUE
)
bootstrap <- bootstrap_cca_uq(
  X = data$X,
  Y = data$Y,
  reference_fit = fit,
  B = n_boot,
  seed = 2025,
  refit_mode = "fixed_lambda",
  progress = TRUE
)

bootstrap_results <- list(
  data = data,
  fit = fit,
  bootstrap = bootstrap
)

artifacts <- write_bootstrap_loading_reports(
  results = bootstrap_results,
  output_dir = output_dir,
  alpha = 0.05,
  features_per_page = features_per_page,
  verbose = TRUE
)
crossfit_variants <- list()
comparison_variants <- list()
crossfit_artifacts <- list()
crossfit_failures <- list()

for (method in variance_shrinkage_methods) {
  variant_dir <- file.path(output_dir, paste0("crossfit_", gsub("[^A-Za-z0-9]+", "_", method)))

  message(sprintf("Running cross-fit inference with %s.", method))
  crossfit_result <- tryCatch(
    crossfit_cca_inference(
      X = data$X,
      Y = data$Y,
      reference_fit = fit,
      K = k_outer,
      seed = 2026,
      alpha = 0.05,
      fit_mode = "fixed_lambda",
      progress = TRUE,
      variance_shrinkage = method
    ),
    error = function(e) e
  )

  if (inherits(crossfit_result, "error")) {
    crossfit_failures[[method]] <- conditionMessage(crossfit_result)
    warning(sprintf("Skipping %s: %s", method, conditionMessage(crossfit_result)), call. = FALSE)
    next
  }

  variant_results <- list(
    data = data,
    fit = fit,
    bootstrap = bootstrap,
    crossfit = crossfit_result,
    comparison = compare_uq_methods(
      reference_fit = fit,
      bootstrap_result = bootstrap,
      crossfit_result = crossfit_result,
      component = 1,
      top_n = 10,
      alpha = 0.05
    )
  )

  crossfit_variants[[method]] <- crossfit_result
  comparison_variants[[method]] <- variant_results$comparison
  crossfit_artifacts[[method]] <- write_crossfit_loading_reports(
    results = variant_results,
    output_dir = variant_dir,
    features_per_page = features_per_page,
    verbose = TRUE
  )
}

successful_methods <- names(crossfit_variants)
if (!length(successful_methods)) {
  stop("Cross-fit inference failed for every requested variance shrinkage method.", call. = FALSE)
}

default_method <- successful_methods[[1]]
variant_summary <- do.call(rbind, lapply(successful_methods, function(method) {
  crossfit_result <- crossfit_variants[[method]]
  data.frame(
    variance_shrinkage = method,
    significant_gene = sum(crossfit_result$x_results$significant, na.rm = TRUE),
    significant_methylation = sum(crossfit_result$y_results$significant, na.rm = TRUE),
    mean_fold_correlation = mean(crossfit_result$fold_correlations, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
variant_summary_csv <- file.path(output_dir, "crossfit_variance_shrinkage_summary.csv")
utils::write.csv(variant_summary, variant_summary_csv, row.names = FALSE)

results <- list(
  data = data,
  fit = fit,
  bootstrap = bootstrap,
  crossfit = crossfit_variants[[default_method]],
  comparison = comparison_variants[[default_method]],
  crossfit_variants = crossfit_variants,
  comparison_variants = comparison_variants,
  crossfit_failures = crossfit_failures,
  default_crossfit_method = default_method,
  variance_shrinkage_methods = variance_shrinkage_methods,
  artifacts = list(
    bootstrap = artifacts,
    crossfit = crossfit_artifacts,
    crossfit_summary_csv = variant_summary_csv
  )
)

print_alcohol_uq_summary(results, top_n = 10)
print(variant_summary)

invisible(results)
