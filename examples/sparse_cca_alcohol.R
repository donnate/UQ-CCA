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

parse_variance_shrinkage_methods <- function(value) {
  methods <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  methods <- methods[nzchar(methods)]
  if (!length(methods)) {
    methods <- "james-stein"
  }

  unique(vapply(methods, match_variance_shrinkage, character(1)))
}

ccar3_path <- Sys.getenv("UQCCA_CCAR3_PATH", unset = "/Users/clairedonnat/Documents/ccar3")
n_boot <- as.integer(Sys.getenv("UQCCA_N_BOOT", unset = "25"))
k_outer <- as.integer(Sys.getenv("UQCCA_K", unset = "5"))
features_per_page <- as.integer(Sys.getenv("UQCCA_FEATURES_PER_PAGE", unset = "40"))
variance_shrinkage_methods <- parse_variance_shrinkage_methods(
  Sys.getenv("UQCCA_VARIANCE_SHRINKAGE", unset = "james-stein,none,mr.mashr")
)
output_dir <- Sys.getenv(
  "UQCCA_OUTPUT_DIR",
  unset = file.path(project_root, "results", "alcohol_uq")
)

api <- get_ccar3_api(
  ccar3_path = ccar3_path,
  prefer_source = TRUE,
  quiet = FALSE
)
data <- load_alcohol_example()
fit <- fit_ecca_cv(
  X = data$X,
  Y = data$Y,
  r = 2,
  lambdas = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2),
  kfolds = 5,
  preprocess_mode = "center",
  ccar3_api = api,
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
  ccar3_api = api,
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
  variant_label <- variance_shrinkage_label(method)
  variant_dir <- file.path(
    output_dir,
    paste0("crossfit_", sanitize_path_label(method))
  )

  message(sprintf("Running cross-fit inference with %s.", variant_label))
  crossfit_result <- tryCatch(
    crossfit_cca_inference(
      X = data$X,
      Y = data$Y,
      reference_fit = fit,
      K = k_outer,
      seed = 2026,
      alpha = 0.05,
      fit_mode = "fixed_lambda",
      ccar3_api = api,
      progress = TRUE,
      variance_shrinkage = method
    ),
    error = function(e) e
  )

  if (inherits(crossfit_result, "error")) {
    crossfit_failures[[method]] <- conditionMessage(crossfit_result)
    warning(
      sprintf("Skipping %s: %s", variant_label, conditionMessage(crossfit_result)),
      call. = FALSE
    )
    next
  }

  variant_results <- list(
    data = data,
    fit = fit,
    bootstrap = bootstrap,
    crossfit = crossfit_result
  )
  variant_results$comparison <- compare_uq_methods(
    reference_fit = fit,
    bootstrap_result = bootstrap,
    crossfit_result = crossfit_result,
    component = 1,
    top_n = 10,
    alpha = 0.05
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
    label = variance_shrinkage_label(method),
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
cat("\nCross-fit shrinkage comparison\n")
print(variant_summary)

cat("\nFull coefficient CI reports\n")
cat(sprintf("Gene coefficient CIs: %s\n", artifacts$gene_csv))
cat(sprintf("Methylation coefficient CIs: %s\n", artifacts$methylation_csv))
cat(sprintf("Gene CI plots: %s\n", artifacts$gene_pdf))
cat(sprintf("Methylation CI plots: %s\n", artifacts$methylation_pdf))
cat(sprintf("Cross-fit comparison summary: %s\n", variant_summary_csv))

for (method in successful_methods) {
  variant_label <- variance_shrinkage_label(method)
  variant_artifacts <- crossfit_artifacts[[method]]
  cat(sprintf("\n%s\n", variant_label))
  cat(sprintf("Cross-fit gene regression CIs: %s\n", variant_artifacts$gene_csv))
  cat(sprintf("Cross-fit methylation regression CIs: %s\n", variant_artifacts$methylation_csv))
  cat(sprintf("Cross-fit gene CI plots: %s\n", variant_artifacts$gene_pdf))
  cat(sprintf("Cross-fit methylation CI plots: %s\n", variant_artifacts$methylation_pdf))
  cat(sprintf("Cross-fit variance diagnostics: %s\n", variant_artifacts$variance_csv))
  cat(sprintf("Cross-fit variance diagnostic plot: %s\n", variant_artifacts$variance_pdf))
  cat(sprintf("Cross-fit alignment diagnostics: %s\n", variant_artifacts$alignment_csv))
  cat(sprintf("Cross-fit preprocess diagnostics: %s\n", variant_artifacts$preprocess_csv))
}

if (length(crossfit_failures)) {
  cat("\nSkipped cross-fit variants\n")
  for (method in names(crossfit_failures)) {
    cat(sprintf(
      "%s: %s\n",
      variance_shrinkage_label(method),
      crossfit_failures[[method]]
    ))
  }
}

if (interactive()) {
  assign("alcohol_uq_results", results, envir = .GlobalEnv)
}

invisible(results)
