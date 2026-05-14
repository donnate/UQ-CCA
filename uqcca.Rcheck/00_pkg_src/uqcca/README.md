# uqcca

`uqcca` packages the sparse CCA uncertainty-quantification workflows in this
repository so they can be shared with collaborators as a conventional R
package.

It includes:

- sparse CCA fitting through `ccar3::ecca` and `ccar3::ecca.cv`
- bootstrap uncertainty summaries for loadings, canonical correlations, and
  subspace stability
- cross-fitted regression-based UQ for conditional and marginal
  feature-score associations
- simulation helpers for benchmarking coverage, power, false discovery, and
  subspace recovery

## Install

From a local checkout:

```r
devtools::install("/path/to/UQ-CCA")
```

## Configure External Dependencies

The fitting functions expect access to `ccar3`. You can either install it as an
R package or point `uqcca` to a local source checkout:

```r
options(uqcca.ccar3_path = "/path/to/ccar3")
```

or

```r
Sys.setenv(UQCCA_CCAR3_PATH = "/path/to/ccar3")
```

The simulation helpers additionally use the external `CCAR3_code` repository:

```r
options(uqcca.ccar3_code_path = "/path/to/CCAR3_code")
Sys.setenv(UQCCA_CCAR3_CODE_PATH = "/path/to/CCAR3_code")
```

## Quick Start

```r
library(uqcca)

alcohol <- load_alcohol_example()

fit <- fit_ecca_cv(
  X = alcohol$X,
  Y = alcohol$Y,
  r = 2,
  lambdas = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2),
  kfolds = 5,
  preprocess_mode = "center"
)

boot <- bootstrap_cca_uq(
  X = alcohol$X,
  Y = alcohol$Y,
  reference_fit = fit,
  B = 25,
  refit_mode = "fixed_lambda"
)

marginal <- crossfit_cca_marginal_inference(
  X = alcohol$X,
  Y = alcohol$Y,
  reference_fit = fit,
  K = 5
)

conditional <- crossfit_cca_conditional_inference(
  X = alcohol$X,
  Y = alcohol$Y,
  reference_fit = fit,
  K = 5
)
```

## High-Level Workflows

The package includes end-to-end wrappers for the two main use cases:

```r
results <- run_alcohol_uq_comparison()
sim <- run_rrr_uq_simulation()
```

To write collaborator-friendly CSV and PDF outputs:

```r
write_bootstrap_loading_reports(results, "results/alcohol_bootstrap")
write_crossfit_loading_reports(results, "results/alcohol_crossfit")
write_rrr_uq_simulation_reports(sim, "results/simulation_rrr_uq")
```

## Installed Scripts

Standalone scripts are installed with the package under:

```r
system.file("scripts", package = "uqcca")
```

They mirror the repository examples:

- `sparse_cca_alcohol.R`
- `simulate_rrr_uq_power.R`
- `simulate_rrr_uq_by_q.R`

## Vignette

After installation, open the vignette with:

```r
browseVignettes("uqcca")
```
