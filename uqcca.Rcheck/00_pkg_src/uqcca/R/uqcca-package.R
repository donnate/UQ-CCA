#' uqcca: Uncertainty Quantification for Sparse CCA
#'
#' `uqcca` provides a packaging layer around sparse canonical correlation
#' analysis workflows based on the external `ccar3` implementation. It includes
#' bootstrap summaries for sparse CCA loadings, cross-fitted regression-based
#' uncertainty quantification, and simulation helpers for benchmarking coverage,
#' false discovery, and subspace recovery.
#'
#' The core fitting functions expect access to either an installed `ccar3`
#' package or a local source checkout configured via
#' `options(uqcca.ccar3_path = "/path/to/ccar3")` or the
#' `UQCCA_CCAR3_PATH` environment variable. The simulation helpers additionally
#' use `options(uqcca.ccar3_code_path = "/path/to/CCAR3_code")` or
#' `UQCCA_CCAR3_CODE_PATH`.
#'
#' @docType package
#' @name uqcca
#' @keywords internal
#' @import grDevices
#' @import methods
#' @import stats
#' @import utils
"_PACKAGE"

utils::globalVariables(
  c(
    "bin_midpoint",
    "false_discovery_proportion",
    "false_positive_rate",
    "feature",
    "label",
    "method",
    "plot_estimate",
    "plot_lower",
    "plot_upper",
    "raw_plot",
    "series",
    "shrunken_plot",
    "side",
    "subspace_distance",
    "value"
  )
)
