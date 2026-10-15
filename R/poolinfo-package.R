#' @details
#' The diagnostics describe a set of site reliabilities. A site's
#' reliability is the proportion of the variance of its estimate that comes
#' from the variation of the true effects; under the normal hierarchical
#' model it is the weight of the site's own estimate, against the mean of
#' the effects, in the site's posterior mean.
#'
#' The main functions:
#' * [reliability_posterior()] computes the posterior of the diagnostics
#'   from estimates and standard errors, by numerical integration over the
#'   between-site standard deviation.
#' * [reliability_profile()] computes them at a given between-site standard
#'   deviation, or at the estimate of a metafor fit.
#' * [reliability_draws()] computes them from posterior draws of another
#'   fit, such as one of brms.
#' * [reliability_design()], [balanced_allocation()] and [solve_design()]
#'   plan a trial.
#' * [compare_profiles()], [report_table()] and [report_text()] compare and
#'   report analyses.
#' * `plot()` and the functions `plot_*()` show them.
#'
#' `vignette("poolinfo")` works through an analysis, and the other
#' vignettes cover the profiles, the posterior summaries, a multisite trial,
#' meta-analysis, small area estimation, planning and the computation.
#'
#' @section Funding:
#' This research was supported by the Institute of Education Sciences,
#' U.S. Department of Education, through Grant R305D240078 to the
#' University of Alabama. The opinions expressed are those of the
#' author and do not represent views of the Institute or the U.S.
#' Department of Education.
#'
#' @keywords internal
#' @importFrom rlang .data
#' @importFrom stats update
#' @importFrom graphics plot
"_PACKAGE"
