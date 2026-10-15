#' Profiles and posteriors of simulated multisite data
#'
#' The data simulated by `multisiteDGP::sim_multisite()` and
#' `multisiteDGP::sim_meta()` carry the true between-site standard deviation
#' of their design. `reliability_profile()` gives the profile at that true
#' value, which an analysis of real data could not know.
#' `reliability_posterior()` fits the posterior from the estimates and the
#' standard errors alone, as an analysis of real data would.
#'
#' @param x A `multisitedgp_data` object.
#' @param ... Must be empty.
#'
#' @details
#' The method reads the columns `site_index`, `tau_j_hat` and `se2_j` and the
#' element `sigma_tau` of the attribute `design`, which are documented in
#' `?multisiteDGP::multisitedgp_data`. The true effects `tau_j` are not used.
#' Selecting rows retains the simulation design. If the design attribute
#' has been removed or has no `sigma_tau`, supply the between-site standard
#' deviation directly with `sigma_tau`, together with `se2 = x$se2_j`.
#'
#' @return A `reliability_profile` object whose source is the design of the
#'   simulation.
#'
#' @examplesIf requireNamespace("multisiteDGP", quietly = TRUE)
#' design <- multisiteDGP::multisitedgp_design(J = 30, sigma_tau = 0.2, seed = 1)
#' d <- multisiteDGP::sim_multisite(design)
#' reliability_profile(d)
#' reliability_posterior(d)
#' @name multisitedgp-methods
NULL

#' @rdname multisitedgp-methods
#' @export
reliability_profile.multisitedgp_data <- function(x, ...) {
  .rp_check_dots_empty(...)
  cols <- .rp_multisitedgp_columns(x)
  design <- attr(x, "design")
  sigma_tau <- if (is.list(design)) design$sigma_tau else NULL
  if (is.null(sigma_tau)) {
    .rp_abort(
      "These simulated data have no known design value of `sigma_tau`.",
      info = "The attribute `design` is absent or has no `sigma_tau`.",
      fix = "Give the standard deviation directly with `sigma_tau`, together with `se2 = x$se2_j`."
    )
  }
  .rp_new_profile(
    se = NULL, se2 = cols$se2, sigma_tau = sigma_tau, site = cols$site, estimate = cols$estimate,
    source = "multisiteDGP design value", call = match.call()
  )
}

.rp_multisitedgp_columns <- function(x) {
  need <- c("site_index", "tau_j_hat", "se2_j")
  missing_cols <- setdiff(need, names(x))
  if (length(missing_cols) > 0L) {
    .rp_abort(
      "These simulated data lack the columns that the method reads.",
      info = sprintf("Missing: %s.", paste(sprintf("`%s`", missing_cols), collapse = ", "))
    )
  }
  list(site = as.character(x[["site_index"]]), estimate = as.numeric(x[["tau_j_hat"]]),
       se2 = as.numeric(x[["se2_j"]]))
}
