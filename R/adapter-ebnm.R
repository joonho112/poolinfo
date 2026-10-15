#' Profiles of ebnm fits with a normal prior
#'
#' `reliability_profile()` of a fit of the ebnm package with a normal prior
#' gives the profile at the standard deviation of the fitted prior.
#'
#' @param x An `ebnm` object fitted with `prior_family = "normal"`.
#' @param ... Must be empty.
#'
#' @details
#' The method reads the elements `data$x`, `data$s` and `fitted_g` of the
#' fit, where `fitted_g` must be a normal distribution (a `normalmix` with
#' one component). ebnm's normal prior has its mean fixed at 0 unless the fit
#' used `mode = "estimate"`; the reliabilities depend only on the standard
#' deviation, but the site table's posterior means use the precision-weighted
#' mean and so differ from ebnm's own when the mode was fixed. A prior of
#' another family is not a normal model and is refused.
#'
#' @return A `reliability_profile` object whose source is the ebnm fit.
#'
#' @examplesIf requireNamespace("ebnm", quietly = TRUE)
#' fit <- ebnm::ebnm(eight_schools$estimate, eight_schools$se, prior_family = "normal",
#'                   mode = "estimate")
#' reliability_profile(fit)
#' @name ebnm-methods
NULL

#' @rdname ebnm-methods
#' @export
reliability_profile.ebnm <- function(x, ...) {
  .rp_check_dots_empty(...)
  g <- x$fitted_g
  normal <- inherits(g, "normalmix") && length(g$sd) == 1L && isTRUE(all.equal(sum(g$pi), 1))
  if (!normal) {
    .rp_abort(
      "Only an ebnm fit with a normal prior has a between-site standard deviation.",
      info = "The fitted prior is not a single normal distribution.",
      fix = "Fit with `prior_family = \"normal\"`.",
      class = "rp_scope_error"
    )
  }
  if (is.null(x$data$x) || is.null(x$data$s)) {
    .rp_abort("This ebnm fit has no data; fit it with the default `output`.")
  }
  .rp_new_profile(
    se = as.numeric(x$data$s), se2 = NULL, sigma_tau = as.numeric(g$sd), site = NULL,
    estimate = as.numeric(x$data$x), source = "ebnm normal prior", call = match.call()
  )
}
