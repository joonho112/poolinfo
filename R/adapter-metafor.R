#' Profiles and posteriors of metafor fits
#'
#' `reliability_profile()` of a random-effects fit of the metafor package
#' gives the profile at its estimate of \eqn{\tau^2}{tau^2}, a plug-in
#' profile; `reliability_posterior()` of the same fit integrates over the
#' between-study standard deviation instead, with the same estimates,
#' sampling variances and moderators.
#'
#' @param x An `rma.uni` object, from `metafor::rma()` or
#'   `metafor::rma.uni()`.
#' @param ... Must be empty.
#'
#' @details
#' The method reads the elements `yi`, `vi`, `tau2`, `method` and `slab`
#' of the fit. With moderators, \eqn{\tau^2}{tau^2} is the residual
#' heterogeneity, the variance of the true effects around the
#' meta-regression, and the reliabilities are computed from it. In metafor's
#' notation
#' \eqn{\tau^2}{tau^2} is what this package writes
#' \eqn{\sigma_\tau^2}{sigma_tau^2}. A fit whose \eqn{\tau^2}{tau^2} is
#' estimated as 0, as often happens with few studies, gives the limit of the
#' profile at `sigma_tau = 0` and a message; the posterior has no such
#' boundary.
#'
#' [map_tau2_interval()] maps the confidence interval of the fit through the
#' diagnostics. Other fits of metafor, such as those of `metafor::rma.mv()`
#' with several variance components, are refused with an error: the
#' diagnostics need one between-study variance and known sampling variances.
#' Fixed-, equal- and common-effect fits (`method = "FE"`, `"EE"` or
#' `"CE"`) are also refused. Refit a random-effects model, or supply the
#' estimates and sampling variances directly for a new analysis.
#'
#' @return A `reliability_profile` object whose source names the method of
#'   the fit.
#'
#' @seealso [map_tau2_interval()], and `vignette("meta-analysis")` for an
#'   example.
#'
#' @examplesIf requireNamespace("metafor", quietly = TRUE)
#' fit <- metafor::rma(yi = estimate, sei = se, slab = school, data = eight_schools)
#' reliability_profile(fit)
#' reliability_posterior(fit)
#' @name metafor-methods
NULL

#' @rdname metafor-methods
#' @export
reliability_profile.rma.uni <- function(x, ...) {
  .rp_check_dots_empty(...)
  cols <- .rp_rma_columns(x)
  .rp_new_profile(
    se = NULL, se2 = cols$se2, sigma_tau = sqrt(max(0, x$tau2)), site = cols$site,
    estimate = cols$estimate, source = sprintf("metafor %s estimate", x$method),
    X = cols$X, call = match.call()
  )
}

#' @rdname metafor-methods
#' @export
reliability_profile.rma <- function(x, ...) {
  .rp_abort_rma_scope(x)
}

# Fits of metafor other than rma.uni.
.rp_abort_rma_scope <- function(x) {
  .rp_abort(
    sprintf("poolinfo reads only `rma.uni` fits of metafor, not %s.", .rp_describe(x)),
    info = "The diagnostics need one between-study variance and known sampling variances, as in `metafor::rma()`.",
    fix = "Fit the model with `metafor::rma()`, or give the estimates and their sampling variances (`se2 = vi`).",
    class = "rp_scope_error"
  )
}

.rp_rma_columns <- function(x) {
  .rp_check_rma_scope(x)
  if (is.null(x$yi) || is.null(x$vi)) {
    .rp_abort("This metafor fit has no estimates or sampling variances.")
  }
  site <- if (!is.null(x$slab)) as.character(x$slab) else NULL
  if (!is.null(site) && anyDuplicated(site) > 0L) site <- NULL
  X <- x$X
  if (!is.null(X) && ncol(X) == 1L && all(X[, 1L] == 1)) X <- NULL
  list(estimate = as.numeric(x$yi), se2 = as.numeric(x$vi), site = site, X = X)
}

.rp_check_rma_scope <- function(x) {
  if (x$method %in% c("FE", "EE", "CE")) {
    .rp_abort(
      sprintf("The metafor method \"%s\" is not a random-effects model.", x$method),
      info = "Profiles and posteriors of metafor fits require a between-site variance component.",
      fix = paste("Refit with a random-effects method, such as `method = \"REML\"`,",
                  "or give the estimates and sampling variances directly for a new analysis."),
      class = "rp_scope_error"
    )
  }
  invisible(x)
}
