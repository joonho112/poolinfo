#' Prior distribution of the between-site standard deviation
#'
#' Specifies the prior of `sigma_tau` for [reliability_posterior()]: the flat
#' prior of Lee (2026), which is the default, or a half-normal, half-t,
#' half-Cauchy or uniform prior.
#'
#' @param type Family of the prior: one of `"flat"` (the default),
#'   `"half_normal"`, `"half_t"`, `"half_cauchy"` or `"uniform"`.
#' @param scale Scale of a half-normal, half-t or half-Cauchy prior, a
#'   positive number in the units of the estimates. It has no default,
#'   because a sensible scale depends on the scale of the effects.
#' @param df Degrees of freedom of a half-t prior, a positive number.
#' @param upper Upper limit of a uniform prior, a positive number.
#'
#' @details
#' Proper SD priors are an experimental extension beyond Lee (2026); their
#' interfaces may change in a later 0.x release.
#'
#' The flat prior is uniform on the mean of the effects and on `sigma_tau`
#' (not on `sigma_tau^2`). Under it the posterior is proper only with at least
#' 3 sites, and the posterior mean of \eqn{\sigma_\tau^{q}}{sigma_tau^q}
#' exists only with more than \eqn{q + 2} sites; so with 3 or 4 sites the
#' posterior mean of `sigma_tau^2` does not exist, although those of the
#' bounded diagnostics do. With a mean model of \eqn{p} columns, the
#' condition is \eqn{J > q + p + 1}{J > q + p + 1}.
#'
#' A prior whose density falls off like \eqn{\sigma^{-k}}{sigma^-k} lowers
#' the bound by \eqn{k}: the moment of order \eqn{q} exists when
#' \eqn{J > q + p + 1 - k}{J > q + p + 1 - k}, with \eqn{k = \nu + 1} for a
#' half-t prior with \eqn{\nu}{nu} degrees of freedom and \eqn{k = 2} for
#' the half-Cauchy. Every moment exists under the half-normal and the
#' uniform priors. Every prior other than the flat one needs at least
#' \eqn{p + 1} sites.
#'
#' With few sites the flat prior can matter, as for any prior on a variance
#' parameter (Gelman, 2006); [prior_sensitivity()] compares the summaries
#' under several priors.
#'
#' @return An object of class `sd_prior`, a list with the components
#'   * `type`: the family.
#'   * `scale`, `df` or `upper`: the parameters that the family uses.
#'   * `tail`: the power \eqn{k} with which the density falls off in its
#'     upper tail, which decides which posterior moments exist: 0 for the
#'     flat prior, `df + 1` for the half-t, 2 for the half-Cauchy and `Inf`
#'     for the half-normal and uniform priors.
#'
#' @seealso [prior_sensitivity()] to compare the summaries under several
#'   priors, and `vignette("posterior-summaries")` for the choice of prior.
#'
#' @references
#' Gelman, A. (2006). Prior distributions for variance parameters in
#' hierarchical models (comment on article by Browne and Draper). *Bayesian
#' Analysis, 1*(3), 515-534. \doi{10.1214/06-BA117A}
#'
#' @template ref-lee2026
#' @examples
#' sd_prior()
#' sd_prior("half_normal", scale = 0.5)
#' sd_prior("half_t", scale = 0.5, df = 3)
#' @export
sd_prior <- function(type = c("flat", "half_normal", "half_t", "half_cauchy", "uniform"),
                     scale = NULL, df = NULL, upper = NULL) {
  type <- .rp_match(type, c("flat", "half_normal", "half_t", "half_cauchy", "uniform"), "type")
  positive <- function(x, arg) {
    .rp_check_number(x, arg, "a single positive number", lower = 0, lower_open = TRUE)
  }
  unused <- function(args) {
    given <- args[!vapply(list(scale = scale, df = df, upper = upper)[args], is.null, TRUE)]
    if (length(given) > 0L) {
      .rp_abort(
        sprintf("A %s prior has no %s.", gsub("_", "-", type),
                paste(sprintf("`%s`", given), collapse = " or ")),
        info = "Leave out the arguments that the family does not have."
      )
    }
  }
  need <- function(x, arg) {
    if (is.null(x)) {
      .rp_abort(
        sprintf("A %s prior needs `%s`.", gsub("_", "-", type), arg),
        info = if (arg == "scale") "The scale is in the units of the estimates and has no default."
      )
    }
    positive(x, arg)
  }
  out <- switch(type,
    flat = {
      unused(c("scale", "df", "upper"))
      list(type = "flat", tail = 0)
    },
    half_normal = {
      unused(c("df", "upper"))
      list(type = type, scale = need(scale, "scale"), tail = Inf)
    },
    half_t = {
      unused("upper")
      df <- need(df, "df")
      list(type = type, scale = need(scale, "scale"), df = df, tail = df + 1)
    },
    half_cauchy = {
      unused(c("df", "upper"))
      list(type = type, scale = need(scale, "scale"), df = 1, tail = 2)
    },
    uniform = {
      unused(c("scale", "df"))
      list(type = type, upper = need(upper, "upper"), tail = Inf)
    }
  )
  structure(out, class = "sd_prior")
}

#' @export
format.sd_prior <- function(x, ...) {
  switch(x$type,
    flat = "flat on the mean and on sigma_tau",
    half_normal = sprintf("half-normal on sigma_tau with scale %s, flat on the mean",
                          .rp_fmt_sig(x$scale)),
    half_t = sprintf("half-t on sigma_tau with scale %s and %s degrees of freedom, flat on the mean",
                     .rp_fmt_sig(x$scale), .rp_fmt_sig(x$df)),
    half_cauchy = sprintf("half-Cauchy on sigma_tau with scale %s, flat on the mean",
                          .rp_fmt_sig(x$scale)),
    uniform = sprintf("uniform on sigma_tau from 0 to %s, flat on the mean", .rp_fmt_sig(x$upper))
  )
}

#' @export
print.sd_prior <- function(x, ...) {
  cat("Prior: ", format(x), "\n", sep = "")
  invisible(x)
}

.rp_check_prior <- function(prior) {
  if (!inherits(prior, "sd_prior")) {
    .rp_abort("`prior` must be made by `sd_prior()`.", info = .rp_you_passed(prior))
  }
  prior
}

# Log density of the prior of sigma, up to a constant, with the value 0 at
# sigma = 0 for the families with a mode there.
.rp_prior_log_density <- function(prior, sigma) {
  switch(prior$type,
    flat = rep(0, length(sigma)),
    half_normal = -0.5 * (sigma / prior$scale)^2,
    half_t = ,
    half_cauchy = -0.5 * (prior$df + 1) * log1p((sigma / prior$scale)^2 / prior$df),
    uniform = ifelse(sigma <= prior$upper, 0, -Inf)
  )
}

# Whether the posterior moment of sigma^q exists, with J sites and p columns
# of the mean model: J > q + p + 1 - k, where sigma^-k is the tail of the
# prior (k = Inf for the families with light tails).
.rp_moment_exists <- function(q, J, p, prior) {
  J > q + p + 1 - prior$tail
}

# The smallest number of sites for which E(sigma^q) exists: J > q + p + 1 -
# tail. Under the half-normal and uniform priors every moment exists.
.rp_moment_min_sites <- function(q, p, prior) {
  as.integer(floor(q + p + 1 - prior$tail) + 1)
}

# The prior as a phrase; with a mean model, the flat prior is on its
# coefficients.
.rp_prior_text <- function(prior, X = NULL) {
  out <- format(prior)
  if (!is.null(X)) out <- sub("the mean", "the coefficients of the mean model", out, fixed = TRUE)
  out
}

# The smallest number of sites for a proper posterior.
.rp_min_sites <- function(p, prior) {
  if (prior$type == "flat") p + 2L else p + 1L
}
