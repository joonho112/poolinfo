# Site reliabilities and the logit scale on which the package computes.
#
# The reliabilities are handled through their logits. With sampling
# variances v_j and a between-site SD sigma, logit S_j = log(sigma^2) - log(v_j)
# exactly, so the logits are formed without computing S_j first, which would
# lose the digits of S_j near 0 or 1. The logits are split into their mean m,
# the logit of the informativeness, and the centered deviations d_j, which
# depend on the v_j only.

#' Site reliabilities
#'
#' Computes the reliability of each site: the proportion of the variance of
#' its estimate that comes from the variation of the true effects across
#' sites. Under the normal hierarchical model with a known mean and a known
#' between-site standard deviation, it is also the weight of the site's own
#' estimate in its posterior mean.
#'
#' @param ... Must be empty. The standard errors and the sampling variances
#'   are named arguments, because an unnamed vector could hold either.
#' @param se Standard errors of the site estimates, a numeric vector of
#'   positive numbers. Give `se` or `se2`, not both.
#' @param se2 Sampling variances of the site estimates, a numeric vector of
#'   positive numbers.
#' @param sigma_tau Standard deviation of the true effects across sites, a
#'   number of at least 0.
#'
#' @details
#' The reliability of site \eqn{j} is
#' \deqn{S_j = \frac{\sigma_\tau^2}{\sigma_\tau^2 + v_j},}{S_j = sigma_tau^2 / (sigma_tau^2 + v_j),}
#' where \eqn{v_j}{v_j} is its sampling variance. It equals
#' \eqn{\mathrm{Var}(\tau_j)/\mathrm{Var}(\hat\tau_j)}{Var(tau_j) / Var(tau_hat_j)}
#' whenever the estimate is unbiased given the effect, without normality or
#' independence across sites (Lee, 2026, Proposition 1).
#'
#' With `sigma_tau = 0` every reliability is 0. The functions of the
#' reliabilities, such as [kl_centroid()], need them above 0, and
#' [reliability_profile()] gives the limits of the diagnostics as
#' `sigma_tau` goes to 0.
#'
#' @return A numeric vector of reliabilities, with the names of `se` or `se2`.
#'
#' @template ref-lee2026
#' @examples
#' site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#' @export
site_reliability <- function(..., se = NULL, se2 = NULL, sigma_tau) {
  .rp_check_named_variances(...)
  if (missing(sigma_tau)) .rp_abort_missing_sigma_tau()
  sigma_tau <- .rp_check_sigma_tau(sigma_tau)
  v <- .rp_variances(se, se2)
  if (sigma_tau == 0) {
    .rp_inform(
      "With `sigma_tau = 0` every reliability is 0.",
      info = "Functions such as `kl_centroid()` need reliabilities above 0; `reliability_profile()` gives the limits of the diagnostics.",
      subclass = "rp_message_sigma_zero"
    )
    return(stats::setNames(rep(0, length(v)), names(v)))
  }
  stats::setNames(.rp_plogis(.rp_logit_from_variances(sigma_tau, v)), names(v))
}

# logit S_j = log(sigma^2) - log(v_j), for sigma > 0.
.rp_logit_from_variances <- function(sigma_tau, v) {
  2 * log(sigma_tau) - log(v)
}

# The logits of the reliabilities split into their mean m and the centered
# deviations d. A second pass makes the deviations sum to zero to rounding.
.rp_center <- function(eta) {
  m <- mean(eta)
  d <- eta - m
  m <- m + mean(d)
  list(m = m, d = eta - m)
}

# The centered logits of a set of reliabilities.
.rp_center_S <- function(S) {
  .rp_center(stats::qlogis(S))
}

# The centered logits of a set of sampling variances: d_j = mean(log v) -
# log(v_j). The logit of I at sigma is 2 log(sigma) - mean(log v).
.rp_center_v <- function(v) {
  lv <- log(v)
  out <- .rp_center(-lv)
  list(log_gm = -out$m, d = out$d)
}

# log(1 + exp(x)) without overflow or loss of precision.
.rp_softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}
