#' Summaries of a profile as a data frame
#'
#' Returns the diagnostics of a reliability profile, or their posterior
#' means and intervals, as a data frame with one row for each quantity.
#'
#' @param x An object of class `reliability_profile`,
#'   `reliability_posterior`, `reliability_draws`, `reliability_design` or
#'   `profile_comparison`.
#' @param level Probability content of the posterior intervals, a number
#'   between 0 and 1. `NULL`, the default, uses the level of the fit. Not
#'   used for a profile at a fixed between-site SD.
#' @param ... Must be empty.
#'
#' @details
#' The quantities are coded as in the rest of the package: `I`
#' (informativeness), `mean` (mean reliability), `harmonic` (harmonic
#' centroid), `H` (dispersion), `V` (logit variance), `J_eff`
#' (concentration count) and `sigma_tau` (between-site standard deviation).
#' A posterior fit or draws add `sigma_tau2` and either `effect_mean`, the
#' mean of the effects, or `coef_<name>` for each coefficient of a mean
#' model. For a profile at a fixed `sigma_tau`, the column `mean` contains
#' the value of each quantity and the other columns are `NA`.
#'
#' @return A data frame with one row per quantity and the columns
#'   * `quantity`: the code of the quantity.
#'   * `mean`, `sd`, `lower`, `upper`: the posterior mean and standard
#'     deviation and the ends of the central posterior interval.
#'   * for draws, `mcse` and `ess`: the Monte Carlo standard error of the
#'     mean and the effective sample size of the draws.
#'   * for a design, `allocation`: the planned allocation, or the balanced one
#'     (equal spending with costs).
#'   * for a comparison, `analysis`, the name of the analysis, and `mcse`.
#'
#' @seealso [site_table()] for the sites, and [report_table()] for a table
#'   for a report.
#'
#' @examples
#' prof <- reliability_profile(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#' profile_summary(prof)
#' @export
profile_summary <- function(x, level = NULL, ...) {
  UseMethod("profile_summary")
}

#' @export
profile_summary.default <- function(x, level = NULL, ...) {
  .rp_abort(
    sprintf("`profile_summary()` has no method for %s.", .rp_describe(x)),
    info = "It takes a profile, a posterior fit, draws, a design or a comparison of these."
  )
}

#' @export
profile_summary.reliability_profile <- function(x, level = NULL, ...) {
  .rp_check_dots_empty(...)
  st <- x$stats
  quantity <- c("I", "mean", "harmonic", "H", "V", "J_eff", "sigma_tau")
  value <- c(st[c("I", "mean", "harmonic", "H", "V", "J_eff")], sigma_tau = x$sigma_tau)
  data.frame(quantity = quantity, mean = unname(value), sd = NA_real_, lower = NA_real_,
             upper = NA_real_)
}

#' Site table of a profile
#'
#' Returns one row for each site: its estimate and standard error, its
#' reliability, its term of the dispersion \eqn{H}, its proportion of the total
#' reliability, and, when the estimates are known, the posterior mean and
#' standard deviation of its effect.
#'
#' @param x An object of class `reliability_profile`,
#'   `reliability_posterior`, `reliability_draws` or `reliability_design`.
#' @inheritParams profile_summary
#'
#' @details
#' For a profile at a fixed `sigma_tau`, the effect's posterior mean and
#' standard deviation are those of the normal model at that `sigma_tau`,
#' with a flat prior on the mean of the effects:
#' \eqn{S_j \hat\tau_j + (1 - S_j)\hat\mu}{S_j tau_hat_j + (1 - S_j) mu_hat} and
#' \eqn{\{\sigma_\tau^2(1 - S_j) + (1 - S_j)^2 / W\}^{1/2}}{sqrt(sigma_tau^2 (1 - S_j) + (1 - S_j)^2 / W)},
#' where \eqn{\hat\mu}{mu_hat} is the precision-weighted mean of the
#' estimates with weights \eqn{w_j = 1/(\sigma_\tau^2 + v_j)}{w_j = 1 / (sigma_tau^2 + v_j)}
#' and \eqn{W = \sum_j w_j}{W = sum(w_j)}.
#' A profile inherited from a fitted mean model uses the site's fitted
#' value \eqn{x_j'\hat\beta}{x_j' beta_hat} in place of the common mean,
#' and \eqn{x_j'(X'\Omega X)^{-1}x_j}{x_j' (X' Omega X)^(-1) x_j} in place of
#' \eqn{1/W}{1 / W}, where \eqn{\Omega}{Omega} is the diagonal matrix of
#' the weights \eqn{w_j}{w_j}.
#'
#' The column `share` is the proportion of the total reliability at site
#' \eqn{j}, \eqn{S_j / \sum_k S_k}{S_j / sum(S_k)}. At a fixed `sigma_tau`
#' the inverse of the sum of the squared proportions is the concentration
#' count. For a posterior fit the columns are posterior means, and that
#' identity is only approximate.
#'
#' @return A data frame with one row per site and the columns
#'   * `site`, `estimate` (when known) and `se`.
#'   * for a design, `n`, `kappa` and `n_balanced`: the planned size, the
#'     variance constant and the size under balance or equal spending.
#'   * `reliability`: the site's reliability, its posterior mean for a
#'     posterior fit or draws, which also have `reliability_lower` and
#'     `reliability_upper`, the ends of its posterior interval.
#'   * `kl_term`: the site's term of \eqn{H}; `share`: its proportion of the
#'     total reliability.
#'   * when the estimates are known, `effect` and `effect_sd`: the
#'     posterior mean and standard deviation of the site's effect, and for a
#'     posterior fit or draws `effect_lower` and `effect_upper`, the ends of
#'     its posterior interval.
#'
#' @seealso [profile_summary()] for the diagnostics, and [plot_pooling()]
#'   and [plot_profile()] for plots of the sites.
#'
#' @examples
#' prof <- reliability_profile(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2,
#'                             estimate = c(0.3, 0.1, -0.1, 0.5))
#' site_table(prof)
#' @export
site_table <- function(x, level = NULL, ...) {
  UseMethod("site_table")
}

#' @export
site_table.default <- function(x, level = NULL, ...) {
  .rp_abort(
    sprintf("`site_table()` has no method for %s.", .rp_describe(x)),
    info = "It takes a profile, a posterior fit, draws or a design."
  )
}

#' @export
site_table.reliability_profile <- function(x, level = NULL, ...) {
  .rp_check_dots_empty(...)
  v <- x$se2
  J <- length(v)
  out <- data.frame(site = x$site)
  if (!is.null(x$estimate)) {
    out$estimate <- x$estimate
  }
  out$se <- sqrt(v)
  out$reliability <- x$S
  if (x$sigma_tau == 0) {
    out$kl_term <- 0
    w <- min(v) / v
    out$share <- w / sum(w)
  } else {
    centered <- .rp_center_v(v)
    m <- 2 * log(x$sigma_tau) - centered$log_gm
    out$kl_term <- .rp_bern_kl(rep(m, J), centered$d)
    log_s <- stats::plogis(m + centered$d, log.p = TRUE)
    r <- exp(log_s - max(log_s))
    out$share <- r / sum(r)
  }
  if (!is.null(x$estimate)) {
    cond <- .rp_conditional_sites(x$sigma_tau, x$estimate, v, x$X)
    out$effect <- cond$mean
    out$effect_sd <- sqrt(cond$var)
  }
  out
}

# Posterior mean and variance of each site's effect at a fixed sigma_tau,
# with the mean of the effects integrated out under a flat prior.
.rp_conditional_sites <- function(sigma_tau, y, v, X = NULL) {
  if (!is.null(X)) {
    cond <- .rp_site_conditional(.rp_kernel(y, v, X), sigma_tau)
    return(list(mean = as.vector(cond$mean), var = as.vector(cond$var)))
  }
  log_a <- 2 * log(sigma_tau)
  log_v <- log(v)
  log_den <- pmax(log_a, log_v) + log1p(exp(-abs(log_a - log_v)))
  log_scale <- min(log_den)
  w <- exp(log_scale - log_den)
  W <- sum(w)
  offset <- y[[1L]]
  yc <- y - offset
  mu <- sum(w * yc) / W
  S <- .rp_plogis(log_a - log_v)
  log_one_minus <- stats::plogis(log_v - log_a, log.p = TRUE)
  mean_uncertainty <- exp(2 * log_one_minus + log_scale - log(W))
  list(mean = offset + mu + S * (yc - mu), var = v * S + mean_uncertainty,
       mu = offset + mu, W = exp(log(W) - log_scale))
}

#' @export
profile_summary.reliability_posterior <- function(x, level = NULL, ...) {
  .rp_check_dots_empty(...)
  level <- .rp_check_level(level, null_ok = TRUE)
  out <- x$summary
  if (!is.null(level) && !isTRUE(all.equal(level, x$level))) {
    kernel <- .rp_fit_kernel(x)
    law <- .rp_fit_law(x, kernel)
    probs <- c((1 - level) / 2, (1 + level) / 2)
    sigma_q <- law$quantile(probs)
    fq <- .rp_functionals_at(kernel, sigma_q)
    branches <- .rp_H_branches(kernel, .rp_turning_points(kernel$d))
    H_q <- .rp_H_quantile(law, kernel, branches, probs)
    coef_q <- matrix(.rp_coef_intervals(law, kernel, probs), nrow = kernel$p)
    for (side in 1:2) {
      col <- c("lower", "upper")[[side]]
      out[[col]][match(c("I", "mean", "harmonic", "J_eff"), out$quantity)] <-
        fq[side, c("I", "mean", "harmonic", "J_eff")]
      out[[col]][out$quantity == "H"] <- H_q[[side]]
      out[[col]][out$quantity == "sigma_tau"] <- sigma_q[[side]]
      out[[col]][out$quantity == "sigma_tau2"] <- sigma_q[[side]]^2
      out[[col]][grepl("^(effect_mean|coef_)", out$quantity)] <- coef_q[, side]
    }
  }
  attr(out, "level") <- if (is.null(level)) x$level else level
  attr(out, "prior") <- .rp_prior_text(x$prior, x$data$X)
  attr(out, "conditioning") <- paste(.rp_conditioning_lines(), collapse = " ")
  out
}

#' @export
site_table.reliability_posterior <- function(x, level = NULL, ...) {
  .rp_check_dots_empty(...)
  level <- .rp_check_level(level, null_ok = TRUE)
  if (is.null(level)) level <- x$level
  kernel <- .rp_fit_kernel(x)
  law <- .rp_fit_law(x, kernel)
  out <- x$sites
  if (!isTRUE(all.equal(level, x$level))) {
    sigma_q <- law$quantile(c((1 - level) / 2, (1 + level) / 2))
    out$reliability_lower <- sigma_q[[1L]]^2 / (sigma_q[[1L]]^2 + kernel$v)
    out$reliability_upper <- sigma_q[[2L]]^2 / (sigma_q[[2L]]^2 + kernel$v)
  }
  nodes <- .rp_positive_nodes(law$second()$rule)
  out$share <- .rp_site_shares(kernel, nodes$sigma, nodes$weight / sum(nodes$weight), x$control$chunk_sites)
  intervals <- .rp_site_intervals(law, kernel, level, x$control$chunk_sites)
  out$effect_lower <- intervals[, 1L]
  out$effect_upper <- intervals[, 2L]
  out <- out[, c("site", "estimate", "se", "reliability", "reliability_lower", "reliability_upper",
                 "kl_term", "share", "effect", "effect_sd", "effect_lower", "effect_upper")]
  attr(out, "level") <- level
  out
}

#' Posterior probability that a diagnostic exceeds a value
#'
#' Computes the posterior probability that the informativeness, the mean
#' reliability, the concentration count or the between-site standard
#' deviation is above a given value.
#'
#' @param x An object of class `reliability_posterior` or
#'   `reliability_draws`.
#' @param threshold The value, a single number. For `which = "I"`, `NULL`
#'   (the default) uses the `threshold` of the fit; for the other quantities
#'   a value must be given.
#' @param which The quantity: one of `"I"` (the default), `"mean"`,
#'   `"J_eff"` or `"sigma_tau"`.
#'
#' @details
#' The four quantities increase with `sigma_tau`, so the probability is
#' \eqn{1 - F(\sigma_c)}{1 - F(sigma_c)}, where \eqn{F} is the posterior
#' distribution function of `sigma_tau` and \eqn{\sigma_c}{sigma_c} the value
#' at which the quantity equals the threshold. For the informativeness and a
#' threshold \eqn{c},
#' \eqn{\sigma_c = \{\mathrm{GM}(v)\, c/(1-c)\}^{1/2}}{sigma_c = sqrt(GM(v) c / (1 - c))}.
#' When all sampling variances are equal, the concentration count is the
#' constant \eqn{J}{J}; its probability of strictly exceeding \eqn{c}{c}
#' is 1 for \eqn{c < J}{c < J} and 0 for \eqn{c \ge J}{c >= J}.
#'
#' Lee (2026) reports the probability that \eqn{I} exceeds 1/2, the value at
#' which the between-site variance equals the geometric mean of the sampling
#' variances. Above it the leading term of the mean reliability minus
#' \eqn{I} is negative (Proposition 2), although the remainder can change
#' the sign of the difference. The value 1/2 is a reference point, not a
#' standard that an analysis must meet.
#'
#' @return A number between 0 and 1. For draws, its attribute `mcse` is the
#'   Monte Carlo standard error.
#'
#' @template ref-lee2026
#' @seealso [plot_functionals()], whose label gives this probability.
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools)
#' prob_above(fit)
#' prob_above(fit, threshold = 0.3)
#' prob_above(fit, threshold = 5, which = "sigma_tau")
#' @export
prob_above <- function(x, threshold = NULL, which = c("I", "mean", "J_eff", "sigma_tau")) {
  UseMethod("prob_above")
}

#' @export
prob_above.default <- function(x, threshold = NULL, which = c("I", "mean", "J_eff", "sigma_tau")) {
  .rp_abort(sprintf("`prob_above()` has no method for %s.", .rp_describe(x)),
            info = "It takes a posterior fit or draws.")
}

#' @export
prob_above.reliability_posterior <- function(x, threshold = NULL,
                                             which = c("I", "mean", "J_eff", "sigma_tau")) {
  which <- .rp_match(which, c("I", "mean", "J_eff", "sigma_tau"), "which")
  kernel <- .rp_fit_kernel(x)
  threshold <- .rp_prob_threshold(threshold, which, x$threshold, kernel$J)
  if (which == "J_eff" && all(kernel$v == kernel$v[[1L]])) {
    return(as.numeric(threshold < kernel$J))
  }
  sigma_c <- .rp_sigma_at_value(kernel, which, threshold)
  if (is.na(sigma_c)) {
    return(if (threshold <= .rp_value_at_zero(kernel, which)) 1 else 0)
  }
  law <- .rp_fit_law(x, kernel)
  law$cdf(sigma_c, lower = FALSE)
}

# The threshold of prob_above(), checked against the range of the quantity.
.rp_prob_threshold <- function(threshold, which, default, J) {
  if (is.null(threshold)) {
    if (which != "I") {
      .rp_abort(
        sprintf("`threshold` is needed for `which = \"%s\"`.", which),
        info = "The default threshold of a fit applies to the informativeness only."
      )
    }
    return(default)
  }
  switch(which,
    I = ,
    mean = .rp_check_level(threshold, "threshold"),
    J_eff = .rp_check_number(threshold, "threshold", sprintf("a single number between 1 and %d", J),
                             lower = 1, upper = J),
    sigma_tau = .rp_check_number(threshold, "threshold", "a single number of at least 0", lower = 0)
  )
}

# The value of an increasing quantity at sigma_tau = 0.
.rp_value_at_zero <- function(kernel, which) {
  switch(which,
    I = ,
    mean = ,
    sigma_tau = 0,
    J_eff = .rp_functionals_at(kernel, 0)[, "J_eff"]
  )
}

# The sigma_tau at which an increasing quantity equals a value, or NA when
# the value is at or below the quantity's value at sigma_tau = 0.
.rp_sigma_at_value <- function(kernel, which, value) {
  if (which == "sigma_tau") {
    return(value)
  }
  if (which == "I") {
    return(exp((stats::qlogis(value) + kernel$log_gm) / 2))
  }
  if (which == "J_eff" && value >= kernel$J) {
    return(Inf)
  }
  if (value <= .rp_value_at_zero(kernel, which)) {
    return(NA_real_)
  }
  f <- function(t) .rp_functionals_at(kernel, exp(t))[, which] - value
  lo <- log(kernel$scale) - 1
  while (f(lo) > 0) lo <- lo - 2
  hi <- log(kernel$scale) + 1
  while (f(hi) < 0) {
    hi <- hi + 2
    if (hi > log(kernel$scale) + 200) return(Inf)
  }
  exp(stats::uniroot(f, c(lo, hi), tol = 1e-13)$root)
}

#' Posterior means against plug-in values
#'
#' Compares the posterior mean of each diagnostic with its value at the
#' posterior mean of \eqn{\sigma_\tau^2}{sigma_tau^2}, the value that a
#' plug-in of that estimate would give.
#'
#' @param x An object of class `reliability_posterior` or
#'   `reliability_draws`.
#'
#' @details
#' The informativeness and every site reliability are concave functions of
#' \eqn{\sigma_\tau^2}{sigma_tau^2}. By Jensen's inequality, the posterior
#' means of the informativeness and of the mean reliability are therefore at
#' most their values at the posterior mean of \eqn{\sigma_\tau^2}{sigma_tau^2},
#' and strictly less unless the posterior of \eqn{\sigma_\tau^2}{sigma_tau^2}
#' is a point mass (Lee, 2026, Proposition 3). The dispersion and the
#' concentration count have no such general order. A plug-in of another
#' estimate, such as the REML estimate, can fall on either side of the
#' posterior mean.
#'
#' The comparison needs the posterior mean of \eqn{\sigma_\tau^2}{sigma_tau^2},
#' which under the flat prior exists only with more than 4 sites (more than
#' \eqn{p + 3} with a mean model of \eqn{p} columns). Without it the
#' values are `NA`, with a message that gives the number of sites needed.
#'
#' @return A data frame with the columns `functional` (`"I"`, `"mean"`,
#'   `"H"` and `"J_eff"`), `posterior_mean`, `at_mean_sigma2` and
#'   `difference` (the second minus the first), and the attribute `sigma2`,
#'   the posterior mean of \eqn{\sigma_\tau^2}{sigma_tau^2}.
#'
#' @template ref-lee2026
#' @seealso `vignette("posterior-summaries")` for the comparison in an
#'   example.
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools)
#' plugin_gap(fit)
#' @export
plugin_gap <- function(x) {
  UseMethod("plugin_gap")
}

#' @export
plugin_gap.default <- function(x) {
  .rp_abort(sprintf("`plugin_gap()` has no method for %s.", .rp_describe(x)),
            info = "It takes a posterior fit or draws.")
}

#' @export
plugin_gap.reliability_posterior <- function(x) {
  names <- c("I", "mean", "H", "J_eff")
  post <- vapply(names, function(q) .rp_summary_value(x$summary, q, "mean"), 0)
  plug <- x$plugin[names]
  if (is.na(x$plugin[["sigma2"]])) {
    .rp_inform(
      sprintf("The posterior mean of sigma_tau^2 does not exist with %d sites under this prior.",
              length(x$data$estimate)),
      info = sprintf("It needs at least %d sites.",
                     .rp_moment_min_sites(2, if (is.null(x$data$X)) 1L else ncol(x$data$X), x$prior)),
      subclass = "rp_message_no_moment"
    )
  }
  out <- data.frame(functional = names, posterior_mean = unname(post),
                    at_mean_sigma2 = unname(plug), difference = unname(plug - post))
  attr(out, "sigma2") <- x$plugin[["sigma2"]]
  out
}

#' Numerical checks of a posterior fit
#'
#' Reports the checks of the integration over the between-site standard
#' deviation: the agreement of the rule with two other rules, the bounds on
#' the omitted upper tail, and the location of the peak of the dispersion
#' found at two resolutions.
#'
#' @param x An object of class `reliability_posterior`.
#'
#' @details
#' For the node rule, the normalizing constant and the posterior means of
#' \eqn{I} and \eqn{H} are compared with those of a rule whose panels are
#' all halved, and with those of a rule whose lowest panel ends at half the
#' value. The differences are relative, and their limit is
#' \eqn{10^{-10}}{1e-10}, or \eqn{64J}{64 J} times the machine epsilon when
#' that is larger, because the log density is a sum of \eqn{J} terms. For
#' the adaptive engine, the check is the error
#' estimate of [stats::integrate()] for the normalizing constant. The bounds
#' on the omitted upper tail are those of Appendix D.3 of Lee (2026),
#' extended to the mean model and to the priors of [sd_prior()]; the bound
#' for \eqn{E(\sigma_\tau^q)}{E(sigma_tau^q)} is shown divided by
#' \eqn{\sqrt{\mathrm{GM}(v)}^q}{sqrt(GM(v))^q}, so that it does not depend on
#' the units. The bounds concern the numerical approximation, not the
#' calibration of the posterior. The last row compares the peak of \eqn{H},
#' where the mean reliability equals \eqn{I}, as located by scans of the
#' logit of \eqn{I} at two resolutions. `vignette("implementation")` has
#' the details.
#'
#' @return A data frame with the columns `check`, `value`, `limit` and
#'   `status` (`"ok"` or `"outside"`).
#'
#' @template ref-lee2026
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools)
#' check_numerics(fit)
#' @export
check_numerics <- function(x) {
  if (!inherits(x, "reliability_posterior")) {
    .rp_abort("`check_numerics()` takes a fit of `reliability_posterior()`.",
              info = .rp_you_passed(x))
  }
  x$checks
}

#' Posterior density of a diagnostic
#'
#' Computes the exact posterior density of a diagnostic or of the
#' between-site standard deviation on a grid of values.
#'
#' @param x An object of class `reliability_posterior`.
#' @param which The quantity: one of `"I"` (the default), `"mean"`, `"H"`,
#'   `"J_eff"` or `"sigma_tau"`.
#' @param n Number of points of the grid, an integer of at least 10.
#'
#' @details
#' The density follows from that of `sigma_tau` by a change of variables:
#' for an increasing quantity \eqn{f}, the density at \eqn{f(\sigma)} is
#' \eqn{g(\sigma) / f'(\sigma)}{g(sigma) / f'(sigma)}. For the dispersion
#' \eqn{H} the terms of the two branches on either side of its largest value
#' are added, with \eqn{dH/d\sigma = 2(\bar S - I)/\sigma}{dH/dsigma = 2 (mean(S_j) - I) / sigma}.
#' Under every prior of [sd_prior()] the posterior density of `sigma_tau`
#' has a finite positive limit at 0. The densities of \eqn{I} and the mean
#' reliability grow without bound at their values for `sigma_tau = 0`,
#' like the inverse square root of the distance from that value; they are
#' still integrable. This behavior also applies to \eqn{H} and the
#' concentration count when the standard errors are unequal. With equal
#' standard errors, their posteriors are point masses at 0 and \eqn{J},
#' respectively.
#'
#' The grid spans the central 99.9% of the posterior, so a plot of the
#' density does not include the extreme tails.
#' With equal sampling variances, \eqn{J_{\mathrm{eff}}}{J_eff} is a point
#' mass at \eqn{J}{J}, has no density, and asking for one is an error. When
#' the sampling variances differ so little (by about 1e-6 relative or less)
#' that the derivative of \eqn{J_{\mathrm{eff}}}{J_eff} cannot be computed
#' reliably, asking for its density is an error too.
#'
#' @return A data frame with the columns `value` and `density`.
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools)
#' head(functional_density(fit, "I", n = 50))
#' @export
functional_density <- function(x, which = c("I", "mean", "H", "J_eff", "sigma_tau"), n = 512) {
  if (!inherits(x, "reliability_posterior")) {
    .rp_abort("`functional_density()` takes a fit of `reliability_posterior()`.",
              info = .rp_you_passed(x))
  }
  which <- .rp_match(which, c("I", "mean", "H", "J_eff", "sigma_tau"), "which")
  n <- .rp_check_whole(n, "n", lower = 10)
  kernel <- .rp_fit_kernel(x)
  if (which == "J_eff" && all(kernel$v == kernel$v[[1L]])) {
    .rp_abort(
      sprintf("J_eff is a point mass at %d and has no ordinary density.", kernel$J),
      info = "The sampling variances are all equal.",
      fix = "Use `prob_above()` for its strict exceedance probability.",
      class = "rp_scope_error"
    )
  }
  law <- .rp_fit_law(x, kernel)
  if (which == "H") {
    return(.rp_H_density(x, kernel, law, n))
  }
  ends <- law$quantile(c(0.0005, 0.9995))
  if (which == "sigma_tau") {
    s <- seq(ends[[1L]], ends[[2L]], length.out = n)
    return(data.frame(value = s, density = law$density(s)))
  }
  f_ends <- .rp_functionals_at(kernel, ends)[, which]
  value <- seq(f_ends[[1L]], f_ends[[2L]], length.out = n)
  if (which == "J_eff" &&
      (diff(f_ends) <= 64 * .Machine$double.eps * max(abs(f_ends)) || anyDuplicated(value))) {
    .rp_abort_concentration_density()
  }
  sigma <- vapply(value, function(z) .rp_sigma_at_value(kernel, which, z), 0)
  deriv <- .rp_functional_derivative(kernel, sigma, which)
  density <- law$density(sigma) / deriv
  if (which == "J_eff" && any(!is.finite(density) | !is.finite(deriv) | deriv <= 0 | density < 0)) {
    .rp_abort_concentration_density()
  }
  data.frame(value = value, density = density)
}

.rp_abort_concentration_density <- function() {
  .rp_abort(
    "The variation of J_eff is too small to resolve its density at machine precision.",
    info = "The concentration count is nonconstant, but the continuous density calculation loses precision.",
    fix = "Use its posterior summary or draws instead of an ordinary density curve.",
    class = "rp_numeric_error"
  )
}

# d f / d sigma for the increasing functionals. With m = 2 log(sigma) -
# log GM(v), dm / dsigma = 2 / sigma, and dS_j / dm = S_j (1 - S_j).
.rp_functional_derivative <- function(kernel, sigma, which) {
  m <- 2 * log(sigma) - kernel$log_gm
  dm <- 2 / sigma
  switch(which,
    I = stats::plogis(m) * stats::plogis(-m) * dm,
    mean = .rp_by_rows(m, kernel$d, function(eta, dd, mm, n) {
      .rp_row_means(stats::plogis(eta) * stats::plogis(-eta), n)
    }) * dm,
    J_eff = .rp_by_rows(m, kernel$d, function(eta, dd, mm, n) {
      # The centered weighted variance avoids subtracting nearly equal
      # terms in a * sum(S^3) - b^2. Scale S before taking powers.
      top <- stats::plogis(mm + max(dd), log.p = TRUE)
      r <- matrix(exp(stats::plogis(eta, log.p = TRUE) - top), nrow = n)
      a <- rowSums(r)
      b <- rowSums(r^2)
      2 * exp(top[seq_len(n)]) * a^2 / b^2 * rowSums(r * (r - b / a)^2)
    }) * dm
  )
}

# The density of H, as the sum over its branches of g(sigma) / |dH/dsigma|.
.rp_H_density <- function(x, kernel, law, n) {
  hb <- x$h_branches
  if (hb$homogeneous || hb$H_max == 0) {
    .rp_abort(
      "H is 0 at every sigma_tau, so it has no density.",
      info = "The standard errors are all equal."
    )
  }
  branches <- .rp_H_branches(kernel, .rp_turning_points(kernel$d))
  ends <- .rp_H_quantile(law, kernel, branches, c(0.0005, 0.9995))
  value <- seq(ends[[1L]], ends[[2L]], length.out = n)
  density <- vapply(value, function(h) {
    total <- 0
    for (k in seq_along(branches$rising)) {
      m_a <- branches$m[[k]]
      m_b <- branches$m[[k + 1L]]
      h_a <- branches$h[[k]]
      h_b <- branches$h[[k + 1L]]
      if (h > min(h_a, h_b) && h < max(h_a, h_b)) {
        m <- .rp_H_inverse(kernel$d, h, m_a, m_b, h_a, h_b)
        s <- .rp_sigma_of_m(kernel, m)
        total <- total + law$density(s) / abs(.rp_gap_at(m, kernel$d) * 2 / s)
      }
    }
    total
  }, 0)
  data.frame(value = value, density = density)
}

#' Draws from the posterior
#'
#' Draws the between-site standard deviation, the mean of the effects and
#' the diagnostics from their posterior, and optionally the site effects.
#'
#' @param x An object of class `reliability_posterior`.
#' @param ndraws Number of draws, a positive integer.
#' @param seed Seed of the random number generator, an integer, or `NULL`
#'   (the default) for a seed from the clock. The seed used is stored in
#'   the attribute `seed`. The state of the session's generator is restored
#'   afterwards, so `set.seed()` before the call does not make the draws
#'   reproducible; give `seed` instead.
#' @param sites If `TRUE`, the draws include the site effects.
#'
#' @details
#' `sigma_tau` is drawn by inverting a table of its distribution function,
#' then the mean of the effects from its normal distribution given
#' `sigma_tau`, and the site effects from their normal distributions given
#' both (Lee, 2026, Appendix D.3). The draws are independent. They are for
#' plots and for summaries that the package does not compute; the summaries
#' of the package itself come from the integrals and have no Monte Carlo
#' error.
#'
#' @return A data frame with one row per draw and the columns `draw`,
#'   `sigma_tau`, `effect_mean` (or one column per coefficient of the mean
#'   model), `I`, `mean`, `H`, `J_eff` and, with `sites = TRUE`, one column
#'   per site, named `effect_` followed by the site's label (`effect_A` for
#'   site A).
#'
#' @template ref-lee2026
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools)
#' draws <- profile_draws(fit, ndraws = 1000, seed = 1)
#' mean(draws$I)
#' @export
profile_draws <- function(x, ndraws = 4000, seed = NULL, sites = FALSE) {
  if (!inherits(x, "reliability_posterior")) {
    .rp_abort("`profile_draws()` takes a fit of `reliability_posterior()`.",
              info = .rp_you_passed(x))
  }
  ndraws <- .rp_check_whole(ndraws, "ndraws")
  sites <- .rp_check_flag(sites, "sites")
  if (is.null(seed)) {
    seed <- .rp_clock_seed()
  } else {
    seed <- .rp_check_whole(seed, "seed", lower = -.Machine$integer.max)
  }
  kernel <- .rp_fit_kernel(x)
  law <- .rp_fit_law(x, kernel)
  out <- .rp_with_seed(seed, .rp_draws(law, kernel, ndraws, sites))
  if (sites) {
    names(out)[grepl("^effect_[0-9]+$", names(out))] <- paste0("effect_", x$data$site)
  }
  attr(out, "seed") <- seed
  out
}

#' @rdname reliability_profile
#' @order 5
#' @export
reliability_profile.reliability_posterior <- function(x, ...) {
  .rp_check_dots_empty(...)
  a <- x$plugin[["sigma2"]]
  if (is.na(a)) {
    .rp_abort(
      "The posterior mean of sigma_tau^2 does not exist for this fit.",
      info = sprintf("Under this prior it needs at least %d sites; there are %d.",
                     .rp_moment_min_sites(2, if (is.null(x$data$X)) 1L else ncol(x$data$X), x$prior),
                     length(x$se2))
    )
  }
  .rp_new_profile(
    se = NULL, se2 = x$data$se2, sigma_tau = sqrt(a), site = x$data$site,
    estimate = x$data$estimate, source = "posterior mean of sigma_tau^2",
    X = x$data$X,
    scale_label = x$labels$scale_label, target_label = x$labels$target_label,
    call = match.call()
  )
}
