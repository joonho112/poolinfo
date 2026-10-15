#' Typical sampling variance of Higgins and Thompson
#'
#' Computes the "typical" within-study variance of Higgins and Thompson
#' (2002), the variance that the \eqn{I^2} statistic of meta-analysis
#' compares with the between-study variance.
#'
#' @inheritParams site_reliability
#'
#' @details
#' With weights \eqn{w_j = 1 / v_j}, \eqn{W = \sum_j w_j} and
#' \eqn{W_2 = \sum_j w_j^2},
#' \deqn{v_{\mathrm{typ}} = \frac{(J-1)W}{W^2 - W_2}.}{v_typ = (J - 1) W / (W^2 - W_2).}
#' For a random-effects model without moderators, metafor's `rma()` reports
#' \eqn{I^2} as the percentage
#' \eqn{100\,\hat\sigma_\tau^2 / (\hat\sigma_\tau^2 + v_{\mathrm{typ}})}{100 sigma_tau_hat^2 / (sigma_tau_hat^2 + v_typ)}
#' at its estimate \eqn{\hat\sigma_\tau^2}{sigma_tau_hat^2} (\eqn{\tau^2}
#' in metafor's notation). The typical variance is
#' at least the harmonic mean of the \eqn{v_j} (Böhning et al., 2017), and
#' in general it differs from their geometric mean, which defines the
#' informativeness. At least 2 sites are needed.
#'
#' @return The typical variance, a positive number.
#'
#' @references
#' Böhning, D., Lerdsuwansri, R., & Holling, H. (2017). Some general points
#' on the \eqn{I^2}{I^2}-measure of heterogeneity in meta-analysis.
#' *Metrika, 80*, 685-695. \doi{10.1007/s00184-017-0622-3}
#'
#' Higgins, J. P. T., & Thompson, S. G. (2002). Quantifying heterogeneity in
#' a meta-analysis. *Statistics in Medicine, 21*, 1539-1558.
#' \doi{10.1002/sim.1186}
#'
#' @seealso [kl_centroid()] for the informativeness, which uses the
#'   geometric mean of the sampling variances instead.
#'
#' @examples
#' se <- c(0.10, 0.15, 0.20, 0.30)
#' typical_variance(se = se)
#' exp(mean(log(se^2)))  # the geometric mean, for comparison
#' @export
typical_variance <- function(..., se = NULL, se2 = NULL) {
  .rp_check_named_variances(...)
  v <- .rp_variances(se, se2, min_length = 2L)
  # The weights are scaled by the largest, min(v) / v, so that the result
  # scales with the variances without overflow or underflow.
  w <- min(v) / v
  J <- length(w)
  # W^2 - W_2 = sum_j w_j (W - w_j), with W - w_j summed from the other
  # terms, so that a dominant weight does not cancel.
  left <- c(0, cumsum(w)[-J])
  right <- rev(c(0, cumsum(rev(w))[-J]))
  multiplier <- (J - 1) * sum(w) / sum(w * (left + right))
  value <- min(v) * multiplier
  if (is.finite(value) && value > 0) return(min(max(v), value))
  # If relative weights underflow, retain their logarithms. With the weights
  # ordered largest first, each prefix contains 1, so its sum is well scaled.
  # The denominator is twice the sum of w_i w_j over i < j.
  log_w <- sort(log(min(v)) - log(v), decreasing = TRUE)
  prefix <- cumsum(exp(log_w))[-J]
  log_pairs <- log_w[-1L] + log(prefix)
  log_pair_sum <- .rp_log_mean_exp(log_pairs) + log(J - 1)
  log_value <- log(min(v)) + log(J - 1) + log(sum(w)) - log(2) - log_pair_sum
  # The typical variance cannot exceed the largest sampling variance.
  min(max(v), exp(log_value))
}

#' Benchmark risk of a reliability profile
#'
#' Computes \eqn{\sigma_\tau^2 (1 - \bar S)}{sigma_tau^2 (1 - mean(S_j))},
#' the average squared error of the best predictors that are affine in each
#' site's own estimate when the mean and the variances are known.
#'
#' @inheritParams site_reliability
#'
#' @details
#' For a known common mean and known \eqn{\sigma_\tau^2} and \eqn{v_j}, the
#' predictor \eqn{\mu + S_j(\hat\tau_j - \mu)}{mu + S_j (tau_hat_j - mu)}
#' has the smallest mean squared error among predictors affine in the
#' site's own estimate, and its error averaged over sites is
#' \eqn{\sigma_\tau^2(1 - \bar S)}{sigma_tau^2 (1 - mean(S_j))}
#' (Lee, 2026, Proposition 2(iii)). Under the normal model with independent
#' sites it is the Bayes risk of the posterior means with known
#' hyperparameters. It is a benchmark at known inputs, not the error of a
#' fitted model, whose hyperparameters are estimated. It is 0 at
#' `sigma_tau = 0` and approaches the mean of the \eqn{v_j} as `sigma_tau`
#' grows.
#'
#' @return The benchmark risk, a number of at least 0, in the squared units
#'   of the estimates.
#'
#' @template ref-lee2026
#' @examples
#' benchmark_risk(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#' @export
benchmark_risk <- function(..., se = NULL, se2 = NULL, sigma_tau) {
  .rp_check_named_variances(...)
  if (missing(sigma_tau)) .rp_abort_missing_sigma_tau()
  sigma_tau <- .rp_check_sigma_tau(sigma_tau)
  v <- .rp_variances(se, se2)
  mean(.rp_parallel_variances(v, sigma_tau))
}

# v sigma^2 / (v + sigma^2), keeping both squares and reciprocal ratios out
# of intermediate arithmetic when they would exceed the floating-point range.
.rp_parallel_variances <- function(v, sigma_tau) {
  if (sigma_tau == 0) return(rep(0, length(v)))
  tau2 <- sigma_tau^2
  if (is.finite(tau2) && tau2 > .Machine$double.xmin) {
    small <- pmin(v, tau2)
    return(small / (1 + small / pmax(v, tau2)))
  }
  lv <- log(v)
  lt <- 2 * log(sigma_tau)
  pmin(v, exp(pmin(lv, lt) - log1p(exp(-abs(lv - lt)))))
}

#' Bregman centroids and dispersions of a reliability profile
#'
#' `bregman_centroid()` computes the centroid of the site reliabilities
#' under one of three Bregman divergences, on either side of the
#' divergence, and `bregman_dispersion()` the average divergence at that
#' centroid.
#'
#' @inheritParams kl_centroid
#' @param potential The convex function that defines the divergence: one of
#'   `"kl"` (the default), `"squared"` or `"burg"`. See Details.
#' @param side Which argument of the divergence the centroid takes: `"first"`
#'   (the default) or `"second"`.
#'
#' @details
#' A strictly convex function \eqn{\phi} defines the divergence
#' \eqn{D_\phi(x, y) = \phi(x) - \phi(y) - \phi'(y)(x - y)}{D(x, y) = phi(x) - phi(y) - phi'(y) (x - y)}.
#' The first-argument centroid minimizes \eqn{\sum_j D_\phi(c, S_j)}{sum_j D(c, S_j)}
#' and is \eqn{(\phi')^{-1}\{J^{-1}\sum_j \phi'(S_j)\}}{the inverse of phi' at mean(phi'(S_j))};
#' the second-argument centroid minimizes \eqn{\sum_j D_\phi(S_j, c)}{sum_j D(S_j, c)}
#' and is the mean reliability for every \eqn{\phi} (Lee, 2026, Table 2 and
#' Appendix B):
#'
#' | `potential` | \eqn{\phi(x)}{phi(x)} | first-argument centroid |
#' |:---|:---|:---|
#' | `"squared"` | \eqn{x^2/2}{x^2 / 2} | mean reliability |
#' | `"kl"` | \eqn{x\log x + (1-x)\log(1-x)}{x log(x) + (1 - x) log(1 - x)} | informativeness \eqn{I} |
#' | `"burg"` | \eqn{-\log x}{-log(x)} | harmonic mean |
#'
#' The dispersion is the minimum average divergence. On the first side of
#' `"kl"` it is the dispersion \eqn{H} of [kl_dispersion()]; on the second
#' side it is the mean Bernoulli divergence of the reliabilities from their
#' mean. The two need not agree: for the reliabilities 0.2 and 0.8 both
#' centroids are 0.5, but the dispersions are \eqn{\log(5/4) \approx 0.2231}{log(5/4), about 0.2231}
#' and about 0.1927. For `"burg"` the dispersions are the logarithm of the
#' ratio of the geometric to the harmonic mean (first side) and of the
#' arithmetic to the geometric mean (second side); for `"squared"` both are
#' half the variance of the reliabilities, with divisor \eqn{J}.
#'
#' @return A number: the centroid, between 0 and 1, or the dispersion, of
#'   at least 0.
#'
#' @template ref-lee2026
#' @examples
#' S <- c(0.2, 0.8)
#' bregman_centroid(S, "kl")
#' bregman_dispersion(S, "kl", side = "first")
#' bregman_dispersion(S, "kl", side = "second")
#'
#' S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#' sapply(c("squared", "kl", "burg"), function(p) bregman_centroid(S, p))
#' @export
bregman_centroid <- function(S, potential = c("kl", "squared", "burg"),
                             side = c("first", "second")) {
  S <- .rp_check_reliabilities(S)
  potential <- .rp_match(potential, c("kl", "squared", "burg"), "potential")
  side <- .rp_match(side, c("first", "second"), "side")
  if (side == "second" || potential == "squared") {
    return(mean(S))
  }
  p <- .rp_center_S(S)
  switch(potential,
    kl = stats::plogis(p$m),
    burg = .rp_harmonic_at(p$m, p$d)
  )
}

#' @rdname bregman_centroid
#' @export
bregman_dispersion <- function(S, potential = c("kl", "squared", "burg"),
                               side = c("first", "second")) {
  S <- .rp_check_reliabilities(S)
  potential <- .rp_match(potential, c("kl", "squared", "burg"), "potential")
  side <- .rp_match(side, c("first", "second"), "side")
  if (potential == "squared") {
    return(mean((S - mean(S))^2) / 2)
  }
  p <- .rp_center_S(S)
  eta <- p$m + p$d
  if (potential == "kl") {
    if (side == "first") {
      return(.rp_H_at(p$m, p$d))
    }
    # mean_j KL(Ber(S_j) || Ber(mean S)): the divergence from base logit
    # eta_j to the logit of the mean.
    return(mean(.rp_bern_kl(eta, stats::qlogis(mean(S)) - eta)))
  }
  # The logarithms of the means are formed from the logits, which keeps
  # their digits when the reliabilities are close to 1.
  log_gm <- mean(stats::plogis(eta, log.p = TRUE))
  if (side == "first") {
    log_hm <- stats::plogis(p$m - .rp_log_mean_exp(-p$d), log.p = TRUE)
    return(max(log_gm - log_hm, 0))
  }
  log_am <- if (mean(S) > 0.5) log1p(-mean(.rp_plogis(-eta))) else log(mean(S))
  max(log_am - log_gm, 0)
}
