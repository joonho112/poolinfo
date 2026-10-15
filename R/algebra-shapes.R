#' Sampling variances with a given informativeness and logit variance
#'
#' Builds sampling variances whose reliability profile at a given
#' between-site standard deviation has exactly the informativeness `I` and
#' the logit variance `V`, in one of four shapes. These are the profiles of
#' the simulation study of Lee (2026, Appendix F).
#'
#' @param J Number of sites, an integer of at least 1.
#' @param I Informativeness of the profile, a number strictly between 0
#'   and 1.
#' @param V Logit variance of the profile, a number of at least 0. With
#'   `J = 1` it must be 0.
#' @param shape Shape of the profile: one of `"even"` (the default),
#'   `"bimodal"`, `"concentrated_precise"` or `"concentrated_noisy"`. See
#'   Details.
#' @param sigma_tau Between-site standard deviation at which the profile has
#'   the informativeness `I`, a positive number.
#'
#' @details
#' The variances are \eqn{v_j = \tilde v\, e^{x_j}}{v_j = v0 exp(x_j)} with
#' \eqn{\tilde v = \sigma_\tau^2 (1 - I) / I}{v0 = sigma_tau^2 (1 - I) / I}
#' and deviations \eqn{x_j} of mean 0 and mean square \eqn{V}, so that the
#' logits of the reliabilities have mean \eqn{\mathrm{logit}(I)}{logit(I)}
#' and variance \eqn{V}. The shapes of the deviations:
#'
#' * `"even"`: equally spaced, a ladder of variances that is even on the log
#'   scale.
#' * `"bimodal"`: two groups of \eqn{\lceil J/2 \rceil}{ceiling(J/2)} precise
#'   and \eqn{\lfloor J/2 \rfloor}{floor(J/2)} noisy sites, which are not
#'   symmetric when \eqn{J} is odd.
#' * `"concentrated_precise"`: \eqn{k = \lceil J/5 \rceil}{k = ceiling(J/5)}
#'   precise sites at \eqn{-\sqrt{V(J-k)/k}}{-sqrt(V (J - k) / k)} and the
#'   other \eqn{J - k} at \eqn{\sqrt{Vk/(J-k)}}{sqrt(V k / (J - k))}, so a
#'   fifth of the sites have reliabilities far above the rest.
#' * `"concentrated_noisy"`: the mirror image, with \eqn{k} noisy sites.
#'
#' At `V = 0` all four shapes give equal variances. The sites come in the
#' order of their deviations as listed: from precise to noisy, except for
#' `"concentrated_noisy"`, whose noisy sites come first.
#'
#' @return A numeric vector of `J` sampling variances.
#'
#' @template ref-lee2026
#' @examples
#' v <- shape_profile(J = 40, I = 0.3, V = 1, shape = "concentrated_precise")
#' S <- site_reliability(se2 = v, sigma_tau = 1)
#' kl_centroid(S)
#' logit_moments(S)
#' @export
shape_profile <- function(J, I, V, shape = c("even", "bimodal", "concentrated_precise",
                                             "concentrated_noisy"),
                          sigma_tau = 1) {
  J <- .rp_check_whole(J, "J")
  I <- .rp_check_number(I, "I", "a single number strictly between 0 and 1", 0, 1,
                        lower_open = TRUE, upper_open = TRUE)
  V <- .rp_check_number(V, "V", "a single number of at least 0", lower = 0)
  shape <- .rp_match(shape, c("even", "bimodal", "concentrated_precise", "concentrated_noisy"),
                     "shape")
  sigma_tau <- .rp_check_number(sigma_tau, "sigma_tau", "a single positive number",
                                lower = 0, lower_open = TRUE)
  if (J == 1L && V > 0) {
    .rp_abort(
      "A single site has a logit variance of 0.",
      info = sprintf("You asked for `J = 1` and `V = %s`.", format(V))
    )
  }
  x <- if (V == 0) rep(0, J) else .rp_shape_deviations(J, V, shape)
  # log(v0) = log(sigma_tau^2) - logit(I).
  v <- exp(2 * log(sigma_tau) - stats::qlogis(I) + x)
  if (any(!is.finite(v) | v == 0)) {
    .rp_abort(
      "Some of the variances are too large or too small to be represented.",
      info = sprintf("Their logarithms range from %s to %s.",
                     format(signif(min(2 * log(sigma_tau) - stats::qlogis(I) + x), 4)),
                     format(signif(max(2 * log(sigma_tau) - stats::qlogis(I) + x), 4))),
      fix = "Use a smaller `V` or a `sigma_tau` closer to 1."
    )
  }
  v
}

# Deviations of mean 0 and mean square V for J >= 2 sites.
.rp_shape_deviations <- function(J, V, shape) {
  if (shape %in% c("concentrated_precise", "concentrated_noisy")) {
    k <- ceiling(J / 5)
    few <- sqrt(V * (J - k) / k)
    rest <- sqrt(V * k / (J - k))
    if (shape == "concentrated_precise") {
      return(c(rep(-few, k), rep(rest, J - k)))
    }
    return(c(rep(few, k), rep(-rest, J - k)))
  }
  base <- switch(shape,
    even = seq(0, 1, length.out = J),
    bimodal = sort(rep(c(-1, 1), length.out = J))
  )
  base <- base - mean(base)
  base / sqrt(mean(base^2)) * sqrt(V)
}
