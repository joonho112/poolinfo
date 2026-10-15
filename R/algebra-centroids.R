#' Informativeness of a reliability profile
#'
#' Computes the informativeness \eqn{I}, the centroid of the site
#' reliabilities under the Bernoulli Kullback-Leibler divergence: the
#' logistic function of the mean logit of the reliabilities.
#'
#' @param S Site reliabilities, a numeric vector of numbers strictly between
#'   0 and 1, as returned by [site_reliability()].
#'
#' @details
#' The informativeness is
#' \deqn{I = \mathrm{logistic}\Big(\frac{1}{J}\sum_j \mathrm{logit}\, S_j\Big),}{I = logistic(mean(logit(S_j))),}
#' the value \eqn{p} that minimizes the average divergence
#' \eqn{J^{-1}\sum_j \mathrm{KL}\{\mathrm{Ber}(p) \| \mathrm{Ber}(S_j)\}}{mean(KL(Ber(p) || Ber(S_j)))}
#' (Lee, 2026, Theorem 1). When the reliabilities come from sampling
#' variances \eqn{v_j} and a between-site standard deviation \eqn{\sigma_\tau},
#' it equals \eqn{\sigma_\tau^2 / \{\sigma_\tau^2 + \mathrm{GM}(v)\}}{sigma_tau^2 / (sigma_tau^2 + GM(v))},
#' where \eqn{\mathrm{GM}(v)}{GM(v)} is the geometric mean of the variances.
#'
#' A reliability of exactly 0 or 1 has no finite logit, and such a profile has
#' no centroid; the function then fails with an error.
#'
#' @return The informativeness, a number between 0 and 1.
#'
#' @seealso [reliability_centroids()] for the harmonic and arithmetic means
#'   beside \eqn{I}, and [kl_dispersion()] for the spread around it.
#'
#' @template ref-lee2026
#' @examples
#' S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#' kl_centroid(S)
#' @export
kl_centroid <- function(S) {
  S <- .rp_check_reliabilities(S)
  .rp_plogis(.rp_center_S(S)$m)
}

#' Three centroids of a reliability profile
#'
#' Computes the harmonic centroid, the informativeness and the mean
#' reliability of a set of site reliabilities.
#'
#' @inheritParams kl_centroid
#'
#' @details
#' The three numbers are the first-argument centroids of the reliabilities
#' under three Bregman divergences (Lee, 2026, Table 2): the harmonic mean
#' (Burg potential), the informativeness \eqn{I} (Bernoulli
#' Kullback-Leibler) and the arithmetic mean \eqn{\bar S}{mean(S_j)}
#' (squared error). The harmonic centroid is never larger than the other
#' two, and equals them only when all the reliabilities are equal. There is no general order between \eqn{I} and the mean reliability:
#' either can be the larger, depending on the level and the skewness of the
#' profile, and the sign of their difference is reported by [local_terms()].
#'
#' @return A named numeric vector with elements `harmonic`, `I` (the
#'   informativeness) and `mean`.
#'
#' @seealso [bregman_centroid()] for the centroids on either side of each
#'   divergence.
#'
#' @template ref-lee2026
#' @examples
#' S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#' reliability_centroids(S)
#' @export
reliability_centroids <- function(S) {
  S <- .rp_check_reliabilities(S)
  p <- .rp_center_S(S)
  I <- .rp_plogis(p$m)
  # The harmonic centroid is at most the other two; the minimum removes the
  # few units in the last place by which rounding can put it above them.
  c(
    harmonic = min(.rp_harmonic_at(p$m, p$d), I, mean(S)),
    I = I,
    mean = mean(S)
  )
}
