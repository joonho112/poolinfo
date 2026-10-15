#' Dispersion of a reliability profile
#'
#' `kl_dispersion()` computes the dispersion \eqn{H} of the site
#' reliabilities around their informativeness \eqn{I}, and
#' `kl_contributions()` the term of each site, whose mean is \eqn{H}.
#'
#' @inheritParams kl_centroid
#'
#' @details
#' The dispersion is the minimum of the average divergence that defines
#' \eqn{I},
#' \deqn{H = \frac{1}{J}\sum_j \mathrm{KL}\{\mathrm{Ber}(I) \| \mathrm{Ber}(S_j)\}
#'   = \frac{1}{J}\sum_j A(m + \Delta_j) - A(m),}{H = mean(KL(Ber(I) || Ber(S_j))) = mean(A(m + Delta_j)) - A(m),}
#' where \eqn{A(t) = \log(1 + e^t)}{A(t) = log(1 + exp(t))}, \eqn{m} is the
#' logit of \eqn{I} and \eqn{\Delta_j}{Delta_j} the centered logit of site
#' \eqn{j} (Lee, 2026, Theorem 2 and Proposition 2). It is 0 only when all
#' the reliabilities are equal.
#'
#' A site's term is the divergence of \eqn{I} from the site's reliability.
#' Removing the site would also move \eqn{I}, so the term is not the change
#' in \eqn{H} without that site.
#'
#' @return `kl_dispersion()` returns \eqn{H}, a number of at least 0.
#'   `kl_contributions()` returns a numeric vector with one term for each
#'   site, with the names of `S`.
#'
#' @seealso [logit_moments()] and [local_terms()] for the logit variance
#'   \eqn{V} and the second-order approximation \eqn{H \approx I(1-I)V/2}{H = I(1 - I) V / 2}.
#'
#' @template ref-lee2026
#' @examples
#' S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#' kl_dispersion(S)
#' kl_contributions(S)
#' mean(kl_contributions(S))
#' @export
kl_dispersion <- function(S) {
  S <- .rp_check_reliabilities(S)
  p <- .rp_center_S(S)
  .rp_H_at(p$m, p$d)
}

#' @rdname kl_dispersion
#' @export
kl_contributions <- function(S) {
  S <- .rp_check_reliabilities(S)
  p <- .rp_center_S(S)
  stats::setNames(.rp_bern_kl(rep(p$m, length(S)), p$d), names(S))
}

#' Logit variance and third absolute moment
#'
#' Computes the variance \eqn{V} and the third absolute moment
#' \eqn{\rho_3}{rho3} of the logits of the site reliabilities.
#'
#' @inheritParams kl_centroid
#'
#' @details
#' With \eqn{\Delta_j = \mathrm{logit}\, S_j - \mathrm{logit}\, I}{Delta_j = logit(S_j) - logit(I)},
#' \deqn{V = \frac{1}{J}\sum_j \Delta_j^2, \qquad \rho_3 = \frac{1}{J}\sum_j |\Delta_j|^3.}{V = mean(Delta_j^2), rho3 = mean(|Delta_j|^3).}
#' Both divide by \eqn{J}. Because
#' \eqn{\mathrm{logit}\, S_j = \log \sigma_\tau^2 - \log v_j}{logit(S_j) = log(sigma_tau^2) - log(v_j)},
#' the deviations, and so \eqn{V} and \eqn{\rho_3}{rho3}, depend on the
#' sampling variances only: they do not change with the between-site
#' standard deviation.
#'
#' @return A named numeric vector with elements `V` and `rho3`.
#'
#' @template ref-lee2026
#' @examples
#' S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#' logit_moments(S)
#'
#' # The same at another sigma_tau: the logit moments do not change.
#' logit_moments(site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 1))
#' @export
logit_moments <- function(S) {
  S <- .rp_check_reliabilities(S)
  .rp_logit_moments(.rp_center_S(S)$d)
}

#' Second-order terms of the dispersion and of the centroid gap
#'
#' Compares the dispersion \eqn{H} and the gap between the mean reliability
#' and the informativeness with their second-order approximations, and
#' bounds the remainders.
#'
#' @inheritParams kl_centroid
#'
#' @details
#' Proposition 2 of Lee (2026) gives, with the logit variance \eqn{V} and
#' the third absolute logit moment \eqn{\rho_3}{rho3} of [logit_moments()],
#' \deqn{H = \tfrac12 I(1-I)V + R_H, \qquad |R_H| \le \frac{\sqrt3}{108}\rho_3,}{H = I(1 - I) V / 2 + R_H, |R_H| <= sqrt(3) rho3 / 108,}
#' \deqn{\bar S - I = \tfrac12 I(1-I)(1-2I)V + R_S, \qquad |R_S| \le \frac{\rho_3}{48}.}{mean(S_j) - I = I(1 - I)(1 - 2I) V / 2 + R_S, |R_S| <= rho3 / 48.}
#' When the leading term is larger than the bound in absolute value, the
#' exact value has the sign of the leading term. For the gap, that sign is
#' the sign of \eqn{1 - 2I}, so the mean reliability is above \eqn{I} when
#' \eqn{I < 1/2}. When the bound is the larger, the remainder can decide the
#' sign, and the gap can have the other sign; the `exact` column always
#' gives the value itself.
#'
#' @return A data frame with two rows, `H` and `gap` (the mean reliability
#'   minus \eqn{I}), and the columns
#'   * `quantity`: `"H"` or `"gap"`.
#'   * `exact`: the value of the quantity.
#'   * `leading`: its second-order term.
#'   * `bound`: the bound on the absolute remainder.
#'   * `holds`: `TRUE` if the absolute leading term is larger than the
#'     bound, `NA` if both are 0, as when all the reliabilities are equal.
#'
#' @template ref-lee2026
#' @examples
#' S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#' local_terms(S)
#'
#' # A skewed profile of 25 sites with I = 0.3 and V = 2 whose mean
#' # reliability is below I, against the sign of the leading term.
#' d <- c(sqrt(48), rep(-sqrt(1 / 12), 24))
#' local_terms(plogis(qlogis(0.3) + d))
#' @export
local_terms <- function(S) {
  S <- .rp_check_reliabilities(S)
  p <- .rp_center_S(S)
  moments <- .rp_logit_moments(p$d)
  local <- .rp_local_at(p$m, moments[["V"]], moments[["rho3"]])
  leading <- c(local$H_leading, local$gap_leading)
  bound <- c(local$H_bound, local$gap_bound)
  data.frame(
    quantity = c("H", "gap"),
    exact = c(.rp_H_at(p$m, p$d), .rp_gap_at(p$m, p$d)),
    leading = leading,
    bound = bound,
    holds = .rp_holds(leading, bound)
  )
}
