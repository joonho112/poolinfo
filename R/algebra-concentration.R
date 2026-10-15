#' Concentration count of a reliability profile
#'
#' Computes the concentration count \eqn{J_{\mathrm{eff}}}{J_eff}, the
#' inverse of the sum of the squared proportions \eqn{S_j / \sum_k S_k}{S_j / sum(S_k)}.
#' It describes how the total reliability is shared among the sites.
#'
#' @inheritParams kl_centroid
#'
#' @details
#' With \eqn{q_j = S_j / \sum_k S_k},
#' \deqn{J_{\mathrm{eff}} = \frac{1}{\sum_j q_j^2} = \frac{(\sum_j S_j)^2}{\sum_j S_j^2}}{J_eff = 1 / sum(q_j^2) = sum(S_j)^2 / sum(S_j^2)}
#' (Lee, 2026, Section 5). It lies between 1 and \eqn{J}, the number of
#' sites, and equals \eqn{J} only when all the reliabilities are equal; it
#' approaches 1 when one site has nearly all of the reliability. It
#' does not change when all the reliabilities are multiplied by a common
#' factor, so profiles that are uniformly low and uniformly high both give
#' \eqn{J}.
#'
#' Its formula is Kish's, with the reliabilities as weights, but the count
#' is not a sample size or a number of degrees of freedom, and it does not
#' measure how well the sites can be ranked. A small value means that a few
#' sites carry much of the total reliability, and [plot_concentration()]
#' shows which.
#'
#' @return The concentration count, a number between 1 and the number of
#'   sites.
#'
#' @template ref-lee2026
#' @examples
#' S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#' concentration_count(S)
#'
#' # One precise site and three noisy ones
#' concentration_count(c(0.95, 0.05, 0.05, 0.05))
#' @export
concentration_count <- function(S) {
  S <- .rp_check_reliabilities(S)
  r <- S / max(S)
  sum(r)^2 / sum(r^2)
}
