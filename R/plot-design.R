#' Design plots
#'
#' Plots what a planned allocation implies: the cost of unequal site sizes,
#' the trade between the number of sites and the precision of the mean at a
#' fixed total, or the site size needed for a target informativeness.
#'
#' @param x An object of class `reliability_design`.
#' @param type The plot: `"imbalance"` (the default), `"sites"` or
#'   `"target"`. See Details.
#' @param ... Must be empty.
#'
#' @details
#' * `"imbalance"`: the ratio of the geometric to the arithmetic mean of the
#'   sizes against their coefficient of variation, for a lognormal population
#'   of sizes (exact) and as the local approximation
#'   \eqn{\exp(-\mathrm{CV}^2/2)}{exp(-CV^2 / 2)}, with the planned
#'   allocation marked (Lee, 2026, Appendix E.2).
#' * `"sites"`: at the planned total \eqn{N}, the informativeness of a
#'   balanced design and the standard error of the mean effect,
#'   \eqn{(\sigma_\tau^2/J + \kappa/N)^{1/2}}{sqrt(sigma_tau^2 / J + kappa / N)},
#'   against the number of sites \eqn{J} (Appendix E.3).
#' * `"target"`: the informativeness of a balanced design against the size
#'   per site, one curve for each `sigma_tau` of the design.
#'
#' The curves use the geometric mean of the design's \eqn{\kappa_j}.
#'
#' @return A ggplot object.
#'
#' @examples
#' n <- round(stats::qlnorm(stats::ppoints(40), 0, 0.47) / 40 * 2000)
#' d <- reliability_design(n, sigma_tau = c(0.15, 0.2, 0.25), kappa = 4)
#' plot_design(d)
#' plot_design(d, type = "sites")
#' plot_design(d, type = "target")
#' @export
plot_design <- function(x, type = c("imbalance", "sites", "target"), ...) {
  .rp_check_dots_empty(...)
  if (!inherits(x, "reliability_design")) {
    .rp_abort("`plot_design()` takes a `reliability_design` object.", info = .rp_you_passed(x))
  }
  type <- .rp_match(type, c("imbalance", "sites", "target"), "type")
  col <- .rp_colors()
  kappa <- exp(mean(log(x$kappa)))
  J <- length(x$n)
  N <- sum(x$n)
  sigma <- x$sigma_tau
  if (type == "imbalance") {
    cv <- seq(0, max(1.5, x$local$cv * 1.2), length.out = 200)
    d <- rbind(data.frame(cv = cv, ratio = imbalance_ratio(cv, "lognormal"), kind = "lognormal population, exact"),
               data.frame(cv = cv, ratio = imbalance_ratio(cv, "local"), kind = "local approximation"))
    here <- data.frame(cv = x$local$cv, ratio = x$efficiency)
    return(
      ggplot2::ggplot(d, ggplot2::aes(.data$cv, .data$ratio, linetype = .data$kind)) +
        ggplot2::geom_line(color = col[["ink"]], linewidth = 0.7) +
        ggplot2::geom_point(data = here, ggplot2::aes(.data$cv, .data$ratio), inherit.aes = FALSE,
                            color = col[["informativeness"]], size = 2.6) +
        ggplot2::scale_linetype_manual(values = c("lognormal population, exact" = "solid",
                                                  "local approximation" = "22")) +
        ggplot2::labs(x = "Coefficient of variation of the site sizes",
                      y = "Geometric / arithmetic mean size",
                      caption = .rp_caption(
                        sprintf("The point is the planned allocation: CV %s, ratio %s.", .rp_fmt(x$local$cv, 2L),
                                .rp_fmt(x$efficiency, 3L)),
                        "The two curves are different statements: an identity for a lognormal population of sizes, and an approximation for sizes close to their mean (Lee, 2026)."
                      )) +
        .rp_theme()
    )
  }
  if (type == "sites") {
    s <- sigma[[1L]]
    Js <- unique(round(seq(2, max(4 * J, 20), length.out = 300)))
    panels <- c("Informativeness I of a balanced design",
                "Standard error of the mean effect, relative to the planned J")
    d <- rbind(
      data.frame(J = Js, value = s^2 / (s^2 + kappa * Js / N), what = panels[[1L]]),
      data.frame(J = Js, value = sqrt(s^2 / Js + kappa / N) / sqrt(s^2 / J + kappa / N), what = panels[[2L]])
    )
    d$what <- factor(d$what, levels = panels)
    return(
      ggplot2::ggplot(d, ggplot2::aes(.data$J, .data$value, color = .data$what)) +
        ggplot2::geom_vline(xintercept = J, color = "grey70") +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::facet_wrap(~ what, ncol = 1, scales = "free_y") +
        ggplot2::scale_color_manual(values = stats::setNames(c(col[["informativeness"]], col[["plugin"]]), panels),
                                    guide = "none") +
        ggplot2::labs(x = sprintf("Number of sites J (N = %s fixed)", format(round(N), big.mark = ",")),
                      y = NULL,
                      caption = .rp_caption(
                        sprintf("At the planned total and sigma_tau = %s, more sites lower the informativeness of a balanced design while the mean effect becomes more precise (Lee, 2026).",
                                .rp_fmt_sig(s)),
                        "The gray line is the planned number of sites."
                      )) +
        .rp_theme()
    )
  }
  nbar <- exp(seq(log(1), log(max(50 * N / J, 10)), length.out = 300))
  d <- do.call(rbind, lapply(sigma, function(s) {
    data.frame(n = nbar, I = s^2 / (s^2 + kappa / nbar), sigma = sprintf("sigma_tau = %s", .rp_fmt_sig(s)))
  }))
  d$sigma <- factor(d$sigma, levels = unique(d$sigma))
  ggplot2::ggplot(d, ggplot2::aes(.data$n, .data$I, color = .data$sigma)) +
    ggplot2::geom_vline(xintercept = N / J, color = "grey70") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_x_log10("Individuals per site (balanced)") +
    ggplot2::scale_color_manual(values = .rp_series_colors(length(sigma))) +
    ggplot2::labs(y = "Informativeness I", caption = .rp_caption(
      "The informativeness of a balanced design depends on the size per site only.",
      "These are scenario calculations over plausible values of sigma_tau, not a calibration. The gray line is the planned mean size."
    )) +
    .rp_theme()
}
