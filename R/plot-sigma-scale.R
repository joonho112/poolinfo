#' Posterior of sigma_tau on the scale of the informativeness
#'
#' Plots the posterior density of the between-site standard deviation with
#' an axis of the informativeness on top, and shades the region where the
#' informativeness exceeds the threshold of the fit.
#'
#' @param x An object of class `reliability_posterior` or
#'   `reliability_draws`.
#' @param prior If `TRUE` (the default), a proper prior is shown too, scaled
#'   to the height of the posterior.
#' @param scale Factor that multiplies `sigma_tau`, a positive number.
#' @param ... Must be empty.
#'
#' @details
#' The informativeness \eqn{I = \sigma_\tau^2/\{\sigma_\tau^2 + \mathrm{GM}(v)\}}{I = sigma_tau^2 / (sigma_tau^2 + GM(v))}
#' increases with `sigma_tau`, so its posterior is the image of that of
#' `sigma_tau`; \eqn{I > c} exactly when
#' \eqn{\sigma_\tau > \{\mathrm{GM}(v)\, c/(1-c)\}^{1/2}}{sigma_tau > sqrt(GM(v) c / (1 - c))},
#' which is \eqn{\mathrm{GM}(v)^{1/2}}{sqrt(GM(v))} for \eqn{c = 1/2}. For
#' draws the density is a kernel estimate.
#'
#' @return A ggplot object.
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools)
#' plot_sigma_scale(fit)
#' @export
plot_sigma_scale <- function(x, prior = TRUE, scale = 1, ...) {
  .rp_check_dots_empty(...)
  if (!inherits(x, c("reliability_posterior", "reliability_draws"))) {
    .rp_abort("`plot_sigma_scale()` takes a posterior fit or draws.", info = .rp_you_passed(x))
  }
  prior <- .rp_check_flag(prior, "prior")
  scale <- .rp_check_scale(scale)
  col <- .rp_colors()
  gm <- exp(mean(log(x$se2)))
  threshold <- x$threshold
  sigma_c <- sqrt(gm * threshold / (1 - threshold))
  if (inherits(x, "reliability_posterior")) {
    d <- functional_density(x, "sigma_tau", n = 512)
  } else {
    kd <- stats::density(x$draws$sigma_tau, n = 512, from = 0)
    d <- data.frame(value = kd$x, density = kd$y)
  }
  d$value <- d$value * scale
  d$density <- d$density / scale
  shade <- d[d$value >= sigma_c * scale, , drop = FALSE]
  p_above <- prob_above(x)
  breaks_I <- c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$value, .data$density)) +
    ggplot2::geom_area(data = shade, fill = col[["informativeness"]], alpha = 0.15) +
    ggplot2::geom_line(color = col[["informativeness"]], linewidth = 0.8)
  pr <- if (inherits(x, "reliability_posterior")) x$prior else NULL
  if (prior && !is.null(pr) && pr$type != "flat") {
    lp <- .rp_prior_log_density(pr, d$value / scale)
    dens_p <- exp(lp - max(lp[is.finite(lp)]))
    dens_p[!is.finite(dens_p)] <- 0
    d$prior <- dens_p / max(dens_p) * max(d$density)
    p <- p + ggplot2::geom_line(data = d, ggplot2::aes(y = .data$prior), color = col[["harmonic"]],
                                linetype = "13", linewidth = 0.6)
  }
  p +
    ggplot2::geom_vline(xintercept = sigma_c * scale, color = "grey60", linetype = "22") +
    ggplot2::annotate("text", x = sigma_c * scale, y = max(d$density) * 0.95, hjust = -0.05, size = 3.1,
                      label = sprintf("I > %s when sigma_tau > %s\nPr(I > %s) = %s", format(threshold),
                                      .rp_fmt_sig(sigma_c * scale), format(threshold),
                                      .rp_fmt(p_above, 3L))) +
    ggplot2::scale_x_continuous(
      "Between-site SD sigma_tau",
      sec.axis = ggplot2::sec_axis(~ (. / scale)^2 / ((. / scale)^2 + gm), name = "Informativeness I",
                                   breaks = breaks_I)
    ) +
    ggplot2::labs(y = "Posterior density", caption = .rp_caption(
      "I increases with sigma_tau, so its posterior is the image of the posterior of sigma_tau (top axis).",
      sprintf("Shaded: I > %s.", format(threshold)),
      if (prior && !is.null(pr) && pr$type != "flat") "Dotted gray: the prior, scaled to the height of the posterior.",
      if (inherits(x, "reliability_draws")) "The density is a kernel estimate of the draws.",
      .rp_source_sentence(x),
      if (scale != 1) sprintf("sigma_tau is multiplied by %s.", format(scale))
    )) +
    .rp_theme()
}
