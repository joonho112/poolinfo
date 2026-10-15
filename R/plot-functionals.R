#' Plot of the posterior density of a diagnostic
#'
#' Plots the exact posterior density of the informativeness, or of another
#' diagnostic, with its posterior mean, its value at the posterior mean of
#' `sigma_tau^2` and its central posterior interval.
#'
#' @param x An object of class `reliability_posterior` or
#'   `reliability_draws`, or a `profile_comparison` of them.
#' @param which The diagnostic: one of `"I"` (the default), `"mean"`, `"H"`,
#'   `"J_eff"` or `"sigma_tau"`.
#' @param level Probability content of the shaded interval, a number between
#'   0 and 1. `NULL`, the default, uses the level of the fit.
#' @param plugin If `TRUE` (the default), a dashed line marks the value at
#'   the posterior mean of `sigma_tau^2`, when that mean exists.
#' @param scale Factor that multiplies `sigma_tau` for `which = "sigma_tau"`,
#'   a positive number.
#' @param ... Must be empty.
#'
#' @details
#' For a posterior fit the density comes from [functional_density()].
#' The density of the informativeness grows without bound as \eqn{I}
#' approaches 0 under every prior offered by [sd_prior()], whereas the
#' posterior density of `sigma_tau` has a finite positive limit at 0.
#' For draws it is a kernel density estimate of the draws.
#' For the informativeness, a label gives the posterior probability that
#' \eqn{I} exceeds the threshold of the fit; the gray line is at one half,
#' where the second-order term of the mean reliability minus \eqn{I}
#' changes sign. See [prob_above()].
#'
#' @return A ggplot object.
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools)
#' plot_functionals(fit)
#' plot_functionals(fit, which = "H")
#' @export
plot_functionals <- function(x, which = c("I", "mean", "H", "J_eff", "sigma_tau"), level = NULL,
                             plugin = TRUE, scale = 1, ...) {
  .rp_check_dots_empty(...)
  which <- .rp_match(which, c("I", "mean", "H", "J_eff", "sigma_tau"), "which")
  level <- .rp_check_level(level, null_ok = TRUE)
  plugin <- .rp_check_flag(plugin, "plugin")
  scale <- .rp_check_scale(scale)
  fits <- if (inherits(x, "profile_comparison")) x$objects else list(x)
  if (is.null(names(fits))) names(fits) <- "fit"
  for (f in fits) {
    if (!inherits(f, c("reliability_posterior", "reliability_draws"))) {
      .rp_abort("`plot_functionals()` takes posterior fits or draws.", info = .rp_you_passed(f))
    }
  }
  pieces <- lapply(names(fits), function(name) .rp_functional_pieces(fits[[name]], which, level, name, scale))
  dens <- do.call(rbind, lapply(pieces, `[[`, "density"))
  marks <- do.call(rbind, lapply(pieces, `[[`, "marks"))
  dens$analysis <- factor(dens$analysis, levels = names(fits))
  marks$analysis <- factor(marks$analysis, levels = names(fits))
  colors <- stats::setNames(.rp_series_colors(length(fits)), names(fits))
  shade <- dens[dens$inside, , drop = FALSE]
  # A transformed exact density can grow without bound at the value for
  # sigma_tau = 0. The vertical axis stops at 1.3 times the largest
  # density beyond the lowest 2% of each grid, and the caption says so.
  bulk <- unlist(lapply(split(dens, dens$analysis), function(d) {
    d$density[d$value >= stats::quantile(d$value, 0.02)]
  }))
  top_all <- max(dens$density[is.finite(dens$density)])
  top <- min(top_all, 1.3 * max(bulk[is.finite(bulk)])) * 1.05
  cut <- top_all > top
  p <- ggplot2::ggplot(dens, ggplot2::aes(.data$value, .data$density, color = .data$analysis))
  if (which == "I") {
    p <- p + ggplot2::geom_vline(xintercept = 0.5, color = "grey85")
  }
  p <- p +
    ggplot2::geom_area(data = shade, ggplot2::aes(fill = .data$analysis), alpha = 0.15, color = NA,
                       position = "identity") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_segment(data = marks, ggplot2::aes(x = .data$mean, xend = .data$mean, y = 0, yend = top,
                                                     color = .data$analysis), linewidth = 0.6)
  has_plugin <- plugin && any(is.finite(marks$plugin))
  if (has_plugin) {
    p <- p + ggplot2::geom_segment(data = marks[is.finite(marks$plugin), , drop = FALSE],
                                   ggplot2::aes(x = .data$plugin, xend = .data$plugin, y = 0, yend = top,
                                                color = .data$analysis),
                                   linetype = "22", linewidth = 0.6)
  }
  if (which == "I") {
    # The label goes on the side away from most of the mass.
    right <- stats::median(marks$mean) < 0.5
    marks$label_y <- top * (0.95 - 0.08 * (seq_len(nrow(marks)) - 1))
    p <- p + ggplot2::geom_text(data = marks, ggplot2::aes(x = if (right) 0.98 else 0.02, y = .data$label_y,
                                                           label = .data$text, color = .data$analysis),
                                hjust = if (right) 1 else 0, size = 3.2, show.legend = FALSE)
  }
  axis <- c(I = "Informativeness I", mean = "Mean reliability", H = "Dispersion H",
            J_eff = "Concentration count J_eff", sigma_tau = "Between-site SD sigma_tau")[[which]]
  p <- p +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = axis, y = "Posterior density") +
    .rp_theme()
  if (length(fits) == 1L) {
    p <- p + ggplot2::guides(color = "none", fill = "none")
  }
  p <- p + ggplot2::coord_cartesian(xlim = if (which %in% c("I", "mean")) c(0, 1), ylim = c(0, top))
  exact <- all(vapply(fits, inherits, TRUE, "reliability_posterior"))
  divergent <- exact && which != "sigma_tau"
  p <- p + ggplot2::labs(caption = .rp_caption(
    if (exact) "Exact posterior density (a change of variables from the posterior of sigma_tau)." else
      "Kernel density estimate of the draws.",
    "Solid line: posterior mean;",
    if (has_plugin) "dashed line: the value at the posterior mean of sigma_tau^2;",
    sprintf("shaded: %s posterior interval.", .rp_percent_label(pieces[[1L]]$level)),
    if (cut && divergent) {
      if (which == "I") "The density grows without bound near I = 0; its vertical axis is cut there." else
        "The density grows without bound near the value at sigma_tau = 0, where the vertical axis is cut."
    } else if (cut) "The vertical axis is cut to show the body of the density.",
    if (which == "I") "The gray line at 1/2 is where the second-order term of the mean reliability minus I changes sign (Lee, 2026).",
    .rp_source_sentence(fits[[1L]])
  ))
  attr(p, ".rp_density_display") <- list(cut = cut, exact = exact, which = which)
  p
}

# The density, the shading and the marks of one fit.
.rp_functional_pieces <- function(x, which, level, name, scale) {
  ps <- profile_summary(x, level = level)
  level <- attr(ps, "level")
  row <- if (which == "sigma_tau") "sigma_tau" else which
  m <- ps$mean[ps$quantity == row]
  lo <- ps$lower[ps$quantity == row]
  hi <- ps$upper[ps$quantity == row]
  if (inherits(x, "reliability_posterior")) {
    d <- functional_density(x, which, n = 512)
  } else {
    k <- .rp_draws_kernel(x$se2)
    f <- .rp_functionals_at(k, x$draws$sigma_tau)
    col <- if (which == "sigma_tau") "sigma" else which
    kd <- stats::density(f[, col], n = 512, from = min(f[, col]), to = max(f[, col]))
    d <- data.frame(value = kd$x, density = kd$y)
  }
  plug <- if (which == "sigma_tau") sqrt(x$plugin[["sigma2"]]) else x$plugin[[which]]
  f <- if (which == "sigma_tau") scale else 1
  d$value <- d$value * f
  d$density <- d$density / f
  d$analysis <- name
  d$inside <- d$value >= lo * f & d$value <= hi * f
  prob <- if (which == "I") prob_above(x) else NA_real_
  marks <- data.frame(analysis = name, mean = m * f, plugin = plug * f,
                      text = if (which == "I") sprintf("%sPr(I > %s) = %s",
                                                       if (name == "fit") "" else paste0(name, ": "),
                                                       format(x$threshold), .rp_fmt(prob, 2L)) else "")
  list(density = d, marks = marks, level = level)
}
