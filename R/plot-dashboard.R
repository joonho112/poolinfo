#' Dashboard of a posterior fit
#'
#' Combines four plots of a posterior fit in the layout of Figure 8 of Lee
#' (2026): (a) the pooling plot, (b) the reliability profile, (c) the
#' posterior density of the informativeness and (d) the informativeness and
#' logit variance plane. `plot()` of a posterior fit shows it.
#'
#' @param x An object of class `reliability_posterior` or
#'   `reliability_draws`.
#' @param level Probability content of the intervals, a number between 0
#'   and 1. `NULL`, the default, uses the level of the fit.
#' @param max_sites Largest number of sites to label, a positive integer.
#' @param scale Factor that multiplies the estimates and the effects in
#'   panel (a), a positive number.
#' @param ... Must be empty.
#'
#' @details
#' Draws without estimates or without draws of the mean have no panel (a).
#' The captions of the panels are gathered below the plot.
#'
#' @return A patchwork object.
#'
#' @template ref-lee2026
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, site = school, data = eight_schools)
#' plot_dashboard(fit)
#' @export
plot_dashboard <- function(x, level = NULL, max_sites = 60, scale = 1, ...) {
  .rp_check_dots_empty(...)
  if (!inherits(x, c("reliability_posterior", "reliability_draws"))) {
    .rp_abort("`plot_dashboard()` takes a posterior fit or draws.", info = .rp_you_passed(x))
  }
  has_sites <- inherits(x, "reliability_posterior") ||
    (!is.null(x$estimate) && !all(is.na(x$draws$mean)))
  strip <- ggplot2::theme(plot.caption = ggplot2::element_blank())
  density_panel <- plot_functionals(x, level = level)
  density_display <- attr(density_panel, ".rp_density_display")
  panels <- list(
    if (has_sites) {
      plot_pooling(x, level = level, max_sites = max_sites, scale = scale) + strip
    },
    suppressMessages(plot_profile(x, level = level, max_sites = max_sites)) + strip,
    density_panel + strip,
    plot_plane(x, level = level) + strip
  )
  panels <- Filter(Negate(is.null), panels)
  pct <- .rp_percent_label(if (is.null(level)) x$level else level)
  mean_model <- inherits(x, "reliability_posterior") && !is.null(x$data$X)
  undefined_mean <- inherits(x, "reliability_posterior") && !mean_model &&
    !is.finite(.rp_summary_value(profile_summary(x), "effect_mean", "mean"))
  has_plugin <- is.finite(x$plugin[["I"]])
  breaks <- .rp_plane_breaks(.rp_plane_top(.rp_summary_value(profile_summary(x), "V", "mean")))
  caption <- .rp_caption(
    if (has_sites) {
      paste(sprintf("(a) Estimates (gray) and posterior means (blue) with %s intervals;", pct),
            if (mean_model) "each posterior mean is pulled toward the site's predicted value." else if (undefined_mean) {
              paste("site posterior means average conditional predictions over sigma_tau;",
                    "the common mean has no posterior mean, so there is no vertical line.")
            } else
              "the vertical line is the mean effect's posterior mean.",
            if (scale != 1) sprintf("Estimates and effects are multiplied by %s.", format(scale)))
    },
    sprintf("(%s) Posterior mean site reliabilities with %s intervals; lines: informativeness I and mean reliability.",
            if (has_sites) "b" else "a", pct),
    sprintf(paste0("(%s) Posterior density of I: solid line, posterior mean",
                   if (has_plugin) "; dashed, the value at the posterior mean of sigma_tau^2" else "",
                   if (!isTRUE(density_display$cut)) "." else if (isTRUE(density_display$exact)) {
                     ". The density grows without bound near I = 0, and its vertical axis is cut there."
                   } else ". The vertical axis is cut to show the body of the density."),
            if (has_sites) "c" else "b"),
    sprintf("(%s) Informativeness and logit variance, with contours of the leading term of H, I (1 - I) V / 2, at %s.",
            if (has_sites) "d" else "c", paste(format(breaks), collapse = ", ")),
    .rp_source_sentence(x), width = 130L
  )
  patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      tag_levels = "a", tag_prefix = "(", tag_suffix = ")", caption = caption,
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, color = "grey30", size = 9))
    )
}
