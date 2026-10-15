#' Pooling plot
#'
#' Plots each site's estimate with its interval beside the posterior mean of
#' its effect with its posterior interval, so that the partial pooling of
#' every site can be seen.
#'
#' @param x An object of class `reliability_posterior`, `reliability_draws`
#'   (with estimates and draws of the mean) or `reliability_profile` (with
#'   estimates).
#' @param order Order of the sites from bottom to top: `"estimate"` (the
#'   default), `"reliability"`, or `"input"` for the order of the input.
#' @param level Probability content of the intervals, a number between 0
#'   and 1. `NULL`, the default, uses the level of the fit (0.95 for a
#'   profile).
#' @param color_by_reliability If `TRUE`, the posterior means are colored by
#'   the site reliability instead of a single color.
#' @param max_sites Largest number of sites to label, a positive integer.
#'   With more sites the axis has no labels.
#' @param scale Factor that multiplies the estimates and the effects, a
#'   positive number; `100` turns proportions into percentage points.
#' @param ... Must be empty.
#'
#' @details
#' The estimate's interval is \eqn{\hat\tau_j \pm z\, \mathrm{se}_j}{tau_hat_j +/- z se_j},
#' with \eqn{z} the standard normal quantile at
#' \eqn{(1 + \mathrm{level})/2}{(1 + level) / 2}. The posterior intervals of a posterior fit are quantiles of the normal
#' mixture over `sigma_tau` of the effect's conditional posterior, so they
#' are not symmetric in general. For a profile at a fixed `sigma_tau`, the
#' posterior is normal, with the mean and standard deviation of
#' [site_table()]. The vertical line is the posterior mean of the effect
#' mean (for a profile, the precision-weighted mean at `sigma_tau`). The
#' horizontal axis names the scale of the effects when the fit records it.
#' With a mean model there is no single mean: each posterior mean is pulled
#' toward the site's own predicted value, and there is no vertical line.
#'
#' @return A ggplot object.
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, site = school, data = eight_schools)
#' plot_pooling(fit)
#' @export
plot_pooling <- function(x, order = c("estimate", "reliability", "input"), level = NULL,
                         color_by_reliability = FALSE, max_sites = 60, scale = 1, ...) {
  .rp_check_dots_empty(...)
  if (!inherits(x, c("reliability_posterior", "reliability_draws", "reliability_profile"))) {
    .rp_abort("`plot_pooling()` takes a posterior fit, draws with estimates, or a profile with estimates.",
              info = .rp_you_passed(x))
  }
  order <- .rp_match(order, c("estimate", "reliability", "input"), "order")
  level <- .rp_check_level(level, null_ok = TRUE)
  color_by_reliability <- .rp_check_flag(color_by_reliability, "color_by_reliability")
  max_sites <- .rp_check_max_sites(max_sites)
  scale <- .rp_check_scale(scale)
  col <- .rp_colors()
  fixed <- inherits(x, "reliability_profile")
  tab <- site_table(x, level = if (fixed) NULL else level)
  if (fixed) level <- if (is.null(level)) 0.95 else level
  if (is.null(level)) level <- attr(tab, "level")
  if (!"estimate" %in% names(tab) || !"effect" %in% names(tab)) {
    .rp_abort(
      "There are no estimates of the sites, or no posterior of their effects, to plot.",
      info = if (inherits(x, "reliability_draws")) "Give `estimate` and `mean` to `reliability_draws()`." else
        "Give `estimate` to `reliability_profile()`."
    )
  }
  z <- stats::qnorm((1 + level) / 2)
  if (fixed) {
    tab$effect_lower <- tab$effect - z * tab$effect_sd
    tab$effect_upper <- tab$effect + z * tab$effect_sd
  }
  center <- if (fixed && !is.null(x$X)) {
    numeric()
  } else if (fixed) {
    .rp_conditional_sites(x$sigma_tau, x$estimate, x$se2)$mu
  } else {
    profile_summary(x)$mean[profile_summary(x)$quantity == "effect_mean"]
  }
  ord <- switch(order,
    estimate = order(tab$estimate),
    reliability = order(tab$reliability),
    input = rev(seq_len(nrow(tab)))
  )
  tab <- tab[ord, , drop = FALSE]
  tab$row <- seq_len(nrow(tab))
  tab$est_lower <- (tab$estimate - z * tab$se) * scale
  tab$est_upper <- (tab$estimate + z * tab$se) * scale
  tab$estimate <- tab$estimate * scale
  tab$effect <- tab$effect * scale
  tab$effect_lower <- tab$effect_lower * scale
  tab$effect_upper <- tab$effect_upper * scale
  pct <- .rp_percent_label(level)
  est_label <- sprintf("estimate (%s interval)", pct)
  post_label <- sprintf("posterior mean (%s posterior interval)", pct)
  has_center <- length(center) == 1L && is.finite(center)
  p <- ggplot2::ggplot(tab, ggplot2::aes(y = .data$row))
  if (has_center) {
    p <- p + ggplot2::geom_vline(xintercept = center * scale, color = "grey60", linewidth = 0.4)
  }
  p <- p +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$est_lower, xmax = .data$est_upper),
                            color = col[["data"]], linewidth = 0.5,
                            position = ggplot2::position_nudge(y = 0.18)) +
    ggplot2::geom_point(ggplot2::aes(x = .data$estimate, shape = est_label), color = col[["data"]],
                        size = 1.7, position = ggplot2::position_nudge(y = 0.18))
  if (color_by_reliability) {
    p <- p +
      ggplot2::geom_linerange(ggplot2::aes(xmin = .data$effect_lower, xmax = .data$effect_upper,
                                           color = .data$reliability),
                              linewidth = 0.6, position = ggplot2::position_nudge(y = -0.18)) +
      ggplot2::geom_point(ggplot2::aes(x = .data$effect, shape = post_label, color = .data$reliability),
                          size = 2, position = ggplot2::position_nudge(y = -0.18)) +
      ggplot2::scale_color_viridis_c("reliability", end = 0.9, limits = c(0, 1)) +
      ggplot2::theme(legend.title = ggplot2::element_text())
  } else {
    p <- p +
      ggplot2::geom_linerange(ggplot2::aes(xmin = .data$effect_lower, xmax = .data$effect_upper),
                              color = col[["informativeness"]], linewidth = 0.6,
                              position = ggplot2::position_nudge(y = -0.18)) +
      ggplot2::geom_point(ggplot2::aes(x = .data$effect, shape = post_label),
                          color = col[["informativeness"]], size = 2,
                          position = ggplot2::position_nudge(y = -0.18))
  }
  shapes <- stats::setNames(c(16, 18), c(est_label, post_label))
  p <- p + ggplot2::scale_shape_manual(values = shapes, breaks = names(shapes)) +
    ggplot2::guides(shape = ggplot2::guide_legend(
      override.aes = list(color = c(col[["data"]], col[["informativeness"]]))
    )) +
    .rp_theme()
  labelled <- nrow(tab) <= max_sites
  if (!labelled) .rp_inform_many_sites(nrow(tab), max_sites)
  p <- p + if (labelled) {
    ggplot2::scale_y_continuous(NULL, breaks = tab$row, labels = tab$site)
  } else {
    ggplot2::scale_y_continuous(NULL, breaks = NULL)
  }
  scale_note <- if (scale != 1) sprintf("Estimates and effects are multiplied by %s.", format(scale)) else ""
  # The axis says what the values are: the scale of the effects, and the
  # factor when they are multiplied.
  unit <- c(x$labels$scale_label, if (scale != 1) paste("\u00d7", format(scale)))
  p + ggplot2::labs(
    x = if (length(unit) == 0L) "Effect" else sprintf("Effect (%s)", paste(unit, collapse = ", ")),
    caption = .rp_caption(
      sprintf("Gray: estimates with %s intervals; blue: posterior means with %s posterior intervals.", pct, pct),
      if (has_center) {
        "The vertical line is the posterior mean of the mean of the effects."
      } else if (inherits(x, "reliability_posterior") && is.null(x$data$X)) {
        paste("Site posterior means average conditional predictions over sigma_tau.",
              "The common mean has no posterior mean, so there is no vertical line.")
      } else {
        "Each posterior mean is pulled toward the site's predicted value from the mean model."
      },
      .rp_source_sentence(x), scale_note
    )
  )
}
