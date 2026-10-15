#' Comparison plot
#'
#' Plots the informativeness, the mean reliability and the dispersion of
#' several analyses as points with their posterior intervals, one row per
#' analysis and one panel per diagnostic.
#'
#' @param x An object of class `profile_comparison`.
#' @param which The diagnostics: any of `"I"`, `"mean"`, `"H"` and `"J_eff"`.
#' @param ... Must be empty.
#'
#' @details
#' The numbers are those of `as.data.frame()` of the comparison. Analyses at
#' a fixed `sigma_tau` have points without intervals.
#'
#' @return A ggplot object.
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools)
#' cmp <- prior_sensitivity(fit, list(`half-normal` = sd_prior("half_normal", scale = 10)))
#' plot_compare(cmp)
#' @export
plot_compare <- function(x, which = c("I", "mean", "H"), ...) {
  .rp_check_dots_empty(...)
  if (!inherits(x, "profile_comparison")) {
    .rp_abort("`plot_compare()` takes a comparison made by `compare_profiles()`.", info = .rp_you_passed(x))
  }
  if (!is.character(which) || length(which) == 0L || !all(which %in% c("I", "mean", "H", "J_eff"))) {
    .rp_abort("`which` must hold some of \"I\", \"mean\", \"H\" and \"J_eff\".", info = .rp_you_passed(which))
  }
  t <- x$table
  labels <- c(I = "informativeness I", mean = "mean reliability", H = "dispersion H",
              J_eff = "concentration count J_eff")
  long <- do.call(rbind, lapply(which, function(w) {
    data.frame(analysis = t$analysis, quantity = labels[[w]], value = t[[w]],
               lower = if (w == "J_eff") NA_real_ else t[[paste0(w, "_lower")]],
               upper = if (w == "J_eff") NA_real_ else t[[paste0(w, "_upper")]])
  }))
  long$quantity <- factor(long$quantity, levels = labels[which])
  long$analysis <- factor(long$analysis, levels = rev(t$analysis))
  col <- .rp_colors()
  p <- ggplot2::ggplot(long, ggplot2::aes(.data$value, .data$analysis))
  if (any(is.finite(long$lower))) {
    p <- p + ggplot2::geom_errorbar(data = long[is.finite(long$lower), , drop = FALSE],
                                    ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                                    width = 0.2, color = "grey55", orientation = "y")
  }
  p + ggplot2::geom_point(color = col[["informativeness"]], size = 2.2) +
    ggplot2::facet_wrap(~ quantity, scales = "free_x", nrow = 1) +
    ggplot2::labs(x = NULL, y = NULL, caption = .rp_caption(
      "Points are posterior means (or values at a fixed sigma_tau); bars are central posterior intervals at the level of each fit.",
      "Each analysis is conditional on its own standard errors."
    )) +
    .rp_theme()
}

# The sorted reliabilities of several analyses in one plot: one line per
# analysis against the fraction of sites.
.rp_plot_compare_profiles <- function(x, ...) {
  .rp_check_dots_empty(...)
  rows <- lapply(names(x$objects), function(name) {
    tab <- site_table(x$objects[[name]])
    S <- sort(tab$reliability)
    data.frame(analysis = name, k = (seq_along(S) - 0.5) / length(S), reliability = S)
  })
  d <- do.call(rbind, rows)
  d$analysis <- factor(d$analysis, levels = names(x$objects))
  t <- x$table
  t$analysis <- factor(t$analysis, levels = names(x$objects))
  colors <- stats::setNames(.rp_series_colors(length(x$objects)), names(x$objects))
  ggplot2::ggplot(d, ggplot2::aes(.data$k, .data$reliability, color = .data$analysis)) +
    ggplot2::geom_hline(data = t, ggplot2::aes(yintercept = .data$I, color = .data$analysis),
                        linewidth = 0.5, alpha = 0.6) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_x_continuous("Fraction of sites, ordered by reliability", limits = c(0, 1)) +
    ggplot2::scale_y_continuous("Site reliability S_j", limits = c(0, 1)) +
    ggplot2::labs(caption = .rp_caption(
      "Sorted site reliabilities of each analysis (posterior means for posterior fits); the horizontal lines are the informativeness of each."
    )) +
    .rp_theme()
}
