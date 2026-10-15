#' Dispersion and centroid gap along the informativeness
#'
#' Plots the dispersion \eqn{H} and the mean reliability minus \eqn{I} as the
#' informativeness moves with the centered logits of the profile held fixed,
#' as they do when only `sigma_tau` changes.
#'
#' @param x An object of class `reliability_profile`,
#'   `reliability_posterior`, `reliability_draws` or `reliability_design`,
#'   or a `profile_comparison` of them.
#' @param leading If `TRUE` (the default), dashed gray lines show the
#'   second-order terms of Proposition 2.
#' @param ... Must be empty.
#'
#' @details
#' Along the path \eqn{dH/dm = \bar S - I}{dH / dm = mean(S_j) - I}, where
#' \eqn{m} is the logit of \eqn{I} (Lee, 2026, Proposition 2).
#' For unequal sampling variances, \eqn{H} is strictly concave in \eqn{I}
#' along the path and is largest at the single point where the gap crosses
#' 0. With equal sampling variances, \eqn{H} and the gap are 0 throughout
#' the path and there is no turning point. The point on
#' each curve marks the profile's own informativeness (for a posterior fit,
#' its posterior median), and the cross marks the turning point.
#'
#' @return A patchwork object with two panels.
#'
#' @examples
#' prof <- reliability_profile(se = eight_schools$se, sigma_tau = 10)
#' plot_gap_path(prof)
#' @export
plot_gap_path <- function(x, leading = TRUE, ...) {
  .rp_check_dots_empty(...)
  leading <- .rp_check_flag(leading, "leading")
  fits <- if (inherits(x, "profile_comparison")) x$objects else list(profile = x)
  col <- .rp_colors()
  paths <- list()
  here <- list()
  for (name in names(fits)) {
    f <- .rp_plot_check(fits[[name]], c("reliability_profile", "reliability_posterior",
                                        "reliability_draws", "reliability_design"), "plot_gap_path")
    path <- profile_path(f, n = 399)
    path$analysis <- name
    paths[[name]] <- path
    sigma <- if (inherits(f, "reliability_profile")) f$sigma_tau else .rp_posterior_median_sigma(f)
    centered <- .rp_center_v(f$se2)
    m_here <- 2 * log(sigma) - centered$log_gm
    turning <- .rp_turning_points(centered$d)
    here[[name]] <- data.frame(
      analysis = name, I = stats::plogis(m_here), H = .rp_H_at(m_here, centered$d),
      gap = .rp_gap_at(m_here, centered$d),
      I_turn = if (turning$homogeneous) NA_real_ else stats::plogis(turning$roots),
      H_turn = if (turning$homogeneous) NA_real_ else .rp_H_at(turning$roots, centered$d)
    )
  }
  path <- do.call(rbind, paths)
  pts <- do.call(rbind, here)
  path$analysis <- factor(path$analysis, levels = names(fits))
  pts$analysis <- factor(pts$analysis, levels = names(fits))
  colors <- stats::setNames(.rp_series_colors(length(fits)), names(fits))
  base <- function(y, lead, ylab) {
    p <- ggplot2::ggplot(path, ggplot2::aes(.data$I, .data[[y]], color = .data$analysis))
    if (y == "gap") p <- p + ggplot2::geom_hline(yintercept = 0, color = "grey60")
    if (leading) {
      p <- p + ggplot2::geom_line(ggplot2::aes(y = .data[[lead]], group = .data$analysis),
                                  color = "grey55", linetype = "22")
    }
    p + ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::labs(x = "Informativeness I", y = ylab) +
      .rp_theme()
  }
  pa <- base("H", "H_leading", "Dispersion H") +
    ggplot2::geom_point(data = pts, ggplot2::aes(.data$I, .data$H), size = 2.2) +
    ggplot2::geom_point(data = pts[is.finite(pts$I_turn), , drop = FALSE],
                        ggplot2::aes(.data$I_turn, .data$H_turn), shape = 4, size = 2.6, stroke = 1)
  pb <- base("gap", "gap_leading", "Mean reliability - I") +
    ggplot2::geom_point(data = pts, ggplot2::aes(.data$I, .data$gap), size = 2.2) +
    ggplot2::geom_point(data = pts[is.finite(pts$I_turn), , drop = FALSE],
                        ggplot2::aes(.data$I_turn, 0), shape = 4, size = 2.6, stroke = 1)
  if (length(fits) == 1L) {
    pa <- pa + ggplot2::guides(color = "none")
    pb <- pb + ggplot2::guides(color = "none")
  }
  where <- if (inherits(fits[[1L]], "reliability_profile")) "the profile's own I" else
    "I at the posterior median of sigma_tau"
  patchwork::wrap_plots(pa, pb, nrow = 1) + patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      caption = .rp_caption(
        "The centered logits of the profile are held fixed while I moves, as they are when only sigma_tau changes.",
        if (leading) "Dashed gray: the second-order term of each curve (Lee, 2026).",
        sprintf("Points: %s; crosses, when present: the turning point of H, where the mean reliability equals I.", where),
        if (any(vapply(fits, function(f) all(f$se2 == f$se2[[1L]]), TRUE))) {
          "With equal sampling variances, H and the gap are 0 throughout the path and there is no turning point."
        },
        width = 120L
      ),
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, color = "grey30", size = 9))
    ) &
    ggplot2::theme(legend.position = "bottom")
}
