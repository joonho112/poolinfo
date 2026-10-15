#' Decomposition of the average divergence
#'
#' Plots the average divergence of a constant reliability \eqn{p} from the
#' site reliabilities, which is smallest at \eqn{p = I} with the value
#' \eqn{H}, and shows at one \eqn{p} how it splits into
#' \eqn{\mathrm{KL}\{\mathrm{Ber}(p) \| \mathrm{Ber}(I)\} + H}{KL(Ber(p) || Ber(I)) + H}.
#'
#' @param x An object of class `reliability_profile`,
#'   `reliability_posterior`, `reliability_draws` or `reliability_design`,
#'   or site reliabilities as a numeric vector.
#' @param at The constant reliability at which the split is shown, a number
#'   strictly between 0 and 1. By default a value well away from \eqn{I}.
#' @param ... Must be empty.
#'
#' @details
#' For every \eqn{p},
#' \eqn{J^{-1}\sum_j \mathrm{KL}\{\mathrm{Ber}(p)\|\mathrm{Ber}(S_j)\} = \mathrm{KL}\{\mathrm{Ber}(p)\|\mathrm{Ber}(I)\} + H}{mean(KL(Ber(p) || Ber(S_j))) = KL(Ber(p) || Ber(I)) + H}
#' (Lee, 2026, Theorem 2). The mean reliability is marked too: it minimizes
#' the average divergence with \eqn{p} in the second argument,
#' \eqn{J^{-1}\sum_j \mathrm{KL}\{\mathrm{Ber}(S_j)\|\mathrm{Ber}(p)\}}{mean(KL(Ber(S_j) || Ber(p)))},
#' not the average plotted here. For a posterior fit the profile is taken at
#' the posterior median of `sigma_tau`.
#'
#' @return A ggplot object.
#'
#' @examples
#' plot_decomposition(c(0.05, 0.10, 0.15, 0.25, 0.92))
#' @export
plot_decomposition <- function(x, at = NULL, ...) {
  .rp_check_dots_empty(...)
  prof <- .rp_teaching_profile(x, "plot_decomposition")
  S <- prof$S
  I <- kl_centroid(S)
  H <- kl_dispersion(S)
  Sbar <- mean(S)
  at <- if (is.null(at)) {
    if (I < 0.5) min(0.9, I + 0.45) else max(0.1, I - 0.45)
  } else {
    .rp_check_number(at, "at", "a single number strictly between 0 and 1", 0, 1,
                     lower_open = TRUE, upper_open = TRUE)
  }
  kl <- function(p, q) p * log(p / q) + (1 - p) * log((1 - p) / (1 - q))
  pp <- seq(0.003, 0.997, length.out = 400)
  objective <- data.frame(p = pp, value = vapply(pp, function(p) mean(kl(p, S)), 0))
  at_value <- mean(kl(at, S))
  sbar_value <- mean(kl(Sbar, S))
  col <- .rp_colors()
  top <- min(max(objective$value), max(1.2 * at_value, 3 * H, 0.2))
  # The two centroids are named in a legend below the plot, because their
  # points are close to each other and to the curve.
  marks <- data.frame(
    p = c(I, Sbar), value = c(H, sbar_value),
    label = c(sprintf("I = %s, the minimum: H = %s", .rp_fmt(I, 3L), .rp_fmt(H, 3L)),
              sprintf("mean reliability = %s", .rp_fmt(Sbar, 3L)))
  )
  marks$label <- factor(marks$label, levels = marks$label)
  ggplot2::ggplot(objective, ggplot2::aes(.data$p, .data$value)) +
    ggplot2::geom_line(color = col[["ink"]], linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = H, color = col[["dispersion"]], linetype = "22") +
    ggplot2::annotate("segment", x = at, xend = at, y = 0, yend = H, color = col[["dispersion"]],
                      linewidth = 1.4) +
    ggplot2::annotate("segment", x = at, xend = at, y = H, yend = at_value,
                      color = col[["informativeness"]], linewidth = 1.4) +
    ggplot2::geom_point(data = marks, ggplot2::aes(color = .data$label, shape = .data$label), size = 2.6) +
    ggplot2::scale_color_manual(NULL, values = stats::setNames(c(col[["informativeness"]], col[["mean"]]),
                                                               levels(marks$label))) +
    ggplot2::scale_shape_manual(NULL, values = stats::setNames(c(16, 17), levels(marks$label))) +
    ggplot2::annotate("text", x = at + 0.015, y = H / 2, label = "H", hjust = 0, size = 3.2,
                      color = col[["dispersion"]]) +
    ggplot2::annotate("text", x = at + 0.015, y = (H + at_value) / 2, label = "KL(Ber(p) || Ber(I))",
                      hjust = 0, size = 3.1, color = col[["informativeness"]]) +
    ggplot2::coord_cartesian(xlim = c(0.02, 0.98), ylim = c(0, top)) +
    ggplot2::labs(
      x = "Constant reliability p", y = "Average KL(Ber(p) || Ber(S_j))",
      caption = .rp_caption(
        sprintf("%d sites. For every p the average divergence is KL(Ber(p) || Ber(I)) plus H (Lee, 2026);",
                length(S)),
        sprintf("the bars show the split at p = %s. The mean reliability minimizes the average of KL(Ber(S_j) || Ber(p)), with p in the second argument, not this average.",
                .rp_fmt(at, 2L)),
        prof$where
      )
    ) +
    .rp_theme() +
    ggplot2::theme(legend.position = "bottom")
}
