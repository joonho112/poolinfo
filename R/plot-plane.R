#' Informativeness and logit variance plane
#'
#' Places one or several analyses by their informativeness and their logit
#' variance, with contours of the leading term of the dispersion and the half
#' of the plane where the leading term of the mean reliability minus \eqn{I}
#' is positive.
#'
#' @param x An object of class `reliability_posterior`,
#'   `reliability_draws`, `reliability_profile` or `profile_comparison`.
#' @param level Probability content of the intervals of \eqn{I}, a number
#'   between 0 and 1. `NULL`, the default, uses the level of the fits.
#' @param contours If `TRUE` (the default), the contours of
#'   \eqn{I(1 - I)V/2} are shown.
#' @param ... Must be empty.
#'
#' @details
#' The dispersion depends on the whole profile, not on \eqn{I} and \eqn{V}
#' alone, so the contours are those of its leading term
#' \eqn{\tfrac12 I(1-I)V}{I (1 - I) V / 2} (Lee, 2026, Proposition 2). Left
#' of \eqn{I = 1/2} the leading term of the mean reliability minus \eqn{I}
#' is positive; whether the exact difference is depends on the remainder,
#' which [local_terms()] bounds.
#'
#' @return A ggplot object.
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools)
#' plot_plane(fit)
#' @export
plot_plane <- function(x, level = NULL, contours = TRUE, ...) {
  .rp_check_dots_empty(...)
  level <- .rp_check_level(level, null_ok = TRUE)
  contours <- .rp_check_flag(contours, "contours")
  fits <- if (inherits(x, "profile_comparison")) x$objects else list(fit = x)
  rows <- lapply(names(fits), function(name) {
    f <- .rp_plot_check(fits[[name]], c("reliability_profile", "reliability_posterior",
                                        "reliability_draws", "reliability_design"), "plot_plane")
    ps <- profile_summary(f, level = level)
    get <- function(q, col) ps[[col]][ps$quantity == q]
    data.frame(analysis = name, I = get("I", "mean"), lower = get("I", "lower"),
               upper = get("I", "upper"), V = get("V", "mean"))
  })
  pts <- do.call(rbind, rows)
  col <- .rp_colors()
  V_top <- .rp_plane_top(pts$V)
  p <- ggplot2::ggplot(pts, ggplot2::aes(.data$I, .data$V)) +
    ggplot2::annotate("rect", xmin = 0, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "grey95") +
    ggplot2::annotate("text", x = 0.25, y = V_top, label = "leading term:\nmean reliability > I",
                      size = 2.8, color = "grey35", vjust = 1, lineheight = 0.9) +
    ggplot2::annotate("text", x = 0.75, y = V_top, label = "leading term:\nmean reliability < I",
                      size = 2.8, color = "grey35", vjust = 1, lineheight = 0.9)
  if (contours) {
    grid <- expand.grid(I = seq(0.005, 0.995, length.out = 150), V = seq(0, V_top, length.out = 150))
    grid$lead <- 0.5 * grid$I * (1 - grid$I) * grid$V
    breaks <- .rp_plane_breaks(V_top)
    p <- p + ggplot2::geom_contour(data = grid, ggplot2::aes(z = .data$lead), breaks = breaks,
                                   color = col[["dispersion"]], linewidth = 0.3, alpha = 0.7)
  }
  if (any(is.finite(pts$lower))) {
    p <- p + ggplot2::geom_errorbar(data = pts[is.finite(pts$lower), , drop = FALSE],
                                    ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                                    width = V_top * 0.02, color = "grey55", orientation = "y")
  }
  p <- p + ggplot2::geom_point(color = col[["informativeness"]], size = 2.4)
  if (length(fits) > 1L) {
    at <- .rp_place_labels(pts$I, pts$V, pts$analysis, V_top, pts$lower, pts$upper)
    pts$label_x <- at$x
    pts$label_y <- at$y
    pts$label_hjust <- at$hjust
    far <- pts[at$far, , drop = FALSE]
    if (nrow(far) > 0L) {
      # A line from the near edge of the label to its point.
      far$line_y <- far$label_y - sign(far$label_y - far$V) * 0.025 * V_top
      p <- p + ggplot2::geom_segment(data = far, ggplot2::aes(x = .data$I, y = .data$V, xend = .data$I,
                                                              yend = .data$line_y),
                                     color = "grey60", linewidth = 0.3)
    }
    p <- p + ggplot2::geom_text(data = pts, ggplot2::aes(x = .data$label_x, y = .data$label_y,
                                                         label = .data$analysis, hjust = .data$label_hjust),
                                size = 3, vjust = 0.5, color = "grey20")
  }
  lvl <- if (any(is.finite(pts$lower))) {
    sprintf(" with its %s posterior interval", .rp_percent_label(if (is.null(level)) {
      attr(profile_summary(fits[[1L]]), "level") %||% 0.95
    } else level))
  } else ""
  p + ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, V_top)) +
    ggplot2::labs(
      x = "Informativeness I", y = "Logit variance V",
      caption = .rp_caption(
        sprintf("Each point is an analysis at its informativeness I%s and its logit variance V.", lvl),
        if (contours) sprintf("Purple lines: the leading term of H, I (1 - I) V / 2, at %s.",
                              paste(format(breaks), collapse = ", ")),
        "Left of I = 1/2 the leading term of the mean reliability minus I is positive.",
        .rp_source_sentence(fits[[1L]])
      )
    ) +
    .rp_theme()
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# The top of the logit variance axis of the plane.
.rp_plane_top <- function(V) max(0.05, max(V) * 1.4)

# The levels of the contours of I (1 - I) V / 2 below the top of the plane.
.rp_plane_breaks <- function(V_top) {
  I <- seq(0.005, 0.995, length.out = 150)
  top <- 0.5 * max(I * (1 - I)) * V_top
  breaks <- pretty(c(0, top), n = 5)
  breaks[breaks > 0 & breaks < top]
}

# Positions of the labels of the analyses in the plane: next to each point,
# the first of a few positions (above, below, then to the right or left of
# the point, then farther up or down) at which the label overlaps neither a
# label already placed nor the point or interval of another analysis, and
# stays inside the plot. Labels are about 0.013 wide per character on the
# axis of I, and 0.05 of the axis of V high. `far` marks the labels placed
# more than one step from their point, which get a line to it.
.rp_place_labels <- function(I, V, labels, V_top, lower = I, upper = I) {
  step <- 0.05 * V_top
  width <- 0.013 * nchar(labels)
  gap <- 0.012
  candidates <- rbind(
    c(0, 1, 0.5), c(0, -1, 0.5), c(1, 1, 0), c(1, -1, 0), c(-1, 1, 1), c(-1, -1, 1),
    c(0, 2, 0.5), c(0, -2, 0.5), c(1, 2, 0), c(-1, 2, 1), c(0, 3, 0.5), c(0, -3, 0.5)
  )
  extent <- function(x, w, hjust) c(x - hjust * w, x + (1 - hjust) * w)
  lower[!is.finite(lower)] <- I[!is.finite(lower)]
  upper[!is.finite(upper)] <- I[!is.finite(upper)]
  out <- data.frame(x = I, y = V + step, hjust = 0.5, far = FALSE)
  placed <- integer()
  for (i in order(V, decreasing = TRUE)) {
    for (k in seq_len(nrow(candidates))) {
      cand <- candidates[k, ]
      x0 <- I[[i]] + cand[[1L]] * gap
      y0 <- V[[i]] + cand[[2L]] * step
      ext <- extent(x0, width[[i]], cand[[3L]])
      if (ext[[1L]] < -0.02 || ext[[2L]] > 1.02) next
      clash <- vapply(placed, function(j) {
        e <- extent(out$x[[j]], width[[j]], out$hjust[[j]])
        ext[[1L]] < e[[2L]] && e[[1L]] < ext[[2L]] && abs(out$y[[j]] - y0) < step
      }, NA)
      hits <- lower < ext[[2L]] + 0.01 & upper > ext[[1L]] - 0.01 & abs(V - y0) < 0.6 * step
      if (!any(clash) && !any(hits)) break
    }
    out[i, c("x", "y", "hjust")] <- c(x0, y0, cand[[3L]])
    out$far[[i]] <- abs(cand[[2L]]) > 1
    placed <- c(placed, i)
  }
  out
}

