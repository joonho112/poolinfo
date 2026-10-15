#' Reliability profile plot
#'
#' Plots the site reliabilities in order, with the informativeness and the
#' mean reliability as horizontal lines; for a posterior fit, with the
#' posterior interval of each reliability and of the informativeness. With
#' `type = "logit"` the profile is shown on the logit scale, where the
#' informativeness is the plain mean.
#'
#' @param x An object of class `reliability_profile`,
#'   `reliability_posterior`, `reliability_draws` or `reliability_design`.
#' @param type The plot: `"reliability"` (the default), the sorted
#'   reliabilities, or `"logit"`, their logits as a strip.
#' @param centroids The centroids shown as lines: any of `"I"` (the
#'   informativeness), `"mean"` and `"harmonic"`.
#' @param level Probability content of the posterior intervals, a number
#'   between 0 and 1. `NULL`, the default, uses the level of the fit.
#' @param intervals If `TRUE` (the default), a posterior fit shows the
#'   interval of each reliability and of the informativeness.
#' @param max_sites Largest number of sites to label, a
#'   positive integer.
#' @param ... Must be empty.
#'
#' @details
#' A corner label gives the logit variance \eqn{V}, the dispersion \eqn{H},
#' the mean reliability minus \eqn{I} and the concentration count; for a
#' posterior fit, their posterior means. The intervals of the reliabilities
#' are the reliabilities at the quantiles of `sigma_tau`, because each
#' reliability increases with it.
#'
#' On the logit scale the profile of a posterior fit is taken at the
#' posterior median of `sigma_tau`, where the informativeness is its
#' posterior median. The logits there are
#' \eqn{\mathrm{logit}\, S_j = \log \sigma_\tau^2 - \log v_j}{logit(S_j) = log(sigma_tau^2) - log(v_j)};
#' their mean is \eqn{\mathrm{logit}(I)}{logit(I)}, and the bracket shows
#' plus and minus \eqn{\sqrt V}{sqrt(V)} around it.
#'
#' @return A ggplot object.
#'
#' @examples
#' prof <- reliability_profile(se = eight_schools$se, sigma_tau = 10, site = eight_schools$school)
#' plot_profile(prof)
#' plot_profile(prof, type = "logit")
#'
#' fit <- reliability_posterior(estimate, se = se, site = school, data = eight_schools)
#' plot_profile(fit)
#' @export
plot_profile <- function(x, type = c("reliability", "logit"),
                         centroids = c("I", "mean"), level = NULL,
                         intervals = TRUE, max_sites = 60, ...) {
  .rp_check_dots_empty(...)
  x <- .rp_plot_check(x, c("reliability_profile", "reliability_posterior", "reliability_draws",
                           "reliability_design"), "plot_profile")
  type <- .rp_match(type, c("reliability", "logit"), "type")
  if (!is.character(centroids) || length(centroids) == 0L ||
      !all(centroids %in% c("I", "mean", "harmonic"))) {
    .rp_abort("`centroids` must hold some of \"I\", \"mean\" and \"harmonic\".",
              info = .rp_you_passed(centroids))
  }
  level <- .rp_check_level(level, null_ok = TRUE)
  intervals <- .rp_check_flag(intervals, "intervals")
  max_sites <- .rp_check_max_sites(max_sites)
  if (type == "logit") {
    return(.rp_plot_profile_logit(x, centroids, max_sites))
  }
  d <- .rp_plot_profile_data(x, level)
  col <- .rp_colors()
  tab <- d$sites
  tab <- tab[order(tab$reliability), , drop = FALSE]
  tab$rank <- seq_len(nrow(tab))
  lines <- data.frame(
    what = c("informativeness I", "mean reliability", "harmonic centroid"),
    key = c("I", "mean", "harmonic"),
    value = c(d$I, d$mean, d$harmonic)
  )
  lines <- lines[lines$key %in% centroids, , drop = FALSE]
  line_types <- c("informativeness I" = "solid", "mean reliability" = "22", "harmonic centroid" = "13")
  line_colors <- c("informativeness I" = col[["informativeness"]], "mean reliability" = col[["mean"]],
                   "harmonic centroid" = col[["harmonic"]])
  p <- ggplot2::ggplot(tab, ggplot2::aes(.data$rank, .data$reliability))
  show_intervals <- intervals && !d$fixed
  if (show_intervals && "I" %in% centroids) {
    p <- p + ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = d$I_lower, ymax = d$I_upper,
                               fill = col[["informativeness"]], alpha = 0.10)
  }
  p <- p +
    ggplot2::geom_hline(data = lines, ggplot2::aes(yintercept = .data$value, linetype = .data$what,
                                                   color = .data$what), linewidth = 0.7)
  if (show_intervals) {
    p <- p + ggplot2::geom_linerange(ggplot2::aes(ymin = .data$reliability_lower,
                                                  ymax = .data$reliability_upper),
                                     color = "grey70", linewidth = 0.5)
  }
  labelled <- d$J <= max_sites
  if (!labelled) .rp_inform_many_sites(d$J, max_sites)
  p <- p +
    ggplot2::geom_point(color = col[["ink"]], size = 1.6) +
    ggplot2::annotate("text", x = -Inf, y = Inf, label = .rp_profile_annotation(d), hjust = -0.03,
                      vjust = 1.2, size = 3, color = "grey20", lineheight = 1) +
    ggplot2::scale_linetype_manual(values = line_types[lines$what], breaks = lines$what) +
    ggplot2::scale_color_manual(values = line_colors[lines$what], breaks = lines$what) +
    ggplot2::scale_y_continuous("Site reliability S_j", limits = c(0, 1.05),
                                breaks = seq(0, 1, 0.25)) +
    .rp_theme()
  p <- if (labelled) {
    p + ggplot2::scale_x_continuous(NULL, breaks = tab$rank, labels = tab$site) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                     panel.grid.major.x = ggplot2::element_blank())
  } else {
    p + ggplot2::scale_x_continuous(if (d$fixed) "Sites, ordered by reliability" else
      "Sites, ordered by posterior mean reliability")
  }
  caption <- c(
    if (d$fixed) "Points are the site reliabilities." else
      "Points are the posterior means of the site reliabilities.",
    if (show_intervals) sprintf("Bars are their %s posterior intervals; the band is that of I.",
                                .rp_percent_label(d$level)),
    "Lines: informativeness I (solid), mean reliability (dashed)",
    if ("harmonic" %in% centroids) "and harmonic centroid (dotted)." else ".",
    .rp_source_sentence(x)
  )
  caption <- sub(" \\.", ".", paste(caption, collapse = " "))
  p + ggplot2::labs(caption = .rp_caption(caption))
}

# The logit scale: a strip of the logits at sigma_tau (for a posterior, at
# its median), the mean logit = logit(I) and +/- sqrt(V).
.rp_plot_profile_logit <- function(x, centroids, max_sites) {
  col <- .rp_colors()
  v <- x$se2
  sigma <- if (inherits(x, "reliability_profile")) x$sigma_tau else .rp_posterior_median_sigma(x)
  if (sigma == 0) {
    .rp_abort("At `sigma_tau = 0` every logit is -Inf, so there is nothing to show on the logit scale.")
  }
  site <- if (inherits(x, "reliability_posterior")) x$data$site else x[["site"]]
  eta <- 2 * log(sigma) - log(v)
  m <- mean(eta)
  V <- mean((eta - m)^2)
  S <- stats::plogis(eta)
  d <- data.frame(site = site, eta = eta)
  marks <- data.frame(
    what = c("informativeness I", "mean reliability", "harmonic centroid"),
    key = c("I", "mean", "harmonic"),
    eta = c(m, stats::qlogis(mean(S)), stats::qlogis(1 / mean(1 / S)))
  )
  marks <- marks[marks$key %in% centroids, , drop = FALSE]
  line_types <- c("informativeness I" = "solid", "mean reliability" = "22", "harmonic centroid" = "13")
  line_colors <- c("informativeness I" = col[["informativeness"]], "mean reliability" = col[["mean"]],
                   "harmonic centroid" = col[["harmonic"]])
  breaks_S <- pretty(stats::plogis(range(eta)), n = 6)
  breaks_S <- breaks_S[breaks_S > 0 & breaks_S < 1]
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$eta, 0)) +
    ggplot2::annotate("segment", x = m - sqrt(V), xend = m + sqrt(V), y = -0.35, yend = -0.35,
                      color = col[["informativeness"]], linewidth = 0.6) +
    ggplot2::annotate("segment", x = m + c(-1, 1) * sqrt(V), xend = m + c(-1, 1) * sqrt(V),
                      y = -0.42, yend = -0.28, color = col[["informativeness"]], linewidth = 0.6) +
    ggplot2::annotate("text", x = m, y = -0.5, label = "mean logit +/- sqrt(V)", size = 3,
                      color = col[["informativeness"]]) +
    ggplot2::geom_vline(data = marks, ggplot2::aes(xintercept = .data$eta, linetype = .data$what,
                                                   color = .data$what), linewidth = 0.7) +
    ggplot2::geom_point(color = col[["ink"]], size = 2, alpha = 0.8) +
    ggplot2::scale_linetype_manual(values = line_types[marks$what], breaks = marks$what) +
    ggplot2::scale_color_manual(values = line_colors[marks$what], breaks = marks$what) +
    ggplot2::scale_x_continuous(
      "logit of the site reliability, log(sigma_tau^2) - log(v_j)",
      sec.axis = ggplot2::sec_axis(~ stats::plogis(.), name = "Site reliability S_j",
                                   breaks = breaks_S)
    ) +
    ggplot2::scale_y_continuous(NULL, breaks = NULL, limits = c(-0.6, 0.6)) +
    .rp_theme()
  if (length(site) <= max_sites) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$site), vjust = -1, size = 2.8,
                                color = "grey30", check_overlap = TRUE)
  }
  where <- if (inherits(x, "reliability_profile")) .rp_source_sentence(x) else
    sprintf("At the posterior median of sigma_tau, %s. %s", .rp_fmt_sig(sigma), .rp_source_sentence(x))
  p + ggplot2::labs(caption = .rp_caption(
    "Points are the logits of the site reliabilities. On this scale the informativeness I is",
    "their plain mean; the bracket spans one logit standard deviation on either side.", where
  ))
}

# The posterior median of sigma_tau of a posterior fit or of draws.
.rp_posterior_median_sigma <- function(x) {
  if (inherits(x, "reliability_posterior")) {
    return(x$sigma_quantiles[[2L]])
  }
  stats::median(x$draws$sigma_tau)
}
