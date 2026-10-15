#' Concentration of the reliability
#'
#' Plots the cumulative proportion of the total reliability against the
#' fraction of sites, the most reliable sites first, with the diagonal of
#' equal proportions.
#'
#' @param x An object of class `reliability_profile`,
#'   `reliability_posterior`, `reliability_draws` or `reliability_design`,
#'   or a `profile_comparison` of them.
#' @param label_top Number of the most reliable sites to name, an integer of
#'   at least 0.
#' @param ... Must be empty.
#'
#' @details
#' The proportions are \eqn{q_j = S_j / \sum_k S_k}{q_j = S_j / sum(S_k)}, and
#' for a fixed profile the inverse of the sum of their squares is the
#' concentration count \eqn{J_{\mathrm{eff}}}{J_eff}. For posterior fits and
#' draws, each step adds the posterior mean of a site's proportion, and the
#' label gives the posterior mean of \eqn{J_{\mathrm{eff}}}{J_eff}, which in
#' general differs from the inverse of the sum of the squared posterior
#' mean proportions. A curve on the diagonal means that every site has the
#' same proportion.
#'
#' @return A ggplot object.
#'
#' @examples
#' prof <- reliability_profile(se = eight_schools$se, sigma_tau = 10,
#'                             site = eight_schools$school)
#' plot_concentration(prof)
#' @export
plot_concentration <- function(x, label_top = 3, ...) {
  .rp_check_dots_empty(...)
  label_top <- .rp_check_whole(label_top, "label_top", lower = 0)
  fits <- if (inherits(x, "profile_comparison")) x$objects else list(profile = x)
  curves <- list()
  labels <- list()
  jeffs <- list()
  for (name in names(fits)) {
    f <- .rp_plot_check(fits[[name]], c("reliability_profile", "reliability_posterior",
                                        "reliability_draws", "reliability_design"), "plot_concentration")
    tab <- site_table(f)
    o <- order(tab$share, decreasing = TRUE)
    q <- tab$share[o]
    J <- length(q)
    # The proportions add up to 1 up to rounding, which can leave the last
    # step a little above the axis.
    cum <- pmin(cumsum(q), 1)
    curves[[name]] <- data.frame(analysis = name, k = c(0, seq_len(J)) / J, share = c(0, cum))
    ps <- profile_summary(f)
    top <- utils::head(o, label_top)
    labels[[name]] <- data.frame(analysis = rep(name, length(top)), k = seq_along(top) / J,
                                 share = cum[seq_along(top)], site = tab$site[top])
    jeffs[[name]] <- data.frame(analysis = name, jeff = ps$mean[ps$quantity == "J_eff"], J = J)
  }
  curve <- do.call(rbind, curves)
  lab <- do.call(rbind, labels)
  curve$analysis <- factor(curve$analysis, levels = names(fits))
  colors <- if (length(fits) == 1L) stats::setNames(.rp_colors()[["concentration"]], names(fits)) else
    stats::setNames(.rp_series_colors(length(fits)), names(fits))
  jeff <- do.call(rbind, jeffs)
  jeff$text <- sprintf("%sJ_eff = %s of %d", if (length(fits) > 1L) paste0(jeff$analysis, ": ") else "",
                       .rp_fmt(jeff$jeff, 1L), jeff$J)
  jeff$y <- 0.3 - 0.08 * (seq_len(nrow(jeff)) - 1)
  p <- ggplot2::ggplot(curve, ggplot2::aes(.data$k, .data$share, color = .data$analysis)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey75") +
    ggplot2::geom_step(linewidth = 0.8, direction = "hv") +
    ggplot2::geom_text(data = jeff, ggplot2::aes(x = 0.45, y = .data$y, label = .data$text),
                       hjust = 0, size = 3.1, show.legend = FALSE)
  if (label_top > 0L && nrow(lab) > 0L) {
    p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(label = .data$site), hjust = -0.2, vjust = 1.2,
                                size = 2.8, show.legend = FALSE)
  }
  p <- p +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_x_continuous("Fraction of sites, most reliable first", limits = c(0, 1)) +
    ggplot2::scale_y_continuous("Cumulative proportion of the reliability", limits = c(0, 1)) +
    .rp_theme()
  if (length(fits) == 1L) p <- p + ggplot2::guides(color = "none")
  posterior <- vapply(fits, inherits, TRUE, c("reliability_posterior", "reliability_draws"))
  fixed_caption <- paste(
    "Each step adds the proportion S_j / sum(S_k) of the next most reliable site.",
    "The concentration count J_eff is the inverse of the sum of the squared proportions."
  )
  posterior_caption <- paste(
    "Each step adds the posterior mean of the proportion S_j / sum(S_k) of the next most reliable site.",
    "The label gives the posterior mean of J_eff, which in general differs from the inverse of the sum of the squared posterior mean proportions."
  )
  p + ggplot2::labs(caption = .rp_caption(
    if (all(posterior)) posterior_caption else if (!any(posterior)) fixed_caption else paste(
      "For fixed profiles and designs:", fixed_caption,
      "For posterior fits and draws:", posterior_caption
    ),
    "On the diagonal every site would have the same proportion.",
    if (label_top > 0L && nrow(lab) > 0L) {
      sprintf("The labels name the %s most reliable %s.", if (label_top == 1L) "" else label_top,
              if (label_top == 1L) "site" else "sites")
    },
    .rp_source_sentence(fits[[1L]])
  ))
}
