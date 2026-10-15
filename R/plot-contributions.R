#' Sites' terms of the dispersion
#'
#' Plots each site's term \eqn{\mathrm{KL}\{\mathrm{Ber}(I) \| \mathrm{Ber}(S_j)\}}{KL(Ber(I) || Ber(S_j))},
#' whose mean is the dispersion \eqn{H}, in order.
#'
#' @param x An object of class `reliability_profile`,
#'   `reliability_posterior`, `reliability_draws` or `reliability_design`.
#' @param max_sites Largest number of sites to label, a positive integer.
#' @param ... Must be empty.
#'
#' @details
#' The terms are those of [kl_contributions()]; for a posterior fit they are
#' posterior means.
#'
#' @return A ggplot object.
#'
#' @examples
#' prof <- reliability_profile(se = eight_schools$se, sigma_tau = 10,
#'                             site = eight_schools$school)
#' plot_contributions(prof)
#' @export
plot_contributions <- function(x, max_sites = 60, ...) {
  .rp_check_dots_empty(...)
  x <- .rp_plot_check(x, c("reliability_profile", "reliability_posterior", "reliability_draws",
                           "reliability_design"), "plot_contributions")
  max_sites <- .rp_check_max_sites(max_sites)
  tab <- site_table(x)
  ps <- profile_summary(x)
  H <- ps$mean[ps$quantity == "H"]
  tab <- tab[order(tab$kl_term), , drop = FALSE]
  tab$rank <- seq_len(nrow(tab))
  col <- .rp_colors()
  labelled <- nrow(tab) <= max_sites
  if (!labelled) .rp_inform_many_sites(nrow(tab), max_sites)
  p <- ggplot2::ggplot(tab, ggplot2::aes(.data$rank, .data$kl_term)) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$rank, y = 0, yend = .data$kl_term),
                          color = col[["dispersion"]], linewidth = 0.7) +
    ggplot2::geom_point(color = col[["dispersion"]], size = 2) +
    ggplot2::geom_hline(yintercept = H, linetype = "22", color = col[["ink"]]) +
    ggplot2::annotate("text", x = 0.6, y = H, label = sprintf("H = %s", .rp_fmt(H, 4L)), vjust = -0.5,
                      hjust = 0, size = 3.2, color = col[["ink"]]) +
    ggplot2::labs(y = "KL(Ber(I) || Ber(S_j))") +
    .rp_theme()
  p <- p + if (labelled) {
    ggplot2::scale_x_continuous(NULL, breaks = tab$rank, labels = tab$site)
  } else {
    ggplot2::scale_x_continuous("Sites, ordered by their term")
  }
  if (labelled) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                            panel.grid.major.x = ggplot2::element_blank())
  }
  p + ggplot2::labs(caption = .rp_caption(
    sprintf("Each site's term KL(Ber(I) || Ber(S_j)), the divergence of I from the site's reliability; the dashed line is their mean, the dispersion H%s.",
            if (inherits(x, "reliability_profile")) "" else " (posterior means)"),
    .rp_source_sentence(x)
  ))
}
