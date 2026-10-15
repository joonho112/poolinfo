#' Geometry of the Kullback-Leibler centroid
#'
#' Shows the geometry behind the informativeness as a centroid of the
#' reliabilities: (a) the
#' Bernoulli negative entropy, a tangent line at each of a few sites, and the
#' vertical gap between the curve and each tangent at \eqn{p = I}, which is
#' that site's term of the dispersion \eqn{H}; (b) the same gaps as bars.
#' Over all the sites, the mean of the gaps is \eqn{H}.
#'
#' @param x An object of class `reliability_profile`,
#'   `reliability_posterior`, `reliability_draws` or `reliability_design`,
#'   or site reliabilities as a numeric vector.
#' @param sites Positions or labels of the sites to show, at most 6. By
#'   default four sites spread over the range of the reliabilities.
#' @param ... Must be empty.
#'
#' @details
#' Site \eqn{j}'s term of \eqn{H},
#' \eqn{\mathrm{KL}\{\mathrm{Ber}(I)\|\mathrm{Ber}(S_j)\}}{KL(Ber(I) || Ber(S_j))}, is the
#' Bregman divergence of
#' \eqn{\phi(p) = p \log p + (1 - p)\log(1 - p)}{phi(p) = p log(p) + (1 - p) log(1 - p)}:
#' the height of \eqn{\phi} at \eqn{I} above its tangent at \eqn{S_j}
#' (Lee, 2026, Section 4). For a posterior fit the profile is taken at the
#' posterior median of `sigma_tau`. Panel (b) shows these sites
#' only; their mean is not \eqn{H} unless all the sites are shown.
#'
#' @return A patchwork object with two panels.
#'
#' @examples
#' plot_geometry(c(0.2, 0.4, 0.6, 0.85))
#' @export
plot_geometry <- function(x, sites = NULL, ...) {
  .rp_check_dots_empty(...)
  prof <- .rp_teaching_profile(x, "plot_geometry")
  S <- prof$S
  J <- length(S)
  chosen <- .rp_choose_sites(sites, prof$site, S)
  I <- kl_centroid(S)
  H <- kl_dispersion(S)
  col <- .rp_colors()
  palette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#7F7F7F")
  phi <- function(p) p * log(p) + (1 - p) * log1p(-p)
  pp <- seq(0.002, 0.998, length.out = 400)
  Sc <- S[chosen]
  lab <- prof$site[chosen]
  tang <- do.call(rbind, lapply(seq_along(Sc), function(j) {
    xs <- seq(max(0.005, min(Sc[j], I) - 0.12), min(0.995, max(Sc[j], I) + 0.12), length.out = 60)
    data.frame(site = lab[[j]], p = xs, y = phi(Sc[j]) + stats::qlogis(Sc[j]) * (xs - Sc[j]))
  }))
  offset <- (seq_along(Sc) - (length(Sc) + 1) / 2) * 0.012
  gaps <- data.frame(site = lab, x = I + offset, top = phi(I),
                     bottom = phi(Sc) + stats::qlogis(Sc) * (I - Sc), S = Sc,
                     kl = kl_contributions(S)[chosen])
  site_colors <- stats::setNames(palette[seq_along(Sc)], lab)
  y_low <- min(c(phi(pp), gaps$bottom)) - 0.05
  pa <- ggplot2::ggplot() +
    ggplot2::geom_line(data = data.frame(p = pp, y = phi(pp)), ggplot2::aes(.data$p, .data$y),
                       color = col[["ink"]], linewidth = 0.8) +
    ggplot2::geom_line(data = tang, ggplot2::aes(.data$p, .data$y, color = .data$site), linewidth = 0.5) +
    ggplot2::geom_point(data = gaps, ggplot2::aes(.data$S, phi(.data$S), color = .data$site), size = 2.2) +
    ggplot2::geom_segment(data = gaps, ggplot2::aes(x = .data$x, xend = .data$x, y = .data$bottom,
                                                    yend = .data$top, color = .data$site), linewidth = 1.2) +
    ggplot2::geom_vline(xintercept = I, color = col[["informativeness"]], linewidth = 0.4, linetype = "22") +
    ggplot2::annotate("text", x = I, y = max(phi(pp)), label = "I", color = col[["informativeness"]],
                      hjust = -0.6, vjust = 1, fontface = "bold") +
    ggplot2::scale_color_manual(values = site_colors, guide = "none") +
    ggplot2::coord_cartesian(ylim = c(max(y_low, -1.2), -0.05), xlim = c(0.02, 0.98)) +
    ggplot2::labs(x = "p", y = "p log(p) + (1 - p) log(1 - p)") +
    .rp_theme()
  gaps$site <- factor(gaps$site, levels = lab)
  pb <- ggplot2::ggplot(gaps, ggplot2::aes(.data$site, .data$kl, fill = .data$site)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("S = %.2f", .data$S)), vjust = -0.4, size = 2.9) +
    ggplot2::scale_fill_manual(values = site_colors, guide = "none") +
    ggplot2::labs(x = "Site", y = "KL(Ber(I) || Ber(S_j))") +
    .rp_theme()
  if (length(chosen) == J) {
    pb <- pb + ggplot2::geom_hline(yintercept = H, linetype = "22", color = col[["ink"]]) +
      ggplot2::annotate("text", x = 0.55, y = H, label = sprintf("H = %s", .rp_fmt(H, 3L)),
                        vjust = -0.5, hjust = 0, size = 3.2)
  }
  patchwork::wrap_plots(pa, pb, nrow = 1) + patchwork::plot_layout(widths = c(1.6, 1)) +
    patchwork::plot_annotation(
      tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
      caption = .rp_caption(
        sprintf("%s with I = %s.", if (length(chosen) == J) sprintf("%d sites", J) else
          sprintf("%d of the %d sites", length(chosen), J), .rp_fmt(I, 3L)),
        "(a) Each colored line is a site's tangent to the curve; the vertical gap between the curve and that tangent at p = I is the site's divergence (the gaps are offset a little so that they do not overlap).",
        sprintf("(b) The gaps%s.", if (length(chosen) == J) "; their mean is H" else
          sprintf(" of the drawn sites; the mean over all sites is H = %s", .rp_fmt(H, 3L))),
        prof$where, width = 120L
      ),
      theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, color = "grey30", size = 9))
    )
}

# The reliabilities of what a teaching plot is given: a vector of
# reliabilities, a profile, or a posterior fit at its median sigma_tau.
.rp_teaching_profile <- function(x, fn) {
  if (is.numeric(x) && !is.object(x)) {
    S <- .rp_check_reliabilities(x, "x")
    site <- if (!is.null(names(S))) names(S) else as.character(seq_along(S))
    return(list(S = unname(S), site = site, where = ""))
  }
  x <- .rp_plot_check(x, c("reliability_profile", "reliability_posterior", "reliability_draws",
                           "reliability_design"), fn)
  if (inherits(x, "reliability_profile")) {
    if (x$sigma_tau == 0) {
      .rp_abort("At `sigma_tau = 0` every reliability is 0, so there is no centroid to show.")
    }
    return(list(S = x$S, site = x$site, where = .rp_source_sentence(x)))
  }
  sigma <- .rp_posterior_median_sigma(x)
  site <- if (inherits(x, "reliability_posterior")) x$data$site else x$site
  list(S = sigma^2 / (sigma^2 + x$se2), site = site,
       where = sprintf("At the posterior median of sigma_tau, %s.", .rp_fmt_sig(sigma)))
}

# Up to 6 sites to draw: those asked for, or 4 spread over the reliabilities.
.rp_choose_sites <- function(sites, labels, S) {
  J <- length(S)
  if (is.null(sites)) {
    o <- order(S)
    at <- unique(round(seq(1, J, length.out = min(4L, J))))
    return(sort(o[at]))
  }
  at <- if (is.character(sites)) match(sites, labels) else sites
  if (anyNA(at) || !is.numeric(at) || any(at < 1 | at > J | at != floor(at)) || length(at) > 6L) {
    .rp_abort("`sites` must name or number at most 6 of the sites.", info = .rp_you_passed(sites))
  }
  unique(as.integer(at))
}
