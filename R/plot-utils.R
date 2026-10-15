# What the plot functions share.
#
# Conventions of every plot_*() function and of the plot() methods:
#
#   1. The result is a ggplot object (a patchwork object for the dashboard and
#      the two-panel plots), so that layers, scales and a theme can be added.
#      Nothing is drawn until it is printed.
#   2. There is no title. A caption says what the lines, bands and points
#      are, with the level of every interval, and where sigma_tau came from:
#      "conditional on the supplied standard errors" for a posterior, the
#      value and its source for a profile at a fixed sigma_tau.
#   3. The numbers that are drawn are those of the accessors: profile_summary(),
#      site_table(), functional_density(), profile_path(). A plot computes no
#      summary of its own.
#   4. Arguments come in this order: the object `x`; what is drawn (`type`,
#      `which`, the sites); the level; what is laid over it; `max_sites`;
#      `scale`; and `...`, which must be empty.
#   5. `scale` multiplies the estimates and the effects (and sigma_tau), so
#      that `scale = 100` turns proportions into percentage points. The
#      reliabilities and the diagnostics have no units.
#   6. Colors have fixed roles (Okabe and Ito): the informativeness is blue
#      and solid, the mean reliability vermilion and dashed, the harmonic
#      centroid gray and dotted, the dispersion and the sites' terms reddish
#      purple, the concentration bluish green, values at the posterior mean of
#      sigma_tau^2 orange and dashed, estimates gray. A difference in color is
#      repeated in the line type or the point shape.
#   7. With more than `max_sites` sites, a plot of sites draws them without
#      labels and says how to see them.
#   8. Axis and legend words are the printed words: "informativeness I",
#      "mean reliability", "dispersion H", "concentration count J_eff".
#   9. No colors of judgment. The line at I = 1/2 is drawn only with the
#      explanation that the leading term of Proposition 2 changes sign there.

.rp_colors <- function() {
  c(
    informativeness = "#0072B2", mean = "#D55E00", harmonic = "#7F7F7F",
    dispersion = "#CC79A7", concentration = "#009E73", plugin = "#E69F00",
    data = "#9E9E9E", ink = "#404040", other1 = "#CC79A7", other2 = "#009E73",
    other3 = "#56B4E9", other4 = "#E69F00"
  )
}

# Colors for several analyses in a fixed order.
.rp_series_colors <- function(n) {
  base <- unname(.rp_colors()[c("informativeness", "other1", "other2", "other3", "other4", "mean")])
  rep_len(base, n)
}

# The theme of the plots. It is added when a plot is built, so a theme that
# the caller adds afterwards replaces it.
.rp_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(hjust = 0, color = "grey30", size = ggplot2::rel(0.85)),
      plot.caption.position = "plot",
      plot.tag = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.95)),
      strip.text = ggplot2::element_text(hjust = 0, face = "bold")
    )
}

# The sentences of a caption, wrapped to the width of a plot.
.rp_caption <- function(..., width = 100L) {
  text <- paste(c(...), collapse = " ")
  paste(strwrap(text, width = width), collapse = "\n")
}

# Where sigma_tau came from, for the caption.
.rp_source_sentence <- function(x) {
  if (inherits(x, "reliability_posterior")) {
    return("Posterior summaries, conditional on the supplied standard errors.")
  }
  if (inherits(x, "reliability_draws")) {
    model <- if (is.null(x$model_label)) "" else paste0(" from ", x$model_label)
    return(sprintf("Posterior draws%s, conditional on the supplied standard errors.", model))
  }
  if (inherits(x, "reliability_profile")) {
    return(sprintf("At sigma_tau = %s (%s).", .rp_fmt_sig(x$sigma_tau), x$source))
  }
  if (inherits(x, "reliability_design")) {
    return(sprintf("Planned allocation at sigma_tau = %s.", .rp_fmt_sig(x$sigma_tau[[1L]])))
  }
  ""
}

.rp_percent_label <- function(level) {
  paste0(format(100 * level, trim = TRUE), "%")
}

# Classes that the plots of a profile take.
.rp_plot_check <- function(x, classes, fn) {
  if (!inherits(x, classes)) {
    .rp_abort(
      sprintf("`%s()` takes %s.", fn, .rp_class_list(classes)),
      info = .rp_you_passed(x)
    )
  }
  if (inherits(x, "reliability_design")) {
    return(x$planned)
  }
  x
}

.rp_class_list <- function(classes) {
  shown <- sprintf("`%s`", classes)
  if (length(shown) == 1L) {
    return(paste("a", shown))
  }
  paste("a", paste(shown[-length(shown)], collapse = ", "), "or", shown[[length(shown)]], "object")
}

.rp_check_scale <- function(scale) {
  .rp_check_number(scale, "scale", "a single positive number", lower = 0, lower_open = TRUE)
}

.rp_check_max_sites <- function(max_sites) {
  .rp_check_whole(max_sites, "max_sites")
}

# Said when a plot leaves the labels of the sites out because there are too
# many of them.
.rp_inform_many_sites <- function(J, max_sites) {
  .rp_inform(
    sprintf("There are %d sites, more than `max_sites` (%d), so they are drawn without labels.",
            J, max_sites),
    info = "Raise `max_sites` to label them.",
    subclass = "rp_message_many_sites"
  )
}

# The numbers of a profile for the plots: one row per site and the
# diagnostics, from the accessors. For a posterior or draws, the site points
# are posterior means and the diagnostics posterior means with intervals.
.rp_plot_profile_data <- function(x, level = NULL) {
  tab <- site_table(x, level = level)
  ps <- profile_summary(x, level = level)
  get <- function(q, col = "mean") ps[[col]][ps$quantity == q]
  fixed <- inherits(x, "reliability_profile")
  list(
    sites = tab,
    J = nrow(tab),
    fixed = fixed,
    I = get("I"), I_lower = if (!fixed) get("I", "lower"), I_upper = if (!fixed) get("I", "upper"),
    mean = get("mean"), harmonic = get("harmonic"), H = get("H"), V = get("V"), J_eff = get("J_eff"),
    level = if (fixed) NULL else attr(tab, "level")
  )
}

.rp_profile_annotation <- function(d) {
  sprintf("V = %s   H = %s\nmean reliability - I = %s   J_eff = %s of %d",
          .rp_fmt(d$V, 2L), .rp_fmt(d$H, 3L), .rp_fmt(d$mean - d$I, 3L), .rp_fmt(d$J_eff, 1L), d$J)
}
