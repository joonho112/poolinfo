# Methods of reliability_posterior objects: printing, the summary, update(),
# and the law of sigma_tau rebuilt from the numbers that a fit stores.

# The kernel and the law of a fit. A fit made with the node rule keeps the
# rule's points and weights; the adaptive engine is rebuilt from the data.
.rp_fit_kernel <- function(x) {
  .rp_kernel(x$data$estimate, x$data$se2, x$data$X, x$prior)
}

.rp_fit_law <- function(x, kernel = .rp_fit_kernel(x)) {
  .rp_law(kernel, x$control$engine, x$control$rel_tol, x$control$abs_tol, rule = x$rule)
}

.rp_model_lines <- function(x) {
  mean_part <- if (is.null(x$data$X)) "mean" else "x'beta"
  c(sprintf("Model: estimate ~ N(effect, se^2), effect ~ N(%s, sigma_tau^2),", mean_part),
    "       sites independent")
}

.rp_conditioning_lines <- function() {
  c("Summaries are conditional on the supplied standard errors.",
    "Each site counts equally in the diagnostics, whatever its size.")
}

.rp_summary_value <- function(summary, quantity, column) {
  summary[[column]][summary$quantity == quantity]
}

# The printed names of the rows of the mean or of the coefficients.
.rp_coef_labels <- function(quantities) {
  ifelse(quantities == "effect_mean", "effect mean", paste("coefficient", sub("^coef_", "", quantities)))
}

.rp_percent <- function(level) {
  paste0(format(100 * level, trim = TRUE), "%")
}

#' @rdname reliability_posterior
#' @order 2
#' @param x,object A fit of `reliability_posterior()`.
#' @param use_cli Whether to print the headings in bold with the cli package,
#'   when it is installed. `NULL`, the default, takes
#'   `getOption("poolinfo.use_cli", TRUE)`.
#' @export
format.reliability_posterior <- function(x, ..., use_cli = NULL) {
  .rp_check_dots_empty(...)
  s <- x$summary
  J <- length(x$data$estimate)
  val <- function(q, col) .rp_summary_value(s, q, col)
  coef_rows <- s$quantity[grepl("^(effect_mean|coef_)", s$quantity)]
  coef_labels <- .rp_coef_labels(coef_rows)
  w <- max(26L, nchar(coef_labels))
  row <- function(label, q, digits, sig = FALSE, tail = "") {
    f <- if (sig) function(z) .rp_fmt_sig(z) else function(z) .rp_fmt(z, digits)
    mean_txt <- if (is.na(val(q, "mean"))) "n.d." else f(val(q, "mean"))
    sprintf("  %-*s %-8s %-8s %-8s%s", w, label, mean_txt, f(val(q, "lower")), f(val(q, "upper")), tail)
  }
  head <- c(
    sprintf("Reliability profile of %s (posterior)", .rp_sites_word(J)),
    .rp_model_lines(x),
    paste0("Prior: ", .rp_prior_text(x$prior, x$data$X)),
    .rp_label_lines(x$labels),
    .rp_conditioning_lines()
  )
  coef_lines <- vapply(seq_along(coef_rows), function(i) row(coef_labels[[i]], coef_rows[[i]], 3L, sig = TRUE), "")
  plug <- x$plugin[["I"]]
  body <- c(
    sprintf("  %-*s %-8s %s", w, "", "mean", paste(.rp_percent(x$level), "interval")),
    row("informativeness I", "I", 3L),
    row("mean reliability", "mean", 3L),
    row("dispersion H", "H", 4L),
    row("concentration count J_eff", "J_eff", 1L, tail = sprintf(" of %d", J)),
    sprintf("  %-*s %-8s %s", w, "logit variance V", .rp_fmt(val("V", "mean"), 3L),
            "(fixed by the standard errors)"),
    row("sigma_tau", "sigma_tau", 3L, sig = TRUE),
    unname(coef_lines)
  )
  extra <- c(
    sprintf("  %-*s %s", w, sprintf("Pr(I > %s)", format(x$threshold)),
            .rp_fmt(x$probability[["I_above"]], 3L)),
    sprintf("  %-*s %s", w, "I at E(sigma_tau^2)",
            if (is.na(plug)) sprintf("not defined: E(sigma_tau^2) does not exist with %d sites", J)
            else .rp_fmt(plug, 3L))
  )
  tail_line <- if (any(x$checks$status != "ok")) {
    "Some numerical checks are outside their limits; see check_numerics()."
  }
  nd_line <- if (anyNA(s$mean[s$quantity == "sigma_tau"])) {
    sprintf("n.d.: the posterior mean does not exist with %d sites under this prior.", J)
  }
  lines <- c(head, "", body, if (!is.null(nd_line)) c("", nd_line), "", extra,
             if (!is.null(tail_line)) c("", tail_line))
  .rp_style_lines(lines, bold = 1L, use_cli = use_cli)
}

#' @rdname reliability_posterior
#' @order 3
#' @export
print.reliability_posterior <- function(x, ..., use_cli = NULL) {
  cat(format(x, ..., use_cli = use_cli), sep = "\n")
  invisible(x)
}

#' @rdname reliability_posterior
#' @order 4
#' @export
summary.reliability_posterior <- function(object, ...) {
  .rp_check_dots_empty(...)
  structure(list(fit = object), class = "summary.reliability_posterior")
}

#' @rdname reliability_posterior
#' @order 5
#' @export
format.summary.reliability_posterior <- function(x, ..., use_cli = NULL) {
  .rp_check_dots_empty(...)
  fit <- x$fit
  s <- fit$summary
  q <- fit$quantiles
  J <- length(fit$data$estimate)
  pct <- function(p) paste0(format(100 * p, trim = TRUE), "%")
  lo_p <- (1 - fit$level) / 2
  hi_p <- (1 + fit$level) / 2
  coef_rows <- s$quantity[grepl("^(effect_mean|coef_)", s$quantity)]
  coef_labels <- .rp_coef_labels(coef_rows)
  w <- max(26L, nchar(coef_labels))
  table_row <- function(label, quantity, f) {
    m <- .rp_summary_value(s, quantity, "mean")
    sdv <- .rp_summary_value(s, quantity, "sd")
    med <- if (quantity %in% q$quantity) q$median[q$quantity == quantity] else NA_real_
    sprintf("  %-*s %-9s %-9s %-9s %-9s %s", w, label,
            if (is.na(m)) "n.d." else f(m), if (is.na(sdv)) "n.d." else f(sdv),
            f(.rp_summary_value(s, quantity, "lower")), f(med), f(.rp_summary_value(s, quantity, "upper")))
  }
  f3 <- function(z) .rp_fmt(z, 3L)
  f4 <- function(z) .rp_fmt(z, 4L)
  f1 <- function(z) .rp_fmt(z, 1L)
  fs <- function(z) .rp_fmt_sig(z)
  header <- sprintf("  %-*s %-9s %-9s %-9s %-9s %s", w, "", "mean", "sd", pct(lo_p), "50%", pct(hi_p))
  diag_rows <- c(
    header,
    table_row("informativeness I", "I", f3),
    table_row("mean reliability", "mean", f3),
    table_row("harmonic centroid", "harmonic", f3),
    table_row("dispersion H", "H", f4),
    table_row("concentration count J_eff", "J_eff", f1),
    sprintf("  %-*s %-9s %s", w, "logit variance V", f3(.rp_summary_value(s, "V", "mean")),
            "(fixed by the standard errors)")
  )
  hyper_rows <- c(
    header,
    table_row("sigma_tau", "sigma_tau", fs),
    table_row("sigma_tau^2", "sigma_tau2", fs),
    vapply(seq_along(coef_rows), function(i) table_row(coef_labels[[i]], coef_rows[[i]], fs), "")
  )
  plug <- fit$plugin
  post_mean <- function(name) .rp_summary_value(s, name, "mean")
  plug_rows <- if (is.na(plug[["sigma2"]])) {
    sprintf("  E(sigma_tau^2) does not exist with %d sites under this prior.", J)
  } else {
    c(
      sprintf("  %-34s %-9s %s", "", "posterior", "at E(sigma_tau^2)"),
      sprintf("  %-34s %-9s %s", "informativeness I", f3(post_mean("I")), f3(plug[["I"]])),
      sprintf("  %-34s %-9s %s", "mean reliability", f3(post_mean("mean")), f3(plug[["mean"]])),
      sprintf("  %-34s %-9s %s", "dispersion H", f4(post_mean("H")), f4(plug[["H"]])),
      sprintf("  %-34s %-9s %s", "concentration count J_eff", f1(post_mean("J_eff")), f1(plug[["J_eff"]]))
    )
  }
  hb <- fit$h_branches
  turning <- if (hb$homogeneous) {
    "The standard errors are all equal, so H is 0 at every sigma_tau."
  } else {
    c(sprintf("H is largest, at %s, where sigma_tau = %s (I = %s).",
              f4(hb$H_max), fs(hb$turning_sigma[[1L]]), f3(stats::plogis(hb$turning_logits[[1L]]))),
      "At that sigma_tau the mean reliability equals I.")
  }
  method <- if (is.null(fit$rule)) "adaptive quadrature" else sprintf("a rule of %d nodes", length(fit$rule$sigma))
  checks_lines <- c(sprintf("Integration over sigma_tau: %s.", method),
                    if (all(fit$checks$status == "ok")) "All numerical checks are within their limits." else
                      "Some numerical checks are outside their limits; see check_numerics().")
  lines <- c(
    sprintf("Reliability profile of %s (posterior)", .rp_sites_word(J)),
    .rp_model_lines(fit), paste0("Prior: ", .rp_prior_text(fit$prior, fit$data$X)), .rp_label_lines(fit$labels),
    .rp_conditioning_lines(), "",
    "Diagnostics", diag_rows, "",
    "Between-site SD and mean", hyper_rows, "",
    sprintf("Pr(I > %s) = %s", format(fit$threshold), f3(fit$probability[["I_above"]])), "",
    "Posterior means and values at the posterior mean of sigma_tau^2", plug_rows, "",
    turning, checks_lines
  )
  bold <- c(1L, which(lines %in% c("Diagnostics", "Between-site SD and mean",
                                   "Posterior means and values at the posterior mean of sigma_tau^2")))
  .rp_style_lines(lines, bold = bold, use_cli = use_cli)
}

#' @rdname reliability_posterior
#' @order 6
#' @export
print.summary.reliability_posterior <- function(x, ..., use_cli = NULL) {
  cat(format(x, ..., use_cli = use_cli), sep = "\n")
  invisible(x)
}

#' @rdname reliability_posterior
#' @order 7
#' @export
update.reliability_posterior <- function(object, ..., prior = NULL, level = NULL, threshold = NULL,
                                         control = NULL) {
  .rp_check_dots_empty(...)
  d <- object$data
  inputs <- list(estimate = d$estimate, se2 = d$se2, site = d$site, X = d$X,
                 mean_model = d$mean_model, omitted = d$omitted)
  .rp_new_posterior(
    inputs = inputs,
    prior = if (is.null(prior)) object$prior else .rp_check_prior(prior),
    level = if (is.null(level)) object$level else .rp_check_level(level),
    threshold = if (is.null(threshold)) object$threshold else .rp_check_level(threshold, "threshold"),
    labels = object$labels,
    control = if (is.null(control)) object$control else .rp_check_control(control),
    call = object$call
  )
}

#' @rdname reliability_posterior
#' @order 8
#' @param row.names,optional Ignored.
#' @export
as.data.frame.reliability_posterior <- function(x, row.names = NULL, optional = FALSE, ...) {
  site_table(x)
}
