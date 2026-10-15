#' Diagnostics of several analyses side by side
#'
#' Collects several posterior fits, sets of draws, profiles or designs, such
#' as the analyses of several outcomes, effect measures or subsets of sites,
#' and tabulates their diagnostics side by side.
#'
#' @param ... For `compare_profiles()`, named objects of class
#'   `reliability_posterior`, `reliability_draws`, `reliability_profile` or
#'   `reliability_design`. For `plot()`, further arguments of the plot
#'   function.
#' @param .list Optionally, a named list of such objects, added after those
#'   in `...`.
#'
#' @details
#' The table has one row per analysis with the informativeness, the mean
#' reliability, the dispersion, the concentration count and the logit
#' variance, and their posterior intervals where there are any. The
#' analyses can be of different kinds, which the column `source` records.
#' Each analysis is conditional on its own standard errors, so analyses of
#' the same sites with different effect measures (a risk difference and a
#' log odds ratio, for example) have different diagnostics.
#'
#' @return An object of class `profile_comparison`, a list with the
#'   components
#'   * `objects`: the named list of analyses.
#'   * `table`: a data frame with one row per analysis and the columns
#'     `analysis`, `source`, `J`, `I`, `I_lower`, `I_upper`, `mean`,
#'     `mean_lower`, `mean_upper`, `H`, `H_lower`, `H_upper`, `J_eff`, `V`,
#'     `prob_I_above` and `sigma_tau` (posterior means, or values at a fixed
#'     `sigma_tau`, with the ends of the posterior intervals).
#'
#'   `print()`, `as.data.frame()` (the table), [profile_summary()],
#'   [report_table()], [report_text()] and `plot()` take it.
#'
#' @seealso [prior_sensitivity()] for the same analysis under several
#'   priors, and `vignette("multisite-trial")` for a comparison of the
#'   analyses of one trial.
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools)
#' prof <- reliability_profile(fit)
#' cmp <- compare_profiles(posterior = fit, `plug-in` = prof)
#' cmp
#' @export
compare_profiles <- function(..., .list = NULL) {
  objects <- c(list(...), .list)
  if (length(objects) == 0L) {
    .rp_abort("Give at least one analysis to compare.")
  }
  names <- names(objects)
  if (is.null(names) || any(!nzchar(names)) || anyDuplicated(names) > 0L) {
    .rp_abort("Every analysis must have a distinct name, as in `compare_profiles(a = fit_a, b = fit_b)`.",
              info = if (is.null(names)) "None has a name." else
                sprintf("The names are %s.", paste(sprintf("\"%s\"", names), collapse = ", ")))
  }
  classes <- c("reliability_posterior", "reliability_draws", "reliability_profile", "reliability_design")
  for (name in names) {
    if (!inherits(objects[[name]], classes)) {
      .rp_abort(sprintf("`%s` is not an analysis of poolinfo.", name),
                info = .rp_you_passed(objects[[name]]))
    }
  }
  structure(list(objects = objects, table = .rp_comparison_table(objects)), class = "profile_comparison")
}

.rp_comparison_table <- function(objects) {
  rows <- lapply(names(objects), function(name) {
    x <- objects[[name]]
    ps <- profile_summary(x)
    if (inherits(x, "reliability_design")) ps <- ps[ps$allocation == "planned", , drop = FALSE]
    get <- function(q, col) {
      out <- ps[[col]][ps$quantity == q]
      if (length(out) == 0L) NA_real_ else out
    }
    source <- if (inherits(x, "reliability_posterior")) "posterior" else
      if (inherits(x, "reliability_draws")) "posterior draws" else
        if (inherits(x, "reliability_design")) "planned design" else x$source
    probability <- if (inherits(x, c("reliability_posterior", "reliability_draws"))) {
      as.numeric(prob_above(x))
    } else NA_real_
    data.frame(
      analysis = name, source = source, J = length(x$se2),
      I = get("I", "mean"), I_lower = get("I", "lower"), I_upper = get("I", "upper"),
      mean = get("mean", "mean"), mean_lower = get("mean", "lower"), mean_upper = get("mean", "upper"),
      H = get("H", "mean"), H_lower = get("H", "lower"), H_upper = get("H", "upper"),
      J_eff = get("J_eff", "mean"), V = get("V", "mean"), prob_I_above = probability,
      sigma_tau = get("sigma_tau", "mean")
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @export
print.profile_comparison <- function(x, ..., use_cli = NULL) {
  cat(format(x, ..., use_cli = use_cli), sep = "\n")
  invisible(x)
}

#' @export
format.profile_comparison <- function(x, ..., use_cli = NULL) {
  .rp_check_dots_empty(...)
  t <- x$table
  interval <- function(m, lo, hi, d) {
    ifelse(is.na(lo), .rp_fmt(m, d), sprintf("%s [%s, %s]", .rp_fmt(m, d), .rp_fmt(lo, d), .rp_fmt(hi, d)))
  }
  rows <- data.frame(
    analysis = t$analysis, J = as.character(t$J),
    I = interval(t$I, t$I_lower, t$I_upper, 3L),
    mean = .rp_fmt(t$mean, 3L),
    H = interval(t$H, t$H_lower, t$H_upper, 4L),
    J_eff = .rp_fmt(t$J_eff, 1L)
  )
  cells <- rbind(c("analysis", "J", "informativeness I", "mean reliability", "dispersion H",
                   "concentration count J_eff"),
                 as.matrix(rows))
  widths <- apply(cells, 2L, function(col) max(nchar(col)))
  width_of <- function(cols) 2L + sum(widths[cols]) + 2L * (length(cols) - 1L)
  # A table wider than 80 characters is split: the centroids, then the
  # dispersion and the concentration count, each after the names; a part
  # that is still too wide is split again, column by column.
  groups <- if (width_of(1:6) <= 80L) list(2:6) else list(2:4, 5:6)
  groups <- unlist(lapply(groups, function(g) {
    out <- list(integer())
    for (k in g) {
      last <- out[[length(out)]]
      if (length(last) > 0L && width_of(c(1L, last, k)) > 80L) out <- c(out, list(k)) else
        out[[length(out)]] <- c(last, k)
    }
    out
  }), recursive = FALSE)
  groups <- lapply(groups, function(g) c(1L, g))
  table <- lapply(groups, function(cols) {
    apply(cells[, cols, drop = FALSE], 1L, function(r) {
      paste0("  ", paste(sprintf("%-*s", widths[cols], r), collapse = "  "))
    })
  })
  posterior <- any(vapply(x$objects, inherits, NA, what = c("reliability_posterior", "reliability_draws")))
  lines <- c(
    sprintf("Comparison of %d analyses", nrow(t)),
    if (posterior) "Intervals are central posterior intervals.",
    "Each analysis is conditional on its own standard errors.",
    unlist(lapply(table, function(block) c("", block)))
  )
  .rp_style_lines(lines, bold = 1L, use_cli = use_cli)
}

#' @export
as.data.frame.profile_comparison <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$table
}

#' @export
profile_summary.profile_comparison <- function(x, level = NULL, ...) {
  .rp_check_dots_empty(...)
  columns <- c("quantity", "mean", "sd", "lower", "upper", "mcse")
  out <- lapply(names(x$objects), function(name) {
    ps <- profile_summary(x$objects[[name]], level = level)
    if ("allocation" %in% names(ps)) ps <- ps[ps$allocation == "planned", , drop = FALSE]
    for (col in setdiff(columns, names(ps))) ps[[col]] <- NA_real_
    cbind(analysis = name, ps[, columns])
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' @rdname compare_profiles
#' @param x A comparison made by `compare_profiles()`.
#' @param type The plot of a comparison: `"functionals"` (the posterior
#'   densities of the informativeness, for posterior fits and draws only),
#'   `"plane"`, `"compare"`, `"profile"`, `"gap"` or `"concentration"`.
#'   `NULL`, the default, gives `"functionals"` when every analysis is a
#'   posterior fit or draws, and `"compare"` otherwise.
#' @export
plot.profile_comparison <- function(x, type = NULL, ...) {
  types <- c("functionals", "plane", "compare", "profile", "gap", "concentration")
  if (is.null(type)) {
    posterior <- all(vapply(x$objects, inherits, NA, what = c("reliability_posterior", "reliability_draws")))
    type <- if (posterior) "functionals" else "compare"
  }
  type <- .rp_match(type, types, "type")
  switch(type,
    functionals = plot_functionals(x, ...),
    plane = plot_plane(x, ...),
    compare = plot_compare(x, ...),
    profile = .rp_plot_compare_profiles(x, ...),
    gap = plot_gap_path(x, ...),
    concentration = plot_concentration(x, ...)
  )
}

#' Posterior summaries under several priors
#'
#' Refits a posterior under each of several priors on the between-site
#' standard deviation and collects the fits in a comparison, to see how much
#' the diagnostics depend on the prior.
#'
#' @param x An object of class `reliability_posterior`.
#' @param priors A named list of priors made by [sd_prior()].
#'
#' @details
#' With few sites the flat prior on the standard deviation can matter
#' (Gelman, 2006): its posterior has a heavy upper tail. A half-normal or
#' half-t prior whose scale is set on the scale of the effects is a common
#' alternative. The comparison shows how far the conclusions move; it does
#' not choose a prior. The scale of a proper prior has no default, because it
#' depends on the units of the estimates.
#'
#' @return An object of class `profile_comparison` whose first analysis is
#'   `x` itself, named after its prior.
#'
#' @references
#' Gelman, A. (2006). Prior distributions for variance parameters in
#' hierarchical models (comment on article by Browne and Draper). *Bayesian
#' Analysis, 1*(3), 515-534. \doi{10.1214/06-BA117A}
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools)
#' prior_sensitivity(fit, list(
#'   `half-normal, scale 10` = sd_prior("half_normal", scale = 10),
#'   `half-t, scale 10, 3 df` = sd_prior("half_t", scale = 10, df = 3)
#' ))
#' @export
prior_sensitivity <- function(x, priors) {
  if (!inherits(x, "reliability_posterior")) {
    .rp_abort("`prior_sensitivity()` takes a fit of `reliability_posterior()`.", info = .rp_you_passed(x))
  }
  if (!is.list(priors) || inherits(priors, "sd_prior") || length(priors) == 0L ||
      is.null(names(priors)) || any(!nzchar(names(priors)))) {
    .rp_abort("`priors` must be a named list of priors made by `sd_prior()`.", info = .rp_you_passed(priors))
  }
  for (name in names(priors)) .rp_check_prior(priors[[name]])
  fits <- lapply(priors, function(pr) stats::update(x, prior = pr))
  base_name <- switch(x$prior$type, flat = "flat", x$prior$type)
  if (base_name %in% names(fits)) base_name <- paste(base_name, "(the fit given)")
  compare_profiles(.list = c(stats::setNames(list(x), base_name), fits))
}
