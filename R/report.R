#' Table of the diagnostics for a report
#'
#' Formats the diagnostics of one analysis or of a comparison as a table in
#' the layout of Table 5 of Lee (2026): the informativeness, the mean
#' reliability and the dispersion with their posterior intervals, the
#' concentration count and the posterior probability that the
#' informativeness exceeds the threshold.
#'
#' @param x An object of class `reliability_posterior`,
#'   `reliability_draws`, `reliability_profile` or `profile_comparison`.
#' @param digits Number of decimals of the informativeness, the mean
#'   reliability and the probability, a positive integer; the dispersion has
#'   one more and the concentration count 1.
#' @param format `"data.frame"` (the default) for a data frame of numbers,
#'   `"markdown"` or `"latex"` for text in the style of the paper.
#'
#' @details
#' In the text formats the leading zeros of the numbers between 0 and 1 are
#' dropped and intervals are in brackets, as in the paper. A table of one
#' analysis has no column for its name; the names of a comparison label its
#' rows. The LaTeX table uses the rules of the booktabs package. Every
#' format ends with a note that states the working model, the prior, the
#' scale and the estimand when the fits record them, the conditioning on the
#' standard errors, and that each site counts equally in the diagnostics.
#'
#' @return For `format = "data.frame"`, a data frame with one row per
#'   analysis and the columns `analysis`, `J`, `I`, `I_lower`, `I_upper`,
#'   `mean`, `H`, `H_lower`, `H_upper`, `J_eff` and `prob_I_above`, and the
#'   note in its attribute `note`. For the text formats, a character string.
#'
#' @template ref-lee2026
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools)
#' report_table(fit)
#' cat(report_table(fit, format = "markdown"))
#' @export
report_table <- function(x, digits = 2, format = c("data.frame", "markdown", "latex")) {
  digits <- .rp_check_whole(digits, "digits")
  format <- .rp_match(format, c("data.frame", "markdown", "latex"), "format")
  single <- !inherits(x, "profile_comparison")
  cmp <- if (!single) x else {
    if (!inherits(x, c("reliability_posterior", "reliability_draws", "reliability_profile"))) {
      .rp_abort("`report_table()` takes a posterior fit, draws, a profile or a comparison.",
                info = .rp_you_passed(x))
    }
    compare_profiles(analysis = x)
  }
  t <- cmp$table
  note <- .rp_report_note(cmp$objects)
  if (format == "data.frame") {
    out <- t[, c("analysis", "J", "I", "I_lower", "I_upper", "mean", "H", "H_lower", "H_upper", "J_eff",
                 "prob_I_above")]
    attr(out, "note") <- note
    return(out)
  }
  num <- function(z, d) {
    out <- formatC(z, format = "f", digits = d)
    out <- sub("^(-?)0\\.", "\\1.", out)
    out[is.na(z)] <- if (format == "latex") "--" else "-"
    out
  }
  interval <- function(m, lo, hi, d) {
    ifelse(is.na(lo), num(m, d), sprintf("%s [%s, %s]", num(m, d), num(lo, d), num(hi, d)))
  }
  threshold <- vapply(cmp$objects, function(o) if (is.null(o$threshold)) NA_real_ else o$threshold, 0)
  prob_head <- if (length(unique(stats::na.omit(threshold))) == 1L) {
    sprintf("Pr(I > %s)", format(stats::na.omit(threshold)[[1L]]))
  } else "Pr(I > threshold)"
  cells <- data.frame(
    Analysis = t$analysis, J = as.character(t$J),
    I = interval(t$I, t$I_lower, t$I_upper, digits),
    mean = num(t$mean, digits),
    H = interval(t$H, t$H_lower, t$H_upper, digits + 1L),
    J_eff = num(t$J_eff, 1L),
    prob = num(t$prob_I_above, digits)
  )
  if (format == "markdown") {
    head <- c("Analysis", "J", "I [95% interval]", "Mean reliability", "H [95% interval]", "J_eff", prob_head)
    head[3] <- .rp_interval_head("I", cmp$objects)
    head[5] <- .rp_interval_head("H", cmp$objects)
    if (single) {
      cells <- cells[, -1L, drop = FALSE]
      head <- head[-1L]
    }
    rows <- apply(cells, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
    return(paste(c(paste0("| ", paste(head, collapse = " | "), " |"),
                   paste0("|", paste(rep(":--", length(head)), collapse = "|"), "|"),
                   rows, "", paste("*Note.*", note)), collapse = "\n"))
  }
  esc <- function(s) gsub("([%&#_])", "\\\\\\1", s)
  head <- c("Analysis", "$J$", sub("I \\[", "$I$ [", .rp_interval_head("I", cmp$objects)), "$\\bar S$",
            sub("H \\[", "$H$ [", .rp_interval_head("H", cmp$objects)), "$J_{\\mathrm{eff}}$",
            sub("^Pr\\(I > ", "$\\\\Pr(I > ", sub("\\)$", ")$", prob_head)))
  # A percent sign starts a comment in LaTeX.
  head <- gsub("%", "\\%", head, fixed = TRUE)
  spec <- "lrlrlrr"
  if (single) {
    cells <- cells[, -1L, drop = FALSE]
    head <- head[-1L]
    spec <- "rlrlrr"
  }
  rows <- apply(cells, 1L, function(r) paste0(paste(esc(r), collapse = " & "), " \\\\"))
  paste(c(
    sprintf("\\begin{tabular}{%s}", spec), "\\toprule",
    paste0(paste(head, collapse = " & "), " \\\\"), "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
    "", paste0("\\par\\smallskip\\noindent\\textit{Note.} ", esc(note))
  ), collapse = "\n")
}

.rp_interval_head <- function(q, objects) {
  levels <- unique(unlist(lapply(objects, function(o) o$level)))
  if (length(levels) == 1L) sprintf("%s [%s interval]", q, .rp_percent_label(levels)) else
    sprintf("%s [interval]", q)
}

# The note of a report: the working model, how sigma_tau was treated, the
# labels and the conditioning. What differs between the analyses of a
# comparison is given for each analysis.
.rp_report_note <- function(objects) {
  # A value shared by every analysis is stated once, in a sentence; values
  # that differ are listed with the names of the analyses.
  per_analysis <- function(values, sentence, lead) {
    values <- unlist(values)
    if (length(values) == 0L) return(NULL)
    if (length(unique(values)) == 1L && length(values) == length(objects)) {
      return(sprintf(sentence, values[[1L]]))
    }
    paste0(lead, " ", paste(sprintf("%s (%s)", .rp_uncapitalize(values), names(values)), collapse = "; "), ".")
  }
  labels <- function(name) Filter(Negate(is.null), lapply(objects, function(o) o$labels[[name]]))
  paste(c(
    "Each estimate is modeled as normal around its site's effect with a known sampling variance, and the site effects as independent and normal across sites.",
    per_analysis(lapply(objects, .rp_sigma_words), "%s.", "By analysis:"),
    per_analysis(labels("scale_label"), "The effects are on the %s scale.", "Scales of the effects:"),
    per_analysis(labels("target_label"), "The estimand at each site is %s.", "Estimands at each site:"),
    "Summaries are conditional on the supplied standard errors, and each site counts equally in the diagnostics, whatever its size."
  ), collapse = " ")
}

.rp_capitalize <- function(text) paste0(toupper(substring(text, 1L, 1L)), substring(text, 2L))
.rp_uncapitalize <- function(text) paste0(tolower(substring(text, 1L, 1L)), substring(text, 2L))

# How an analysis treats sigma_tau, in words for a report: the priors of a
# posterior fit, the source of draws, or the value of a profile.
.rp_sigma_words <- function(x) {
  if (inherits(x, "reliability_posterior")) {
    return(paste(if (x$prior$type == "flat") "Flat priors on" else "A flat prior on", .rp_prior_words(x)))
  }
  if (inherits(x, "reliability_draws")) {
    return(sprintf("Posterior draws of the between-site standard deviation%s",
                   if (is.null(x$model_label)) "" else paste0(" from ", x$model_label)))
  }
  if (inherits(x, "reliability_profile")) {
    return(sprintf("Values at a between-site standard deviation of %s%s", .rp_fmt_sig(x$sigma_tau),
                   .rp_source_words(x$source)))
  }
  NULL
}

# The parts of a prior in words: the mean or the coefficients, then the
# between-site standard deviation.
.rp_prior_words <- function(x) {
  prior <- x$prior
  mean_part <- if (is.null(x$data$X)) "the mean" else "the coefficients of the mean model"
  sd_part <- switch(prior$type,
    flat = "the between-site standard deviation",
    half_normal = sprintf("a half-normal prior with scale %s on the between-site standard deviation",
                          .rp_fmt_sig(prior$scale)),
    half_t = sprintf("a half-t prior with scale %s and %s degrees of freedom on the between-site standard deviation",
                     .rp_fmt_sig(prior$scale), .rp_fmt_sig(prior$df)),
    half_cauchy = sprintf("a half-Cauchy prior with scale %s on the between-site standard deviation",
                          .rp_fmt_sig(prior$scale)),
    uniform = sprintf("a uniform prior from 0 to %s on the between-site standard deviation", .rp_fmt_sig(prior$upper))
  )
  paste(mean_part, "and", sd_part)
}

# The source of the sigma_tau of a profile, in words, as a parenthesis.
.rp_source_words <- function(source) {
  if (identical(source, "supplied")) return("")
  if (grepl("sigma_tau^2", source, fixed = TRUE)) {
    source <- paste("the square root of the", gsub("sigma_tau^2", "the between-site variance", source, fixed = TRUE))
  }
  sprintf(" (%s)", source)
}

#' Sentences for a report
#'
#' Writes the diagnostics of an analysis as sentences for the results
#' section of a paper, to be edited into the paper's own prose.
#'
#' @param x An object of class `reliability_posterior`,
#'   `reliability_draws`, `reliability_profile`, `reliability_design` or
#'   `profile_comparison`.
#' @param outcome Optional name of the outcome, a string, as in
#'   `"course passing"`.
#' @param digits Number of decimals of the informativeness, the mean
#'   reliability and the probability, a positive integer; the dispersion has
#'   one more.
#'
#' @details
#' The sentences open with the scale of the effects and the estimand at each
#' site when the fit records them (`scale_label` and `target_label`). For a
#' posterior fit or draws, the next sentence states the working model and
#' the prior; for a profile, the between-site standard deviation and where
#' it came from; and for a design, the planned allocation and its
#' comparison with balance or equal spending. When the threshold is 1/2, a
#' sentence says what that value of the informativeness means: the
#' between-site variance equals the geometric mean of the sampling
#' variances. The last sentence states that the summaries are conditional
#' on the supplied standard errors. The sentences give numbers and
#' intervals without describing them as low or high.
#'
#' @return A character vector with one sentence per element; for a
#'   comparison, the sentences of each analysis in turn.
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, data = eight_schools,
#'                              scale_label = "SAT-V point",
#'                              target_label = "the effect of coaching on the SAT-V")
#' report_text(fit, outcome = "the SAT-V")
#' @export
report_text <- function(x, outcome = NULL, digits = 2) {
  outcome <- .rp_check_string(outcome, "outcome")
  digits <- .rp_check_whole(digits, "digits")
  if (inherits(x, "profile_comparison")) {
    return(unlist(lapply(names(x$objects), function(name) {
      s <- report_text(x$objects[[name]], outcome = outcome, digits = digits)
      c(sprintf("%s:", name), s)
    }), use.names = FALSE))
  }
  f <- function(z, d = digits) .rp_fmt(z, d)
  of_outcome <- if (is.null(outcome)) "" else paste0(" for ", outcome)
  if (inherits(x, "reliability_design")) {
    p <- x$planned$stats
    b <- x$balanced$stats
    return(c(
      sprintf("%s of %s individuals across %d sites would have an informativeness of %s at a between-site standard deviation of %s, against %s for %s of the same %s.",
              if (is.null(outcome)) "The planned allocation" else sprintf("For %s, the planned allocation", outcome),
              format(round(sum(x$n)), big.mark = ","), length(x$n), f(p[["I"]]), .rp_fmt_sig(x$sigma_tau[[1L]]),
              f(b[["I"]]), if (is.null(x$cost)) "a balanced allocation" else "equal spending",
              if (is.null(x$cost)) "total" else "budget"),
      sprintf("The mean reliability would be %s (%s %s), and the concentration count %s (at most %d, the number of sites).",
              f(p[["mean"]]), f(b[["mean"]]), if (is.null(x$cost)) "balanced" else "with equal spending",
              .rp_fmt(p[["J_eff"]], 1L), length(x$n)),
      if (is.null(x$cost)) {
        "At a fixed number of sites and total, balance maximizes the informativeness and the mean reliability when every site has the same variance constant and the same cost per individual."
      } else {
        "With costs per individual that differ, equal spending maximizes the informativeness for the same budget."
      },
      "These values are computed at the planning values of the between-site standard deviation and the sampling variances."
    ))
  }
  labels <- x$labels
  scale_sentence <- if (!is.null(labels$scale_label)) sprintf("The effects are on the %s scale.", labels$scale_label)
  target_sentence <- if (!is.null(labels$target_label)) sprintf("The estimand at each site is %s.", labels$target_label)
  ps <- profile_summary(x)
  get <- function(q, col = "mean") ps[[col]][ps$quantity == q]
  J <- length(x$se2)
  if (inherits(x, "reliability_profile")) {
    return(c(
      scale_sentence, target_sentence,
      sprintf("At a between-site standard deviation of %s%s, the informativeness of the %d sites%s was %s and the mean reliability %s.",
              .rp_fmt_sig(x$sigma_tau), .rp_source_words(x$source), J, of_outcome, f(get("I")), f(get("mean"))),
      sprintf("The dispersion of the site reliabilities was H = %s, with logit variance V = %s, and the concentration count was %s (at most %d, the number of sites).",
              f(get("H"), digits + 1L), f(get("V")), .rp_fmt(get("J_eff"), 1L), J),
      "These values are conditional on the supplied standard errors and on the between-site standard deviation."
    ))
  }
  pct <- .rp_percent_label(attr(ps, "level"))
  model <- if (inherits(x, "reliability_posterior")) {
    paste("Under the normal working model with", .rp_uncapitalize(.rp_sigma_words(x)))
  } else {
    sprintf("With posterior draws of the between-site standard deviation%s",
            if (is.null(x$model_label)) "" else paste0(" from ", x$model_label))
  }
  c(
    scale_sentence, target_sentence,
    sprintf("%s, the informativeness of the %d sites%s was %s (%s posterior interval %s to %s), and the posterior probability that it exceeded %s was %s.",
            model, J, of_outcome, f(get("I")), pct, f(get("I", "lower")), f(get("I", "upper")),
            format(x$threshold), f(as.numeric(prob_above(x)))),
    if (x$threshold == 0.5) {
      "At an informativeness of 0.5, the between-site variance equals the geometric mean of the sampling variances."
    },
    sprintf("The mean reliability was %s (%s to %s).", f(get("mean")), f(get("mean", "lower")), f(get("mean", "upper"))),
    sprintf("The dispersion of the site reliabilities was H = %s (%s to %s), with logit variance V = %s, and the concentration count was %s (at most %d, the number of sites).",
            f(get("H"), digits + 1L), f(get("H", "lower"), digits + 1L), f(get("H", "upper"), digits + 1L),
            f(get("V")), .rp_fmt(get("J_eff"), 1L), J),
    "These summaries are conditional on the supplied standard errors, and each site counts equally in them, whatever its size."
  )
}
