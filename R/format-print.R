# Printed output. Every object has a format() method that returns its lines
# and a print() method that writes them. The words are those of the
# documentation: informativeness I, mean reliability, dispersion H, logit
# variance V, concentration count J_eff, between-site SD sigma_tau. No
# verdicts are printed.

# Numbers as text. Bounded summaries have 3 decimals, H 4, J_eff 1; values on
# the scale of the effects have 3 significant digits.
.rp_fmt <- function(x, digits) {
  out <- formatC(x, format = "f", digits = digits)
  out[is.na(x)] <- "NA"
  out
}

.rp_fmt_sig <- function(x, digits = 3L) {
  out <- vapply(x, function(z) {
    if (!is.finite(z)) return(format(z))
    s <- format(signif(z, digits), trim = TRUE)
    # A value with no more digits than shown, as a value that the user gave
    # (sigma_tau = 0.2, scale 10), is printed as it is. A rounded value keeps
    # the trailing zeros of its significant digits: -0.130, not -0.13.
    if (z == 0 || signif(z, digits) == z || grepl("e", s, fixed = TRUE)) return(s)
    shown <- sub("^0+", "", gsub(".", "", sub("^-", "", s), fixed = TRUE))
    missing <- digits - nchar(shown)
    if (missing > 0L) {
      if (!grepl(".", s, fixed = TRUE)) s <- paste0(s, ".")
      s <- paste0(s, strrep("0", missing))
    }
    s
  }, "")
  out[is.na(x)] <- "NA"
  out
}

# Whether to style the output with cli.
.rp_use_cli <- function(use_cli = NULL) {
  if (!is.null(use_cli) &&
      (!is.logical(use_cli) || length(use_cli) != 1L || is.na(use_cli))) {
    .rp_abort("`use_cli` must be NULL, TRUE or FALSE.", info = .rp_you_passed(use_cli))
  }
  if (is.null(use_cli)) {
    use_cli <- getOption("poolinfo.use_cli", TRUE)
  }
  isTRUE(use_cli) && requireNamespace("cli", quietly = TRUE)
}

# Bold headings when cli is used. Trailing blanks are dropped.
.rp_style_lines <- function(lines, bold = integer(), use_cli = NULL) {
  lines <- sub("[[:space:]]+$", "", lines)
  if (.rp_use_cli(use_cli) && length(bold) > 0L) {
    lines[bold] <- cli::style_bold(lines[bold])
  }
  lines
}

.rp_sites_word <- function(J) {
  if (J == 1L) "1 site" else sprintf("%d sites", J)
}

# Label and value in two columns.
.rp_rows <- function(labels, values, width = 28L) {
  sprintf("  %-*s %s", width, labels, values)
}

# Wraps a sentence at its spaces into lines of at most `width` characters,
# with `exdent` spaces before every line after the first. Unlike strwrap(), it
# keeps the spacing of the text. A no-break space (\u00a0) joins words that
# must stay on one line, such as "sigma_tau = 0.2", and prints as a space.
.rp_wrap <- function(text, width = 80L, exdent = 0L) {
  words <- strsplit(text, " ", fixed = TRUE)[[1L]]
  if (length(words) == 0L) return(text)
  lines <- character()
  line <- words[[1L]]
  for (word in words[-1L]) {
    if (nchar(line) + 1L + nchar(word) <= width) {
      line <- paste(line, word)
    } else {
      lines <- c(lines, line)
      line <- paste0(strrep(" ", exdent), word)
    }
  }
  gsub("\u00a0", " ", c(lines, line), fixed = TRUE)
}

# Joins the words of a phrase with no-break spaces for .rp_wrap().
.rp_nobreak <- function(text) gsub(" ", "\u00a0", text, fixed = TRUE)

# The heading of a profile; the source of sigma_tau goes to a second line
# when the two do not fit on one.
.rp_profile_heading <- function(x) {
  main <- sprintf("Reliability profile of %s at sigma_tau = %s", .rp_sites_word(length(x$se2)),
                  .rp_fmt_sig(x$sigma_tau))
  source <- sprintf("(%s)", x$source)
  if (nchar(main) + 1L + nchar(source) <= 80L) paste(main, source) else c(main, paste0("  ", source))
}

.rp_label_lines <- function(labels) {
  hang <- function(prefix, text) .rp_wrap(paste0(prefix, text), exdent = nchar(prefix))
  c(
    if (!is.null(labels$scale_label)) hang("Scale of the effects: ", labels$scale_label),
    if (!is.null(labels$target_label)) hang("Site estimates: ", labels$target_label)
  )
}

.rp_jeff_text <- function(jeff, J, limit = FALSE) {
  sprintf("%s of %d%s", .rp_fmt(jeff, 1L), J, if (limit) " (limit as sigma_tau goes to 0)" else "")
}

#' @rdname reliability_profile
#' @order 7
#' @param use_cli Whether to print the headings in bold with the cli package,
#'   when it is installed. `NULL`, the default, takes
#'   `getOption("poolinfo.use_cli", TRUE)`.
#' @export
format.reliability_profile <- function(x, ..., use_cli = NULL) {
  .rp_check_dots_empty(...)
  st <- x$stats
  J <- length(x$se2)
  zero <- x$sigma_tau == 0
  heading <- .rp_profile_heading(x)
  head <- c(heading, .rp_label_lines(x$labels))
  body <- .rp_rows(
    c("informativeness I", "mean reliability", "dispersion H", "logit variance V",
      "concentration count J_eff"),
    c(.rp_fmt(st[["I"]], 3L), .rp_fmt(st[["mean"]], 3L), .rp_fmt(st[["H"]], 4L),
      .rp_fmt(st[["V"]], 3L), .rp_jeff_text(st[["J_eff"]], J, zero))
  )
  range_line <- if (zero) {
    "At sigma_tau = 0 every reliability is 0."
  } else if (J == 1L) {
    sprintf("The site's reliability is %s.", .rp_fmt(x$S, 2L))
  } else {
    sprintf("Site reliabilities range from %s to %s.", .rp_fmt(min(x$S), 2L), .rp_fmt(max(x$S), 2L))
  }
  lines <- c(head, "", body, "", range_line)
  .rp_style_lines(lines, bold = seq_along(heading), use_cli = use_cli)
}

#' @rdname reliability_profile
#' @order 8
#' @export
print.reliability_profile <- function(x, ..., use_cli = NULL) {
  cat(format(x, ..., use_cli = use_cli), sep = "\n")
  invisible(x)
}

#' @rdname reliability_profile
#' @order 9
#' @export
summary.reliability_profile <- function(object, ...) {
  .rp_check_dots_empty(...)
  st <- object$stats
  centroids <- data.frame(
    centroid = c("harmonic", "I", "mean"),
    value = unname(st[c("harmonic", "I", "mean")]),
    definition = c("harmonic mean", "logit mean", "arithmetic mean")
  )
  related <- data.frame(
    quantity = c("benchmark_risk", "trace", "v_typ", "R_typ"),
    value = unname(st[c("benchmark_risk", "trace", "v_typ", "R_typ")])
  )
  structure(list(profile = object, centroids = centroids, local = object$local, related = related),
            class = "summary.reliability_profile")
}

#' @rdname reliability_profile
#' @order 10
#' @export
format.summary.reliability_profile <- function(x, ..., use_cli = NULL) {
  .rp_check_dots_empty(...)
  prof <- x$profile
  st <- prof$stats
  J <- length(prof$se2)
  zero <- prof$sigma_tau == 0
  heading <- .rp_profile_heading(prof)
  head <- c(heading, .rp_label_lines(prof$labels))
  centroids <- .rp_rows(
    c("harmonic centroid", "informativeness I", "mean reliability"),
    sprintf("%-8s (%s)", .rp_fmt(x$centroids$value, 3L), x$centroids$definition)
  )
  spread <- .rp_rows(
    c("dispersion H", "logit variance V", "concentration count J_eff"),
    c(.rp_fmt(st[["H"]], 4L), .rp_fmt(st[["V"]], 3L), .rp_jeff_text(st[["J_eff"]], J, zero))
  )
  lt <- x$local
  decided <- ifelse(is.na(lt$holds), "not defined", ifelse(lt$holds, "yes", "no"))
  # H is never negative, so the question of its sign does not arise.
  decided[lt$quantity == "H"] <- "-"
  g <- function(z) {
    out <- formatC(z, digits = 3L, format = "g", flag = "#")
    out <- sub("[.]$", "", out)
    ifelse(z < 0, out, paste0(" ", out))
  }
  local <- c(
    sprintf("  %-22s %-10s %-10s %-10s %s", "", "", "", "", "leading term"),
    sprintf("  %-22s %-10s %-10s %-10s %s", "", " exact", " leading", " bound", "decides the sign"),
    sprintf("  %-22s %-10s %-10s %-10s %s", c("H", "mean reliability - I"), g(lt$exact),
            g(lt$leading), g(lt$bound), decided)
  )
  rel <- x$related$value
  related <- .rp_rows(
    c("benchmark risk sigma_tau^2 (1 - mean S_j)", "sum of reliabilities",
      "typical variance (Higgins-Thompson)"),
    c(.rp_fmt_sig(rel[[1L]]), .rp_fmt(rel[[2L]], 1L),
      if (is.na(rel[[3L]])) "NA (needs 2 sites)" else
        sprintf("%-8s R_typ %s", .rp_fmt_sig(rel[[3L]]), .rp_fmt(rel[[4L]], 3L))),
    width = 44L
  )
  lines <- c(head, "", "Centroids", centroids, "", "Dispersion and concentration", spread, "",
             "Second-order terms (Lee, 2026)", local, "",
             "Related summaries", related)
  bold <- c(seq_along(heading), which(lines %in% c("Centroids", "Dispersion and concentration",
                                   "Second-order terms (Lee, 2026)",
                                   "Related summaries")))
  .rp_style_lines(lines, bold = bold, use_cli = use_cli)
}

#' @rdname reliability_profile
#' @order 11
#' @export
print.summary.reliability_profile <- function(x, ..., use_cli = NULL) {
  cat(format(x, ..., use_cli = use_cli), sep = "\n")
  invisible(x)
}

#' @rdname reliability_profile
#' @order 12
#' @param row.names,optional Ignored.
#' @export
as.data.frame.reliability_profile <- function(x, row.names = NULL, optional = FALSE, ...) {
  site_table(x)
}
