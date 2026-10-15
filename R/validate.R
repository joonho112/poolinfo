# Checks of the arguments. Each returns the value it checked, possibly in a
# normalized form, or stops with an error of class rp_input_error.

.rp_is_number <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

# A single finite number in a range. `rule` says the range in words for the
# error, for example "a single number between 0 and 1".
.rp_check_number <- function(x, arg, rule, lower = -Inf, upper = Inf,
                             lower_open = FALSE, upper_open = FALSE,
                             null_ok = FALSE) {
  if (is.null(x) && null_ok) {
    return(NULL)
  }
  ok <- .rp_is_number(x) &&
    (if (lower_open) x > lower else x >= lower) &&
    (if (upper_open) x < upper else x <= upper)
  if (!ok) {
    .rp_abort(
      sprintf("`%s` must be %s%s.", arg, rule, if (null_ok) ", or NULL" else ""),
      info = .rp_you_passed(x)
    )
  }
  as.numeric(x)
}

.rp_check_level <- function(x, arg = "level", null_ok = FALSE) {
  .rp_check_number(x, arg, "a single number between 0 and 1", 0, 1,
                   lower_open = TRUE, upper_open = TRUE, null_ok = null_ok)
}

.rp_check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .rp_abort(sprintf("`%s` must be TRUE or FALSE.", arg), info = .rp_you_passed(x))
  }
  x
}

.rp_check_whole <- function(x, arg, lower = 1, upper = Inf, null_ok = FALSE) {
  if (is.null(x) && null_ok) {
    return(NULL)
  }
  limit <- min(upper, .Machine$integer.max)
  if (!.rp_is_number(x) || x != floor(x) || x < lower || x > limit) {
    rule <- if (is.finite(upper)) {
      sprintf("an integer between %s and %s", format(lower), format(upper))
    } else {
      sprintf("an integer of at least %s", format(lower))
    }
    .rp_abort(sprintf("`%s` must be %s.", arg, rule), info = .rp_you_passed(x))
  }
  as.integer(x)
}

.rp_check_string <- function(x, arg, null_ok = TRUE) {
  if (is.null(x) && null_ok) {
    return(NULL)
  }
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    .rp_abort(
      sprintf("`%s` must be a single string%s.", arg, if (null_ok) ", or NULL" else ""),
      info = .rp_you_passed(x)
    )
  }
  x
}

# One of a set of strings. The whole set, as in the usage of a function,
# means the first one.
.rp_match <- function(x, choices, arg) {
  if (identical(x, choices)) {
    return(choices[[1L]])
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% choices) {
    .rp_abort(
      sprintf("`%s` must be one of %s.", arg,
              paste(sprintf("\"%s\"", choices), collapse = ", ")),
      info = .rp_you_passed(x)
    )
  }
  x
}

# The dots of a method that takes no further arguments.
.rp_check_dots_empty <- function(...) {
  if (...length() > 0L) {
    # The names are read without evaluating the arguments, which can be
    # column names that exist only in `data`.
    names <- ...names()
    named <- if (is.null(names)) character() else names[!is.na(names) & nzchar(names)]
    unnamed <- ...length() - length(named)
    shown <- c(
      if (length(named) > 0L) paste(sprintf("`%s`", named), collapse = ", "),
      if (unnamed > 0L) sprintf("%d unnamed argument%s", unnamed, if (unnamed > 1L) "s" else "")
    )
    .rp_abort(
      "This function takes no further arguments.",
      info = sprintf("Not used: %s. Check the spelling of the argument names.", paste(shown, collapse = " and "))
    )
  }
  invisible(TRUE)
}

# Standard errors and sampling variances are named arguments. An unnamed
# argument in their place could hold either, so it is refused.
.rp_check_named_variances <- function(...) {
  names <- ...names()
  if (...length() > 0L && (is.null(names) || any(is.na(names) | !nzchar(names)))) {
    .rp_abort(
      "An unnamed argument after the first is ambiguous.",
      info = "It could hold standard errors or sampling variances.",
      fix = "Name it: `se = ` for standard errors or `se2 = ` for sampling variances."
    )
  }
  .rp_check_dots_empty(...)
}

# A vector of finite numbers. `positive` asks for numbers above 0. `arg` is
# the name of an argument, shown in backticks, or a phrase such as "the
# sampling variances of the fit", shown as it is.
.rp_check_vector <- function(x, arg, positive = FALSE, min_length = 1L) {
  subject <- if (grepl(" ", arg, fixed = TRUE)) {
    paste0(toupper(substr(arg, 1L, 1L)), substr(arg, 2L, nchar(arg)))
  } else {
    sprintf("`%s`", arg)
  }
  if (!is.numeric(x) || (length(dim(x)) > 1L && min(dim(x)) > 1L)) {
    .rp_abort(sprintf("%s must be a numeric vector.", subject), info = .rp_you_passed(x))
  }
  x <- as.vector(x, mode = "double")
  if (length(x) < min_length) {
    .rp_abort(
      sprintf("%s must hold at least %d number%s.", subject, min_length,
              if (min_length > 1L) "s" else ""),
      info = .rp_you_passed(x)
    )
  }
  bad <- which(!is.finite(x))
  if (length(bad) > 0L) {
    .rp_abort(
      sprintf("%s must hold finite numbers.", subject),
      info = sprintf("Position%s %s %s missing or not finite.",
                     if (length(bad) > 1L) "s" else "", .rp_positions(bad),
                     if (length(bad) > 1L) "are" else "is")
    )
  }
  if (positive) {
    bad <- which(x <= 0)
    if (length(bad) > 0L) {
      what <- if (all(x[bad] == 0)) "0" else if (all(x[bad] < 0)) "negative" else "0 or negative"
      .rp_abort(
        sprintf("%s must hold positive numbers.", subject),
        info = paste0(
          sprintf("Position%s %s %s %s.", if (length(bad) > 1L) "s" else "",
                  .rp_positions(bad), if (length(bad) > 1L) "are" else "is", what),
          if (any(x[bad] == 0)) {
            " A site without sampling error has reliability 1, whose logit is not finite."
          }
        )
      )
    }
  }
  x
}

# The sampling variances from the standard errors `se` or from the variances
# `se2`: exactly one of the two. Names are kept.
.rp_variances <- function(se, se2, min_length = 1L) {
  if (!is.null(se) && !is.null(se2)) {
    .rp_abort(
      "`se` and `se2` were both supplied.",
      info = "Give the standard errors or the sampling variances, not both."
    )
  }
  if (is.null(se) && is.null(se2)) {
    .rp_abort(
      "The standard errors are missing.",
      info = "Give the standard errors as `se` or the sampling variances as `se2`."
    )
  }
  if (!is.null(se)) {
    v <- .rp_check_vector(se, "se", positive = TRUE, min_length = min_length)^2
    bad <- which(v == 0 | !is.finite(v))
    if (length(bad) > 0L) {
      .rp_abort(
        "The squares of `se` must be positive finite numbers.",
        info = sprintf("At position%s %s the square is %s.", if (length(bad) > 1L) "s" else "",
                       .rp_positions(bad), if (all(v[bad] == 0)) "0" else "not finite"),
        fix = "Rescale the estimates and their standard errors."
      )
    }
    names(v) <- names(se)
  } else {
    v <- .rp_check_vector(se2, "se2", positive = TRUE, min_length = min_length)
    names(v) <- names(se2)
  }
  v
}

# Reliabilities: finite numbers strictly between 0 and 1.
.rp_check_reliabilities <- function(S, arg = "S") {
  values <- .rp_check_vector(S, arg)
  bad <- which(values <= 0 | values >= 1)
  if (length(bad) > 0L) {
    .rp_abort(
      sprintf("`%s` must hold reliabilities strictly between 0 and 1.", arg),
      info = sprintf("Position%s %s %s outside the open interval (0, 1).",
                     if (length(bad) > 1L) "s" else "", .rp_positions(bad),
                     if (length(bad) > 1L) "are" else "is"),
      fix = if (any(values[bad] %in% c(0, 1))) {
        "A reliability of 0 or 1 has no finite logit, so the centroids are not defined."
      }
    )
  }
  names(values) <- names(S)
  values
}

.rp_check_sigma_tau <- function(sigma_tau, arg = "sigma_tau") {
  .rp_check_number(sigma_tau, arg, "a single number of at least 0", lower = 0)
}

.rp_abort_missing_sigma_tau <- function() {
  .rp_abort(
    "`sigma_tau` is missing.",
    info = "Give the standard deviation of the true effects across sites."
  )
}

# Labels of the sites: NULL, or one distinct label for each of `J` sites.
.rp_check_site <- function(site, J, arg = "site") {
  if (is.null(site)) {
    return(NULL)
  }
  if (!is.atomic(site) || length(site) != J) {
    .rp_abort(
      sprintf("`%s` must hold one label for each of the %d sites.", arg, J),
      info = .rp_you_passed(site)
    )
  }
  site <- as.character(site)
  if (anyNA(site)) {
    .rp_abort(
      sprintf("`%s` must not be missing.", arg),
      info = sprintf("Position%s %s %s missing.",
                     if (sum(is.na(site)) > 1L) "s" else "",
                     .rp_positions(which(is.na(site))),
                     if (sum(is.na(site)) > 1L) "are" else "is")
    )
  }
  if (anyDuplicated(site) > 0L) {
    dup <- unique(site[duplicated(site)])
    .rp_abort(
      sprintf("`%s` must hold distinct labels.", arg),
      info = sprintf("Repeated: %s.",
                     paste(sprintf("\"%s\"", utils::head(dup, 5L)), collapse = ", "))
    )
  }
  site
}

# The labels of `J` sites: those given, the names of a vector, or 1 to J.
.rp_site_labels <- function(site, J, names = NULL) {
  if (!is.null(site)) {
    return(site)
  }
  if (!is.null(names) && !anyNA(names) && all(nzchar(names)) && anyDuplicated(names) == 0L) {
    return(as.character(names))
  }
  as.character(seq_len(J))
}
