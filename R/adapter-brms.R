# Read the normal hierarchical model from brms metadata before converting draws.
.rp_brms_parts <- function(x) {
  scope <- function(info) {
    .rp_abort(
      "This brms fit is outside the supported normal hierarchical model.",
      info = info,
      fix = "Use a Gaussian identity-link response with known standard errors and one independent site intercept.",
      class = "rp_scope_error"
    )
  }
  f <- x$formula$formula
  if (!inherits(f, "formula") || length(f) != 3L ||
      length(x$formula$pforms) || length(x$formula$pfix) ||
      isTRUE(attr(f, "nl")) || length(attr(f, "autocor")) || length(x$autocor)) {
    scope("Only a univariate linear response without distributional or autocorrelation terms is supported.")
  }
  family <- if (!is.null(x$family)) x$family else x$formula$family
  if (!is.list(family) || !identical(family$family, "gaussian") ||
      !identical(family$link, "identity")) {
    scope("The response family must be Gaussian with an identity link.")
  }
  lhs <- f[[2L]]
  if (!is.call(lhs) || !identical(lhs[[1L]], as.name("|")) ||
      length(lhs) != 3L || !is.symbol(lhs[[2L]])) {
    scope("The response must be a direct data column with `se(se_var)`.")
  }
  se_call <- lhs[[3L]]
  if (!is.call(se_call) || !identical(se_call[[1L]], as.name("se")) ||
      !length(se_call) %in% c(2L, 3L) || !is.symbol(se_call[[2L]])) {
    scope("Known standard errors must be a direct data column in `se(se_var, sigma = FALSE)`.")
  }
  se_args <- as.list(se_call)[-1L]
  se_names <- names(se_args)
  if (!is.null(se_names) && nzchar(se_names[[1L]]) && se_names[[1L]] != "x") {
    scope("The first argument of `se()` must name the standard-error column.")
  }
  if (length(se_args) == 2L &&
      (!identical(se_args[[2L]], FALSE) ||
       (!is.null(se_names) && nzchar(se_names[[2L]]) && se_names[[2L]] != "sigma"))) {
    scope("`se(..., sigma = TRUE)` adds residual variance; only omitted or FALSE `sigma` is supported.")
  }
  unsupported_term <- function(expr) {
    if (!is.call(expr)) return(FALSE)
    head <- if (is.symbol(expr[[1L]])) as.character(expr[[1L]]) else ""
    if (!nzchar(head) && is.call(expr[[1L]]) &&
        identical(expr[[1L]][[1L]], as.name("::"))) {
      head <- as.character(expr[[1L]][[3L]])
    }
    if (head %in% c("arma", "ar", "ma", "cosy", "unstr", "sar", "car", "fcor",
                    "s", "t2", "gp", "mi", "me")) return(TRUE)
    any(vapply(as.list(expr)[-1L], unsupported_term, logical(1L)))
  }
  if (unsupported_term(f[[3L]])) {
    scope("Autocorrelation, smooths, Gaussian processes and latent measurement-error terms are not supported.")
  }
  groups <- list()
  has_group <- function(expr) {
    if (!is.call(expr)) return(FALSE)
    if (identical(expr[[1L]], as.name("|")) || identical(expr[[1L]], as.name("||"))) return(TRUE)
    any(vapply(as.list(expr)[-1L], has_group, logical(1L)))
  }
  strip_group <- function(expr) {
    if (!is.call(expr)) return(expr)
    head <- if (is.symbol(expr[[1L]])) as.character(expr[[1L]]) else ""
    if (head == "(") return(strip_group(expr[[2L]]))
    if (head == "+" && length(expr) == 3L) {
      left <- strip_group(expr[[2L]])
      right <- strip_group(expr[[3L]])
      if (is.null(left)) return(right)
      if (is.null(right)) return(left)
      return(call("+", left, right))
    }
    if (head %in% c("|", "||")) {
      if (length(expr) != 3L || !identical(expr[[2L]], 1) || !is.symbol(expr[[3L]])) {
        scope("Only a group-level intercept `(1 | site)` with a direct grouping column is supported.")
      }
      groups[[length(groups) + 1L]] <<- as.character(expr[[3L]])
      return(NULL)
    }
    if (has_group(expr)) {
      scope("Group-level terms must enter additively as one independent `(1 | site)` intercept.")
    }
    expr
  }
  fixed <- strip_group(f[[3L]])
  if (length(groups) != 1L) {
    scope("Exactly one independent group-level intercept is supported; additional grouping components are not supported.")
  }
  if (is.null(fixed)) fixed <- 1
  fixed_formula <- stats::as.formula(call("~", fixed), env = environment(f))
  fixed_terms <- tryCatch(stats::terms(fixed_formula, data = x$data), error = function(e) NULL)
  if (is.null(fixed_terms)) scope("The population mean formula could not be read.")
  labels <- attr(fixed_terms, "term.labels")
  has_intercept <- identical(attr(fixed_terms, "intercept"), 1L)
  # brms writes a common mean as `0 + Intercept` too.
  intercept_term <- !has_intercept && identical(labels, "Intercept")
  if (!has_intercept && !length(labels)) {
    scope("The formula has no intercept or coefficient for the mean of the effects; the model needs an intercept.")
  }
  intercept_only <- (has_intercept && !length(labels) || intercept_term) &&
    !length(attr(fixed_terms, "offset"))
  if (identical(attr(x$prior, "sample_prior"), "only")) {
    scope("The fit was made with `sample_prior = \"only\"`, so its draws are from the prior.")
  }
  response <- as.character(lhs[[2L]])
  se <- as.character(se_call[[2L]])
  group <- groups[[1L]]
  data <- x$data
  if (!is.data.frame(data) || !all(c(response, se, group) %in% names(data))) {
    scope("The response, standard-error and grouping columns must be present in the fit's data.")
  }
  if (!is.numeric(data[[response]]) || any(!is.finite(data[[response]])) ||
      !is.numeric(data[[se]]) || any(!is.finite(data[[se]])) || any(data[[se]] <= 0)) {
    scope("The response must be finite and the known standard errors must be finite and positive.")
  }
  if (anyNA(data[[group]]) || anyDuplicated(data[[group]])) {
    scope("Each site must have exactly one observation and a nonmissing group label.")
  }
  list(response = response, se = se, group = group,
       sd_variable = sprintf("sd_%s__Intercept", group),
       mean_variable = if (intercept_only) "b_Intercept" else NULL)
}
