#' Reliability profile at a fixed between-site SD
#'
#' Computes the site reliabilities and the diagnostics of their profile at a
#' given between-site standard deviation: a planning value, the true value of
#' a simulation, or an estimate from another fit.
#'
#' @param x Usually left out, with the standard errors given as `se` or
#'   `se2`. It can also be a data frame, in which `se`, `se2`, `site` and
#'   `estimate` name columns; a posterior fit or posterior draws, for the
#'   profile at the posterior mean of \eqn{\sigma_\tau^2}{sigma_tau^2}; or a
#'   fit of metafor, ebnm or multisiteDGP (see [metafor-methods],
#'   [ebnm-methods] and [multisitedgp-methods]).
#' @param ... Must be empty.
#' @param se Standard errors of the site estimates, a numeric vector of
#'   positive numbers. With a data frame `x`, an unquoted column name. Give
#'   `se` or `se2`, not both.
#' @param se2 Sampling variances of the site estimates, as `se`.
#' @param sigma_tau Standard deviation of the true effects across sites, a
#'   number of at least 0.
#' @param site Optional labels of the sites, one for each site and all
#'   distinct. With a data frame `x`, an unquoted column name.
#' @param estimate Optional estimates of the sites, a numeric vector. When
#'   they are given, the site table also has the posterior mean and
#'   standard deviation of each site's effect at `sigma_tau`.
#' @param data Optional data frame in which `se`, `se2`, `site` and
#'   `estimate` are found, as in [reliability_posterior()]; the same as
#'   giving the data frame as `x`.
#' @param scale_label,target_label Optional strings that name the scale of
#'   the effects, such as `"risk difference"` (read as "on the risk difference
#'   scale"), and the estimand at each site, such as `"the block-adjusted
#'   intention-to-treat contrast for course passing"`. They appear in the
#'   printed output and in [report_text()].
#'
#' @details
#' With sampling variances \eqn{v_j} and the between-site standard
#' deviation \eqn{\sigma_\tau}{sigma_tau}, the reliability of site \eqn{j}
#' is \eqn{S_j = \sigma_\tau^2/(\sigma_\tau^2 + v_j)}{S_j = sigma_tau^2 / (sigma_tau^2 + v_j)}.
#' The profile's diagnostics are the informativeness \eqn{I}
#' ([kl_centroid()]), the mean reliability, the harmonic centroid, the
#' dispersion \eqn{H} ([kl_dispersion()]), the logit variance \eqn{V}
#' ([logit_moments()]), the concentration count
#' \eqn{J_{\mathrm{eff}}}{J_eff} ([concentration_count()]) and the terms of
#' [local_terms()]. `summary()` adds the benchmark risk, the sum of the
#' reliabilities and the typical variance of [typical_variance()].
#'
#' This is the profile at one value of `sigma_tau`. When `sigma_tau` is
#' estimated from the same data, [reliability_posterior()] averages the
#' diagnostics over its posterior distribution instead, which is what we
#' recommend for reporting.
#' Profiles obtained from a posterior or a metafor fit with moderators
#' retain its mean-model matrix. Their site effects pool toward the
#' site's fitted mean, with uncertainty in the mean-model coefficients
#' integrated out under a flat prior.
#'
#' At `sigma_tau = 0` every reliability is 0 and the profile is a limit:
#' \eqn{I}, the mean reliability and \eqn{H} are 0, and the concentration
#' count is its limit \eqn{(\sum_j v_j^{-1})^2 / \sum_j v_j^{-2}}{sum(1/v)^2 / sum(1/v^2)},
#' and a message is shown.
#'
#' @return An object of class `reliability_profile`, a list with components
#'   * `se2`, `site` and `estimate` (or `NULL`): the data.
#'   * `X`: the mean-model matrix inherited from a fitted model, or `NULL`
#'     for a common mean.
#'   * `sigma_tau` and `source`: the between-site SD and where it came from.
#'   * `S`: the reliabilities.
#'   * `stats`: a named vector of the diagnostics: `I`, `mean`, `harmonic`,
#'     `H`, `V`, `rho3` (the third absolute moment of the centered logits),
#'     `J_eff`, `gap` (the mean reliability minus \eqn{I}), `benchmark_risk`
#'     (see [benchmark_risk()]), `trace` (the sum of the reliabilities),
#'     `v_typ` (the typical variance of [typical_variance()]) and `R_typ`
#'     (the reliability at that variance).
#'   * `local`: the table of [local_terms()].
#'   * `labels`: the scale and estimand labels.
#'   * `call` and `version`: the call and the version of the package that
#'     made the profile.
#'
#'   It has methods for `print()`, `summary()`, `as.data.frame()` (the site
#'   table), `plot()` and `update()`.
#'
#' @seealso [profile_summary()] and [site_table()] for the numbers as data
#'   frames, [reliability_posterior()] for the posterior over `sigma_tau`,
#'   and `vignette("reliability-profiles")` for what the diagnostics
#'   describe.
#'
#' @template ref-lee2026
#' @examples
#' se <- c(0.10, 0.12, 0.15, 0.20, 0.25, 0.30)
#' prof <- reliability_profile(se = se, sigma_tau = 0.2)
#' prof
#' summary(prof)
#'
#' # The same from a data frame, with the columns named
#' dat <- data.frame(school = letters[1:6], y = c(0.3, 0.1, 0.25, -0.05, 0.4, 0),
#'                   s = se)
#' reliability_profile(dat, se = s, estimate = y, site = school, sigma_tau = 0.2)
#'
#' # The profile at another sigma_tau
#' update(prof, sigma_tau = 0.1)
#' @order 1
#' @export
reliability_profile <- function(x = NULL, ...) {
  # Dispatch on x itself: without it, UseMethod() would evaluate the first
  # argument given, which can be a column name of `data`.
  UseMethod("reliability_profile", x)
}

#' @rdname reliability_profile
#' @order 2
#' @export
reliability_profile.default <- function(x = NULL, ..., se = NULL, se2 = NULL, sigma_tau,
                                        site = NULL, estimate = NULL, data = NULL,
                                        scale_label = NULL, target_label = NULL) {
  if (!is.null(x)) {
    if (is.numeric(x)) {
      .rp_abort(
        "An unnamed vector of standard errors is ambiguous.",
        info = "It could hold standard errors or sampling variances.",
        fix = "Name it: `se = ` for standard errors or `se2 = ` for sampling variances."
      )
    }
    .rp_abort(
      sprintf("There is no method of `reliability_profile()` for %s.", .rp_describe(x)),
      info = "Give the standard errors as `se` or `se2`, or a data frame with those columns."
    )
  }
  .rp_check_dots_empty(...)
  if (missing(sigma_tau)) .rp_abort_missing_sigma_tau()
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      .rp_abort("`data` must be a data frame.", info = .rp_you_passed(data))
    }
    cols <- .rp_eval_columns(data, list(
      se = rlang::enquo(se), se2 = rlang::enquo(se2), site = rlang::enquo(site),
      estimate = rlang::enquo(estimate)
    ))
    se <- cols$se
    se2 <- cols$se2
    site <- cols$site
    estimate <- cols$estimate
  }
  .rp_new_profile(
    se = se, se2 = se2, sigma_tau = sigma_tau, site = site, estimate = estimate,
    source = "supplied", scale_label = scale_label, target_label = target_label,
    call = match.call()
  )
}

#' @rdname reliability_profile
#' @order 3
#' @export
reliability_profile.data.frame <- function(x, ..., se = NULL, se2 = NULL, sigma_tau,
                                           site = NULL, estimate = NULL, scale_label = NULL,
                                           target_label = NULL) {
  .rp_check_dots_empty(...)
  if (missing(sigma_tau)) .rp_abort_missing_sigma_tau()
  cols <- .rp_eval_columns(x, list(
    se = rlang::enquo(se), se2 = rlang::enquo(se2), site = rlang::enquo(site),
    estimate = rlang::enquo(estimate)
  ))
  .rp_new_profile(
    se = cols$se, se2 = cols$se2, sigma_tau = sigma_tau, site = cols$site,
    estimate = cols$estimate, source = "supplied", scale_label = scale_label,
    target_label = target_label, call = match.call()
  )
}

# Builds a reliability_profile from checked or unchecked parts.
.rp_new_profile <- function(se, se2, sigma_tau, site, estimate, source, scale_label = NULL,
                            target_label = NULL, call = NULL, X = NULL) {
  v <- .rp_variances(se, se2)
  J <- length(v)
  sigma_tau <- .rp_check_sigma_tau(sigma_tau)
  site <- .rp_site_labels(.rp_check_site(site, J), J, names(v))
  if (!is.null(estimate)) {
    estimate <- .rp_check_vector(estimate, "estimate")
    if (length(estimate) != J) {
      .rp_abort(
        sprintf("`estimate` must hold one number for each of the %d sites.", J),
        info = .rp_you_passed(estimate)
      )
    }
  }
  labels <- list(
    scale_label = .rp_check_string(scale_label, "scale_label"),
    target_label = .rp_check_string(target_label, "target_label")
  )
  if (sigma_tau == 0) {
    .rp_inform(
      "At `sigma_tau = 0` the profile is a limit: every reliability is 0.",
      info = "The concentration count is its limit as sigma_tau goes to 0.",
      subclass = "rp_message_sigma_zero"
    )
  }
  fixed <- .rp_fixed_profile(unname(v), sigma_tau)
  out <- list(
    se2 = unname(v), site = site, estimate = if (!is.null(estimate)) unname(estimate),
    sigma_tau = sigma_tau, source = source, X = X,
    S = fixed$S, stats = fixed$stats, local = fixed$local,
    labels = labels, call = call, version = .rp_version()
  )
  .rp_validate_profile(structure(out, class = "reliability_profile"))
}

# The reliabilities and diagnostics of sampling variances v at sigma_tau,
# with the limits at sigma_tau = 0.
.rp_fixed_profile <- function(v, sigma_tau) {
  J <- length(v)
  centered <- .rp_center_v(v)
  d <- centered$d
  moments <- .rp_logit_moments(d)
  v_typ <- if (J >= 2L) typical_variance(se2 = v) else NA_real_
  if (sigma_tau == 0) {
    r <- min(v) / v
    stats <- c(I = 0, mean = 0, harmonic = 0, H = 0, V = moments[["V"]],
               rho3 = moments[["rho3"]], J_eff = sum(r)^2 / sum(r^2), gap = 0,
               benchmark_risk = 0, trace = 0, v_typ = v_typ, R_typ = if (J >= 2L) 0 else NA_real_)
    m <- -Inf
    S <- rep(0, J)
    local <- data.frame(quantity = c("H", "gap"), exact = c(0, 0), leading = c(0, 0),
                        bound = c(sqrt(3) / 108, 1 / 48) * moments[["rho3"]])
  } else {
    m <- 2 * log(sigma_tau) - centered$log_gm
    S <- .rp_plogis(m + d)
    stats <- .rp_profile_stats(m, d)
    stats <- c(stats,
               benchmark_risk = mean(.rp_parallel_variances(v, sigma_tau)),
               trace = sum(S),
               v_typ = v_typ,
               R_typ = if (J >= 2L) .rp_plogis(2 * log(sigma_tau) - log(v_typ)) else NA_real_)
    lt <- .rp_local_at(m, moments[["V"]], moments[["rho3"]])
    local <- data.frame(quantity = c("H", "gap"), exact = stats[c("H", "gap")],
                        leading = c(lt$H_leading, lt$gap_leading),
                        bound = c(lt$H_bound, lt$gap_bound))
  }
  local$holds <- .rp_holds(local$leading, local$bound)
  rownames(local) <- NULL
  list(S = S, stats = stats, local = local, m = m, d = d)
}

.rp_validate_profile <- function(x) {
  parts <- c("se2", "site", "sigma_tau", "source", "S", "stats", "local", "labels")
  missing_parts <- setdiff(parts, names(x))
  if (length(missing_parts) > 0L) {
    .rp_abort(
      "This is not a complete `reliability_profile` object.",
      info = sprintf("Missing: %s.", paste(missing_parts, collapse = ", "))
    )
  }
  J <- length(x$se2)
  if (length(x$S) != J || length(x$site) != J ||
      (!is.null(x$estimate) && length(x$estimate) != J) ||
      (!is.null(x$X) && nrow(x$X) != J)) {
    .rp_abort("The parts of this `reliability_profile` object do not have one value per site.")
  }
  x
}

#' @rdname reliability_profile
#' @order 4
#' @param object A profile made by `reliability_profile()`.
#' @export
update.reliability_profile <- function(object, ..., sigma_tau = NULL, scale_label = NULL,
                                       target_label = NULL) {
  .rp_check_dots_empty(...)
  labels <- object$labels
  .rp_new_profile(
    se = NULL, se2 = object$se2,
    sigma_tau = if (is.null(sigma_tau)) object$sigma_tau else sigma_tau,
    site = object$site, estimate = object$estimate, X = object$X,
    source = if (is.null(sigma_tau)) object$source else "supplied",
    scale_label = if (is.null(scale_label)) labels$scale_label else scale_label,
    target_label = if (is.null(target_label)) labels$target_label else target_label,
    call = object$call
  )
}

# Columns of a data frame named in the call, as for metafor's `data`: a name
# is looked up among the columns first and then in the calling environment;
# `.env$name` reaches past the columns.
.rp_eval_columns <- function(data, quos) {
  lapply(stats::setNames(names(quos), names(quos)), function(arg) {
    quo <- quos[[arg]]
    if (rlang::quo_is_null(quo)) {
      return(NULL)
    }
    tryCatch(
      rlang::eval_tidy(quo, data = data),
      error = function(err) {
        if (grepl("not found", conditionMessage(err), fixed = TRUE)) {
          .rp_abort(
            sprintf("`%s` could not be found in the data.", arg),
            info = sprintf("`%s = %s` is neither a column of the data nor an object.",
                           arg, rlang::as_label(quo)),
            fix = sprintf("The columns are %s.", .rp_column_list(names(data))),
            parent = err
          )
        }
        .rp_abort(sprintf("`%s` could not be evaluated in the data.", arg), parent = err)
      }
    )
  })
}

.rp_column_list <- function(columns, max_shown = 8L) {
  shown <- sprintf("`%s`", utils::head(columns, max_shown))
  if (length(columns) > max_shown) {
    return(paste0(paste(shown, collapse = ", "), sprintf(" and %d more", length(columns) - max_shown)))
  }
  paste(shown, collapse = ", ")
}

.rp_version <- function() {
  as.character(utils::packageVersion("poolinfo"))
}
