#' Confidence interval for the between-site variance mapped through the diagnostics
#'
#' Maps a confidence interval for the between-site variance, such as the
#' Q-profile interval of a metafor fit, through the informativeness, the mean
#' reliability, the concentration count and the dispersion. In metafor's
#' notation the between-site variance is \eqn{\tau^2}{tau^2}.
#'
#' @param x An optional `rma.uni` fit of the metafor package, from which the
#'   sampling variances, the estimate of \eqn{\tau^2}{tau^2} and its
#'   Q-profile confidence interval are computed.
#' @param se,se2 Standard errors or sampling variances of the site
#'   estimates, when `x` is not given. Give one of them.
#' @param interval The confidence interval for the between-site variance, two
#'   numbers of at least 0, optionally with the point estimate, as in
#'   `c(estimate = 0.02, lower = 0.005, upper = 0.08)`. It is needed without
#'   `x`; with `x` it replaces the computed Q-profile interval.
#' @param level The level of the Q-profile interval computed for `x`,
#'   a number between 0 and 1, when `interval` is not given.
#'
#' @details
#' Confidence-interval mapping is an experimental extension beyond Lee
#' (2026); its interface may change in a later 0.x release.
#'
#' The informativeness and the mean reliability increase with
#' \eqn{\tau^2}{tau^2}, so the interval between their values at the ends of
#' a confidence interval \eqn{C} has the coverage of \eqn{C}. So does the
#' interval of the concentration count when the sampling variances differ;
#' when they are all equal the count is \eqn{J} whatever \eqn{\tau^2}{tau^2}
#' is, and its interval is the single value \eqn{J}. The dispersion
#' \eqn{H} first rises
#' and then falls; its interval is the image \eqn{[\min_C H, \max_C H]}{[min H over C, max H over C]},
#' which covers the true value at least as often as \eqn{C} covers
#' \eqn{\tau^2}{tau^2}. These are frequentist intervals, conditional on the
#' sampling variances like the rest of the package; they are not posterior
#' intervals.
#'
#' For a metafor fit the interval is the Q-profile interval for the model
#' of the fit, moderators included. poolinfo computes it for every
#' estimator, also for `GENQ` and `GENQM` fits, whose own `confint()` uses
#' another method, and the result does not depend on the units of the
#' effects. A fit without a between-site variance (`method = "FE"`, `"EE"`
#' or `"CE"`) is an error.
#'
#' @return A data frame with one row for each of `"I"`, `"mean"`, `"J_eff"`
#'   and `"H"`, and the columns `functional`, `estimate` (at the point
#'   estimate of \eqn{\tau^2}{tau^2}, when it is known), `lower` and `upper`.
#'   Its attribute `interval` is the interval for \eqn{\tau^2}{tau^2} that
#'   was mapped, and its attribute `method` describes where that interval
#'   came from.
#'
#' @seealso [metafor-methods] for profiles and posteriors of metafor fits,
#'   and `vignette("meta-analysis")` for an example.
#'
#' @examples
#' map_tau2_interval(se = eight_schools$se, interval = c(estimate = 20, lower = 0, upper = 500))
#' @export
map_tau2_interval <- function(x = NULL, se = NULL, se2 = NULL, interval = NULL, level = 0.95) {
  if (!is.null(x)) {
    if (!inherits(x, "rma.uni")) {
      .rp_abort("`x` must be an `rma.uni` fit of the metafor package, or NULL.",
                info = .rp_you_passed(x))
    }
    if (!requireNamespace("metafor", quietly = TRUE)) {
      .rp_abort("Mapping a metafor fit needs the metafor package.")
    }
    .rp_check_rma_scope(x)
    v <- .rp_check_vector(x$vi, "the sampling variances of the fit", positive = TRUE)
    if (is.null(interval)) {
      random <- .rp_rma_tau2_interval(x, .rp_check_level(level))
      interval <- c(estimate = x$tau2, lower = random[["ci.lb"]], upper = random[["ci.ub"]])
      method <- sprintf("%s%% Q-profile confidence interval for tau^2 of the %s fit", format(100 * level), x$method)
    } else {
      method <- "the given confidence interval for tau^2"
    }
  } else {
    v <- .rp_variances(se, se2)
    method <- "the given confidence interval for tau^2"
  }
  interval <- .rp_check_tau2_interval(interval)
  centered <- .rp_center_v(v)
  d <- centered$d
  at <- function(tau2) {
    if (tau2 == 0) {
      r <- exp(d - max(d))
      return(c(I = 0, mean = 0, J_eff = sum(r)^2 / sum(r^2), H = 0))
    }
    st <- .rp_profile_stats(log(tau2) - centered$log_gm, d)
    st[c("I", "mean", "J_eff", "H")]
  }
  lo <- at(interval[["lower"]])
  hi <- at(interval[["upper"]])
  est <- if (!is.na(interval[["estimate"]])) at(interval[["estimate"]]) else rep(NA_real_, 4L)
  # The largest H over the interval: at the turning point if it is inside.
  turning <- .rp_turning_points(d)
  H_hi <- max(lo[["H"]], hi[["H"]])
  if (!turning$homogeneous) {
    tau2_turn <- exp(turning$roots + centered$log_gm)
    if (tau2_turn > interval[["lower"]] && tau2_turn < interval[["upper"]]) {
      H_hi <- .rp_H_at(turning$roots, d)
    }
  }
  out <- data.frame(
    functional = c("I", "mean", "J_eff", "H"),
    estimate = unname(est),
    lower = c(lo[["I"]], lo[["mean"]], lo[["J_eff"]], min(lo[["H"]], hi[["H"]])),
    upper = c(hi[["I"]], hi[["mean"]], hi[["J_eff"]], H_hi)
  )
  attr(out, "interval") <- interval
  attr(out, "method") <- method
  out
}

# Solve generalized-Q endpoints after removing the units of the variance.
.rp_rma_tau2_interval <- function(x, level) {
  v <- as.numeric(x$vi)
  scale2 <- exp(mean(log(v)))
  v <- v / scale2
  X <- as.matrix(x$X)
  y <- as.numeric(x$yi) / sqrt(scale2)
  df <- length(y) - ncol(X)
  if (df <= 0L || !is.finite(scale2) || scale2 <= 0) {
    .rp_abort("The Q-profile interval needs positive residual degrees of freedom and a finite variance scale.",
              class = "rp_scope_error")
  }
  Q <- function(tau2) {
    inv_sd <- 1 / sqrt(v + tau2)
    residual <- qr.resid(qr(X * inv_sd), y * inv_sd)
    sum(residual^2)
  }
  at_zero <- Q(0)
  critical <- stats::qchisq(c((1 + level) / 2, (1 - level) / 2), df = df)
  # Q falls from Q(0) towards 0 as tau^2 grows. Each end is bracketed by
  # powers of 10 and found on the log scale, so that an end far below the
  # typical sampling variance has the same relative precision as one far
  # above it.
  endpoint <- function(target) {
    if (at_zero <= target) return(0)
    upper <- max(1, x$tau2 / scale2)
    for (i in seq_len(400L)) {
      if (!is.finite(upper) || Q(upper) <= target) break
      upper <- upper * 10
    }
    lower <- upper / 10
    for (i in seq_len(700L)) {
      if (lower < 1e-300 || Q(lower) > target) break
      lower <- lower / 10
    }
    if (!is.finite(upper) || Q(upper) > target || !(Q(lower) > target)) {
      .rp_abort("The ends of the Q-profile interval could not be bracketed.",
                fix = "Give the interval with `interval`.", class = "rp_numeric_error")
    }
    exp(stats::uniroot(function(t) Q(exp(t)) - target, log(c(lower, upper)),
                       tol = 1e-12, maxiter = 1000L)$root)
  }
  stats::setNames(vapply(critical, endpoint, 0) * scale2, c("ci.lb", "ci.ub"))
}

.rp_check_tau2_interval <- function(interval) {
  if (is.null(interval)) {
    .rp_abort("`interval` is missing.", info = "Give the confidence interval for tau^2, or a metafor fit.")
  }
  if (!is.numeric(interval) || !length(interval) %in% c(2L, 3L) || anyNA(interval) ||
      any(interval < 0) || any(!is.finite(interval))) {
    .rp_abort("`interval` must hold two or three finite numbers of at least 0.",
              info = .rp_you_passed(interval))
  }
  if (length(interval) == 2L) {
    if (is.null(names(interval))) {
      return(c(estimate = NA_real_, lower = min(interval), upper = max(interval)))
    }
    if (!setequal(names(interval), c("lower", "upper"))) {
      .rp_abort("Two named values of `interval` must be named `lower` and `upper`.",
                info = .rp_you_passed(interval))
    }
    out <- c(estimate = NA_real_, lower = interval[["lower"]], upper = interval[["upper"]])
  } else {
    if (is.null(names(interval)) || !setequal(names(interval), c("estimate", "lower", "upper"))) {
      .rp_abort("Three values of `interval` must be named `estimate`, `lower` and `upper`.",
                info = .rp_you_passed(interval))
    }
    out <- interval[c("estimate", "lower", "upper")]
  }
  if (out[["lower"]] > out[["upper"]]) {
    .rp_abort("The lower end of `interval` is above its upper end.")
  }
  out
}
