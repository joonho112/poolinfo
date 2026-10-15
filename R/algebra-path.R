#' Dispersion and centroid gap along a path of informativeness
#'
#' Keeps the centered logits of a reliability profile fixed and moves its
#' informativeness \eqn{I}, as a change of the between-site standard
#' deviation does when the sampling variances stay fixed. Returns the
#' dispersion, the mean reliability and their second-order terms at each
#' point of the path.
#'
#' @param x Site reliabilities, a numeric vector of numbers strictly between
#'   0 and 1, or an object of class `reliability_profile`,
#'   `reliability_posterior`, `reliability_draws` or `reliability_design`.
#' @param I Values of the informativeness at which to evaluate the path, a
#'   numeric vector of numbers strictly between 0 and 1. By default `n`
#'   equally spaced values.
#' @param n Number of points of the default path, an integer of at least 2.
#'
#' @details
#' Write \eqn{m} for the logit of \eqn{I} and \eqn{\Delta_j}{Delta_j} for the
#' centered logits of the profile. Along the path \eqn{\Delta_j}{Delta_j}
#' stay fixed and \eqn{m} moves, and Proposition 2 of Lee (2026) gives
#' \eqn{\partial H / \partial m = \bar S - I}{dH/dm = mean(S_j) - I}. So the
#' dispersion rises from 0 while the mean reliability is above \eqn{I},
#' falls back to 0 after it, and is largest where the gap changes sign.
#'
#' For an object, the path also gives the between-site standard deviation
#' at each point,
#' \eqn{\sigma_\tau = \{\mathrm{GM}(v)\, I/(1-I)\}^{1/2}}{sigma_tau = sqrt(GM(v) I / (1 - I))},
#' where \eqn{\mathrm{GM}(v)}{GM(v)} is the geometric mean of the sampling
#' variances. For reliabilities alone the scale is unknown, and that column
#' is absent.
#'
#' @return A data frame with one row for each value of `I` and the columns
#'   * `I` and `logit_I`.
#'   * `sigma_tau`, for an object.
#'   * `H`, `mean` and `gap`: the dispersion, the mean reliability and the
#'     mean reliability minus `I`.
#'   * `H_leading`, `gap_leading`, `H_bound` and `gap_bound`: the
#'     second-order terms and the bounds on the remainders, as in
#'     [local_terms()].
#'
#' @template ref-lee2026
#' @examples
#' S <- site_reliability(se = c(0.10, 0.15, 0.20, 0.30), sigma_tau = 0.2)
#' path <- profile_path(S, n = 9)
#' path[, c("I", "H", "gap")]
#' @export
profile_path <- function(x, I = NULL, n = 201) {
  source <- .rp_path_source(x)
  if (is.null(I)) {
    n <- .rp_check_whole(n, "n", lower = 2)
    I <- seq_len(n) / (n + 1)
  } else {
    I <- .rp_check_vector(I, "I")
    if (any(I <= 0 | I >= 1)) {
      .rp_abort(
        "`I` must hold numbers strictly between 0 and 1.",
        info = .rp_you_passed(I)
      )
    }
  }
  m <- stats::qlogis(I)
  d <- source$d
  moments <- .rp_logit_moments(d)
  local <- .rp_local_at(m, moments[["V"]], moments[["rho3"]])
  out <- data.frame(I = I, logit_I = m)
  if (!is.null(source$log_gm)) {
    out$sigma_tau <- exp((m + source$log_gm) / 2)
  }
  out$H <- .rp_H_at(m, d)
  out$mean <- .rp_mean_at(m, d)
  out$gap <- .rp_gap_at(m, d)
  out$H_leading <- local$H_leading
  out$gap_leading <- local$gap_leading
  out$H_bound <- rep(local$H_bound, length(m))
  out$gap_bound <- rep(local$gap_bound, length(m))
  out
}

# The centered logits of what profile_path() is given, and the log of the
# geometric mean of the sampling variances when they are known.
.rp_path_source <- function(x) {
  v <- .rp_object_se2(x)
  if (!is.null(v)) {
    centered <- .rp_center_v(v)
    return(list(d = centered$d, log_gm = centered$log_gm))
  }
  if (is.numeric(x) && !is.object(x)) {
    S <- .rp_check_reliabilities(x, "x")
    return(list(d = .rp_center_S(S)$d, log_gm = NULL))
  }
  .rp_abort(
    "`x` must be site reliabilities, or a profile, a posterior fit, draws or a design.",
    info = .rp_you_passed(x)
  )
}

# The sampling variances of an object of the package, checked, or NULL for
# anything else.
.rp_object_se2 <- function(x) {
  classes <- c("reliability_profile", "reliability_posterior", "reliability_draws",
               "reliability_design")
  if (!inherits(x, classes)) {
    return(NULL)
  }
  if (!is.list(x) || is.null(x$se2)) {
    .rp_abort(sprintf("This %s object has no sampling variances.", class(x)[[1L]]),
              info = "It may have been changed after it was made; make it again.")
  }
  .rp_check_vector(x$se2, "the object's sampling variances", positive = TRUE)
}
