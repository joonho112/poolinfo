# The adaptive engine: the computation of the paper (Appendix D.3) with
# stats::integrate(). The density is divided by its largest value, and each
# integral over sigma is split at sqrt(GM(v)), as in the paper, because the
# first rule that integrate() applies to a long range can miss a posterior
# concentrated in a small part of it. It is also split across the stretch
# where the posterior is: a posterior concentrated far above sqrt(GM(v)), as
# with many sites and an informativeness close to 1, is otherwise missed the
# same way. Every
# expectation is a separate integral, which makes the engine slow for many
# sites; it is kept as the reference against which the node rule is
# checked. The second moments about a mean, for the standard deviations,
# come from a node rule (see .rp_law()).
#
# The integrals run over u = sigma / c, with c a rough value of the
# integral of the density divided by its largest value (its width), so that
# the integral of the density over u is near 1; and each integrand is divided
# by its largest value on a set of points. The absolute tolerance of
# integrate() is then relative to the posterior and to the function,
# whatever the units of the estimates. The integrals are kept in u (the
# factor c cancels from every ratio), and they end at the upper end of the
# node rule, beyond which the bounds on the omitted tail meet their limits.

.rp_adaptive_law <- function(kernel, rel_tol = 1e-10, abs_tol = 1e-12) {
  B <- .rp_upper_end(kernel)
  # Where the posterior is, on the log scale: the stretch of t = log(sigma)
  # where the density in t is above exp(-70) of its largest value, from a
  # wide grid.
  t_grid <- seq(log(.rp_sigma_floor(kernel)) - 22, log(B), length.out = 2001L)
  lk <- .rp_log_kernel(kernel, exp(t_grid))
  lk_t <- lk + t_grid
  inside <- which(lk_t > max(lk_t) - 70)
  window <- t_grid[range(inside)]
  mode <- exp(t_grid[[which.max(lk_t)]])
  shift <- max(.rp_log_kernel(kernel, 0), lk)
  c0 <- sum(exp(lk_t - shift)) * (t_grid[[2L]] - t_grid[[1L]])
  dens <- function(s) {
    out <- exp(.rp_log_kernel(kernel, s) - shift)
    out[!is.finite(out)] <- 0
    out
  }
  probe <- exp(seq(window[[1L]], window[[2L]], length.out = 61L))
  probe_dens <- dens(probe)
  # The largest value of the integrand of `fun` at the probe points, by which
  # the integrand is divided.
  size_of <- function(fun) {
    m <- max(abs(fun(probe)) * probe_dens)
    if (is.finite(m) && m > 0) m else 1
  }
  # integrate() reports a roundoff error where the integrand is nearly 0, or
  # where the integrand's own rounding is above rel_tol, as for the squared
  # deviation of a quantity that hardly varies from its mean. Its result is
  # then taken when its error estimate is below 1e-6 of the value (the
  # normalizing constant is checked against 1e-10 separately). Otherwise, and
  # after any other failure, the stretch is halved (on the log scale above
  # 0) and each half integrated, up to six times.
  quad <- function(fun, lo, hi, size, depth = 0L) {
    if (hi <= lo) {
      return(c(value = 0, error = 0))
    }
    f <- function(u) {
      s <- c0 * u
      d <- dens(s)
      ans <- d * fun(s) / size
      ans[d == 0] <- 0
      ans
    }
    z <- stats::integrate(f, lo / c0, hi / c0, subdivisions = 2000L, rel.tol = rel_tol,
                          abs.tol = abs_tol, stop.on.error = FALSE)
    rounding <- grepl("roundoff|bad integrand", z$message)
    if (identical(z$message, "OK") ||
        (rounding && is.finite(z$abs.error) && z$abs.error <= max(abs_tol, 1e-6 * abs(z$value)))) {
      return(c(value = z$value, error = z$abs.error) * size)
    }
    if (depth < 6L) {
      mid <- if (lo > 0) sqrt(lo * hi) else hi / 2
      return(quad(fun, lo, mid, size, depth + 1L) + quad(fun, mid, hi, size, depth + 1L))
    }
    .rp_abort(
      "The adaptive engine could not integrate the posterior of sigma_tau.",
      info = sprintf("stats::integrate() reported: %s.", z$message),
      fix = "Use `posterior_control(engine = \"nodes\")`.",
      class = "rp_numeric_error", subclass = "rp_integrate_error"
    )
  }
  # The integrals are split at sqrt(GM(v)), as in the paper, and at nine
  # points evenly spaced on the log scale across that stretch.
  cuts <- c(kernel$scale, exp(seq(window[[1L]], window[[2L]], length.out = 9)))
  cuts <- sort(unique(cuts[cuts > 0 & cuts < B]))
  breaks <- c(0, cuts, B)
  n_pieces <- length(breaks) - 1L
  one <- function(s) rep(1, length(s))
  one_size <- size_of(one)
  pieces <- vapply(seq_len(n_pieces), function(i) quad(one, breaks[[i]], breaks[[i + 1L]], one_size),
                   numeric(2))
  Z <- sum(pieces["value", ])
  expect <- function(fun) {
    size <- size_of(fun)
    sum(vapply(seq_len(n_pieces), function(i) {
      quad(fun, breaks[[i]], breaks[[i + 1L]], size)[["value"]]
    }, numeric(1))) / Z
  }
  below <- c(0, cumsum(pieces["value", ]))
  above <- c(rev(cumsum(rev(pieces["value", ]))), 0)
  # The probability above s is summed from the pieces above, so that a small
  # upper tail keeps its digits.
  cdf <- function(s, lower = TRUE) {
    vapply(s, function(x) {
      if (is.na(x)) return(NA_real_)
      if (x <= 0) return(if (lower) 0 else 1)
      if (x >= B) return(if (lower) 1 else 0)
      i <- findInterval(x, breaks)
      part <- if (lower) {
        below[[i]] + quad(one, breaks[[i]], x, one_size)[["value"]]
      } else {
        above[[i + 1L]] + quad(one, x, breaks[[i + 1L]], one_size)[["value"]]
      }
      min(1, max(0, part / Z))
    }, numeric(1))
  }
  # Quantiles: roots of the CDF on the log scale, with a relative tolerance,
  # so that a quantile far below the mode keeps its digits.
  quantile <- function(p) {
    vapply(p, function(pp) {
      if (pp <= 0) return(0)
      if (pp >= 1) return(B)
      lo <- mode
      for (i in seq_len(60L)) {
        if (cdf(lo) < pp) break
        lo <- lo * 1e-4
      }
      hi <- mode
      while (hi < B && cdf(hi) < pp) hi <- min(B, hi * 1e4)
      exp(stats::uniroot(function(t) cdf(exp(t)) - pp, log(c(lo, hi)), tol = 1e-12)$root)
    }, numeric(1))
  }
  list(
    engine = "adaptive", kernel = kernel,
    expect = expect, cdf = cdf, quantile = quantile,
    density = function(s) dens(s) / (Z * c0),
    log_Z = shift + log(Z) + log(c0), upper = B,
    checks = c(normalization_error = sum(pieces["error", ]) / Z)
  )
}
