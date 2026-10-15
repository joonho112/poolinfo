# Quantiles of the diagnostics. The informativeness, the mean reliability,
# the harmonic centroid and the concentration count increase with sigma_tau,
# so their quantiles are those of sigma_tau mapped through them. The
# dispersion H does not: it rises from 0, reaches a largest value where the
# mean reliability equals I, and falls back to 0. Its distribution is the sum
# over the branches on which it is monotone (Appendix D.3 of the paper).

# The turning point of H in m = logit(I): the zero of the gap mean(S) - I.
# Along a path with the deviations d fixed, d^2 H / dI^2 =
# -mean((S_j - I)^2) / (I (1 - I))^2 < 0, so H is strictly concave in I and
# the gap, which has the sign of dH/dI, changes sign exactly once: it is
# positive as m goes to -Inf and negative as m goes to Inf. The root is
# found in a bracket where the gap's sign is known, and two scans at two
# resolutions check that each finds exactly that one change of sign.
.rp_turning_points <- function(d) {
  V <- mean(d^2)
  if (V < .RP_V_HOMOGENEOUS) {
    return(list(roots = numeric(), homogeneous = TRUE, agree = TRUE, difference = 0))
  }
  lo <- min(-40, -max(d) - 10)
  hi <- max(40, -min(d) + 10)
  gap <- function(m) .rp_gap_at(m, d)
  root <- stats::uniroot(gap, c(lo, hi), tol = 1e-13)$root
  scan <- function(n) {
    m <- seq(lo, hi, length.out = n)
    g <- gap(m)
    at <- which(g[-n] * g[-1L] < 0 | g[-n] == 0)
    vapply(at, function(i) {
      if (g[[i]] == 0) return(m[[i]])
      stats::uniroot(gap, m[c(i, i + 1L)], tol = 1e-13)$root
    }, numeric(1))
  }
  n1 <- max(161L, ceiling(2 * (hi - lo)) + 1L)
  roots1 <- scan(n1)
  roots2 <- scan(2L * n1 - 1L)
  agree <- length(roots1) == 1L && length(roots2) == 1L
  difference <- if (agree) max(abs(c(roots1, roots2) - root)) else Inf
  list(roots = root, homogeneous = FALSE, agree = agree && difference < 1e-7,
       difference = difference)
}

# The branches of H as a function of sigma: their ends in m (with -Inf and
# Inf for sigma = 0 and sigma = Inf), the values of H at the ends, and
# whether H rises on each.
.rp_H_branches <- function(kernel, turning) {
  ends <- c(-Inf, turning$roots, Inf)
  h <- c(0, .rp_H_at(turning$roots, kernel$d), 0)
  list(m = ends, h = h, rising = diff(h) > 0, H_max = max(h))
}

# The m on a monotone branch where H equals h: Newton's method with the exact
# derivative dH/dm = mean(S) - I, kept inside a bracket that bisection
# shrinks whenever a Newton step would leave it. The outer branches are
# bracketed by stepping away from the turning point.
.rp_H_inverse <- function(d, h, m_a, m_b, h_a, h_b) {
  f <- function(m) .rp_H_at(m, d) - h
  if (is.infinite(m_a)) {
    step <- 1
    m_a <- m_b - step
    while (f(m_a) > 0 && step < 1e5) {
      step <- 2 * step
      m_a <- m_b - step
    }
  }
  if (is.infinite(m_b)) {
    step <- 1
    m_b <- m_a + step
    while (f(m_b) > 0 && step < 1e5) {
      step <- 2 * step
      m_b <- m_a + step
    }
  }
  .rp_newton(f, function(m) .rp_gap_at(m, d), m_a, m_b, tol = 1e-13)
}

# Safeguarded Newton's method for a root of f in [a, b], where f changes
# sign; `df` is the derivative.
.rp_newton <- function(f, df, a, b, tol = 1e-13, max_iter = 100L) {
  fa <- f(a)
  fb <- f(b)
  if (fa == 0) return(a)
  if (fb == 0) return(b)
  if (sign(fa) == sign(fb)) {
    return(stats::uniroot(f, c(a, b), tol = tol, extendInt = "yes")$root)
  }
  x <- (a + b) / 2
  for (i in seq_len(max_iter)) {
    fx <- f(x)
    if (fx == 0) return(x)
    if (sign(fx) == sign(fa)) {
      a <- x
      fa <- fx
    } else {
      b <- x
    }
    step <- fx / df(x)
    nx <- x - step
    if (!is.finite(nx) || nx <= min(a, b) || nx >= max(a, b)) {
      nx <- (a + b) / 2
    }
    if (abs(nx - x) <= tol * max(1, abs(x))) return(nx)
    x <- nx
  }
  x
}

.rp_sigma_of_m <- function(kernel, m) {
  exp((m + kernel$log_gm) / 2)
}

# P(H <= h) for a vector h.
.rp_H_cdf <- function(law, kernel, branches, h) {
  vapply(h, function(x) {
    if (x <= 0) return(0)
    if (x >= branches$H_max) return(1)
    total <- 0
    for (k in seq_along(branches$rising)) {
      m_a <- branches$m[[k]]
      m_b <- branches$m[[k + 1L]]
      h_a <- branches$h[[k]]
      h_b <- branches$h[[k + 1L]]
      s_a <- if (is.infinite(m_a)) 0 else .rp_sigma_of_m(kernel, m_a)
      s_b <- if (is.infinite(m_b)) Inf else .rp_sigma_of_m(kernel, m_b)
      F_a <- if (s_a == 0) 0 else law$cdf(s_a)
      F_b <- if (is.infinite(s_b)) 1 else law$cdf(s_b)
      if (x >= max(h_a, h_b)) {
        total <- total + F_b - F_a
      } else if (x > min(h_a, h_b)) {
        root <- .rp_sigma_of_m(kernel, .rp_H_inverse(kernel$d, x, m_a, m_b, h_a, h_b))
        F_r <- law$cdf(root)
        total <- total + if (branches$rising[[k]]) F_r - F_a else F_b - F_r
      }
    }
    min(1, max(0, total))
  }, numeric(1))
}

.rp_H_quantile <- function(law, kernel, branches, p) {
  vapply(p, function(pp) {
    if (branches$H_max == 0 || pp <= 0) return(0)
    if (pp >= 1) return(branches$H_max)
    stats::uniroot(function(h) .rp_H_cdf(law, kernel, branches, h) - pp,
                   c(0, branches$H_max), tol = 1e-13 * branches$H_max)$root
  }, numeric(1))
}
