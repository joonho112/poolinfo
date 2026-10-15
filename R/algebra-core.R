# The diagnostics as functions of m, the logit of the informativeness, for
# a fixed vector d of centered logits. Each function takes a vector m, one
# value for each point of a path or of a rule of integration, and returns one
# value for each element of m. The functions for users call them with a
# single m; the posterior engine calls them with its nodes.
#
# With m and d, the logit of site j's reliability is m + d_j. Moving m with
# d fixed is what a change of sigma_tau does when the sampling variances are
# fixed: m = 2 log(sigma_tau) - mean(log v).
#
# In floating point, the deviations d sum to a small number rather than to
# 0, of the order of the rounding of m. The gap mean(S) - I is first-order in
# that number, so it is corrected for it; the other summaries are not
# sensitive to it.

# The logistic function, without the underflow of plogis() below -709.78,
# where R returns 0 rather than the subnormal value.
.rp_plogis <- function(x) {
  out <- stats::plogis(x)
  low <- !is.na(x) & x < -700
  out[low] <- exp(stats::plogis(x[low], log.p = TRUE))
  out
}

# s(m) (1 - s(m)), without underflow.
.rp_dlogis <- function(m) {
  exp(stats::plogis(m, log.p = TRUE) + stats::plogis(-m, log.p = TRUE))
}

# KL(Ber(s(m)) || Ber(s(m + delta))), elementwise, where s is the logistic
# function: the Bregman divergence A(m + delta) - A(m) - s(m) delta of the
# softplus A. The base probability p is reflected to at most 1/2. For
# |delta| <= 1e-3 the Taylor series of A to the fifth order is used, whose
# relative error is about delta^4 / 360; up to delta = 700,
# log1p(p expm1(delta)) - p delta, which cannot overflow because p <= 1/2;
# beyond, the logarithm of p expm1(delta) is formed first.
.rp_bern_kl <- function(m, delta) {
  out <- rep(NaN, length(m))
  ok <- !is.na(m) & !is.na(delta)
  delta <- delta[ok] * (1 - 2 * (m[ok] > 0))
  m <- -abs(m[ok])
  p <- .rp_plogis(m)
  r <- numeric(length(m))
  small <- abs(delta) <= 1e-3
  if (any(small)) {
    ps <- p[small]
    ds <- delta[small]
    q <- 1 - 2 * ps
    r[small] <- ps * (1 - ps) * ds^2 *
      (1 / 2 + ds * (q / 6 + ds * ((1 - 6 * ps + 6 * ps^2) / 24 +
                                     ds * q * (1 - 12 * ps + 12 * ps^2) / 120)))
  }
  mid <- !small & delta <= 700
  r[mid] <- log1p(p[mid] * expm1(delta[mid])) - p[mid] * delta[mid]
  big <- !small & delta > 700
  if (any(big)) {
    z <- stats::plogis(m[big], log.p = TRUE) + delta[big] + log(-expm1(-delta[big]))
    r[big] <- .rp_softplus(z) - p[big] * delta[big]
  }
  r[is.infinite(delta)] <- Inf
  out[ok] <- pmax(r, 0)
  out
}

# The terms S_j - I of the gap, elementwise, from
# s(m + d) - s(m) = -s(m + d) (1 - s(m)) expm1(-d), which does not subtract
# two numbers close to 1. Where expm1(-d) would overflow (d < -700) the term
# is formed on the log scale.
.rp_gap_terms <- function(m, d) {
  out <- numeric(length(m))
  far <- !is.na(d) & d < -700
  near <- !far
  out[near] <- -.rp_plogis(-m[near]) * .rp_plogis(m[near] + d[near]) * expm1(-d[near])
  if (any(far)) {
    out[far] <- -exp(stats::plogis(-m[far], log.p = TRUE) +
                       stats::plogis(m[far] + d[far], log.p = TRUE) - d[far] + log(-expm1(d[far])))
  }
  out
}

# Applies `fun(eta, d, m, n)` to blocks of the rows of the matrix of logits
# m_i + d_j, so that no block has more than `max_cells` cells. `fun` gets
# the logits, the deviations and m as vectors in column-major order, and n,
# the number of rows of the block; it returns one value for each row.
.rp_by_rows <- function(m, d, fun, max_cells = 5e5) {
  J <- length(d)
  out <- numeric(length(m))
  if (length(m) == 0L) {
    return(out)
  }
  rows <- max(1L, floor(max_cells / J))
  for (start in seq.int(1L, length(m), by = rows)) {
    at <- start:min(length(m), start + rows - 1L)
    n <- length(at)
    mm <- rep(m[at], times = J)
    dd <- rep(d, each = n)
    out[at] <- fun(mm + dd, dd, mm, n)
  }
  out
}

.rp_row_means <- function(x, n) {
  rowMeans(matrix(x, nrow = n))
}

# Mean reliability.
.rp_mean_at <- function(m, d) {
  .rp_by_rows(m, d, function(eta, dd, mm, n) .rp_row_means(.rp_plogis(eta), n))
}

# Mean reliability minus I, from the terms above, less s'(m) mean(d), the
# first-order effect of the deviations not summing exactly to 0.
.rp_gap_at <- function(m, d) {
  .rp_by_rows(m, d, function(eta, dd, mm, n) .rp_row_means(.rp_gap_terms(mm, dd), n)) -
    .rp_dlogis(m) * mean(d)
}

# Dispersion H, the mean of the site terms KL(Ber(I) || Ber(S_j)).
.rp_H_at <- function(m, d) {
  .rp_by_rows(m, d, function(eta, dd, mm, n) .rp_row_means(.rp_bern_kl(mm, dd), n))
}

# Concentration count (sum S)^2 / sum S^2, with the reliabilities divided by
# the largest one on the log scale, so that it is defined when they are all
# too small to be represented. The largest is that of the largest d. At
# m = -Inf it is the limit (sum e^d)^2 / sum e^2d.
.rp_jeff_at <- function(m, d) {
  top <- max(d)
  out <- .rp_by_rows(m, d, function(eta, dd, mm, n) {
    r <- matrix(exp(stats::plogis(eta, log.p = TRUE) - stats::plogis(mm + top, log.p = TRUE)),
                nrow = n)
    rowSums(r)^2 / rowSums(r^2)
  })
  low <- !is.na(m) & m == -Inf
  if (any(low)) {
    r <- exp(d - top)
    out[low] <- sum(r)^2 / sum(r^2)
  }
  out
}

# Harmonic centroid 1 / mean(1 / S_j) = s(m - c) with c = log(mean(exp(-d))).
.rp_harmonic_at <- function(m, d) {
  .rp_plogis(m - .rp_log_mean_exp(-d))
}

.rp_log_mean_exp <- function(x) {
  top <- max(x)
  top + log(mean(exp(x - top)))
}

# The logit moments of the deviations: V and rho3.
.rp_logit_moments <- function(d) {
  c(V = mean(d^2), rho3 = mean(abs(d)^3))
}

# The leading terms of Proposition 2 and the bounds on their remainders.
.rp_local_at <- function(m, V, rho3) {
  i1i <- .rp_dlogis(m)
  list(
    H_leading = 0.5 * i1i * V,
    gap_leading = 0.5 * i1i * (-tanh(m / 2)) * V,
    H_bound = sqrt(3) / 108 * rho3,
    gap_bound = rho3 / 48
  )
}

# Whether a leading term decides the sign of the exact value: it is larger
# than the bound on the remainder. NA where both are 0, as for a profile whose
# reliabilities are all equal.
.rp_holds <- function(leading, bound) {
  out <- abs(leading) > bound
  out[leading == 0 & bound == 0] <- NA
  out
}

# All the summaries of one profile, from m and d. The harmonic centroid is
# at most both I and the mean reliability; rounding can put it a few units
# in the last place above them, which the minimum removes.
.rp_profile_stats <- function(m, d) {
  moments <- .rp_logit_moments(d)
  I <- .rp_plogis(m)
  mean_s <- .rp_mean_at(m, d)
  c(
    I = I,
    mean = mean_s,
    harmonic = min(.rp_harmonic_at(m, d), I, mean_s),
    H = .rp_H_at(m, d),
    V = moments[["V"]],
    rho3 = moments[["rho3"]],
    J_eff = .rp_jeff_at(m, d),
    gap = .rp_gap_at(m, d)
  )
}
