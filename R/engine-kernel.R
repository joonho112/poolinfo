# The marginal posterior of sigma_tau, up to a constant. With w_j = 1 /
# (sigma^2 + v_j), W = diag(w), the mean of the effects x_j' beta integrated
# out under a flat prior, and a prior p(sigma) on sigma,
#
#   log k(sigma) = -1/2 log|X'WX| - 1/2 sum_j log(sigma^2 + v_j)
#                  - 1/2 (y - X beta_sigma)' W (y - X beta_sigma) + log p(sigma),
#
# with beta_sigma the weighted least-squares fit (Lee, 2026, Appendix D; for
# the mean model, the same integral over beta). Without a mean model X is a
# column of ones, X'WX = W and beta_sigma the weighted mean.

.rp_kernel <- function(y, v, X = NULL, prior = sd_prior()) {
  J <- length(y)
  p <- if (is.null(X)) 1L else ncol(X)
  centered <- .rp_center_v(v)
  # Remove a fixed fitted value in the span of X before any weighted fit.
  # In particular, subtracting y[1] exactly preserves represented contrasts
  # when a large common location dwarfs the standard errors.
  if (is.null(X)) {
    coef_offset <- y[[1L]]
    fitted_offset <- rep(y[[1L]], J)
  } else {
    constant <- which(vapply(seq_len(p), function(j) {
      X[1L, j] != 0 && all(X[, j] == X[1L, j])
    }, logical(1)))
    if (length(constant)) {
      coef_offset <- numeric(p)
      coef_offset[constant[[1L]]] <- y[[1L]] / X[1L, constant[[1L]]]
      fitted_offset <- rep(y[[1L]], J)
    } else {
      # A no-intercept basis can still span a constant (for example,
      # columns 1 + x and 1 - x). Use a row-pivoted square system to find
      # that representation and accept it only when X b equals 1 exactly
      # in floating arithmetic; nearly constant fits must retain their model.
      rows <- qr(t(X), LAPACK = TRUE)$pivot[seq_len(p)]
      constant_coef <- tryCatch(solve(X[rows, , drop = FALSE], rep(1, p)),
                                error = function(e) NULL)
      if (!is.null(constant_coef) && all(is.finite(constant_coef)) &&
          all(as.vector(X %*% constant_coef) == 1)) {
        coef_offset <- y[[1L]] * constant_coef
        fitted_offset <- rep(y[[1L]], J)
      } else {
        coef_offset <- qr.coef(qr(X), y)
        fitted_offset <- as.vector(X %*% coef_offset)
      }
    }
  }
  y_work <- y - fitted_offset
  if (!all(is.finite(y_work^2))) {
    .rp_abort("The estimates are too far apart for double precision: their squared differences overflow.",
              fix = "Rescale the estimates and the standard errors together.",
              class = "rp_numeric_error")
  }
  if (!is.null(X)) {
    # Then remove the least-squares fit of what is left, which can still be
    # large compared with the standard errors (a steep slope, say), so that
    # the weighted fits at every sigma start from small residuals.
    rest <- qr.coef(qr(X), y_work)
    y_work <- y_work - as.vector(X %*% rest)
    coef_offset <- coef_offset + rest
    fitted_offset <- fitted_offset + as.vector(X %*% rest)
  }
  log_det_xtx <- if (is.null(X)) log(J) else {
    as.numeric(determinant(crossprod(X), logarithm = TRUE)$modulus)
  }
  list(y = y, y_work = y_work, coef_offset = coef_offset,
       fitted_offset = fitted_offset,
       v = v, X = X, prior = prior, J = J, p = p, log_gm = centered$log_gm,
       d = centered$d, scale = exp(centered$log_gm / 2), v_max = max(v),
       log_det_xtx = log_det_xtx)
}

# log k(sigma) for a vector of sigma, in blocks of rows.
.rp_log_kernel <- function(kernel, sigma, max_cells = 5e5) {
  y <- kernel$y_work
  v <- kernel$v
  J <- kernel$J
  out <- numeric(length(sigma))
  # Where sigma^2 overflows the kernel is 0: it falls as sigma^(p - J).
  ok <- is.finite(sigma^2) & sigma >= 0
  out[!ok] <- -Inf
  idx <- which(ok)
  if (length(idx) > 0L) {
    if (is.null(kernel$X)) {
      rows <- max(1L, floor(max_cells / J))
      for (start in seq.int(1L, length(idx), by = rows)) {
        at <- idx[start:min(length(idx), start + rows - 1L)]
        den <- outer(sigma[at]^2, v, "+")
        w <- 1 / den
        W <- rowSums(w)
        mu <- as.vector(w %*% y) / W
        r <- matrix(y, nrow = length(at), ncol = J, byrow = TRUE) - mu
        out[at] <- -0.5 * log(W) - 0.5 * rowSums(log(den)) - 0.5 * rowSums(w * r^2)
      }
    } else {
      out[idx] <- vapply(sigma[idx], function(s) .rp_log_kernel_x(kernel, s), numeric(1))
    }
    out[idx] <- out[idx] + .rp_prior_log_density(kernel$prior, sigma[idx])
  }
  out[is.nan(out)] <- -Inf
  out
}

# One value of the kernel with a mean model, from the QR decomposition of
# W^(1/2) X.
.rp_log_kernel_x <- function(kernel, s) {
  den <- s^2 + kernel$v
  sw <- 1 / sqrt(den)
  qx <- qr(kernel$X * sw)
  r <- qr.resid(qx, kernel$y_work * sw)
  log_det <- 2 * sum(log(abs(diag(qr.R(qx)))))
  -0.5 * log_det - 0.5 * sum(log(den)) - 0.5 * sum(r^2)
}

# The weighted least-squares fit at sigma and the diagonal of
# X (X'WX)^-1 X', for the moments of the site effects with a mean model.
.rp_mean_fit <- function(kernel, s) {
  den <- s^2 + kernel$v
  w <- 1 / den
  if (is.null(kernel$X)) {
    W <- sum(w)
    beta <- sum(w * kernel$y_work) / W
    return(list(fitted = rep(beta, kernel$J) + kernel$fitted_offset,
                fitted_centered = rep(beta, kernel$J), leverage = rep(1 / W, kernel$J),
                coef = beta + kernel$coef_offset, coef_centered = beta))
  }
  sw <- sqrt(w)
  qx <- qr(kernel$X * sw)
  coef <- qr.coef(qx, kernel$y_work * sw)
  Q <- qr.Q(qx)
  # x_j' (X'WX)^-1 x_j = h_jj / w_j, with h_jj the leverage of W^(1/2) X.
  fitted <- as.vector(kernel$X %*% coef)
  list(fitted = fitted + kernel$fitted_offset, fitted_centered = fitted,
       leverage = rowSums(Q^2) / w,
       coef = coef + kernel$coef_offset, coef_centered = coef)
}

# An upper bound on the integral of sigma^q k(sigma) over [B, Inf), on the
# log scale and before normalization. For sigma >= B,
#   k(sigma) <= |X'X|^(-1/2) (1 + v_max / B^2)^(p/2) sigma^(p - J) p(sigma),
# because exp(-Q/2) <= 1, prod (sigma^2 + v_j)^(-1/2) <= sigma^-J and
# X'WX >= X'X / (sigma^2 + v_max). For p = 1 this is the bound of Appendix D.3
# of the paper. The prior is bounded on [B, Inf) by its tail.
.rp_log_tail_bound <- function(kernel, B, q = 0) {
  prior <- kernel$prior
  p <- kernel$p
  J <- kernel$J
  log_c <- -0.5 * kernel$log_det_xtx + 0.5 * p * log1p(kernel$v_max / B^2)
  c_pow <- q + p - J
  switch(prior$type,
    flat = {
      if (c_pow >= -1) return(Inf)
      log_c + (c_pow + 1) * log(B) - log(-(c_pow + 1))
    },
    half_t = ,
    half_cauchy = {
      nu <- prior$df
      e <- c_pow - nu - 1
      if (e >= -1) return(Inf)
      log_c + 0.5 * (nu + 1) * log(nu * prior$scale^2) + (e + 1) * log(B) - log(-(e + 1))
    },
    half_normal = {
      s <- prior$scale
      if (c_pow <= 0) {
        log_c + c_pow * log(B) + log(s) + 0.5 * log(2 * pi) +
          stats::pnorm(B / s, lower.tail = FALSE, log.p = TRUE)
      } else {
        k <- (c_pow + 1) / 2
        log_c + k * log(2 * s^2) - log(2) + lgamma(k) +
          stats::pgamma(B^2 / (2 * s^2), k, lower.tail = FALSE, log.p = TRUE)
      }
    },
    uniform = if (B >= prior$upper) -Inf else Inf
  )
}

# The orders q of the moments of sigma whose tails the rules bound: the mass
# and the moments that exist among 1, 2 and 4.
.rp_tail_orders <- function(kernel) {
  q <- c(0, 1, 2, 4)
  q[vapply(q, .rp_moment_exists, logical(1), J = kernel$J, p = kernel$p, prior = kernel$prior)]
}

# The limits that the normalized tail bounds must meet, on the log scale:
# 1e-12 for the mass, 1e-10 scale^q for the moment of order q. On the
# natural scale scale^q underflows or overflows for extreme units.
.rp_log_tail_limit <- function(kernel, q) {
  if (q == 0) log(1e-12) else log(1e-10) + q * log(kernel$scale)
}

# A starting point for the upper end of sigma: well beyond the spread of the
# estimates and of their standard errors.
.rp_upper_guess <- function(kernel) {
  y <- kernel$y_work
  v <- kernel$v
  spread <- if (length(y) > 1L) c(2 * stats::sd(y), 3 * sqrt(max(0, stats::var(y) - mean(v)))) else 0
  guess <- max(c(spread, 3 * sqrt(stats::median(v)), kernel$scale), na.rm = TRUE) * exp(2)
  if (kernel$prior$type == "uniform") min(guess, kernel$prior$upper) else guess
}
