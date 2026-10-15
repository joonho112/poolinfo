# Draws from the posterior, for plots and for summaries that users compute
# themselves; the summaries that the package reports come from the
# integrals. sigma_tau is drawn by inverting a table of its CDF, built on a
# grid by the trapezoid rule and checked against the exact CDF; the mean (or
# the coefficients) and the site effects are then drawn from their normal
# conditional distributions (Appendix D.3 of the paper).

# A table of the CDF of sigma on a grid of t = log(sigma) between its
# quantiles 1e-12 and 1 - 1e-12, built by the trapezoid rule and pinned to the
# exact CDF at both ends. The grid is doubled until the inverse is within
# `max_error` of the exact CDF at a set of probabilities.
.rp_inverse_cdf <- function(law, kernel, max_error = 2e-8) {
  ends <- law$quantile(c(1e-12, 1 - 1e-12))
  t_lo <- log(ends[[1L]])
  t_hi <- log(ends[[2L]])
  F_ends <- law$cdf(ends)
  checked <- sort(unique(c(1e-6, 1e-4, 1e-3, seq(0.005, 0.995, length.out = 41), 0.999, 0.9999,
                           1 - 1e-6)))
  n <- 4097L
  repeat {
    t <- seq(t_lo, t_hi, length.out = n)
    d <- law$density(exp(t)) * exp(t)
    cum <- c(0, cumsum(diff(t) * (d[-n] + d[-1L]) / 2))
    cdf <- F_ends[[1L]] + cum / cum[[n]] * (F_ends[[2L]] - F_ends[[1L]])
    keep <- c(TRUE, diff(cdf) > 0)
    cc <- cdf[keep]
    tt <- t[keep]
    inverse <- function(p) exp(stats::approx(cc, tt, xout = p, rule = 2, ties = "ordered")$y)
    error <- max(abs(law$cdf(inverse(checked)) - checked))
    if (error <= max_error || n >= 262145L) break
    n <- 2L * n - 1L
  }
  list(inverse = inverse, error = error, n = n)
}

# Draws of sigma, of the mean model's coefficients and of the diagnostics,
# and optionally of the site effects.
.rp_draws <- function(law, kernel, ndraws, sites = FALSE) {
  table <- .rp_inverse_cdf(law, kernel)
  sigma <- table$inverse(stats::runif(ndraws))
  p <- kernel$p
  coef <- matrix(NA_real_, ndraws, p)
  fun <- .rp_functionals_at(kernel, sigma)
  out <- data.frame(draw = seq_len(ndraws), sigma_tau = sigma)
  fits <- lapply(sigma, function(s) .rp_mean_fit(kernel, s))
  for (i in seq_len(ndraws)) {
    if (is.null(kernel$X)) {
      coef[i, ] <- stats::rnorm(1L, fits[[i]]$coef_centered, sqrt(fits[[i]]$leverage[[1L]]))
    } else {
      w <- 1 / (sigma[[i]]^2 + kernel$v)
      cov <- solve(crossprod(kernel$X * sqrt(w)))
      coef[i, ] <- fits[[i]]$coef_centered + as.vector(t(chol(cov)) %*% stats::rnorm(p))
    }
  }
  reported_coef <- sweep(coef, 2L, kernel$coef_offset, "+")
  if (is.null(kernel$X)) {
    out$effect_mean <- reported_coef[, 1L]
  } else {
    for (k in seq_len(p)) out[[paste0("coef_", colnames(kernel$X)[[k]])]] <- reported_coef[, k]
  }
  out$I <- fun[, "I"]
  out$mean <- fun[, "mean"]
  out$H <- fun[, "H"]
  out$J_eff <- fun[, "J_eff"]
  if (sites) {
    mu <- if (is.null(kernel$X)) {
      matrix(coef[, 1L], ndraws, kernel$J)
    } else {
      coef %*% t(kernel$X)
    }
    eta <- outer(2 * log(sigma), -log(kernel$v), "+")
    S <- stats::plogis(eta)
    one_minus <- stats::plogis(-eta)
    y <- matrix(kernel$y_work, ndraws, kernel$J, byrow = TRUE)
    center <- S * y + one_minus * mu
    sd <- sqrt(matrix(kernel$v, ndraws, kernel$J, byrow = TRUE) * S)
    effects <- center + sd * matrix(stats::rnorm(ndraws * kernel$J), ndraws, kernel$J)
    effects <- sweep(effects, 2L, kernel$fitted_offset, "+")
    colnames(effects) <- paste0("effect_", seq_len(kernel$J))
    out <- cbind(out, as.data.frame(effects))
  }
  attr(out, "inverse_cdf_error") <- table$error
  out
}

# Runs `expr` with the seed set, and puts back the state of the generator of
# the calling session afterwards.
.rp_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  set.seed(seed)
  force(expr)
}

# A seed from the clock, for when none is given.
.rp_clock_seed <- function() {
  as.integer((as.numeric(Sys.time()) * 1000) %% 2147483647)
}
