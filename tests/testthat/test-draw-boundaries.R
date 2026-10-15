test_that("zero SD draws use inverse-variance site shares and degenerate effects", {
  x <- reliability_draws(c(0, 0), se = c(0.5, 1), estimate = c(1, 2), mean = c(4, 4))
  sites <- site_table(x)
  expect_equal(sites$share, c(0.8, 0.2))
  expect_equal(sites$reliability, c(0, 0))
  expect_equal(sites$effect, c(4, 4))
  expect_equal(sites$effect_sd, c(0, 0))
  expect_equal(sites$effect_lower, c(4, 4))
  expect_equal(sites$effect_upper, c(4, 4))
  expect_true(all(is.finite(as.matrix(sites[-1]))))
})

test_that("zero SD components with different means retain their discrete quantiles", {
  x <- reliability_draws(c(0, 0), se = c(0.5, 1), estimate = c(1, 2),
                          mean = c(1, 3), level = 0.5)
  sites <- site_table(x)
  expect_equal(sites$effect, c(2, 2))
  expect_equal(sites$effect_sd, c(1, 1))
  expect_equal(sites$effect_lower, c(1, 1))
  expect_equal(sites$effect_upper, c(3, 3))
  center <- matrix(c(1, 1, 3, 7), ncol = 1)
  expect_equal(.rp_mixture_quantile(c(0.1, 0.2, 0.3, 0.4), center,
                                    matrix(0, 4, 1), 0.3), 1)
  expect_equal(.rp_mixture_quantile(c(0.1, 0.2, 0.3, 0.4), center,
                                    matrix(0, 4, 1), 0.31), 3)
})

test_that("mixed zero and positive SD draws use atom-aware mixture intervals", {
  sigma <- c(0, 1)
  se <- c(0.5, 1)
  y <- c(1, 2)
  mu <- c(4, 5)
  x <- reliability_draws(sigma, se = se, estimate = y, mean = mu)
  sites <- site_table(x)
  positive_share <- (1 / (1 + se^2)) / sum(1 / (1 + se^2))
  expect_equal(sites$share, (c(0.8, 0.2) + positive_share) / 2)
  expect_equal(sum(sites$share), 1)
  for (j in seq_along(se)) {
    reliability <- 1 / (1 + se[[j]]^2)
    center <- 5 + reliability * (y[[j]] - 5)
    sd <- sqrt(se[[j]]^2 * reliability)
    cdf <- function(z, strict = FALSE) {
      0.5 * stats::pnorm((z - center) / sd) + 0.5 * if (strict) z > 4 else z >= 4
    }
    oracle <- function(p) {
      if (cdf(4, TRUE) <= p && cdf(4) >= p) return(4)
      stats::uniroot(function(z) cdf(z) - p, c(-10, 15), tol = 1e-11)$root
    }
    expect_equal(sites$effect_lower[[j]], oracle(0.025), tolerance = 1e-9)
    expect_equal(sites$effect_upper[[j]], oracle(0.975), tolerance = 1e-9)
    for (p in c(0.025, 0.975)) {
      q <- if (p < 0.5) sites$effect_lower[[j]] else sites$effect_upper[[j]]
      expect_lte(cdf(q, TRUE), p + 1e-10)
      expect_gte(cdf(q), p - 1e-10)
    }
  }
})

test_that("site shares at zero agree with their positive-SD limit across chunks", {
  kernel <- .rp_draws_kernel(c(0.25, 1, 4, 9))
  zero <- .rp_site_shares(kernel, c(0, 0), c(0.25, 0.75), chunk_sites = 1)
  expected <- (1 / kernel$v) / sum(1 / kernel$v)
  expect_equal(zero, expected, tolerance = 1e-13)
  expect_equal(.rp_site_shares(kernel, c(0, 1e-100), c(0.5, 0.5), chunk_sites = 2), expected,
               tolerance = 1e-13)
})

test_that("a distant point mass does not coarsen the quantiles of the rest of a mixture", {
  cases <- list(c(scale = 1, offset = 0), c(scale = 0.25, offset = 1e6),
                c(scale = 1e-6, offset = 0), c(scale = 1e-4, offset = 1e8))
  for (case in cases) {
    scale <- case[["scale"]]
    offset <- case[["offset"]]
    for (direction in c(-1, 1)) {
      x <- reliability_draws(c(0, scale), se2 = scale^2, estimate = offset,
                              mean = c(offset + direction * 1e12 * scale, offset))
      sites <- site_table(x)
      p <- if (direction > 0) 0.025 else 0.975
      q <- if (direction > 0) sites$effect_lower else sites$effect_upper
      reference <- stats::qnorm(if (direction > 0) 0.05 else 0.95) / sqrt(2)
      standardized <- (q - offset) / scale
      resolution <- 8 * .Machine$double.eps * max(1, abs(offset)) / scale
      expect_lte(abs(standardized - reference), 1e-9 + resolution)
      F <- 0.5 * stats::pnorm(standardized * sqrt(2)) +
        0.5 * (direction < 0)
      expect_lte(abs(F - p), 1e-10 + resolution / sqrt(pi))
      atom_end <- if (direction > 0) sites$effect_upper else sites$effect_lower
      expect_identical(atom_end, offset + direction * 1e12 * scale)
    }
  }
})

test_that("high-reliability site draws retain their conditional sampling variance", {
  y <- c(-2, -1, 0, 1, 2) * 1e9
  fit <- reliability_posterior(y, se = rep(1, 5))
  draws <- profile_draws(fit, ndraws = 4000, seed = 1, sites = TRUE)
  effects <- draws[, paste0("effect_", seq_along(y))]
  empirical <- vapply(effects, stats::sd, 0)
  expect_equal(site_table(fit)$effect_sd, rep(1, 5), tolerance = 1e-7)
  # For a normal SD estimate, its Monte Carlo SE is about 1/sqrt(2*n).
  # A fixed bound of 0.08 exceeds six such standard errors at n = 4000.
  expect_true(all(abs(empirical - 1) < 0.08))
  expect_true(all(empirical > 0.9))
})

test_that("ordinary site-effect draws reproduce analytic posterior variances", {
  fit <- reliability_posterior(estimate, se = se, data = eight_schools)
  draws <- profile_draws(fit, ndraws = 20000, seed = 13, sites = TRUE)
  empirical <- vapply(draws[, paste0("effect_", 1:8)], stats::sd, 0)
  expect_true(all(abs(empirical / site_table(fit)$effect_sd - 1) < 0.04))
})

test_that("site draws keep their variance when the reliabilities round to 1", {
  fit <- reliability_posterior(c(-2, -1, 0, 1, 2) * 1e9, se = rep(1, 5))
  draws <- profile_draws(fit, ndraws = 2000, seed = 1, sites = TRUE)
  effects <- as.matrix(draws[, grep("^effect_[0-9]", names(draws))])
  # With S rounded to 1 each site's posterior SD is its standard error.
  expect_equal(site_table(fit)$effect_sd, rep(1, 5), tolerance = 1e-8)
  expect_true(all(abs(apply(effects, 2L, stats::sd) - 1) < 0.1))
})

test_that("intervals from external draws keep the variance of very reliable sites", {
  sigma <- rep(c(1, 1e9), each = 200)
  se <- c(0.5, 1)
  y <- c(1, 2)
  sites <- site_table(reliability_draws(sigma, se = se, estimate = y, mean = rep(0, 400)))
  for (j in 1:2) {
    S <- sigma^2 / (sigma^2 + se[[j]]^2)
    center <- S * y[[j]]
    sd <- sqrt(se[[j]]^2 * S)
    cdf <- function(z) mean(stats::pnorm((z - center) / sd))
    end <- function(p) stats::uniroot(function(z) cdf(z) - p, c(-20, 20), tol = 1e-12)$root
    expect_equal(sites$effect_lower[[j]], end(0.025), tolerance = 1e-8)
    expect_equal(sites$effect_upper[[j]], end(0.975), tolerance = 1e-8)
  }
})

test_that("discrete quantiles at p = k/n pick the k-th smallest value", {
  withr::local_seed(100)
  vals <- stats::rnorm(100)
  for (k in c(1, 10, 37, 99)) {
    q <- .rp_mixture_quantile(rep(1 / 100, 100), matrix(vals, ncol = 1), matrix(0, 100, 1), k / 100)
    expect_identical(q, sort(vals)[[k]])
  }
})
