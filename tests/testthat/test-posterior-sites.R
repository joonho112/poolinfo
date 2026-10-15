# Site summaries against a brute-force integral over sigma_tau and the mean,
# on a fine grid, for a small data set.
brute_force_sites <- function(y, v) {
  scale <- sqrt(exp(mean(log(v))))
  # A grid in sigma near 0, where the density is positive under the flat
  # prior, and a grid in log(sigma) above it.
  # Trapezoid weights on both grids.
  trap <- function(n, h) h * c(0.5, rep(1, n - 2), 0.5)
  s0 <- seq(0, scale * 1e-4, length.out = 401)
  t <- seq(log(scale * 1e-4), log(scale) + 12, length.out = 40001)
  s <- c(s0, exp(t))
  k <- .rp_kernel(y, v)
  lk <- .rp_log_kernel(k, s)
  dx <- c(trap(401, diff(s0)[[1]]), exp(t) * trap(40001, diff(t)[[1]]))
  w <- exp(lk - max(lk)) * dx
  w <- w / sum(w)
  cc <- .rp_site_conditional(k, s)
  pm <- colSums(w * cc$mean)
  pv <- colSums(w * (cc$var + cc$mean^2)) - pm^2
  list(pm = pm, psd = sqrt(pv), S = colSums(w * cc$S), w = w, cc = cc)
}

test_that("site posterior means and SDs match a brute-force integral", {
  y <- c(0.30, 0.10, 0.25, -0.05, 0.40, 0.00, 0.22, 0.15)
  v <- c(0.01, 0.015, 0.02, 0.03, 0.05, 0.06, 0.012, 0.025)
  fit <- reliability_posterior(y, se2 = v)
  bf <- brute_force_sites(y, v)
  # The grid integral is accurate to about 1e-9 here.
  expect_lt(max(abs(fit$sites$effect - bf$pm)), 1e-8)
  expect_lt(max(abs(fit$sites$effect_sd - bf$psd)), 1e-8)
  expect_lt(max(abs(fit$sites$reliability - bf$S)), 1e-8)
})

test_that("site intervals are the quantiles of the mixture over sigma_tau", {
  y <- c(0.30, 0.10, 0.25, -0.05, 0.40, 0.00)
  v <- c(0.01, 0.015, 0.02, 0.03, 0.05, 0.06)
  fit <- reliability_posterior(y, se2 = v)
  tab <- site_table(fit)
  bf <- brute_force_sites(y, v)
  for (j in seq_along(y)) {
    F <- function(x) sum(bf$w * stats::pnorm((x - bf$cc$mean[, j]) / sqrt(bf$cc$var[, j])))
    lo <- stats::uniroot(function(x) F(x) - 0.025, c(-5, 5), tol = 1e-12)$root
    hi <- stats::uniroot(function(x) F(x) - 0.975, c(-5, 5), tol = 1e-12)$root
    expect_lt(abs(tab$effect_lower[[j]] - lo), 1e-7)
    expect_lt(abs(tab$effect_upper[[j]] - hi), 1e-7)
  }
  expect_true(all(tab$reliability_lower < tab$reliability & tab$reliability < tab$reliability_upper))
  expect_equal(sum(tab$share), 1, tolerance = 1e-12)
  expect_equal(attr(tab, "level"), 0.95)
  tab80 <- site_table(fit, level = 0.8)
  expect_true(all(tab80$effect_lower > tab$effect_lower))
  expect_identical(as.data.frame(fit), tab)
})

test_that("the reliability intervals are the reliabilities at the quantiles of sigma_tau", {
  fit <- reliability_posterior(estimate, se = se, data = eight_schools)
  q <- fit$sigma_quantiles
  v <- eight_schools$se^2
  expect_equal(fit$sites$reliability_lower, q[[1]]^2 / (q[[1]]^2 + v))
  expect_equal(fit$sites$reliability_upper, q[[3]]^2 / (q[[3]]^2 + v))
})

test_that("a mean model fits the effects around a regression line", {
  withr::local_seed(91)
  x <- seq(-1, 1, length.out = 12)
  v <- runif(12, 0.01, 0.05)
  y <- 0.2 + 0.3 * x + rnorm(12, 0, sqrt(0.02 + v))
  fit <- reliability_posterior(y, se2 = v, mean_model = ~ x)
  expect_true(all(c("coef_(Intercept)", "coef_x") %in% fit$summary$quantity))
  expect_true(all(fit$checks$status == "ok"))
  # The kernel with covariates equals a brute-force integral over beta up to
  # a constant (Appendix D with p columns).
  k <- .rp_kernel(y, v, cbind(1, x))
  brute <- function(s) {
    f <- function(b0, b1) {
      mu <- b0 + b1 * x
      exp(sum(stats::dnorm(y, mu, sqrt(s^2 + v), log = TRUE)) + 0)
    }
    g <- expand.grid(b0 = seq(-2, 2, length.out = 401), b1 = seq(-3, 3, length.out = 401))
    vals <- mapply(f, g$b0, g$b1)
    log(sum(vals) * (4 / 400) * (6 / 400))
  }
  s <- c(0.05, 0.1, 0.2)
  diffs <- vapply(s, brute, 0) - .rp_log_kernel(k, s)
  expect_lt(diff(range(diffs)), 1e-4)
})
