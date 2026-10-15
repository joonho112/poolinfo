test_that("sd_prior() builds the five families and refuses what they lack", {
  expect_identical(sd_prior()$type, "flat")
  expect_identical(sd_prior("half_normal", scale = 2)$tail, Inf)
  expect_identical(sd_prior("half_t", scale = 1, df = 3)$tail, 4)
  expect_identical(sd_prior("half_cauchy", scale = 1)$tail, 2)
  expect_identical(sd_prior("uniform", upper = 5)$upper, 5)
  expect_rp_error(sd_prior("half_normal"), regexp = "needs `scale`")
  expect_rp_error(sd_prior("half_t", scale = 1), regexp = "needs `df`")
  expect_rp_error(sd_prior("flat", scale = 1), regexp = "has no `scale`")
  expect_rp_error(sd_prior("uniform", upper = -1), regexp = "positive")
  expect_rp_error(sd_prior("gamma"), regexp = "one of")
  expect_output(print(sd_prior("half_t", scale = 0.5, df = 3)),
                "half-t on sigma_tau with scale 0.5 and 3 degrees of freedom")
})

# The posterior under a prior by brute force on a fine grid.
brute_posterior <- function(y, v, prior, fun) {
  scale <- sqrt(exp(mean(log(v))))
  trap <- function(n, h) h * c(0.5, rep(1, n - 2), 0.5)
  s0 <- seq(0, scale * 1e-4, length.out = 401)
  # The grid ends at the upper limit of a uniform prior, where its density
  # jumps to 0.
  top <- if (prior$type == "uniform") log(prior$upper) else log(scale) + 14
  t <- seq(log(scale * 1e-4), top, length.out = 60001)
  s <- c(s0, exp(t))
  k <- .rp_kernel(y, v, prior = prior)
  lk <- .rp_log_kernel(k, s)
  w <- exp(lk - max(lk)) * c(trap(401, diff(s0)[[1]]), exp(t) * trap(60001, diff(t)[[1]]))
  sum(w * fun(s, k)) / sum(w)
}

test_that("posteriors under proper priors match a brute-force integral", {
  withr::local_seed(111)
  v <- runif(6, 0.02, 0.2)
  y <- rnorm(6, 0.1, sqrt(0.05 + v))
  I_at <- function(s, k) .rp_functionals_at(k, s)[, "I"]
  for (prior in list(sd_prior("half_normal", scale = 0.3), sd_prior("half_t", scale = 0.3, df = 3),
                     sd_prior("half_cauchy", scale = 0.3), sd_prior("uniform", upper = 1))) {
    fit <- reliability_posterior(y, se2 = v, prior = prior)
    expect_lt(abs(fit$summary$mean[1] - brute_posterior(y, v, prior, I_at)), 1e-8)
    expect_true(all(fit$checks$status == "ok"))
    a <- reliability_posterior(y, se2 = v, prior = prior, control = posterior_control(engine = "adaptive"))
    expect_lt(max(abs(a$summary$mean[1:5] - fit$summary$mean[1:5])), tol_bounded)
  }
})

test_that("moments exist according to the tail of the prior", {
  y <- c(0.1, 0.5, -0.2)
  se <- c(0.2, 0.3, 0.25)
  flat <- reliability_posterior(y, se = se)
  expect_true(is.na(flat$summary$mean[flat$summary$quantity == "sigma_tau2"]))
  cauchy <- reliability_posterior(y, se = se, prior = sd_prior("half_cauchy", scale = 0.3))
  # J > q + p + 1 - k: 3 > 2 + 1 + 1 - 2.
  expect_false(is.na(cauchy$summary$mean[cauchy$summary$quantity == "sigma_tau2"]))
  # The fourth moment needs J > 4 + 1 + 1 - 2 = 4.
  expect_true(is.na(cauchy$summary$sd[cauchy$summary$quantity == "sigma_tau2"]))
  normal <- reliability_posterior(y, se = se, prior = sd_prior("half_normal", scale = 0.3))
  expect_false(anyNA(normal$summary$sd))
})

test_that("proper priors allow 2 sites, with a message", {
  expect_message(fit <- reliability_posterior(c(0.1, 0.4), se = c(0.2, 0.2),
                                              prior = sd_prior("half_normal", scale = 0.3)),
                 class = "rp_message_few_sites")
  expect_true(all(fit$checks$status == "ok"))
  expect_rp_error(reliability_posterior(0.1, se = 0.2, prior = sd_prior("half_normal", scale = 0.3)),
                  regexp = "At least 2 sites")
})

test_that("the bounds on the tail are above the tail itself", {
  withr::local_seed(112)
  v <- runif(5, 0.02, 0.2)
  y <- rnorm(5, 0, sqrt(0.05 + v))
  for (prior in list(sd_prior(), sd_prior("half_t", scale = 0.3, df = 2),
                     sd_prior("half_normal", scale = 0.3))) {
    k <- .rp_kernel(y, v, prior = prior)
    for (B in c(1, 3, 10)) {
      for (q in .rp_tail_orders(k)) {
        f <- function(s) s^q * exp(.rp_log_kernel(k, s))
        tail <- stats::integrate(f, B, Inf, rel.tol = 1e-10)$value
        expect_gte(exp(.rp_log_tail_bound(k, B, q)), tail * (1 - 1e-6))
      }
    }
  }
})

test_that("the first five OCRS sites give the stored posterior means of I under three priors", {
  g <- ocrs_analyses()$pooled
  y <- g$y[1:5]
  v <- g$v[1:5]
  I_mean <- function(prior) reliability_posterior(y, se2 = v, prior = prior)$summary$mean[1]
  expect_equal(round(I_mean(sd_prior()), 4), 0.6449)
  expect_equal(round(I_mean(sd_prior("half_normal", scale = 0.5)), 4), 0.6123)
  expect_equal(round(I_mean(sd_prior("half_t", scale = 0.5, df = 3)), 4), 0.6083)
})
