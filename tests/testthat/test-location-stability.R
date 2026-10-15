test_that("a large common location preserves small represented contrasts", {
  y <- 1e8 + c(-8, -5, -1, 1, 4, 9) * 1e-4
  se <- c(1, 1.2, 1.6, 1.4, 1, 1.7) * 1e-4
  yc <- y - y[1]
  for (engine in c("nodes", "adaptive")) {
    control <- posterior_control(engine = engine)
    for (moderators in c(FALSE, TRUE)) {
      d <- data.frame(y = y, yc = yc, se = se, x = c(-2, -1, 0, 0, 1, 2))
      form <- if (moderators) ~ x else ~ 1
      shifted <- reliability_posterior(y, se = se, data = d, mean_model = form, control = control)
      centered <- reliability_posterior(yc, se = se, data = d, mean_model = form, control = control)
      expect_equal(shifted$summary[1:8, ], centered$summary[1:8, ], tolerance = 1e-11)
      expect_equal(shifted$sites$effect_sd / centered$sites$effect_sd, rep(1, length(y)), tolerance = 1e-11)
      expect_lte(max(abs((shifted$sites$effect - y[1]) - centered$sites$effect)),
                 2 * .Machine$double.eps * abs(y[1]))
      i <- grepl("^(effect_mean|coef_)", shifted$summary$quantity)
      expect_equal(shifted$summary$sd[i] / centered$summary$sd[i], rep(1, sum(i)), tolerance = 1e-11)
      expect_true(all(shifted$checks$status == "ok"))
      if (!moderators) expect_equal(shifted$summary$mean[1], 0.96405246977799, tolerance = 1e-11)
    }
  }
})

test_that("centering a no-intercept fit stays within its model span", {
  x <- c(0.5, 1, 1.5, 2, 2.5, 3)
  X <- cbind(x = x)
  y <- 0.7 + sin(x)
  v <- c(0.2, 0.3, 0.4, 0.2, 0.5, 0.6)
  kernel <- .rp_kernel(y, v, X)
  for (sigma in c(0, 0.2, 2)) {
    w <- 1 / (sigma^2 + v)
    beta <- sum(w * x * y) / sum(w * x^2)
    residual <- y - x * beta
    reference <- -0.5 * log(sum(w * x^2)) - 0.5 * sum(log(sigma^2 + v)) - 0.5 * sum(w * residual^2)
    expect_equal(.rp_log_kernel(kernel, sigma), reference, tolerance = 1e-13)
    expect_equal(unname(.rp_mean_fit(kernel, sigma)$coef), beta, tolerance = 1e-13)
  }
  # A constant shift is outside this one-column mean model.
  expect_gt(abs(diff(.rp_log_kernel(.rp_kernel(y + 10, v, X), c(0.2, 2))) -
                  diff(.rp_log_kernel(kernel, c(0.2, 2)))), 1)
})

test_that("an exact constant in a no-intercept basis preserves large-location contrasts", {
  y <- 1e8 + c(-8, -5, -1, 1, 4, 9) * 1e-4
  x <- seq(-1, 1, length.out = length(y))
  d <- data.frame(y = y, yc = y - y[[1L]], u = 1 + x, v = 1 - x,
                  se = c(1, 1.2, 1.6, 1.4, 1, 1.7) * 1e-4)
  X <- stats::model.matrix(~ 0 + u + v, d)
  k <- .rp_kernel(y, d$se^2, X)
  # The offset lies in the span of X and leaves residuals on the scale of
  # the standard errors, while offset plus residuals give back the estimates.
  expect_equal(as.vector(X %*% k$coef_offset), k$fitted_offset, tolerance = 1e-15)
  expect_lt(max(abs(k$y_work)), 1e-2)
  expect_lte(max(abs(k$y_work + k$fitted_offset - y)), 2 * .Machine$double.eps * y[[1L]])
  for (engine in c("nodes", "adaptive")) {
    control <- posterior_control(engine = engine)
    shifted <- reliability_posterior(y, se = se, mean_model = ~ 0 + u + v,
                                      data = d, control = control)
    centered <- reliability_posterior(yc, se = se, mean_model = ~ 0 + u + v,
                                       data = d, control = control)
    expect_equal(shifted$summary[1:8, ], centered$summary[1:8, ], tolerance = 1e-11)
    expect_equal(shifted$sites$effect_sd / centered$sites$effect_sd,
                 rep(1, length(y)), tolerance = 1e-11)
    rows <- grepl("^coef_", shifted$summary$quantity)
    expect_equal(shifted$summary$sd[rows] / centered$summary$sd[rows], rep(1, 2), tolerance = 1e-11)
    expect_lte(max(abs(shifted$summary$mean[rows] - y[[1L]] / 2 -
                          centered$summary$mean[rows])), 2 * .Machine$double.eps * y[[1L]])
    expect_lte(max(abs(shifted$sites$effect - y[[1L]] - centered$sites$effect)),
                 2 * .Machine$double.eps * y[[1L]])
  }
  # A basis that only nearly spans a constant must use its original model.
  near_X <- cbind(u = 1 + x, v = 1 - x + 1e-7 * x^2)
  near_k <- .rp_kernel(y, d$se^2, near_X)
  expect_equal(near_k$coef_offset, qr.coef(qr(near_X), y), tolerance = 1e-13)
  expect_equal(near_k$fitted_offset, as.vector(near_X %*% qr.coef(qr(near_X), y)), tolerance = 1e-13)
  expect_false(all(near_k$fitted_offset == y[[1L]]))
})

test_that("a moment that exists can still have too heavy a tail to bound", {
  for (df in c(0.05, 0.01)) {
    prior <- sd_prior("half_t", scale = 1, df = df)
    expect_true(.rp_moment_exists(1, J = 2, p = 1, prior = prior))
    for (engine in c("nodes", "adaptive")) {
      err <- suppressMessages(expect_rp_error(
        reliability_posterior(c(-1, 1), se = c(1, 1.2), prior = prior,
                              control = posterior_control(engine = engine)),
        "rp_numeric_error", regexp = "moment of order 1"))
      expect_match(conditionMessage(err), "B^(-", fixed = TRUE)
      expect_match(conditionMessage(err), "The moment exists", fixed = TRUE)
      expect_match(conditionMessage(err), "lighter tail", fixed = TRUE)
    }
  }
})

test_that("a steep fitted mean in the estimates leaves the fit unchanged", {
  withr::local_seed(1)
  x <- stats::rnorm(12)
  se <- stats::runif(12, 0.2, 0.7)
  y <- 0.2 + 0.4 * x + stats::rnorm(12, 0, sqrt(0.15 + se^2))
  d <- data.frame(x = x, y = y, steep = y + 1e6 * x)
  base <- reliability_posterior(y, se = se, data = d, mean_model = ~ x)
  steep <- reliability_posterior(steep, se = se, data = d, mean_model = ~ x)
  rows <- c("I", "mean", "H", "J_eff", "sigma_tau")
  expect_equal(steep$summary$mean[steep$summary$quantity %in% rows],
               base$summary$mean[base$summary$quantity %in% rows], tolerance = 1e-9)
  expect_equal(steep$sites$effect - 1e6 * x, base$sites$effect, tolerance = 1e-6)
})

test_that("coefficient means move with the estimates and agree with a direct integral", {
  x <- c(-1.6, -1, -0.7, -0.2, 0, 0.3, 0.6, 0.9, 1.3, 1.8)
  v <- c(0.2, 0.3, 0.25, 0.4, 0.3, 0.2, 0.35, 0.3, 0.25, 0.4)
  y <- 0.4 + 0.5 * x + c(0.3, -0.5, 0.2, 0.6, -0.4, 0.1, -0.2, 0.5, -0.6, 0.2)
  d <- data.frame(y = y, x = x, v = v, y2 = y + 3 - 2 * x)
  coef_mean <- function(f) f$summary$mean[grepl("^coef_", f$summary$quantity)]
  fit <- reliability_posterior(y, se2 = v, data = d, mean_model = ~ x)
  moved <- reliability_posterior(y2, se2 = v, data = d, mean_model = ~ x)
  expect_equal(coef_mean(moved), coef_mean(fit) + c(3, -2), tolerance = 1e-9)
  # The posterior mean of the coefficients is the average over sigma of the
  # weighted least-squares fit, under the restricted likelihood.
  X <- cbind(1, x)
  wls <- function(s) {
    w <- 1 / (s^2 + v)
    A <- crossprod(X * sqrt(w))
    b <- solve(A, crossprod(X, w * y))
    list(b = b, log_k = -0.5 * as.numeric(determinant(A)$modulus) - 0.5 * sum(log(s^2 + v)) -
           0.5 * sum(w * (y - X %*% b)^2))
  }
  top <- max(vapply(seq(0.01, 5, length.out = 200), function(s) wls(s)$log_k, 0))
  dens <- function(s, k = 0) vapply(s, function(si) {
    f <- wls(si)
    exp(f$log_k - top) * if (k == 0) 1 else f$b[[k]]
  }, 0)
  norm <- stats::integrate(dens, 0, Inf, rel.tol = 1e-10)$value
  direct <- vapply(1:2, function(k) stats::integrate(dens, 0, Inf, k = k, rel.tol = 1e-10)$value / norm, 0)
  expect_equal(coef_mean(fit), direct, tolerance = 1e-7)
})

test_that("estimates too far apart for double precision give a clear error", {
  expect_rp_error(reliability_posterior(c(1e300, -1, 0, 1, 2), se = rep(1, 5)),
                  "rp_numeric_error", regexp = "too far apart")
})
