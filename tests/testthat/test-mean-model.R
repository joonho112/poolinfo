test_that("a mean model needs more sites than columns", {
  x <- 1:3
  expect_rp_error(reliability_posterior(c(0.1, 0.3, 0.2), se = rep(0.1, 3), mean_model = ~ x),
                  regexp = "at least 4 sites")
  expect_rp_error(reliability_posterior(c(0.1, 0.3, 0.2, 0.4), se = rep(0.1, 4),
                                        mean_model = ~ x2, data = data.frame(x2 = c(1, 1, 1, 1))),
                  regexp = "not of full rank")
  expect_rp_error(reliability_posterior(1:4, se = rep(1, 4), mean_model = ~ 0), regexp = "no column|`~ 0`")
  expect_rp_error(reliability_posterior(1:4, se = rep(1, 4), mean_model = "x"), regexp = "one-sided formula")
})

test_that("an intercept-only model is the common mean", {
  a <- reliability_posterior(estimate, se = se, data = eight_schools)
  b <- reliability_posterior(estimate, se = se, data = eight_schools, mean_model = ~ 1)
  expect_identical(a$summary, b$summary)
})

test_that("site moments with a mean model follow from a joint normal with a diffuse beta", {
  withr::local_seed(121)
  J <- 7
  x <- seq(-1, 1, length.out = J)
  X <- cbind(1, x)
  v <- runif(J, 0.01, 0.05)
  y <- 0.2 + 0.4 * x + rnorm(J, 0, sqrt(0.03 + v))
  k <- .rp_kernel(y, v, X)
  s <- 0.15
  cond <- .rp_site_conditional(k, s)
  # tau = X beta + u, y = tau + e, beta ~ N(0, c I) with c large.
  c0 <- 1e8
  cov_tau <- c0 * X %*% t(X) + diag(s^2, J)
  cov_y <- cov_tau + diag(v)
  post_mean <- as.vector(cov_tau %*% solve(cov_y, y))
  post_var <- diag(cov_tau - cov_tau %*% solve(cov_y, cov_tau))
  expect_equal(as.vector(cond$mean), post_mean, tolerance = 1e-6)
  expect_equal(as.vector(cond$var), post_var, tolerance = 1e-6)
})

test_that("the moments of sigma_tau need more sites with more columns", {
  withr::local_seed(122)
  x <- rnorm(5)
  fit <- reliability_posterior(rnorm(5), se = rep(0.5, 5) * runif(5, 0.8, 1.2), mean_model = ~ x)
  # p = 2: E(sigma^2) needs J > 2 + 2 + 1 = 5.
  expect_true(is.na(fit$summary$mean[fit$summary$quantity == "sigma_tau2"]))
  expect_false(is.na(fit$summary$mean[fit$summary$quantity == "sigma_tau"]))
  expect_true(all(fit$checks$status == "ok"))
})

test_that("the variables of a mean model come from the data", {
  withr::local_seed(123)
  d <- data.frame(y = rnorm(10), s = runif(10, 0.2, 0.5), z = rnorm(10), g = rep(c("a", "b"), 5))
  fit <- reliability_posterior(y, se = s, data = d, mean_model = ~ z + g)
  expect_true(all(c("coef_(Intercept)", "coef_z", "coef_gb") %in% fit$summary$quantity))
  a <- reliability_posterior(y, se = s, data = d, mean_model = ~ z,
                             control = posterior_control(engine = "adaptive"))
  n <- reliability_posterior(y, se = s, data = d, mean_model = ~ z)
  expect_lt(max(abs(a$summary$mean[1:5] - n$summary$mean[1:5])), tol_bounded)
  expect_lt(max(abs(a$sites$effect - n$sites$effect)), 1e-8)
})

test_that("the kernel is metafor's restricted likelihood, for p = 1, 2 and 3", {
  skip_if_not_installed("metafor")
  withr::local_seed(124)
  J <- 25
  d <- data.frame(x1 = rnorm(J), x2 = runif(J))
  d$vi <- runif(J, 0.01, 0.08)
  d$yi <- 0.1 + 0.3 * d$x1 - 0.2 * d$x2 + rnorm(J, 0, sqrt(0.04 + d$vi))
  forms <- list(~ 1, ~ x1, ~ x1 + x2)
  for (f in forms) {
    X <- if (identical(f, ~ 1)) NULL else stats::model.matrix(f, d)
    k <- .rp_kernel(d$yi, d$vi, X)
    fit <- metafor::rma(yi, vi, mods = f, data = d, method = "REML",
                        control = list(threshold = 1e-12, maxiter = 1000))
    # The maximum of the kernel over sigma_tau is at the REML estimate.
    opt <- stats::optimize(function(s) -.rp_log_kernel(k, s), c(0, 2), tol = 1e-12)$minimum
    if (fit$tau2 > 0) {
      expect_equal(opt^2, fit$tau2, tolerance = 1e-6)
    } else {
      expect_lt(opt^2, 1e-8)
    }
    # Differences of the kernel are differences of the restricted log likelihood.
    ll <- function(t2) {
      as.numeric(stats::logLik(metafor::rma(yi, vi, mods = f, data = d, method = "REML", tau2 = t2)))
    }
    t2 <- c(0.01, 0.05, 0.2)
    expect_equal(diff(.rp_log_kernel(k, sqrt(t2))), diff(vapply(t2, ll, 0)), tolerance = 1e-8)
  }
})

test_that("a mean model is named in the prior line, and long names keep the columns", {
  d <- data.frame(y = eight_schools$estimate, se = eight_schools$se, group_size = c(1, 0, 1, 0, 1, 0, 1, 1))
  fit <- reliability_posterior(y, se = se, data = d, mean_model = ~group_size)
  out <- format(fit, use_cli = FALSE)
  expect_true(any(grepl("Prior: flat on the coefficients of the mean model and on sigma_tau", out, fixed = TRUE)))
  # The labels are padded to 26 characters, one more than the longest
  # standard label, "concentration count J_eff".
  rows <- out[grepl("^  (informativeness I|coefficient|sigma_tau) ", out)]
  expect_length(rows, 4L)
  expect_true(all(substr(rows, 29, 29) == " " & substr(rows, 30, 30) != " "))
  long <- format(summary(fit), use_cli = FALSE)
  rows <- c(long[grep("^  informativeness I", long)[[1L]]], long[grepl("^  coefficient ", long)])
  expect_true(all(substr(rows, 29, 29) == " " & substr(rows, 30, 30) != " "))
  # A longer name of a coefficient widens the column.
  names(d)[3] <- "number_of_tutoring_sessions"
  fit <- reliability_posterior(y, se = se, data = d, mean_model = ~number_of_tutoring_sessions)
  out <- format(fit, use_cli = FALSE)
  rows <- out[grepl("^  (informativeness I|coefficient|sigma_tau) ", out)]
  w <- nchar("coefficient number_of_tutoring_sessions")
  expect_true(all(substr(rows, w + 3L, w + 3L) == " " & substr(rows, w + 4L, w + 4L) != " "))
})

test_that("the sites needed for E(sigma_tau^2) count the columns of the mean model", {
  d <- data.frame(y = c(0.1, 0.5, -0.4, 0.2, 0.3), se = c(0.2, 0.3, 0.45, 0.1, 0.2), x = c(0, 1, 0, 1, 0))
  fit <- reliability_posterior(y, se = se, data = d, mean_model = ~x)
  # J > q + p + 1 with q = 2 and p = 2: at least 6 sites.
  expect_message(plugin_gap(fit), "It needs at least 6 sites")
  expect_rp_error(reliability_profile(fit), regexp = "needs at least 6 sites; there are 5")
})
