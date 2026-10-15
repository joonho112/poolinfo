test_that("the ends of a tau^2 interval map to the ends of the monotone diagnostics", {
  se <- eight_schools$se
  m <- map_tau2_interval(se = se, interval = c(estimate = 30, lower = 4, upper = 200))
  expect_identical(m$functional, c("I", "mean", "J_eff", "H"))
  at <- function(t2) reliability_profile(se = se, sigma_tau = sqrt(t2))$stats
  for (row in 1:3) {
    q <- c("I", "mean", "J_eff")[[row]]
    expect_equal(m$lower[[row]], at(4)[[q]], tolerance = 1e-12)
    expect_equal(m$upper[[row]], at(200)[[q]], tolerance = 1e-12)
    expect_equal(m$estimate[[row]], at(30)[[q]], tolerance = 1e-12)
  }
})

test_that("the interval of H is the image of the interval, with its peak inside", {
  se <- eight_schools$se
  v <- se^2
  m <- map_tau2_interval(se = se, interval = c(1, 10000))
  centered <- .rp_center_v(v)
  turning <- .rp_turning_points(centered$d)
  expect_equal(m$upper[[4]], .rp_H_at(turning$roots, centered$d), tolerance = 1e-14)
  narrow <- map_tau2_interval(se = se, interval = c(1, 2))
  expect_equal(narrow$upper[[4]], max(reliability_profile(se = se, sigma_tau = 1)$stats[["H"]],
                                      reliability_profile(se = se, sigma_tau = sqrt(2))$stats[["H"]]))
  expect_true(all(is.na(narrow$estimate)))
})

test_that("a metafor fit gives its Q-profile interval mapped", {
  skip_if_not_installed("metafor")
  withr::local_seed(141)
  vi <- runif(20, 0.01, 0.1)
  yi <- rnorm(20, 0.3, sqrt(0.04 + vi))
  fit <- metafor::rma(yi = yi, vi = vi)
  m <- map_tau2_interval(fit)
  ci <- stats::confint(fit, type = "qp", control = list(tol = 1e-12))$random
  expect_equal(m$lower[[1]], reliability_profile(se2 = vi, sigma_tau = sqrt(ci["tau^2", "ci.lb"]))$stats[["I"]])
  expect_equal(m$estimate[[1]], reliability_profile(se2 = vi, sigma_tau = sqrt(fit$tau2))$stats[["I"]])
  expect_match(attr(m, "method"), "Q-profile")
  expect_identical(names(attr(m, "interval")), c("estimate", "lower", "upper"))
  expect_rp_error(map_tau2_interval(list()), regexp = "rma.uni")
})

test_that("the upper end of a metafor interval is found beyond metafor's search limit", {
  skip_if_not_installed("metafor")
  fit <- metafor::rma(yi = eight_schools$estimate, sei = eight_schools$se, method = "REML")
  # confint() stops at tau2.max = 100 by default and marks the end as a bound.
  expect_identical(stats::confint(fit)$ub.sign, ">")
  wide <- stats::confint(fit, control = list(tau2.max = 1e5))$random["tau^2", "ci.ub"]
  expect_gt(wide, 270)
  m <- map_tau2_interval(fit)
  expect_equal(m$upper[[1]], reliability_profile(se = eight_schools$se, sigma_tau = sqrt(wide))$stats[["I"]],
               tolerance = 1e-6)
})

test_that("a fit without a between-site variance is refused", {
  skip_if_not_installed("metafor")
  fit <- metafor::rma(yi = eight_schools$estimate, sei = eight_schools$se, method = "FE")
  expect_rp_error(map_tau2_interval(fit), class = "rp_scope_error", regexp = "not a random-effects model")
})

test_that("map_tau2_interval() checks the interval", {
  expect_rp_error(map_tau2_interval(se = c(1, 2)), regexp = "`interval` is missing")
  expect_rp_error(map_tau2_interval(se = c(1, 2), interval = c(-1, 2)), regexp = "at least 0")
  expect_rp_error(map_tau2_interval(se = c(1, 2), interval = c(1, 2, 3)), regexp = "named")
  expect_rp_error(map_tau2_interval(se = c(1, 2), interval = c(lower = 5, upper = 1)), regexp = "above its upper end")
  expect_rp_error(map_tau2_interval(se = c(1, 2), interval = c(a = 1, b = 5)), regexp = "`lower` and `upper`")
  ends <- map_tau2_interval(se = c(1, 2), interval = c(upper = 5, lower = 1))
  expect_equal(ends$lower[[1]], map_tau2_interval(se = c(1, 2), interval = c(1, 5))$lower[[1]])
})

test_that("the REML interval of the OCRS passing outcome maps to its stored values", {
  skip_if_not_installed("metafor")
  g <- ocrs_analyses()$pooled
  fit <- metafor::rma(yi = g$y, vi = g$v, method = "REML")
  m <- map_tau2_interval(fit)
  expect_equal(round(fit$tau2, 5), 0.03213)
  ci <- stats::confint(fit)$random
  expect_equal(round(ci["tau^2", c("ci.lb", "ci.ub")], 4), c(ci.lb = 0.0142, ci.ub = 0.0639))
  expect_equal(round(m$estimate[[1]], 3), 0.663)
  expect_equal(round(c(m$lower[[1]], m$upper[[1]]), 3), c(0.465, 0.796))
  expect_equal(round(c(m$lower[[4]], m$upper[[4]]), 4), c(0.0196, 0.0289))
  expect_equal(round(fit$I2 / 100, 3), 0.686)
})

# A direct generalized-least-squares computation of the Q-profile interval,
# for comparison. Its inputs are on the scale of the example.
q_profile_reference <- function(y, v, X, level = 0.95) {
  Q <- function(tau2) {
    W <- diag(1 / (v + tau2))
    beta <- solve(crossprod(X, W %*% X), crossprod(X, W %*% y))
    residual <- y - drop(X %*% beta)
    drop(crossprod(residual, W %*% residual))
  }
  df <- length(y) - ncol(X)
  critical <- stats::qchisq(c((1 + level) / 2, (1 - level) / 2), df)
  vapply(critical, function(target) {
    if (Q(0) <= target) return(0)
    upper <- max(1, stats::var(y))
    while (Q(upper) > target) upper <- upper * 10
    stats::uniroot(function(tau2) Q(tau2) - target, c(0, upper), tol = 1e-11)$root
  }, 0)
}

test_that("Q-profile endpoints and mapped intervals preserve effect units", {
  skip_if_not_installed("metafor")
  y <- eight_schools$estimate
  v <- eight_schools$se^2
  expected <- q_profile_reference(y, v, matrix(1, length(y), 1))
  results <- lapply(c(1, 1e-4, 1e-6), function(a) {
    fit <- metafor::rma(yi = y * a, vi = v * a^2, method = "REML")
    mapped <- map_tau2_interval(fit)
    ends <- unname(attr(mapped, "interval")[c("lower", "upper")])
    expect_equal(ends / a^2, expected, tolerance = 1e-9)
    expect_equal(mapped$upper[1], 0.64736618, tolerance = 1e-7)
    mapped[c("lower", "upper")]
  })
  expect_equal(results[[2]], results[[1]], tolerance = 1e-10)
  expect_equal(results[[3]], results[[1]], tolerance = 1e-10)
})

test_that("Q-profile intervals use the full meta-regression mean model", {
  skip_if_not_installed("metafor")
  x <- seq(-1, 1, length.out = 12)
  y <- 0.3 + 1.8 * x + c(-0.6, 0.9, -0.7, 1, -0.3, 0.5, -0.9, 0.4, -0.5, 0.8, -0.4, 0.2)
  v <- seq(0.03, 0.12, length.out = 12)
  for (X in list(cbind(1, x), matrix(x, ncol = 1))) {
    expected <- q_profile_reference(y, v, X, level = 0.9)
    fit <- metafor::rma(yi = y, vi = v, mods = X, intercept = FALSE, method = "REML")
    mapped <- map_tau2_interval(fit, level = 0.9)
    expect_equal(unname(attr(mapped, "interval")[c("lower", "upper")]), expected, tolerance = 1e-9)
    # Positive endpoints satisfy the chi-squared equations independently.
    for (a in c(1e-4, 1e-6)) {
      scaled <- metafor::rma(yi = y * a, vi = v * a^2, mods = X,
                              intercept = FALSE, method = "REML")
      ends <- attr(map_tau2_interval(scaled, level = 0.9), "interval")
      expect_equal(unname(ends[c("lower", "upper")]) / a^2, expected, tolerance = 1e-9)
    }
  }
})

test_that("generalized-Q estimators are explicitly mapped with Q-profile intervals", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("CompQuadForm")
  y <- eight_schools$estimate
  v <- eight_schools$se^2
  expected <- q_profile_reference(y, v, matrix(1, length(y), 1))
  for (method in c("GENQ", "GENQM")) {
    fit <- metafor::rma(yi = y, vi = v, method = method, weights = 1 / v)
    mapped <- map_tau2_interval(fit)
    expect_match(attr(mapped, "method"), "Q-profile")
    expect_equal(unname(attr(mapped, "interval")[c("lower", "upper")]), expected, tolerance = 1e-9)
    qp <- stats::confint(fit, type = "qp", control = list(tau2.max = 1e5, tol = 1e-9))
    expect_equal(unname(qp$random["tau^2", c("ci.lb", "ci.ub")]), expected, tolerance = 1e-8)
    genq <- stats::confint(fit, type = "genq")
    expect_gt(abs(mapped$upper[1] - reliability_profile(se2 = v,
                sigma_tau = sqrt(genq$random["tau^2", "ci.ub"]))$stats[["I"]]), 0.01)
  }
})

test_that("homogeneous mapped concentration intervals are a certain singleton", {
  mapped <- map_tau2_interval(se = rep(1, 4), interval = c(1, 2))
  row <- mapped[mapped$functional == "J_eff", ]
  expect_equal(row$lower, 4)
  expect_equal(row$upper, 4)
  expect_equal(reliability_profile(se = rep(1, 4), sigma_tau = 3)$stats[["J_eff"]], 4)
})

test_that("Q-profile ends are found when the sampling variances span many orders of magnitude", {
  skip_if_not_installed("metafor")
  withr::local_seed(57)
  vi <- 10^stats::runif(40, -7, 7)
  yi <- stats::rnorm(40, 0, sqrt(vi + 1e-6))
  fit <- suppressWarnings(metafor::rma(yi, vi, method = "DL"))
  ends <- attr(map_tau2_interval(fit, level = 0.9), "interval")
  ci <- suppressWarnings(stats::confint(fit, type = "QP", level = 90,
    control = list(tol = 1e-14, maxiter = 2000, tau2.max = 1e20)))$random
  expect_equal(unname(ends[c("lower", "upper")]), unname(ci["tau^2", c("ci.lb", "ci.ub")]), tolerance = 1e-6)
})
