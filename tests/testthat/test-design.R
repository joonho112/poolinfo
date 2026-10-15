# A design example: 40 sizes from a lognormal grid
# with a coefficient of variation of 0.5, scaled to a total of 2,000.
example_sizes <- function() {
  n <- stats::qlnorm(stats::ppoints(40), 0, sqrt(log(1.25)))
  round(n / sum(n) * 2000)
}

test_that("neyman_variance() is kappa / n with the kappa of multisiteDGP", {
  expect_equal(neyman_variance(c(20, 50)), 4 / c(20, 50))
  expect_equal(neyman_variance(100, p = 0.6, R2 = 0.4, var_outcome = 2),
               (1 / 0.6 + 1 / 0.4) * 0.6 * 2 / 100)
  expect_equal(neyman_variance(c(10, 20), p = c(0.5, 0.25)), c(4 / 10, (4 + 4 / 3) / 20))
  expect_rp_error(neyman_variance(10, p = 1), regexp = "strictly between 0 and 1")
  expect_rp_error(neyman_variance(c(10, 20, 30), p = c(0.5, 0.5)), regexp = "one value for each")
  skip_if_not_installed("multisiteDGP")
  expect_equal(neyman_variance(50, p = 0.3, R2 = 0.2, var_outcome = 1.5),
               multisiteDGP::compute_kappa(p = 0.3, R2 = 0.2, var_outcome = 1.5) / 50)
})

test_that("the design example has the expected profile", {
  d <- reliability_design(example_sizes(), sigma_tau = 0.2, kappa = 4)
  expect_s3_class(d, "reliability_design")
  expect_identical(range(d$n), c(16, 129))
  expect_equal(round(d$local$cv, 2), 0.48)
  expect_equal(round(exp(mean(log(d$n))), 1), 44.9)
  expect_equal(round(d$planned$stats[["I"]], 3), 0.310)
  expect_equal(round(d$balanced$stats[["I"]], 3), 0.333)
  expect_equal(round(d$planned$stats[["J_eff"]], 1), 36.6)
  expect_identical(d$balanced$stats[["H"]], 0)
  expect_equal(d$efficiency, exp(mean(log(d$n))) / 50)
  # The largest size is 2.6 times the mean, beyond the reach of the bound.
  expect_true(is.na(d$local$bound))
  expect_snapshot(print(d, use_cli = FALSE))
})

test_that("balance maximizes I and the mean reliability (Proposition 4)", {
  withr::local_seed(151)
  for (i in 1:2000) {
    J <- sample(2:30, 1)
    n <- stats::rexp(J) + 0.01
    n <- n / sum(n) * 100
    sigma <- exp(stats::runif(1, -2, 1))
    kappa <- exp(stats::runif(1, -1, 2))
    planned <- reliability_profile(se2 = kappa / n, sigma_tau = sigma)$stats
    balanced <- reliability_profile(se2 = rep(kappa * J / 100, J), sigma_tau = sigma)$stats
    expect_lte(planned[["I"]], balanced[["I"]] + 1e-14)
    expect_lte(planned[["mean"]], balanced[["mean"]] + 1e-14)
    expect_gte(planned[["benchmark_risk"]], balanced[["benchmark_risk"]] - 1e-14)
  }
})

test_that("the local approximation of the geometric mean is within its bound (Appendix E.2)", {
  withr::local_seed(152)
  for (i in 1:500) {
    J <- sample(2:50, 1)
    r <- stats::runif(1, 0.01, 0.9)
    u <- stats::runif(J, -r, r)
    u <- u - mean(u)
    if (max(abs(u)) >= 1) next
    n <- 100 * (1 + u)
    a <- .rp_size_approximation(n)
    expect_lte(abs(a$exact - a$leading), a$bound * (1 + 1e-12))
  }
})

test_that("unequal costs are compared with equal spending", {
  n <- c(20, 40, 60, 80)
  cost <- c(1, 1, 2, 2)
  d <- reliability_design(n, sigma_tau = 0.3, kappa = 4, cost = cost)
  B <- sum(cost * n)
  expect_equal(d$balanced_n, B / (4 * cost))
  expect_identical(d$balanced$source, "equal spending")
})

test_that("balanced_allocation() gives equal sizes, the bounded solution and equal spending", {
  expect_equal(as.vector(balanced_allocation(J = 4, N = 100, integer = FALSE)), rep(25, 4))
  n <- balanced_allocation(J = 7, N = 100)
  expect_identical(sum(n), 100)
  expect_lte(diff(range(n)), 1)
  b <- balanced_allocation(J = 5, N = 100, upper = c(10, 50, 50, 50, 50), integer = FALSE)
  expect_equal(as.vector(b), c(10, 22.5, 22.5, 22.5, 22.5))
  b2 <- balanced_allocation(J = 3, N = 30, lower = c(15, 1, 1), integer = FALSE)
  expect_equal(as.vector(b2), c(15, 7.5, 7.5))
  s <- balanced_allocation(J = 4, cost = c(1, 1, 2, 4), budget = 400, integer = FALSE)
  expect_equal(as.vector(s), 400 / (4 * c(1, 1, 2, 4)))
  si <- balanced_allocation(J = 4, cost = c(1, 1, 2, 4), budget = 403)
  expect_lte(attr(si, "spent"), 403)
  expect_true(all(si == round(si)))
  expect_rp_error(balanced_allocation(J = 3, N = 30, cost = 1, budget = 30), regexp = "not both")
  expect_rp_error(balanced_allocation(J = 3), regexp = "Give the total")
  expect_rp_error(balanced_allocation(J = 3, N = 10, lower = 5), regexp = "lower bounds")
  expect_rp_error(balanced_allocation(J = 3, N = 100, upper = 10), regexp = "upper bounds")
})

test_that("the integer allocation is the most even one within the bounds", {
  withr::local_seed(153)
  for (i in 1:200) {
    J <- sample(2:15, 1)
    N <- sample(J:(20 * J), 1)
    n <- balanced_allocation(J = J, N = N)
    expect_equal(sum(n), N)
    expect_lte(diff(range(n)), 1)
    # No transfer of one individual raises the product of the sizes.
    a <- which.max(n)
    b <- which.min(n)
    if (n[[a]] - n[[b]] >= 2) fail("a transfer would raise the product")
  }
})

test_that("integer allocations check feasibility before rounding", {
  expect_rp_error(balanced_allocation(2, N = 3, lower = 1.1), regexp = "lower bounds")
  expect_rp_error(balanced_allocation(2, N = 3, lower = 1.2, upper = 1.8),
                  regexp = "no positive integer")
  expect_rp_error(balanced_allocation(3, N = 2, lower = 0), regexp = "lower bounds")
  expect_rp_error(balanced_allocation(2, N = 3.5), regexp = "must be an integer")
  expect_rp_error(balanced_allocation(2, cost = c(1, 3), budget = 7, lower = 1.1),
                  regexp = "lower bounds")
  n <- balanced_allocation(2, N = 5, lower = c(1.1, 0), upper = c(3.9, 2.9))
  expect_equal(as.vector(n), c(3, 2))
  expect_equal(sum(n), 5)
  expect_true(all(n >= c(1.1, 0) & n <= c(3.9, 2.9)))
  expect_true(all(n >= 1 & n == floor(n)))
  # Fractional sizes remain available on the continuous path.
  expect_equal(as.vector(balanced_allocation(2, N = 3.5, lower = 1.1, integer = FALSE)),
               c(1.75, 1.75))
})

test_that("a budget is a ceiling even when upper bounds prevent full spending", {
  for (integer in c(TRUE, FALSE)) {
    n <- balanced_allocation(2, cost = c(1, 3), budget = 100,
                             upper = c(3, 2), integer = integer)
    expect_equal(as.vector(n), c(3, 2))
    expect_equal(attr(n, "spent"), 9)
    expect_lte(sum(c(1, 3) * n), 100)
  }
  expect_rp_error(balanced_allocation(2, N = 100, upper = c(3, 2)), regexp = "upper bounds")
})

test_that("equal-cost integer solutions agree with independent enumeration", {
  # Exhaust all positive bounded allocations; the enumeration does not use the
  # continuous solution or the rounding algorithm.
  candidates <- expand.grid(n1 = 1:6, n2 = 2:5, n3 = 1:4)
  candidates <- as.matrix(candidates)
  for (N in 4:15) {
    feasible <- candidates[rowSums(candidates) == N, , drop = FALSE]
    n <- balanced_allocation(3, N = N, lower = c(0.1, 1.1, 0),
                             upper = c(6.9, 5.4, 4.8))
    optimum <- max(apply(feasible, 1, prod))
    expect_equal(prod(n), optimum)
    expect_equal(sum(n), N)
    expect_true(all(n >= c(1, 2, 1) & n <= c(6, 5, 4)))
    # Equal non-unit costs give the same optimum under a budget ceiling.
    b <- balanced_allocation(3, cost = 2, budget = 2 * N + 0.5,
                             lower = c(0.1, 1.1, 0), upper = c(6.9, 5.4, 4.8))
    expect_equal(prod(b), optimum)
    expect_lte(attr(b, "spent"), 2 * N + 0.5)
  }
})

test_that("unequal-cost integer rounding is a feasible heuristic", {
  cost <- c(1, 3)
  n <- balanced_allocation(2, cost = cost, budget = 10)
  candidates <- as.matrix(expand.grid(n1 = 1:10, n2 = 1:3))
  feasible <- candidates[as.vector(candidates %*% cost) <= 10, , drop = FALSE]
  products <- apply(feasible, 1, prod)
  expect_equal(as.vector(n), c(7, 1))
  expect_equal(max(products), 8)
  expect_lt(prod(n), max(products))
  expect_true(all(n >= 1 & n == floor(n)))
  expect_lte(sum(cost * n), 10)
})

test_that("solve_design() inverts the balanced informativeness", {
  s <- solve_design(target_I = 0.5, sigma_tau = c(0.1, 0.2), kappa = 4, J = 40)
  expect_equal(s$n_per_site, 4 / c(0.01, 0.04))
  expect_equal(0.04 / (0.04 + 4 / s$n_per_site[[2]]), 0.5)
  t <- solve_design(target_I = 0.5, sigma_tau = 0.2, kappa = 4, N = 4000)
  expect_identical(t$J, 40)
  expect_gte(t$I, 0.5)
  expect_lt(0.04 / (0.04 + 4 * 41 / 4000), 0.5)
  expect_true(is.na(solve_design(0.99, 0.01, 4, N = 10)$J))
  expect_rp_error(solve_design(0.5, 0.2, 4), regexp = "not both and not neither")
})

test_that("design scenarios and inverse sizes preserve representable extreme results", {
  d <- reliability_design(c(20, 40), sigma_tau = c(1e200, 1e-160), kappa = 4)
  expect_identical(d$scenarios$I_planned[[1]], 1)
  expect_identical(d$scenarios$I_balanced[[1]], 1)
  expect_true(all(is.finite(as.matrix(d$scenarios))))
  expect_gt(d$scenarios$I_planned[[2]], 0)
  expect_equal(d$scenarios$I_planned[[2]] / exp(2 * log(1e-160) - mean(log(c(4 / 20, 4 / 40)))),
               1, tolerance = 1e-12)
  s <- solve_design(0.5, sigma_tau = 1e200, kappa = 1e300, J = 2)
  expect_equal(s$n_per_site / 1e-100, 1, tolerance = 1e-12)
  expect_equal(s$N / 2e-100, 1, tolerance = 1e-12)
  t <- solve_design(0.5, sigma_tau = 1e200, kappa = 1e300, N = 2e-100)
  expect_identical(t$J, 2)
  expect_equal(t$n_per_site / 1e-100, 1, tolerance = 1e-12)
  expect_equal(t$I, 0.5, tolerance = 1e-12)
  low <- solve_design(0.5, sigma_tau = 1e-200, kappa = 1e-300, J = 2)
  expect_equal(low$n_per_site / 1e100, 1, tolerance = 1e-12)
  expect_equal(low$N / 2e100, 1, tolerance = 1e-12)
  expect_rp_error(solve_design(0.5, 1e-300, 1e300, J = 2), regexp = "numeric range")
})

test_that("imbalance_ratio() gives the two documented values", {
  expect_equal(round(imbalance_ratio(0.75), 2), 0.80)
  expect_equal(round(imbalance_ratio(0.75, "local"), 3), 0.755)
  expect_rp_error(imbalance_ratio(-1), regexp = "at least 0")
})

test_that("a design works with the path, the profile and the site table", {
  d <- reliability_design(c(20, 40, 60), sigma_tau = 0.3, kappa = 4)
  expect_s3_class(reliability_profile(d), "reliability_profile")
  expect_true("sigma_tau" %in% names(profile_path(d, n = 5)))
  tab <- site_table(d)
  expect_named(tab, c("site", "n", "kappa", "n_balanced", "se", "reliability", "kl_term", "share"))
  ps <- profile_summary(d)
  expect_identical(unique(ps$allocation), c("planned", "balanced"))
})

test_that("a single site gets the whole total or what the budget buys", {
  expect_equal(as.numeric(balanced_allocation(1, N = 25)), 25)
  n <- balanced_allocation(1, cost = 0.57, budget = 27.9)
  expect_equal(as.numeric(n), floor(27.9 / 0.57))
  expect_lte(attr(n, "spent"), 27.9)
})
