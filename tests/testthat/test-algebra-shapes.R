test_that("shaped profiles have exactly the asked I and V", {
  withr::local_seed(71)
  shapes <- c("even", "bimodal", "concentrated_precise", "concentrated_noisy")
  for (i in 1:200) {
    J <- sample(2:100, 1)
    I <- runif(1, 0.01, 0.99)
    V <- runif(1, 0, 4)
    sigma_tau <- exp(runif(1, -3, 1))
    shape <- sample(shapes, 1)
    v <- shape_profile(J, I, V, shape, sigma_tau = sigma_tau)
    S <- site_reliability(se2 = v, sigma_tau = sigma_tau)
    expect_length(v, J)
    expect_equal(kl_centroid(S), I, tolerance = 1e-12)
    expect_equal(logit_moments(S)[["V"]], V, tolerance = 1e-11)
    gm <- exp(mean(log(v)))
    expect_equal(gm, sigma_tau^2 * (1 - I) / I, tolerance = 1e-12)
  }
})

test_that("the concentrated shapes have ceiling(J / 5) sites in the small group", {
  for (J in c(2, 6, 7, 11, 25, 50, 100)) {
    k <- ceiling(J / 5)
    v <- shape_profile(J, 0.3, 1, "concentrated_precise")
    expect_identical(sum(v < min(v) * 1.0001), as.integer(k))
    v <- shape_profile(J, 0.3, 1, "concentrated_noisy")
    expect_identical(sum(v > max(v) * 0.9999), as.integer(k))
  }
})

test_that("the bimodal shape has ceiling(J / 2) precise sites", {
  for (J in c(2, 7, 30, 31)) {
    v <- shape_profile(J, 0.5, 1, "bimodal")
    expect_identical(sum(v < exp(mean(log(v)))), as.integer(ceiling(J / 2)))
  }
})

test_that("total precision rises with V as in Appendix F of the paper", {
  # At V = 2 and J = 50 the paper reports the factor mean(exp(-x)): 2.35
  # (even), 2.18 (bimodal), 3.78 (concentrated precise), 1.63 (concentrated
  # noisy).
  factor <- function(shape) {
    v <- shape_profile(50, 0.3, 2, shape)
    mean(exp(mean(log(v)) - log(v)))
  }
  expect_equal(round(factor("even"), 2), 2.35)
  expect_equal(round(factor("bimodal"), 2), 2.18)
  expect_equal(round(factor("concentrated_precise"), 2), 3.78)
  expect_equal(round(factor("concentrated_noisy"), 2), 1.63)
})

test_that("V = 0 gives equal variances, and one site needs V = 0", {
  for (shape in c("even", "bimodal", "concentrated_precise", "concentrated_noisy")) {
    v <- shape_profile(10, 0.4, 0, shape, sigma_tau = 0.5)
    expect_equal(v, rep(0.25 * 0.6 / 0.4, 10), tolerance = 1e-14)
  }
  expect_equal(shape_profile(1, 0.4, 0), 1.5, tolerance = 1e-15)
  expect_rp_error(shape_profile(1, 0.4, 1), regexp = "single site")
})

test_that("shape_profile() checks its arguments", {
  expect_rp_error(shape_profile(0, 0.3, 1), regexp = "at least 1")
  expect_rp_error(shape_profile(10, 1, 1), regexp = "strictly between 0 and 1")
  expect_rp_error(shape_profile(10, 0.3, -1), regexp = "at least 0")
  expect_rp_error(shape_profile(10, 0.3, 1, "skewed"), regexp = "one of")
  expect_rp_error(shape_profile(10, 0.3, 1, sigma_tau = 0), regexp = "positive")
})
