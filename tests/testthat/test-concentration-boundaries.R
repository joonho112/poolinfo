test_that("a homogeneous concentration count has strict point-mass probabilities", {
  fit <- reliability_posterior(c(-2, -1, 0, 1, 2), se = rep(1, 5))
  expect_identical(prob_above(fit, 5, "J_eff"), 0)
  expect_identical(prob_above(fit, 4.9, "J_eff"), 1)
  expect_identical(prob_above(fit, 1, "J_eff"), 1)
  expect_rp_error(functional_density(fit, "J_eff", n = 10), "rp_scope_error",
                  regexp = "point mass")
})

test_that("an unresolved nonconstant concentration density explains its precision limit", {
  v <- 1 + 1e-7 * c(-2, -1, 0, 1, 2)
  fit <- reliability_posterior(c(-2, -1, 0, 1, 2), se2 = v)
  expect_gt(diff(range(fit$data$se2)), 0)
  expect_rp_error(functional_density(fit, "J_eff", n = 10), "rp_numeric_error",
                  regexp = "precision")
  expect_identical(prob_above(fit, 5, "J_eff"), 0)
})

test_that("nearby nonconstant concentration profiles keep a continuous density", {
  for (spread in c(.05, .005)) {
    v <- 1 + spread * c(-2, -1, 0, 1, 2)
    fit <- reliability_posterior(c(-2, -1, 0, 1, 2), se2 = v)
    density <- functional_density(fit, "J_eff", n = 100)
    expect_true(all(is.finite(density$density) & density$density >= 0))
    expect_true(all(diff(density$value) > 0))
    expect_identical(prob_above(fit, 5, "J_eff"), 0)
  }
})

test_that("the concentration derivative agrees with an independent finite difference", {
  v <- c(.3, .8, 1.2, 2, 4)
  kernel <- .rp_kernel(seq_along(v), v)
  sigma <- c(.1, .4, 1, 3)
  jeff <- function(s) {
    S <- s^2 / (s^2 + v)
    sum(S)^2 / sum(S^2)
  }
  step <- sigma * 1e-4
  reference <- vapply(seq_along(sigma), function(i) {
    (jeff(sigma[[i]] + step[[i]]) - jeff(sigma[[i]] - step[[i]])) / (2 * step[[i]])
  }, 0)
  expect_equal(.rp_functional_derivative(kernel, sigma, "J_eff"), reference, tolerance = 1e-7)
})

test_that("the concentration derivative keeps its positive tiny weighted variance", {
  v <- 1 + 1e-7 * c(-2, -1, 0, 1, 2)
  kernel <- .rp_kernel(seq_along(v), v)
  sigma <- 1
  S <- sigma^2 / (sigma^2 + v)
  pairs <- utils::combn(seq_along(S), 2)
  pair_sum <- sum(S[pairs[1, ]] * S[pairs[2, ]] * (S[pairs[1, ]] - S[pairs[2, ]])^2)
  reference <- 4 * sum(S) * pair_sum / (sigma * sum(S^2)^2)
  actual <- .rp_functional_derivative(kernel, sigma, "J_eff")
  expect_gt(actual, 0)
  expect_equal(actual / reference, 1, tolerance = 1e-8)
})
