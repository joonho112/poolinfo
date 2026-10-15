test_that("the concentration count lies in (1, J] and is J only for a constant profile", {
  withr::local_seed(41)
  for (i in 1:300) {
    p <- random_profile()
    jeff <- concentration_count(p$S)
    expect_equal(jeff, naive_jeff(p$S), tolerance = 1e-13)
    expect_gt(jeff, 1)
    expect_lte(jeff, p$J * (1 + 1e-14))
  }
  expect_equal(concentration_count(rep(0.3, 7)), 7)
  expect_equal(concentration_count(0.3), 1)
})

test_that("the concentration count does not change with a common factor", {
  S <- c(0.05, 0.1, 0.2, 0.4)
  expect_equal(concentration_count(S), concentration_count(S * 2), tolerance = 1e-15)
  expect_equal(concentration_count(S), concentration_count(S / 1e6), tolerance = 1e-15)
})

test_that("uniformly low and uniformly high profiles both give J", {
  expect_equal(concentration_count(rep(0.01, 10)), 10)
  expect_equal(concentration_count(rep(0.99, 10)), 10)
})

test_that("the count rises with sigma_tau for fixed variances", {
  withr::local_seed(42)
  for (i in 1:100) {
    v <- exp(runif(sample(2:40, 1), -3, 3))
    centered <- .rp_center_v(v)
    m <- seq(-20, 20, length.out = 200)
    jeff <- .rp_jeff_at(m, centered$d)
    expect_true(all(diff(jeff) >= -1e-12 * jeff[-1]))
  }
})

test_that("the count tends to (sum 1/v)^2 / sum 1/v^2 as sigma_tau goes to 0", {
  v <- c(0.01, 0.03, 0.1, 0.5)
  centered <- .rp_center_v(v)
  limit <- sum(1 / v)^2 / sum(1 / v^2)
  expect_equal(.rp_jeff_at(-60, centered$d), limit, tolerance = 1e-12)
  # Far below the smallest double, where every S_j is 0 as a number.
  expect_equal(.rp_jeff_at(-1e4, centered$d), limit, tolerance = 1e-12)
})
