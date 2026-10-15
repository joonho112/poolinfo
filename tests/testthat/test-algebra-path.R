test_that("the path passes through the profile at its own I", {
  withr::local_seed(61)
  S <- random_profile(J = 15)$S
  I <- kl_centroid(S)
  path <- profile_path(S, I = I)
  expect_equal(path$H, kl_dispersion(S), tolerance = 1e-14)
  expect_equal(path$mean, mean(S), tolerance = 1e-14)
  expect_equal(path$gap, mean(S) - I, tolerance = 1e-14)
  lt <- local_terms(S)
  expect_equal(c(path$H_leading, path$gap_leading), lt$leading, tolerance = 1e-14)
  expect_equal(c(path$H_bound, path$gap_bound), lt$bound, tolerance = 1e-15)
})

test_that("dH/dm is the gap along the path (Proposition 2(i))", {
  withr::local_seed(62)
  for (i in 1:50) {
    S <- random_profile(log_ratio = runif(1, 0.5, 5))$S
    I <- runif(5, 0.02, 0.98)
    m <- stats::qlogis(I)
    h <- 1e-5
    up <- profile_path(S, I = stats::plogis(m + h))$H
    down <- profile_path(S, I = stats::plogis(m - h))$H
    expect_equal((up - down) / (2 * h), profile_path(S, I = I)$gap, tolerance = 1e-7)
  }
})

test_that("the default path has n points strictly inside (0, 1)", {
  path <- profile_path(c(0.2, 0.5, 0.7), n = 9)
  expect_identical(nrow(path), 9L)
  expect_equal(path$I, (1:9) / 10)
  expect_named(path, c("I", "logit_I", "H", "mean", "gap", "H_leading", "gap_leading",
                       "H_bound", "gap_bound"))
})

test_that("H rises while the gap is positive and falls after, with one turning point", {
  withr::local_seed(63)
  m <- seq(-40, 40, length.out = 1601)
  for (i in 1:300) {
    J <- sample(2:60, 1)
    d <- rnorm(J) * exp(runif(1, -3, 1.5))
    if (runif(1) < 0.5) d <- sign(d) * abs(d)^runif(1, 1, 3)
    d <- d - mean(d)
    if (mean(d^2) < 1e-14) next
    gap <- .rp_gap_at(m, d)
    # Exactly one change of sign, from positive to negative.
    signs <- sign(gap[gap != 0])
    expect_identical(sum(diff(signs) != 0), 1L)
    expect_gt(signs[[1L]], 0)
    expect_lt(signs[[length(signs)]], 0)
  }
})

test_that("the gap keeps its sign where both means are close to 1", {
  d <- c(-1e-4, 0, 1e-4)
  m <- 35
  gap <- .rp_gap_at(m, d)
  # mean(S) - I = -(1 - I) mean(expm1(-d)) + O(exp(-2m)) < 0.
  expected <- -stats::plogis(-m) * mean(stats::plogis(m + d) * expm1(-d))
  expect_lt(gap, 0)
  expect_equal(gap, expected, tolerance = 1e-12)
  # The plain difference has lost every digit.
  expect_identical(mean(stats::plogis(m + d)) - stats::plogis(m), 0)
})

test_that("profile_path() checks its arguments", {
  expect_rp_error(profile_path("a"), regexp = "site reliabilities, or a profile")
  expect_rp_error(profile_path(c(0.2, 0.5), I = c(0.2, 1)), regexp = "strictly between 0 and 1")
  expect_rp_error(profile_path(c(0.2, 0.5), n = 1), regexp = "at least 2")
  # n is not used when I is given.
  expect_identical(nrow(profile_path(c(0.2, 0.5), I = 0.5, n = 1)), 1L)
  expect_rp_error(profile_path(c(0.2, 1.2)), regexp = "strictly between 0 and 1")
})
