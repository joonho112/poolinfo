test_that("I is the logit mean and the geometric-mean form (Theorem 1)", {
  withr::local_seed(21)
  for (i in 1:200) {
    p <- random_profile()
    gm <- exp(mean(log(p$v)))
    I_gm <- p$sigma_tau^2 / (p$sigma_tau^2 + gm)
    expect_equal(kl_centroid(p$S), I_gm, tolerance = tol_identity)
  }
})

test_that("the harmonic centroid is at most I and at most the mean", {
  withr::local_seed(22)
  for (i in 1:300) {
    p <- random_profile()
    cc <- reliability_centroids(p$S)
    expect_named(cc, c("harmonic", "I", "mean"))
    expect_lte(cc[["harmonic"]], cc[["I"]] + 1e-15)
    expect_lte(cc[["harmonic"]], cc[["mean"]] + 1e-15)
    expect_equal(cc[["harmonic"]], 1 / mean(1 / p$S), tolerance = tol_identity)
    expect_equal(cc[["mean"]], mean(p$S), tolerance = 1e-15)
  }
})

test_that("I and the mean reliability have no general order", {
  below <- reliability_centroids(profile_from_logits(stats::qlogis(0.3), paper_counterexample_d))
  above <- reliability_centroids(profile_from_logits(stats::qlogis(0.2), c(-1, 0, 1)))
  expect_lt(below[["mean"]], below[["I"]])
  expect_gt(above[["mean"]], above[["I"]])
})

test_that("the three centroids agree for a constant profile", {
  cc <- reliability_centroids(rep(0.37, 5))
  expect_equal(unname(cc), rep(0.37, 3), tolerance = 1e-15)
})

test_that("Bregman centroids follow Table 2 of the paper", {
  withr::local_seed(23)
  p <- random_profile(J = 12)
  S <- p$S
  expect_equal(bregman_centroid(S), kl_centroid(S))
  expect_equal(bregman_centroid(S, "kl", "first"), kl_centroid(S))
  expect_equal(bregman_centroid(S, "squared"), mean(S))
  expect_equal(bregman_centroid(S, "burg"), 1 / mean(1 / S), tolerance = tol_identity)
  for (pot in c("kl", "squared", "burg")) {
    expect_identical(bregman_centroid(S, pot, "second"), mean(S))
  }
  expect_rp_error(bregman_centroid(S, "entropy"), regexp = "one of")
  expect_rp_error(bregman_centroid(S, side = "left"), regexp = "one of")
})

test_that("a first-argument centroid minimizes the average divergence", {
  withr::local_seed(24)
  S <- random_profile(J = 9)$S
  phi <- list(
    kl = function(x) x * log(x) + (1 - x) * log(1 - x),
    squared = function(x) x^2 / 2,
    burg = function(x) -log(x)
  )
  dphi <- list(kl = stats::qlogis, squared = function(x) x, burg = function(x) -1 / x)
  for (pot in names(phi)) {
    D <- function(x, y) phi[[pot]](x) - phi[[pot]](y) - dphi[[pot]](y) * (x - y)
    first <- stats::optimize(function(c) mean(D(c, S)), c(1e-9, 1 - 1e-9), tol = 1e-12)$minimum
    second <- stats::optimize(function(c) mean(D(S, c)), c(1e-9, 1 - 1e-9), tol = 1e-12)$minimum
    expect_equal(bregman_centroid(S, pot, "first"), first, tolerance = 1e-7)
    expect_equal(bregman_centroid(S, pot, "second"), second, tolerance = 1e-7)
    expect_equal(bregman_dispersion(S, pot, "first"), mean(D(first, S)), tolerance = 1e-9)
    expect_equal(bregman_dispersion(S, pot, "second"), mean(D(S, second)), tolerance = 1e-9)
  }
})

test_that("centroid functions refuse reliabilities of 0 or 1", {
  for (fn in list(kl_centroid, reliability_centroids, kl_dispersion, kl_contributions,
                  logit_moments, concentration_count, local_terms, bregman_centroid,
                  bregman_dispersion)) {
    expect_rp_error(fn(c(0.2, 1)), regexp = "strictly between 0 and 1")
    expect_rp_error(fn(c(0, 0.5)), regexp = "strictly between 0 and 1")
    expect_rp_error(fn("a"), regexp = "numeric vector")
  }
})
