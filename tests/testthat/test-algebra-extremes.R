# Profiles at the edges of what doubles can hold: the logit of I far from 0,
# sampling variances that differ by many orders of magnitude, and profiles
# that are nearly constant. Every summary must be finite and in its range.

summaries_at <- function(m, d) {
  c(.rp_profile_stats(m, d), harmonic_check = .rp_harmonic_at(m, d))
}

test_that("summaries are finite and in range far into both tails", {
  withr::local_seed(81)
  for (i in 1:200) {
    J <- sample(2:200, 1)
    d <- rnorm(J) * exp(runif(1, -8, 3))
    d <- d - mean(d)
    for (m in c(-700, -60, -5, 0, 5, 60, 700)) {
      s <- .rp_profile_stats(m, d)
      expect_true(all(is.finite(s)), info = paste("m =", m))
      expect_true(s[["H"]] >= 0)
      expect_true(s[["J_eff"]] >= 1 - 1e-12 && s[["J_eff"]] <= J * (1 + 1e-12))
      expect_true(all(s[c("I", "mean", "harmonic")] >= 0 & s[c("I", "mean", "harmonic")] <= 1))
      expect_lte(s[["harmonic"]], min(s[["I"]], s[["mean"]]) + 1e-15)
    }
  }
})

test_that("H agrees with a high-precision value in the tails", {
  # For m -> -inf, KL(Ber(I) || Ber(S_j)) ~ I (e^d - 1 - d), so
  # H ~ I mean(e^d - 1 - d); for m -> +inf the mirror image.
  d <- c(-1.5, -0.5, 0.2, 1.8)
  d <- d - mean(d)
  for (m in c(-40, -100, -300)) {
    expected <- stats::plogis(m) * mean(expm1(d) - d)
    expect_gt(expected, 0)
    expect_equal(.rp_H_at(m, d) / expected, 1, tolerance = 1e-10)
    expect_equal(.rp_H_at(-m, -d) / expected, 1, tolerance = 1e-10)
  }
})

test_that("H of a nearly constant profile is I(1 - I) V / 2 to high relative precision", {
  for (scale in c(1e-3, 1e-5, 1e-7)) {
    d <- c(-1, -0.2, 0.4, 0.8) * scale
    d <- d - mean(d)
    V <- mean(d^2)
    for (m in c(-3, 0, 3)) {
      I <- stats::plogis(m)
      # The remainder is at most sqrt(3) / 108 rho3, of relative order scale.
      expected <- 0.5 * I * (1 - I) * V
      expect_gt(expected, 0)
      expect_equal(.rp_H_at(m, d) / expected, 1, tolerance = 10 * scale)
    }
  }
})

test_that("variance ratios up to 1e16 give finite summaries", {
  for (ratio in 10^c(4, 8, 12, 16)) {
    v <- c(1, sqrt(ratio), ratio)
    for (sigma_tau in c(1e-6, 1, 1e6, 1e12)) {
      S <- site_reliability(se2 = v, sigma_tau = sigma_tau)
      centered <- .rp_center_v(v)
      m <- 2 * log(sigma_tau) - centered$log_gm
      s <- .rp_profile_stats(m, centered$d)
      expect_true(all(is.finite(s)))
      expect_equal(s[["I"]], sigma_tau^2 / (sigma_tau^2 + exp(mean(log(v)))), tolerance = 1e-12)
    }
  }
})

test_that("the site terms keep their relative precision when I is close to 1", {
  d <- c(-2, 0, 2)
  m <- 30
  terms <- .rp_bern_kl(rep(m, 3), d)
  # KL(Ber(s(m)) || Ber(s(m + d))) for m -> inf is (1 - I)(e^-d - 1 + d).
  expected <- stats::plogis(-m) * (exp(-d) - 1 + d)
  nonzero <- expected > 0
  expect_true(any(nonzero))
  expect_equal(terms[nonzero] / expected[nonzero], rep(1, sum(nonzero)), tolerance = 1e-10)
  expect_identical(terms[!nonzero], rep(0, sum(!nonzero)))
})

test_that("a site term is infinite only when a reliability is 0 or 1", {
  expect_identical(.rp_bern_kl(0, -Inf), Inf)
  expect_identical(.rp_bern_kl(0, Inf), Inf)
  expect_true(is.finite(.rp_bern_kl(0, 700)))
  expect_true(is.finite(.rp_bern_kl(0, -700)))
})
