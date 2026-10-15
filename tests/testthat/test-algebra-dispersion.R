test_that("H is the mean Bernoulli divergence from I, and the mean of the site terms", {
  withr::local_seed(31)
  for (i in 1:200) {
    p <- random_profile(log_ratio = runif(1, 0, 4), sigma_tau = exp(runif(1, -2, 1)))
    H <- kl_dispersion(p$S)
    # The naive formula loses relative precision when H is small, so the two
    # are compared on the absolute scale.
    expect_lt(abs(H - naive_H(p$S)), 1e-12)
    expect_equal(mean(kl_contributions(p$S)), H, tolerance = 1e-14)
    expect_true(all(kl_contributions(p$S) >= 0))
  }
})

test_that("H is log(5/4) for the reliabilities 0.2 and 0.8, and the right minimum differs", {
  expect_equal(kl_dispersion(c(0.2, 0.8)), log(5 / 4), tolerance = 1e-15)
  expect_equal(bregman_dispersion(c(0.2, 0.8), "kl", "first"), log(5 / 4), tolerance = 1e-15)
  right <- bregman_dispersion(c(0.2, 0.8), "kl", "second")
  expect_equal(right, (0.2 * log(0.2) + 0.8 * log(0.8)) - log(0.5), tolerance = 1e-14)
  expect_equal(round(right, 6), 0.192745)
  expect_equal(kl_centroid(c(0.2, 0.8)), 0.5)
  expect_equal(mean(c(0.2, 0.8)), 0.5)
})

test_that("the divergence identity of Theorem 2 holds for every p", {
  withr::local_seed(32)
  for (i in 1:100) {
    S <- random_profile(log_ratio = runif(1, 0, 4))$S
    I <- kl_centroid(S)
    H <- kl_dispersion(S)
    p <- runif(1, 0.01, 0.99)
    expect_equal(mean(naive_kl(p, S)), naive_kl(p, I) + H, tolerance = 1e-10)
  }
})

test_that("site terms are given with the names of the reliabilities", {
  S <- c(a = 0.2, b = 0.5, c = 0.9)
  expect_named(kl_contributions(S), c("a", "b", "c"))
  expect_null(names(kl_contributions(unname(S))))
})

test_that("the logit moments divide by J and do not depend on sigma_tau", {
  v <- c(0.01, 0.02, 0.05, 0.1)
  m1 <- logit_moments(site_reliability(se2 = v, sigma_tau = 0.1))
  m2 <- logit_moments(site_reliability(se2 = v, sigma_tau = 3))
  expect_equal(m1, m2, tolerance = 1e-12)
  d <- mean(log(v)) - log(v)
  expect_equal(m1, c(V = mean(d^2), rho3 = mean(abs(d)^3)), tolerance = 1e-12)
})

test_that("the leading terms and bounds are those of Proposition 2", {
  withr::local_seed(33)
  S <- random_profile(J = 20)$S
  lt <- local_terms(S)
  expect_identical(lt$quantity, c("H", "gap"))
  I <- kl_centroid(S)
  mom <- logit_moments(S)
  expect_equal(lt$leading, c(0.5 * I * (1 - I) * mom[["V"]],
                             0.5 * I * (1 - I) * (1 - 2 * I) * mom[["V"]]), tolerance = 1e-14)
  expect_equal(lt$bound, c(sqrt(3) / 108, 1 / 48) * mom[["rho3"]], tolerance = 1e-15)
  expect_equal(lt$exact, c(kl_dispersion(S), mean(S) - I), tolerance = 1e-12)
})

test_that("the remainders are within the bounds of Proposition 2(ii)", {
  withr::local_seed(34)
  for (i in 1:500) {
    J <- sample(2:80, 1)
    d <- rnorm(J) * runif(1, 0.01, 3)
    if (runif(1) < 0.3) d <- d^3
    d <- d - mean(d)
    # Keep every logit within about 35 of 0, where the reliabilities are
    # still distinguishable from 0 and 1 as doubles.
    d <- d * min(1, 25 / max(abs(d)))
    m <- max(-10, min(10, rnorm(1, 0, 3)))
    lt <- local_terms(profile_from_logits(m, d))
    remainder <- abs(lt$exact - lt$leading)
    expect_true(all(remainder <= lt$bound * (1 + 1e-9) + 1e-14))
    # Where the leading term decides the sign, the exact value has that sign.
    decided <- which(lt$holds)
    expect_true(all(sign(lt$exact[decided]) == sign(lt$leading[decided])))
  }
})

test_that("the skewed example of the paper reverses the sign of the leading term (Appendix C)", {
  S <- profile_from_logits(stats::qlogis(0.3), paper_counterexample_d)
  expect_equal(kl_centroid(S), 0.3, tolerance = 1e-14)
  expect_equal(logit_moments(S)[["V"]], 2, tolerance = 1e-13)
  expect_equal(round(mean(S), 10), 0.2732468357)
  lt <- local_terms(S)
  expect_lt(lt$exact[[2L]], 0)
  expect_gt(lt$leading[[2L]], 0)
  expect_false(lt$holds[[2L]])
  mirror <- profile_from_logits(-stats::qlogis(0.3), -paper_counterexample_d)
  expect_equal(kl_centroid(mirror), 0.7, tolerance = 1e-14)
  expect_equal(round(mean(mirror), 10), 0.7267531643)
})

test_that("a constant profile has no dispersion and undecided signs", {
  lt <- local_terms(rep(0.4, 6))
  expect_identical(lt$exact, c(0, 0))
  expect_identical(lt$leading, c(0, 0))
  expect_identical(lt$holds, c(NA, NA))
  expect_identical(kl_dispersion(0.4), 0)
  expect_identical(unname(logit_moments(0.4)), c(0, 0))
})
