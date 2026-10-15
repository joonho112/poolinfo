se6 <- c(0.10, 0.12, 0.15, 0.20, 0.25, 0.30)
y6 <- c(0.30, 0.10, 0.25, -0.05, 0.40, 0.00)

test_that("a profile holds the reliabilities and the diagnostics of the algebra", {
  prof <- reliability_profile(se = se6, sigma_tau = 0.2)
  expect_s3_class(prof, "reliability_profile")
  S <- site_reliability(se = se6, sigma_tau = 0.2)
  expect_equal(prof$S, unname(S), tolerance = 1e-15)
  st <- prof$stats
  expect_equal(st[["I"]], kl_centroid(S), tolerance = 1e-14)
  expect_equal(st[["mean"]], mean(S), tolerance = 1e-14)
  expect_equal(st[["harmonic"]], reliability_centroids(S)[["harmonic"]], tolerance = 1e-14)
  expect_equal(st[["H"]], kl_dispersion(S), tolerance = 1e-14)
  expect_equal(st[c("V", "rho3")], logit_moments(S), tolerance = 1e-12)
  expect_equal(st[["J_eff"]], concentration_count(S), tolerance = 1e-13)
  expect_equal(st[["gap"]], mean(S) - kl_centroid(S), tolerance = 1e-14)
  expect_equal(st[["benchmark_risk"]], benchmark_risk(se = se6, sigma_tau = 0.2), tolerance = 1e-14)
  expect_equal(st[["trace"]], sum(S), tolerance = 1e-14)
  expect_equal(st[["v_typ"]], typical_variance(se = se6))
  expect_equal(st[["R_typ"]], 0.04 / (0.04 + typical_variance(se = se6)))
  expect_equal(prof$local, local_terms(S), tolerance = 1e-13)
  expect_identical(prof$source, "supplied")
  expect_identical(prof$site, as.character(1:6))
})

test_that("se and se2 give the same profile, and names become site labels", {
  a <- reliability_profile(se = se6, sigma_tau = 0.2)
  b <- reliability_profile(se2 = se6^2, sigma_tau = 0.2)
  expect_equal(a$stats, b$stats, tolerance = 1e-15)
  named <- reliability_profile(se = stats::setNames(se6, letters[1:6]), sigma_tau = 0.2)
  expect_identical(named$site, letters[1:6])
  labelled <- reliability_profile(se = se6, sigma_tau = 0.2, site = LETTERS[1:6])
  expect_identical(labelled$site, LETTERS[1:6])
})

test_that("an unnamed vector is refused, as are objects without a method", {
  err <- expect_rp_error(reliability_profile(se6, sigma_tau = 0.2), regexp = "ambiguous")
  expect_identical(error_call_name(err), "reliability_profile")
  expect_rp_error(reliability_profile(list(1, 2), sigma_tau = 0.2), regexp = "no method")
  expect_rp_error(reliability_profile(se = se6), regexp = "`sigma_tau` is missing")
  expect_rp_error(reliability_profile(se = se6, sigma_tau = 0.2, level = 0.9),
                  regexp = "no further arguments")
  expect_rp_error(reliability_profile(se = se6, sigma_tau = 0.2, estimate = 1:3),
                  regexp = "one number for each of the 6 sites")
  expect_rp_error(reliability_profile(se = se6, sigma_tau = 0.2, site = rep("a", 6)),
                  regexp = "distinct")
})

test_that("a data frame's columns are found by name, before objects of the same name", {
  dat <- data.frame(school = letters[1:6], y = y6, s = se6)
  prof <- reliability_profile(dat, se = s, estimate = y, site = school, sigma_tau = 0.2)
  expect_identical(prof$site, letters[1:6])
  expect_identical(prof$estimate, y6)
  ref <- reliability_profile(se = se6, estimate = y6, site = letters[1:6], sigma_tau = 0.2)
  expect_equal(prof$stats, ref$stats)
  s <- rep(1, 6)
  expect_equal(reliability_profile(dat, se = s, sigma_tau = 0.2)$se2, se6^2)
  expect_equal(reliability_profile(dat, se = .env$s, sigma_tau = 0.2)$se2, rep(1, 6))
  expect_equal(reliability_profile(dat, se2 = s^2, sigma_tau = 0.2)$se2, se6^2)
  err <- expect_rp_error(reliability_profile(dat, se = sd_col, sigma_tau = 0.2),
                         regexp = "could not be found")
  expect_match(conditionMessage(err), "`school`, `y`, `s`", fixed = TRUE)
})

test_that("one site gives a profile without dispersion", {
  prof <- reliability_profile(se = 0.2, sigma_tau = 0.2)
  expect_equal(prof$stats[["I"]], 0.5)
  expect_identical(prof$stats[["H"]], 0)
  expect_identical(prof$stats[["J_eff"]], 1)
  expect_true(is.na(prof$stats[["v_typ"]]))
  expect_identical(prof$local$holds, c(NA, NA))
})

test_that("sigma_tau = 0 gives the limits and a message", {
  expect_message(prof <- reliability_profile(se = se6, sigma_tau = 0, estimate = y6),
                 class = "rp_message_sigma_zero")
  st <- prof$stats
  expect_identical(unname(st[c("I", "mean", "harmonic", "H", "gap", "benchmark_risk", "trace")]),
                   rep(0, 7))
  w <- 1 / se6^2
  expect_equal(st[["J_eff"]], sum(w)^2 / sum(w^2))
  # The limit agrees with a profile at a tiny sigma_tau.
  tiny <- reliability_profile(se = se6, sigma_tau = 1e-9)
  expect_equal(tiny$stats[["J_eff"]], st[["J_eff"]], tolerance = 1e-12)
  expect_lt(abs(tiny$stats[["I"]]), 1e-15)
  # Complete pooling: every effect's posterior mean is the weighted mean.
  tab <- site_table(prof)
  expect_equal(tab$effect, rep(sum(w * y6) / sum(w), 6))
  expect_equal(tab$effect_sd, rep(sqrt(1 / sum(w)), 6))
  expect_equal(tab$share, w / sum(w))
})

test_that("update() recomputes the profile at another sigma_tau", {
  prof <- reliability_profile(se = se6, sigma_tau = 0.2, scale_label = "effect size")
  up <- update(prof, sigma_tau = 0.1)
  expect_equal(up$stats, reliability_profile(se = se6, sigma_tau = 0.1, scale_label = "effect size")$stats)
  expect_identical(up$labels$scale_label, "effect size")
  expect_identical(update(prof)$stats, prof$stats)
  expect_rp_error(update(prof, sigma = 0.1), regexp = "no further arguments")
})

test_that("the site table has one row per site and its columns add up", {
  prof <- reliability_profile(se = se6, sigma_tau = 0.2, estimate = y6)
  tab <- site_table(prof)
  expect_named(tab, c("site", "estimate", "se", "reliability", "kl_term", "share",
                      "effect", "effect_sd"))
  expect_identical(nrow(tab), 6L)
  expect_equal(sum(tab$share), 1)
  expect_equal(mean(tab$kl_term), prof$stats[["H"]], tolerance = 1e-14)
  expect_equal(1 / sum(tab$share^2), prof$stats[["J_eff"]], tolerance = 1e-12)
  # Posterior mean and variance at the fixed sigma_tau, by direct formulas.
  a <- 0.04
  w <- 1 / (a + se6^2)
  mu <- sum(w * y6) / sum(w)
  S <- a / (a + se6^2)
  expect_equal(tab$effect, S * y6 + (1 - S) * mu, tolerance = 1e-14)
  expect_equal(tab$effect_sd^2, a * (1 - S) + (1 - S)^2 / sum(w), tolerance = 1e-14)
  expect_identical(as.data.frame(prof), tab)
  expect_named(site_table(reliability_profile(se = se6, sigma_tau = 0.2)),
               c("site", "se", "reliability", "kl_term", "share"))
})

test_that("profile_summary() gives one row per quantity", {
  prof <- reliability_profile(se = se6, sigma_tau = 0.2)
  ps <- profile_summary(prof)
  expect_identical(ps$quantity, c("I", "mean", "harmonic", "H", "V", "J_eff", "sigma_tau"))
  expect_equal(ps$mean[ps$quantity == "I"], prof$stats[["I"]])
  expect_equal(ps$mean[ps$quantity == "sigma_tau"], 0.2)
  expect_true(all(is.na(ps$sd)))
  expect_rp_error(profile_summary(1), regexp = "no method")
  expect_rp_error(site_table("a"), regexp = "no method")
})

test_that("the printed profile and its summary keep their words", {
  prof <- reliability_profile(se = se6, sigma_tau = 0.2, scale_label = "standardized mean difference")
  expect_snapshot(print(prof, use_cli = FALSE))
  expect_snapshot(print(summary(prof), use_cli = FALSE))
  expect_snapshot(print(reliability_profile(se = 0.2, sigma_tau = 0.3), use_cli = FALSE))
  expect_snapshot(print(suppressMessages(reliability_profile(se = se6, sigma_tau = 0)), use_cli = FALSE))
})

test_that("the profile of the OCRS passing outcome at sigma_tau = 0.188 matches its reference values", {
  g <- ocrs_analyses()$pooled
  p <- ocrs_profile_reference()
  prof <- reliability_profile(se2 = g$v, sigma_tau = 0.188)
  st <- prof$stats
  expect_equal(st[["I"]], p[["I"]], tolerance = 1e-12)
  expect_equal(st[["mean"]], p[["Sbar"]], tolerance = 1e-12)
  expect_equal(st[["harmonic"]], p[["HM"]], tolerance = 1e-12)
  expect_equal(st[["H"]], p[["H"]], tolerance = 1e-12)
  expect_equal(st[["V"]], p[["V"]], tolerance = 1e-12)
  expect_equal(st[["rho3"]], p[["rho3"]], tolerance = 1e-12)
  expect_equal(st[["J_eff"]], p[["Jeff"]], tolerance = 1e-12)
  expect_equal(st[["gap"]], p[["gap"]], tolerance = 1e-12)
  expect_equal(st[["benchmark_risk"]], p[["benchmark_risk"]], tolerance = 1e-12)
  expect_equal(st[["trace"]], p[["trace"]], tolerance = 1e-12)
  expect_equal(st[["v_typ"]], p[["v_typ"]], tolerance = 1e-12)
  expect_equal(st[["R_typ"]], p[["R_typ"]], tolerance = 1e-12)
  expect_equal(prof$local$leading, unname(p[c("lead_H", "lead_gap")]), tolerance = 1e-12)
  expect_equal(prof$local$bound, unname(p[c("bound_H", "bound_gap")]), tolerance = 1e-12)
  expect_equal(range(prof$S), unname(p[c("S_min", "S_max")]), tolerance = 1e-12)
})
