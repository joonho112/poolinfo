test_that("simulated data give the profile at the true SD and a posterior from the estimates", {
  skip_if_not_installed("multisiteDGP")
  design <- multisiteDGP::multisitedgp_design(J = 30, sigma_tau = 0.2, seed = 1)
  d <- multisiteDGP::sim_multisite(design)
  prof <- reliability_profile(d)
  expect_identical(prof$source, "multisiteDGP design value")
  expect_identical(prof$sigma_tau, 0.2)
  expect_equal(prof$se2, d$se2_j)
  expect_equal(prof$estimate, d$tau_j_hat)
  fit <- reliability_posterior(d)
  ref <- reliability_posterior(d$tau_j_hat, se2 = d$se2_j, site = d$site_index)
  expect_equal(fit$summary, ref$summary)
  expect_rp_error(reliability_posterior(d, se = d$se_j), regexp = "leave out")
})

test_that("the informativeness agrees with multisiteDGP's compute_I() and mean_shrinkage()", {
  skip_if_not_installed("multisiteDGP")
  withr::local_seed(171)
  for (i in 1:50) {
    v <- exp(runif(sample(2:40, 1), -6, 0))
    s <- exp(runif(1, -3, 0))
    expect_equal(kl_centroid(site_reliability(se2 = v, sigma_tau = s)), multisiteDGP::compute_I(v, s),
                 tolerance = 1e-12)
  }
  d <- multisiteDGP::sim_multisite(multisiteDGP::multisitedgp_design(J = 25, sigma_tau = 0.15, seed = 2))
  shrink <- tryCatch(multisiteDGP::mean_shrinkage(d), error = function(e) NULL)
  skip_if(is.null(shrink), "mean_shrinkage() of this multisiteDGP version does not take simulated data")
  expect_equal(reliability_profile(d)$stats[["mean"]], as.numeric(shrink), tolerance = 1e-12)
})

test_that("a metafor fit gives a plug-in profile, a posterior and its moderators", {
  skip_if_not_installed("metafor")
  fit <- metafor::rma(yi = estimate, sei = se, slab = school, data = eight_schools)
  expect_message(prof <- reliability_profile(fit), class = "rp_message_sigma_zero")
  expect_match(prof$source, "metafor REML estimate")
  expect_identical(prof$site, LETTERS[1:8])
  post <- reliability_posterior(fit)
  ref <- reliability_posterior(estimate, se = se, site = school, data = eight_schools)
  expect_equal(post$summary, ref$summary)
  withr::local_seed(172)
  dd <- data.frame(yi = rnorm(15), vi = runif(15, 0.05, 0.3), x = rnorm(15))
  fm <- metafor::rma(yi, vi, mods = ~ x, data = dd)
  pm <- reliability_posterior(fm)
  expect_true(all(c("coef_intrcpt", "coef_x") %in% pm$summary$quantity))
  direct <- reliability_posterior(yi, se2 = vi, data = dd, mean_model = ~ x)
  expect_equal(pm$summary$mean[1:7], direct$summary$mean[1:7], tolerance = 1e-12)
  pf <- reliability_profile(fm)
  expect_equal(pf$sigma_tau^2, fm$tau2)
})

test_that("an ebnm fit with a normal prior gives the profile at its SD", {
  skip_if_not_installed("ebnm")
  withr::local_seed(173)
  x <- rnorm(30, 0.3, 1.2)
  s <- runif(30, 0.5, 1)
  fit <- ebnm::ebnm(x, s, prior_family = "normal", mode = "estimate")
  prof <- reliability_profile(fit)
  expect_identical(prof$source, "ebnm normal prior")
  expect_equal(prof$sigma_tau, fit$fitted_g$sd)
  # With the mean estimated, ebnm's posterior means are the normal model's.
  expect_lt(max(abs(site_table(prof)$effect - fit$posterior$mean)), 1e-4)
  pn <- ebnm::ebnm(x, s, prior_family = "point_normal")
  expect_rp_error(reliability_profile(pn), "rp_scope_error")
})

test_that("a brms formula is read for its parts", {
  fake <- structure(list(
    formula = list(formula = estimate | se(se) ~ 1 + (1 | school)),
    family = stats::gaussian(),
    data = data.frame(estimate = c(0, 1), se = c(1, 2), school = c("a", "b"))
  ), class = "brmsfit")
  parts <- .rp_brms_parts(fake)
  expect_identical(parts$response, "estimate")
  expect_identical(parts$se, "se")
  expect_identical(parts$group, "school")
  expect_identical(parts$sd_variable, "sd_school__Intercept")
  expect_identical(parts$mean_variable, "b_Intercept")
  other <- structure(list(formula = list(formula = y ~ x + (1 | g))), class = "brmsfit")
  expect_rp_error(.rp_brms_parts(other), "rp_scope_error")
  slope <- structure(list(
    formula = list(formula = y | se(s) ~ 1 + x + (1 | g)),
    family = stats::gaussian(),
    data = data.frame(y = c(0, 1), s = c(1, 2), x = c(-1, 1), g = c("a", "b"))
  ), class = "brmsfit")
  expect_null(.rp_brms_parts(slope)$mean_variable)
})

test_that("fits of metafor other than rma.uni are refused with a scope error", {
  skip_if_not_installed("metafor")
  d <- data.frame(yi = eight_schools$estimate, vi = eight_schools$se^2, study = seq_len(8))
  fit <- metafor::rma.mv(yi, vi, random = ~ 1 | study, data = d)
  expect_rp_error(reliability_posterior(fit), class = "rp_scope_error", regexp = "rma.uni")
  expect_rp_error(reliability_profile(fit), class = "rp_scope_error", regexp = "rma.uni")
})
