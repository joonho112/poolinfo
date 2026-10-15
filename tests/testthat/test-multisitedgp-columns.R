test_that("a missing simulation design identifies the absent SD and a usable alternative", {
  x <- structure(data.frame(site_index = 1:4, tau_j_hat = c(-0.2, 0.1, 0.3, 0.5),
                            se2_j = c(0.1, 0.2, 0.15, 0.1)), class = c("multisitedgp_data", "data.frame"))
  for (design in list(NULL, list(J = 4))) {
    attr(x, "design") <- design
    err <- expect_rp_error(reliability_profile(x), regexp = "absent or has no `sigma_tau`")
    expect_match(conditionMessage(err), "se2 = x\\$se2_j", fixed = FALSE)
    expect_false(grepl("Selecting rows", conditionMessage(err), fixed = TRUE))
    expect_false(grepl("sigma_tau = )", conditionMessage(err), fixed = TRUE))
  }
  direct <- reliability_profile(se2 = x$se2_j, sigma_tau = 0.2)
  expect_s3_class(direct, "reliability_profile")
})

test_that("row selection of simulated data preserves its known design SD", {
  skip_if_not_installed("multisiteDGP")
  design <- multisiteDGP::multisitedgp_design(J = 12, sigma_tau = 0.2, seed = 7)
  x <- multisiteDGP::sim_multisite(design)
  subset <- suppressWarnings(x[1:8, ])
  expect_s3_class(subset, "multisitedgp_data")
  expect_equal(attr(subset, "design")$sigma_tau, 0.2)
  profile <- reliability_profile(subset)
  direct <- reliability_profile(se2 = subset$se2_j, sigma_tau = 0.2)
  expect_equal(profile$stats, direct$stats)
  expect_equal(reliability_posterior(subset)$summary,
               reliability_posterior(subset$tau_j_hat, se2 = subset$se2_j)$summary)
})
