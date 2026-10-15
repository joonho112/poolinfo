# Reference posterior summaries computed for Lee (2026), recomputed by the
# package: the six OCRS analyses, with inputs from the data set ocrs, and six
# simulated data sets, in cases/.

# Compares a fit with a named vector of stored reference summaries, with the tolerances of the tolerance table.
expect_matches_reference <- function(fit, ref, v, label) {
  s <- fit$summary
  q <- fit$quantiles
  get <- function(quantity, col) s[[col]][s$quantity == quantity]
  getq <- function(quantity, col) q[[col]][q$quantity == quantity]
  scale_tol <- tol_effect_scale(v)
  bounded <- list(I = "I", Sbar = "mean", H = "H")
  for (name in names(bounded)) {
    quantity <- bounded[[name]]
    expect_lt(abs(get(quantity, "mean") - ref[[paste0(name, ".mean")]]), tol_bounded, label = paste(label, name))
    if (paste0(name, ".q025") %in% names(ref)) {
      expect_lt(abs(getq(quantity, "lower") - ref[[paste0(name, ".q025")]]), tol_bounded)
      expect_lt(abs(getq(quantity, "median") - ref[[paste0(name, ".q500")]]), tol_bounded)
      expect_lt(abs(getq(quantity, "upper") - ref[[paste0(name, ".q975")]]), tol_bounded)
    }
  }
  expect_lt(rel_diff(get("J_eff", "mean"), ref[["Jeff.mean"]]), tol_jeff_rel)
  if ("Jeff.q025" %in% names(ref)) {
    expect_lt(rel_diff(getq("J_eff", "lower"), ref[["Jeff.q025"]]), tol_jeff_rel)
    expect_lt(rel_diff(getq("J_eff", "upper"), ref[["Jeff.q975"]]), tol_jeff_rel)
  }
  expect_lt(abs(get("sigma_tau", "mean") - ref[["sigma.mean"]]), scale_tol)
  expect_lt(abs(get("effect_mean", "mean") - ref[["tau.mean"]]), scale_tol)
  expect_lt(rel_diff(get("sigma_tau2", "mean"), ref[["a.mean"]]), tol_sigma2_rel)
  expect_lt(abs(fit$plugin[["I"]] - ref[["I_plugin"]]), tol_bounded)
  expect_lt(abs(fit$probability[["I_above"]] - ref[["P_I_gt_half"]]), tol_bounded)
  expect_lt(abs(get("V", "mean") - ref[["V"]]), 1e-12)
  if ("sigma.q025" %in% names(ref)) {
    expect_lt(abs(getq("sigma_tau", "lower") - ref[["sigma.q025"]]), scale_tol)
    expect_lt(abs(getq("sigma_tau", "upper") - ref[["sigma.q975"]]), scale_tol)
    expect_lt(rel_diff(getq("sigma_tau2", "median"), ref[["a.q500"]]), tol_sigma2_rel)
  }
  expect_true(all(fit$checks$status == "ok"), label = paste(label, "numerical checks"))
}

test_that("both engines reproduce the reference summaries of the OCRS analyses", {
  analyses <- ocrs_analyses()
  for (id in names(analyses)) {
    g <- analyses[[id]]
    for (engine in c("nodes", "adaptive")) {
      fit <- reliability_posterior(g$y, se2 = g$v, control = posterior_control(engine = engine))
      expect_matches_reference(fit, g$summary, g$v, paste(id, engine))
      expect_equal(fit$h_branches$turning_logits, g$H_critical_logits, tolerance = tol_turning)
    }
  }
})

test_that("both engines reproduce the reference summaries of six simulated data sets", {
  cases <- read_boundary_cases()
  for (i in seq_along(cases)) {
    cc <- cases[[i]]
    for (engine in c("nodes", "adaptive")) {
      fit <- reliability_posterior(cc$y, se2 = cc$v, control = posterior_control(engine = engine))
      expect_matches_reference(fit, cc$reference, cc$v, sprintf("case %d %s", i, engine))
    }
  }
})

test_that("the profile at the posterior mean of sigma_tau^2 is the plug-in profile", {
  g <- ocrs_analyses()$pooled
  fit <- reliability_posterior(g$y, se2 = g$v)
  prof <- reliability_profile(fit)
  expect_identical(prof$source, "posterior mean of sigma_tau^2")
  expect_equal(prof$sigma_tau^2, g$summary[["a.mean"]], tolerance = tol_sigma2_rel)
  expect_equal(prof$stats[["I"]], g$summary[["I_plugin"]], tolerance = tol_bounded)
})
