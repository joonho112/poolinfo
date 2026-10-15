test_that("a comparison tabulates each analysis from its accessors", {
  o <- plot_objects()
  cmp <- compare_profiles(posterior = o$fit, draws = o$draws, `plug-in` = reliability_profile(o$fit))
  expect_s3_class(cmp, "profile_comparison")
  t <- as.data.frame(cmp)
  expect_identical(t$analysis, c("posterior", "draws", "plug-in"))
  expect_identical(t$source, c("posterior", "posterior draws", "posterior mean of sigma_tau^2"))
  ps <- profile_summary(o$fit)
  expect_equal(t$I[[1]], ps$mean[ps$quantity == "I"])
  expect_equal(t$H_upper[[1]], ps$upper[ps$quantity == "H"])
  expect_true(is.na(t$I_lower[[3]]))
  expect_equal(t$prob_I_above[[1]], prob_above(o$fit))
  long <- profile_summary(cmp)
  expect_identical(unique(long$analysis), t$analysis)
  expect_snapshot(print(cmp, use_cli = FALSE))
})

test_that("compare_profiles() needs distinct names and analyses of poolinfo", {
  o <- plot_objects()
  expect_rp_error(compare_profiles(o$fit), regexp = "distinct name")
  expect_rp_error(compare_profiles(a = o$fit, a = o$prof), regexp = "distinct name")
  expect_rp_error(compare_profiles(a = 1), regexp = "not an analysis")
  expect_rp_error(compare_profiles(), regexp = "at least one")
  expect_identical(names(compare_profiles(.list = list(x = o$fit))$objects), "x")
})

test_that("prior_sensitivity() refits under each prior", {
  o <- plot_objects()
  cmp <- prior_sensitivity(o$fit, list(hn = sd_prior("half_normal", scale = 10)))
  expect_identical(names(cmp$objects), c("flat", "hn"))
  expect_identical(cmp$objects$hn$prior$type, "half_normal")
  expect_equal(cmp$objects$hn$summary,
               update(o$fit, prior = sd_prior("half_normal", scale = 10))$summary)
  expect_rp_error(prior_sensitivity(o$fit, sd_prior("half_normal", scale = 1)), regexp = "named list")
  expect_rp_error(prior_sensitivity(o$prof, list(a = sd_prior())), regexp = "reliability_posterior")
})

test_that("report_table() follows the paper's table and carries the note", {
  o <- plot_objects()
  t <- report_table(o$fit)
  expect_named(t, c("analysis", "J", "I", "I_lower", "I_upper", "mean", "H", "H_lower", "H_upper",
                    "J_eff", "prob_I_above"))
  expect_match(attr(t, "note"), "conditional on the supplied standard errors")
  md <- report_table(o$fit, format = "markdown")
  ps <- profile_summary(o$fit)
  I <- ps$mean[ps$quantity == "I"]
  expect_match(md, sub("^0", "", formatC(I, format = "f", digits = 2)), fixed = TRUE)
  expect_match(md, "*Note.*", fixed = TRUE)
  rows <- strsplit(md, "\n", fixed = TRUE)[[1]][3]
  expect_false(grepl(" 0\\.", rows))
  # One analysis has no column for its name.
  expect_false(grepl("Analysis", md, fixed = TRUE))
  tex <- report_table(o$fit, format = "latex")
  expect_match(tex, "\\toprule", fixed = TRUE)
  # Every percent sign is escaped for LaTeX.
  expect_false(grepl("[^\\\\]%", tex))
  expect_match(tex, "95\\%", fixed = TRUE)
  expect_match(tex, "\\textit{Note.}", fixed = TRUE)
  expect_match(tex, "conditional on the supplied standard errors")
})

test_that("report_text() always ends with the conditioning and never judges", {
  o <- plot_objects()
  for (x in list(o$fit, o$draws, o$prof, o$design)) {
    s <- report_text(x)
    expect_type(s, "character")
    expect_match(s[[length(s)]], "conditional on the supplied standard errors|computed at the planning values")
    expect_false(any(grepl("\\b(low|high|adequate|poor|good|weak|strong)\\b", s)))
  }
  s <- report_text(o$fit, outcome = "the SAT-V")
  ps <- profile_summary(o$fit)
  expect_match(s[[1]], .rp_fmt(ps$mean[ps$quantity == "I"], 2L), fixed = TRUE)
  expect_match(s[[1]], "flat priors on the mean and the between-site standard deviation")
  hn <- update(o$fit, prior = sd_prior("half_normal", scale = 10))
  expect_match(report_text(hn)[[1]], "half-normal prior with scale 10")
  labelled <- reliability_posterior(estimate, se = se, data = eight_schools, scale_label = "SAT-V point",
                                    target_label = "the effect of coaching on the SAT-V")
  s <- report_text(labelled)
  expect_true(any(grepl("The effects are on the SAT-V point scale.", s, fixed = TRUE)))
  expect_true(any(grepl("The estimand at each site is the effect of coaching on the SAT-V.", s, fixed = TRUE)))
  cmp <- compare_profiles(a = o$fit, b = o$prof)
  expect_identical(report_text(cmp)[[1]], "a:")
})

test_that("the example profiles are rebuilt by shape_profile()", {
  expect_identical(levels(example_profiles$profile),
                   c("homogeneous", "even", "bimodal", "concentrated_precise", "concentrated_noisy"))
  for (name in levels(example_profiles$profile)) {
    d <- example_profiles[example_profiles$profile == name, ]
    expect_equal(kl_centroid(d$reliability), 0.3, tolerance = 1e-12)
    expect_equal(logit_moments(d$reliability)[["V"]], if (name == "homogeneous") 0 else 1.5, tolerance = 1e-12)
    shape <- if (name == "homogeneous") "even" else name
    expect_equal(d$se2, shape_profile(40, 0.3, if (name == "homogeneous") 0 else 1.5, shape, sigma_tau = 0.2),
                 tolerance = 1e-14)
  }
  # The values of the paper's Figure 2.
  even <- example_profiles[example_profiles$profile == "even", "reliability"]
  conc <- example_profiles[example_profiles$profile == "concentrated_precise", "reliability"]
  expect_equal(round(kl_dispersion(even), 3), 0.149)
  expect_equal(round(kl_dispersion(conc), 3), 0.168)
  expect_equal(round(concentration_count(conc), 1), 24.1)
})

test_that("eight_schools has the eight estimates of Rubin (1981)", {
  expect_identical(eight_schools$estimate, c(28, 8, -3, 7, -1, 1, 18, 12))
  expect_identical(eight_schools$se, c(15, 10, 16, 11, 9, 11, 10, 18))
})

test_that("ocrs has the 32 sites of the analytic sample", {
  expect_identical(dim(ocrs), c(32L, 18L))
  expect_identical(ocrs$site, sprintf("S%02d", 1:32))
  expect_identical(as.vector(table(ocrs$subject)), c(13L, 19L))
  expect_identical(length(unique(ocrs$school)), 24L)
  expect_identical(sum(table(ocrs$school) == 2L), 8L)
  expect_identical(ocrs$n, ocrs$n_treatment + ocrs$n_control)
  expect_identical(c(range(ocrs$n), stats::median(ocrs$n)), c(24, 104, 49))
  expect_identical(sum(ocrs$n), 1737L)
  expect_identical(sum(ocrs$blocks), 144L)
  expect_gte(min(ocrs$n_treatment, ocrs$n_control), 8L)
  se <- ocrs[grep("^se", names(ocrs))]
  expect_true(all(vapply(se, function(s) all(s > 0), NA)))
  expect_false(any(ocrs$zero_cell))
})

test_that("the comparison of the OCRS analyses gives the paper's Table 5", {
  g <- ocrs_analyses()
  fits <- lapply(g, function(z) reliability_posterior(z$y, se2 = z$v))
  cmp <- compare_profiles(.list = fits)
  t <- as.data.frame(cmp)
  for (id in names(g)) {
    row <- t[t$analysis == id, ]
    ref <- g[[id]]$summary
    expect_lt(abs(row$I - ref[["I.mean"]]), tol_bounded)
    expect_lt(abs(row$I_lower - ref[["I.q025"]]), tol_bounded)
    expect_lt(abs(row$H - ref[["H.mean"]]), tol_bounded)
    expect_lt(abs(row$prob_I_above - ref[["P_I_gt_half"]]), tol_bounded)
  }
  md <- report_table(cmp, format = "markdown")
  expect_match(md, "| pooled | 32 | .67 [.50, .81] | .66 | .025 [.019, .029] | 31.2 | .97 |", fixed = TRUE)
  expect_match(md, "| credits | 32 | .32 [.02, .60] |", fixed = TRUE)
})

test_that("printed output fits in 80 columns, with long labels and names", {
  withr::local_options(poolinfo.use_cli = FALSE)
  target <- "the block-adjusted intention-to-treat contrast for course passing in each school"
  fit <- reliability_posterior(estimate, se = se, site = site, data = ocrs, mean_model = ~ subject,
                               scale_label = "risk difference", target_label = target)
  prof <- reliability_profile(fit)
  withr::local_seed(91)
  draws <- reliability_draws(sqrt(stats::rgamma(400, 20, 20 / 0.035)), se = ocrs$se,
                             model_label = "a brms fit with a half-normal prior on sd(Intercept)")
  shown <- list(
    fit, summary(fit), prof, summary(prof), draws,
    compare_profiles(`pooled, adjusted for covariates` = fit, `plug-in at E(sigma_tau^2)` = prof),
    reliability_design(n = c(rep(20, 500), rep(2000, 500)), sigma_tau = 0.25,
                       kappa = c(rep(0.25, 500), rep(4, 500)))
  )
  for (x in shown) {
    out <- utils::capture.output(print(x))
    expect_lte(max(nchar(out)), 80L, label = class(x)[[1L]])
  }
  expect_identical(.rp_wrap("a b c d", width = 3L), c("a", "b c", "d"))
})

test_that("a design with costs is reported against equal spending", {
  d <- reliability_design(rep(50, 40), sigma_tau = 0.15, kappa = 4, cost = rep(c(1, 2), each = 20))
  text <- report_text(d, digits = 3)
  expect_match(text[[2]], "with equal spending", fixed = TRUE)
  expect_false(any(grepl("balanced", text, fixed = TRUE)))
})

test_that("a comparison of profiles plots without a type", {
  p1 <- reliability_profile(se = eight_schools$se, sigma_tau = 5)
  p2 <- reliability_profile(se = eight_schools$se, sigma_tau = 10)
  expect_s3_class(plot(compare_profiles(a = p1, b = p2)), "ggplot")
})

test_that("draws carry the labels of the scale and the estimand", {
  withr::local_options(poolinfo.use_cli = FALSE)
  withr::local_seed(93)
  rd <- reliability_draws(abs(stats::rnorm(500, 8, 4)), se = eight_schools$se,
                          scale_label = "SAT-V point", target_label = "the effect of coaching")
  out <- format(rd)
  expect_true(any(out == "Scale of the effects: SAT-V point"))
  expect_true(any(out == "Site estimates: the effect of coaching"))
  expect_match(report_text(rd)[[1]], "SAT-V point scale", fixed = TRUE)
  expect_identical(reliability_profile(rd)$labels$scale_label, "SAT-V point")
})
