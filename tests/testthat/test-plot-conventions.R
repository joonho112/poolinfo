plot_function_names <- function() {
  sort(grep("^plot_", getNamespaceExports("poolinfo"), value = TRUE))
}

test_that("the arguments of every plot function come in the same order", {
  expect_identical(
    plot_function_names(),
    c("plot_compare", "plot_concentration", "plot_contributions", "plot_dashboard", "plot_decomposition", "plot_design",
      "plot_functionals", "plot_gap_path", "plot_geometry", "plot_plane", "plot_pooling", "plot_profile",
      "plot_sigma_scale")
  )
  for (name in plot_function_names()) {
    args <- names(formals(get(name, asNamespace("poolinfo"))))
    expect_identical(args[[1L]], "x", info = name)
    expect_identical(args[[length(args)]], "...", info = name)
    if ("scale" %in% args) expect_identical(which(args == "scale"), length(args) - 1L, info = name)
    if ("max_sites" %in% args) {
      expect_identical(which(args == "max_sites"),
                       if ("scale" %in% args) which(args == "scale") - 1L else length(args) - 1L, info = name)
    }
    if (all(c("type", "level") %in% args)) expect_lt(which(args == "type"), which(args == "level"))
    if (all(c("which", "level") %in% args)) expect_lt(which(args == "which"), which(args == "level"))
  }
})

test_that("every plot is a ggplot without a title and with a caption", {
  o <- plot_objects()
  plots <- list(
    profile = plot_profile(o$fit), profile_fixed = plot_profile(o$prof),
    logit = plot_profile(o$prof, type = "logit"), pooling = plot_pooling(o$fit),
    pooling_fixed = plot_pooling(o$prof), functionals = plot_functionals(o$fit),
    functionals_draws = plot_functionals(o$draws), plane = plot_plane(o$fit),
    concentration = plot_concentration(o$prof), contributions = plot_contributions(o$fit),
    decomposition = plot_decomposition(o$prof), sigma = plot_sigma_scale(o$fit),
    design = plot_design(o$design), design_sites = plot_design(o$design, type = "sites"),
    design_target = plot_design(o$design, type = "target")
  )
  for (name in names(plots)) {
    p <- plots[[name]]
    expect_s3_class(p, "ggplot")
    expect_null(p$labels$title, label = name)
    expect_true(nzchar(p$labels$caption), label = name)
    expect_silent(ggplot2::ggplot_build(p))
  }
  for (p in list(plot_dashboard(o$fit), plot_gap_path(o$prof), plot_geometry(o$prof))) {
    expect_s3_class(p, "patchwork")
    expect_true(nzchar(p$patches$annotation$caption))
  }
})

test_that("captions say where sigma_tau came from", {
  o <- plot_objects()
  expect_match(plot_profile(o$fit)$labels$caption, "conditional on the supplied standard errors")
  expect_match(plot_profile(o$prof)$labels$caption, "sigma_tau = 10 \\(supplied\\)")
  expect_match(plot_functionals(o$draws)$labels$caption, "exact draws")
  expect_match(plot_functionals(o$draws)$labels$caption, "Kernel density")
})

test_that("plot() draws the plot that `type` names", {
  o <- plot_objects()
  same <- function(a, b) expect_identical(ggplot2::ggplot_build(a)$data, ggplot2::ggplot_build(b)$data)
  same(plot(o$prof), plot_profile(o$prof))
  same(plot(o$prof, type = "contributions"), plot_contributions(o$prof))
  same(plot(o$fit, type = "pooling"), plot_pooling(o$fit))
  same(plot(o$fit, type = "sigma"), plot_sigma_scale(o$fit))
  same(plot(o$design), plot_design(o$design))
  same(plot(o$design, type = "profile"), plot_profile(o$design))
  expect_s3_class(plot(o$fit), "patchwork")
  expect_s3_class(plot(o$draws), "patchwork")
  expect_rp_error(plot(o$fit, type = "bars"), regexp = "one of")
})

test_that("plots of many sites leave the labels out and say so", {
  withr::local_seed(161)
  v <- runif(80, 0.01, 0.1)
  fit <- reliability_posterior(rnorm(80, 0, sqrt(0.03 + v)), se2 = v)
  expect_message(plot_profile(fit), class = "rp_message_many_sites")
  expect_message(plot_pooling(fit), class = "rp_message_many_sites")
  expect_message(plot_contributions(fit), class = "rp_message_many_sites")
  expect_no_message(plot_profile(fit, max_sites = 100))
})

test_that("plot functions refuse what they cannot draw", {
  o <- plot_objects()
  expect_rp_error(plot_profile(1), regexp = "takes")
  expect_rp_error(plot_pooling(reliability_profile(se = c(1, 2), sigma_tau = 1)), regexp = "no estimates")
  expect_rp_error(plot_functionals(o$prof), regexp = "posterior fits or draws")
  expect_rp_error(plot_design(o$fit), regexp = "reliability_design")
  expect_rp_error(plot_profile(o$fit, centroids = "median"), regexp = "centroids")
  expect_rp_error(plot_geometry(o$prof, sites = 1:7), regexp = "at most 6")
  expect_rp_error(plot_pooling(o$fit, scale = 0), regexp = "positive")
})
