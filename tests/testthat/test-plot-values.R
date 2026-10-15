# The numbers that the plots draw are those of the accessors.

test_that("the profile plot draws the reliabilities and the centroids", {
  o <- plot_objects()
  p <- plot_profile(o$fit)
  pts <- layer_data_of(p, "GeomPoint")
  tab <- site_table(o$fit)
  expect_equal(sort(pts$y), sort(tab$reliability))
  lines <- layer_data_of(p, "GeomHline")
  ps <- profile_summary(o$fit)
  expect_equal(sort(lines$yintercept), sort(ps$mean[ps$quantity %in% c("I", "mean")]))
  ranges <- layer_data_of(p, "GeomLinerange")
  expect_equal(sort(ranges$ymin), sort(tab$reliability_lower))
  expect_equal(sort(ranges$ymax), sort(tab$reliability_upper))
  band <- layer_data_of(p, "GeomRect")
  expect_equal(c(band$ymin, band$ymax), unlist(ps[ps$quantity == "I", c("lower", "upper")]), ignore_attr = TRUE)
})

test_that("the pooling plot draws the site table's posterior means and intervals", {
  o <- plot_objects()
  p <- plot_pooling(o$fit, scale = 2)
  tab <- site_table(o$fit)
  built <- ggplot2::ggplot_build(p)
  points <- built$data[vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), TRUE)]
  expect_equal(sort(points[[1]]$x), sort(2 * tab$estimate))
  expect_equal(sort(points[[2]]$x), sort(2 * tab$effect))
  ranges <- built$data[vapply(p$layers, function(l) inherits(l$geom, "GeomLinerange"), TRUE)]
  expect_equal(sort(ranges[[2]]$xmin), sort(2 * tab$effect_lower))
  vline <- layer_data_of(p, "GeomVline")
  ps <- profile_summary(o$fit)
  expect_equal(vline$xintercept, 2 * ps$mean[ps$quantity == "effect_mean"])
})

test_that("the density plot draws functional_density() and the summaries", {
  o <- plot_objects()
  p <- plot_functionals(o$fit, which = "H")
  line <- layer_data_of(p, "GeomLine")
  d <- functional_density(o$fit, "H", n = 512)
  expect_equal(line$x, d$value)
  expect_equal(line$y, d$density)
  segs <- layer_data_of(p, "GeomSegment")
  ps <- profile_summary(o$fit)
  expect_equal(segs$x[[1]], ps$mean[ps$quantity == "H"])
})

test_that("the plane plot places each analysis at its I and V", {
  o <- plot_objects()
  pts <- layer_data_of(plot_plane(o$fit), "GeomPoint")
  ps <- profile_summary(o$fit)
  expect_equal(pts$x, ps$mean[ps$quantity == "I"])
  expect_equal(pts$y, ps$mean[ps$quantity == "V"])
})

test_that("the gap path plot draws profile_path()", {
  o <- plot_objects()
  p <- plot_gap_path(o$prof)
  path <- profile_path(o$prof, n = 399)
  built <- ggplot2::ggplot_build(p[[1]])
  lines <- built$data[vapply(p[[1]]$layers, function(l) inherits(l$geom, "GeomLine"), TRUE)]
  expect_equal(lines[[2]]$y, path$H)
  expect_equal(lines[[1]]$y, path$H_leading)
})

test_that("the concentration plot's steps add up the shares", {
  o <- plot_objects()
  step <- layer_data_of(plot_concentration(o$prof, label_top = 0), "GeomStep")
  q <- sort(site_table(o$prof)$share, decreasing = TRUE)
  expect_equal(step$y, c(0, cumsum(q)))
  # With the 32 OCRS sites the proportions add up to a little more than 1;
  # the curve still ends at 1.
  fit <- reliability_posterior(estimate, se = se, data = ocrs)
  expect_no_warning(step <- layer_data_of(plot_concentration(fit, label_top = 0), "GeomStep"))
  expect_identical(step$y[[nrow(step)]], 1)
})

test_that("the contributions plot draws the sites' terms and H", {
  o <- plot_objects()
  p <- plot_contributions(o$prof)
  pts <- layer_data_of(p, "GeomPoint")
  expect_equal(sort(pts$y), sort(site_table(o$prof)$kl_term))
  expect_equal(layer_data_of(p, "GeomHline")$yintercept, o$prof$stats[["H"]])
})

test_that("the sigma plot shades I above the threshold", {
  o <- plot_objects()
  area <- layer_data_of(plot_sigma_scale(o$fit), "GeomArea")
  sigma_c <- sqrt(exp(mean(log(eight_schools$se^2))))
  expect_true(all(area$x >= sigma_c - 1e-12))
})

test_that("scale multiplies the effects and sigma_tau and leaves the reliabilities", {
  o <- plot_objects()
  a <- layer_data_of(plot_sigma_scale(o$fit), "GeomLine")
  b <- layer_data_of(plot_sigma_scale(o$fit, scale = 100), "GeomLine")
  expect_equal(b$x, 100 * a$x)
  expect_equal(b$y, a$y / 100)
})

test_that("the pooling plot of a fit with a mean model draws no single mean", {
  d <- data.frame(y = eight_schools$estimate, se = eight_schools$se, x = c(1, 0, 1, 0, 1, 0, 1, 1))
  p <- plot_pooling(reliability_posterior(y, se = se, data = d, mean_model = ~x))
  expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomVline"), NA)))
  expect_match(p$labels$caption, "predicted value", fixed = TRUE)
})
