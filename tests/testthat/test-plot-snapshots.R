# Pictures of the plots, compared with stored ones. A new version of
# ggplot2 can change a picture without changing what it shows, so the
# comparison is made only where POOLINFO_SNAPSHOTS is "true", as in the macOS
# job of the check workflow, and never on CRAN. After a change that is meant,
# review the new pictures with testthat::snapshot_review().

test_that("the plots look as stored", {
  pictures <- c("profile of a posterior", "profile at a fixed sigma_tau", "profile on the logit scale",
                "pooling", "density of I", "plane", "dashboard", "gap path", "concentration",
                "contributions", "geometry", "decomposition", "sigma scale", "design imbalance")
  # Every stored picture is announced before the comparison can be skipped,
  # so that testthat keeps the files when it does not run.
  for (title in pictures) {
    announce_snapshot_file(name = paste0(gsub("^-|-$", "", gsub("[^a-z0-9]+", "-", tolower(title))), ".svg"))
  }
  skip_on_cran()
  skip_if_not_installed("vdiffr")
  skip_if_not(identical(Sys.getenv("POOLINFO_SNAPSHOTS"), "true"), "POOLINFO_SNAPSHOTS is not \"true\"")
  o <- plot_objects()
  plots <- list(
    plot_profile(o$fit), plot_profile(o$prof), plot_profile(o$prof, type = "logit"),
    plot_pooling(o$fit), plot_functionals(o$fit), plot_plane(o$fit), plot_dashboard(o$fit),
    plot_gap_path(o$prof), plot_concentration(o$prof), plot_contributions(o$prof),
    plot_geometry(c(0.2, 0.4, 0.6, 0.85)), plot_decomposition(c(0.05, 0.10, 0.15, 0.25, 0.92)),
    plot_sigma_scale(o$fit), plot_design(o$design)
  )
  for (i in seq_along(pictures)) vdiffr::expect_doppelganger(pictures[[i]], plots[[i]])
})
