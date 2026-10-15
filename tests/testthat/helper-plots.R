# Objects for the plot tests, made once.
plot_objects <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      fit <- reliability_posterior(estimate, se = se, site = school, data = eight_schools)
      prof <- reliability_profile(se = eight_schools$se, sigma_tau = 10, site = eight_schools$school,
                                  estimate = eight_schools$estimate)
      n <- round(stats::qlnorm(stats::ppoints(40), 0, sqrt(log(1.25))) / 40 * 2000)
      design <- reliability_design(n, sigma_tau = c(0.15, 0.2, 0.25), kappa = 4)
      d <- profile_draws(fit, ndraws = 2000, seed = 1)
      draws <- reliability_draws(d$sigma_tau, se = eight_schools$se, site = eight_schools$school,
                                 estimate = eight_schools$estimate, mean = d$effect_mean,
                                 model_label = "exact draws")
      cache <<- list(fit = fit, prof = prof, design = design, draws = draws)
    }
    cache
  }
})

# The data of a layer of a built plot.
layer_data_of <- function(plot, geom) {
  built <- ggplot2::ggplot_build(plot)
  for (i in seq_along(plot$layers)) {
    if (inherits(plot$layers[[i]]$geom, geom)) return(built$data[[i]])
  }
  NULL
}
