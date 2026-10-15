#' Plots of poolinfo objects
#'
#' `plot()` shows one of the package's plots of a profile, a posterior fit,
#' draws or a design. `type` chooses the plot, and most types call the plot
#' function of the same name, so `plot(fit, type = "pooling")` is
#' `plot_pooling(fit)`.
#'
#' @param x A `reliability_profile`, `reliability_posterior`,
#'   `reliability_draws` or `reliability_design` object.
#' @param type The plot:
#'   * for a profile, `"profile"` (the default), `"logit"`, `"geometry"`,
#'     `"decomposition"`, `"gap"`, `"concentration"` or `"contributions"`.
#'   * for a posterior fit or draws, `"dashboard"` (the default),
#'     `"profile"`, `"logit"`, `"pooling"`, `"functionals"`, `"plane"`,
#'     `"sigma"`, `"gap"`, `"concentration"` or `"contributions"`.
#'   * for a design, `"imbalance"` (the default), `"sites"`, `"target"` or
#'     `"profile"`.
#'
#'   `"logit"` is [plot_profile()] with `type = "logit"`, `"sigma"` is
#'   [plot_sigma_scale()], `"gap"` is [plot_gap_path()], and `"imbalance"`,
#'   `"sites"` and `"target"` are the types of [plot_design()]; `"profile"`
#'   of a design shows the planned profile.
#' @param ... Further arguments of the plot function.
#'
#' @return A ggplot or patchwork object.
#'
#' @seealso [plot_profile()], [plot_dashboard()], [plot_design()], and
#'   [compare_profiles()] for the plots of a comparison.
#'
#' @examples
#' fit <- reliability_posterior(estimate, se = se, site = school, data = eight_schools)
#' plot(fit)
#' plot(fit, type = "functionals")
#' @name plot.poolinfo
NULL

.rp_plot_dispatch <- function(x, type, ...) {
  switch(type,
    profile = plot_profile(x, ...),
    logit = plot_profile(x, type = "logit", ...),
    geometry = plot_geometry(x, ...),
    decomposition = plot_decomposition(x, ...),
    gap = plot_gap_path(x, ...),
    concentration = plot_concentration(x, ...),
    contributions = plot_contributions(x, ...),
    dashboard = plot_dashboard(x, ...),
    pooling = plot_pooling(x, ...),
    functionals = plot_functionals(x, ...),
    plane = plot_plane(x, ...),
    sigma = plot_sigma_scale(x, ...),
    imbalance = plot_design(x, type = "imbalance", ...),
    sites = plot_design(x, type = "sites", ...),
    target = plot_design(x, type = "target", ...)
  )
}

#' @rdname plot.poolinfo
#' @export
plot.reliability_profile <- function(x, type = c("profile", "logit", "geometry", "decomposition",
                                                 "gap", "concentration", "contributions"), ...) {
  .rp_plot_dispatch(x, .rp_match(type, c("profile", "logit", "geometry", "decomposition", "gap",
                                         "concentration", "contributions"), "type"), ...)
}

.rp_posterior_plot_types <- c("dashboard", "profile", "logit", "pooling", "functionals", "plane",
                              "sigma", "gap", "concentration", "contributions")

#' @rdname plot.poolinfo
#' @export
plot.reliability_posterior <- function(x, type = c("dashboard", "profile", "logit", "pooling",
                                                   "functionals", "plane", "sigma", "gap",
                                                   "concentration", "contributions"), ...) {
  .rp_plot_dispatch(x, .rp_match(type, .rp_posterior_plot_types, "type"), ...)
}

#' @rdname plot.poolinfo
#' @export
plot.reliability_draws <- function(x, type = c("dashboard", "profile", "logit", "pooling",
                                               "functionals", "plane", "sigma", "gap",
                                               "concentration", "contributions"), ...) {
  .rp_plot_dispatch(x, .rp_match(type, .rp_posterior_plot_types, "type"), ...)
}

#' @rdname plot.poolinfo
#' @export
plot.reliability_design <- function(x, type = c("imbalance", "sites", "target", "profile"), ...) {
  type <- .rp_match(type, c("imbalance", "sites", "target", "profile"), "type")
  if (type == "profile") {
    return(plot_profile(x, ...))
  }
  .rp_plot_dispatch(x, type, ...)
}
