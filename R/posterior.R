#' Settings of the posterior computation
#'
#' Chooses the engine that integrates over the between-site standard
#' deviation and sets its tolerances.
#'
#' @param engine The engine: `"nodes"` (the default) builds one quadrature
#'   rule and reuses it for every summary of a fit; `"adaptive"` computes
#'   each posterior mean, probability and quantile with [stats::integrate()],
#'   as Lee (2026, Appendix D.3) describes.
#' @param rel_tol,abs_tol Relative and absolute tolerances of
#'   [stats::integrate()] for the adaptive engine, positive numbers. The
#'   absolute tolerance applies to integrands scaled to the size of the
#'   posterior and of the integrated quantity, so it does not depend on the
#'   units of the estimates.
#' @param chunk_sites Number of sites computed at a time for the site
#'   summaries, a positive integer. Smaller values use less memory.
#'
#' @details
#' Both engines integrate the same posterior, and their summaries agree to
#' many more digits than are printed. The node rule is the faster: one rule
#' of Gauss-Legendre panels serves every summary of a fit, while the
#' adaptive engine integrates each summary, and each site, separately. The
#' adaptive engine takes the posterior standard deviations, and the
#' intervals of the site effects and of the coefficients, from a node rule.
#' `vignette("implementation")` describes the rule, its checks and the
#' comparison of the engines.
#' With a prior whose tail is very heavy, such as a half-t prior with
#' far fewer than one degree of freedom, a moment can exist and still have
#' a tail too slow for the rule to bound. Both engines then stop with an
#' error that names the moment. A moment that does not exist is `NA`
#' instead. `vignette("implementation")` gives the limits.
#'
#' @return A list of class `posterior_control` with the components
#'   * `engine`: `"nodes"` or `"adaptive"`.
#'   * `rel_tol`, `abs_tol`: the tolerances of the adaptive engine.
#'   * `chunk_sites`: the number of sites computed at a time.
#'
#' @seealso [reliability_posterior()], whose argument `control` takes these
#'   settings, and [check_numerics()] for the checks of a fit.
#'
#' @examples
#' posterior_control()
#' posterior_control(engine = "adaptive")
#' @export
posterior_control <- function(engine = c("nodes", "adaptive"), rel_tol = 1e-10, abs_tol = 1e-12,
                              chunk_sites = 2000) {
  engine <- .rp_match(engine, c("nodes", "adaptive"), "engine")
  structure(
    list(
      engine = engine,
      rel_tol = .rp_check_number(rel_tol, "rel_tol", "a single positive number", 0, lower_open = TRUE),
      abs_tol = .rp_check_number(abs_tol, "abs_tol", "a single positive number", 0, lower_open = TRUE),
      chunk_sites = .rp_check_whole(chunk_sites, "chunk_sites")
    ),
    class = "posterior_control"
  )
}

#' @export
print.posterior_control <- function(x, ...) {
  cat(sprintf("Posterior computation: %s engine", x$engine),
      if (x$engine == "adaptive") sprintf(", rel_tol %s, abs_tol %s", x$rel_tol, x$abs_tol),
      "\n", sep = "")
  invisible(x)
}

#' Posterior summaries of the reliability profile
#'
#' Computes the posterior distribution of the informativeness, the mean
#' reliability, the dispersion and the concentration count of a set of site
#' estimates with standard errors, under the normal hierarchical model, and
#' the posterior means and intervals of the site effects.
#'
#' @param estimate The site estimates: a numeric vector, or a column of
#'   `data`. It can also be a random-effects fit of metafor (`rma.uni`) or a
#'   data set simulated by multisiteDGP, which carry the estimates and their
#'   sampling variances (see [metafor-methods] and [multisitedgp-methods]).
#' @param se,se2 The standard errors of the estimates, or their sampling
#'   variances (the squared standard errors): positive numbers, or a column
#'   of `data`. Give one of the two.
#' @param site Optional labels of the sites, all distinct, used in the
#'   tables and plots.
#' @param data A data frame in which `estimate`, `se`, `se2`, `site` and the
#'   variables of `mean_model` are found. A column takes precedence over a
#'   variable of the same name in the calling environment; write `.env$x`
#'   for the variable.
#' @param mean_model A one-sided formula for the mean of the effects, as in a
#'   meta-regression; `~ 1`, the default, is a common mean.
#' @param prior The prior of the between-site standard deviation, made by
#'   [sd_prior()]; flat by default. In `update()`, `NULL` keeps the prior of
#'   the fit, as it does for `level`, `threshold` and `control`.
#' @param level Probability content of the posterior intervals, a number
#'   between 0 and 1.
#' @param threshold The value of the informativeness whose posterior
#'   probability of being exceeded is reported, between 0 and 1. See
#'   [prob_above()].
#' @param scale_label,target_label Optional strings that name the scale of
#'   the effects, such as `"risk difference"` (read as "on the risk difference
#'   scale"), and the estimand at each site, such as `"the block-adjusted
#'   intention-to-treat contrast for course passing"`. They appear in the
#'   printed output and in [report_text()].
#' @param na_action `"fail"` (the default) stops at a site with a missing
#'   value; `"omit"` leaves such sites out, with a message that names them.
#' @param control Settings of the numerical integration, from
#'   [posterior_control()].
#' @param ... Must be empty. The arguments after `estimate` are named,
#'   because an unnamed vector could hold standard errors or sampling
#'   variances.
#'
#' @details
#' The working model is
#' \deqn{\hat\tau_j \mid \tau_j \sim N(\tau_j, v_j), \qquad \tau_j \sim N(x_j^\top\beta, \sigma_\tau^2),}{tau_hat_j | tau_j ~ N(tau_j, v_j), tau_j ~ N(x_j' beta, sigma_tau^2),}
#' independently over sites, with the sampling variances \eqn{v_j} taken as
#' known. The coefficients \eqn{\beta} (by default the common mean) have a
#' flat prior and are integrated out, which leaves a posterior for
#' `sigma_tau` alone (Lee, 2026, Appendix D). The default prior of
#' `sigma_tau` is flat too, as in Lee (2026). Each summary is a
#' one-dimensional integral over this posterior, computed numerically to
#' about 10 digits, so there is no Monte Carlo error.
#'
#' The diagnostics are functions of `sigma_tau` and the \eqn{v_j}. The
#' informativeness \eqn{I}, the mean reliability, the harmonic centroid and
#' the concentration count increase with `sigma_tau`. The posterior
#' intervals are central, with probability \eqn{(1 - \mathrm{level})/2}{(1 - level) / 2}
#' below and above them, so for these four quantities they are the images
#' of the interval of `sigma_tau`. The dispersion \eqn{H} first rises and
#' then falls, and its quantiles take both branches into account. The logit
#' variance \eqn{V} depends on the \eqn{v_j} only and has no posterior
#' distribution.
#'
#' The summaries are conditional on the supplied standard errors, which are
#' treated as exact, and every site counts the same in them, whatever its
#' size. The last two lines of the printed header state both.
#'
#' Under the flat prior the posterior is proper with at least \eqn{p + 2}
#' sites, where \eqn{p} is the number of columns of the mean model (3 sites
#' for a common mean). The posterior mean of `sigma_tau` needs at least
#' \eqn{p + 3} sites and that of `sigma_tau^2` at least \eqn{p + 4}. A
#' posterior standard deviation needs the next moment as well: at least
#' \eqn{p + 4} sites for `sigma_tau` and the mean, and \eqn{p + 6} for
#' `sigma_tau^2`. The posterior mean of the common mean or a mean-model
#' coefficient requires \eqn{E(\sigma_\tau)}{E(sigma_tau)} to exist, hence
#' at least \eqn{p + 3} sites under the flat prior. With fewer sites its
#' mean is not defined, although a central posterior interval is available.
#' Undefined moments are `NA`, and the printed
#' fit marks them as not defined. [sd_prior()] gives the conditions under
#' the other priors. The paper's
#' results are for a common mean and the flat prior. Mean models with
#' covariates and proper SD priors are experimental extensions beyond Lee
#' (2026); their interfaces may change in a later 0.x release.
#' `vignette("implementation")` describes their computation.
#'
#' @return `reliability_posterior()` and `update()` return an object of class
#'   `reliability_posterior`, a list whose main components are
#'   * `summary`: a data frame with the posterior mean, standard deviation
#'     and interval of each diagnostic, of `sigma_tau`, `sigma_tau^2` and of
#'     the mean or the coefficients.
#'   * `quantiles`: the lower end, median and upper end of the same
#'     quantities.
#'   * `probability`: the posterior probability that \eqn{I} exceeds
#'     `threshold`.
#'   * `plugin`: the diagnostics at the posterior mean of `sigma_tau^2`.
#'   * `sites`: the posterior means of the site reliabilities and site
#'     effects.
#'   * `checks`: the checks of the numerical integration.
#'   * `data`, `prior`, `control`, `level`, `threshold`, `labels` and `call`:
#'     the inputs.
#'
#'   The accessors [profile_summary()], [site_table()], [prob_above()],
#'   [plugin_gap()], [functional_density()], [profile_draws()] and
#'   [check_numerics()] return these numbers as data frames and vectors, and
#'   `plot()` shows [plot_dashboard()]. `summary()` returns an object that
#'   prints the longer summary, and `as.data.frame()` returns the site table.
#'
#' @seealso [reliability_profile()] for a profile at a fixed `sigma_tau`;
#'   [reliability_draws()] for posterior draws from other software;
#'   [compare_profiles()] and [prior_sensitivity()] for several analyses;
#'   `vignette("poolinfo")` for a worked analysis and
#'   `vignette("posterior-summaries")` for the choices behind the summaries.
#'
#' @template ref-lee2026
#' @order 1
#' @examples
#' fit <- reliability_posterior(estimate, se = se, site = school, data = eight_schools)
#' fit
#' profile_summary(fit)
#' update(fit, prior = sd_prior("half_normal", scale = 10))
#' @export
reliability_posterior <- function(estimate, ..., se = NULL, se2 = NULL, site = NULL, data = NULL,
                                  mean_model = ~1, prior = sd_prior("flat"), level = 0.95,
                                  threshold = 0.5, scale_label = NULL, target_label = NULL,
                                  na_action = c("fail", "omit"), control = posterior_control()) {
  .rp_check_named_variances(...)
  if (missing(estimate)) {
    .rp_abort("`estimate` is missing.", info = "Give the estimates of the sites.")
  }
  # A fit of metafor or data simulated by multisiteDGP in the place of the
  # estimates: their estimates, sampling variances and (for metafor)
  # moderators. `estimate` is looked at only without `data`, where looking
  # at it cannot break the lookup of a column name.
  if (is.null(data) && inherits(estimate, "rma") && !inherits(estimate, "rma.uni")) .rp_abort_rma_scope(estimate)
  from_object <- is.null(data) && inherits(estimate, c("rma.uni", "multisitedgp_data"))
  if (from_object) {
    source <- estimate
    if (!is.null(se) || !is.null(se2)) {
      .rp_abort(sprintf("With %s, leave out `se` and `se2`.", .rp_describe(source)),
                info = "They are read from the object.")
    }
    cols <- if (inherits(source, "rma.uni")) .rp_rma_columns(source) else .rp_multisitedgp_columns(source)
    if (inherits(source, "rma.uni") && !is.null(cols$X) && identical(mean_model, ~1)) {
      X <- cols$X
      mean_model <- stats::as.formula(paste("~ 0 +", paste(sprintf("`%s`", colnames(X)), collapse = " + ")))
      data <- as.data.frame(X, check.names = FALSE)
    }
    estimate <- cols$estimate
    se2 <- cols$se2
    if (is.null(site)) site <- cols$site
  }
  if (!is.null(data) && !from_object) {
    if (!is.data.frame(data)) {
      .rp_abort("`data` must be a data frame.", info = .rp_you_passed(data))
    }
    cols <- .rp_eval_columns(data, list(
      estimate = rlang::enquo(estimate), se = rlang::enquo(se), se2 = rlang::enquo(se2),
      site = rlang::enquo(site)
    ))
    estimate <- cols$estimate
    se <- cols$se
    se2 <- cols$se2
    site <- cols$site
  }
  inputs <- .rp_posterior_inputs(estimate, se, se2, site, data, mean_model,
                                 .rp_match(na_action, c("fail", "omit"), "na_action"))
  .rp_new_posterior(
    inputs = inputs, prior = .rp_check_prior(prior), level = .rp_check_level(level),
    threshold = .rp_check_level(threshold, "threshold"),
    labels = list(scale_label = .rp_check_string(scale_label, "scale_label"),
                  target_label = .rp_check_string(target_label, "target_label")),
    control = .rp_check_control(control), call = match.call()
  )
}

.rp_check_control <- function(control) {
  if (!inherits(control, "posterior_control")) {
    .rp_abort("`control` must be made by `posterior_control()`.", info = .rp_you_passed(control))
  }
  control
}

# The estimates, variances, labels and design matrix, with missing values
# handled and the mean model checked.
.rp_posterior_inputs <- function(estimate, se, se2, site, data, mean_model, na_action) {
  for (arg in c("se", "se2")) {
    x <- if (arg == "se") se else se2
    if (is.matrix(x) && nrow(x) == ncol(x) && nrow(x) > 1L) {
      .rp_abort(
        "Correlated estimates are outside the model of this package.",
        info = sprintf("`%s` is a %d x %d matrix.", arg, nrow(x), ncol(x)),
        fix = "Give one standard error or sampling variance per site, for independent estimates.",
        class = "rp_scope_error"
      )
    }
  }
  if (!is.numeric(estimate) || (length(dim(estimate)) > 1L && min(dim(estimate)) > 1L)) {
    .rp_abort("`estimate` must be a numeric vector.", info = .rp_you_passed(estimate))
  }
  J <- length(estimate)
  given <- if (!is.null(se)) se else se2
  if (!is.null(se) && !is.null(se2)) {
    .rp_variances(se, se2)
  }
  if (is.null(given)) {
    .rp_variances(NULL, NULL)
  }
  if (length(given) != J) {
    .rp_abort(
      sprintf("`%s` must have one value for each of the %d estimates.", if (!is.null(se)) "se" else "se2", J),
      info = .rp_you_passed(given)
    )
  }
  if (!is.null(site) && length(site) != J) {
    .rp_abort(sprintf("`site` must hold one label for each of the %d estimates.", J),
              info = .rp_you_passed(site))
  }
  X <- .rp_model_matrix(mean_model, data, J)
  missing_row <- is.na(estimate) | is.na(given) | (if (!is.null(site)) is.na(site) else FALSE)
  if (!is.null(X)) missing_row <- missing_row | rowSums(is.na(X)) > 0
  labels <- .rp_site_labels(NULL, J, names(estimate))
  if (!is.null(site)) labels[!is.na(site)] <- as.character(site[!is.na(site)])
  if (any(missing_row)) {
    if (na_action == "fail") {
      .rp_abort(
        sprintf("%d site%s %s a missing value.", sum(missing_row), if (sum(missing_row) > 1L) "s" else "",
                if (sum(missing_row) > 1L) "have" else "has"),
        info = sprintf("Position%s %s.", if (sum(missing_row) > 1L) "s" else "",
                       .rp_positions(which(missing_row))),
        fix = "Use `na_action = \"omit\"` to leave those sites out."
      )
    }
    .rp_inform(
      sprintf("Left out %d site%s with a missing value: %s.", sum(missing_row),
              if (sum(missing_row) > 1L) "s" else "",
              paste(utils::head(labels[missing_row], 10L), collapse = ", ")),
      subclass = "rp_message_omitted"
    )
    keep <- !missing_row
    estimate <- estimate[keep]
    if (!is.null(se)) se <- se[keep] else se2 <- se2[keep]
    if (!is.null(site)) site <- site[keep]
    if (!is.null(X)) X <- X[keep, , drop = FALSE]
    labels <- labels[keep]
  }
  estimate <- .rp_check_vector(estimate, "estimate")
  v <- .rp_variances(se, se2)
  if (!is.null(site)) {
    labels <- .rp_check_site(site, length(estimate))
  }
  if (!is.null(X)) {
    rank <- qr(X)$rank
    if (rank < ncol(X)) {
      .rp_abort(
        "The design matrix of `mean_model` is not of full rank.",
        info = sprintf("It has %d columns and rank %d.", ncol(X), rank),
        fix = "Drop the columns that repeat others."
      )
    }
  }
  list(estimate = unname(estimate), se2 = unname(v), site = labels, X = X,
       mean_model = mean_model, omitted = sum(missing_row))
}

# The design matrix of a one-sided formula, or NULL for a common mean.
.rp_model_matrix <- function(mean_model, data, J) {
  if (!inherits(mean_model, "formula") || length(mean_model) != 2L) {
    .rp_abort("`mean_model` must be a one-sided formula such as `~ 1` or `~ x`.",
              info = .rp_you_passed(mean_model))
  }
  terms <- stats::terms(mean_model)
  if (length(attr(terms, "term.labels")) == 0L) {
    if (attr(terms, "intercept") == 1L) {
      return(NULL)
    }
    .rp_abort("`mean_model` must have at least one column; `~ 0` has none.")
  }
  frame <- tryCatch(
    stats::model.frame(mean_model, data = if (is.null(data)) environment(mean_model) else data,
                       na.action = stats::na.pass),
    error = function(err) {
      .rp_abort("The variables of `mean_model` could not be found.", parent = err)
    }
  )
  if (nrow(frame) != J) {
    .rp_abort(
      sprintf("The variables of `mean_model` must have one value for each of the %d sites.", J),
      info = sprintf("They have %d.", nrow(frame))
    )
  }
  stats::model.matrix(terms, frame)
}

# Fits the posterior and builds the object.
.rp_new_posterior <- function(inputs, prior, level, threshold, labels, control, call) {
  J <- length(inputs$estimate)
  p <- if (is.null(inputs$X)) 1L else ncol(inputs$X)
  need <- .rp_min_sites(p, prior)
  if (J < need) {
    if (prior$type == "flat") {
      .rp_abort(
        sprintf("The flat prior on sigma_tau needs at least %d sites; there %s %d.", need,
                if (J == 1L) "is" else "are", J),
        info = sprintf("The posterior is proper only when J > %d (Lee, 2026).", p + 1L),
        fix = "Use a proper prior, for example `prior = sd_prior(\"half_normal\", scale = 0.5)`."
      )
    }
    .rp_abort(
      sprintf("At least %d sites are needed; there %s %d.", need, if (J == 1L) "is" else "are", J),
      info = "With as many sites as columns of the mean model, the posterior equals the prior."
    )
  }
  if (J < 3L) {
    .rp_inform(
      sprintf("With %d sites the posterior of sigma_tau is close to its prior.", J),
      subclass = "rp_message_few_sites"
    )
  }
  if (control$engine == "adaptive" && J > 500L) {
    .rp_inform(
      sprintf("The adaptive engine computes each summary of the %d sites as a separate integral, which takes a while.", J),
      info = "`posterior_control(engine = \"nodes\")` is faster, and the two agree within 1e-9 on I, the mean reliability and H.",
      subclass = "rp_message_slow_engine"
    )
  }
  kernel <- .rp_kernel(inputs$estimate, inputs$se2, inputs$X, prior)
  law <- .rp_law(kernel, control$engine, control$rel_tol, control$abs_tol)
  turning <- .rp_turning_points(kernel$d)
  if (!turning$agree && !turning$homogeneous) {
    .rp_abort(
      "The turning points of H were not found consistently.",
      info = sprintf("Two scans differ by %s on the logit scale.", format(turning$difference)),
      class = "rp_numeric_error", subclass = "rp_turning_error"
    )
  }
  branches <- .rp_H_branches(kernel, turning)
  fitted <- .rp_posterior_summaries(kernel, law, branches, level, threshold)
  sites <- .rp_site_moments(law, kernel, control$chunk_sites)
  sigma_q <- fitted$sigma_quantiles
  S_lower <- sigma_q[[1L]]^2 / (sigma_q[[1L]]^2 + kernel$v)
  S_upper <- sigma_q[[3L]]^2 / (sigma_q[[3L]]^2 + kernel$v)
  checks <- .rp_numeric_checks(kernel, law, turning)
  out <- list(
    data = list(estimate = inputs$estimate, se2 = inputs$se2, site = inputs$site, X = inputs$X,
                mean_model = inputs$mean_model, omitted = inputs$omitted),
    se2 = inputs$se2,
    prior = prior, control = control,
    rule = if (law$engine == "nodes") law$rule,
    summary = fitted$summary,
    quantiles = fitted$quantiles,
    sigma_quantiles = sigma_q,
    probability = fitted$probability,
    plugin = fitted$plugin,
    sites = data.frame(site = inputs$site, estimate = inputs$estimate, se = sqrt(inputs$se2),
                       reliability = sites$reliability, reliability_lower = S_lower,
                       reliability_upper = S_upper, kl_term = sites$kl_term,
                       effect = sites$effect, effect_sd = sites$effect_sd),
    h_branches = list(turning_logits = turning$roots,
                      turning_sigma = .rp_sigma_of_m(kernel, turning$roots),
                      H_max = branches$H_max, homogeneous = turning$homogeneous),
    moments = fitted$moments,
    checks = checks,
    level = level, threshold = threshold, labels = labels, call = call, version = .rp_version()
  )
  out <- structure(out, class = "reliability_posterior")
  if (any(checks$status != "ok")) {
    .rp_warn(
      "Some checks of the numerical integration are outside their limits.",
      info = "See `check_numerics()`.",
      subclass = "rp_warning_numerics"
    )
  }
  out
}

# Posterior means, standard deviations and quantiles of the diagnostics,
# the probability above the threshold, and the values at E(sigma_tau^2).
.rp_posterior_summaries <- function(kernel, law, branches, level, threshold) {
  J <- kernel$J
  p <- kernel$p
  prior <- kernel$prior
  exists <- function(q) .rp_moment_exists(q, J, p, prior)
  # Only the moments that exist are integrated: the adaptive engine cannot
  # integrate the others, and the node rule would return the integral over
  # its finite range instead.
  bounded <- c("I", "mean", "harmonic", "H", "J_eff")
  # Means first, then the second moments about them, which keeps the
  # standard deviation of a quantity that hardly varies.
  means <- law$expect_matrix(function(s) .rp_functionals_at(kernel, s)[, bounded, drop = FALSE],
                             length(bounded))
  names(means) <- bounded
  second <- law$second()
  centered <- second$expect_matrix(function(s) {
    f <- .rp_functionals_at(kernel, s)[, bounded, drop = FALSE]
    sweep(f, 2L, means)^2
  }, length(bounded))
  names(centered) <- paste0(bounded, "_var")
  # The moments of sigma are taken in units of c = sqrt(GM(v)), u = sigma / c,
  # and converted, so that sigma^4 stays in range whatever the units. A
  # moment of order 2 exists only if that of order 1 does.
  c0 <- kernel$scale
  sigma_moments <- numeric()
  sd_sigma <- c(sigma = NA_real_, sigma2 = NA_real_)
  if (exists(1)) {
    u1 <- law$expect(function(s) s / c0)
    sigma_moments[["sigma"]] <- c0 * u1
  }
  if (exists(2)) {
    u2 <- law$expect(function(s) (s / c0)^2)
    u_var <- second$expect(function(s) (s / c0 - u1)^2)
    sigma_moments[["sigma2"]] <- c0^2 * u2
    sigma_moments[["sigma_var"]] <- c0^2 * u_var
    sd_sigma[["sigma"]] <- c0 * sqrt(max(0, u_var))
  }
  if (exists(4)) {
    u2_var <- second$expect(function(s) ((s / c0)^2 - u2)^2)
    sigma_moments[["sigma2_var"]] <- c0^4 * u2_var
    sd_sigma[["sigma2"]] <- c0^2 * sqrt(max(0, u2_var))
  }
  if (branches$H_max == 0) {
    # All the sampling variances are equal to rounding: H is 0 at every
    # sigma_tau, and the integral of its rounding error is not reported.
    means[["H"]] <- 0
    centered[["H_var"]] <- 0
  }
  moments <- c(means, centered, sigma_moments)
  mean_of <- function(name) moments[[name]]
  sd_of <- function(name) {
    if (name %in% names(sd_sigma)) sd_sigma[[name]] else sqrt(max(0, moments[[paste0(name, "_var")]]))
  }
  probs <- c((1 - level) / 2, 0.5, (1 + level) / 2)
  sigma_q <- law$quantile(probs)
  fq <- .rp_functionals_at(kernel, sigma_q)
  H_q <- .rp_H_quantile(law, kernel, branches, probs)
  # The mean of the effects, or the coefficients of the mean model; then
  # their variances, as the mean of the conditional variance and of the
  # squared deviation of the conditional mean from the mean, which keeps the
  # digits of a coefficient far from 0.
  coef_mean <- if (exists(1)) {
    law$expect_matrix(function(s) .rp_coef_at(kernel, s, centered = TRUE)[, seq_len(p), drop = FALSE], p)
  } else {
    rep(NA_real_, p)
  }
  coef_sd <- if (exists(2)) {
    sqrt(pmax(0, second$expect_matrix(function(s) {
      cf <- .rp_coef_at(kernel, s, centered = TRUE)
      sweep(cf[, seq_len(p), drop = FALSE], 2L, coef_mean)^2 + cf[, p + seq_len(p), drop = FALSE]
    }, p)))
  } else {
    rep(NA_real_, p)
  }
  coef_mean <- coef_mean + kernel$coef_offset
  coef_names <- if (is.null(kernel$X)) "effect_mean" else paste0("coef_", colnames(kernel$X))
  coef_q <- matrix(.rp_coef_intervals(law, kernel, probs), nrow = p)
  rows <- c("I", "mean", "harmonic", "H", "J_eff")
  summary <- data.frame(
    quantity = c(rows, "V", "sigma_tau", "sigma_tau2", coef_names),
    mean = c(vapply(rows, mean_of, 0), mean(kernel$d^2),
             if (exists(1)) mean_of("sigma") else NA_real_,
             if (exists(2)) mean_of("sigma2") else NA_real_,
             coef_mean),
    sd = c(vapply(rows, sd_of, 0), 0,
           if (exists(2)) sd_of("sigma") else NA_real_,
           if (exists(4)) sd_of("sigma2") else NA_real_,
           coef_sd),
    lower = c(fq[1L, c("I", "mean", "harmonic")], H_q[[1L]], fq[1L, "J_eff"], mean(kernel$d^2),
              sigma_q[[1L]], sigma_q[[1L]]^2, coef_q[, 1L]),
    upper = c(fq[3L, c("I", "mean", "harmonic")], H_q[[3L]], fq[3L, "J_eff"], mean(kernel$d^2),
              sigma_q[[3L]], sigma_q[[3L]]^2, coef_q[, 3L])
  )
  rownames(summary) <- NULL
  # A quantity that hardly varies with sigma_tau can have its ends a unit in
  # the last place out of order.
  ends <- cbind(summary$lower, summary$upper)
  summary$lower <- pmin(ends[, 1L], ends[, 2L])
  summary$upper <- pmax(ends[, 1L], ends[, 2L])
  quantiles <- data.frame(
    quantity = c(rows, "sigma_tau", "sigma_tau2", coef_names),
    lower = c(fq[1L, c("I", "mean", "harmonic")], H_q[[1L]], fq[1L, "J_eff"], sigma_q[[1L]],
              sigma_q[[1L]]^2, coef_q[, 1L]),
    median = c(fq[2L, c("I", "mean", "harmonic")], H_q[[2L]], fq[2L, "J_eff"], sigma_q[[2L]],
               sigma_q[[2L]]^2, coef_q[, 2L]),
    upper = c(fq[3L, c("I", "mean", "harmonic")], H_q[[3L]], fq[3L, "J_eff"], sigma_q[[3L]],
              sigma_q[[3L]]^2, coef_q[, 3L])
  )
  rownames(quantiles) <- NULL
  sigma_c <- sqrt(exp(kernel$log_gm) * threshold / (1 - threshold))
  probability <- c(I_above = law$cdf(sigma_c, lower = FALSE))
  plugin <- if (exists(2)) {
    a <- mean_of("sigma2")
    m <- log(a) - kernel$log_gm
    st <- .rp_profile_stats(m, kernel$d)
    c(sigma2 = a, st[c("I", "mean", "H", "J_eff")])
  } else {
    c(sigma2 = NA_real_, I = NA_real_, mean = NA_real_, H = NA_real_, J_eff = NA_real_)
  }
  moment_table <- data.frame(
    moment = c("E(sigma_tau)", "E(sigma_tau^2)", "Var(sigma_tau^2)"),
    exists = c(exists(1), exists(2), exists(4)),
    value = c(if (exists(1)) mean_of("sigma") else NA_real_,
              if (exists(2)) mean_of("sigma2") else NA_real_,
              if (exists(4)) moments[["sigma2_var"]] else NA_real_)
  )
  list(summary = summary, quantiles = quantiles, sigma_quantiles = sigma_q,
       probability = probability, plugin = plugin, moments = moment_table)
}

# Central intervals of the mean (or of each coefficient): quantiles of the
# normal mixture over sigma of N(beta_sigma, (X'WX)^-1).
.rp_coef_intervals <- function(law, kernel, probs) {
  nodes <- .rp_site_nodes(law$second())
  p <- kernel$p
  cf <- .rp_coef_at(kernel, nodes$sigma, centered = TRUE)
  center <- cf[, seq_len(p), drop = FALSE]
  sd <- sqrt(cf[, p + seq_len(p), drop = FALSE])
  vapply(probs, function(pp) .rp_mixture_quantile(nodes$weight, center, sd, pp) + kernel$coef_offset,
         numeric(p))
}
