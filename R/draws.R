#' Diagnostics from posterior draws of the between-site SD
#'
#' Computes the diagnostics of the reliability profile at each posterior
#' draw of the between-site standard deviation from another fit, such as a
#' model fitted with Stan or brms, and summarizes them with their Monte Carlo
#' error.
#'
#' @param x Posterior draws of the between-site standard deviation: a
#'   numeric vector, an object that the posterior package can convert to
#'   draws (from Stan, CmdStan or brms, among others), or a fit of brms, from
#'   which the standard errors, the estimates and the sites are taken too.
#' @param se,se2 Standard errors or sampling variances of the site
#'   estimates, as in [reliability_profile()]. Give one of them.
#' @param site Optional labels of the sites.
#' @param estimate Optional estimates of the sites. With `mean`, the site
#'   table has the posterior means and intervals of the site effects.
#' @param mean Optional draws of the mean of the effects, in the same order
#'   as the draws of `x`: a numeric vector, or with a draws object the name
#'   of its variable.
#' @param variable Name of the variable of the standard deviation in a draws
#'   object, a string. `NULL`, the default, takes the one variable whose name
#'   starts with `sd_` or is `sigma_tau`, `tau` or `sigma`; when there is
#'   not exactly one such variable, `variable` must be given.
#' @param model_label A description of the model and the prior that gave the
#'   draws, a string, printed with the summaries.
#' @param scale_label,target_label Optional strings that name the scale of
#'   the effects and the estimand at each site, as in
#'   [reliability_posterior()]. They are printed and used by [report_text()]
#'   and [report_table()].
#' @param level Probability content of the posterior intervals, a number
#'   between 0 and 1.
#' @param threshold Value of the informativeness whose posterior probability
#'   of being exceeded is reported, a number between 0 and 1.
#' @param data Optional data frame in which `se`, `se2`, `site` and
#'   `estimate` are found, as in [reliability_posterior()].
#' @param ... Must be empty. The arguments after `x` are named, because an
#'   unnamed vector could hold standard errors or sampling variances.
#'
#' @details
#' External-draw summaries are an experimental extension beyond Lee (2026);
#' their interfaces may change in a later 0.x release.
#'
#' The diagnostics depend on the draws only through the between-site
#' standard deviation, so they apply to any posterior algorithm whose draws
#' target the stated model and prior (Lee, 2026, Appendix D.3). Nothing is
#' refitted: the model and the prior are those of the fit that gave the
#' draws, and `model_label` names them for the printed output.
#'
#' The posterior means are draw averages and come with their Monte Carlo
#' standard errors. With the posterior package installed, the effective
#' sample sizes of the draws (`posterior::ess_mean()`), which account for the
#' autocorrelation of each chain, give the standard errors; otherwise the
#' draws are treated as independent, with a message. For the
#' integrals of the normal model with a flat or proper prior on the SD,
#' [reliability_posterior()] has no Monte Carlo error.
#'
#' A univariate Gaussian identity-link fit of brms of the form
#' `y | se(sei, sigma = FALSE) ~ 1 + (1 | group)` is read directly.
#' `sigma = FALSE` may be omitted. The response, known standard errors and
#' grouping variable must be direct data columns, with one observation per
#' group. The standard deviation comes from `sd_<group>__Intercept` and the
#' common mean from `b_Intercept`. Extra residual variance, other families,
#' additional group-level components, random slopes, distributional formulas,
#' autocorrelation, smooths, Gaussian processes and latent measurement-error
#' terms are not supported.
#'
#' Gaussian models with fixed moderators can supply variance-based
#' diagnostics and site reliabilities. Site effect estimates are omitted:
#' a single common-mean draw cannot represent their site-specific population
#' predictions. Giving `mean` for such a fit is an error.
#'
#' The site summaries use the draws of the mean rather than the flat-prior
#' formula: the posterior mean of site \eqn{j}'s effect averages
#' \eqn{\mu + S_j(\hat\tau_j - \mu)}{mu + S_j (tau_hat_j - mu)} over the
#' draws, and its interval is a quantile of the normal mixture over them.
#'
#' @return An object of class `reliability_draws`, a list whose main
#'   components are
#'   * `summary`: a data frame with the mean, standard deviation, interval,
#'     Monte Carlo standard error (`mcse`) and
#'     effective sample size of the draws (`ess`) of each diagnostic, of
#'     `sigma_tau` and `sigma_tau^2`, and of the mean when its draws are given.
#'   * `draws`: the draws of `sigma_tau` and of the mean.
#'   * `probability` and `probability_mcse`: the proportion of draws with
#'     \eqn{I} above `threshold`, and its Monte Carlo standard error.
#'   * `plugin`: the diagnostics at the mean of the draws of `sigma_tau^2`.
#'   * `se2`, `site`, `estimate`, `level`, `threshold`, `model_label` and
#'     `call`: the inputs.
#'
#'   [profile_summary()], [site_table()], [prob_above()], [plugin_gap()] and
#'   `plot()` take it.
#' @seealso [reliability_posterior()], whose integrals have no Monte Carlo
#'   error, and `vignette("implementation")` for a brms example.
#'
#' @template ref-lee2026
#' @examples
#' # Draws of sigma_tau from some fit; here simulated for the example.
#' set.seed(1)
#' sigma_draws <- abs(rnorm(4000, 8, 4))
#' rd <- reliability_draws(sigma_draws, se = eight_schools$se,
#'                         model_label = "an example")
#' rd
#' @export
reliability_draws <- function(x, ..., se = NULL, se2 = NULL, site = NULL, estimate = NULL,
                              data = NULL, mean = NULL, variable = NULL, model_label = NULL,
                              scale_label = NULL, target_label = NULL, level = 0.95,
                              threshold = 0.5) {
  .rp_check_named_variances(...)
  if (!is.null(data)) {
    if (inherits(x, "brmsfit")) {
      .rp_abort("With a brms fit, leave out `data`.", info = "The data are read from the fit.")
    }
    if (!is.data.frame(data)) {
      .rp_abort("`data` must be a data frame.", info = .rp_you_passed(data))
    }
    cols <- .rp_eval_columns(data, list(
      se = rlang::enquo(se), se2 = rlang::enquo(se2), site = rlang::enquo(site),
      estimate = rlang::enquo(estimate)
    ))
    se <- cols$se
    se2 <- cols$se2
    site <- cols$site
    estimate <- cols$estimate
  }
  if (inherits(x, "brmsfit")) {
    parts <- .rp_brms_parts(x)
    data <- x$data
    if (is.null(se) && is.null(se2)) se <- data[[parts$se]]
    if (is.null(parts$mean_variable)) {
      if (!is.null(mean)) {
        .rp_abort(
          "A brms fit with fixed moderators supplies variance-based diagnostics only.",
          info = "A common `mean` draw cannot represent site-specific population predictions.",
          class = "rp_scope_error"
        )
      }
      estimate <- NULL
    } else if (is.null(estimate)) {
      estimate <- data[[parts$response]]
    }
    if (is.null(site)) site <- as.character(data[[parts$group]])
    if (is.null(variable)) variable <- parts$sd_variable
    if (is.null(mean) && !is.null(parts$mean_variable)) mean <- parts$mean_variable
    if (is.null(model_label)) model_label <- "brms"
  }
  extracted <- .rp_extract_draws(x, variable, mean)
  sigma <- extracted$sigma
  mean_draws <- extracted$mean
  if (!is.null(mean_draws) && length(mean_draws) != length(sigma)) {
    .rp_abort(
      "`mean` must have one draw for each draw of the standard deviation.",
      info = sprintf("There are %d draws of the standard deviation and %d of the mean.",
                     length(sigma), length(mean_draws))
    )
  }
  v <- .rp_variances(se, se2)
  J <- length(v)
  site <- .rp_site_labels(.rp_check_site(site, J), J, names(v))
  if (!is.null(estimate)) {
    estimate <- .rp_check_vector(estimate, "estimate")
    if (length(estimate) != J) {
      .rp_abort(sprintf("`estimate` must hold one number for each of the %d sites.", J),
                info = .rp_you_passed(estimate))
    }
  }
  level <- .rp_check_level(level)
  threshold <- .rp_check_level(threshold, "threshold")
  model_label <- .rp_check_string(model_label, "model_label")
  kernel <- .rp_draws_kernel(v)
  f <- .rp_functionals_at(kernel, sigma)
  # With equal sampling variances (to rounding) H is 0 at every sigma_tau, as
  # in a posterior fit; its computed values are rounding.
  if (mean(kernel$d^2) < .RP_V_HOMOGENEOUS) f[, "H"] <- 0
  chains <- extracted$chains
  ess <- .rp_ess(f[, c("I", "mean", "harmonic", "H", "J_eff", "sigma", "sigma2"), drop = FALSE], chains)
  if (!extracted$has_chains && !requireNamespace("posterior", quietly = TRUE)) {
    .rp_inform_once(
      "draws_independent",
      "The Monte Carlo errors treat the draws as independent.",
      info = "Install the posterior package to account for their autocorrelation.",
      subclass = "rp_message_independent_draws"
    )
  }
  probs <- c((1 - level) / 2, 0.5, (1 + level) / 2)
  rows <- c("I", "mean", "harmonic", "H", "J_eff", "sigma", "sigma2")
  quantity <- c("I", "mean", "harmonic", "H", "J_eff", "sigma_tau", "sigma_tau2")
  qs <- t(apply(f[, rows, drop = FALSE], 2L, stats::quantile, probs = probs, names = FALSE))
  means <- colMeans(f[, rows, drop = FALSE])
  sds <- apply(f[, rows, drop = FALSE], 2L, stats::sd)
  summary <- data.frame(
    quantity = c(quantity[1:5], "V", quantity[6:7]),
    mean = c(means[1:5], mean(kernel$d^2), means[6:7]),
    sd = c(sds[1:5], 0, sds[6:7]),
    lower = c(qs[1:5, 1L], mean(kernel$d^2), qs[6:7, 1L]),
    upper = c(qs[1:5, 3L], mean(kernel$d^2), qs[6:7, 3L]),
    mcse = c(sds[1:5] / sqrt(ess[1:5]), 0, sds[6:7] / sqrt(ess[6:7])),
    ess = c(ess[1:5], NA_real_, ess[6:7])
  )
  if (!is.null(mean_draws)) {
    m_ess <- .rp_ess(matrix(mean_draws, ncol = 1L), chains)
    summary <- rbind(summary, data.frame(
      quantity = "effect_mean", mean = base::mean(mean_draws), sd = stats::sd(mean_draws),
      lower = stats::quantile(mean_draws, probs[[1L]], names = FALSE),
      upper = stats::quantile(mean_draws, probs[[3L]], names = FALSE),
      mcse = stats::sd(mean_draws) / sqrt(m_ess), ess = m_ess
    ))
  }
  rownames(summary) <- NULL
  medians <- c(qs[1:5, 2L], qs[6:7, 2L])
  above <- f[, "I"] > threshold
  a2 <- base::mean(sigma^2)
  plug <- .rp_profile_stats(log(a2) - kernel$log_gm, kernel$d)
  out <- list(
    se2 = unname(v), site = site, estimate = if (!is.null(estimate)) unname(estimate),
    draws = data.frame(sigma_tau = sigma, mean = if (!is.null(mean_draws)) mean_draws else NA_real_),
    chains = chains,
    summary = summary,
    medians = stats::setNames(medians, c(quantity[1:5], quantity[6:7])),
    probability = c(I_above = base::mean(above)),
    probability_mcse = sqrt(base::mean(above) * (1 - base::mean(above)) / .rp_ess(matrix(as.numeric(above), ncol = 1L), chains)),
    plugin = c(sigma2 = a2, plug[c("I", "mean", "H", "J_eff")]),
    model_label = model_label, level = level, threshold = threshold,
    labels = list(scale_label = .rp_check_string(scale_label, "scale_label"),
                  target_label = .rp_check_string(target_label, "target_label")),
    call = match.call(), version = .rp_version()
  )
  structure(out, class = "reliability_draws")
}

# A kernel with only what the diagnostics need: the sampling variances.
.rp_draws_kernel <- function(v) {
  centered <- .rp_center_v(v)
  list(v = v, J = length(v), log_gm = centered$log_gm, d = centered$d,
       scale = exp(centered$log_gm / 2))
}

# The draws of the SD and of the mean, from a vector or a draws object, with
# the chain of each draw when it is known.
.rp_extract_draws <- function(x, variable, mean) {
  if (is.numeric(x) && is.null(dim(x)) && !inherits(x, "draws")) {
    sigma <- .rp_check_sd_draws(.rp_check_vector(x, "x"))
    if (!is.null(mean) && !is.numeric(mean)) {
      .rp_abort("With a vector of draws, `mean` must be a vector of draws too.",
                info = .rp_you_passed(mean))
    }
    return(list(sigma = sigma, mean = if (!is.null(mean)) .rp_check_vector(mean, "mean"),
                chains = rep(1L, length(sigma)), has_chains = FALSE))
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    .rp_abort(
      sprintf("Draws in %s need the posterior package.", .rp_describe(x)),
      fix = "Install posterior, or give the draws of the standard deviation as a numeric vector."
    )
  }
  draws <- tryCatch(
    posterior::as_draws_df(x),
    error = function(err) {
      .rp_abort(sprintf("%s could not be turned into draws.", .rp_describe(x)), parent = err)
    }
  )
  names_all <- posterior::variables(draws)
  if (is.null(variable)) {
    candidates <- grep("^(sd_|sigma_tau|tau$|sigma$)", names_all, value = TRUE)
    if (length(candidates) != 1L) {
      .rp_abort(
        "Name the variable of the between-site standard deviation with `variable`.",
        info = if (length(candidates) > 1L) {
          sprintf("Candidates: %s.", paste(sprintf("`%s`", utils::head(candidates, 8L)), collapse = ", "))
        } else {
          "None of the variables looks like one."
        }
      )
    }
    variable <- candidates
  }
  variable <- .rp_check_string(variable, "variable", null_ok = FALSE)
  if (!variable %in% names_all) {
    .rp_abort(sprintf("The draws have no variable `%s`.", variable),
              info = sprintf("They have %s.", .rp_column_list(names_all)))
  }
  sigma <- as.numeric(draws[[variable]])
  mean_draws <- NULL
  if (!is.null(mean)) {
    if (is.character(mean)) {
      if (!mean %in% names_all) {
        .rp_abort(sprintf("The draws have no variable `%s`.", mean))
      }
      mean_draws <- as.numeric(draws[[mean]])
    } else {
      mean_draws <- .rp_check_vector(mean, "mean")
    }
  }
  sigma <- .rp_check_sd_draws(.rp_check_vector(sigma, sprintf("the draws of `%s`", variable)))
  list(sigma = sigma, mean = mean_draws, chains = as.integer(draws[[".chain"]]), has_chains = TRUE)
}

# Draws of a standard deviation: at least two, none negative.
.rp_check_sd_draws <- function(sigma) {
  if (length(sigma) < 2L) {
    .rp_abort("At least 2 draws of the standard deviation are needed.",
              info = sprintf("There %s %d.", if (length(sigma) == 1L) "is" else "are", length(sigma)))
  }
  if (any(sigma < 0)) {
    .rp_abort("Draws of a standard deviation cannot be negative.",
              info = sprintf("%d of the %d draws are negative.", sum(sigma < 0), length(sigma)))
  }
  sigma
}

# Effective sample sizes of the columns of f, from posterior::ess_mean() when
# posterior is installed, and the number of draws otherwise.
.rp_ess <- function(f, chains) {
  n <- nrow(f)
  if (!requireNamespace("posterior", quietly = TRUE)) {
    return(rep(n, ncol(f)))
  }
  n_chains <- length(unique(chains))
  per_chain <- n %/% n_chains
  apply(f, 2L, function(col) {
    if (!is.finite(stats::sd(col)) || stats::sd(col) == 0) return(n)
    if (n_chains > 1L && all(table(chains) == per_chain)) {
      m <- matrix(col[order(chains)], nrow = per_chain)
      return(posterior::ess_mean(m))
    }
    posterior::ess_mean(col)
  })
}

#' @export
format.reliability_draws <- function(x, ..., use_cli = NULL) {
  .rp_check_dots_empty(...)
  s <- x$summary
  J <- length(x$se2)
  source <- if (is.null(x$model_label)) "posterior draws" else paste("posterior draws from", x$model_label)
  heading <- .rp_wrap(sprintf("Reliability profile of %s (%s)", .rp_sites_word(J), source), exdent = 2L)
  val <- function(q, col) s[[col]][s$quantity == q]
  row <- function(label, q, f, tail = "") {
    sprintf("  %-26s %-8s %-10s %-8s %-8s%s", label, f(val(q, "mean")),
            sprintf("(%s)", .rp_fmt_sig(val(q, "mcse"), 2L)), f(val(q, "lower")), f(val(q, "upper")), tail)
  }
  f3 <- function(z) .rp_fmt(z, 3L)
  lines <- c(
    heading,
    sprintf("%d draws; Monte Carlo standard errors of the means in parentheses.", nrow(x$draws)),
    .rp_label_lines(x$labels),
    .rp_conditioning_lines(), "",
    sprintf("  %-26s %-19s %s", "", "mean", paste(.rp_percent(x$level), "interval")),
    row("informativeness I", "I", f3),
    row("mean reliability", "mean", f3),
    row("dispersion H", "H", function(z) .rp_fmt(z, 4L)),
    row("concentration count J_eff", "J_eff", function(z) .rp_fmt(z, 1L), sprintf(" of %d", J)),
    sprintf("  %-26s %-19s %s", "logit variance V", f3(val("V", "mean")), "(fixed by the standard errors)"),
    row("sigma_tau", "sigma_tau", function(z) .rp_fmt_sig(z)),
    if ("effect_mean" %in% s$quantity) row("effect mean", "effect_mean", function(z) .rp_fmt_sig(z)),
    "",
    sprintf("  %-26s %s (%s)", sprintf("Pr(I > %s)", format(x$threshold)),
            f3(x$probability[["I_above"]]), .rp_fmt_sig(x$probability_mcse, 2L)),
    sprintf("  %-26s %s", "I at E(sigma_tau^2)", f3(x$plugin[["I"]]))
  )
  .rp_style_lines(lines, bold = seq_along(heading), use_cli = use_cli)
}

#' @export
print.reliability_draws <- function(x, ..., use_cli = NULL) {
  cat(format(x, ..., use_cli = use_cli), sep = "\n")
  invisible(x)
}

#' @export
profile_summary.reliability_draws <- function(x, level = NULL, ...) {
  .rp_check_dots_empty(...)
  level <- .rp_check_level(level, null_ok = TRUE)
  out <- x$summary
  if (!is.null(level) && !isTRUE(all.equal(level, x$level))) {
    kernel <- .rp_draws_kernel(x$se2)
    f <- .rp_functionals_at(kernel, x$draws$sigma_tau)
    rows <- c(I = "I", mean = "mean", harmonic = "harmonic", H = "H", J_eff = "J_eff",
              sigma_tau = "sigma", sigma_tau2 = "sigma2")
    probs <- c((1 - level) / 2, (1 + level) / 2)
    for (q in names(rows)) {
      ends <- stats::quantile(f[, rows[[q]]], probs, names = FALSE)
      out$lower[out$quantity == q] <- ends[[1L]]
      out$upper[out$quantity == q] <- ends[[2L]]
    }
    if ("effect_mean" %in% out$quantity) {
      ends <- stats::quantile(x$draws$mean, probs, names = FALSE)
      out$lower[out$quantity == "effect_mean"] <- ends[[1L]]
      out$upper[out$quantity == "effect_mean"] <- ends[[2L]]
    }
  }
  attr(out, "level") <- if (is.null(level)) x$level else level
  attr(out, "model") <- x$model_label
  attr(out, "conditioning") <- paste(.rp_conditioning_lines(), collapse = " ")
  out
}

#' @export
site_table.reliability_draws <- function(x, level = NULL, ...) {
  .rp_check_dots_empty(...)
  level <- .rp_check_level(level, null_ok = TRUE)
  if (is.null(level)) level <- x$level
  v <- x$se2
  J <- length(v)
  sigma <- x$draws$sigma_tau
  n <- length(sigma)
  kernel <- .rp_draws_kernel(v)
  probs <- c((1 - level) / 2, (1 + level) / 2)
  q_sigma <- stats::quantile(sigma, probs, names = FALSE)
  out <- data.frame(site = x$site)
  if (!is.null(x$estimate)) out$estimate <- x$estimate
  out$se <- sqrt(v)
  reliability <- kl <- share <- numeric(J)
  for (start in seq.int(1L, J, by = 2000L)) {
    cols <- start:min(J, start + 1999L)
    S <- .rp_plogis(outer(2 * log(sigma), -log(v[cols]), "+"))
    reliability[cols] <- colMeans(S)
    kl[cols] <- colMeans(.rp_site_kl(kernel, sigma, cols))
  }
  out$reliability <- reliability
  out$reliability_lower <- .rp_plogis(2 * log(q_sigma[[1L]]) - log(v))
  out$reliability_upper <- .rp_plogis(2 * log(q_sigma[[2L]]) - log(v))
  out$kl_term <- kl
  out$share <- .rp_site_shares(kernel, sigma, rep(1 / n, n))
  if (!is.null(x$estimate) && !all(is.na(x$draws$mean))) {
    tau <- x$draws$mean
    y <- x$estimate
    em <- es <- lo <- hi <- numeric(J)
    for (start in seq.int(1L, J, by = 500L)) {
      cols <- start:min(J, start + 499L)
      eta <- outer(2 * log(sigma), -log(v[cols]), "+")
      S <- .rp_plogis(eta)
      center <- S * matrix(y[cols], n, length(cols), byrow = TRUE) + .rp_plogis(-eta) * tau
      sd <- sqrt(matrix(v[cols], n, length(cols), byrow = TRUE) * S)
      em[cols] <- colMeans(center)
      es[cols] <- sqrt(colMeans(sd^2) + apply(center, 2L, stats::var) * (n - 1) / n)
      w <- rep(1 / n, n)
      lo[cols] <- .rp_mixture_quantile(w, center, sd, probs[[1L]])
      hi[cols] <- .rp_mixture_quantile(w, center, sd, probs[[2L]])
    }
    out$effect <- em
    out$effect_sd <- es
    out$effect_lower <- lo
    out$effect_upper <- hi
  }
  attr(out, "level") <- level
  out
}

#' @export
prob_above.reliability_draws <- function(x, threshold = NULL,
                                         which = c("I", "mean", "J_eff", "sigma_tau")) {
  which <- .rp_match(which, c("I", "mean", "J_eff", "sigma_tau"), "which")
  kernel <- .rp_draws_kernel(x$se2)
  threshold <- .rp_prob_threshold(threshold, which, x$threshold, kernel$J)
  f <- .rp_functionals_at(kernel, x$draws$sigma_tau)
  col <- if (which == "sigma_tau") "sigma" else which
  above <- f[, col] > threshold
  p <- mean(above)
  attr(p, "mcse") <- sqrt(p * (1 - p) / .rp_ess(matrix(as.numeric(above), ncol = 1L), x$chains))
  p
}

#' @export
plugin_gap.reliability_draws <- function(x) {
  names <- c("I", "mean", "H", "J_eff")
  post <- vapply(names, function(q) x$summary$mean[x$summary$quantity == q], 0)
  plug <- x$plugin[names]
  out <- data.frame(functional = names, posterior_mean = unname(post),
                    at_mean_sigma2 = unname(plug), difference = unname(plug - post))
  attr(out, "sigma2") <- x$plugin[["sigma2"]]
  out
}

#' @rdname reliability_profile
#' @order 6
#' @export
reliability_profile.reliability_draws <- function(x, ...) {
  .rp_check_dots_empty(...)
  .rp_new_profile(
    se = NULL, se2 = x$se2, sigma_tau = sqrt(x$plugin[["sigma2"]]), site = x$site,
    estimate = x$estimate, source = "mean of the draws of sigma_tau^2",
    scale_label = x$labels$scale_label, target_label = x$labels$target_label, call = match.call()
  )
}
