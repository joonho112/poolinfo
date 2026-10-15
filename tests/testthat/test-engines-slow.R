# Slow comparisons, run only with POOLINFO_SLOW=true: the two engines on the
# six OCRS analyses and on 500 random inputs, and the node rule on 400
# extreme inputs. They take several minutes.

skip_slow <- function() {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("POOLINFO_SLOW"), "true"),
              "Set POOLINFO_SLOW=true to run the slow comparisons.")
}

# The largest differences between the engines, each divided by its
# tolerance, so that a value of at most 1 passes.
engine_differences <- function(y, v) {
  a <- reliability_posterior(y, se2 = v, control = posterior_control(engine = "adaptive"))
  n <- reliability_posterior(y, se2 = v, control = posterior_control(engine = "nodes"))
  get <- function(fit, q, col) fit$summary[[col]][fit$summary$quantity == q]
  cols <- c("mean", "lower", "upper")
  both <- function(q) vapply(cols, function(col) c(get(a, q, col), get(n, q, col)), c(0, 0))
  bounded <- max(abs(unlist(lapply(c("I", "mean", "harmonic", "H"), function(q) {
    b <- both(q)
    b[1, ] - b[2, ]
  }))), abs(a$probability - n$probability), abs(a$sites$reliability - n$sites$reliability))
  jeff <- both("J_eff")
  scale_pairs <- rbind(t(both("sigma_tau")), t(both("effect_mean")),
                       cbind(a$sites$effect, n$sites$effect),
                       cbind(a$sites$effect_sd, n$sites$effect_sd))
  scale_tol <- vapply(scale_pairs[, 2L], function(value) tol_effect_scale(v, value), 0)
  s2 <- both("sigma_tau2")
  c(
    bounded = bounded / tol_bounded,
    jeff = max(rel_diff(jeff[1, ], jeff[2, ])) / tol_jeff_rel,
    scale = max(abs(scale_pairs[, 1L] - scale_pairs[, 2L]) / scale_tol, na.rm = TRUE),
    sigma2 = max(rel_diff(s2[1, ], s2[2, ]), abs(a$plugin[["I"]] - n$plugin[["I"]]), na.rm = TRUE) /
      tol_sigma2_rel,
    checks = as.numeric(!all(a$checks$status == "ok") || !all(n$checks$status == "ok"))
  )
}

test_that("the two engines agree on the OCRS analyses and on 500 random inputs", {
  skip_slow()
  cases <- lapply(ocrs_analyses(), function(z) list(y = z$y, v = z$v))
  withr::local_seed(20261014)
  for (i in 1:500) {
    J <- sample(c(3:12, 15, 20, 30, 50, 80, 120, 200), 1)
    v <- exp(stats::rnorm(J, log(stats::runif(1, 0.001, 10)), stats::runif(1, 0.05, 2.5)))
    sigma <- sqrt(exp(mean(log(v)))) * exp(stats::rnorm(1, 0, 1.5))
    y <- stats::rnorm(J, stats::rnorm(1), sqrt(sigma^2 + v))
    if (stats::runif(1) < 0.1) y[1] <- y[1] + 10 * sqrt(max(v))
    cases[[length(cases) + 1L]] <- list(y = y, v = v)
  }
  d <- vapply(cases, function(z) engine_differences(z$y, z$v), numeric(5))
  expect_lte(max(d["bounded", ]), 1)
  expect_lte(max(d["jeff", ]), 1)
  expect_lte(max(d["scale", ]), 1)
  expect_lte(max(d["sigma2", ]), 1)
  expect_identical(sum(d["checks", ]), 0)
})

test_that("the node rule gives finite, ordered summaries for 400 extreme inputs", {
  skip_slow()
  withr::local_seed(20261015)
  for (i in 1:400) {
    kind <- sample(c("ratio", "outlier", "equal", "tiny", "huge", "big_J"), 1,
                   prob = c(0.3, 0.2, 0.15, 0.1, 0.1, 0.15))
    J <- if (kind == "big_J") sample(c(1000L, 3000L, 10000L), 1, prob = c(0.6, 0.3, 0.1)) else
      sample(c(3:10, 15, 25, 50, 100, 300), 1)
    base <- exp(stats::runif(1, -12, 8))
    v <- switch(kind,
      ratio = base * exp(stats::runif(J, 0, log(10^stats::runif(1, 0, 8)))),
      equal = base * (1 + stats::runif(J, -1, 1) * 10^stats::runif(1, -12, -4)),
      base * exp(stats::rnorm(J, 0, stats::runif(1, 0.1, 2))))
    sigma <- sqrt(exp(mean(log(v)))) * exp(stats::rnorm(1, 0, 2))
    y <- stats::rnorm(J, stats::rnorm(1, 0, 10 * sqrt(base)), sqrt(sigma^2 + v))
    if (kind == "outlier") y[sample(J, max(1, J %/% 20))] <- y[1] + 50 * sqrt(max(v)) * sample(c(-1, 1), 1)
    if (kind == "tiny") y <- stats::rnorm(J, 0, sqrt(v))
    if (kind == "huge") y <- y * 1e3
    fit <- reliability_posterior(y, se2 = v)
    s <- fit$summary
    label <- sprintf("case %d (%s, %d sites)", i, kind, J)
    expect_true(all(is.finite(s$mean[s$quantity %in% c("I", "mean", "harmonic", "H", "J_eff", "V")])) &&
                  all(is.finite(s$lower)) && all(is.finite(s$upper)), label = label)
    expect_true(all(s$lower <= s$upper), label = label)
    expect_true(fit$probability[["I_above"]] >= 0 && fit$probability[["I_above"]] <= 1, label = label)
    expect_true(all(is.finite(as.matrix(fit$sites[, -1L]))), label = label)
    expect_true(all(fit$checks$status == "ok"), label = label)
  }
})
