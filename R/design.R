#' Sampling variances of a planned trial
#'
#' Computes the planning approximation of the sampling variance of each site's
#' impact estimate in an individually randomized multisite trial.
#'
#' @param n Planned numbers of individuals in the sites, a numeric vector of
#'   positive numbers.
#' @param p Proportion of each site's individuals assigned to treatment,
#'   numbers strictly between 0 and 1, one for all sites or one per site.
#' @param R2 Proportion of the outcome variance explained by covariates,
#'   numbers in \eqn{[0, 1)}, one for all sites or one per site.
#' @param var_outcome Variance of the outcome within each arm of a site
#'   (\eqn{s_j^2}{s_j^2} in equation 2 of Lee, 2026), positive
#'   numbers, one for all sites or one per site.
#'
#' @details
#' The variance of a difference in means with covariate adjustment is
#' approximately
#' \deqn{v_j = \Big(\frac{1}{p_j} + \frac{1}{1 - p_j}\Big)\frac{s_j^2 (1 - R^2)}{n_j} = \frac{\kappa_j}{n_j}}{v_j = (1/p_j + 1/(1 - p_j)) s_j^2 (1 - R^2) / n_j = kappa_j / n_j}
#' (Lee, 2026, equation 2). The argument names are those of
#' `multisiteDGP::compute_kappa()`. Designs with blocks, clusters or other
#' estimators need their own variance formula; [reliability_design()] also
#' takes the constants \eqn{\kappa_j} directly.
#'
#' @return A numeric vector of sampling variances, one per site.
#'
#' @template ref-lee2026
#' @examples
#' neyman_variance(c(20, 50, 100))
#' neyman_variance(c(20, 50, 100), p = 0.6, R2 = 0.4)
#' @export
neyman_variance <- function(n, p = 0.5, R2 = 0, var_outcome = 1) {
  n <- .rp_check_vector(n, "n", positive = TRUE)
  .rp_kappa(length(n), p, R2, var_outcome) / n
}

# The constants kappa_j = (1/p + 1/(1 - p)) (1 - R2) var_outcome, recycled to
# J sites.
.rp_kappa <- function(J, p, R2, var_outcome) {
  recycle <- function(x, arg, ok, rule) {
    x <- .rp_check_vector(x, arg)
    if (!length(x) %in% c(1L, J)) {
      .rp_abort(sprintf("`%s` must have one value or one value for each of the %d sites.", arg, J),
                info = .rp_you_passed(x))
    }
    if (!all(ok(x))) {
      .rp_abort(sprintf("`%s` must hold %s.", arg, rule), info = .rp_you_passed(x))
    }
    rep_len(x, J)
  }
  p <- recycle(p, "p", function(x) x > 0 & x < 1, "numbers strictly between 0 and 1")
  R2 <- recycle(R2, "R2", function(x) x >= 0 & x < 1, "numbers of at least 0 and below 1")
  var_outcome <- recycle(var_outcome, "var_outcome", function(x) x > 0, "positive numbers")
  (1 / p + 1 / (1 - p)) * (1 - R2) * var_outcome
}

#' Reliability profile of a planned allocation
#'
#' Computes the reliability profile of a planned allocation of individuals to
#' sites and compares it with the balanced allocation of the same total, or
#' with equal spending when the costs per individual differ.
#'
#' @param n Planned numbers of individuals in the sites, a numeric vector of
#'   positive numbers.
#' @param sigma_tau Planned between-site standard deviation of the effects, a
#'   positive number, or a numeric vector of them for a range of scenarios.
#'   The printed comparison uses the first.
#' @param kappa Optional constants \eqn{\kappa_j} with \eqn{v_j = \kappa_j / n_j},
#'   one for all sites or one per site. By default they come from `p`, `R2`
#'   and `var_outcome`, as in [neyman_variance()].
#' @inheritParams neyman_variance
#' @param cost Optional costs per individual in each site, positive numbers.
#'   The comparison is then with equal spending of the same budget.
#'
#' @details
#' With a fixed number of sites \eqn{J}, a fixed total \eqn{N} and a common
#' \eqn{\kappa}, the balanced allocation \eqn{n_j = N/J} maximizes the
#' informativeness and the mean reliability, and minimizes the benchmark
#' risk; there \eqn{I = \sigma_\tau^2 / (\sigma_\tau^2 + \kappa J / N)}{I = sigma_tau^2 / (sigma_tau^2 + kappa J / N)}
#' (Lee, 2026, Proposition 4). The informativeness depends on the sizes
#' through their geometric mean \eqn{\tilde n}{n_tilde}: with
#' \eqn{u_j = n_j / \bar n - 1}{u_j = n_j / n_bar - 1} and
#' \eqn{\mathrm{CV}_J^2 = J^{-1}\sum_j u_j^2}{CV_J^2 = mean(u_j^2)},
#' \eqn{\log(\tilde n / \bar n) = -\mathrm{CV}_J^2/2 + R}{log(n_tilde / n_bar) = -CV_J^2 / 2 + R}
#' with \eqn{|R| \le J^{-1}\sum_j |u_j|^3 / \{3(1 - r)^3\}}{|R| <= mean(|u_j|^3) / (3 (1 - r)^3)}
#' when every \eqn{|u_j| \le r < 1}{|u_j| <= r < 1} (Appendix E.2). The
#' component `local` contains this approximation and its bound.
#'
#' With unequal costs \eqn{c_j}, the allocation that maximizes the
#' informativeness for a budget \eqn{B} spends the same on every site,
#' \eqn{n_j = B / (J c_j)}{n_j = B / (J c_j)} (Appendix E.3). With
#' constants \eqn{\kappa_j} that differ, balance still maximizes the
#' informativeness, because \eqn{\mathrm{GM}(v) = \mathrm{GM}(\kappa)/\mathrm{GM}(n)}{GM(v) = GM(kappa) / GM(n)},
#' but it need not make the reliabilities equal or maximize the mean
#' reliability (Appendix E.3).
#'
#' These results concern the site profile at known inputs: the
#' informativeness, the mean reliability and the benchmark risk, not ranks
#' or the distribution of the effects.
#'
#' @return An object of class `reliability_design`, a list whose main
#'   components are
#'   * `n`, `kappa`, `se2`, `sigma_tau`, `cost`: the inputs, with the
#'     sampling variances `se2`.
#'   * `planned`, `balanced`: `reliability_profile` objects of the planned
#'     allocation and of the balanced one (or equal spending) at the first
#'     `sigma_tau`.
#'   * `balanced_n`: the sizes of the balanced allocation or of equal
#'     spending.
#'   * `efficiency`: the ratio of the geometric to the arithmetic mean of the
#'     planned sizes.
#'   * `local`: the approximation of Appendix E.2, a list with the
#'     coefficient of variation of the sizes `cv`, the logarithm of the
#'     efficiency `exact`, its approximation `leading`
#'     (\eqn{-\mathrm{CV}^2/2}{-CV^2 / 2}) and the bound on the remainder
#'     `bound` (`NA` when a size is twice the mean or more).
#'   * `scenarios`: a data frame with the informativeness of both
#'     allocations at each `sigma_tau`.
#'
#'   `print()`, [profile_summary()] (with a column `allocation`),
#'   [site_table()] (with the columns `n`, `kappa` and `n_balanced`),
#'   [reliability_profile()] (the planned profile), [report_text()] and
#'   [plot_design()] take it.
#'
#' @seealso [balanced_allocation()] and [solve_design()] for allocations and
#'   sizes, and `vignette("planning")` for a worked plan.
#'
#' @template ref-lee2026
#' @examples
#' n <- round(exp(qnorm(seq(0.5 / 40, 1 - 0.5 / 40, length.out = 40), 3.8, 0.47)))
#' design <- reliability_design(n, sigma_tau = 0.2, kappa = 4)
#' design
#'
#' # The informativeness of both allocations over a range of sigma_tau
#' reliability_design(n, sigma_tau = c(0.1, 0.2, 0.3), kappa = 4)$scenarios
#' @export
reliability_design <- function(n, sigma_tau, kappa = NULL, p = 0.5, R2 = 0, var_outcome = 1,
                               cost = NULL) {
  n <- .rp_check_vector(n, "n", positive = TRUE)
  J <- length(n)
  if (missing(sigma_tau)) .rp_abort_missing_sigma_tau()
  sigma_tau <- .rp_check_vector(sigma_tau, "sigma_tau", positive = TRUE)
  kappa <- if (is.null(kappa)) {
    .rp_kappa(J, p, R2, var_outcome)
  } else {
    k <- .rp_check_vector(kappa, "kappa", positive = TRUE)
    if (!length(k) %in% c(1L, J)) {
      .rp_abort(sprintf("`kappa` must have one value or one value for each of the %d sites.", J),
                info = .rp_you_passed(kappa))
    }
    rep_len(k, J)
  }
  if (!is.null(cost)) {
    cost <- .rp_check_vector(cost, "cost", positive = TRUE)
    if (!length(cost) %in% c(1L, J)) {
      .rp_abort(sprintf("`cost` must have one value or one value for each of the %d sites.", J),
                info = .rp_you_passed(cost))
    }
    cost <- rep_len(cost, J)
  }
  balanced_n <- if (is.null(cost)) rep(sum(n) / J, J) else sum(cost * n) / (J * cost)
  se2 <- kappa / n
  planned <- .rp_new_profile(NULL, se2, sigma_tau[[1L]], NULL, NULL, source = "planned")
  balanced <- .rp_new_profile(NULL, kappa / balanced_n, sigma_tau[[1L]], NULL, NULL,
                              source = if (is.null(cost)) "balanced" else "equal spending")
  I_at <- function(v, s) .rp_plogis(2 * log(s) - mean(log(v)))
  scenarios <- data.frame(
    sigma_tau = sigma_tau,
    I_planned = vapply(sigma_tau, function(s) I_at(se2, s), 0),
    I_balanced = vapply(sigma_tau, function(s) I_at(kappa / balanced_n, s), 0)
  )
  out <- list(
    n = n, kappa = kappa, se2 = se2, sigma_tau = sigma_tau, cost = cost,
    balanced_n = balanced_n, planned = planned, balanced = balanced,
    efficiency = exp(mean(log(n))) / mean(n),
    local = .rp_size_approximation(n),
    scenarios = scenarios,
    call = match.call(), version = .rp_version()
  )
  structure(out, class = "reliability_design")
}

# The local approximation of Appendix E.2: log(GM(n) / mean(n)) = -CV_J^2 / 2
# + R, with |R| <= mean(|u|^3) / (3 (1 - r)^3) when every |u_j| <= r < 1.
.rp_size_approximation <- function(n) {
  u <- n / mean(n) - 1
  cv2 <- mean(u^2)
  r <- max(abs(u))
  exact <- mean(log(n)) - log(mean(n))
  # The bound applies only when every size is within a factor (1 - r, 1 + r)
  # of the mean with r < 1.
  bound <- if (r < 1) mean(abs(u)^3) / (3 * (1 - r)^3) else NA_real_
  list(cv = sqrt(cv2), exact = exact, leading = -cv2 / 2, bound = bound)
}

#' @export
format.reliability_design <- function(x, ..., use_cli = NULL) {
  .rp_check_dots_empty(...)
  J <- length(x$n)
  N <- sum(x$n)
  kappa_txt <- if (length(unique(x$kappa)) == 1L) {
    sprintf("kappa = %s", .rp_fmt_sig(x$kappa[[1L]]))
  } else {
    sprintf("kappa from %s to %s", .rp_fmt_sig(min(x$kappa)), .rp_fmt_sig(max(x$kappa)))
  }
  count <- function(z) format(round(z), big.mark = ",")
  heading <- .rp_wrap(sprintf("Planned allocation: %s sites, %s individuals, %s, %s", count(J), count(N),
                              .rp_nobreak(kappa_txt),
                              .rp_nobreak(paste("sigma_tau =", .rp_fmt_sig(x$sigma_tau[[1L]])))),
                      exdent = 20L)
  head <- c(
    heading,
    sprintf("Site sizes from %s to %s (CV %s)", count(min(x$n)), count(max(x$n)), .rp_fmt(x$local$cv, 2L))
  )
  p <- x$planned$stats
  b <- x$balanced$stats
  other <- if (is.null(x$cost)) "balanced" else "equal spending"
  gm <- function(n) exp(mean(log(n)))
  rows <- c(
    sprintf("  %-26s %-10s %s", "", "planned", other),
    sprintf("  %-26s %-10s %s", "informativeness I", .rp_fmt(p[["I"]], 3L), .rp_fmt(b[["I"]], 3L)),
    sprintf("  %-26s %-10s %s", "mean reliability", .rp_fmt(p[["mean"]], 3L), .rp_fmt(b[["mean"]], 3L)),
    sprintf("  %-26s %-10s %s", "dispersion H", .rp_fmt(p[["H"]], 4L), .rp_fmt(b[["H"]], 4L)),
    sprintf("  %-26s %-10s %s", "concentration count J_eff", .rp_fmt(p[["J_eff"]], 1L), .rp_fmt(b[["J_eff"]], 1L)),
    sprintf("  %-26s %-10s %s", "geometric mean size", .rp_fmt(gm(x$n), 1L), .rp_fmt(gm(x$balanced_n), 1L))
  )
  note <- if (is.null(x$cost)) {
    paste("Balance maximizes I and the mean reliability at fixed J and N when every site has",
          "the same kappa and the same cost per individual", .rp_nobreak("(Lee, 2026)."))
  } else {
    paste("With costs per individual that differ, equal spending maximizes I for the same",
          "budget", .rp_nobreak("(Lee, 2026)."))
  }
  lines <- c(head, "", rows, "", .rp_wrap(note))
  .rp_style_lines(lines, bold = seq_along(heading), use_cli = use_cli)
}

#' @export
print.reliability_design <- function(x, ..., use_cli = NULL) {
  cat(format(x, ..., use_cli = use_cli), sep = "\n")
  invisible(x)
}

#' @export
profile_summary.reliability_design <- function(x, level = NULL, ...) {
  .rp_check_dots_empty(...)
  a <- profile_summary(x$planned)
  b <- profile_summary(x$balanced)
  out <- rbind(cbind(allocation = "planned", a), cbind(allocation = if (is.null(x$cost)) "balanced" else "equal spending", b))
  rownames(out) <- NULL
  out
}

#' @export
site_table.reliability_design <- function(x, level = NULL, ...) {
  .rp_check_dots_empty(...)
  out <- site_table(x$planned)
  cbind(out[, "site", drop = FALSE], n = x$n, kappa = x$kappa, n_balanced = x$balanced_n,
        out[, setdiff(names(out), "site"), drop = FALSE])
}

#' @export
as.data.frame.reliability_design <- function(x, row.names = NULL, optional = FALSE, ...) {
  site_table(x)
}

#' @export
reliability_profile.reliability_design <- function(x, ...) {
  .rp_check_dots_empty(...)
  x$planned
}

#' Site allocations under size and budget constraints
#'
#' Computes allocations under size and budget constraints. The continuous
#' allocation maximizes the informativeness. With integer sizes and equal
#' costs the allocation is optimal; with unequal costs the function uses a
#' greedy rounding rule that need not find the integer optimum.
#'
#' @param J Number of sites, a positive integer.
#' @param N Total number of individuals, a positive number. Give `N`, for
#'   equal costs per individual, or `cost` and `budget`, not both. It must be
#'   an integer when `integer = TRUE` and is then allocated in full.
#' @param cost Costs per individual in each site, positive numbers, one for
#'   all sites or one per site.
#' @param budget Maximum budget, a positive number, in the units of `cost`.
#' @param lower,upper Optional bounds on the site sizes, numbers, one for all
#'   sites or one per site.
#' @param integer If `TRUE` (the default), the sizes are integers.
#'
#' @details
#' The informativeness increases with the product of the site sizes when the
#' constants \eqn{\kappa_j} are fixed, so the continuous allocation maximizes
#' \eqn{\sum_j \log n_j}{sum(log(n_j))} (Lee, 2026, Proposition 4 and
#' Appendix E.3):
#'
#' * with equal costs and no bounds, every site gets \eqn{N/J}.
#' * with bounds \eqn{L_j \le n_j \le U_j}{L_j <= n_j <= U_j},
#'   \eqn{n_j = \min\{U_j, \max(L_j, t)\}}{n_j = min(U_j, max(L_j, t))}
#'   with the level \eqn{t} that uses up \eqn{N}.
#' * with costs \eqn{c_j} and a budget \eqn{B}, every site gets the same
#'   spending, \eqn{n_j = B / (J c_j)}{n_j = B / (J c_j)}, clipped to the
#'   bounds in the same way.
#'
#' With `integer = TRUE`, the continuous solution is rounded down and the
#' remaining individuals (or budget) go one at a time to the sites where an
#' individual adds most to \eqn{\sum_j \log n_j}{sum(log(n_j))} per unit of
#' cost. With equal costs the sizes then differ by at most 1 wherever the
#' bounds allow; with unequal costs some of the budget can be left over.
#' This greedy rule is optimal with equal costs, but is a heuristic with
#' unequal costs. For example, with costs `c(1, 3)` and budget 10 it returns
#' `c(7, 1)`, whereas `c(4, 2)` has a larger product of site sizes at the same
#' cost. Bounds are intersected with the positive integers before an integer
#' allocation is computed; an infeasible request gives an error. A budget
#' may also be left unused when every site has reached its upper bound.
#'
#' @return A numeric vector of `J` site sizes, with the attributes `total`
#'   (the number of individuals) and, with costs, `spent`.
#'
#' @template ref-lee2026
#' @seealso [reliability_design()] to compare a planned allocation with
#'   these, and `vignette("planning")`.
#'
#' @examples
#' balanced_allocation(J = 7, N = 100)
#' balanced_allocation(J = 5, N = 100, upper = c(10, 50, 50, 50, 50))
#' balanced_allocation(J = 4, cost = c(1, 1, 2, 4), budget = 400)
#' @export
balanced_allocation <- function(J, N = NULL, cost = NULL, budget = NULL, lower = NULL,
                                upper = NULL, integer = TRUE) {
  J <- .rp_check_whole(J, "J")
  integer <- .rp_check_flag(integer, "integer")
  if (!is.null(N) && (!is.null(cost) || !is.null(budget))) {
    .rp_abort("Give `N`, or `cost` and `budget`, not both.",
              info = "`N` is a total number of individuals at equal costs; `budget` is spent at the costs `cost`.")
  }
  if (is.null(N) && (is.null(cost) || is.null(budget))) {
    .rp_abort("Give the total `N`, or both `cost` and `budget`.")
  }
  per_site <- function(x, arg, default) {
    if (is.null(x)) return(rep(default, J))
    x <- .rp_check_vector(x, arg)
    if (!length(x) %in% c(1L, J)) {
      .rp_abort(sprintf("`%s` must have one value or one value for each of the %d sites.", arg, J),
                info = .rp_you_passed(x))
    }
    rep_len(x, J)
  }
  if (!is.null(N)) {
    N <- .rp_check_number(N, "N", "a single positive number", lower = 0, lower_open = TRUE)
    if (integer && N != floor(N)) {
      .rp_abort("`N` must be an integer when `integer = TRUE`.", info = .rp_you_passed(N))
    }
    c_j <- rep(1, J)
    B <- N
  } else {
    c_j <- per_site(cost, "cost", 1)
    if (any(c_j <= 0)) .rp_abort("`cost` must hold positive numbers.", info = .rp_you_passed(cost))
    B <- .rp_check_number(budget, "budget", "a single positive number", lower = 0, lower_open = TRUE)
  }
  L <- per_site(lower, "lower", if (integer) 1 else 0)
  U <- per_site(upper, "upper", Inf)
  if (any(L < 0) || any(U <= 0) || any(L > U)) {
    .rp_abort("The bounds must satisfy 0 <= lower <= upper with upper > 0.")
  }
  if (integer) {
    L <- pmax(1, ceiling(L))
    U <- floor(U)
    if (any(L > U)) {
      .rp_abort("The bounds allow no positive integer size for at least one site.")
    }
  }
  if (sum(c_j * L) > B * (1 + 1e-12)) {
    .rp_abort("The lower bounds already need more than the total.",
              info = sprintf("They need %s; the total is %s.", .rp_fmt_sig(sum(c_j * L), 6L), .rp_fmt_sig(B, 6L)))
  }
  capacity <- sum(c_j * U)
  if (!is.null(N) && capacity < B * (1 - 1e-12)) {
    .rp_abort("The upper bounds cannot take the whole total.",
              info = sprintf("They take at most %s; the total is %s.", .rp_fmt_sig(sum(c_j * U), 6L), .rp_fmt_sig(B, 6L)))
  }
  allocation_total <- min(B, capacity)
  # Equal spending t per site, clipped to the bounds: n_j = clip(t / c_j).
  spend <- function(t) sum(c_j * pmin(U, pmax(L, t / c_j)))
  t <- if (all(is.infinite(U)) && all(L == 0)) allocation_total / J else {
    hi <- max(c_j * pmin(U, allocation_total))
    if (!is.finite(hi)) hi <- allocation_total
    # spend() rises to the total at the root; doubling the upper end keeps
    # the root inside when it falls on the end, as with a single site.
    stats::uniroot(function(z) spend(z) - allocation_total,
                   c(0, 2 * max(hi, allocation_total)), tol = 1e-12 * allocation_total)$root
  }
  n <- pmin(U, pmax(L, t / c_j))
  if (integer) {
    n <- .rp_integer_allocation(n, c_j, allocation_total, L, U)
  }
  if (any(!is.finite(n)) || any(n < L | n > U) ||
      (integer && any(n < 1 | n != floor(n))) ||
      sum(c_j * n) > B * (1 + 1e-12) ||
      (!is.null(N) && (if (integer) sum(n) != N else abs(sum(n) - N) > 1e-10 * N))) {
    .rp_abort("The allocation could not be represented within the supplied constraints.",
              fix = "Check the total, costs and bounds, or rescale the costs and budget.")
  }
  attr(n, "total") <- sum(n)
  if (!is.null(cost)) attr(n, "spent") <- sum(c_j * n)
  n
}

# From the continuous solution down to whole numbers, then one individual at
# a time to the site with the largest gain in log size per unit of cost
# among those within their bounds and the budget.
.rp_integer_allocation <- function(n, c_j, B, L, U) {
  out <- pmin(U, pmax(L, floor(n + 1e-9)))
  left <- B - sum(c_j * out)
  repeat {
    can <- out + 1 > out & out + 1 <= U & c_j <= left + 1e-12 * B
    if (!any(can)) break
    gain <- ifelse(can, log1p(1 / pmax(out, 1e-12)) / c_j, -Inf)
    j <- which.max(gain)
    out[[j]] <- out[[j]] + 1
    left <- left - c_j[[j]]
  }
  out
}

#' Sizes that reach a target informativeness
#'
#' Computes, for a balanced design, the site size needed for a target
#' informativeness with a given number of sites, or the largest number of
#' sites that reaches the target with a given total.
#'
#' @param target_I Target informativeness, a number strictly between 0 and
#'   1.
#' @param sigma_tau Planned between-site standard deviation, a positive
#'   number, or a numeric vector of them for a range of scenarios.
#' @param kappa Constant \eqn{\kappa} of the sampling variances
#'   \eqn{v_j = \kappa / n_j}{v_j = kappa / n_j}, a positive number (see
#'   [neyman_variance()]).
#' @param J,N The number of sites or the total number of individuals: give
#'   one of them.
#'
#' @details
#' Under balance \eqn{I = \sigma_\tau^2 / (\sigma_\tau^2 + \kappa / \bar n)}{I = sigma_tau^2 / (sigma_tau^2 + kappa / n_bar)},
#' so the target needs \eqn{\bar n = (\kappa / \sigma_\tau^2)\, I/(1 - I)}{n_bar = (kappa / sigma_tau^2) I / (1 - I)}
#' individuals per site, and with a total \eqn{N} at most
#' \eqn{J = N \sigma_\tau^2 (1 - I) / (\kappa I)}{J = N sigma_tau^2 (1 - I) / (kappa I)}
#' sites, rounded down. Under balance the informativeness depends on the
#' size per site, not on the number of sites: with `J` given, the target
#' fixes the size per site, and with `N` given, it fixes how many sites the
#' total can be divided into. Because `sigma_tau` is unknown at the planning
#' stage, a vector of plausible values gives one row per scenario.
#'
#' At a fixed total, more sites lower the informativeness but make the mean
#' effect more precise: its variance under balance is
#' \eqn{\sigma_\tau^2/J + \kappa/N}{sigma_tau^2 / J + kappa / N} (Appendix
#' E.3). The two are different design questions.
#'
#' @return A data frame with one row per value of `sigma_tau` and the columns
#'   * `target_I`, `sigma_tau`, `kappa`: the inputs.
#'   * `J`, `n_per_site`, `N`: the design, with `n_per_site` not rounded.
#'   * `I`, with `N` given: the informativeness reached with the number of
#'     sites rounded down. `J`, `n_per_site` and `I` are `NA` when the total
#'     cannot reach the target even in a single site.
#'
#' @template ref-lee2026
#' @seealso [reliability_design()] and [plot_design()], and
#'   `vignette("planning")`.
#'
#' @examples
#' solve_design(target_I = 0.5, sigma_tau = c(0.1, 0.2, 0.3), kappa = 4, J = 40)
#' solve_design(target_I = 0.5, sigma_tau = 0.2, kappa = 4, N = 4000)
#' @export
solve_design <- function(target_I, sigma_tau, kappa, J = NULL, N = NULL) {
  target_I <- .rp_check_number(target_I, "target_I", "a single number strictly between 0 and 1",
                               0, 1, lower_open = TRUE, upper_open = TRUE)
  sigma_tau <- .rp_check_vector(sigma_tau, "sigma_tau", positive = TRUE)
  kappa <- .rp_check_number(kappa, "kappa", "a single positive number", 0, lower_open = TRUE)
  if (is.null(J) == is.null(N)) {
    .rp_abort("Give the number of sites `J` or the total `N`, not both and not neither.")
  }
  log_nbar <- log(kappa) - 2 * log(sigma_tau) + stats::qlogis(target_I)
  check_sizes <- function(x) {
    if (any(!is.finite(x) | x <= 0)) {
      .rp_abort("The requested design sizes are outside the positive finite numeric range.",
                fix = "Revise the target, planned standard deviation or variance constant.")
    }
    x
  }
  if (!is.null(J)) {
    J <- .rp_check_whole(J, "J")
    nbar <- check_sizes(exp(log_nbar))
    return(data.frame(target_I = target_I, sigma_tau = sigma_tau, kappa = kappa, J = J,
                      n_per_site = nbar, N = check_sizes(exp(log(J) + log_nbar))))
  }
  N <- .rp_check_number(N, "N", "a single positive number", 0, lower_open = TRUE)
  J_max <- floor(exp(log(N) - log_nbar) * (1 + 1e-12))
  if (any(!is.finite(J_max))) {
    .rp_abort("The number of sites is outside the finite numeric range.",
              fix = "Revise the target, total or planned standard deviation.")
  }
  n_per_site <- rep(NA_real_, length(J_max))
  feasible <- J_max >= 1
  n_per_site[feasible] <- check_sizes(N / J_max[feasible])
  reached <- ifelse(J_max >= 1,
                     .rp_plogis(log(n_per_site) - log_nbar + stats::qlogis(target_I)), NA_real_)
  data.frame(target_I = target_I, sigma_tau = sigma_tau, kappa = kappa,
             J = ifelse(J_max >= 1, J_max, NA_real_), n_per_site = n_per_site,
             N = N, I = reached)
}

#' Ratio of the geometric to the arithmetic mean site size
#'
#' Gives the ratio of the geometric to the arithmetic mean of the site sizes,
#' the factor by which unequal sizes shrink the size that enters the
#' informativeness, as a function of their coefficient of variation.
#'
#' @param cv Coefficients of variation of the site sizes, numbers of at least
#'   0.
#' @param type `"lognormal"` (the default) for the exact ratio in a population
#'   of sizes (a lognormal population), \eqn{(1 + \mathrm{CV}^2)^{-1/2}}{(1 + CV^2)^(-1/2)};
#'   `"local"` for the approximation \eqn{\exp(-\mathrm{CV}^2/2)}{exp(-CV^2 / 2)}
#'   of a finite set of sizes close to their mean.
#'
#' @details
#' The two are different statements (Lee, 2026, Appendix E.2). The
#' lognormal ratio is an identity for a population; it does not give the
#' geometric mean of a finite set of sizes from their coefficient of
#' variation. The local approximation applies to a finite set, with an
#' explicit remainder, when the sizes are close to their mean, and needs the
#' coefficient of variation with divisor \eqn{J}, not \eqn{J - 1}.
#' [reliability_design()] reports the exact ratio of a planned allocation.
#' At a coefficient of variation of 0.75 the lognormal ratio is 0.80 and the
#' local approximation 0.755.
#'
#' @return A numeric vector of ratios, one per value of `cv`.
#'
#' @template ref-lee2026
#' @examples
#' imbalance_ratio(c(0.25, 0.5, 0.75))
#' imbalance_ratio(c(0.25, 0.5, 0.75), type = "local")
#' @export
imbalance_ratio <- function(cv, type = c("lognormal", "local")) {
  cv <- .rp_check_vector(cv, "cv")
  if (any(cv < 0)) .rp_abort("`cv` must hold numbers of at least 0.", info = .rp_you_passed(cv))
  type <- .rp_match(type, c("lognormal", "local"), "type")
  if (type == "lognormal") (1 + cv^2)^(-1 / 2) else exp(-cv^2 / 2)
}
