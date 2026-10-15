# Checks of the numerical integration, kept in the fit and shown by
# check_numerics(). The node rule is compared with a rule whose panels are
# all halved (which shows whether the panels resolve the density) and with a
# rule whose lower end sigma_lo is halved (which shows the mass below the
# window on the log scale, which refining the panels cannot see). The
# omitted upper tail is bounded analytically. Both scans for the turning
# points of H must agree.
#
# The comparisons of two rules are relative. Their limit is 1e-10, or the
# rounding noise of log k, a sum of J terms, when that is larger: below it
# two rules cannot agree.

.rp_numeric_checks <- function(kernel, law, turning) {
  rows <- list()
  add <- function(check, value, limit) {
    rows[[length(rows) + 1L]] <<- data.frame(check = check, value = value, limit = limit)
  }
  limit <- max(1e-10, 64 * kernel$J * .Machine$double.eps)
  key <- function(rule) {
    nodes <- .rp_positive_nodes(rule)
    m <- 2 * log(nodes$sigma) - kernel$log_gm
    c(sum(nodes$weight * .rp_plogis(m)),
      if (!turning$homogeneous) sum(nodes$weight * .rp_H_at(m, kernel$d)))
  }
  relative <- function(a, b) max(abs(a - b) / pmax(abs(b), .Machine$double.xmin))
  if (law$engine == "nodes") {
    rule <- law$rule
    base <- key(rule)
    refined <- .rp_rule_from_edges(kernel, .rp_halve_edges(rule$edges), rule$sigma_lo, rule$upper, rule$n)
    lower_lo <- rule$sigma_lo / 2
    lower <- .rp_rule_from_edges(kernel, c(log(lower_lo), rule$edges), lower_lo, rule$upper, rule$n)
    add("normalizing constant, rule against halved panels (relative)",
        abs(expm1(refined$log_Z - law$rule$log_Z)), limit)
    add("E(I) and E(H), rule against halved panels (relative)", relative(key(refined), base), limit)
    add("normalizing constant, lower end moved down (relative)",
        abs(expm1(lower$log_Z - law$rule$log_Z)), limit)
    add("E(I) and E(H), lower end moved down (relative)", relative(key(lower), base), limit)
  } else {
    add("normalizing constant, error estimate of integrate() (relative)",
        law$checks[["normalization_error"]], limit)
  }
  # The bounds are relative to sqrt(GM(v))^q and compared on the log scale,
  # which keeps them in range whatever the units.
  for (q in .rp_tail_orders(kernel)) {
    log_bound <- .rp_log_tail_bound(kernel, law$upper, q) - law$log_Z - q * log(kernel$scale)
    label <- if (q == 0) "mass beyond the upper end (bound)" else
      sprintf("E(sigma_tau^%d) beyond the upper end, over sqrt(GM(v))^%d (bound)", q, q)
    add(label, exp(log_bound), exp(.rp_log_tail_limit(kernel, q) - q * log(kernel$scale)))
  }
  add("turning points of H at two resolutions (logit scale)",
      if (turning$homogeneous) 0 else turning$difference, 1e-7)
  out <- do.call(rbind, rows)
  out$status <- ifelse(is.finite(out$value) & out$value <= out$limit, "ok", "outside")
  out
}
