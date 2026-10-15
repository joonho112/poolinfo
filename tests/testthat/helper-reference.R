# A brute-force posterior of sigma_tau for the tests: 20-point
# Gauss-Legendre rules on many equal panels, in sigma on [0, lo] and in
# log(sigma) on [lo, hi]. It shares only the kernel with the package; the
# rule is independent of the package's rules.
brute_rule <- function(y, v, prior = sd_prior(), X = NULL, lo, hi, panels = 4000L, zero_panels = 100L) {
  g <- .rp_gauss_legendre(20L)
  rule_on <- function(edges) {
    half <- diff(edges) / 2
    mid <- (edges[-1L] + edges[-length(edges)]) / 2
    list(x = as.vector(outer(g$x, half) + rep(mid, each = 20L)), w = as.vector(outer(g$w, half)))
  }
  zero <- rule_on(seq(0, lo, length.out = zero_panels + 1L))
  logs <- rule_on(seq(log(lo), log(hi), length.out = panels + 1L))
  sigma <- c(zero$x, exp(logs$x))
  k <- .rp_kernel(y, v, X, prior)
  lk <- .rp_log_kernel(k, sigma)
  w <- exp(lk - max(lk)) * c(zero$w, logs$w * exp(logs$x))
  list(sigma = sigma, weight = w / sum(w), kernel = k, log_z = max(lk) + log(sum(w)))
}

# The CDF of sigma_tau at q, by integrate() on pieces evenly spaced in
# log(sigma) below q, with the normalizing constant of a brute-force rule.
brute_cdf <- function(rule, q) {
  f <- function(s) exp(.rp_log_kernel(rule$kernel, s) - rule$log_z)
  edges <- c(0, exp(seq(log(q) - 40, log(q), length.out = 400)))
  sum(vapply(seq_len(length(edges) - 1L), function(i) {
    stats::integrate(f, edges[[i]], edges[[i + 1L]], rel.tol = 1e-12, abs.tol = 0,
                     stop.on.error = FALSE)$value
  }, numeric(1)))
}

brute_mean <- function(rule, fun) sum(rule$weight * fun(rule$sigma))

summary_value <- function(fit, quantity, column = "mean") {
  fit$summary[[column]][fit$summary$quantity == quantity]
}
