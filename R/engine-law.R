# The posterior law of sigma_tau in one form for both engines: a list with
# the functions expect(), cdf(), quantile() and density(). It is rebuilt
# from the numbers stored in a fit whenever a summary needs it, so that the
# fit itself holds no functions.
#
# second() is the law for the second moments about a mean, from which the
# standard deviations come. For the node rule it is the law itself. The
# adaptive engine uses a node rule for them: the squared deviation of a
# quantity that hardly varies is mostly rounding, which integrate() cannot
# bring to its relative tolerance, while a sum over nodes is not affected.

.rp_law <- function(kernel, engine = c("nodes", "adaptive"), rel_tol = 1e-10, abs_tol = 1e-12,
                    rule = NULL) {
  engine <- match.arg(engine)
  if (engine == "adaptive") {
    law <- .rp_adaptive_law(kernel, rel_tol = rel_tol, abs_tol = abs_tol)
    law$expect_matrix <- function(fun, k) {
      vapply(seq_len(k), function(j) law$expect(function(s) as.matrix(fun(s))[, j]), numeric(1))
    }
    nodes <- NULL
    law$second <- function() {
      if (is.null(nodes)) nodes <<- .rp_law(kernel, "nodes")
      nodes
    }
    return(law)
  }
  if (is.null(rule)) {
    rule <- .rp_build_node_rule(kernel)
  }
  # The points of weight 0 (where the kernel underflows, or where sigma^2
  # overflows) are left out, so that a value that is not finite there does
  # not spoil a sum.
  nodes <- .rp_positive_nodes(rule)
  expect_matrix <- function(fun, k) colSums(nodes$weight * as.matrix(fun(nodes$sigma)))
  law <- list(
    engine = "nodes", kernel = kernel, rule = rule,
    expect = function(fun) sum(nodes$weight * fun(nodes$sigma)),
    expect_matrix = expect_matrix,
    cdf = function(s, lower = TRUE) .rp_rule_cdf(rule, s, lower),
    quantile = function(p) .rp_rule_quantile(rule, p),
    density = function(s) exp(.rp_log_kernel(kernel, s) - rule$log_Z),
    log_Z = rule$log_Z, upper = rule$upper
  )
  law$second <- function() law
  law
}

# The diagnostics at each sigma of a vector: a matrix with the columns I,
# mean, harmonic, H, J_eff, sigma and sigma2. At sigma = 0 they take their
# limits.
.rp_functionals_at <- function(kernel, sigma) {
  d <- kernel$d
  m <- 2 * log(sigma) - kernel$log_gm
  I <- .rp_plogis(m)
  mean_s <- .rp_mean_at(m, d)
  cbind(
    I = I,
    mean = mean_s,
    harmonic = pmin(.rp_harmonic_at(m, d), I, mean_s),
    H = .rp_H_at(m, d),
    J_eff = .rp_jeff_at(m, d),
    sigma = sigma,
    sigma2 = sigma^2
  )
}

.rp_functional_names <- c("I", "mean", "harmonic", "H", "J_eff", "sigma", "sigma2")

# The functionals that are increasing in sigma, whose quantiles are those of
# sigma mapped through them.
.rp_monotone_functionals <- c("I", "mean", "harmonic", "J_eff", "sigma", "sigma2")

# The mean model at each sigma: the weighted least-squares coefficients and
# the diagonal of (X'WX)^-1, as a matrix with 2p columns.
.rp_coef_at <- function(kernel, sigma, centered = FALSE) {
  p <- kernel$p
  t(vapply(sigma, function(s) {
    if (is.null(kernel$X)) {
      w <- 1 / (s^2 + kernel$v)
      W <- sum(w)
      beta <- sum(w * kernel$y_work) / W
      return(c(beta + if (centered) 0 else kernel$coef_offset, 1 / W))
    }
    w <- 1 / (s^2 + kernel$v)
    sw <- sqrt(w)
    qx <- qr(kernel$X * sw)
    coef <- qr.coef(qx, kernel$y_work * sw)
    if (!centered) coef <- coef + kernel$coef_offset
    rinv <- backsolve(qr.R(qx), diag(p))
    # qr.coef() returns the coefficients in the order of the columns of X;
    # R is in the order of the pivoting.
    c(coef, rowSums(rinv^2)[order(qx$pivot)])
  }, numeric(2L * p)))
}
