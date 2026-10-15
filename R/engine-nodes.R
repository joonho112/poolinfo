# The node rule: one set of points and weights on sigma_tau that every
# posterior expectation reuses.
#
# The rule has two parts. Above sigma_lo it is made of Gauss-Legendre panels
# on t = log(sigma), with the Jacobian d sigma = sigma dt; the panels are
# split where a 16-point and a 32-point rule disagree. Below sigma_lo it is
# one Gauss-Legendre panel in sigma itself: under the flat prior the density
# is positive at 0, so the mass below sigma_lo is about g(0) sigma_lo, which
# a window on the log scale alone would miss. sigma_lo is exp(-8) times the
# smallest of the smallest standard error, the scale of the prior and the
# upper limit of a uniform prior, so that the kernel and the prior hardly
# change below it. The upper end is where the bounds on the omitted tail of
# Appendix D.3 of the paper fall below 1e-12 for the mass and 1e-10 scale^q
# for the moments that exist.

# Gauss-Legendre points and weights on [-1, 1] (Golub-Welsch).
.rp_gauss_legendre <- function(n) {
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i * i - 1)
  A <- matrix(0, n, n)
  A[cbind(i, i + 1L)] <- b
  A[cbind(i + 1L, i)] <- b
  e <- eigen(A, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1L, o]^2)
}

.rp_gl_cache <- new.env(parent = emptyenv())

.rp_gl <- function(n) {
  key <- as.character(n)
  if (is.null(.rp_gl_cache[[key]])) {
    assign(key, .rp_gauss_legendre(n), envir = .rp_gl_cache)
  }
  .rp_gl_cache[[key]]
}

# Points and weights of n-point rules on the panels [lo_k, hi_k] of t =
# log(sigma), in sigma, with the Jacobian in the weights. Panel k's points
# are in the k-th block of n.
.rp_panel_points <- function(lo, hi, n) {
  g <- .rp_gl(n)
  half <- (hi - lo) / 2
  mid <- (hi + lo) / 2
  t <- as.vector(t(outer(half, g$x) + mid))
  w <- as.vector(t(outer(half, g$w)))
  list(sigma = exp(t), weight = w * exp(t), panel = rep(seq_along(lo), each = n))
}

# Points and weights of an n-point rule on [0, upper] in sigma.
.rp_zero_points <- function(upper, n) {
  g <- .rp_gl(n)
  list(sigma = upper / 2 * (g$x + 1), weight = upper / 2 * g$w)
}

# Integrals of exp(log k - shift) times 1 and times (sigma / scale)^2 over
# each panel, with n-point rules.
.rp_panel_integrals <- function(kernel, lo, hi, n, shift, with_moment) {
  pts <- .rp_panel_points(lo, hi, n)
  dens <- exp(.rp_log_kernel(kernel, pts$sigma) - shift) * pts$weight
  mass <- rowsum(dens, pts$panel, reorder = FALSE)[, 1L]
  moment <- if (with_moment) rowsum(dens * (pts$sigma / kernel$scale)^2, pts$panel, reorder = FALSE)[, 1L]
  list(mass = mass, moment = moment, max_log = max(.rp_log_kernel(kernel, pts$sigma)))
}

# The lower end sigma_lo of the panels on the log scale. Below it, sigma^2
# is below exp(-16) of every sampling variance and the prior changes by less
# than about exp(-16): one panel in sigma resolves the density there.
.rp_sigma_floor <- function(kernel) {
  prior <- kernel$prior
  sizes <- c(sqrt(min(kernel$v)),
             if (prior$type %in% c("half_normal", "half_t", "half_cauchy")) prior$scale,
             if (prior$type == "uniform") prior$upper)
  exp(-8) * min(sizes)
}

# Builds the rule. `refine` halves every panel once more (for the check of
# the rule against itself) and `lower_factor` moves sigma_lo down (for the
# check of the mass below it).
.rp_build_node_rule <- function(kernel, n = 32L, refine = 0L, lower_factor = 1,
                                max_panels = 20000L) {
  sigma_lo <- .rp_sigma_floor(kernel) * lower_factor
  t_lo <- log(sigma_lo)
  upper <- .rp_upper_end(kernel)
  t_hi <- log(upper)
  edges <- seq(t_lo, t_hi, length.out = max(16L, ceiling((t_hi - t_lo) / 0.5)) + 1L)
  with_moment <- .rp_moment_exists(2, kernel$J, kernel$p, kernel$prior)
  noise <- 64 * kernel$J * .Machine$double.eps
  # A first pass on a coarse grid gives the shift that keeps exp() in range.
  shift <- max(.rp_log_kernel(kernel, c(0, exp(seq(t_lo, t_hi, length.out = 513)))))
  lo <- edges[-length(edges)]
  hi <- edges[-1L]
  done <- rep(FALSE, length(lo))
  mass <- moment <- rep(NA_real_, length(lo))
  for (iter in seq_len(60L)) {
    todo <- which(!done)
    if (length(todo) == 0L) break
    coarse <- .rp_panel_integrals(kernel, lo[todo], hi[todo], n %/% 2L, shift, with_moment)
    fine <- .rp_panel_integrals(kernel, lo[todo], hi[todo], n, shift, with_moment)
    if (fine$max_log > shift) {
      # A higher mode than the first pass saw: rescale to it and redo.
      ratio <- exp(shift - fine$max_log)
      mass <- mass * ratio
      moment <- moment * ratio
      shift <- fine$max_log
      next
    }
    mass[todo] <- fine$mass
    if (with_moment) moment[todo] <- fine$moment
    # A panel is resolved when its two rules agree to 1e-14 of the total, or
    # to the rounding noise of the kernel, a sum of J terms, relative to the
    # panel's own integral: below that, splitting cannot help.
    bad <- abs(fine$mass - coarse$mass) > pmax(1e-14 * sum(mass, na.rm = TRUE), noise * fine$mass)
    if (with_moment) {
      bad <- bad | abs(fine$moment - coarse$moment) >
        pmax(1e-14 * sum(moment, na.rm = TRUE), noise * fine$moment)
    }
    done[todo[!bad]] <- TRUE
    if (!any(bad)) break
    split <- todo[bad]
    mid <- (lo[split] + hi[split]) / 2
    lo <- c(lo[-split], lo[split], mid)
    hi <- c(hi[-split], mid, hi[split])
    done <- c(done[-split], rep(FALSE, 2L * length(split)))
    mass <- c(mass[-split], rep(NA_real_, 2L * length(split)))
    moment <- c(moment[-split], rep(NA_real_, 2L * length(split)))
    o <- order(lo)
    lo <- lo[o]
    hi <- hi[o]
    done <- done[o]
    mass <- mass[o]
    moment <- moment[o]
    if (length(lo) > max_panels) {
      .rp_abort(
        "The node rule needs more panels than its limit.",
        info = sprintf("It has %d panels after %d rounds of splitting.", length(lo), iter),
        fix = paste("This should not happen; please report the estimates and standard errors at",
                    "https://github.com/joonho112/poolinfo/issues."),
        class = "rp_numeric_error", subclass = "rp_rule_error"
      )
    }
  }
  if (!all(done)) {
    .rp_abort(
      "The panels of the node rule did not settle.",
      fix = paste("This should not happen; please report the estimates and standard errors at",
                  "https://github.com/joonho112/poolinfo/issues."),
      class = "rp_numeric_error", subclass = "rp_rule_error"
    )
  }
  edges <- c(lo, hi[length(hi)])
  for (r in seq_len(refine)) {
    edges <- .rp_halve_edges(edges)
  }
  .rp_rule_from_edges(kernel, edges, sigma_lo, upper, n)
}

.rp_halve_edges <- function(edges) {
  sort(c(edges, (edges[-1L] + edges[-length(edges)]) / 2))
}

# The rule on given panels: the points, the normalized weights, the mass of
# each panel and the Legendre coefficients of the density in each panel.
.rp_rule_from_edges <- function(kernel, edges, sigma_lo, upper, n = 32L) {
  pts <- .rp_panel_points(edges[-length(edges)], edges[-1L], n)
  zero <- .rp_zero_points(sigma_lo, n)
  sigma <- c(zero$sigma, pts$sigma)
  log_k <- .rp_log_kernel(kernel, sigma)
  shift <- max(log_k)
  raw <- exp(log_k - shift) * c(zero$weight, pts$weight)
  Z <- sum(raw)
  panel <- c(rep(0L, n), pts$panel)
  # The density at the points of each panel, in the panel's own variable
  # (sigma for the lowest panel, t = log(sigma) for the others), as the
  # coefficients of its Legendre series: the CDF inside a panel is the
  # integral of that series.
  dens <- exp(log_k - shift) / Z
  dens_var <- c(dens[panel == 0L], dens[panel > 0L] * sigma[panel > 0L])
  coef <- .rp_legendre_coef(matrix(dens_var, nrow = n), .rp_gl(n))
  list(
    engine = "nodes", sigma = sigma, weight = raw / Z,
    panel = panel, edges = edges, sigma_lo = sigma_lo, upper = upper,
    panel_mass = c(sum(raw[panel == 0L]), rowsum(raw[panel > 0L], panel[panel > 0L], reorder = FALSE)[, 1L]) / Z,
    coef = coef, log_Z = shift + log(Z), n = n
  )
}

# Legendre coefficients c_0, ..., c_{n-1} of the polynomial that takes the
# values f (one column per panel) at the n Gauss-Legendre points of [-1, 1].
.rp_legendre_coef <- function(f, gl) {
  P <- .rp_legendre_values(gl$x, nrow(f))
  k <- seq_len(nrow(f)) - 1L
  ((2 * k + 1) / 2) * crossprod(P * gl$w, f)
}

# The Legendre polynomials P_0, ..., P_{n-1} at the points x: a matrix with
# one row per point.
.rp_legendre_values <- function(x, n) {
  P <- matrix(0, length(x), n + 1L)
  P[, 1L] <- 1
  if (n >= 1L) P[, 2L] <- x
  if (n >= 2L) {
    for (k in 1:(n - 1L)) {
      P[, k + 2L] <- ((2 * k + 1) * x * P[, k + 1L] - k * P[, k]) / (k + 1)
    }
  }
  P[, seq_len(n), drop = FALSE]
}

# The integral from -1 to x of the Legendre series with coefficients `coef`
# (one column per panel; `which` gives the panel of each x). With
# P_{k+1}' - P_{k-1}' = (2k + 1) P_k, the integral of P_k from -1 is
# (P_{k+1}(x) - P_{k-1}(x)) / (2k + 1), and x + 1 for P_0.
.rp_legendre_integral <- function(coef, which, x) {
  n <- nrow(coef)
  P <- .rp_legendre_values(x, n + 1L)
  ints <- matrix(0, length(x), n)
  ints[, 1L] <- x + 1
  for (k in seq_len(n - 1L)) {
    ints[, k + 1L] <- (P[, k + 2L] - P[, k]) / (2 * k + 1)
  }
  rowSums(ints * t(coef[, which, drop = FALSE]))
}

# The upper end of sigma: doubled from a starting point until the bounds on
# the omitted tail meet their limits, with the normalizing constant over
# [0, upper] (smaller than the whole, so the bound is conservative).
.rp_upper_end <- function(kernel) {
  if (kernel$prior$type == "uniform") {
    return(kernel$prior$upper)
  }
  B <- .rp_upper_guess(kernel)
  orders <- .rp_tail_orders(kernel)
  log_z <- .rp_log_mass_below(kernel, B)
  for (i in seq_len(400L)) {
    bounds <- vapply(orders, function(q) .rp_log_tail_bound(kernel, B, q), numeric(1)) - log_z
    limits <- vapply(orders, function(q) .rp_log_tail_limit(kernel, q), numeric(1))
    if (all(bounds <= limits)) {
      return(B)
    }
    B <- 2 * B
    if (!is.finite(B)) break
  }
  failed <- which(!is.finite(bounds) | bounds > limits)
  q <- if (length(failed)) max(orders[failed]) else max(orders)
  exponent <- kernel$J - kernel$p + kernel$prior$tail - q - 1
  .rp_abort(
    sprintf("The posterior of sigma_tau has too heavy a tail to bound its moment of order %g.", q),
    info = paste0(if (is.finite(exponent)) sprintf("The omitted tail falls off only like B^(-%g). ", exponent),
                  "The moment exists, but no upper end of the integral is far enough."),
    fix = "Use a prior with a lighter tail, such as a half-t prior with more degrees of freedom.",
    class = "rp_numeric_error", subclass = "rp_rule_error"
  )
}

# log of the integral of k over [0, B], from a coarse rule; used only to
# make the tail bounds relative.
.rp_log_mass_below <- function(kernel, B) {
  lo <- log(.rp_sigma_floor(kernel))
  edges <- seq(lo, log(B), length.out = max(65L, ceiling((log(B) - lo) / 0.5) + 1L))
  pts <- .rp_panel_points(edges[-length(edges)], edges[-1L], 16L)
  zero <- .rp_zero_points(exp(lo), 16L)
  log_k <- .rp_log_kernel(kernel, c(zero$sigma, pts$sigma))
  shift <- max(log_k)
  shift + log(sum(exp(log_k - shift) * c(zero$weight, pts$weight)))
}

# The CDF of sigma at s under a node rule: the mass of the panels below s
# and the integral of the Legendre series of the panel that holds s. With
# `lower = FALSE`, the probability above s, summed from the panels above, so
# that a small upper tail keeps its digits.
.rp_rule_cdf <- function(rule, s, lower = TRUE) {
  out <- numeric(length(s))
  cum <- cumsum(rule$panel_mass)
  above <- rev(cumsum(rev(rule$panel_mass)))
  for (i in seq_along(s)) {
    x <- s[[i]]
    if (is.na(x)) {
      out[[i]] <- NA_real_
      next
    }
    if (x <= 0) {
      out[[i]] <- if (lower) 0 else 1
      next
    }
    if (x >= rule$upper) {
      out[[i]] <- if (lower) 1 else 0
      next
    }
    if (x <= rule$sigma_lo) {
      k <- 1L
      u <- 2 * x / rule$sigma_lo - 1
      part <- rule$sigma_lo / 2 * .rp_legendre_integral(rule$coef, 1L, u)
    } else {
      t <- log(x)
      j <- min(findInterval(t, rule$edges, rightmost.closed = TRUE), length(rule$edges) - 1L)
      k <- j + 1L
      a <- rule$edges[[j]]
      b <- rule$edges[[j + 1L]]
      part <- (b - a) / 2 * .rp_legendre_integral(rule$coef, k, (2 * t - a - b) / (b - a))
    }
    part <- min(max(part, 0), rule$panel_mass[[k]])
    out[[i]] <- if (lower) {
      min(1, max(0, (if (k > 1L) cum[[k - 1L]] else 0) + part))
    } else {
      min(1, max(0, (if (k < length(above)) above[[k + 1L]] else 0) + rule$panel_mass[[k]] - part))
    }
  }
  out
}

# Quantiles of sigma under a node rule: the panel from the cumulative masses,
# then a root of the CDF inside it.
.rp_rule_quantile <- function(rule, p) {
  cum <- cumsum(rule$panel_mass)
  vapply(p, function(pp) {
    if (pp <= 0) return(0)
    if (pp >= 1) return(rule$upper)
    k <- which(cum >= pp)[1L]
    if (is.na(k)) return(rule$upper)
    in_t <- function(t) .rp_rule_cdf(rule, exp(t)) - pp
    if (k == 1L) {
      # In the lowest panel the CDF is about g(0) sigma: the root is found on
      # the log scale, with a lower end below it.
      lo <- rule$sigma_lo * min(1, pp / cum[[1L]]) / 16
      for (i in seq_len(60L)) {
        if (.rp_rule_cdf(rule, lo) < pp) break
        lo <- lo / 16
      }
      return(exp(stats::uniroot(in_t, log(c(lo, rule$sigma_lo)), tol = 1e-14)$root))
    }
    exp(stats::uniroot(in_t, rule$edges[c(k - 1L, k)], tol = 1e-14)$root)
  }, numeric(1))
}
