# Posterior moments and intervals of the site effects. Given sigma, with the
# mean of the effects (or the coefficients of the mean model) integrated out
# under a flat prior (Appendix A.5 of the paper),
#
#   E(tau_j | sigma, y)   = S_j y_j + (1 - S_j) x_j' beta_sigma,
#   Var(tau_j | sigma, y) = sigma^2 (1 - S_j) + (1 - S_j)^2 x_j' (X'WX)^-1 x_j,
#
# and the posterior moments average these over sigma. Without a mean model
# x_j' beta_sigma is the weighted mean and x_j' (X'WX)^-1 x_j = 1 / W.

# Conditional means and variances at a vector of sigma, for the sites in
# `cols`: two matrices with one row per sigma.
.rp_site_conditional <- function(kernel, sigma, cols = seq_len(kernel$J), centered = FALSE) {
  v <- kernel$v
  y <- kernel$y_work
  a <- sigma^2
  den <- outer(a, v[cols], "+")
  S <- a / den
  one_minus <- matrix(v[cols], nrow = length(sigma), ncol = length(cols), byrow = TRUE) / den
  if (is.null(kernel$X)) {
    w <- 1 / outer(a, v, "+")
    W <- rowSums(w)
    mu <- as.vector(w %*% y) / W
    fitted <- matrix(mu, nrow = length(sigma), ncol = length(cols))
    lev <- matrix(1 / W, nrow = length(sigma), ncol = length(cols))
  } else {
    fits <- lapply(sigma, function(s) .rp_mean_fit(kernel, s))
    fitted <- do.call(rbind, lapply(fits, function(f) f$fitted_centered[cols]))
    lev <- do.call(rbind, lapply(fits, function(f) f$leverage[cols]))
  }
  yy <- matrix(y[cols], nrow = length(sigma), ncol = length(cols), byrow = TRUE)
  vv <- matrix(v[cols], nrow = length(sigma), ncol = length(cols), byrow = TRUE)
  center <- fitted + S * (yy - fitted)
  if (!centered) center <- sweep(center, 2L, kernel$fitted_offset[cols], "+")
  list(mean = center, var = vv * S + one_minus^2 * lev, S = S)
}

# The nodes of a rule whose weight is above 1e-16 of the largest, for the
# quantiles of the normal mixtures of the site effects and coefficients:
# the others change a probability by less than about 1e-16 times their
# number.
.rp_site_nodes <- function(law) {
  keep <- law$rule$weight > 1e-16 * max(law$rule$weight)
  list(sigma = law$rule$sigma[keep], weight = law$rule$weight[keep] / sum(law$rule$weight[keep]))
}

# The nodes of a rule with positive weight, for posterior means. Every one
# is kept: a quantity such as I is tiny where the weight is largest when the
# posterior of sigma_tau is far below sqrt(GM(v)), and its mean then comes
# from nodes of small weight.
.rp_positive_nodes <- function(rule) {
  keep <- rule$weight > 0
  list(sigma = rule$sigma[keep], weight = rule$weight[keep])
}

# The site's term KL(Ber(I) || Ber(S_j)) of the dispersion at each sigma, for
# the sites in `cols`.
.rp_site_kl <- function(kernel, sigma, cols) {
  m <- 2 * log(sigma) - kernel$log_gm
  matrix(.rp_bern_kl(rep(m, times = length(cols)), rep(kernel$d[cols], each = length(sigma))),
         nrow = length(sigma))
}

# Posterior means and standard deviations of the site effects, and the
# posterior means of the reliabilities and of the sites' terms of H. The
# variance is the mean of the conditional variance and of the squared
# deviation of the conditional mean from the posterior mean, which keeps the
# digits of effects far from 0.
.rp_site_moments <- function(law, kernel, chunk_sites = 2000L) {
  J <- kernel$J
  if (law$engine == "adaptive") {
    second <- law$second()
    one <- function(j) {
      m1 <- law$expect(function(s) .rp_site_conditional(kernel, s, j, centered = TRUE)$mean[, 1L])
      m2 <- second$expect(function(s) {
        cc <- .rp_site_conditional(kernel, s, j, centered = TRUE)
        cc$var[, 1L] + (cc$mean[, 1L] - m1)^2
      })
      r <- law$expect(function(s) s^2 / (s^2 + kernel$v[[j]]))
      kl <- law$expect(function(s) .rp_site_kl(kernel, s, j)[, 1L])
      c(m1 + kernel$fitted_offset[[j]], m2, r, kl)
    }
    out <- vapply(seq_len(J), one, numeric(4))
    return(data.frame(effect = out[1L, ], effect_sd = sqrt(pmax(0, out[2L, ])),
                      reliability = out[3L, ], kl_term = out[4L, ]))
  }
  nodes <- .rp_positive_nodes(law$rule)
  effect_mean <- effect_sd <- reliability <- kl_term <- numeric(J)
  for (start in seq.int(1L, J, by = chunk_sites)) {
    cols <- start:min(J, start + chunk_sites - 1L)
    cc <- .rp_site_conditional(kernel, nodes$sigma, cols, centered = TRUE)
    m1 <- colSums(nodes$weight * cc$mean)
    m2 <- colSums(nodes$weight * (cc$var + sweep(cc$mean, 2L, m1)^2))
    effect_mean[cols] <- m1 + kernel$fitted_offset[cols]
    effect_sd[cols] <- sqrt(pmax(0, m2))
    reliability[cols] <- colSums(nodes$weight * cc$S)
    kl_term[cols] <- colSums(nodes$weight * .rp_site_kl(kernel, nodes$sigma, cols))
  }
  data.frame(effect = effect_mean, effect_sd = effect_sd, reliability = reliability,
             kl_term = kl_term)
}

# Central posterior intervals of the site effects: quantiles of the normal
# mixture over sigma, found for all sites at once by a safeguarded Newton's
# method on the mixture's CDF, whose derivative is the mixture's density.
# The adaptive engine uses the nodes of a node rule here.
.rp_site_intervals <- function(law, kernel, level, chunk_sites = 2000L) {
  nodes <- .rp_site_nodes(law$second())
  J <- kernel$J
  probs <- c((1 - level) / 2, (1 + level) / 2)
  out <- matrix(NA_real_, J, 2L)
  for (start in seq.int(1L, J, by = chunk_sites)) {
    cols <- start:min(J, start + chunk_sites - 1L)
    cc <- .rp_site_conditional(kernel, nodes$sigma, cols, centered = TRUE)
    for (k in 1:2) {
      out[cols, k] <- .rp_mixture_quantile(nodes$weight, cc$mean, sqrt(cc$var), probs[[k]]) +
        kernel$fitted_offset[cols]
    }
  }
  out
}

# The p-quantile of each column's normal mixture with weights w (one per
# row), means `center` and standard deviations `sd`. The iterations start at
# the normal quantile with the mixture's mean and standard deviation, and
# stop when a step is below 1e-12 of that standard deviation (or a few units
# in the last place of the quantile). Zero-SD components are point masses:
# their CDF jumps are checked first, and the remaining quantiles reduce to
# the continuous mixture with its probability adjusted for lower atoms.
.rp_mixture_quantile <- function(w, center, sd, p) {
  atomic_cols <- which(colSums(sd == 0) > 0L)
  if (length(atomic_cols)) {
    out <- numeric(ncol(center))
    continuous_cols <- setdiff(seq_len(ncol(center)), atomic_cols)
    if (length(continuous_cols)) {
      out[continuous_cols] <- .rp_mixture_quantile(
        w, center[, continuous_cols, drop = FALSE], sd[, continuous_cols, drop = FALSE], p
      )
    }
    for (j in atomic_cols) {
      atom <- sd[, j] == 0
      values <- center[, j]
      # The sums of weights carry rounding, so p = k/n can be missed by a
      # unit in the last place; the comparisons allow for it.
      slack <- 64 * .Machine$double.eps
      if (all(atom)) {
        ordered <- order(values)
        out[[j]] <- values[ordered][which(cumsum(w[ordered]) >= p - slack)[[1L]]]
        next
      }
      continuous_cdf <- function(x) sum(w[!atom] * stats::pnorm((x - values[!atom]) / sd[!atom, j]))
      atoms <- sort(unique(values[atom]))
      left <- vapply(atoms, function(x) continuous_cdf(x) + sum(w[atom & values < x]), 0)
      right <- vapply(atoms, function(x) continuous_cdf(x) + sum(w[atom & values <= x]), 0)
      hit <- which(left <= p + slack & right >= p - slack)
      if (length(hit)) {
        out[[j]] <- atoms[[hit[[1L]]]]
        next
      }
      # Between adjacent atoms their CDF contribution is constant. Remove
      # that mass before solving the continuous mixture, so a distant atom
      # cannot set the quantile's tolerance or search interval.
      below <- atoms[right < p]
      mass_below <- if (length(below)) sum(w[atom & values <= max(below)]) else 0
      continuous_weight <- sum(w[!atom])
      anchor <- values[which(!atom)[[1L]]]
      out[[j]] <- .rp_mixture_quantile(
        w[!atom] / continuous_weight,
        matrix(values[!atom] - anchor, ncol = 1L), sd[!atom, j, drop = FALSE],
        (p - mass_below) / continuous_weight
      )[[1L]] + anchor
    }
    return(out)
  }
  n <- nrow(center)
  lo <- apply(center - 12 * sd, 2L, min)
  hi <- apply(center + 12 * sd, 2L, max)
  mix_mean <- colSums(w * center)
  mix_sd <- sqrt(colSums(w * (sd^2 + sweep(center, 2L, mix_mean)^2)))
  tol <- pmax(1e-12 * mix_sd, 4 * .Machine$double.eps * abs(mix_mean))
  x <- mix_mean + stats::qnorm(p) * mix_sd
  x <- pmin(pmax(x, lo), hi)
  for (iter in seq_len(100L)) {
    z <- (matrix(x, n, length(x), byrow = TRUE) - center) / sd
    F <- colSums(w * stats::pnorm(z)) - p
    dens <- colSums(w * stats::dnorm(z) / sd)
    below <- F < 0
    lo[below] <- x[below]
    hi[!below] <- x[!below]
    nx <- x - F / dens
    bad <- !is.finite(nx) | nx <= lo | nx >= hi
    nx[bad] <- (lo[bad] + hi[bad]) / 2
    done <- abs(nx - x) <= pmax(tol, 4 * .Machine$double.eps * abs(nx))
    x <- nx
    if (all(done)) break
  }
  x
}

# Posterior means of the sites' shares S_j / sum_k S_k, over weighted values
# of sigma, in blocks of sites. Each S is divided by the largest of its row
# on the log scale (the largest is that of the largest d), so that nothing
# underflows; a first pass adds up the rows, a second averages the shares.
# At sigma = 0 the normalized inverse variances give the limiting shares.
.rp_site_shares <- function(kernel, sigma, weight, chunk_sites = 2000L) {
  J <- kernel$J
  d <- kernel$d
  m <- 2 * log(sigma) - kernel$log_gm
  top <- stats::plogis(m + max(d), log.p = TRUE)
  blocks <- split(seq_len(J), ceiling(seq_len(J) / chunk_sites))
  scaled <- function(cols) {
    out <- matrix(exp(d[cols] - max(d)), nrow = length(m), ncol = length(cols), byrow = TRUE)
    positive <- sigma > 0
    if (any(positive)) {
      out[positive, ] <- exp(stats::plogis(outer(m[positive], d[cols], "+"), log.p = TRUE) - top[positive])
    }
    out
  }
  total <- numeric(length(m))
  for (cols in blocks) total <- total + rowSums(scaled(cols))
  share <- numeric(J)
  for (cols in blocks) share[cols] <- colSums(weight * scaled(cols) / total)
  share
}
