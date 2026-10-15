# Random reliability profiles for the property tests. Each profile is given by
# sampling variances and a between-site SD, so that the reliabilities and
# their logits can be formed in several ways and compared.
random_profile <- function(J = sample(2:60, 1), log_ratio = runif(1, 0, 8),
                           sigma_tau = exp(runif(1, -4, 2))) {
  v <- exp(runif(J, 0, log_ratio) + runif(1, -6, 2))
  list(v = v, sigma_tau = sigma_tau, J = J,
       S = sigma_tau^2 / (sigma_tau^2 + v))
}

# Profiles whose logits are given directly: m is the logit of I and d the
# centered deviations.
profile_from_logits <- function(m, d) {
  d <- d - mean(d)
  stats::plogis(m + d)
}

# The paper's skewed example (Appendix C): J = 25, I = 0.3, V = 2.
paper_counterexample_d <- c(sqrt(48), rep(-sqrt(1 / 12), 24))

# Naive formulas, accurate for moderate inputs, to compare with.
naive_I <- function(S) stats::plogis(mean(stats::qlogis(S)))
naive_kl <- function(p, q) p * log(p / q) + (1 - p) * log((1 - p) / (1 - q))
naive_H <- function(S) mean(naive_kl(naive_I(S), S))
naive_jeff <- function(S) sum(S)^2 / sum(S^2)
