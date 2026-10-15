# Tolerances of the tests, by kind of quantity. The posterior tolerances are
# those of the table in vignette("implementation"); the algebra tolerance is
# for identities that hold exactly in exact arithmetic.

# Identities of the profile algebra, computed in closed form.
tol_identity <- 1e-12

# Bounded summaries: I, mean reliability, H, S_j and probabilities (absolute).
tol_bounded <- 1e-9

# The concentration count J_eff (relative).
tol_jeff_rel <- 1e-9

# Quantities on the scale of the effects: sigma_tau, the effect mean, site
# posterior means and standard deviations, and quantiles of sigma_tau. The
# scale is the square root of the geometric mean of the sampling variances;
# a value far beyond it, as an upper quantile of sigma_tau with 3 sites, is
# compared relative to its own size.
tol_effect_scale <- function(v, value = 0) 1e-9 * max(sqrt(exp(mean(log(v)))), abs(value))

# E(sigma_tau^2) and its quantiles (relative); the tail of sigma_tau^2 is
# heavy. I at E(sigma_tau^2) is bounded and is compared absolutely, with
# tol_bounded against stored values and with this number between engines.
tol_sigma2_rel <- 1e-8

# Relative error of the normalizing constant, and bounds on the omitted tail.
tol_normalization <- 1e-10
tol_tail_mass <- 1e-12
tol_tail_moment <- function(v, q) 1e-10 * exp(mean(log(v)))^(q / 2)

# Turning points of H found at two resolutions (logit scale).
tol_turning <- 1e-7

# Relative difference between two values.
rel_diff <- function(x, y) abs(x - y) / pmax(abs(y), .Machine$double.xmin)
