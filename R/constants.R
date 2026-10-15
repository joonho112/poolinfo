# Numerical constants of the package, kept in one place. The tolerances of
# the posterior computations are the defaults of posterior_control().

# A profile whose logit variance V is below this bound is treated as one
# whose reliabilities are all equal: its dispersion H is 0 and has no
# turning point. At this V, H is at most V / 8 plus a remainder, about 1e-15.
.RP_V_HOMOGENEOUS <- 1e-14
