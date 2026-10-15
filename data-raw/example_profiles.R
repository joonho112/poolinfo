# Five profiles of 40 sites with the informativeness I = 0.3 at
# sigma_tau = 0.2, built by shape_profile(): the profiles of Figure 2 of
# Lee (2026), equal variances and the even and concentrated-precise shapes at
# V = 1.5, and the bimodal and concentrated-noisy shapes at the same V.
#
#   Rscript data-raw/example_profiles.R   (from the package root)
pkgload::load_all(".", quiet = TRUE)
shapes <- c(homogeneous = "even", even = "even", bimodal = "bimodal",
            concentrated_precise = "concentrated_precise", concentrated_noisy = "concentrated_noisy")
example_profiles <- do.call(rbind, lapply(names(shapes), function(name) {
  V <- if (name == "homogeneous") 0 else 1.5
  v <- shape_profile(J = 40, I = 0.3, V = V, shape = shapes[[name]], sigma_tau = 0.2)
  data.frame(profile = name, site = seq_len(40), se2 = v, se = sqrt(v), sigma_tau = 0.2,
             reliability = 0.04 / (0.04 + v))
}))
example_profiles$profile <- factor(example_profiles$profile, levels = names(shapes))
save(example_profiles, file = "data/example_profiles.rda", compress = "xz", version = 2)
