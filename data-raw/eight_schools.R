# The eight schools of Rubin (1981): estimated effects of coaching programs
# on the SAT-V in eight high schools, with their standard errors, as given in
# Table 5.2 of Gelman, Carlin, Stern, Dunson, Vehtari and Rubin (2013).
#
#   Rscript data-raw/eight_schools.R   (from the package root)
eight_schools <- data.frame(
  school = LETTERS[1:8],
  estimate = c(28, 8, -3, 7, -1, 1, 18, 12),
  se = c(15, 10, 16, 11, 9, 11, 10, 18)
)
save(eight_schools, file = "data/eight_schools.rda", compress = "xz", version = 2)
