# Site estimates of the Online Credit Recovery Study (Rickles et al., 2024)
# for the 32 school-by-subject sites analyzed by Lee (2026).
#
# Within each site:
#   - the effect on course passing is the coefficient of the assignment
#     indicator in an ordinary least squares regression on the indicator and
#     the randomization blocks, with its HC2 standard error;
#   - the adjusted effect adds the grade 8 state test score in the subject,
#     the grade 9 GPA (both with means imputed) and their indicators of
#     missing values;
#   - the effect on credits by the end of year 4 is the passing regression
#     with the credits in the subject as the outcome;
#   - the unadjusted effect is the difference of the passing proportions,
#     with the usual standard error;
#   - the log odds ratio of passing has the Woolf standard error, with 0.5
#     added to each cell of the site's table where a cell is 0.
#
# The public-use file is not part of the package. It is available from
# openICPSR, project 202181 (https://www.openicpsr.org/openicpsr/project/202181),
# as OCRS_publicuse_impactanalysisdata.csv. Point OCRS_DATA at it and run the
# script from the package root:
#
#   OCRS_DATA=/path/to/OCRS_publicuse_impactanalysisdata.csv Rscript data-raw/ocrs.R
#
# The sites are labeled S01 to S32, Algebra 1 first, each subject in the
# order of the study's school identifiers; the schools are numbered 1 to 24
# in the same order. The study's identifiers are not kept.
path <- Sys.getenv("OCRS_DATA")
if (!nzchar(path) || !file.exists(path)) {
  stop("Set OCRS_DATA to the path of OCRS_publicuse_impactanalysisdata.csv.", call. = FALSE)
}
d <- utils::read.csv(path, stringsAsFactors = FALSE)
d$subject <- ifelse(d$crssub == "A1", "Algebra 1", "English 9")
d$site <- paste(d$studyschids, d$subject, sep = " x ")
d$credits_y4 <- ifelse(d$subject == "Algebra 1", d$cby_y4_math, d$cby_y4_eng)
d$pre <- ifelse(d$subject == "Algebra 1", d$state_math_imp, d$state_ela_imp)
needed <- c("studyschids", "rablock", "ragrp", "scm_passcrs", "credits_y4", "pre", "state_miss",
            "gpa9th_imp", "gpa9th_miss")
stopifnot(!anyNA(d[needed]), all(d$ragrp %in% 0:1), all(d$scm_passcrs %in% 0:1))

# OLS with HC2 standard errors; the coefficient of `ragrp`.
hc2_ols <- function(X, y, coef_name = "ragrp") {
  qr_x <- qr(X)
  keep <- qr_x$pivot[seq_len(qr_x$rank)]
  Xk <- X[, keep, drop = FALSE]
  XtXi <- chol2inv(chol(crossprod(Xk)))
  b <- XtXi %*% crossprod(Xk, y)
  e <- as.numeric(y - Xk %*% b)
  h <- rowSums((Xk %*% XtXi) * Xk)
  meat <- crossprod(Xk * (e^2 / pmax(1 - h, 1e-8)), Xk)
  V <- XtXi %*% meat %*% XtXi
  i <- match(coef_name, colnames(Xk))
  c(estimate = b[i], se = sqrt(V[i, i]))
}
site_effect <- function(dd, outcome, covars = NULL) {
  X <- cbind(`(Intercept)` = 1, ragrp = dd$ragrp)
  if (length(unique(dd$rablock)) > 1L) {
    X <- cbind(X, stats::model.matrix(~ factor(rablock), data = dd)[, -1, drop = FALSE])
  }
  if (!is.null(covars)) X <- cbind(X, as.matrix(dd[, covars, drop = FALSE]))
  hc2_ols(X, dd[[outcome]])
}
site_row <- function(dd) {
  pass <- site_effect(dd, "scm_passcrs")
  adjusted <- site_effect(dd, "scm_passcrs", covars = c("pre", "state_miss", "gpa9th_imp", "gpa9th_miss"))
  credits <- site_effect(dd, "credits_y4")
  y1 <- dd$scm_passcrs[dd$ragrp == 1]
  y0 <- dd$scm_passcrs[dd$ragrp == 0]
  cells <- c(sum(y1), sum(1 - y1), sum(y0), sum(1 - y0))
  correction <- if (min(cells) == 0) 0.5 else 0
  cells <- cells + correction
  data.frame(
    key = dd$site[[1L]], study_school = dd$studyschids[[1L]], subject = dd$subject[[1L]],
    n = nrow(dd), n_treatment = length(y1), n_control = length(y0),
    blocks = length(unique(dd$rablock)),
    estimate = pass[["estimate"]], se = pass[["se"]],
    estimate_adjusted = adjusted[["estimate"]], se_adjusted = adjusted[["se"]],
    estimate_credits = credits[["estimate"]], se_credits = credits[["se"]],
    estimate_unadjusted = mean(y1) - mean(y0),
    se_unadjusted = sqrt(stats::var(y1) / length(y1) + stats::var(y0) / length(y0)),
    log_odds_ratio = log(cells[[1L]] / cells[[2L]]) - log(cells[[3L]] / cells[[4L]]),
    se_log_odds_ratio = sqrt(sum(1 / cells)),
    zero_cell = correction > 0
  )
}
sites <- do.call(rbind, lapply(split(d, d$site), site_row))
sites <- sites[order(sites$subject, sites$study_school), ]
rownames(sites) <- NULL

# The analytic sample of the paper: 24 schools and 2 subjects give 32 sites
# with both arms, 24 to 104 students (median 49), at least 8 per arm and 144
# randomization blocks in all.
stopifnot(
  nrow(sites) == 32L, sum(sites$subject == "Algebra 1") == 13L, length(unique(sites$study_school)) == 24L,
  min(sites$n) == 24L, max(sites$n) == 104L, stats::median(sites$n) == 49,
  min(sites$n_treatment, sites$n_control) >= 8L, sum(sites$blocks) == 144L, sum(sites$n) == nrow(d),
  all(c(sites$se, sites$se_adjusted, sites$se_credits, sites$se_unadjusted, sites$se_log_odds_ratio) > 0)
)

schools <- sort(unique(sites$study_school))
ocrs <- data.frame(
  site = sprintf("S%02d", seq_len(nrow(sites))),
  school = match(sites$study_school, schools),
  sites[, c("subject", "n", "n_treatment", "n_control", "blocks", "estimate", "se", "estimate_adjusted",
            "se_adjusted", "estimate_credits", "se_credits", "estimate_unadjusted", "se_unadjusted",
            "log_odds_ratio", "se_log_odds_ratio", "zero_cell")]
)
integer_columns <- c("school", "n", "n_treatment", "n_control", "blocks")
ocrs[integer_columns] <- lapply(ocrs[integer_columns], as.integer)

save(ocrs, file = "data/ocrs.rda", compress = "xz", version = 2)
