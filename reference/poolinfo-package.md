# poolinfo: Reliability-Profile Diagnostics for Partial Pooling

Summarizes how much a set of estimates with known standard errors, such
as the site estimates of a multisite trial, the studies of a
meta-analysis or the direct estimates of small areas, can say about the
individual units under a normal hierarchical model. Computes the site
reliabilities and their informativeness, mean reliability, dispersion
and concentration, with posterior means and intervals from numerical
integration over the between-site standard deviation, design
calculations, plots made with 'ggplot2', and tables and sentences for
reporting. The diagnostics follow Lee (2026). Supported by the Institute
of Education Sciences through Grant R305D240078 to the University of
Alabama.

## Details

The diagnostics describe a set of site reliabilities. A site's
reliability is the proportion of the variance of its estimate that comes
from the variation of the true effects; under the normal hierarchical
model it is the weight of the site's own estimate, against the mean of
the effects, in the site's posterior mean.

The main functions:

- [`reliability_posterior()`](https://joonho112.github.io/poolinfo/reference/reliability_posterior.md)
  computes the posterior of the diagnostics from estimates and standard
  errors, by numerical integration over the between-site standard
  deviation.

- [`reliability_profile()`](https://joonho112.github.io/poolinfo/reference/reliability_profile.md)
  computes them at a given between-site standard deviation, or at the
  estimate of a metafor fit.

- [`reliability_draws()`](https://joonho112.github.io/poolinfo/reference/reliability_draws.md)
  computes them from posterior draws of another fit, such as one of
  brms.

- [`reliability_design()`](https://joonho112.github.io/poolinfo/reference/reliability_design.md),
  [`balanced_allocation()`](https://joonho112.github.io/poolinfo/reference/balanced_allocation.md)
  and
  [`solve_design()`](https://joonho112.github.io/poolinfo/reference/solve_design.md)
  plan a trial.

- [`compare_profiles()`](https://joonho112.github.io/poolinfo/reference/compare_profiles.md),
  [`report_table()`](https://joonho112.github.io/poolinfo/reference/report_table.md)
  and
  [`report_text()`](https://joonho112.github.io/poolinfo/reference/report_text.md)
  compare and report analyses.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and the
  functions `plot_*()` show them.

[`vignette("poolinfo")`](https://joonho112.github.io/poolinfo/articles/poolinfo.md)
works through an analysis, and the other vignettes cover the profiles,
the posterior summaries, a multisite trial, meta-analysis, small area
estimation, planning and the computation.

## Funding

This research was supported by the Institute of Education Sciences, U.S.
Department of Education, through Grant R305D240078 to the University of
Alabama. The opinions expressed are those of the author and do not
represent views of the Institute or the U.S. Department of Education.

## See also

Useful links:

- <https://github.com/joonho112/poolinfo>

- <https://joonho112.github.io/poolinfo/>

- Report bugs at <https://github.com/joonho112/poolinfo/issues>

## Author

**Maintainer**: JoonHo Lee <jlee296@ua.edu>
([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]

Authors:

- JoonHo Lee <jlee296@ua.edu>
  ([ORCID](https://orcid.org/0009-0006-4019-8703)) \[copyright holder\]
