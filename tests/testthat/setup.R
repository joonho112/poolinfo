# A partial match of `$` on a list (x$site finding x$sites) is a warning in
# the tests, so that it cannot pass unnoticed.
withr::local_options(warnPartialMatchDollar = TRUE, .local_envir = teardown_env())
