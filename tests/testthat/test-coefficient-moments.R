test_that("coefficient moments require the absolute moments of the marginal law", {
  for (engine in c("nodes", "adaptive")) {
    for (p in 1:2) {
      for (extra in 2:4) {
        J <- p + extra
        d <- data.frame(y = sin(seq_len(J)), se = seq(0.7, 1.2, length.out = J),
                        x = seq(-1, 1, length.out = J))
        fit <- reliability_posterior(y, se = se, data = d,
          mean_model = if (p == 1) ~ 1 else ~ x,
          control = posterior_control(engine = engine))
        tab <- fit$summary[grepl("^(effect_mean|coef_)", fit$summary$quantity), ]
        expect_identical(is.na(tab$mean), rep(extra == 2, p))
        expect_identical(is.na(tab$sd), rep(extra < 4, p))
        expect_true(all(is.finite(tab$lower) & is.finite(tab$upper)))
        expect_true(all(is.finite(fit$sites$effect) & is.finite(fit$sites$effect_sd)))
        expect_true(all(fit$checks$status == "ok"))
        if (extra == 2) {
          lines <- format(fit, use_cli = FALSE)
          effect_lines <- lines[grepl("^  (effect mean|coefficient)", lines)]
          expect_length(effect_lines, p)
          expect_true(all(grepl("n.d.", effect_lines, fixed = TRUE)))
        }
        # A proper light-tailed prior restores both moments at the same J.
        if (extra == 2) {
          proper <- update(fit, prior = sd_prior("half_normal", scale = 1))
          tab <- proper$summary[grepl("^(effect_mean|coef_)", proper$summary$quantity), ]
          expect_true(all(is.finite(tab$mean) & is.finite(tab$sd)))
        }
      }
    }
  }
})
