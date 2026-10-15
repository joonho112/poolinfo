# a posterior fit has its parts and prints its words

    Code
      print(fit, use_cli = FALSE)
    Output
      Reliability profile of 8 sites (posterior)
      Model: estimate ~ N(effect, se^2), effect ~ N(mean, sigma_tau^2),
             sites independent
      Prior: flat on the mean and on sigma_tau
      Summaries are conditional on the supplied standard errors.
      Each site counts equally in the diagnostics, whatever its size.
      
                                   mean     95% interval
        informativeness I          0.225    0.000    0.745
        mean reliability           0.230    0.000    0.735
        dispersion H               0.0143   0.0000   0.0283
        concentration count J_eff  7.2      6.8      7.9      of 8
        logit variance V           0.230    (fixed by the standard errors)
        sigma_tau                  6.58     0.243    20.8
        effect mean                7.93     -2.09    18.2
      
        Pr(I > 0.5)                0.138
        I at E(sigma_tau^2)        0.338

---

    Code
      print(summary(fit), use_cli = FALSE)
    Output
      Reliability profile of 8 sites (posterior)
      Model: estimate ~ N(effect, se^2), effect ~ N(mean, sigma_tau^2),
             sites independent
      Prior: flat on the mean and on sigma_tau
      Summaries are conditional on the supplied standard errors.
      Each site counts equally in the diagnostics, whatever its size.
      
      Diagnostics
                                   mean      sd        2.5%      50%       97.5%
        informativeness I          0.225     0.217     0.000     0.157     0.745
        mean reliability           0.230     0.215     0.000     0.167     0.735
        harmonic centroid          0.210     0.208     0.000     0.142     0.722
        dispersion H               0.0143    0.0101    0.0000    0.0144    0.0283
        concentration count J_eff  7.2       0.3       6.8       7.1       7.9
        logit variance V           0.230     (fixed by the standard errors)
      
      Between-site SD and mean
                                   mean      sd        2.5%      50%       97.5%
        sigma_tau                  6.58      5.65      0.243     5.24      20.8
        sigma_tau^2                75.2      162       0.0589    27.4      431
        effect mean                7.93      5.18      -2.09     7.89      18.2
      
      Pr(I > 0.5) = 0.138
      
      Posterior means and values at the posterior mean of sigma_tau^2
                                           posterior at E(sigma_tau^2)
        informativeness I                  0.225     0.338
        mean reliability                   0.230     0.346
        dispersion H                       0.0143    0.0249
        concentration count J_eff          7.2       7.4
      
      H is largest, at 0.0284, where sigma_tau = 12.6 (I = 0.517).
      At that sigma_tau the mean reliability equals I.
      Integration over sigma_tau: a rule of 1408 nodes.
      All numerical checks are within their limits.

