# draws of sigma_tau give the diagnostics at each draw

    Code
      print(rd, use_cli = FALSE)
    Output
      Reliability profile of 8 sites (posterior draws from a test)
      3000 draws; Monte Carlo standard errors of the means in parentheses.
      Summaries are conditional on the supplied standard errors.
      Each site counts equally in the diagnostics, whatever its size.
      
                                   mean                95% interval
        informativeness I          0.297    (0.0033)   0.005    0.630
        mean reliability           0.303    (0.0032)   0.006    0.624
        dispersion H               0.0199   (0.00016)  0.0006   0.0284
        concentration count J_eff  7.3      (0.0049)   6.8      7.8      of 8
        logit variance V           0.230               (fixed by the standard errors)
        sigma_tau                  7.95     (0.072)    0.895    15.8
      
        Pr(I > 0.5)                0.149 (0.0065)
        I at E(sigma_tau^2)        0.347

