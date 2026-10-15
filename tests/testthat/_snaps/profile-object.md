# the printed profile and its summary keep their words

    Code
      print(prof, use_cli = FALSE)
    Output
      Reliability profile of 6 sites at sigma_tau = 0.2 (supplied)
      Scale of the effects: standardized mean difference
      
        informativeness I            0.571
        mean reliability             0.562
        dispersion H                 0.0720
        logit variance V             0.609
        concentration count J_eff    5.5 of 6
      
      Site reliabilities range from 0.31 to 0.80.

---

    Code
      print(summary(prof), use_cli = FALSE)
    Output
      Reliability profile of 6 sites at sigma_tau = 0.2 (supplied)
      Scale of the effects: standardized mean difference
      
      Centroids
        harmonic centroid            0.501    (harmonic mean)
        informativeness I            0.571    (logit mean)
        mean reliability             0.562    (arithmetic mean)
      
      Dispersion and concentration
        dispersion H                 0.0720
        logit variance V             0.609
        concentration count J_eff    5.5 of 6
      
      Second-order terms (Lee, 2026)
                                                                leading term
                                exact      leading    bound     decides the sign
        H                       0.0720     0.0746     0.00933   -
        mean reliability - I   -0.00922   -0.0107     0.0121    no
      
      Related summaries
        benchmark risk sigma_tau^2 (1 - mean S_j)    0.0175
        sum of reliabilities                         3.4
        typical variance (Higgins-Thompson)          0.0251   R_typ 0.614

---

    Code
      print(reliability_profile(se = 0.2, sigma_tau = 0.3), use_cli = FALSE)
    Output
      Reliability profile of 1 site at sigma_tau = 0.3 (supplied)
      
        informativeness I            0.692
        mean reliability             0.692
        dispersion H                 0.0000
        logit variance V             0.000
        concentration count J_eff    1.0 of 1
      
      The site's reliability is 0.69.

---

    Code
      print(suppressMessages(reliability_profile(se = se6, sigma_tau = 0)), use_cli = FALSE)
    Output
      Reliability profile of 6 sites at sigma_tau = 0 (supplied)
      
        informativeness I            0.000
        mean reliability             0.000
        dispersion H                 0.0000
        logit variance V             0.609
        concentration count J_eff    4.0 of 6 (limit as sigma_tau goes to 0)
      
      At sigma_tau = 0 every reliability is 0.

