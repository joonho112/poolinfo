# the design example has the expected profile

    Code
      print(d, use_cli = FALSE)
    Output
      Planned allocation: 40 sites, 2,000 individuals, kappa = 4, sigma_tau = 0.2
      Site sizes from 16 to 129 (CV 0.48)
      
                                   planned    balanced
        informativeness I          0.310      0.333
        mean reliability           0.318      0.333
        dispersion H               0.0227     0.0000
        concentration count J_eff  36.6       40.0
        geometric mean size        44.9       50.0
      
      Balance maximizes I and the mean reliability at fixed J and N when every site
      has the same kappa and the same cost per individual (Lee, 2026).

