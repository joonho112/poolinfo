# a comparison tabulates each analysis from its accessors

    Code
      print(cmp, use_cli = FALSE)
    Output
      Comparison of 3 analyses
      Intervals are central posterior intervals.
      Each analysis is conditional on its own standard errors.
      
        analysis   J  informativeness I     mean reliability
        posterior  8  0.225 [0.000, 0.745]  0.230
        draws      8  0.224 [0.000, 0.748]  0.228
        plug-in    8  0.338                 0.346
      
        analysis   dispersion H             concentration count J_eff
        posterior  0.0143 [0.0000, 0.0283]  7.2
        draws      0.0141 [0.0000, 0.0283]  7.2
        plug-in    0.0249                   7.4

