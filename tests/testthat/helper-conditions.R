# An error of the package, of a given class. Returns the condition so that a
# test can look at its parts.
expect_rp_error <- function(expr, class = "rp_input_error", regexp = NULL) {
  err <- rlang::catch_cnd(expr, classes = "error")
  testthat::expect_s3_class(err, "rp_error")
  testthat::expect_s3_class(err, class)
  if (!is.null(regexp)) {
    testthat::expect_match(conditionMessage(err), regexp)
  }
  invisible(err)
}

# The name of the function that an error names.
error_call_name <- function(err) {
  if (is.null(err$call)) NA_character_ else as.character(err$call[[1L]])
}
