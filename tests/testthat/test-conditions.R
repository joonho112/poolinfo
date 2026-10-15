test_that("errors carry the class of their kind and rp_error", {
  for (kind in c("rp_input_error", "rp_numeric_error", "rp_scope_error")) {
    err <- expect_rp_error(.rp_abort("Something is wrong.", class = kind), kind)
    expect_s3_class(err, "rlang_error")
  }
  err <- expect_rp_error(
    .rp_abort("Something is wrong.", class = "rp_numeric_error", subclass = "rp_tail_error"),
    "rp_tail_error"
  )
  expect_identical(class(err)[1:3], c("rp_tail_error", "rp_numeric_error", "rp_error"))
})

test_that("an unknown class of error is refused", {
  expect_error(.rp_abort("x", class = "rp_other_error"), "Unknown error class")
})

test_that("an error has a header, an information line and a remedy", {
  err <- expect_rp_error(
    .rp_abort("What is wrong.", info = "What was received.", fix = "What to do.")
  )
  lines <- strsplit(conditionMessage(err), "\n", fixed = TRUE)[[1L]]
  expect_length(lines, 3L)
  expect_identical(lines[[1L]], "What is wrong.")
  expect_match(lines[[2L]], "What was received.$")
  expect_match(lines[[3L]], "What to do.$")
})

test_that("an error without information or remedy has one line", {
  err <- expect_rp_error(.rp_abort("Only this."))
  expect_length(strsplit(conditionMessage(err), "\n", fixed = TRUE)[[1L]], 1L)
})

test_that("an error names the outermost function of the package that was called", {
  err <- expect_rp_error(.rp_check_level(2))
  expect_identical(error_call_name(err), ".rp_check_level")
  outer_fn <- function() .rp_check_level(2)
  err <- expect_rp_error(outer_fn())
  expect_identical(error_call_name(err), ".rp_check_level")
})

test_that("warnings and messages have their classes", {
  w <- rlang::catch_cnd(.rp_warn("Careful.", subclass = "rp_warn_test"), classes = "warning")
  expect_s3_class(w, "rp_warn_test")
  expect_s3_class(w, "rp_warning")
  m <- rlang::catch_cnd(.rp_inform("Note.", subclass = "rp_note_test"), classes = "message")
  expect_s3_class(m, "rp_note_test")
  expect_s3_class(m, "rp_message")
})

test_that("a message for once in a session is shown once", {
  key <- "test_once_key"
  withr::defer(rm(list = key, envir = .rp_session))
  expect_message(shown <- .rp_inform_once(key, "First time."), class = "rp_message")
  expect_true(shown)
  expect_no_message(shown <- .rp_inform_once(key, "First time."))
  expect_false(shown)
})

test_that("values are described for the information line", {
  expect_identical(.rp_describe(NULL), "NULL")
  expect_identical(.rp_describe(2), "2")
  expect_identical(.rp_describe("a"), "\"a\"")
  expect_identical(.rp_describe(c(1, 2)), "a numeric vector of length 2")
  expect_identical(.rp_describe(c(TRUE, FALSE)), "a logical vector of length 2")
  expect_identical(.rp_describe(1:3), "an integer vector of length 3")
  expect_identical(.rp_describe(data.frame(a = 1:2)), "a data frame with 2 rows")
  expect_identical(.rp_describe(matrix(0, 2, 3)), "a 2 x 3 matrix")
  expect_identical(.rp_describe(list(1, 2)), "a list of length 2")
  expect_identical(.rp_describe(mean), "a function")
  expect_identical(.rp_describe(structure(list(), class = "foo")), "an object of class foo")
  expect_identical(.rp_positions(3), "3")
  expect_identical(.rp_positions(c(3, 7)), "3 and 7")
  expect_identical(.rp_positions(c(3, 7, 9)), "3, 7 and 9")
  expect_identical(.rp_positions(1:8), "1, 2, 3, 4, 5 and 3 more")
})
