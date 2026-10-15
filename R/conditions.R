# Conditions of the package. Every error goes through .rp_abort(), every
# warning through .rp_warn() and every message through .rp_inform().
#
# An error has up to three lines: what is wrong, what was received or which
# rule applies ("i"), and a remedy that the first two lines do not make
# obvious (">"). Its class says what kind of problem it is:
#
#   rp_input_error    an argument is not one that the function accepts
#   rp_numeric_error  a computation could not reach its stated accuracy
#   rp_scope_error    the input asks for something outside the model of the
#                     package, such as correlated estimates
#
# All three also have the class rp_error. Warnings have the class rp_warning
# and messages rp_message, each with a more specific class in front of it.
# The tests check the classes, not the wording.

# The frame to name in an error or a warning: the outermost call into this
# package, so that the message names the function that the user called and
# not an internal helper.
.rp_user_frame <- function() {
  namespace <- topenv(environment(.rp_user_frame))
  for (i in seq_len(sys.nframe())) {
    fn <- sys.function(i)
    if (is.function(fn) && identical(topenv(environment(fn)), namespace)) {
      return(sys.frame(i))
    }
  }
  NULL
}

.rp_error_classes <- c("rp_input_error", "rp_numeric_error", "rp_scope_error")

.rp_abort <- function(message,
                      info = NULL,
                      fix = NULL,
                      class = "rp_input_error",
                      subclass = NULL,
                      ...,
                      parent = NULL,
                      call = .rp_user_frame()) {
  if (!class %in% .rp_error_classes) {
    rlang::abort(sprintf("Unknown error class of poolinfo: %s.", class))
  }
  rlang::abort(
    c(message, i = info, ">" = fix),
    class = c(subclass, class, "rp_error"),
    ...,
    parent = parent,
    call = call
  )
}

.rp_warn <- function(message,
                     info = NULL,
                     subclass = NULL,
                     ...,
                     call = .rp_user_frame()) {
  rlang::warn(
    c(message, i = info),
    class = c(subclass, "rp_warning"),
    ...,
    call = call
  )
}

.rp_inform <- function(message, info = NULL, subclass = NULL, ...) {
  rlang::inform(c(message, i = info), class = c(subclass, "rp_message"), ...)
}

# Messages that are shown once in a session. The keys of the ones already
# shown are kept here.
.rp_session <- new.env(parent = emptyenv())

.rp_inform_once <- function(key, message, info = NULL, subclass = NULL) {
  if (isTRUE(.rp_session[[key]])) {
    return(invisible(FALSE))
  }
  assign(key, TRUE, envir = .rp_session)
  .rp_inform(c(message, " " = "This message is shown once in a session."),
             info = info, subclass = subclass)
  invisible(TRUE)
}

# A short description of a value, for the line of an error that says what
# was received.
.rp_describe <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }
  if (is.function(x)) {
    return("a function")
  }
  if (is.data.frame(x)) {
    return(sprintf("a data frame with %d rows", nrow(x)))
  }
  if (is.factor(x)) {
    return(sprintf("a factor of length %d", length(x)))
  }
  if (is.atomic(x) && is.null(dim(x))) {
    if (length(x) != 1L) {
      kind <- if (is.double(x)) "numeric" else class(x)[[1L]]
      article <- if (grepl("^[aeiou]", kind)) "an" else "a"
      return(sprintf("%s %s vector of length %d", article, kind, length(x)))
    }
    if (is.character(x) && !is.na(x)) {
      shown <- if (nchar(x) > 40L) paste0(substr(x, 1L, 40L), "...") else x
      return(sprintf("\"%s\"", shown))
    }
    return(format(x))
  }
  if (is.matrix(x)) {
    return(sprintf("a %d x %d matrix", nrow(x), ncol(x)))
  }
  if (is.list(x) && !is.object(x)) {
    return(sprintf("a list of length %d", length(x)))
  }
  sprintf("an object of class %s", paste(class(x), collapse = "/"))
}

.rp_you_passed <- function(x) {
  sprintf("You passed %s.", .rp_describe(x))
}

# Positions as text: "3", "3 and 7", "3, 7 and 9", or the first few of many.
.rp_positions <- function(at, max_shown = 5L) {
  at <- as.integer(at)
  if (length(at) > max_shown) {
    return(sprintf("%s and %d more", paste(at[seq_len(max_shown)], collapse = ", "),
                   length(at) - max_shown))
  }
  if (length(at) == 1L) {
    return(as.character(at))
  }
  paste(paste(at[-length(at)], collapse = ", "), "and", at[[length(at)]])
}
