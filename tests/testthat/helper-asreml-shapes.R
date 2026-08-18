# helper-asreml-shapes.R -- recorded shapes of ASReml-R objects the
# na_action normaliser has to recognise.
#
# flexyBayes never depends on asreml, not even in Suggests: it is
# commercial, licensed per seat, and unavailable to most readers of this
# suite. The normaliser therefore detects an incoming `na.method()` value
# by SHAPE (a list carrying `x` and `y` character slots), never by class,
# and the shapes below are the recorded oracle it is tested against. The
# tests that consume this file run with asreml absent.
#
# Provenance of every literal in this file
# ----------------------------------------
#   Recorded : 2026-08-17, on the licensed local installation.
#   Package  : asreml 4.2.0.392 (utils::packageVersion("asreml")).
#   Commands : asreml::na.method(y = "include", x = "fail")
#              asreml::na.method()
#              asreml::na.method(y = "include")
#              asreml::na.method(y = "omit")
#              asreml::na.method(y = "fail")
#              asreml::na.method(y = "include", x = "omit")
#              asreml::na.method(y = "include", x = "include")
#              formals(asreml::na.method)
#   Method   : each value captured with dput(), transcribed verbatim
#              below. Nothing here is reconstructed from documentation.
#
# Two properties of the recorded values drive the normaliser, and neither
# is guessable from the ASReml-R manual:
#
#   1. The returned object is a PLAIN list. class() is "list", the only
#      attribute is `names`, and the names come back in the order
#      c("x", "y") -- x first, whatever order the caller supplied.
#
#   2. na.method() does NOT reduce an unsupplied argument to its default
#      scalar. It returns that argument's whole default vector, so
#      na.method(y = "include") carries x = c("fail", "include", "omit"),
#      length 3. A normaliser that assumes each slot is a scalar will
#      read the wrong policy off a partially specified call. The
#      effective value is the first element, matching match.arg().
#
# The single-argument default vectors, as recorded from
# formals(asreml::na.method):
#   y : c("include", "omit", "fail")   -> default "include"
#   x : c("fail", "include", "omit")   -> default "fail"


# .asreml_na_method_recorded() --- one recorded na.method() return value
#
# `case` selects which recorded call to return. `"explicit"` is the fully
# specified call an ASReml user writes on a field trial with missing
# plots and is the primary fixture. The remaining cases exist so the
# normaliser is tested against partial specification, where the
# unsupplied slot arrives as a length-3 vector rather than a scalar.
#
# @returns The list asreml returned for that call, transcribed verbatim
#   from dput() output.
# @noRd
# @keywords internal
.asreml_na_method_recorded <- function(
  case = c("explicit", "default", "y_include", "y_omit", "y_fail",
           "x_omit", "x_include")
) {
  case <- match.arg(case)
  switch(
    case,
    # asreml::na.method(y = "include", x = "fail")
    explicit = list(x = "fail", y = "include"),
    # asreml::na.method()
    default = list(
      x = c("fail", "include", "omit"),
      y = c("include", "omit", "fail")
    ),
    # asreml::na.method(y = "include")
    y_include = list(x = c("fail", "include", "omit"), y = "include"),
    # asreml::na.method(y = "omit")
    y_omit = list(x = c("fail", "include", "omit"), y = "omit"),
    # asreml::na.method(y = "fail")
    y_fail = list(x = c("fail", "include", "omit"), y = "fail"),
    # asreml::na.method(y = "include", x = "omit")
    x_omit = list(x = "omit", y = "include"),
    # asreml::na.method(y = "include", x = "include")
    x_include = list(x = "include", y = "include")
  )
}


# .asreml_na_method_defaults() --- the recorded argument defaults
#
# Transcribed from formals(asreml::na.method) on the same run. The first
# element of each vector is ASReml's own default for that argument, which
# is what a normaliser must take when the slot arrives unreduced.
#
# @returns A named list with character vectors `y` and `x`.
# @noRd
# @keywords internal
.asreml_na_method_defaults <- function() {
  list(
    y = c("include", "omit", "fail"),
    x = c("fail", "include", "omit")
  )
}


# .asreml_na_method_version() --- the asreml build the shapes came from
#
# Recorded so a future re-grounding run can say which version it is
# re-checking against rather than diffing against an unnamed baseline.
#
# @returns A single string.
# @noRd
# @keywords internal
.asreml_na_method_version <- function() {
  "4.2.0.392"
}
