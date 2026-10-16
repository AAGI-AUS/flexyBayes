# Deprecation conditions, in base R.
#
# A deprecated export, argument or keyword warns with a condition of class
# `flexybayes_deprecated`. The class vector also keeps
# `lifecycle_warning_deprecated`, the class these warnings carried before
# 0.13.0, so a handler written for it still fires, and base R's
# `deprecatedWarning`, the class `.Deprecated()` signals. A defunct argument
# stops with base R's `defunctError` and keeps `lifecycle_error_deprecated`
# for the same reason. Every call warns: there is no once-per-session latch.
# A caller running compatibility checks may explicitly set
# `options(flexyBayes.deprecation_warnings = FALSE)`. Ordinary sessions default
# to warnings on.

.fb_warn_deprecated <- function(message) {
  if (identical(getOption("flexyBayes.deprecation_warnings", TRUE), FALSE)) {
    return(invisible(NULL))
  }
  warning(warningCondition(
    paste(message, collapse = "\n"),
    class = c(
      "flexybayes_deprecated",
      "lifecycle_warning_deprecated",
      "deprecatedWarning"
    ),
    call = NULL
  ))
}

.fb_stop_defunct <- function(message) {
  stop(errorCondition(
    paste(message, collapse = "\n"),
    class = c("lifecycle_error_deprecated", "defunctError"),
    call = NULL
  ))
}

# The warning a deprecated export raises: the export, the release that
# deprecated it, and what to call.
.fb_warn_deprecated_export <- function(name, use) {
  .fb_warn_deprecated(c(
    paste0(
      "`",
      name,
      "()` is deprecated as of flexyBayes 0.13.0 and will be ",
      "removed in a later release."
    ),
    use
  ))
}
