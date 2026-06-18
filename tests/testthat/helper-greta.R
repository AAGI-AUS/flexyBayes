# Backend-readiness gate for the greta engine.
#
# `skip_if_not_installed("greta")` only confirms the R package is present.
# greta additionally needs a working Python / TensorFlow stack at run time,
# which a fresh, sandboxed, or CI environment frequently lacks even when the R
# package installs fine. A test that reaches real greta execution must skip
# (not fail) when the backend is installed-but-unusable; otherwise a reviewer
# running on such a machine sees a wall of errors of the form
# "the expected python packages are not available".
#
# `.fb_greta_usable()` probes actual usability via `tensorflow::tf_version()`
# (the reachability of the Python TF stack greta requires) rather than building
# a greta graph, so it carries no greta global state and cannot perturb the
# emit-state isolation the other helpers rely on. The probe is cached for the
# session: it runs once, on the first gate that needs it.

.fb_greta_usable <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) {
      return(cached)
    }
    usable <- requireNamespace("greta", quietly = TRUE) &&
      requireNamespace("tensorflow", quietly = TRUE) &&
      isTRUE(tryCatch(
        !is.null(tensorflow::tf_version()),
        error = function(e) FALSE
      ))
    cached <<- usable
    usable
  }
})

# Skip a test that reaches real greta execution unless the greta backend is
# both installed and usable. Use this everywhere a test fits through greta.
skip_if_greta_backend_unusable <- function() {
  testthat::skip_if_not_installed("greta")
  # greta's Python / TensorFlow stack is intentionally absent on the CI matrix.
  # Skip on CI *before* the usability probe: on Windows, calling
  # `tensorflow::tf_version()` triggers reticulate's uv-based auto-provisioning
  # (a large, slow TensorFlow download) that fails the job non-deterministically
  # rather than returning FALSE the way it does on Linux / macOS. Skipping first
  # keeps the probe -- and the download -- off CI entirely.
  testthat::skip_on_ci()
  if (!.fb_greta_usable()) {
    testthat::skip("greta backend unusable (Python / TensorFlow stack unavailable)")
  }
}

# Canonical project-wide gate for greta-dependent tests. Skips on CRAN and on
# the standard CI matrix (where greta's Python stack is intentionally absent),
# and otherwise skips cleanly when greta is installed-but-unusable. Replaces the
# per-file definitions that previously rolled their own (and only checked
# install state).
skip_if_no_greta <- function() {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  skip_if_greta_backend_unusable()
}
