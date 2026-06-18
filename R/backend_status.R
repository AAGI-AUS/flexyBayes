# fb_backend_status() -- user-facing backend-readiness diagnostic.
#
# Reports, for each inference backend flexyBayes can route to, whether the
# backend's package is installed and whether it is actually usable in this
# session. greta needs a working Python / TensorFlow stack at run time, so
# "installed" and "usable" differ for it; the other backends are usable as
# soon as their package is present. The function is read-only -- it probes
# availability and never builds a model graph or starts a fit.

# --- internal probe ------------------------------------------------- #

# Absorb R-level stdout, messages, and warnings from an invasive probe so a
# readiness check never leaks raw infrastructure noise into the console. Note:
# output written by a *subprocess* at the OS file-descriptor level (for example
# a Python launcher) cannot be captured from R; the `deep = FALSE` path exists
# precisely so such a probe is never triggered.
.fb_capture_quiet <- function(expr) {
  out_con <- textConnection("fb_quiet_out", open = "w", local = TRUE)
  sink(out_con)
  sink(out_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(out_con)
  }, add = TRUE)
  suppressWarnings(suppressMessages(force(expr)))
}

# greta's run-time usability is the R package PLUS a reachable Python /
# TensorFlow stack. This mirrors the gate the test suite and dispatch use, so
# the user sees the same verdict the package will act on. Probing via
# tensorflow::tf_version() carries no greta global state, but it does initialise
# the Python / TensorFlow stack, which can emit low-level loader noise -- so the
# probe is wrapped in output capture, and `deep = FALSE` skips it entirely
# (returning NA, "installed but not probed").
.greta_backend_usable <- function(deep = TRUE) {
  if (!requireNamespace("greta", quietly = TRUE) ||
      !requireNamespace("tensorflow", quietly = TRUE)) {
    return(FALSE)
  }
  if (!isTRUE(deep)) {
    return(NA)
  }
  isTRUE(tryCatch(
    .fb_capture_quiet(!is.null(tensorflow::tf_version())),
    error = function(e) FALSE
  ))
}

# --- exported diagnostic -------------------------------------------- #

#' Report inference-backend readiness
#'
#' Checks which inference backends flexyBayes can route to are installed and
#' usable in the current session, returning a small table you can inspect
#' before fitting. The check is read-only: it probes package availability and,
#' for greta, the reachability of its Python / TensorFlow stack, without
#' building a model or starting a fit.
#'
#' `installed` records whether the backend's R package is present. `usable`
#' records whether the backend can actually run a fit now -- for greta this
#' additionally requires a working Python / TensorFlow stack, so a greta that
#' is `installed` but not `usable` needs `greta::install_greta_deps()`. The
#' dormant opt-in `gretaR` engine is reported separately by [gretaR_status()].
#'
#' @return A data frame of class `fb_backend_status` with one row per backend
#'   and the columns `backend`, `installed` (logical), `usable` (logical), and
#'   `note` (a human-readable status, including the install command when a
#'   backend is absent). A `print` method renders it as a readiness table.
#'
#' @seealso [flexybayes()] for the universal entry, the `fb_greta()` /
#'   `fb_inla()` / `fb_brms()` single-engine pins, and [gretaR_status()] for
#'   the dormant gretaR slot.
#'
#' @examples
#' \donttest{
#' # Probing greta initialises its Python / TensorFlow stack, so the first
#' # call can be slow; wrapped in \donttest{} for that reason.
#' fb_backend_status()
#' }
#'
#' @param deep Logical; if `TRUE` (default) greta's readiness is probed by
#'   initialising its Python / TensorFlow stack (the probe's output is captured,
#'   never leaked). If `FALSE`, greta is reported as installed without touching
#'   Python -- a fast, non-invasive check whose `usable` value is `NA`
#'   ("installed, not probed").
#' @export
fb_backend_status <- function(deep = TRUE) {
  greta_inst <- requireNamespace("greta", quietly = TRUE)
  greta_use <- if (greta_inst) .greta_backend_usable(deep = deep) else FALSE
  inla_inst <- requireNamespace("INLA", quietly = TRUE)
  brms_inst <- requireNamespace("brms", quietly = TRUE)

  note_greta <- if (!greta_inst) {
    paste0(
      "not installed: install.packages('greta', repos = ",
      "c('https://greta-dev.r-universe.dev', getOption('repos')))"
    )
  } else if (is.na(greta_use)) {
    "installed; Python / TensorFlow stack not probed (call with deep = TRUE)"
  } else if (!greta_use) {
    paste0(
      "installed, but the Python / TensorFlow stack is unavailable: ",
      "run greta::install_greta_deps()"
    )
  } else {
    "ready (MCMC)"
  }

  note_inla <- if (inla_inst) {
    "ready (approximate inference: integrated nested Laplace)"
  } else {
    paste0(
      "not installed: install.packages('INLA', repos = ",
      "c(getOption('repos'), ",
      "INLA = 'https://inla.r-inla-download.org/R/stable'))"
    )
  }

  note_brms <- if (brms_inst) {
    "ready (MCMC via Stan; first-call compile)"
  } else {
    "not installed: install.packages('brms')"
  }

  out <- data.frame(
    backend = c("greta", "INLA", "brms"),
    installed = c(greta_inst, inla_inst, brms_inst),
    usable = c(greta_use, inla_inst, brms_inst),
    note = c(note_greta, note_inla, note_brms),
    stringsAsFactors = FALSE
  )
  structure(out, class = c("fb_backend_status", "data.frame"))
}

#' @export
print.fb_backend_status <- function(x, ...) {
  cat("flexyBayes backend readiness\n")
  cat(strrep("-", 64L), "\n", sep = "")
  for (i in seq_len(nrow(x))) {
    mark <- if (isTRUE(x$usable[i])) {
      "ok"
    } else if (isTRUE(x$installed[i])) {
      "!!"
    } else {
      "--"
    }
    cat(sprintf("  [%s] %-7s %s\n", mark, x$backend[i], x$note[i]))
  }
  active <- x$backend %in% c("greta", "INLA", "brms")
  if (!any(x$usable[active])) {
    cat(
      "\n  No active inference backend is usable -- install at least one ",
      "of the above before fitting.\n",
      sep = ""
    )
  }
  invisible(x)
}
