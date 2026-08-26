# fb_backend_status() -- user-facing backend-readiness diagnostic.
#
# Reports, for each inference backend flexyBayes can route to, whether the
# backend's package is installed and whether it is actually usable in this
# session. The function is read-only -- it probes availability and never
# builds a model graph or starts a fit.

# --- exported diagnostic -------------------------------------------- #

#' Report inference-backend readiness
#'
#' Checks which inference backends flexyBayes can route to are installed and
#' usable in the current session, returning a small table you can inspect
#' before fitting. The check is read-only: it probes package availability
#' without building a model or starting a fit.
#'
#' `installed` records whether the backend's R package is present. `usable`
#' records whether the backend can actually run a fit now.
#'
#' @return A data frame of class `fb_backend_status` with one row per backend
#'   and the columns `backend`, `installed` (logical), `usable` (logical), and
#'   `note` (a human-readable status, including the install command when a
#'   backend is absent). A `print` method renders it as a readiness table.
#'
#' @seealso [flexybayes()] for the universal entry, and the `fb_inla()` /
#'   `fb_brms()` single-engine pins.
#'
#' @examples
#' fb_backend_status()
#'
#' @param deep Logical; currently unused (reserved for a future backend whose
#'   readiness probe is expensive enough to need an opt-out). Default `TRUE`.
#' @export
fb_backend_status <- function(deep = TRUE) {
  inla_inst <- requireNamespace("INLA", quietly = TRUE)
  brms_inst <- requireNamespace("brms", quietly = TRUE)

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
    backend = c("INLA", "brms"),
    installed = c(inla_inst, brms_inst),
    usable = c(inla_inst, brms_inst),
    note = c(note_inla, note_brms),
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
  if (!any(x$usable)) {
    cat(
      "\n  No active inference backend is usable -- install at least one ",
      "of the above before fitting.\n",
      sep = ""
    )
  }
  invisible(x)
}
