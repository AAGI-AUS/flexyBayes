# gretaR slot scaffold -- v0.2 scaffolding for the planned v0.3
# activation of the gretaR backend.
#
# Design contract: registered in the canonical-name registry
# (`R/canonical_names.R`), in lgm_gate() (`R/lgm_gate.R`), and in
# the shared dispatch helper (`R/dispatch.R`), but not advertised
# on any user-facing `backend` argument until activation. The slot
# becomes a live backend in v0.3 when (a) the gretaR package is
# public, (b) it has cleared the audit gate, and (c) the activation
# boolean `options(flexyBayes.gretaR_activated = TRUE)` is set.
#
# Until then, every internal path that would have dispatched to
# gretaR raises a structured `gretaR_dormant_refusal` error naming
# the dormancy reason + the activation procedure. v0.3 activation
# is a four-touch-point flip: widen the match.arg sets on the
# user-facing entries, swap the registry stub for the real mapper,
# replace the dispatch-helper branch body with the real emitter,
# and update the lgm_gate capability flag value.
#
# The single source of truth for activation state is the option
# `flexyBayes.gretaR_activated`, queried here only.

# ---------------------------------------------------------------- #
# Activation state introspection                                    #
# ---------------------------------------------------------------- #

# Resolve the current dormancy reason for the gretaR slot. Returns
# one of the documented sentinels; consumed by the dispatch-helper
# branch, the canonical-name mapper stub, and the gretaR_status()
# helper.
.gretaR_dormancy_reason <- function() {
  activated <- isTRUE(getOption("flexyBayes.gretaR_activated", FALSE))
  installed <- nzchar(system.file(package = "gretaR"))
  if (!installed) {
    return("gretaR_not_installed")
  }
  if (!activated) {
    return("slot_provisioned_not_activated")
  }
  # Reachable only when the activation boolean is TRUE but no audit
  # mechanism has been wired yet. v0.2 has no audit mechanism, so
  # the activation boolean alone moves the slot from
  # "slot_provisioned_not_activated" to "gretaR_not_audit_clean"
  # until the audit-status mechanism lands in v0.3.
  "gretaR_not_audit_clean"
}

# Capability flag appended by lgm_gate() to an LGM-compatible IR.
# v0.2 always returns the dormant flag because no audit mechanism
# has shipped; v0.3 will flip to the eligible flag when both the
# activation boolean and the audit mechanism agree.
.gretaR_slot_capability <- function() {
  if (identical(.gretaR_dormancy_reason(), "gretaR_dispatch_eligible")) {
    "gretaR_dispatch_eligible"
  } else {
    "gretaR_slot_dormant"
  }
}

# Standard activation procedure -- documented in one place so the
# dispatch refusal, the canonical-name stub, and gretaR_status() all
# surface the same instructions.
.gretaR_activation_procedure <- function() {
  c(
    "Install the gretaR package when it goes public on CRAN.",
    "Verify the local gretaR build against the audit checklist.",
    "Run options(flexyBayes.gretaR_activated = TRUE) in the session."
  )
}


# ---------------------------------------------------------------- #
# Canonical-name mapper stub                                        #
# ---------------------------------------------------------------- #
# Registered into .canonical_mappers under the key "gretaR" at
# package load (R/canonical_names.R bottom block). Returns the
# scaffold-dormant sentinel so triangulate() short-circuits cleanly
# against a (currently impossible) gretaR fit.

.mapper_gretaR_stub <- function(fit, fb_terms) {
  list(
    map = character(0),
    transform = list(),
    source = "scaffold_dormant",
    dormancy_reason = .gretaR_dormancy_reason()
  )
}


# ---------------------------------------------------------------- #
# Exported helper                                                   #
# ---------------------------------------------------------------- #

#' Introspect the gretaR backend slot
#'
#' Returns the activation state of the gretaR backend slot plus the
#' procedure to activate it. The gretaR slot is provisioned at v0.2
#' but not yet a live backend; activation lands in v0.3 when the
#' gretaR package is publicly available and has cleared the audit
#' gate. Until then, this helper is the canonical discoverability
#' surface --- users who see "gretaR" referenced in the backend
#' matrix or release notes can call `gretaR_status()` to inspect why
#' the slot is dormant and how to wake it.
#'
#' @return A list with components:
#'   \describe{
#'     \item{`activated`}{Logical: `TRUE` if
#'       `options(flexyBayes.gretaR_activated)` is set AND a future
#'       audit mechanism has cleared. v0.2: always `FALSE`.}
#'     \item{`gretaR_installed`}{Logical: result of
#'       `nzchar(system.file(package = "gretaR"))`.}
#'     \item{`audit_clean`}{Logical or `NA`: audit-status indicator.
#'       v0.2: always `NA` (no audit-status mechanism shipped).}
#'     \item{`dormancy_reason`}{Character: one of
#'       `"slot_provisioned_not_activated"`, `"gretaR_not_installed"`,
#'       `"gretaR_not_audit_clean"`, or `"gretaR_dispatch_eligible"`
#'       (the last meaning the slot is fully active).}
#'     \item{`activation_procedure`}{Character vector: numbered steps
#'       to activate the slot, in order.}
#'   }
#'
#' @examples
#' gs <- gretaR_status()
#' gs$activated         # v0.2: FALSE
#' gs$dormancy_reason   # v0.2: "slot_provisioned_not_activated" or
#'                      #       "gretaR_not_installed"
#' cat(gs$activation_procedure, sep = "\n")
#'
#' @export
gretaR_status <- function() {
  installed <- nzchar(system.file(package = "gretaR"))
  reason <- .gretaR_dormancy_reason()
  # v0.2: the slot is fully active only when the dormancy_reason
  # returns the eligibility sentinel. No code path produces that
  # value at v0.2 (no audit mechanism); the boolean exists to keep
  # the activation contract self-documenting.
  fully_active <- identical(reason, "gretaR_dispatch_eligible")
  list(
    activated = fully_active,
    gretaR_installed = installed,
    audit_clean = NA,
    dormancy_reason = reason,
    activation_procedure = .gretaR_activation_procedure()
  )
}


# ---------------------------------------------------------------- #
# Dispatch-helper refusal constructor                               #
# ---------------------------------------------------------------- #
# Called from .dispatch_backend() (R/dispatch.R) when the internal
# value backend = "gretaR" reaches the dispatch branch. The branch
# is unreachable from user-facing entries at v0.2 (match.arg
# rejects); reachable only from internal callers + tests. v0.3
# activation replaces the branch body with the real emit_gretaR()
# call; this helper is removed at that time.

.gretaR_dormant_refusal <- function(call = NULL) {
  reason <- .gretaR_dormancy_reason()
  proc <- .gretaR_activation_procedure()
  msg <- paste0(
    "The gretaR backend slot is currently dormant (reason: ",
    reason,
    "). Recommended fallback: pass backend = \"greta\" ",
    "for the same model. Activation procedure:\n",
    paste0("  ", seq_along(proc), ". ", proc, collapse = "\n")
  )
  cond <- structure(
    class = c("gretaR_dormant_refusal", "error", "condition"),
    list(
      message = msg,
      call = call,
      dormancy_reason = reason,
      activation_procedure = proc
    )
  )
  stop(cond)
}


# ---------------------------------------------------------------- #
# Package-load registration                                         #
# ---------------------------------------------------------------- #
# Sourced after R/canonical_names.R (alphabetic 'g' > 'c'), so the
# .canonical_mappers registry + register_canonical_mapper() are
# defined by the time this line runs. v0.3 activation: replace
# .mapper_gretaR_stub with .mapper_gretaR_real in one place here.

# Activated: the real mapper (defined in R/emit_gretaR.R, sourced first,
# alphabetic 'e' < 'g'). gretaR draws are already canonically named by the
# worker's model_from_arrays(names=), so the mapper is a near-identity.
register_canonical_mapper("gretaR", .mapper_gretaR_real)
