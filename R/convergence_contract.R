# Convergence contract --- one tri-state verdict, read by every consumer
#
# Closes audit A25. Before this file there were two consumers of the
# convergence evidence and they made opposite mistakes on the same absent
# field:
#
#   run_one.R:82   isTRUE(s$converge$ok %||% TRUE)
#                  `ok` never existed, and NULL %||% TRUE is TRUE, so this
#                  was not a weak measurement -- it was a CONSTANT. No input
#                  produced FALSE, and all eight banked fixtures carry a
#                  field that measured nothing.
#
#   test-metamorphic-scale-invariance.R:203
#                  expect_false(isTRUE(summary(fit)$converge$ok))
#                  isTRUE(NULL) is FALSE, so this passed for every fit --- a
#                  clean one, a failed one, and one with no diagnostics at
#                  all were indistinguishable to it.
#
# There is a third, quieter version of the same error inside the package:
# .fb_inla_diagnostics() computed ok = isTRUE(nc$pass), and nc$pass is
# deliberately NA when no numerical confirm ran. isTRUE(NA) is FALSE, so a
# fit whose confirmation was ABSENT was reported as having FAILED. Absent
# and failed are different answers and only one of them is true.
#
# The contract here is tri-state and says so:
#
#   TRUE   evidence exists and all of it passes
#   FALSE  some evidence says the fit did not converge
#   NA     the evidence needed to decide is not there
#
# Missing evidence never becomes TRUE, and it never becomes FALSE either.


# --- evidence items -------------------------------------------------- #

# Each check returns "pass", "fail" or "absent" and a reason when it fails.
# Keeping the vote and the reason together is what lets the verdict say
# WHICH evidence failed rather than only that something did.
.fb_conv_item <- function(state, reason = NULL) {
  list(state = state, reason = reason)
}

.fb_conv_check_mode_status <- function(ms) {
  if (is.null(ms) || length(ms) != 1L || !is.finite(ms)) {
    return(.fb_conv_item("absent"))
  }
  if (ms == 0) {
    return(.fb_conv_item("pass"))
  }
  .fb_conv_item("fail", paste0("INLA mode status is ", ms, ", not 0"))
}

.fb_conv_check_mlik <- function(mlik) {
  if (is.null(mlik) || length(mlik) != 1L || is.na(mlik)) {
    return(.fb_conv_item("absent"))
  }
  if (is.finite(mlik)) {
    return(.fb_conv_item("pass"))
  }
  .fb_conv_item("fail", "the marginal log-likelihood is not finite")
}

.fb_conv_check_confirm <- function(pass, reasons = character(0)) {
  if (is.null(pass) || length(pass) != 1L || is.na(pass)) {
    return(.fb_conv_item("absent"))
  }
  if (isTRUE(pass)) {
    return(.fb_conv_item("pass"))
  }
  .fb_conv_item("fail", paste0(
    "the numerical confirmation failed",
    if (length(reasons) > 0L) paste0(": ", paste(reasons, collapse = "; ")) else ""
  ))
}


# --- combining ------------------------------------------------------- #

# .fb_conv_combine() --- the tri-state rule, in one place.
#
# `required` names the items that must be PRESENT for a TRUE to be
# defensible. A fit whose mode status is unknown has not been shown to
# converge, whatever else passed, so the verdict is NA rather than TRUE.
# Any failure anywhere is decisive: FALSE wins over a missing requirement,
# because evidence of failure is still evidence.
.fb_conv_combine <- function(items, required = character(0)) {
  states <- vapply(items, function(i) i$state, character(1L))
  reasons <- unlist(lapply(items, function(i) i$reason), use.names = FALSE)
  reasons <- reasons[!is.na(reasons)]

  if (any(states == "fail")) {
    return(list(ok = FALSE, reasons = as.character(reasons)))
  }

  missing_required <- intersect(names(states)[states == "absent"], required)
  if (length(missing_required) > 0L) {
    return(list(ok = NA, reasons = paste0(
      "no evidence for: ", paste(missing_required, collapse = ", ")
    )))
  }

  if (all(states == "absent")) {
    return(list(ok = NA, reasons = "this fit records no convergence diagnostics"))
  }

  list(ok = TRUE, reasons = character(0))
}


# --- per-engine verdicts ---------------------------------------------- #

.fb_conv_verdict_inla <- function(fit) {
  ms <- tryCatch(fit$inla$mode$mode.status, error = function(e) NULL)
  mlik <- tryCatch(
    if (is.null(fit$inla$mlik)) NULL else fit$inla$mlik[1L, 1L],
    error = function(e) NULL
  )
  nc <- fit$num_check
  # The aggregated route builds its own object and skips the post-fit
  # confirm, so read it off the raw INLA object rather than calling the
  # evidence absent when it is merely elsewhere.
  if (is.null(nc) && !is.null(fit$inla)) {
    nc <- tryCatch(.lgm_check_numerical(fit$inla), error = function(e) NULL)
  }

  items <- list(
    mode_status = .fb_conv_check_mode_status(ms),
    mlik = .fb_conv_check_mlik(mlik),
    numerical_confirm = .fb_conv_check_confirm(nc$pass, nc$reasons %||% character(0))
  )

  # mode_status is required: without it there is no basis for TRUE. The
  # numerical confirm is not, because the aggregated route legitimately has
  # none and a fit is not unverified merely for taking that route.
  v <- .fb_conv_combine(items, required = "mode_status")
  c(
    list(engine = "inla"), v,
    list(evidence = list(
      mode_status = if (is.null(ms)) NA_real_ else as.numeric(ms),
      mlik = if (is.null(mlik)) NA_real_ else as.numeric(mlik),
      numerical_confirm = if (is.null(nc$pass)) NA else nc$pass,
      max_rhat = NA_real_, min_ess_bulk = NA_real_, n_divergent = NA_integer_
    ))
  )
}

.fb_conv_verdict_brms <- function(fit) {
  d <- .fb_brms_diagnostics(fit)
  # .fb_brms_diagnostics() already applies the package's own thresholds and
  # already reports NA where the sampler diagnostics are absent, so the
  # verdict adopts its answer rather than re-deriving one with a second set
  # of thresholds that could drift from the first.
  list(
    engine = "brms",
    ok = d$ok,
    reasons = d$reasons,
    evidence = list(
      mode_status = NA_real_, mlik = NA_real_, numerical_confirm = NA,
      max_rhat = d$max_rhat,
      min_ess_bulk = d$min_ess_bulk,
      n_divergent = d$n_divergent
    )
  )
}


# --- the one function every consumer calls ---------------------------- #

# .fb_convergence_verdict() --- the authoritative tri-state answer.
#
# summary(), triangulate(), the ASReml oracle runner, the live-backend tests
# and fb_select()'s fold gate all read this and nothing else. A second
# implementation anywhere is how A25 happened.
.fb_convergence_verdict <- function(fit) {
  if (inherits(fit, "flexybayes_inla")) {
    return(.fb_conv_verdict_inla(fit))
  }
  if (inherits(fit, "flexybayes_brms")) {
    return(.fb_conv_verdict_brms(fit))
  }
  list(
    engine = tryCatch(.triangulate_source(fit), error = function(e) "unknown"),
    ok = NA,
    reasons = "this fit records no convergence diagnostics",
    evidence = list(
      mode_status = NA_real_, mlik = NA_real_, numerical_confirm = NA,
      max_rhat = NA_real_, min_ess_bulk = NA_real_, n_divergent = NA_integer_
    )
  )
}
