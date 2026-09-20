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
    if (length(reasons) > 0L) {
      paste0(": ", paste(reasons, collapse = "; "))
    } else {
      ""
    }
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
    return(list(
      ok = NA,
      reasons = "this fit records no convergence diagnostics"
    ))
  }

  list(ok = TRUE, reasons = character(0))
}


# --- per-engine verdicts ---------------------------------------------- #

# The marginal likelihood a verdict reports: the fit's evidence record where
# it has one (R/model_manifest.R), under which an aggregated fit has none,
# and the engine's value otherwise. The mode check above still reads the
# engine's value, which is finite on every route.
.fb_conv_reported_mlik <- function(fit, mlik) {
  evidence <- .fb_manifest(fit)[["evidence"]]
  if (!is.null(evidence)) {
    if (identical(evidence[["status"]], "available")) {
      return(evidence[["value"]])
    }
    return(NA_real_)
  }
  if (is.null(mlik)) NA_real_ else as.numeric(mlik)
}

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
    numerical_confirm = .fb_conv_check_confirm(
      nc$pass,
      nc$reasons %||% character(0)
    )
  )

  # mode_status is required: without it there is no basis for TRUE. The
  # numerical confirm is not, because the aggregated route legitimately has
  # none and a fit is not unverified merely for taking that route.
  v <- .fb_conv_combine(items, required = "mode_status")
  c(
    list(engine = "inla"), v,
    list(evidence = list(
      mode_status = if (is.null(ms)) NA_real_ else as.numeric(ms),
      mlik = .fb_conv_reported_mlik(fit, mlik),
      numerical_confirm = if (is.null(nc$pass)) NA else nc$pass,
      max_rhat = NA_real_, min_ess_bulk = NA_real_, n_divergent = NA_integer_,
      # The sampler items, carried as absent so the evidence list has one
      # shape on both engines and a consumer can read a field without
      # first asking which engine answered.
      min_ess_tail = NA_real_, min_ebfmi = NA_real_,
      n_max_treedepth = NA_integer_, max_treedepth = NA_real_
    ))
  )
}

# --- the summary-slot diagnostics on a brms fit ------------------------ #

# Thresholds. R-hat 1.01 and bulk ESS 400 are the rank-normalised
# workflow defaults (Vehtari et al. 2021); they are heuristics with no
# error rate.
.FB_CONV_RHAT_MAX <- 1.01
.FB_CONV_ESS_MIN <- 400

# Parameters excluded from the sampler diagnostics: Stan's own bookkeeping
# quantities are not model parameters, and their mixing is not what the
# verdict is about.
.FB_CONV_DIAG_EXCLUDE <- c("lp__", "lprior", "lp")

# Drop the bookkeeping quantities from a named diagnostic vector. An
# unnamed vector is returned whole: subsetting on `names(v) %in% ...`
# when `names(v)` is NULL empties the vector, which would silently turn
# a failing fit into one with no diagnostics at all.
.fb_drop_bookkeeping <- function(v) {
  if (is.null(v) || is.null(names(v))) {
    return(v)
  }
  v[!(names(v) %in% .FB_CONV_DIAG_EXCLUDE)]
}

# Name of the element at `idx`, or a generic label when the vector has no
# names (a fit whose diagnostics arrived without variable labels).
.fb_diag_label <- function(v, idx) {
  nms <- names(v)
  if (is.null(nms) || is.na(nms[[idx]])) "a model parameter" else nms[[idx]]
}

# max R-hat, bulk ESS and the divergence count from the summary slot
# .brms_convergence() wrote at fit time, against the thresholds above:
# ok (TRUE / FALSE / NA when the fit carries none), the reasons a FALSE
# was reached, and the numbers behind them.
.fb_brms_diagnostics <- function(fit) {
  conv <- fit$extras$convergence
  out <- list(
    engine = "brms",
    ok = NA,
    reasons = character(0),
    max_rhat = NA_real_,
    min_ess_bulk = NA_real_,
    n_divergent = NA_integer_
  )
  if (is.null(conv)) {
    out$reasons <- "this fit records no sampler diagnostics"
    return(out)
  }

  rhat <- tryCatch(
    conv$gelman$psrf[, "Point est.", drop = FALSE],
    error = function(e) numeric(0)
  )
  if (is.matrix(rhat)) {
    rhat <- stats::setNames(as.numeric(rhat[, 1L]), rownames(rhat))
  }
  rhat <- .fb_drop_bookkeeping(rhat)
  rhat <- rhat[is.finite(rhat)]
  ess <- .fb_drop_bookkeeping(conv$n_eff)
  ess <- ess[is.finite(ess)]
  ndiv <- conv$n_divergent

  reasons <- character(0)
  ok <- TRUE
  if (length(rhat) == 0L && length(ess) == 0L) {
    out$reasons <- "this fit records no sampler diagnostics"
    return(out)
  }
  if (length(rhat) > 0L) {
    out$max_rhat <- max(rhat)
    if (out$max_rhat > .FB_CONV_RHAT_MAX) {
      ok <- FALSE
      reasons <- c(reasons, sprintf(
        "max R-hat %.4f on `%s` exceeds %.2f",
        out$max_rhat,
        .fb_diag_label(rhat, which.max(rhat)),
        .FB_CONV_RHAT_MAX
      ))
    }
  }
  if (length(ess) > 0L) {
    out$min_ess_bulk <- min(ess)
    if (out$min_ess_bulk < .FB_CONV_ESS_MIN) {
      ok <- FALSE
      reasons <- c(reasons, sprintf(
        "min bulk ESS %.0f on `%s` is below %d",
        out$min_ess_bulk,
        .fb_diag_label(ess, which.min(ess)),
        as.integer(.FB_CONV_ESS_MIN)
      ))
    }
  }
  if (!is.null(ndiv) && is.finite(ndiv)) {
    out$n_divergent <- as.integer(ndiv)
    if (out$n_divergent > 0L) {
      ok <- FALSE
      reasons <- c(reasons, sprintf(
        "%d divergent transition%s",
        out$n_divergent,
        if (out$n_divergent == 1L) "" else "s"
      ))
    }
  }
  out$ok <- ok
  out$reasons <- reasons
  out
}

# --- the NUTS evidence a summary table does not carry ------------------ #
#
# max R-hat, bulk ESS and the divergence count reach the verdict through
# the summary slot .brms_convergence() wrote at fit time. Three further
# items are properties of the sampler's own trajectory and are read off
# the live brmsfit instead: treedepth saturation, E-BFMI, and tail ESS.
# They are read here rather than recomputed from the summary because the
# summary does not carry them at all -- an absent field is what audit A25
# was about, and inventing one from what is present is the same mistake in
# a new place.

# Stan's guidance, and the reason the three items are not treated alike:
# a divergence, an E-BFMI below 0.2 and a starved tail are evidence the
# posterior was not explored correctly; saturating max treedepth is
# evidence the sampler was slow, not that the answer is wrong. The first
# two join the failure vote; treedepth is recorded and reported as
# evidence, and does not flip the verdict.
.FB_CONV_EBFMI_MIN <- 0.2

# E-BFMI on one chain: the mean squared successive difference in energy
# against the energy variance (Betancourt 2016, eq. 30).
.fb_conv_ebfmi <- function(energy) {
  energy <- energy[is.finite(energy)]
  if (length(energy) < 2L) {
    return(NA_real_)
  }
  denom <- sum((energy - mean(energy))^2)
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  sum(diff(energy)^2) / denom
}

# The sampler's own treedepth ceiling, read from the arguments the fit was
# run with. Saturation is `treedepth == max_depth`, so the ceiling has to
# come from the run rather than from Stan's documented default: a fit with
# `control = list(max_treedepth = 12)` saturates at 12 and a hard-coded 10
# would count twelve-deep trajectories that never hit the ceiling.
.fb_conv_brms_max_treedepth <- function(brmsfit) {
  sa <- tryCatch(brmsfit$fit@stan_args, error = function(e) NULL)
  if (!is.list(sa) || length(sa) == 0L) {
    return(NA_real_)
  }
  md <- unlist(
    lapply(sa, function(a) a$max_depth %||% a$control$max_treedepth),
    use.names = FALSE
  )
  md <- suppressWarnings(as.numeric(md))
  md <- md[is.finite(md)]
  if (length(md) == 0L) {
    return(NA_real_)
  }
  max(md)
}

.fb_conv_brms_sampler_evidence <- function(fit) {
  out <- list(
    n_divergent = NA_integer_,
    n_max_treedepth = NA_integer_,
    max_treedepth = NA_real_,
    min_ebfmi = NA_real_,
    min_ess_tail = NA_real_
  )
  brmsfit <- fit$brms
  if (is.null(brmsfit) || !requireNamespace("brms", quietly = TRUE)) {
    return(out)
  }

  np <- tryCatch(brms::nuts_params(brmsfit), error = function(e) NULL)
  if (
    is.data.frame(np) &&
      all(c("Parameter", "Value", "Chain") %in% names(np))
  ) {
    par <- as.character(np$Parameter)
    div <- np$Value[par == "divergent__"]
    if (length(div) > 0L) {
      out$n_divergent <- as.integer(sum(div > 0, na.rm = TRUE))
    }
    td <- np$Value[par == "treedepth__"]
    if (length(td) > 0L) {
      ceiling_td <- .fb_conv_brms_max_treedepth(brmsfit)
      out$max_treedepth <- max(td, na.rm = TRUE)
      if (is.finite(ceiling_td)) {
        out$n_max_treedepth <- as.integer(sum(td >= ceiling_td, na.rm = TRUE))
      }
    }
    energy <- np[par == "energy__", , drop = FALSE]
    if (nrow(energy) > 0L) {
      per_chain <- vapply(
        split(energy$Value, energy$Chain),
        .fb_conv_ebfmi,
        numeric(1L)
      )
      per_chain <- per_chain[is.finite(per_chain)]
      if (length(per_chain) > 0L) {
        out$min_ebfmi <- min(per_chain)
      }
    }
  }

  if (requireNamespace("posterior", quietly = TRUE)) {
    tail_ess <- tryCatch(
      posterior::summarise_draws(brmsfit, "ess_tail"),
      error = function(e) NULL
    )
    if (is.data.frame(tail_ess) && "ess_tail" %in% names(tail_ess)) {
      v <- stats::setNames(
        as.numeric(tail_ess$ess_tail),
        as.character(tail_ess$variable)
      )
      v <- .fb_drop_bookkeeping(v)
      v <- v[is.finite(v)]
      if (length(v) > 0L) {
        out$min_ess_tail <- min(v)
      }
    }
  }

  out
}

.fb_conv_verdict_brms <- function(fit) {
  d <- .fb_brms_diagnostics(fit)
  # .fb_brms_diagnostics() already applies the package's own thresholds and
  # already reports NA where the sampler diagnostics are absent, so the
  # verdict adopts its answer rather than re-deriving one with a second set
  # of thresholds that could drift from the first.
  s <- .fb_conv_brms_sampler_evidence(fit)

  ok <- d$ok
  reasons <- d$reasons
  fail <- function(reason) {
    ok <<- FALSE
    reasons <<- c(reasons, reason)
  }

  # A divergence count the summary slot did not carry, taken from the
  # trajectory itself. Where both are present the summary's own count has
  # already voted, so this only fills a gap.
  if (
    !isTRUE(is.finite(d$n_divergent)) &&
      is.finite(s$n_divergent) &&
      s$n_divergent > 0L
  ) {
    fail(sprintf(
      "%d divergent transition%s",
      s$n_divergent,
      if (s$n_divergent == 1L) "" else "s"
    ))
  }
  if (is.finite(s$min_ebfmi) && s$min_ebfmi < .FB_CONV_EBFMI_MIN) {
    fail(sprintf(
      "E-BFMI %.3f on the worst chain is below %.1f",
      s$min_ebfmi,
      .FB_CONV_EBFMI_MIN
    ))
  }
  if (is.finite(s$min_ess_tail) && s$min_ess_tail < .FB_CONV_ESS_MIN) {
    fail(sprintf(
      "min tail ESS %.0f is below %d",
      s$min_ess_tail,
      as.integer(.FB_CONV_ESS_MIN)
    ))
  }
  # Missing evidence never becomes TRUE: an NA from the summary layer
  # stays NA when the trajectory evidence merely fails to contradict it.

  list(
    engine = "brms",
    ok = ok,
    reasons = as.character(reasons),
    evidence = list(
      mode_status = NA_real_, mlik = NA_real_, numerical_confirm = NA,
      max_rhat = d$max_rhat,
      min_ess_bulk = d$min_ess_bulk,
      n_divergent = if (is.finite(d$n_divergent)) {
        d$n_divergent
      } else {
        s$n_divergent
      },
      min_ess_tail = s$min_ess_tail,
      min_ebfmi = s$min_ebfmi,
      n_max_treedepth = s$n_max_treedepth,
      max_treedepth = s$max_treedepth
    )
  )
}


# --- the one function every consumer calls ---------------------------- #

# The engine a fit came from, named for a verdict: "brms" or "inla" by
# class, otherwise the class vector itself, so an object of no known engine
# is reported as what it is.
.fb_conv_engine_label <- function(fit) {
  if (inherits(fit, "flexybayes_brms")) {
    return("brms")
  }
  if (inherits(fit, "flexybayes_inla")) {
    return("inla")
  }
  paste(class(fit), collapse = "/")
}

# .fb_convergence_verdict() --- the authoritative tri-state answer.
#
# summary(), the ASReml oracle runner, the live-backend tests and
# fb_select()'s fold gate all read this and nothing else. A second
# implementation anywhere is how A25 happened.
.fb_convergence_verdict <- function(fit) {
  if (inherits(fit, "flexybayes_inla")) {
    return(.fb_conv_verdict_inla(fit))
  }
  if (inherits(fit, "flexybayes_brms")) {
    return(.fb_conv_verdict_brms(fit))
  }
  list(
    engine = tryCatch(
      .fb_conv_engine_label(fit),
      error = function(e) "unknown"
    ),
    ok = NA,
    reasons = "this fit records no convergence diagnostics",
    evidence = list(
      mode_status = NA_real_, mlik = NA_real_, numerical_confirm = NA,
      max_rhat = NA_real_, min_ess_bulk = NA_real_, n_divergent = NA_integer_,
      min_ess_tail = NA_real_, min_ebfmi = NA_real_,
      n_max_treedepth = NA_integer_, max_treedepth = NA_real_
    )
  )
}
