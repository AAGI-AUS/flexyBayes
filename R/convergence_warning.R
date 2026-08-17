# Fit-time convergence warning.
#
# The print and summary methods show an Rhat / effective-size badge, but
# a user who assigns a fit and moves straight to downstream work never
# sees it. The package's standing rule is to surface an anomaly the
# moment it arises rather than wait to be asked, so a fit whose sampler
# has not mixed emits a warning() as it is handed back. The warning reads
# the same Gelman psrf + effective-size diagnostics the badge uses, so
# the two surfaces never disagree.
#
# It is a no-op when no MCMC convergence information is attached -- the
# INLA / Laplace path is deterministic and carries no Rhat, and code /
# plan objects are not fits -- and it is suppressible (for intentionally
# short illustrative fits) via
# options(flexyBayes.silence_convergence_warning = TRUE).
.fb_warn_poor_convergence <- function(
  fit,
  rhat_threshold = 1.1,
  ess_floor = 100
) {
  if (isTRUE(getOption("flexyBayes.silence_convergence_warning", FALSE))) {
    return(invisible(fit))
  }
  if (!inherits(fit, "flexybayes")) {
    return(invisible(fit))
  }
  conv <- fit$extras$convergence
  if (is.null(conv)) {
    return(invisible(fit))
  }

  rhat <- tryCatch(
    conv$gelman$psrf[, "Point est."],
    error = function(e) numeric(0)
  )
  rhat <- rhat[is.finite(rhat)]
  max_rhat <- if (length(rhat)) max(rhat) else NA_real_
  n_over <- if (length(rhat)) sum(rhat >= rhat_threshold) else 0L

  ess <- conv$n_eff
  ess <- ess[is.finite(ess)]
  min_ess <- if (length(ess)) min(ess) else NA_real_
  low_ess <- is.finite(min_ess) && min_ess < ess_floor

  if (n_over == 0L && !low_ess) {
    return(invisible(fit))
  }

  parts <- character(0)
  if (n_over > 0L) {
    parts <- c(parts, sprintf(
      "%d parameter%s with Rhat >= %.2f (max %.2f)",
      n_over,
      if (n_over == 1L) "" else "s",
      rhat_threshold,
      max_rhat
    ))
  }
  if (low_ess) {
    parts <- c(parts, sprintf("min effective sample size %.0f", min_ess))
  }

  # Two structured-term notes, and they are different diagnoses.
  #
  # A factor-analytic term's raw loadings are identified only up to
  # rotation and sign, so their per-entry Rhat is meaningless and
  # inflates this count -- fb_structured_cov() reconstructs the
  # identified covariance and is the right place to look.
  #
  # An unstructured term is NOT that case, and saying so was worse than
  # saying nothing. brms parameterises us(f):g by standard deviations
  # and a Cholesky correlation factor, which is identified; the note
  # invited a user to dismiss a real mixing failure as a labelling
  # artefact, and pointed at a function that abstains for us() terms.
  # What actually fails is the split between the covariance diagonal and
  # the residual, and it fails when there is one observation per cell:
  # measured on agridat::yan.winterwheat (18 genotypes x 9 environments,
  # one plot per cell) the 45 covariance parameters converge while
  # `sigma` reaches Rhat 1.13 with a bulk effective size of 36, and
  # doubling the budget makes it worse rather than better.
  rt <- fit$extras$parse_info$random %||% list()
  has_type <- function(types) {
    any(vapply(rt, function(t) (t$type %||% "") %in% types, logical(1)))
  }
  struct_note <- ""
  if (has_type("fa_gxe")) {
    struct_note <- paste0(
      struct_note,
      " Note: this model has a factor-analytic term whose raw loadings ",
      "are non-identified (rotation/sign), so their Rhat is expected to ",
      "be high; consult fb_structured_cov() for the Rhat of the ",
      "identified covariance."
    )
  }
  if (has_type("us_gxe")) {
    struct_note <- paste0(
      struct_note,
      " Note: this model has an unstructured us() term. With one ",
      "observation per cell of the outer factor the residual variance ",
      "is confounded with the diagonal of that covariance, and the ",
      "sampler walks the ridge between them -- divergences and a low ",
      "effective size on `sigma` are the usual symptom, and more ",
      "iterations do not fix it. Check which parameters are failing in ",
      "summary(): if they are `sigma` and `lp__` rather than the ",
      "covariance entries, the remedy is replication within cell or an ",
      "informative prior on the residual, not a longer chain. Raising ",
      "`control = list(adapt_delta = 0.95)` will lower the divergence ",
      "count on this model without identifying the split, so it hides ",
      "the symptom rather than treating it."
    )
  }

  warning(
    "flexyBayes: the sampler may not have converged -- ",
    paste(parts, collapse = "; "),
    ". Treat the posterior with caution: increase `warmup` / ",
    "`n_samples`, raise `control = list(adapt_delta = 0.95)` where the ",
    "run reported divergent transitions, simplify the model, or supply ",
    "a more informative prior. Inspect the full diagnostics with ",
    "summary().",
    struct_note,
    " Silence this warning via ",
    "options(flexyBayes.silence_convergence_warning = TRUE).",
    call. = FALSE
  )
  invisible(fit)
}
