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

  # Factor-analytic / unstructured loadings are identified only up to
  # rotation and sign, so their per-entry Rhat is meaningless and inflates
  # this count. Point the user at the identified-quantity diagnostic.
  rt <- fit$extras$parse_info$random %||% list()
  has_struct <- any(vapply(
    rt,
    function(t) (t$type %||% "") %in% c("fa_gxe", "us_gxe"),
    logical(1)
  ))
  struct_note <- if (has_struct) {
    paste0(
      " Note: this model has a factor-analytic / unstructured term whose ",
      "raw loadings are non-identified (rotation/sign), so their Rhat is ",
      "expected to be high; consult fb_structured_cov() for the Rhat of ",
      "the identified covariance."
    )
  } else {
    ""
  }

  warning(
    "flexyBayes: the sampler may not have converged -- ",
    paste(parts, collapse = "; "),
    ". Treat the posterior with caution: increase `warmup` / ",
    "`n_samples`, simplify the model, or supply a more informative ",
    "prior. Inspect the full diagnostics with summary().",
    struct_note,
    " Silence this warning via ",
    "options(flexyBayes.silence_convergence_warning = TRUE).",
    call. = FALSE
  )
  invisible(fit)
}
