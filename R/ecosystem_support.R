# Shared helpers for the downstream-ecosystem integrations
# (emmeans, marginaleffects). These reconcile the two backends' fixed-
# effect parametrisations onto a single design-matrix contract:
#
#   * greta fits carry an over-parameterised (cell-means + intercept)
#     fixed-effect basis -- coef() names include every factor level
#     (e.g. (Intercept), fa, fb, fc). The design matrix must use
#     all-levels coding to match, and emmeans handles the resulting
#     rank deficiency through a non-estimability basis.
#   * INLA fits carry a treatment-contrast basis (e.g. (Intercept),
#     fb, fc) read from summary.fixed -- full rank.
#
# Rather than hard-code per backend, the design matrix is built to match
# the names of the fit's own coef() vector, detecting per factor whether
# the reference level is present (all-levels coding) or absent
# (treatment coding). Structures whose model matrix cannot be reconciled
# with the coefficient names (e.g. factor interactions on the over-
# parameterised greta basis) are refused rather than silently mis-mapped.

# Data frame a fit was trained on (greta keeps it on the glm shim;
# INLA keeps it at the top level).
.fb_fit_data <- function(object) {
  object$data %||% object$glm$data
}

# Per-factor `xlev` (level vocabulary) from the fit data, so a reference
# grid is coded against the fit-time factor levels.
.fb_xlev <- function(trms, data) {
  vars <- all.vars(trms)
  fac <- vars[vapply(
    vars,
    function(v) {
      is.factor(data[[v]]) || is.character(data[[v]])
    },
    logical(1L)
  )]
  stats::setNames(lapply(fac, function(v) levels(as.factor(data[[v]]))), fac)
}

# Decide, per factor in `trms`, whether the fit's coefficient names use
# all-levels coding (reference-level column present) or treatment coding
# (reference dropped). Returns a `contrasts.arg` list naming only the
# factors that need all-levels coding; treatment-coded factors are left
# to model.matrix()'s default.
.fb_detect_contrasts <- function(trms, target_names, data) {
  vars <- all.vars(trms)
  contr <- list()
  for (v in vars) {
    col <- data[[v]]
    if (is.null(col) || !(is.factor(col) || is.character(col))) {
      next
    }
    lev <- levels(as.factor(col))
    ref_col <- paste0(v, lev[1L])
    if (ref_col %in% target_names) {
      # All-levels (cell-means) coding. Use a level-labelled identity so
      # model.matrix() names the columns "<v><level>" (e.g. fa, fb, fc)
      # to match the fit's coefficient names -- contr.treatment(n,
      # contrasts = FALSE) would label them by position (f1, f2, f3).
      cm <- diag(length(lev))
      dimnames(cm) <- list(lev, lev)
      contr[[v]] <- cm
    }
  }
  contr
}

# Build the fixed-effect design matrix on `grid`, coded to match the
# coefficient names `target_names`, then align column order to them.
# Refuses (structured message) when the reconciliation is not exact.
.fb_fixef_model_matrix <- function(trms, grid, target_names, data) {
  contr <- .fb_detect_contrasts(trms, target_names, data)
  xlev <- .fb_xlev(trms, data)
  X <- tryCatch(
    stats::model.matrix(
      trms,
      data = grid,
      contrasts.arg = if (length(contr)) contr else NULL,
      xlev = xlev
    ),
    error = function(e) {
      stop(
        "Could not build the fixed-effect design matrix for this ",
        "fit: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  if (!setequal(colnames(X), target_names)) {
    stop(
      "The downstream-ecosystem integration (emmeans / ",
      "marginaleffects) could not reconcile this model's fixed-effect ",
      "design matrix (columns: ",
      paste(colnames(X), collapse = ", "),
      ") with the fitted coefficients (",
      paste(target_names, collapse = ", "),
      "). This happens for fixed-effect structures ",
      "whose model-matrix columns do not line up one-to-one with the ",
      "fitted coefficient basis. Work with the posterior draws via ",
      "fb_as_draws_simple() for this model.",
      call. = FALSE
    )
  }
  X[, target_names, drop = FALSE]
}

# Fixed-effect terms (response deleted) for either backend, from the
# captured fixed-effect formula.
.fb_fixef_terms <- function(object) {
  stats::delete.response(stats::terms(stats::formula(object)))
}
