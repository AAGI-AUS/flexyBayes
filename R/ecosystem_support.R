# Shared helpers for the downstream-ecosystem integrations
# (emmeans, marginaleffects). These reconcile the active backends' fixed-
# effect parametrisations onto a single design-matrix contract:
#
#   * brms and INLA fits both carry a treatment-contrast basis
#     (e.g. (Intercept), fb, fc) -- full rank. A since-withdrawn native
#     engine (see NEWS.md, 0.9.3) carried an over-parameterised
#     (cell-means + intercept) fixed-effect basis instead -- coef()
#     names included every factor level (e.g. (Intercept), fa, fb, fc),
#     needing all-levels coding to match, with emmeans handling the
#     resulting rank deficiency through a non-estimability basis.
#
# Rather than hard-code per backend, the design matrix is built to match
# the names of the fit's own coef() vector, detecting per factor whether
# the reference level is present (all-levels coding) or absent
# (treatment coding). This keeps the reconciliation general -- both
# active engines are full rank today, but the per-factor detection is a
# defensive contract, not an assumption baked in for two backends.
# Structures whose model matrix cannot be reconciled with the
# coefficient names (e.g. factor interactions on an over-parameterised
# basis) are refused rather than silently mis-mapped.

# Data frame a fit was trained on (brms keeps it on the glm shim;
# INLA keeps it at the top level).
#
# C4/FS-26: an INLA fit's `$data` is the level-legalised copy
# emit_inla() built the fit against (`$level_labels` carries the map).
# coef.flexybayes_inla() restores the user's own coefficient labels, so
# every consumer here (recover_data()/emm_basis() for predict(classify
# = ), marginaleffects' get_data()/get_predict(), plot(), the missing-
# cell summary) needs the SAME labels to build a design matrix or
# reference grid that lines up with those coefficient names --
# .fb_fixef_model_matrix()'s own reconciliation check below exists to
# catch exactly that kind of mismatch. .inla_restore_data_original_
# levels() is a no-op for any fit `$level_labels` never touched
# (every non-INLA engine, and an INLA fit with no non-syntactic level).
.fb_fit_data <- function(object) {
  data <- object$data %||% object$glm$data
  .inla_restore_data_original_levels(object$level_labels, data)
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
