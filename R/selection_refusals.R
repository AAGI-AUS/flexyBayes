# Selection refusals --- registry population for the selection module
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md section 10.2.
#
# Codes are registered as their throwing sites land, never in advance. A
# registered code with no triggering test is a defect by the same standard
# the module applies to everything else: a refusal nobody has seen fire is
# not a refusal, it is a string.

.populate_refusal_registry_selection <- function() {
  .register_selection_refusal(
    reason_code = "selection_covariance_not_buildable",
    description = paste0(
      "fb_select(): no admitted covariance builder applies to the base ",
      "fit, or the fit does not report what a builder needs."
    ),
    message_template = "fb_select(): %s",
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_covariance_not_positive_definite",
    description = paste0(
      "fb_select(): a covariance block has no Cholesky factor. No jitter ",
      "is applied, because a silent repair here is misspecification."
    ),
    message_template =
      "fb_select(): block %d of the covariance is not positive definite.",
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_dense_reference_too_large",
    description = paste0(
      "fb_select(): the dense covariance reference was asked for above its ",
      "size cap. It is an oracle for small cases, not a production path."
    ),
    message_template = paste0(
      "fb_select(): the dense covariance reference is %d by %d, above the cap ",
      "of %d."
    ),
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_dense_covariance_too_large",
    description = paste0(
      "fb_select(): a fold covariance would be formed as a dense n x n ",
      "matrix above the bound set by ",
      "options(flexyBayes.selection_dense_gb = ): the REML plug-in forms it ",
      "densely, and a base-model covariance that is not block-diagonal by ",
      "any grouping of its random terms cannot be factored block by block."
    ),
    message_template = paste0(
      "fb_select(): the fold covariance over these %d rows would be formed as ",
      "a dense %d x %d matrix of %s GiB, above the bound of %s GiB."
    ),
    registered_in_adr = "PLAN_AND_UNIT_BRIEFS_0.14.0_2026-09-19 (U-COV-MODEL)",
    plan_field = "selection_plan",
    since_version = "0.14.0"
  )

  .register_selection_refusal(
    reason_code = "selection_reduction_term_not_in_fit",
    description = paste0(
      "fb_select(): the random term the reduction was asked for is not in ",
      "the base fit, or its levels cannot be recovered from the data."
    ),
    message_template = "fb_select(): the base fit carries no random term `%s`.",
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_reduction_level_mismatch",
    description = paste0(
      "fb_select(): a random term reports a level that is absent from the ",
      "data, so effects and covariates cannot be matched by name."
    ),
    message_template = paste0(
      "fb_select(): random term `%s` reports level(s) absent from the data: ",
      "%s."
    ),
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_candidates_not_constant_within_term",
    description = paste0(
      "fb_select(): a candidate is not admissible for a term's reduction -- ",
      "it is not constant within the term's levels, or it does not vary ",
      "across them."
    ),
    message_template =
      "fb_select(): these candidates cannot enter the reduction for `%s`: %s.",
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_base_shape_unreadable",
    description = paste0(
      "fb_select(): the base fit does not carry the parsed model and data ",
      "the selection builders read, as on the aggregated route."
    ),
    message_template = paste0(
      "fb_select(): this base fit does not carry the parsed model and data ",
      "the reduction reads."
    ),
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_reduction_absorption_unsupported",
    description = paste0(
      "fb_select(): absorbing the reduction term as fixed leaves further ",
      "random terms, which needs a mixed refit this release does not do."
    ),
    message_template = paste0(
      "fb_select(): absorbing `%s` as a fixed term leaves %d further random ",
      "term(s)."
    ),
    registered_in_adr = "SPEC v3 section 4.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  # `selection_penalty_not_implemented` was registered at 0.10.0.9000, when
  # `penalty = NULL` was the only implemented corner. It is withdrawn at
  # 0.11.0 because the penalty is implemented and nothing throws it any
  # more: a registered code with no throwing site is a string, and the
  # registry is the count that public documents defer to.

  .register_selection_refusal(
    reason_code = "selection_base_not_allowlisted",
    description = paste0(
      "fb_select(): the base model's structure is not on the selection ",
      "allowlist, or the object is not a fitted flexyBayes model."
    ),
    message_template = paste0(
      "fb_select(): the base carries %s, which is not on the selection ",
      "allowlist."
    ),
    registered_in_adr = "SPEC v3 section 10.1",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_selection_refusal(
    reason_code = "selection_family_not_supported",
    description = paste0(
      "fb_select(): the whitening identity is exact only for a Gaussian ",
      "response with an identity link; every other family refuses."
    ),
    message_template = paste0(
      "fb_select(): the whitening identity is exact only for a Gaussian ",
      "response with an identity link, and this fit is %s/%s."
    ),
    registered_in_adr = "SPEC v3 section 5.1",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  invisible(NULL)
}

# The loss-surface figure's refusals, raised in R/selection_loss_surface.R.
.populate_refusal_registry_selection_figures <- function() {
  # A suggested-package guard: the loss surface is a three-dimensional
  # plotly display, and plotly is a suggested package, not an import.
  .register_selection_refusal(
    reason_code = "loss_surface_requires_plotly",
    description = paste0(
      "plot(type = \"loss_surface\") on an fb_select result draws a ",
      "three-dimensional surface with plotly, a suggested package that is ",
      "not installed."
    ),
    message_template = paste0(
      "plot(type = \"loss_surface\") draws a three-dimensional surface with ",
      "plotly, which is not installed. Install it with ",
      "install.packages(\"plotly\")."
    ),
    registered_in_adr = "0.14.0 figures",
    plan_field = NA_character_,
    since_version = "0.14.0"
  )

  # One fb_select() call walks the path on s at one elastic-net mixing
  # value, so the loss surface over s and alpha is drawn from several
  # results. They are drawn together only when their losses were scored on
  # one score contract, on one grid and chosen by one rule; otherwise two
  # rows of the surface would not be the same quantity.
  .register_selection_refusal(
    reason_code = "selection_loss_surface_undefined",
    description = paste0(
      "plot(type = \"loss_surface\") was given optimiser results that do ",
      "not span a surface: fewer than two elastic-net mixing values, a ",
      "mixing value given twice, an entry with no path on ",
      "s = lambda / lambda_max, or results that differ in score contract, ",
      "grid on s, decision rule, penalised columns or covariance source."
    ),
    message_template = paste0(
      "plot(type = \"loss_surface\") draws the out-of-fold loss of ",
      "optimiser results at two or more elastic-net mixing values, scored ",
      "on one score contract: %s."
    ),
    registered_in_adr = "0.14.0 figures",
    plan_field = NA_character_,
    since_version = "0.14.0"
  )

  invisible(NULL)
}

# The selection layer's refusal codes form the "flexySelect" partition of the
# registry. They are registered through the exported fb_register_refusals(),
# in one table, the call the flexySelect package makes from its .onLoad().
.populate_refusal_registry_flexyselect <- function() {
  .selection_refusal_rows$rows <- list()
  on.exit(.selection_refusal_rows$rows <- NULL)

  .populate_refusal_registry_selection()
  .populate_refusal_registry_selection_design()
  .populate_refusal_registry_selection_weights()
  .populate_refusal_registry_selection_resamples()
  .populate_refusal_registry_penalty()
  .populate_refusal_registry_selection_projection()
  .populate_refusal_registry_selection_bayes_cv()
  .populate_refusal_registry_selection_glmnet()
  .populate_refusal_registry_selection_verb()
  .populate_refusal_registry_selection_parallel()
  .populate_refusal_registry_selection_api()
  .populate_refusal_registry_selection_figures()

  rows <- .selection_refusal_rows$rows
  column <- function(name) vapply(rows, `[[`, character(1L), name)
  fb_register_refusals(
    data.frame(
      reason_code = column("reason_code"),
      description = column("description"),
      message_template = column("message_template"),
      registered_in_adr = column("registered_in_adr"),
      plan_field = column("plan_field"),
      since_version = column("since_version"),
      stringsAsFactors = FALSE
    ),
    package = "flexySelect"
  )
  invisible(NULL)
}

.selection_refusal_rows <- new.env(parent = emptyenv())

# Collects one selection code for .populate_refusal_registry_flexyselect().
.register_selection_refusal <- function(
  reason_code,
  description,
  message_template,
  registered_in_adr,
  plan_field = NA_character_,
  since_version
) {
  rows <- .selection_refusal_rows$rows
  rows[[length(rows) + 1L]] <- list(
    reason_code = reason_code,
    description = description,
    message_template = message_template,
    registered_in_adr = registered_in_adr,
    plan_field = plan_field,
    since_version = since_version
  )
  .selection_refusal_rows$rows <- rows
  invisible(NULL)
}
