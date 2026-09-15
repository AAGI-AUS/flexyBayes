# Selection resampling refusals --- registry population for the fold,
# prediction and scoring machinery.
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md sections 3.3, 4.1, 4.2,
# 5.4, 9 and 10.2.
#
# Codes register as their throwing sites land, never in advance. A registered
# code with no triggering test is a defect: a refusal nobody has seen fire is
# not a refusal, it is a string. Each entry below is thrown by
# R/selection_resamples.R and triggered by a test in
# tests/testthat/test-selection-resamples.R.
#
# Two of the twelve codes in section 10.2 belong to this tranche --
# `selection_seed_required` and `selection_fold_level_absent`. The remaining
# eight are the shape refusals the fold, prediction and scoring surfaces need
# and the contract left to the implementation, each named for the thing the
# caller got wrong rather than for the function that noticed.

.populate_refusal_registry_selection_resamples <- function() {
  .register_refusal(
    reason_code = "selection_seed_required",
    description = paste0(
      "fb_select(): a stochastic resampling step was asked for without a ",
      "seed. No seed is invented, because an invented seed is an ",
      "unreproducible fold assignment that looks reproducible."
    ),
    message_template = "fb_select(): `repeats = %d` randomises the fold assignment and no `seed` was given.",
    registered_in_adr = "SPEC v3 section 3.1",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_fold_level_absent",
    description = paste0(
      "fb_select(): an assessment fold carries a factor level that is ",
      "absent from its training rows, so the level's coefficient is not ",
      "estimable and the fold cannot be scored."
    ),
    message_template = "fb_select(): fold %s carries level(s) of `%s` absent from training: %s.",
    registered_in_adr = "SPEC v3 section 3.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_resamples_underspecified",
    description = paste0(
      "fb_resamples(): a grouping unit was named without a resampling ",
      "scheme. A unit alone is not a specification -- it says what may not ",
      "be split, not how the splitting is done."
    ),
    message_template = "fb_resamples(): `unit` was given without a scheme; supply `v`, `repeats` or `fold_id`.",
    registered_in_adr = "SPEC v3 section 3.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_fold_id_invalid",
    description = paste0(
      "fb_resamples(): the explicit `fold_id` does not describe a usable ",
      "partition -- wrong length, missing values, fewer than two folds, or ",
      "a fold boundary drawn through the middle of a grouping unit."
    ),
    message_template = "fb_resamples(): the explicit `fold_id` is not a usable partition: %s.",
    registered_in_adr = "SPEC v3 section 3.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_fold_variable_not_in_data",
    description = paste0(
      "fb_select(): a variable named by `unit`, `strata` or the scoring ",
      "cluster is not a column of the data the base model was fitted to."
    ),
    message_template = "fb_select(): `%s` names column(s) absent from the data: %s.",
    registered_in_adr = "SPEC v3 section 3.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_strata_not_constant_within_unit",
    description = paste0(
      "fb_select(): a stratifying variable varies inside a grouping unit, ",
      "so the unit belongs to two strata at once and cannot be dealt into ",
      "a fold without splitting it."
    ),
    message_template = "fb_select(): `%s` is not constant within unit(s): %s.",
    registered_in_adr = "SPEC v3 section 3.3",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_fold_not_usable",
    description = paste0(
      "fb_select(): a fold has no training rows, no assessment rows, or ",
      "fewer training rows than its design can estimate, so nothing can be ",
      "learned on it or nothing can be scored on it."
    ),
    message_template = "fb_select(): fold %s has %d training and %d assessment rows.",
    registered_in_adr = "SPEC v3 section 4.2",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_fold_covariance_not_factorisable",
    description = paste0(
      "fb_select(): the training-block covariance of a fold has no Cholesky ",
      "factor. No jitter is applied: a repair here is silent ",
      "misspecification carried into every score the fold contributes."
    ),
    message_template = "fb_select(): the training covariance of fold %s is not positive definite.",
    registered_in_adr = "SPEC v3 section 5.4",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  .register_refusal(
    reason_code = "selection_leakage_detected",
    description = paste0(
      "fb_select(): a quantity the substrate declares as learned on the ",
      "training rows moved when only the assessment rows changed, so it was ",
      "not learned on the training rows."
    ),
    message_template = "fb_select(): fold %s learned these quantities from its assessment rows: %s.",
    registered_in_adr = "SPEC v3 section 4.1",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  # Thrown by R/selection_resamples.R on the benchmark substrate and by
  # R/selection_bayes_cv.R and R/fb_select.R on the Bayesian route; triggered
  # in tests/testthat/test-selection-fold-isolation.R.
  .register_refusal(
    reason_code = "selection_formula_transform_not_fold_local",
    description = paste0(
      "fb_select(): a formula variable of the base or the candidates takes ",
      "a value on a row that depends on other rows, so a fold's training ",
      "rows would read its assessment rows. The benchmark route learns ",
      "scale(), poly(), ",
      "splines::ns() and splines::bs() on each fold's training rows and ",
      "refuses any other such transform; the Bayesian route, whose engines ",
      "evaluate the formula over the rows they are handed, refuses every ",
      "one. The refusal names the variables before the fold is fitted."
    ),
    message_template = "fb_select(): on fold %s these formula variables read other rows: %s.",
    registered_in_adr = "SPEC v3 section 4.1",
    plan_field = "selection_plan",
    since_version = "0.13.0"
  )

  .register_refusal(
    reason_code = "selection_scoring_unit_undeclared",
    description = paste0(
      "fb_select(): the standard-error unit for the reported score cannot ",
      "be formed -- no independent cluster was declared, or fewer than two ",
      "of them were realised. One-standard-error selection is not defined ",
      "without it."
    ),
    message_template = "fb_select(): the scoring unit `%s` yields %d independent unit(s); at least two are needed.",
    registered_in_adr = "SPEC v3 section 9",
    plan_field = "selection_plan",
    since_version = "0.10.0.9000"
  )

  invisible(NULL)
}
