# fb_select() --- the public verb
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 2, 3.1, 5, 6, 8
# and 9.
#
# The partition principle, which is why this function exists in this shape.
# The base model, the extended base model and a penalised selection are ONE
# computation under three partitions of the design:
#
#   base model M0      everything protected, no candidates, no penalty
#   extended base M1   base plus named candidates, still no penalty
#   selection          base protected, candidate block penalised
#
# Everything expensive is shared -- the folds, the prediction estimand, the
# scoring unit and the certificate. What differs is where the vertical line
# falls in the design matrix.
#
# Four axes, not one `engine` argument. A result records what was fitted
# (`workflow`), what fitted it (`provider`), under what prior (`prior`) and
# how the answer was read off (`decision`), because conflating them is how a
# penalised-likelihood mode comes to be reported under a Bayesian label.
#
#   workflow   bayes | penalized_optimizer | unpenalized_benchmark
#   provider   brms | inla | glmnet | gls_reml_plugin
#   prior      the joint prior IR, or "none"
#   decision   full | projection | posterior_threshold | optimizer_support
#
# Three things bind every surface this function writes to:
#
#   1. Posterior draws under a continuous shrinkage prior contain no exact
#      zeros. "Selected" is always a named decision, never a property of the
#      draws.
#   2. Every cross-provider comparison declares one of three equivalence
#      levels. The words "same as glmnet" never appear unqualified.
#   3. Every certificate field that names an estimand, a provider, a
#      convention or a route is copied from the object that did the work,
#      never from an argument. The public argument that used to name the
#      covariance estimand said `plugin_mean_covariance` while the folds used
#      `reml_plugin_covariance`; the field is now read from the fold.


# --- the gate -------------------------------------------------------- #

# .fb_select_gate() --- the axes answered before anything runs.
#
# A base is admitted only when every axis passes, and the refusal states the
# measurement behind the conclusion. `us()` stays off the allowlist
# whichever engine fitted the base: INLA already refuses an estimated
# correlation matrix at the base fit itself
# (`correlated_covariance_requires_brms`, R/lgm_gate.R), on the same
# measurement named below, so the only base that can reach here carrying one
# is a brms fit, and no fold-refit test shows the selection substrate's
# whitening rotation reproduces it. `fa()` fits on INLA only and is marked
# experimental there; the same absence of a fold-refit test keeps it off the
# allowlist too. A correlated random slope on INLA stays out on the same
# measurement taken on a second route, and only on that engine: the brms
# base is admitted.
.fb_select_gate <- function(base) {
  if (!inherits(base, "flexybayes")) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_not_allowlisted",
      message = "fb_select(): `base` must be a fitted flexyBayes model."
    ))
  }
  if (is.null(base$fb) || is.null(base$data)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_shape_unreadable",
      message = paste0(
        "fb_select(): this base fit does not carry the parsed model and the ",
        "data the selection substrate reads. Fits on the aggregated route, ",
        "which an intercept-only fixed model takes, have a different shape."
      )
    ))
  }

  fam <- base$fb$family %||% NA_character_
  link <- base$fb$link %||% NA_character_
  if (!identical(fam, "gaussian") || !identical(link, "identity")) {
    stop(.fb_refusal_condition(
      reason_code = "selection_family_not_supported",
      message = paste0(
        "fb_select(): the whitening identity is exact only for a Gaussian ",
        "response with an identity link, and this fit is ", fam, "/", link,
        ". For other families the observation variance depends on the mean, ",
        "no single V survives integrating the random effects, and an IRLS ",
        "formulation changes the weights as the coefficients change. A ",
        "penalised-likelihood route is a separate programme with its own ",
        "validation."
      )
    ))
  }

  # Structures whose covariance is not admitted, by name and with the
  # measurement that put them there. Read from .FB_GXE_PAIR_TYPE
  # (R/parse_formula.R), the same table .classify_interaction() writes
  # "us_gxe" and "fa_gxe" from: one source for both readers, so a name added
  # there is caught here automatically.
  types <- vapply(
    base$fb$random_terms %||% list(),
    function(t) t$type %||% NA_character_, character(1L)
  )
  us_bad <- intersect(types, c("us", .FB_GXE_PAIR_TYPE[["us"]]))
  if (length(us_bad) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_not_allowlisted",
      message = paste0(
        "fb_select(): the base carries ", paste(us_bad, collapse = ", "),
        ", which is not on the selection allowlist. An estimated ",
        "correlation matrix fits on brms only -- INLA already refuses it at ",
        "the base fit, citing a genetic correlation of 0.142 against ",
        "ASReml's 0.930 on agridat::besag.met, where rescaling the response ",
        "by ten does not reproduce it -- and no fold-refit test shows the ",
        "selection substrate's whitening rotation reproduces a brms fit's ",
        "covariance for this structure, so it stays off the allowlist on ",
        "either route."
      )
    ))
  }
  fa_bad <- intersect(types, c("fa", .FB_GXE_PAIR_TYPE[["fa"]]))
  if (length(fa_bad) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_not_allowlisted",
      message = paste0(
        "fb_select(): the base carries ", paste(fa_bad, collapse = ", "),
        ", which is not on the selection allowlist. A factor-analytic ",
        "covariance fits on INLA only, and the package's own README already ",
        "marks it experimental there: repeat fits of one model have ",
        "disagreed. The selection substrate refits the base on every fold, ",
        "and no test shows its whitening rotation is stable for this ",
        "structure, so it stays off the allowlist."
      )
    ))
  }

  # A correlated random slope carries the same defect on a second route, and
  # it is engine-specific rather than structural: brms writes `(x | g)`
  # natively and stays admitted, while the INLA lowering does not reproduce
  # the covariance a rotation would be built on.
  if (
    inherits(base, "flexybayes_inla") && "simple_slope_cor" %in% types
  ) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_not_allowlisted",
      message = paste0(
        "fb_select(): the base carries a correlated random slope `(x | g)` ",
        "fitted on INLA, which is not on the selection allowlist. On ",
        "lme4::sleepstudy that route returns a slope variance of 46.4 at ",
        "the posterior mode against 35.08 in ASReml and 35.07 in lme4, and ",
        "multiplying the response by ten -- an exact transformation under ",
        "which a correlation is invariant -- moves the intercept-slope ",
        "correlation from 0.02 to 0.969 while REML holds at 0.066. The ",
        "selection substrate whitens by the base's own covariance, so that ",
        "block would set the rotation every score is read through. Refit ",
        "the base with backend = \"brms\", which writes this term natively ",
        "and is admitted here."
      )
    ))
  }

  invisible(TRUE)
}


# --- which provider, and is it here ------------------------------------- #

# .fb_select_engine_of() --- the engine that fitted the base, from its class.
#
# The class of the fitted object is the oracle: `call_info$backend` holds the
# routing POLICY the call asked for ("auto"), which is not the engine the
# request landed on.
.fb_select_engine_of <- function(base) {
  if (inherits(base, "flexybayes_brms")) {
    return("brms")
  }
  if (inherits(base, "flexybayes_inla")) {
    return("inla")
  }
  stop(.fb_refusal_condition(
    reason_code = "selection_base_not_allowlisted",
    message = paste0(
      "fb_select(): `base` must be a fitted flexybayes_brms or ",
      "flexybayes_inla model; this one is <",
      paste(class(base), collapse = ", "),
      ">."
    )
  ))
}

# .fb_select_required_providers() --- the packages a request needs, by name.
#
# The Bayesian route with a joint prior is brms and nothing else: the
# reference model is re-emitted on brms whatever engine fitted the base, and
# no public call hands a joint shrinkage prior to INLA. `projpred` joins
# the set only when the projection decision was asked for, so a
# `decision = "full"` run on a machine without it still works.
.fb_select_required_providers <- function(engine, penalty, decision) {
  if (is.null(penalty)) {
    return(if (identical(engine, "brms")) "brms" else "INLA")
  }
  packages <- "brms"
  if (identical(decision$method %||% "projection", "projection")) {
    packages <- c(packages, "projpred")
  }
  packages
}

.fb_select_missing_providers <- function(packages) {
  packages[!vapply(packages, .fb_sel_provider_installed, logical(1L))]
}

.fb_select_install_line <- function(package) {
  if (identical(package, "INLA")) {
    return(paste0(
      "install.packages(\"INLA\", repos = c(getOption(\"repos\"), ",
      "INLA = \"https://inla.r-inla-download.org/R/stable\"))"
    ))
  }
  paste0("install.packages(\"", package, "\")")
}


# --- the decision, resolved against the structure gate --------------------- #

# .fb_select_resolve_decision() --- turn `auto` into the decision that can
# actually run, and refuse an explicit one that cannot.
#
# The gate is the section 6.3 matrix evaluated for THIS model. A row that has
# not been executed refuses exactly as a row that was executed and came back
# red, so "green" here means every feature the model carries has a run behind
# it. `auto` reads that verdict; an explicit `projection` does not get to,
# and still refuses by name at the projection itself with the feature that
# failed.
.fb_select_resolve_decision <- function(decision, structure_gate) {
  requested <- decision$requested_method %||% decision$method
  if (!identical(decision$method, "auto")) {
    decision$requested_method <- requested
    return(list(
      decision = decision,
      resolved_from = NA_character_,
      gate_verdict = .fb_select_gate_verdict(structure_gate),
      reason = NA_character_
    ))
  }
  verdict <- .fb_select_gate_verdict(structure_gate)
  green <- identical(verdict, "supported")
  decision$method <- if (green) "projection" else "full"
  decision$requested_method <- "auto"
  reason <- if (is.na(verdict)) {
    paste0(
      "the projection structure gate could not be read for this base, so ",
      "the decision falls to `full`, which needs no projectable structure"
    )
  } else if (green) {
    paste0(
      "every structural feature of the reference model has an executed ",
      "projection gate behind it, so the decision is the projection"
    )
  } else {
    paste0(
      "the projection structure gate is red (", verdict,
      "), so the decision is `full`: the reference posterior under the ",
      "shrinkage prior, with no selection. Naming decision = ",
      "\"projection\" explicitly refuses instead, with the feature that ",
      "failed"
    )
  }
  list(
    decision = decision,
    resolved_from = "auto",
    gate_verdict = verdict,
    reason = reason
  )
}


# The one-line reading of the gate table: "supported", or the features that
# are not, or NA where no gate could be built.
.fb_select_gate_verdict <- function(structure_gate) {
  if (!is.data.frame(structure_gate) || nrow(structure_gate) == 0L) {
    return(NA_character_)
  }
  bad <- structure_gate$feature[structure_gate$verdict != "supported"]
  if (length(bad) == 0L) {
    return("supported")
  }
  paste0("refuses: ", paste(bad, collapse = ", "))
}


# --- the information a prior scale is calibrated against ------------------ #

# .fb_select_calibration_n() --- the n a prior scale is stated against.
#
# The Piironen and Vehtari calibration of a global-local prior reads
# `tau0 = (p0 / (p - p0)) * sigma / sqrt(n)`, and the elastic net's
# canonical objective carries `n` in front of its penalty. Both were given
# the ROW count. In a clustered design that overstates the information about
# a coefficient by roughly the cluster size -- a hundred plots in ten sites
# are not a hundred independent pieces of evidence about a site-level
# covariate -- and the prior over-shrinks by a factor nobody sees.
#
# The number of resampling units is the calibration this module can state
# from what the caller declared. It is not the effective sample size, and it
# does not claim to be: it is the count of the blocks the folds treat as
# independent, which is the same claim the standard errors are aggregated
# under. Where no unit was declared, the rows ARE the units and the row
# count is correct rather than a fallback -- and the basis says which.
.fb_select_calibration_n <- function(resamples, data) {
  unit <- if (is.null(resamples)) character(0) else resamples$unit
  unit <- as.character(unit %||% character(0))
  if (length(unit) == 0L || !all(unit %in% names(data))) {
    return(list(
      n = nrow(data),
      basis = "rows",
      unit = character(0),
      n_rows = nrow(data)
    ))
  }
  list(
    n = length(unique(.fb_sel_key(data, unit))),
    basis = "resampling_units",
    unit = unit,
    n_rows = nrow(data)
  )
}


# .fb_select_calibration_normalise() --- accept either the resolved record
# or a bare count, so a caller building a certificate by hand (which the
# falsifiers do) is not obliged to construct the whole list.
.fb_select_calibration_normalise <- function(calibration, base) {
  n_rows <- nrow(base$data %||% data.frame())
  if (is.null(calibration)) {
    return(list(
      n = n_rows, basis = "rows", unit = character(0), n_rows = n_rows
    ))
  }
  if (is.list(calibration)) {
    return(list(
      n = calibration$n %||% n_rows,
      basis = calibration$basis %||% "rows",
      unit = as.character(calibration$unit %||% character(0)),
      n_rows = calibration$n_rows %||% n_rows
    ))
  }
  list(
    n = calibration, basis = "rows", unit = character(0), n_rows = n_rows
  )
}


# --- the penalty pipeline ------------------------------------------------ #

# .fb_select_fit_candidate_transform() --- learn the candidate transform and
# compile only the representation the requested workflow executes.
#
# It is the fit half of a fit/apply pair. The reference posterior calls it on
# the whole sample; the out-of-fold score calls it on each fold's analysis
# rows, so every constant it learns -- the centres and scales, the
# coefficient map and `calibration_n` -- and, on the Bayesian route, the
# compiled prior and `n_train` -- is a function of the rows handed to it and
# of nothing else.
#
# Four steps, in the order the contract fixes them, and each of them is a
# function a sibling unit owns: resolve the targets against the model matrix
# (`.fb_pen_resolve`), learn the scaling constants (`.fb_pen_scaling`),
# decide which of those constants the reference fit can actually apply, and
# compile the workflow's representation (`.fb_pen_compile_brms` for Bayes;
# a compiler-shaped, prior-free record for the optimiser).
#
# The third step is the one that needs saying. The scaling is declared on the
# candidate DESIGN columns, and the reference fit is a re-emission of the
# base FORMULA on a data frame -- so a design column that is a factor
# contrast or an interaction product has no data column to standardise. Those
# columns are left on their own scale, with `centre = 0` and `scale = 1`, so
# that the back-transformation is the identity on exactly the columns the fit
# did not move, and the certificate records the departure by name. Silently
# dividing the coefficient of an unstandardised column by a standard
# deviation would be a wrong answer with a right-looking provenance.
#' Resolve a penalty against the rows given
#'
#' @param data The rows every constant is learned on.
#' @param candidates,penalty,protected_labels As for `fb_select()`.
#' @param calibration A `.fb_select_calibration_n()` record for `data`.
#' @param workflow The workflow whose representation is being prepared.
#' @return The resolved blocks, the workflow representation, the column map
#'   and the constants `.fb_select_apply_candidate_transform()` reads.
#' @noRd
#' @keywords internal
.fb_select_fit_candidate_transform <- function(data,
                                               candidates,
                                               penalty,
                                               protected_labels,
                                               calibration,
                                               workflow = c(
                                                 "bayes",
                                                 "penalized_optimizer"
                                               )) {
  workflow <- match.arg(workflow)
  candidate_labels <- attr(
    stats::terms(candidates, data = data),
    "term.labels"
  )
  fixed_rhs <- c(protected_labels, candidate_labels)
  map_formula <- stats::as.formula(paste0(
    "~ ",
    if (length(fixed_rhs) > 0L) paste(fixed_rhs, collapse = " + ") else "1"
  ))
  map <- .fb_sel_term_map(map_formula, data)
  part <- .fb_sel_partition_design(map, protected = protected_labels)

  blocks <- .fb_pen_resolve(penalty, map, protected = protected_labels)
  blocks <- .fb_pen_scaling(blocks, part$Xp, train_rows = seq_len(nrow(data)))

  # Which candidate design columns ARE a plain numeric data column, and so
  # can be standardised in the data the reference fit is emitted on.
  candidate_names <- colnames(part$Xp)
  scalable <- vapply(
    candidate_names,
    function(nm) {
      nm %in% names(data) && is.numeric(data[[nm]]) &&
        nm %in% map$term_labels &&
        length(map$term_map[[nm]]) == 1L
    },
    logical(1L)
  )
  left_alone <- character(0)
  for (i in seq_along(blocks)) {
    b <- blocks[[i]]
    for (k in seq_along(b$coef_names)) {
      column <- candidate_names[b$columns[[k]]]
      if (!isTRUE(scalable[[column]])) {
        if (!identical(b$scale[[k]], 1) || !identical(b$centre[[k]], 0)) {
          left_alone <- c(left_alone, column)
        }
        b$centre[[k]] <- 0
        b$scale[[k]] <- 1
      }
    }
    blocks[[i]] <- b
  }

  # The information a prior scale is calibrated against. In clustered data
  # the row count overstates it -- a hundred plots in ten sites carry
  # nothing like a hundred independent pieces of evidence about a
  # site-level covariate -- so the calibration is the number of resampling
  # units where one was declared, and the row count only where none was.
  blocks <- .fb_pen_calibrate_n(
    blocks,
    calibration_n = calibration$n,
    n_rows = nrow(data)
  )

  # The class-wide families cannot be told which coefficient is on which
  # scale: brms puts one global-local prior on the whole `b` vector, so an
  # unstandardised 0/1 contrast column sits under the same local scale as a
  # numeric already at unit standard deviation. That refuses rather than
  # being certified as an approximation.
  class_wide_unscaled <- if (length(left_alone) == 0L) {
    character(0)
  } else {
    unique(unlist(lapply(blocks, function(b) {
      if (!identical(b$route, "class_wide")) {
        return(character(0))
      }
      intersect(candidate_names[unlist(b$columns)], left_alone)
    }), use.names = FALSE))
  }
  if (length(class_wide_unscaled) > 0L) {
    families <- paste(unique(vapply(
      Filter(function(b) identical(b$route, "class_wide"), blocks),
      function(b) b$family, character(1L)
    )), collapse = ", ")
    stop(.fb_refusal_condition(
      reason_code = "penalty_candidate_columns_unscaled",
      message = paste0(
        "fb_penalty(): ", families, " is a class-wide global-local prior -- ",
        "brms puts one of them on every population-level coefficient -- and ",
        "these candidate design columns were left on their own scale: ",
        paste(class_wide_unscaled, collapse = ", "),
        ". They are factor contrasts or interaction products with no data ",
        "column to centre and divide, so the reference fit, which is a ",
        "re-emission of the formula rather than a design-matrix fit, cannot ",
        "standardise them. A shared local-scale prior then treats a 0/1 ",
        "dummy as comparable to a numeric already at unit standard ",
        "deviation, and the shrinkage it applies is a property of the ",
        "contrast coding rather than of the science. Two ways forward: use ",
        "ridge(), laplace() or elastic_net(), which put one prior row on ",
        "each coefficient and can therefore state a scale per column; or ",
        "restrict the candidate block to numeric columns, which the ",
        "re-emission can standardise."
      ),
      columns = class_wide_unscaled,
      # Not `family`: it partially matches `family_class`, and the value
      # would become a condition class rather than a field.
      penalty_family = families
    ))
  }

  coef_names_all <- .fb_pen_brms_coef_names(
    colnames(map$X)[map$assign != 0L]
  )
  compiled <- if (identical(workflow, "bayes")) {
    .fb_pen_compile_brms(
      blocks,
      coef_names_all,
      n_train = calibration$n
    )
  } else {
    .fb_select_compile_optimizer(blocks, coef_names_all)
  }

  # The constants the apply step reads: one row per candidate column in a
  # block, and whether the apply step moves it.
  constants <- do.call(rbind, c(
    list(data.frame(
      column = character(0), block = integer(0), centre = numeric(0),
      scale = numeric(0), standardised = logical(0),
      stringsAsFactors = FALSE
    )),
    lapply(blocks, function(b) {
      columns <- candidate_names[b$columns]
      data.frame(
        column = columns,
        block = rep(b$id, length(columns)),
        centre = as.numeric(b$centre),
        scale = as.numeric(b$scale),
        standardised = vapply(
          seq_along(columns),
          function(k) {
            isTRUE(scalable[[columns[[k]]]]) && !identical(b$scale[[k]], 1)
          },
          logical(1L)
        ),
        stringsAsFactors = FALSE
      )
    })
  ))

  scope_all_fixed <- any(vapply(
    blocks,
    function(b) identical(b$hyperparameters$scope %||% "candidates",
                          "all_fixed"),
    logical(1L)
  ))

  list(
    map = map,
    partition = part,
    blocks = blocks,
    compiled = compiled,
    candidate_labels = candidate_labels,
    candidate_names = candidate_names,
    coef_names_all = coef_names_all,
    calibration = calibration,
    n_train = calibration$n,
    n_rows = nrow(data),
    constants = constants,
    unscaled_columns = unique(left_alone),
    scope_all_fixed = scope_all_fixed,
    column_map = .fb_select_column_map(part, blocks, protected_labels, map)
  )
}


#' Put a learned candidate transform on a data frame
#'
#' The apply half of the fit/apply pair: centres and scales the standardised
#' candidate columns by the stored constants and leaves every other column
#' as it was. It reads the constants and nothing about the rows it is given.
#'
#' @param data A data frame carrying the candidate columns.
#' @param transform A `.fb_select_fit_candidate_transform()` result.
#' @return `data`, with the standardised columns transformed.
#' @noRd
#' @keywords internal
.fb_select_apply_candidate_transform <- function(data, transform) {
  constants <- transform$constants
  for (i in which(constants$standardised)) {
    column <- constants$column[[i]]
    data[[column]] <-
      (data[[column]] - constants$centre[[i]]) / constants$scale[[i]]
  }
  data
}


#' Resolve a penalty for the whole-sample reference posterior
#'
#' The fit step on every row, the apply step on the same rows, and the
#' approximations the resolution commits to. They are settled before
#' anything is fitted, so `plan = TRUE` reports the list the certificate
#' will carry.
#'
#' @param base A normalised base fit carrying `$data`.
#' @param candidates,penalty,protected_labels As for `fb_select()`.
#' @param calibration A `.fb_select_calibration_n()` record, or `NULL`.
#' @param workflow The workflow whose representation is being prepared.
#' @return The fit step's result plus `reference_data` and `approximations`.
#' @noRd
#' @keywords internal
.fb_select_resolve_penalty <- function(base,
                                       candidates,
                                       penalty,
                                       protected_labels,
                                       calibration = NULL,
                                       workflow = c(
                                         "bayes",
                                         "penalized_optimizer"
                                       )) {
  workflow <- match.arg(workflow)
  data <- base$data
  calibration <- calibration %||% .fb_select_calibration_n(NULL, data)
  resolved <- .fb_select_fit_candidate_transform(
    data, candidates, penalty, protected_labels, calibration, workflow
  )
  resolved$reference_data <- .fb_select_apply_candidate_transform(
    data, resolved
  )

  approximations <- character(0)
  if (length(resolved$unscaled_columns) > 0L) {
    approximations <- c(approximations, paste0(
      "candidate_columns_not_standardised: ",
      paste(resolved$unscaled_columns, collapse = ", "),
      " are factor contrasts or interaction products with no data column ",
      "to centre and scale, so their prior scale is on their own scale"
    ))
  }
  if (!identical(calibration$basis, "rows")) {
    approximations <- c(approximations, paste0(
      "calibration_n_at_the_resampling_unit: the prior scales are ",
      "calibrated against ", calibration$n, " units of ",
      paste(calibration$unit, collapse = ", "), " rather than against ",
      nrow(data), " rows, because rows inside a cluster are not ",
      "independent pieces of evidence about a coefficient"
    ))
  }
  if (isTRUE(resolved$scope_all_fixed)) {
    approximations <- c(approximations, "protected_under_global_shrinkage")
  }
  resolved$approximations <- approximations
  resolved
}


#' Record an optimiser transform without pretending it is a brms prior
#'
#' The optimiser consumes resolved penalty blocks directly. Compiling those
#' blocks through brms is both unnecessary and wrong for an `optimizer_only`
#' family: grouped penalties are proximal objectives, not densities. This
#' record keeps the common transform shape used by plans and certificates,
#' but has no prior rows or Stan variables because neither is executed.
#'
#' @param blocks Resolved penalty blocks.
#' @param coef_names_all Population-level coefficient names in design order.
#' @return A compiler-shaped record for the optimiser route.
#' @noRd
#' @keywords internal
.fb_select_compile_optimizer <- function(blocks, coef_names_all) {
  .fb_pen_assert_targeting(blocks, coef_names_all)
  lapply(blocks, function(block) {
    .fb_sel_glmnet_assert_family(block$family)
  })

  rows <- lapply(blocks, function(block) {
    idx <- match(block$coef_names, coef_names_all)
    if (anyNA(idx)) {
      stop(
        ".fb_select_compile_optimizer(): resolved coefficient(s) are not ",
        "in the fitted design: ",
        paste(block$coef_names[is.na(idx)], collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    data.frame(
      coefficient = block$coef_names,
      b_index = as.integer(idx),
      block = rep(block$id, length(idx)),
      family = rep(block$family, length(idx)),
      route = rep("penalized_optimizer", length(idx)),
      hyperparameter_mode = rep(block$hyperparameter_mode, length(idx)),
      stringsAsFactors = FALSE
    )
  })
  param_map <- do.call(rbind, rows)
  param_map <- param_map[order(param_map$b_index), , drop = FALSE]
  rownames(param_map) <- NULL

  list(
    prior = NULL,
    stanvars = NULL,
    param_map = param_map,
    route = "penalized_optimizer"
  )
}


#' The fold-local candidate transform an out-of-fold score learns
#'
#' Returns the `fold_transform` that `.fb_sel_bayes_cv()` calls with a
#' fold's analysis rows. It re-resolves the penalty on those rows and returns
#' the fit step's result, the apply step and the record the certificate
#' carries.
#'
#' @param candidates,penalty,protected_labels,resamples As for `fb_select()`.
#' @return A function of the analysis rows.
#' @noRd
#' @keywords internal
.fb_select_fold_transform <- function(candidates,
                                      penalty,
                                      protected_labels,
                                      resamples) {
  function(analysis_data) {
    calibration <- .fb_select_calibration_n(resamples, analysis_data)
    transform <- .fb_select_fit_candidate_transform(
      analysis_data, candidates, penalty, protected_labels, calibration
    )
    list(
      transform = transform,
      apply = function(frame) {
        .fb_select_apply_candidate_transform(frame, transform)
      },
      record = .fb_select_transform_record(transform)
    )
  }
}


#' What a fold's transform learned, in the form the certificate stores
#'
#' @param transform A `.fb_select_fit_candidate_transform()` result.
#' @return A list of the constants, the calibration, `n_train`, the row
#'   count and the compiled prior rows as text.
#' @noRd
#' @keywords internal
.fb_select_transform_record <- function(transform) {
  prior <- transform$compiled$prior
  list(
    constants = transform$constants,
    calibration_n = transform$calibration$n,
    calibration_basis = transform$calibration$basis,
    n_train = transform$n_train,
    n_rows = transform$n_rows,
    coef_names = transform$coef_names_all,
    prior = if (is.null(prior)) {
      character(0)
    } else {
      paste0(
        prior$class,
        ifelse(nzchar(prior$coef), paste0(":", prior$coef), ""),
        " ~ ", prior$prior
      )
    }
  )
}


#' The fold-local preprocessing of an out-of-fold run, flattened
#'
#' @param cv A `.fb_sel_bayes_cv()` result, or any list.
#' @return `NULL` where no fold-local transform ran; otherwise
#'   `list(learned_on, constants, folds)`, the constants one row per fold
#'   and column, the folds one row per fold.
#' @noRd
#' @keywords internal
.fb_select_fold_preprocessing <- function(cv) {
  by_fold <- cv$preprocessing$by_fold
  if (!isTRUE(cv$preprocessing$fold_local_transform) ||
        length(by_fold) == 0L) {
    return(NULL)
  }
  labels <- names(by_fold)
  list(
    learned_on = cv$preprocessing$learned_on,
    constants = do.call(rbind, lapply(labels, function(label) {
      cbind(
        fold = rep(label, nrow(by_fold[[label]]$constants)),
        by_fold[[label]]$constants,
        stringsAsFactors = FALSE
      )
    })),
    folds = do.call(rbind, lapply(labels, function(label) {
      r <- by_fold[[label]]
      data.frame(
        fold = label,
        n_rows = r$n_rows,
        calibration_n = r$calibration_n,
        calibration_basis = r$calibration_basis,
        n_train = r$n_train,
        prior = paste(r$prior, collapse = "; "),
        stringsAsFactors = FALSE
      )
    }))
  )
}

# .fb_select_column_map() --- one row per design column, with its role, its
# block, the scaling constants that were actually applied to it, and the
# FORMULA TERM the column came from.
#
# The term is not decoration. A decision over a formula chooses terms (`f`,
# `x1:f`); a design carries contrast columns (`fB`, `x1:fB`). Matching a
# chosen term against a column name reports every factor and interaction
# column as dropped however the decision went, so the term travels beside
# the column and `tidy()` reads retention off it.
.fb_select_column_map <- function(part, blocks, protected_labels, map = NULL) {
  candidate_names <- colnames(part$Xp)
  block_of <- rep(NA_integer_, length(candidate_names))
  family_of <- rep(NA_character_, length(candidate_names))
  centre <- rep(0, length(candidate_names))
  scale <- rep(1, length(candidate_names))
  for (b in blocks) {
    for (k in seq_along(b$columns)) {
      j <- b$columns[[k]]
      block_of[[j]] <- b$id
      family_of[[j]] <- b$family
      centre[[j]] <- b$centre[[k]]
      scale[[j]] <- b$scale[[k]]
    }
  }
  protected_names <- colnames(part$X0)
  # The term each design column belongs to, read off the assignment map the
  # resolver already built rather than parsed back out of the column name.
  term_of <- function(columns, fallback) {
    if (is.null(map) || is.null(map$columns)) {
      return(as.character(fallback))
    }
    as.character(map$columns$term[columns])
  }
  data.frame(
    column = c(protected_names, candidate_names),
    term = c(
      term_of(part$protected_columns, protected_names),
      term_of(part$candidate_columns, candidate_names)
    ),
    role = c(
      rep("protected", length(protected_names)),
      rep("candidate", length(candidate_names))
    ),
    block = c(rep(NA_integer_, length(protected_names)), block_of),
    family = c(rep("none", length(protected_names)), family_of),
    centre = c(rep(0, length(protected_names)), centre),
    scale = c(rep(1, length(protected_names)), scale),
    stringsAsFactors = FALSE
  )
}


# --- the certificate ----------------------------------------------------- #

# How the prior's hyperparameters were treated, and how the selection reused
# the data. Ridge, Laplace and elastic-net reference fits are Bayesian
# conditional on penalty values the caller fixed; they do not learn the
# global penalty strength. The global-local families sample their scales
# with the coefficients. An optimiser tunes its lambda on the folds. A
# single "Bayesian" label cannot tell these apart, so the certificate does.

#' How one prior block's hyperparameters were treated
#'
#' @param block A resolved prior block.
#' @param workflow The resolved workflow.
#' @return `"fixed_user"`, `"calibrated_rule"`, `"tuned_external"` or
#'   `"hierarchical"`; `NA` for an explicitly unpenalised block.
#' @noRd
#' @keywords internal
.fb_select_block_treatment <- function(block, workflow) {
  if (identical(block$route %||% "", "none")) {
    return(NA_character_)
  }
  if (identical(workflow, "penalized_optimizer")) {
    return("tuned_external")
  }
  switch(
    block$hyperparameter_mode %||% "",
    fixed = "fixed_user",
    learned = "hierarchical",
    NA_character_
  )
}


#' The hyperparameter treatments of a run's prior blocks
#'
#' @param blocks The resolved prior blocks, or `NULL`.
#' @param workflow The resolved workflow.
#' @return The distinct treatments in block order; empty where no joint prior
#'   was declared.
#' @noRd
#' @keywords internal
.fb_select_hyperparameter_treatment <- function(blocks, workflow) {
  treatments <- vapply(
    blocks %||% list(),
    .fb_select_block_treatment,
    character(1L),
    workflow = workflow
  )
  unique(treatments[!is.na(treatments)])
}


#' How the selection reused the data it is scored on
#'
#' @param axes The resolved axes.
#' @param projection The projection result, or `NULL`.
#' @return `list(level, reason)`, the level one of `"none"`,
#'   `"decision_reads_reference_posterior"`, `"decision_on_scored_folds"` or
#'   `"tuning_on_scored_folds"`.
#' @noRd
#' @keywords internal
.fb_select_data_reuse <- function(axes, projection = NULL) {
  reuse <- function(level, reason) list(level = level, reason = reason)
  if (identical(axes$workflow, "penalized_optimizer")) {
    return(reuse("tuning_on_scored_folds", paste0(
      "the optimiser's tuning point, s = lambda / lambda_max with each ",
      "fold's lambda_max read on its training rows (lambda itself for ridge, ",
      "which has no finite lambda_max), is chosen on the out-of-fold curve ",
      "over these folds and the chosen model's score is read off the same ",
      "curve, so it is not a nested estimate"
    )))
  }
  if (!identical(axes$workflow, "bayes") || identical(axes$prior, "none")) {
    return(reuse(
      "none", "every candidate is protected, so nothing is chosen from the data"
    ))
  }
  switch(
    axes$decision,
    projection = reuse("decision_on_scored_folds", paste0(
      "the submodel size is chosen on the elpd curve over these folds",
      if (isTRUE(projection$settings$validate_search)) {
        ", with the search repeated inside each fold"
      } else {
        paste0(
          ", and the search was not repeated inside the folds, so that ",
          "curve is optimistic"
        )
      },
      if (identical(
        projection$settings$decision_stage_search,
        "each_fold_training_rows"
      )) {
        paste0(
          "; the decision-stage score searches the terms again on each ",
          "fold's training rows at that size"
        )
      } else {
        ""
      },
      "; a decision-stage score on the same folds is not a nested estimate"
    )),
    posterior_threshold = reuse("decision_reads_reference_posterior", paste0(
      "the retained set is a fixed rule applied to the one posterior fitted ",
      "on every row; no fold score enters the choice"
    )),
    reuse("none", paste0(
      "decision = \"full\" chooses nothing: the answer is the one reference ",
      "posterior, whose hyperparameters were declared or learned inside it"
    ))
  )
}


# .fb_select_certificate() --- what was done, and every place it was
# approximate.
#
# Contract section 9.2. This is the package's strongest existing idea --
# stated inability -- carried into selection. A certificate that recorded
# only the happy path would be decoration; the fields that earn it are
# `approximations` and `refusals`, which say where the answer is weaker than
# it looks.
#
# Every field naming an estimand, a provider, a convention or a route is read
# off the object that produced it. The `cv` argument is the fold object and
# is the only thing consulted for the covariance estimand.
.fb_select_certificate <- function(base,
                                   axes,
                                   cv,
                                   resamples,
                                   candidates,
                                   blocks = NULL,
                                   adequacy = NULL,
                                   projection = NULL,
                                   benchmark = NULL,
                                   fallback_fired = FALSE,
                                   approximations = character(0),
                                   refusals = character(0),
                                   calibration = NULL,
                                   provider_versions = list(),
                                   column_map = NULL,
                                   reference_engine = NA_character_,
                                   reference_predictive = NULL,
                                   decision_predictive = NULL,
                                   structure_gate = NULL,
                                   decision_resolution = NULL,
                                   fold_preprocessing = NULL,
                                   candidate_space = NULL) {
  conv <- .fb_convergence_verdict(base)
  calibration <- .fb_select_calibration_normalise(calibration, base)
  not_assembled <- paste0(
    "no out-of-fold comparison was assembled for this certificate."
  )
  table_metric <- .fb_select_benchmark_metric(cv$metric %||% "elpd")
  reference_predictive <- reference_predictive %||%
    .fb_select_stage_table(
      "reference", list(), table_metric, unavailable = not_assembled
    )
  decision_predictive <- decision_predictive %||%
    .fb_select_stage_table(
      "decision", list(), table_metric, unavailable = not_assembled
    )
  bs <- try(.fb_sel_block_structure(base), silent = TRUE)

  block_record <- lapply(blocks %||% list(), function(b) {
    list(
      id = b$id,
      family = b$family,
      target = b$target,
      coef_names = b$coef_names,
      hyperparameters = b$hyperparameters,
      hyperparameter_mode = b$hyperparameter_mode,
      hyperparameter_treatment = .fb_select_block_treatment(b, axes$workflow),
      route = b$route,
      scaling = b$scaling,
      centre = b$centre,
      scale = b$scale
    )
  })

  structure(
    list(
      workflow = axes$workflow,
      provider = axes$provider,
      prior = axes$prior,
      decision = axes$decision,
      # The decision the caller asked for beside the one that ran, and the
      # gate that decided between them. `auto` is the default, so the
      # resolution has to be visible rather than inferable.
      decision_requested = decision_resolution$requested %||% axes$decision,
      decision_resolved_from = decision_resolution$resolved_from %||%
        NA_character_,
      decision_resolution_reason = decision_resolution$reason %||%
        NA_character_,
      structure_gate = structure_gate,
      structure_verdict = .fb_select_gate_verdict(structure_gate),
      base = list(
        engine = base$fb$source %||% NA_character_,
        family = base$fb$family,
        link = base$fb$link,
        n_observed = nrow(base$data),
        convergence = conv$ok,
        convergence_reasons = conv$reasons
      ),
      covariance = list(
        # Read from the fold object, never from an argument. A run in which
        # the folds disagreed reports two estimands rather than one of them.
        estimand = cv$covariance_estimand %||% NA_character_,
        # The builder is a property of the base's own covariance, not of the
        # request, so it is derived from the fit rather than declared.
        # A builder assembles a plugged-in V; on a route where the engine
        # integrates over the covariance there is no V to build and the
        # field says so rather than naming a builder that never ran.
        builder = if (identical(cv$covariance_estimand, "integrated")) {
          "none: covariance integrated by the engine"
        } else if (inherits(bs, "try-error") || !isTRUE(bs$applicable)) {
          "random_term_reduction"
        } else {
          "block_cholesky"
        },
        integrated = identical(cv$covariance_estimand, "integrated")
      ),
      design = list(
        protected = attr(cv, "protected") %||% character(0),
        candidates = if (is.null(candidates)) {
          character(0)
        } else {
          attr(stats::terms(candidates, data = base$data), "term.labels")
        },
        column_map = column_map,
        # The declared candidate columns, the rank rule and its tolerance,
        # every column a fold dropped with its reason, and per column the
        # folds in which it was estimable. NULL where no fold was compiled.
        candidate_space = .fb_sel_candidate_space_record(candidate_space)
      ),
      resampling = list(
        unit = resamples$unit,
        v = resamples$v,
        repeats = resamples$repeats,
        n_folds = length(cv$folds %||% list()),
        prediction = cv$prediction,
        # The projection scores at the resampling unit inside the search
        # itself, so where one ran its own record is the authority; the
        # out-of-fold table's unit is the answer on every other route. The
        # projection comes first because a `benchmark = TRUE` projection run
        # now carries both, and the two name the unit differently -- the
        # search by the grouping variable, the table by the kind of unit.
        scoring_unit = projection$scoring$unit %||%
          cv$aggregated$se_unit[[1L]] %||% NA_character_,
        seed = cv$seed,
        metric = cv$metric
      ),
      penalty = if (length(block_record) == 0L) {
        "none: every candidate protected"
      } else {
        paste(
          vapply(block_record, function(b) b$family, character(1L)),
          collapse = " + "
        )
      },
      prior_blocks = block_record,
      # Whether the penalty strength was fixed by the caller, set by a rule,
      # tuned on the folds or learned inside the posterior, and how the
      # selection reused the data it is scored on.
      hyperparameter_treatment = .fb_select_hyperparameter_treatment(
        blocks, axes$workflow
      ),
      selection_data_reuse = .fb_select_data_reuse(axes, projection),
      # Which engine fitted the REFERENCE model, which is not always the
      # engine that fitted the base: a joint shrinkage prior is always
      # emitted on brms, whatever the base was fitted by, because the parsed
      # model is engine-independent and INLA's latent field is Gaussian
      # given the hyperparameters while a scale mixture is not.
      reference_engine = reference_engine,
      reference = if (is.null(adequacy)) {
        NULL
      } else {
        list(
          verdict = adequacy$verdict,
          reasons = adequacy$reasons,
          ppc = adequacy$ppc,
          convergence = adequacy$convergence,
          convergence_reasons = adequacy$convergence_reasons
        )
      },
      projection = if (is.null(projection)) NULL else projection$settings,
      benchmark = if (is.null(benchmark)) {
        NULL
      } else {
        lapply(benchmark, function(arm) {
          list(
            workflow = arm$workflow,
            provider = arm$provider,
            equivalence_level = arm$equivalence_level$level,
            equivalence_reason = arm$equivalence_level$reason,
            covariance_estimand = arm$covariance_estimand,
            provider_version = arm$provider_version,
            provider_lambda_convention = arm$provider_lambda_convention,
            # The optimiser's tuning report: the coordinate the point was
            # chosen on, the chosen point on it, and the lambda that point is
            # on each fold. NULL on an arm with no tuning.
            tuning_scale = arm$selection$tuning_scale,
            selected_s = arm$selection$selected_s,
            selected_lambda_canonical = arm$selection$selected_lambda_canonical,
            # Per penalised column: the folds, the folds that estimated it,
            # and the folds whose support held it at the chosen point.
            candidate_availability = arm$selection$candidate_availability
          )
        })
      },
      # `mode_oracle` is the Tier-1 claim and is empty unless the covariance
      # and the residual scale were both held fixed, which is never true of a
      # sampled route. The two predictive tables are the out-of-fold
      # comparison, one per model stage, and each is written only when every
      # arm in it carries the primary arm's score contract. Reading a glmnet
      # support as a check on a posterior is the category error the split
      # exists to prevent.
      mode_oracle = .fb_select_mode_oracle(benchmark %||% list()),
      reference_predictive = reference_predictive,
      decision_predictive = decision_predictive,
      provider_versions = provider_versions,
      # The information a prior scale was calibrated against. It is on the
      # certificate rather than implied because in clustered data the row
      # count overstates it, and a reader has to be able to see by how much:
      # the basis says whether the number counts resampling units or rows,
      # and `n_rows` sits beside it so the two can be compared without
      # going back to the data.
      calibration_n = calibration$n,
      calibration_basis = calibration$basis,
      calibration_unit = calibration$unit,
      n_rows = calibration$n_rows,
      # The constants each fold's out-of-fold score was computed under,
      # learned on that fold's analysis rows: centres and scales per column,
      # and per fold the calibration, `n_train` and the compiled prior. NULL
      # where no fold-local transform ran.
      fold_preprocessing = fold_preprocessing,
      fallback_fired = isTRUE(fallback_fired),
      # Named because a reader must be able to see what was not exact.
      approximations = as.character(approximations),
      refusals = as.character(refusals),
      contract = "SPEC_bayesian_selection_v5_2026-09-12"
    ),
    class = "fb_select_certificate"
  )
}


# --- the primary arm's own out-of-fold score ------------------------------ #

# .fb_select_reference_cv() --- the unit-level out-of-fold table for the
# REFERENCE model, on the folds the rest of the selection uses.
#
# Why it has to exist. The `penalty = NULL` route has always written a
# unit-level score table -- elpd, rmse and nlpd per resampling unit, under a
# named prediction target, on named folds -- and the penalised route wrote
# none. A `glmnet` arm was then placed beside a Bayesian answer under the
# label `posterior_predictive_comparable`, whose own definition is "scored on
# the same folds, target, loss, unit and candidate design". Nothing had been
# scored on the Bayesian side at all, so the label described a comparison
# that had not happened. The two are now scored the same way or the label is
# not printed.
#
# Why it is `benchmark = TRUE` only. It costs one reference-model refit per
# fold per model -- the protected baseline and the reference itself -- on top
# of the projection's own refits. That is the price of a comparison, and it
# is paid only by a caller who asked for one; the certificate says so where
# it was not paid.
#
# The machinery is `.fb_sel_bayes_cv()`'s: the same folds, the same
# relabelling of held-out grouping levels, the same predictive draws and the
# same aggregation. Two things differ. The candidate block arrives under a
# joint prior, and that prior -- with the standardising constants, the
# coefficient map, `calibration_n` and `n_train` it is stated against -- is
# re-resolved on each fold's analysis rows rather than taken from the
# whole-sample reference fit. A scale learned on every row would let an
# assessment row's covariate values set the training fold's shrinkage.
#
# Where the decision is a projection, `decision_search` carries the size the
# whole-sample decision chose and the search's settings. Each fold's
# reference is searched again to that size, projected onto the terms its own
# search chose, and scored on the same assessment rows: the decision stage's
# out-of-fold score. Projecting every fold onto the whole-sample choice would
# score a choice each fold's assessment rows helped make.
.fb_select_reference_cv <- function(base,
                                    candidates,
                                    penalty,
                                    protected_labels,
                                    resamples,
                                    prediction,
                                    metric,
                                    engine,
                                    seed,
                                    control = NULL,
                                    decision_search = NULL) {
  refit <- function(fold_data, fixed, fit_seed, fold_state) {
    extra <- setdiff(
      attr(stats::terms(fixed, data = fold_data), "term.labels"),
      protected_labels
    )
    if (length(extra) == 0L) {
      # The protected baseline. There is no candidate coefficient for the
      # compiled prior to name, so the fold is fitted without it -- which is
      # what makes the pair an increment over a baseline rather than two
      # runs of one model. It is fitted on brms like the reference, so the
      # two halves of the increment come from one engine.
      return(.fb_sel_bayes_refit(
        base, fold_data, fixed, "brms",
        engine = engine, seed = fit_seed
      ))
    }
    # The fold's own prior and column map, resolved on its analysis rows.
    transform <- fold_state$transform
    args <- list(
      base = base,
      candidates = stats::as.formula(
        paste0("~ ", paste(extra, collapse = " + "))
      ),
      prior = transform$compiled$prior,
      stanvars = transform$compiled$stanvars,
      column_map = transform$column_map,
      engine = engine,
      seed = fit_seed,
      data = fold_data
    )
    if (!is.null(control)) {
      args$control <- control
    }
    do.call(.fb_sel_reference_fit, args)
  }

  derived <- NULL
  if (!is.null(decision_search) && decision_search$size > 0L) {
    derived <- list(
      projected_submodel = function(fits, inputs, fit_seed) {
        terms <- .fb_sel_fold_search(fits$full, decision_search, fit_seed)
        predictive <- .fb_sel_projected_predictive(
          fits$full,
          newdata = inputs$assess_data,
          predictor_terms = terms,
          response = base$fb$response,
          seed = fit_seed
        )
        predictive$terms <- terms
        predictive
      }
    )
  }

  .fb_sel_bayes_cv(
    base = base,
    candidates = candidates,
    resamples = resamples,
    prediction = prediction,
    metric = metric,
    engine = engine,
    seed = seed,
    # The data as supplied. Each fold standardises its candidate columns with
    # its own constants, so nothing is transformed before the folds are cut.
    data = base$data,
    refit_fn = refit,
    # The provider is brms whatever fitted the base: a joint shrinkage prior
    # is always emitted there.
    provider = "brms",
    # The base's own recorded wall time, so the cost estimate costs no fit.
    probe_seconds = .fb_select_base_seconds(base),
    fold_transform = .fb_select_fold_transform(
      candidates, penalty, protected_labels, resamples
    ),
    derived_models = derived
  )
}


# --- the two stanzas that replace one overloaded label --------------------- #

# .fb_select_mode_oracle() --- the Tier-1 stanza, and why it is usually empty.
#
# A mode identity is defined on ONE substrate: the whitened design with the
# covariance plugged in and the residual scale held at unity. An arm reaches
# this stanza only when it solved the canonical objective on that substrate
# AND the comparison against the primary answer was made under the same two
# conditions. A sampled Bayesian route integrates the covariance and samples
# `sigma`, so the second never holds there and the stanza is empty on every
# sampled run.
#
# It is a stanza rather than an absent field because "no mode oracle ran" is
# the thing a reader has to be told: the arm's own internal fidelity to its
# objective is recorded beside it, and reading that as a statement about the
# posterior is the category error this split exists to prevent.
.fb_select_mode_oracle <- function(arms) {
  empty <- list(
    available = FALSE,
    level = "map_exact_fixed_covariance",
    arms = character(0),
    reason = paste0(
      "no arm carries a mode oracle: the identity is defined only with the ",
      "covariance plugged in and the residual scale held fixed, and the ",
      "Bayesian route integrates the first and samples the second. The ",
      "exact-limit bridges test that identity on the fixed-covariance ",
      "substrate; a posterior is not on it."
    )
  )
  if (length(arms) == 0L) {
    return(empty)
  }
  qualifies <- vapply(
    arms,
    function(arm) {
      own <- arm$provider_equivalence_level
      comparison <- arm$equivalence_level
      isTRUE(identical(own$level, "map_exact_fixed_covariance")) &&
        isTRUE(comparison$covariance_fixed) &&
        isTRUE(comparison$sigma_fixed)
    },
    logical(1L)
  )
  if (!any(qualifies)) {
    return(empty)
  }
  list(
    available = TRUE,
    level = "map_exact_fixed_covariance",
    arms = names(arms)[qualifies],
    reason = paste0(
      "the covariance and the residual scale were held fixed on both sides, ",
      "so the posterior mode is the minimiser of the named penalised ",
      "objective and a numerical oracle gates the claim."
    )
  )
}


# --- the score contract and the two stage tables -------------------------- #

# An out-of-fold comparison is a claim that the arms scored the same thing:
# the same units in the same folds, the same rows, the same prediction
# target, response, metric and unit, the same candidate design, the same
# preprocessing boundary and the same model stage. Each arm carries that
# claim as a score contract read off the fold object it was scored on, and
# an arm is placed beside the primary only when the two contracts agree
# field for field. The fold LABELS are not the partition: two assignments of
# units to folds called 1..K share a label set and nothing else.
#
# Two tables, because two different questions were sharing one. The
# reference stage places the Bayesian reference beside full-design
# comparators. The decision stage places the model each procedure chose
# beside the others' choices. A stage with no score says so; nothing is
# substituted into it.

# The canonical encoding the score contract hashes. Every value is written
# with an explicit type tag, explicit lengths and little-endian byte order,
# so one input gives one byte stream under any R version and on any
# platform. `serialize()` is not used: its header carries the version of R
# that wrote it. The digest is MD5, through tools::md5sum() over a temporary
# file (DEC-HASH-01), which needs nothing beyond base R and keeps the R floor.
# MD5 identifies the bytes here; it is not used for security, and the
# contract makes no claim against a constructed collision.
.FB_SCORE_CONTRACT_VERSION <- 2L


#' Write one value of the score contract's canonical encoding
#'
#' One ASCII tag byte, then the value, with every integer a 32-bit and every
#' double a 64-bit little-endian number. `S` character: the count, the byte
#' length of each element after `enc2utf8()` (`-1` for `NA`), then the bytes
#' of the elements that are not `NA`. `I` integer and `L` logical: the count,
#' then the values (`NA` is the smallest 32-bit integer). `D` double: the
#' count, an `I` record of the positions holding `NA`, an `I` record of the
#' positions holding `NaN`, then the values with `NA`, `NaN` and `-0` written
#' as `+0`, because the bit pattern of a `NaN` depends on the platform. `F`
#' factor: an `S` record of the levels, then an `I` record of the codes. A
#' matrix is written as its values in column order.
#'
#' @param con An open binary connection.
#' @param x An atomic vector or factor.
#' @return `NULL`, invisibly.
#' @noRd
#' @keywords internal
.fb_select_write_value <- function(con, x) {
  int32 <- function(v) {
    writeBin(as.integer(v), con, size = 4L, endian = "little")
  }
  tag <- function(t) writeBin(charToRaw(t), con)
  if (length(x) > .Machine$integer.max) {
    stop(
      "the score contract's encoding counts in 32-bit integers, and this ",
      "value has ", length(x), " elements.",
      call. = FALSE
    )
  }
  if (is.factor(x)) {
    tag("F")
    .fb_select_write_value(con, levels(x))
    .fb_select_write_value(con, as.integer(unclass(x)))
  } else if (is.character(x)) {
    x <- enc2utf8(x)
    absent <- is.na(x)
    bytes <- nchar(x, type = "bytes")
    bytes[absent] <- -1L
    tag("S")
    int32(length(x))
    int32(bytes)
    writeBin(charToRaw(paste(x[!absent], collapse = "")), con)
  } else if (is.logical(x) || is.integer(x)) {
    tag(if (is.logical(x)) "L" else "I")
    int32(length(x))
    int32(x)
  } else if (is.double(x)) {
    values <- as.double(x)
    nan <- is.nan(values)
    tag("D")
    int32(length(values))
    .fb_select_write_value(con, which(is.na(values) & !nan))
    .fb_select_write_value(con, which(nan))
    values[is.na(values) | values == 0] <- 0
    writeBin(values, con, size = 8L, endian = "little")
  } else {
    stop(
      "the score contract has no canonical encoding for a value of type `",
      typeof(x), "`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Write one named field of the score contract's canonical encoding
#'
#' @param con An open binary connection.
#' @param name The field's name, written as a one-element `S` record.
#' @param x The field's value.
#' @return `NULL`, invisibly.
#' @noRd
#' @keywords internal
.fb_select_write_field <- function(con, name, x) {
  .fb_select_write_value(con, as.character(name))
  .fb_select_write_value(con, x)
}


#' The MD5 digest of one part of the score contract
#'
#' The part is written behind a header naming it and the encoding version to
#' a temporary file in `tempdir()`, hashed by `tools::md5sum()`, and the file
#' is removed on exit. MD5 identifies the bytes; it is not used for security.
#'
#' @param label The part's name.
#' @param write A function of one open binary connection that writes the
#'   part's fields.
#' @return 32 lower-case hexadecimal characters.
#' @noRd
#' @keywords internal
.fb_select_digest <- function(label, write) {
  path <- tempfile("fb_score_contract_", tmpdir = tempdir(), fileext = ".bin")
  on.exit(unlink(path), add = TRUE)
  con <- file(path, open = "wb")
  con_open <- TRUE
  on.exit(if (con_open) close(con), add = TRUE, after = FALSE)
  .fb_select_write_field(con, "label", label)
  .fb_select_write_field(con, "version", .FB_SCORE_CONTRACT_VERSION)
  write(con)
  close(con)
  con_open <- FALSE
  digest <- unname(tools::md5sum(path))
  if (length(digest) != 1L || is.na(digest)) {
    stop(
      "the score contract's digest of `", label, "` could not be read back ",
      "from its temporary file.",
      call. = FALSE
    )
  }
  digest
}


#' The score contract of one arm's out-of-fold table
#'
#' Read off the fold object the arm was scored on. The unit-to-fold map is
#' canonical -- folds are numbered in the order their units first appear
#' when the units are sorted -- so one partition under two labellings has
#' one contract, and two partitions under one labelling do not.
#'
#' The response and the design enter by their values, each digested through
#' the canonical encoding (`.fb_select_write_value()`): the response resolved
#' through the fit's own resolver, with the offsets the benchmark scorer
#' subtracts, and every design column over all rows before fold-local
#' preprocessing. A rescaled covariate or a transformed response under an
#' unchanged name gives a different contract. The hash is recomputed from the
#' stored fields by `.fb_select_contract_hash()`, which is how a stage table
#' verifies it.
#'
#' @param folds An `<fb_sel_folds>`.
#' @param data The data the folds index.
#' @param response The response as written: a column name or an expression
#'   over the columns.
#' @param prediction The prediction target.
#' @param metric The metric the comparison is scored at.
#' @param design The design matrix the arm was scored on (`cv$design`): the
#'   protected and candidate columns over all rows, named.
#' @param learned_on Where the arm's fold-learned quantities were learned.
#' @param stage `"reference"` or `"decision"`.
#' @param scoring_unit The unit the arm aggregated at; `NULL` reads the
#'   folds' resampling unit, or `"row"` where none was declared.
#' @param offsets The model's offset expressions.
#' @param candidate_space The arm's `.fb_sel_candidate_space()`, whose
#'   per-fold dropped columns are a field, or `NULL` where no fold dropped a
#'   column.
#' @return An `fb_select_score_contract`.
#' @noRd
#' @keywords internal
.fb_select_score_contract <- function(folds,
                                      data,
                                      response,
                                      prediction,
                                      metric,
                                      design,
                                      learned_on,
                                      stage,
                                      scoring_unit = NULL,
                                      offsets = character(0),
                                      candidate_space = NULL) {
  # The fold of every scored row, per repeat -----
  n <- attr(folds, "n_rows") %||% nrow(data)
  unit_key <- attr(folds, "unit_key") %||% as.character(seq_len(n))
  repeats <- sort(unique(vapply(
    folds, function(s) as.integer(s$repeat_id), integer(1L)
  )))
  row_fold <- matrix(NA_integer_, nrow = n, ncol = length(repeats))
  for (i in seq_along(folds)) {
    split <- folds[[i]]
    row_fold[split$assess_rows, match(as.integer(split$repeat_id), repeats)] <-
      i
  }

  # Canonical fold numbers, by the sorted units -----
  scored <- rowSums(!is.na(row_fold)) > 0L
  units <- sort(unique(unit_key[scored]), method = "radix")
  first_row <- match(units, unit_key)
  canonical_rows <- row_fold
  seen_by_repeat <- list()
  unit_to_fold <- do.call(rbind, lapply(seq_along(repeats), function(r) {
    raw <- row_fold[first_row, r]
    seen <- unique(raw[!is.na(raw)])
    seen_by_repeat[[r]] <<- seen
    canonical_rows[, r] <<- match(row_fold[, r], seen)
    data.frame(
      repeat_id = rep(repeats[[r]], length(units)),
      unit = units,
      fold = match(raw, seen),
      stringsAsFactors = FALSE
    )
  }))

  # The values the arm was scored on, digested -----
  if (!is.matrix(design) || is.null(colnames(design)) ||
        nrow(design) != nrow(data)) {
    stop(
      "the score contract needs the design matrix the arm was scored on, ",
      "with named columns and one row per row of the data.",
      call. = FALSE
    )
  }
  y <- .fb_response_values(list(response = response), data)
  if (is.null(y)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_variable_not_in_data",
      message = paste0(
        "fb_select(): the response `", response, "` is neither a column of ",
        "the data nor an expression that evaluates over its rows, so the ",
        "score contract has no response values to identify."
      )
    ))
  }
  # Both scorers read a numeric response through as.numeric().
  if (is.numeric(y) || is.logical(y)) {
    y <- as.double(y)
  }
  offsets <- sort(as.character(offsets), method = "radix")
  response_values <- .fb_select_digest("response", function(con) {
    .fb_select_write_field(con, "expression", as.character(response))
    .fb_select_write_field(con, "values", y)
    .fb_select_write_field(con, "offsets", offsets)
    if (length(offsets) > 0L) {
      .fb_select_write_field(
        con, "offset_values", as.double(.fb_offset_total(offsets, data))
      )
    }
  })
  # Columns in radix order of their names, so the column order of the matrix
  # is not part of the contract and the row order within a column is.
  columns <- order(colnames(design), method = "radix")
  design_values <- .fb_select_digest("design", function(con) {
    .fb_select_write_field(con, "dim", dim(design))
    .fb_select_write_field(con, "columns", colnames(design)[columns])
    for (j in columns) {
      .fb_select_write_field(
        con, colnames(design)[[j]], as.double(design[, j])
      )
    }
  })

  # The columns each fold dropped, by canonical fold -----
  dropped <- candidate_space$dropped
  space_rows <- data.frame(
    repeat_id = integer(0), fold = integer(0), column = character(0),
    stringsAsFactors = FALSE
  )
  if (!is.null(dropped) && nrow(dropped) > 0L) {
    labels <- vapply(folds, function(s) s$label, character(1L))
    index <- match(dropped$fold, labels)
    r <- match(
      vapply(index, function(i) as.integer(folds[[i]]$repeat_id), integer(1L)),
      repeats
    )
    space_rows <- data.frame(
      repeat_id = repeats[r],
      fold = vapply(
        seq_along(index),
        function(k) match(index[[k]], seen_by_repeat[[r[[k]]]]),
        integer(1L)
      ),
      column = as.character(dropped$column),
      stringsAsFactors = FALSE
    )
    space_rows <- space_rows[order(
      space_rows$repeat_id, space_rows$fold, space_rows$column,
      method = "radix"
    ), , drop = FALSE]
    rownames(space_rows) <- NULL
  }

  # The remaining fields, and the hash over all of them -----
  if (is.null(scoring_unit)) {
    scoring_unit <- attr(folds, "unit")
    if (length(scoring_unit) == 0L) {
      scoring_unit <- "row"
    }
  }
  fields <- list(
    version = .FB_SCORE_CONTRACT_VERSION,
    unit_to_fold = unit_to_fold,
    row_fold = canonical_rows,
    scoring_unit = as.character(scoring_unit),
    prediction = as.character(prediction),
    response = list(
      name = as.character(response),
      offsets = offsets,
      values = response_values
    ),
    metric = as.character(metric),
    design_columns = sort(
      unique(as.character(colnames(design))), method = "radix"
    ),
    design_values = design_values,
    learned_on = as.character(learned_on),
    stage = as.character(stage),
    candidate_space = space_rows
  )
  structure(
    c(fields, list(units = units, hash = .fb_select_contract_hash(fields))),
    class = "fb_select_score_contract"
  )
}


#' The hash of a score contract, recomputed from its fields
#'
#' The fold maps and the scalar fields are written through the canonical
#' encoding beside the two value digests, so the hash depends on every value
#' and is recomputed without the data.
#'
#' @param contract A score contract, or the list of its fields.
#' @return 32 lower-case hexadecimal characters.
#' @noRd
#' @keywords internal
.fb_select_contract_hash <- function(contract) {
  .fb_select_digest("contract", function(con) {
    map <- contract$unit_to_fold
    .fb_select_write_field(con, "repeat_id", as.integer(map$repeat_id))
    .fb_select_write_field(con, "unit", as.character(map$unit))
    .fb_select_write_field(con, "fold", as.integer(map$fold))
    .fb_select_write_field(
      con, "row_fold_dim", as.integer(dim(contract$row_fold))
    )
    .fb_select_write_field(con, "row_fold", as.integer(contract$row_fold))
    for (field in c("scoring_unit", "prediction")) {
      .fb_select_write_field(con, field, as.character(contract[[field]]))
    }
    .fb_select_write_field(
      con, "response_name", as.character(contract$response$name)
    )
    .fb_select_write_field(
      con, "response_offsets", as.character(contract$response$offsets)
    )
    .fb_select_write_field(
      con, "response_values", as.character(contract$response$values)
    )
    for (field in c(
      "metric", "design_columns", "design_values", "learned_on", "stage"
    )) {
      .fb_select_write_field(con, field, as.character(contract[[field]]))
    }
    # Each fold's dropped candidate columns, by canonical fold. Written by
    # column and by type, as the unit map above is: a contract that dropped
    # a column is not the contract that kept it.
    space <- contract$candidate_space
    .fb_select_write_field(
      con, "candidate_space_repeat_id", as.integer(space$repeat_id)
    )
    .fb_select_write_field(
      con, "candidate_space_fold", as.integer(space$fold)
    )
    .fb_select_write_field(
      con, "candidate_space_column", as.character(space$column)
    )
  })
}


#' Verify a score contract before it is compared, or refuse by name
#'
#' A contract is verified when it carries this release's version, an MD5
#' hash, and a hash equal to the one recomputed from its fields. A contract
#' written before score contracts identified values carries no version field
#' and a polynomial checksum over names, and neither can be checked against
#' the values.
#'
#' @param contract The arm's score contract.
#' @param arm The arm's name, for the message.
#' @return `TRUE`, invisibly; otherwise a
#'   `selection_score_contract_unverifiable` refusal.
#' @noRd
#' @keywords internal
.fb_select_verify_contract <- function(contract, arm) {
  hash <- contract[["hash"]]
  problem <- if (!inherits(contract, "fb_select_score_contract")) {
    "it is not a score contract"
  } else if (is.null(contract[["version"]])) {
    paste0(
      "it carries no version field, so it was written before score ",
      "contracts identified values"
    )
  } else if (!identical(contract[["version"]], .FB_SCORE_CONTRACT_VERSION)) {
    paste0(
      "it is version ", paste(format(contract[["version"]]), collapse = " "),
      ", and this release verifies version ", .FB_SCORE_CONTRACT_VERSION
    )
  } else if (!is.character(hash) || length(hash) != 1L ||
               !grepl("^[0-9a-f]{32}$", hash)) {
    "its hash is not an MD5 digest of the canonical encoding"
  } else if (!identical(.fb_select_contract_hash(contract), hash)) {
    "its hash does not match its fields"
  } else {
    NULL
  }
  if (!is.null(problem)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_score_contract_unverifiable",
      message = paste0(
        "fb_select(): the score contract of the `", arm, "` arm cannot be ",
        "verified: ", problem, ". Its out-of-fold score is not placed ",
        "beside another arm's."
      )
    ))
  }
  invisible(TRUE)
}


#' The fields on which two score contracts disagree
#'
#' @param a,b Two `fb_select_score_contract` objects.
#' @return The names of the fields that are not identical.
#' @noRd
#' @keywords internal
.fb_select_contract_differences <- function(a, b) {
  fields <- c(
    "unit_to_fold", "row_fold", "scoring_unit", "prediction", "response",
    "metric", "design_columns", "design_values", "learned_on", "stage",
    "candidate_space"
  )
  fields[!vapply(
    fields, function(f) identical(a[[f]], b[[f]]), logical(1L)
  )]
}


#' One stage's out-of-fold comparison, or the reason there is none
#'
#' Each entry is `list(arm, model, clusters, contract, label)` and the first
#' is the primary arm. An arm enters the table only when its contract agrees
#' with the primary's field for field and the units its table scored are the
#' units its contract names; every other arm is listed with the reason. A
#' contract whose hash cannot be verified refuses with
#' `selection_score_contract_unverifiable`.
#' Scores are re-aggregated from the per-unit tables through one function.
#'
#' @param stage `"reference"` or `"decision"`.
#' @param entries The arms, primary first.
#' @param metric The metric every arm is scored at.
#' @param unavailable A reason the stage has no comparison, or `NULL`.
#' @return `list(available, stage, metric, arms, excluded, contract, reason)`.
#' @noRd
#' @keywords internal
.fb_select_stage_table <- function(stage,
                                   entries,
                                   metric,
                                   unavailable = NULL) {
  absent <- function(reason, excluded = NULL) {
    list(
      available = FALSE, stage = stage, metric = metric, arms = NULL,
      excluded = excluded, contract = NULL, reason = reason
    )
  }
  if (!is.null(unavailable)) {
    return(absent(unavailable))
  }
  if (length(entries) == 0L) {
    return(absent(paste0(
      "no arm carries an out-of-fold score at the ", stage, " stage."
    )))
  }

  # Admit each arm against the primary's contract -----
  primary <- entries[[1L]]
  kept <- list()
  excluded <- list()
  primary_admitted <- FALSE
  for (i in seq_along(entries)) {
    entry <- entries[[i]]
    clusters <- if (is.data.frame(entry$clusters)) {
      as.data.frame(entry$clusters)
    } else {
      NULL
    }
    if (!is.null(clusters)) {
      clusters <- clusters[clusters$model == entry$model, , drop = FALSE]
    }
    # An arm that carries a score table is compared only once its contract's
    # hash is verified against its fields; the primary is verified first. An
    # arm with no table is excluded below for that reason.
    if (!is.null(clusters) && nrow(clusters) > 0L) {
      .fb_select_verify_contract(entry$contract, entry$arm)
    }
    differing <- .fb_select_contract_differences(
      primary$contract, entry$contract
    )
    reason <- NULL
    if (is.null(clusters) || nrow(clusters) == 0L) {
      reason <- paste0(
        "it carries no out-of-fold score table for `", entry$model, "`"
      )
    } else if (!metric %in% names(clusters)) {
      reason <- paste0("its per-unit table carries no `", metric, "` column")
    } else if (length(differing) > 0L) {
      reason <- paste0(
        "its score contract differs from the primary arm's in ",
        paste(differing, collapse = ", ")
      )
    } else if (!setequal(
      unique(as.character(clusters$cluster)), entry$contract$units
    )) {
      reason <- paste0(
        "the units its table scored are not the units its score contract ",
        "names"
      )
    } else if (length(entry$contract$units) < 2L) {
      reason <- paste0(
        "its score contract names fewer than two units, so no standard ",
        "error exists"
      )
    }
    if (!is.null(reason)) {
      excluded[[length(excluded) + 1L]] <- data.frame(
        arm = entry$arm, model = entry$model, reason = reason,
        stringsAsFactors = FALSE
      )
      next
    }
    primary_admitted <- primary_admitted || i == 1L
    summary <- .fb_sel_aggregate_bayes_scores(
      clusters, metric = metric, unit = "cluster"
    )
    kept[[length(kept) + 1L]] <- data.frame(
      arm = entry$arm,
      model = entry$model,
      stage = stage,
      score = as.numeric(summary$score),
      se = as.numeric(summary$se),
      n_units = as.integer(summary$n_units),
      unit = paste(entry$contract$scoring_unit, collapse = "+"),
      contract_hash = entry$contract$hash,
      label = entry$label %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }
  excluded <- if (length(excluded) == 0L) NULL else do.call(rbind, excluded)
  named <- function(x) {
    paste(paste0(x$arm, " (", x$reason, ")"), collapse = "; ")
  }

  # The primary anchors the table, and a table needs a second arm -----
  if (!primary_admitted) {
    return(absent(paste0(
      "the primary arm cannot anchor the ", stage, " stage: ",
      named(excluded[1L, , drop = FALSE]), "."
    ), excluded))
  }
  if (length(kept) < 2L) {
    return(absent(paste0(
      "only the ", primary$arm, " arm carries a score under this contract ",
      "at the ", stage, " stage, so there is nothing to place beside it",
      if (is.null(excluded)) "" else paste0(". Excluded: ", named(excluded)),
      "."
    ), excluded))
  }
  list(
    available = TRUE,
    stage = stage,
    metric = metric,
    arms = do.call(rbind, kept),
    excluded = excluded,
    contract = primary$contract,
    reason = paste0(
      "every arm's score is the out-of-fold `", metric, "` aggregated at ",
      "the resampling unit under one score contract (",
      primary$contract$hash, "): one unit-to-fold map and set of scored ",
      "rows, one prediction target, response, scoring unit and candidate ",
      "design, one preprocessing boundary and one model stage. ",
      if (identical(stage, "reference")) {
        paste0(
          "Every arm is the full candidate design: the Bayesian reference ",
          "beside full-design comparators. "
        )
      } else {
        paste0(
          "Every arm is the model its own procedure chose, and each choice ",
          "was made on these folds, so every score carries that choice's ",
          "optimism and none is a nested estimate. "
        )
      },
      "Similar predictive risk is what this licenses; identical ",
      "coefficients or supports are not.",
      if (is.null(excluded)) "" else paste0(" Excluded: ", named(excluded), ".")
    )
  )
}


#' The reference and decision stage tables of one run
#'
#' @param base The normalised base fit.
#' @param axes The resolved axes.
#' @param cv The primary arm's out-of-fold object.
#' @param benchmark The benchmark object, or `NULL`.
#' @param metric The metric the tables are scored at.
#' @return `list(reference, decision)`, each a `.fb_select_stage_table()`.
#' @noRd
#' @keywords internal
.fb_select_stage_tables <- function(base, axes, cv, benchmark, metric) {
  arms <- benchmark$arms %||% list()
  contract_of <- function(source, stage, scoring_unit = NULL) {
    # A source that recorded no design scored nothing and carries no
    # contract; an arm that reaches a comparison without one refuses there.
    if (is.null(source$design)) {
      return(NULL)
    }
    .fb_select_score_contract(
      folds = source$folds,
      data = base$data,
      response = base$fb$response,
      prediction = source$prediction,
      metric = metric,
      design = source$design,
      learned_on = source$preprocessing$learned_on %||% NA_character_,
      stage = stage,
      scoring_unit = scoring_unit,
      offsets = base$fb$offsets %||% character(0),
      candidate_space = source$candidate_space
    )
  }
  bayes <- identical(axes$workflow, "bayes")
  bayes_scored <- bayes && is.data.frame(cv$clusters) &&
    nrow(cv$clusters) > 0L
  no_primary <- if (bayes && !bayes_scored) {
    paste0(
      "the primary arm carries no out-of-fold score table, so there is ",
      "nothing to place beside a benchmark. `benchmark = TRUE` scores the ",
      "reference model on the benchmark's own folds and is what makes this ",
      "comparison available."
    )
  } else if (length(arms) == 0L) {
    "no benchmark arm ran, so the primary arm's table stands alone."
  } else {
    NULL
  }

  # The reference stage: the full candidate design on every arm -----
  reference <- list()
  if (bayes_scored) {
    reference[[length(reference) + 1L]] <- list(
      arm = "bayes", model = "full", clusters = cv$clusters,
      contract = contract_of(cv, "reference"),
      label = "the reference posterior over the full candidate design"
    )
  }
  if (!is.null(arms$gls)) {
    reference[[length(reference) + 1L]] <- list(
      arm = "gls", model = "full", clusters = benchmark$cv$clusters,
      contract = contract_of(benchmark$cv, "reference"),
      label = paste0(
        "generalised least squares on a plugged-in fold covariance over ",
        "the full candidate design"
      )
    )
  }
  if (!is.null(arms$glmnet$clusters_by_model$fixed_1)) {
    reference[[length(reference) + 1L]] <- list(
      arm = "glmnet", model = "fixed_1",
      clusters = arms$glmnet$clusters_by_model$fixed_1,
      contract = contract_of(
        benchmark$cv, "reference", arms$glmnet$scoring_unit
      ),
      label = paste0(
        "the penalised optimiser at the lambda the prior's hyperparameters ",
        "map to, with no selection step"
      )
    )
  }

  # The decision stage: the model each procedure chose -----
  decision <- list()
  decision_unavailable <- NULL
  if (bayes && identical(axes$prior, "none")) {
    decision_unavailable <- paste0(
      "every candidate is protected on this route, so no selection was ",
      "made and there is no decision stage to score."
    )
  } else if (bayes && identical(axes$decision, "posterior_threshold")) {
    decision_unavailable <- paste0(
      "a posterior threshold retains coefficients read off one posterior ",
      "and defines no predictive distribution of its own, so it has no ",
      "out-of-fold score; nothing is put in its place."
    )
  } else if (bayes_scored && identical(axes$decision, "projection")) {
    if ("projected_submodel" %in% cv$clusters$model) {
      decision[[1L]] <- list(
        arm = "bayes", model = "projected_submodel", clusters = cv$clusters,
        contract = contract_of(cv, "decision"),
        label = paste0(
          "the submodel at the chosen size, its terms searched again on each ",
          "fold's training rows, projected from that fold's reference and ",
          "scored on that fold's assessment rows"
        )
      )
    } else {
      decision_unavailable <- paste0(
        "the projected submodel was not scored out of fold, so the decision ",
        "stage has no Bayesian score; the reference is not put in its place."
      )
    }
  } else if (bayes_scored && identical(axes$decision, "full")) {
    decision[[1L]] <- list(
      arm = "bayes", model = "full", clusters = cv$clusters,
      contract = contract_of(cv, "decision"),
      label = paste0(
        "decision = \"full\" keeps every candidate, so the chosen model is ",
        "the reference posterior itself"
      )
    )
  }
  if (is.null(decision_unavailable) && isTRUE(arms$glmnet$selection$ran)) {
    chosen <- arms$glmnet$selection$selected_model
    decision[[length(decision) + 1L]] <- list(
      arm = "glmnet", model = chosen,
      clusters = arms$glmnet$clusters_by_model[[chosen]],
      contract = contract_of(
        benchmark$cv, "decision", arms$glmnet$scoring_unit
      ),
      label = paste0(
        "the optimiser's tuning point on ",
        arms$glmnet$selection$tuning_scale %||% "lambda",
        ", chosen on these folds"
      )
    )
  }

  list(
    reference = .fb_select_stage_table(
      "reference", reference, metric, unavailable = no_primary
    ),
    decision = .fb_select_stage_table(
      "decision", decision, metric,
      unavailable = decision_unavailable %||% no_primary
    )
  )
}


# --- the benchmark arms --------------------------------------------------- #

# .fb_select_benchmark() --- the two comparison arms, on the same folds.
#
# Neither is labelled Bayesian and neither can be. The GLS arm fits REML
# variance components independently in each fold and plugs them in; the
# `glmnet` arm solves a penalised objective at a mode. Both are scored on the
# same folds, the same prediction target and the same candidate design as the
# primary answer, which is what makes the comparison a comparison.
#
# The equivalence level stored on each arm is the level of the comparison
# AGAINST THE PRIMARY ANSWER, not the arm's internal fidelity to its own
# objective: the covariance is integrated and the residual scale is free on
# the Bayesian side, so a mode identity cannot hold whatever the provider can
# express. The provider's own internal level travels beside it.
.fb_select_benchmark <- function(base,
                                 candidates,
                                 protected,
                                 resamples,
                                 prediction,
                                 seed,
                                 metric,
                                 prior_family,
                                 blocks = NULL,
                                 decision = NULL,
                                 shared_predictive_table = FALSE) {
  arms <- list()
  refusals <- character(0)

  gls <- .fb_sel_cv(
    data = base$data,
    response = base$fb$response,
    protected = protected,
    candidates = candidates,
    random = .fb_sel_benchmark_random(base$fb),
    resamples = resamples,
    seed = seed,
    prediction = prediction,
    metric = metric,
    # The base's offset, subtracted from the response in each fold before
    # whitening; the penalised optimiser reads the same substrates.
    offsets = base$fb$offsets %||% character(0)
  )
  arms$gls <- list(
    workflow = gls$workflow,
    provider = gls$provider,
    fully_bayesian = FALSE,
    scores = gls$aggregated,
    covariance_estimand = gls$covariance_estimand,
    provider_version = NA_character_,
    provider_lambda_convention = NULL,
    equivalence_level = .fb_sel_equivalence_level(
      prior_family = prior_family,
      provider = "gls_reml_plugin",
      covariance_fixed = FALSE,
      sigma_fixed = FALSE,
      shared_predictive_table = shared_predictive_table
    ),
    metric = gls$metric
  )

  if (!is.null(blocks) && length(blocks) > 0L) {
    families <- unique(vapply(blocks, function(b) b$family, character(1L)))
    grouped_provider <- length(families) > 0L &&
      all(families %in% .fb_sel_grouped_families())
    if (!grouped_provider && !.fb_sel_provider_installed("glmnet")) {
      refusals <- c(refusals, paste0(
        "selection_provider_absent: the `glmnet` benchmark arm was ",
        "requested and `glmnet` is not installed. The Bayesian answer ",
        "stands; no substitute optimiser was chosen."
      ))
    } else {
      # The scoring unit is the RESAMPLING unit, as it is on every other
      # arm. Left at its own default the optimiser path scores at the fold,
      # so a side-by-side table would have put a twelve-site standard error
      # beside a two-fold one and called them one loss.
      cluster_key <- if (length(resamples$unit %||% character(0)) > 0L) {
        .fb_sel_key(base$data, resamples$unit)
      } else {
        NULL
      }
      path <- tryCatch(
        .fb_sel_glmnet_path(
          substrate = gls$substrates,
          blocks = blocks,
          folds = gls$folds,
          decision = list(
            rule = decision$rule %||% "one_se",
            metric = if (identical(metric, "mse")) "mse" else metric,
            tolerance = decision$tolerance,
            benchmark = TRUE
          ),
          seed = seed,
          cluster = cluster_key
        ),
        # A refusal IS an error condition (class `flexybayes_refusal` sits
        # above `error`), so one handler catches both and the reason code is
        # read off it where there is one. The benchmark arm is a comparison,
        # not the answer: its refusal is carried on the certificate rather
        # than raised over the top of a Bayesian result that stands.
        error = function(e) e
      )
      if (inherits(path, "condition")) {
        refusals <- c(refusals, paste0(
          (path$reason_code %||% "error"), ": ", conditionMessage(path)
        ))
      } else {
        arms$glmnet <- list(
          workflow = path$workflow,
          provider = path$provider,
          fully_bayesian = path$fully_bayesian,
          family = path$family,
          selection = path$selection,
          coefficients = path$coefficients_original_scale,
          support = path$support_by_resolution,
          kkt = path$kkt_by_lambda,
          covariance_estimand = path$covariance_estimand,
          provider_version = path$provider_version,
          provider_lambda_convention = path$provider_lambda_convention,
          provider_equivalence_level = path$equivalence_level,
          # The per-unit tables of the path's models, and the unit they were
          # aggregated at, so a stage table re-aggregates the model it names.
          clusters_by_model = path$clusters_by_model,
          assessment_predictions = path$fitted_values,
          scoring_unit = if (is.null(cluster_key)) "fold" else resamples$unit,
          equivalence_level = .fb_sel_equivalence_level(
            prior_family = prior_family,
            provider = path$provider,
            covariance_fixed = FALSE,
            sigma_fixed = FALSE,
            shared_predictive_table = shared_predictive_table
          )
        )
      }
    }
  }

  list(arms = arms, refusals = refusals, folds = gls$folds, cv = gls)
}


# --- the verb ---------------------------------------------------------- #

#' Bayesian selection over a fitted flexyBayes model
#'
#' Fits a Bayesian reference model in which a declared block of candidate
#' fixed effects carries a joint shrinkage prior and the covariance structure
#' is integrated over rather than plugged in; selects a parsimonious submodel
#' by projection predictive inference where the structure permits it, and
#' refuses by name where it does not; and, on request, triangulates the
#' result against an external penalised optimiser on the same folds.
#'
#' @details
#' Every result records four independent fields rather than one `engine`
#' argument, because conflating the inferential workflow with the software
#' that ran it is how a penalised-likelihood mode comes to be reported under
#' a Bayesian label:
#'
#' \describe{
#'   \item{`workflow`}{`"bayes"`, `"penalized_optimizer"` or
#'     `"unpenalized_benchmark"`.}
#'   \item{`provider`}{`"brms"`, `"inla"`, `"glmnet"`, `"prox"` or
#'     `"gls_reml_plugin"`; `"prox"` is the package's deterministic
#'     proximal solver for grouped penalties.}
#'   \item{`prior`}{the joint prior, or `"none"`.}
#'   \item{`decision`}{`"full"`, `"projection"`, `"posterior_threshold"` or
#'     `"optimizer_support"`.}
#' }
#'
#' Posterior draws under a continuous shrinkage prior contain no exact zeros,
#' so "selected" is always a named decision and never a property of the
#' draws, and a posterior-threshold retention is not an inclusion
#' probability. Every comparison against an external optimiser carries one of
#' three equivalence levels -- `map_exact_fixed_covariance`,
#' `posterior_predictive_comparable` or `no_identity_expected` -- and none of
#' them licenses the sentence "the same as glmnet". The mode identity a
#' ridge or Laplace prior has with a penalised objective holds only with the
#' covariance plugged in and the residual scale held fixed; the Bayesian
#' route integrates the covariance and samples `sigma`, so no arm of it
#' reaches `map_exact_fixed_covariance`.
#'
#' Where the decision is `"projection"`, the search and the predictive
#' evaluation run on CLUSTERED DRAWS rather than on the whole posterior --
#' 20 clusters for the search, 100 for the evaluation. Projecting a
#' multilevel submodel is a fit per draw, so projecting every draw at every
#' search step inside every fold does not finish on a problem of the size
#' this verb is for. Both counts are recorded on the result and on its
#' certificate and printed on the object's first line, so a reader can see
#' what the elpd column was computed from.
#'
#' `penalty = NULL` keeps its earlier meaning: every candidate is protected
#' and the fit is the extended base model. Under `workflow = "bayes"` it
#' refits the base's own engine inside every assessment fold and scores
#' held-out clusters by posterior predictive draws, which is the genuinely
#' Bayesian no-penalty reference and the correct baseline for every shrinkage
#' family. The closed-form generalised-least-squares route on a plugged-in
#' fold covariance remains available through `benchmark = TRUE` and is
#' labelled `unpenalized_benchmark` wherever it appears.
#'
#' Fold refits run one at a time by default.
#' `options(flexyBayes.selection_workers = n)` runs `n` of them at once on a
#' platform that forks, and divides the engine threads between the workers so
#' that workers times threads does not exceed the machine. Nothing the run
#' records changes: the fold allocation, each fold's candidate space, the
#' score contract, the certificate and the selected support are identical to
#' the serial run's, and the worker count is deliberately absent from the
#' certificate, which describes the estimate while the worker count describes
#' the machine that produced it. Engine numbers move as they move between any
#' two repeated runs on an engine whose posterior sampler is not
#' bit-reproducible. Set `n` against MEMORY: every concurrent fold holds its
#' own rows and its own engine process at the same time, so the count a
#' machine supports follows from its memory, and usually sits well below its
#' core count. A worker the operating system kills under that pressure
#' refuses by name, and the aggregation keeps every fold or none.
#'
#' `fallback = "none"` is the default. A Bayesian request whose provider is
#' absent refuses with `selection_provider_absent`, naming the package and
#' the install line. `fallback = "glmnet"` is a deliberate choice; when it
#' fires, `workflow` becomes `penalized_optimizer`, the object prints the
#' substitution on its first line, and the certificate records
#' `fallback_fired = TRUE`. There is no configuration that produces a
#' penalised-likelihood answer under a Bayesian label.
#'
#' Each fold's candidate space is fixed from the covariates on its training
#' rows before any fold is fitted. The protected columns, then the candidate
#' columns in the order `candidates` writes them, are factorised by
#' `qr()` with tolerance `1e-7`, every column but the intercept centred on
#' the training rows. A protected column aliased there refuses with
#' `selection_fold_protected_unidentified`, and a factor level present only in
#' a fold's assessment rows refuses with `selection_fold_level_absent`, both
#' before any engine runs. A candidate column aliased there is dropped for that
#' fold, as when a fold's training trials leave fewer degrees of freedom than
#' its trial-level candidates need: the fold's fits and scores use the columns
#' it keeps, every coefficient it reports keeps the declared columns in order
#' with `NA` for a dropped one, and the certificate's
#' `design$candidate_space` names each drop with its reason and counts, for
#' each candidate column, the folds in which it was estimable. The `glmnet`
#' arm's `selection$candidate_availability` places that count beside the
#' number of folds whose support holds the column at the chosen point. The
#' Bayesian route drops whole terms and refuses a term aliased only in part, and
#' the penalised reference's out-of-fold score refuses any drop, both with
#' `selection_fold_candidate_drop_not_representable`.
#'
#' @param base A fitted flexyBayes model. Its fixed effects are the protected
#'   block and are never penalised.
#' @param candidates A one-sided formula over columns of the base fit's data.
#'   A transform whose value on a row depends on other rows, in the
#'   candidates or in the base's fixed part, is learned inside each fold or
#'   refused: the benchmark workflows learn `scale()`, `poly()`,
#'   `splines::ns()` and `splines::bs()` on each fold's training rows and
#'   refuse any other such transform, and the Bayesian route's fold refits
#'   refuse every one, with `selection_formula_transform_not_fold_local`.
#' @param penalty `NULL` (every candidate protected) or an [fb_penalty()]
#'   declaring the joint shrinkage prior over the candidate block.
#' @param workflow `"bayes"` (default), `"penalized_optimizer"` or
#'   `"unpenalized_benchmark"`. The third asks for the closed-form
#'   generalised-least-squares arm on its own -- the same computation
#'   `benchmark = TRUE` adds beside a Bayesian answer -- for a caller who
#'   wants the plugged-in-covariance comparison and not the posterior. It is
#'   never labelled Bayesian, on the object, in the print method or on the
#'   certificate.
#' @param fallback `"none"` (default) or `"glmnet"`. Fires only when the
#'   Bayesian provider is absent, never after a failure.
#' @param benchmark `TRUE` also runs the unpenalised generalised-least-squares
#'   arm and, where a penalty was given, the penalised-optimiser arm on the
#'   same folds (`glmnet` for ridge, Laplace and elastic net; `prox` for
#'   grouped penalties) -- and scores the reference model itself out of fold,
#'   so the arms have
#'   one loss at one unit to be placed beside. That last part is what the
#'   certificate's `reference_predictive` and `decision_predictive` tables
#'   report -- the full-design models at one stage, the chosen models at the
#'   other, each written only when every arm shares one score contract,
#'   which identifies the folds, the response and every design column by
#'   their values --
#'   and what the equivalence level `posterior_predictive_comparable`
#'   requires; without it the certificate says so by name rather than
#'   claiming a comparison nobody made. It costs one reference refit per
#'   fold per model, and each fold's preprocessing and prior are learned on
#'   that fold's analysis rows. A projection's submodel is scored with terms
#'   searched again on each fold's reference fit, at the size the
#'   whole-sample decision chose. The optimiser arm's tuning point is reported
#'   as `s = lambda / lambda_max`, each fold's `lambda_max` read on its
#'   training rows, beside the lambda that point is on each fold.
#' @param resamples An [fb_resamples()] specification.
#' @param prediction `"new_cluster"` (default) or `"known_cluster"`. A
#'   held-out cluster's own effect is not estimable from the training rows,
#'   so scoring it conditionally would reward the model for a random effect
#'   it was handed. Under `"new_cluster"` every column a random term's effect
#'   is indexed by is relabelled on the assessment rows; a `vm()`, `ped()` or
#'   autoregressive term cannot take a new level and refuses with
#'   `selection_new_cluster_not_representable`. Each assessment row's level
#'   `L` becomes `L::fb_heldout`, checked row by row before any fold is
#'   fitted; a fold in which a training row already carries such a label
#'   refuses with the same code.
#' @param decision An [fb_decision()]: the method, metric, rule and search.
#' @param engine `NULL`, or a named list setting `chains`, `n_samples`,
#'   `warmup` and `adapt_delta` for the brms reference fit. Deprecated
#'   [fb_engine()] objects remain accepted during the compatibility period.
#' @param seed Required. The sampler, the posterior predictive check and the
#'   projection search all consume the random-number stream, and a default
#'   seed is an unreproducible run that looks reproducible.
#' @param budget Seconds. The route refuses before the first fit when the
#'   estimated cost exceeds it.
#' @param plan `TRUE` returns the resolved axes, the prior blocks, the
#'   projection structure gate and a cost estimate without fitting anything.
#' @param ... Reserved. `metric` and `covariance` are accepted here for one
#'   release with a deprecation notice; `metric` now lives in
#'   [fb_decision()], and the covariance summary is a property of the
#'   benchmark workflow and is read off the folds rather than declared.
#'
#' @returns An object of class `fb_select` carrying the four resolved axes,
#'   the reference fit, the column map, the decision, the out-of-fold scores,
#'   the increment over the protected model, the diagnostics, any benchmark
#'   arms, the folds and the certificate. With `plan = TRUE`, an object of
#'   class `fb_select_plan`.
#'
#' @seealso [fb_penalty()] for the prior, [fb_resamples()] for the folds,
#'   [fb_decision()] for the rule.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE) &&
#'   requireNamespace("brms", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(
#'     site = factor(rep(1:8, each = 10)),
#'     x1 = rnorm(80),
#'     x2 = rnorm(80)
#'   )
#'   d$y <- 1 + 0.5 * d$x1 + rnorm(8)[d$site] + rnorm(80)
#'   base <- flexybayes(y ~ 1 + (1 | site), data = d, backend = "inla",
#'                      verbose = FALSE)
#'   # The plan resolves the axes, the prior blocks and the structure gate,
#'   # and estimates the cost, without fitting the reference model.
#'   fb_select(
#'     base,
#'     candidates = ~ x1 + x2,
#'     penalty = fb_penalty(term(everything()) ~ laplace(scale = 1)),
#'     seed = 1L,
#'     plan = TRUE
#'   )
#' }
#' }
#' @export
fb_select <- function(base,
                      candidates,
                      penalty = NULL,
                      workflow = c("bayes", "penalized_optimizer",
                                   "unpenalized_benchmark"),
                      fallback = c("none", "glmnet"),
                      benchmark = FALSE,
                      resamples = fb_resamples(),
                      prediction = c("new_cluster", "known_cluster"),
                      decision = fb_decision(),
                      engine = NULL,
                      seed = NULL,
                      budget = NULL,
                      plan = FALSE,
                      ...) {
  workflow <- match.arg(workflow)
  fallback <- match.arg(fallback)
  prediction <- match.arg(prediction)
  decision <- .fb_select_check_decision(decision)
  dots <- .fb_select_absorb_deprecated(list(...))
  # Read here so that a request for process isolation the session cannot
  # honour refuses before any fold runs (R/process_isolation.R).
  .fb_isolation_active("fb_select()")

  if (!is.null(engine) && !inherits(engine, "fb_engine")) {
    if (!is.list(engine) || is.null(names(engine)) || any(!nzchar(names(engine)))) {
      stop("`engine` must be NULL or a named list of brms sampler settings.", call. = FALSE)
    }
    unknown_engine_opts <- setdiff(names(engine), .fb_engine_opt_names)
    if (length(unknown_engine_opts)) {
      stop(
        "`engine` contains unrecognised setting",
        if (length(unknown_engine_opts) == 1L) " " else "s ",
        paste(shQuote(unknown_engine_opts), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  if (!is.logical(benchmark) || length(benchmark) != 1L || is.na(benchmark)) {
    stop("`benchmark` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(plan) || length(plan) != 1L || is.na(plan)) {
    stop("`plan` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(penalty) && !inherits(penalty, "fb_penalty")) {
    stop(.fb_refusal_condition(
      reason_code = "selection_penalty_spec_malformed",
      message = paste0(
        "fb_select(): `penalty` must be NULL or an <fb_penalty>. A plain ",
        "list of tuning arguments is the earlier release's shape; the ",
        "constructor validates at construction and names what it refuses."
      )
    ))
  }

  # The gate reads the fit as it stands, before any normalisation: a fit
  # that carries neither the parsed model nor the data must refuse under the
  # gate's own name rather than under a shim's. The normalisation that
  # follows is a compatibility step for fits made before the brms emit
  # started writing `$fb` and `$data`, and for the `$extras$fb_terms`
  # spelling the shared dispatch path uses.
  .fb_select_gate(base)
  base <- .fb_sel_bayes_normalise_base(base)
  engine_name <- .fb_select_engine_of(base)

  protected_labels <- .fb_sel_protected_labels(base$fb)
  protected <- stats::as.formula(paste0(
    "~ ",
    if (length(protected_labels) > 0L) {
      paste(protected_labels, collapse = " + ")
    } else {
      "1"
    }
  ))

  metric <- decision$metric

  # The projection structure gate, before anything runs. It is the first
  # thing the plan and the printed result report, and it is what resolves
  # `decision = "auto"`: projection where the reference model's structure is
  # one the provider can project, full where it is not. The default used to
  # be projection outright, so a genomic or a G-by-E model refused only
  # after a reference model had been sampled -- and most of what this
  # package fits is red at that gate.
  structure_gate <- tryCatch(
    .fb_sel_structure_gate(base$fb),
    error = function(e) NULL
  )
  resolved_decision <- .fb_select_resolve_decision(decision, structure_gate)
  decision <- resolved_decision$decision
  decision_resolution <- list(
    requested = resolved_decision$decision$requested_method,
    resolved_from = resolved_decision$resolved_from,
    gate_verdict = resolved_decision$gate_verdict,
    reason = resolved_decision$reason
  )

  # The decision and its metric are one request. A pair the projection cannot
  # choose on refuses here, on the resolved method, before the provider check,
  # the budget and the first fit -- and before the `plan = TRUE` branch, so
  # the plan returns the same verdict the run would.
  if (identical(workflow, "bayes") && !is.null(penalty)) {
    .fb_sel_require_projection_metric(
      decision$method,
      decision$metric,
      requested_method = decision$requested_method
    )
  }

  # The information the prior scales and the certificate are stated against.
  # Resolved once, before any route branches, so every route reports the
  # same number.
  calibration <- .fb_select_calibration_n(resamples, base$data)

  # --- provider availability, and the one place a substitution may happen -- #
  fallback_fired <- FALSE
  requested_workflow <- workflow
  if (identical(workflow, "bayes")) {
    needed <- .fb_select_required_providers(engine_name, penalty, decision)
    missing <- .fb_select_missing_providers(needed)
    if (length(missing) > 0L) {
      if (identical(fallback, "none")) {
        stop(.fb_refusal_condition(
          reason_code = if (is.null(penalty)) {
            "selection_provider_absent"
          } else {
            "selection_fallback_not_authorised"
          },
          message = paste0(
            "fb_select(): the Bayesian workflow needs ",
            paste(missing, collapse = ", "),
            ", which ", if (length(missing) > 1L) "are" else "is",
            " not installed, and `fallback = \"none\"` is the default. ",
            "Nothing was substituted: a penalised optimiser is a different ",
            "estimator and reporting it under a Bayesian label is the one ",
            "thing this argument exists to prevent. Install with ",
            paste(
              vapply(missing, .fb_select_install_line, character(1L)),
              collapse = " ; "
            ),
            ", or pass fallback = \"glmnet\" to run the penalised optimiser ",
            "under its own name."
          ),
          missing = missing
        ))
      }
      if (is.null(penalty)) {
        stop(.fb_refusal_condition(
          reason_code = "selection_provider_absent",
          message = paste0(
            "fb_select(): the Bayesian workflow needs ",
            paste(missing, collapse = ", "),
            ", which is not installed, and `fallback = \"glmnet\"` cannot ",
            "stand in for `penalty = NULL`: there is no penalised problem ",
            "for an optimiser to solve. Install ",
            paste(
              vapply(missing, .fb_select_install_line, character(1L)),
              collapse = " ; "
            ),
            ", or ask for the unpenalised benchmark with benchmark = TRUE."
          ),
          missing = missing
        ))
      }
      if (!.fb_sel_provider_installed("glmnet")) {
        stop(.fb_refusal_condition(
          reason_code = "selection_provider_absent",
          message = paste0(
            "fb_select(): the Bayesian provider (",
            paste(missing, collapse = ", "),
            ") is absent and `fallback = \"glmnet\"` was asked for, but ",
            "`glmnet` is not installed either. Nothing was substituted."
          ),
          missing = c(missing, "glmnet")
        ))
      }
      workflow <- "penalized_optimizer"
      fallback_fired <- TRUE
    }
  }

  if (identical(workflow, "penalized_optimizer") && is.null(penalty)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_penalty_spec_malformed",
      message = paste0(
        "fb_select(): workflow = \"penalized_optimizer\" needs a penalty. ",
        "With every candidate protected there is no penalised objective for ",
        "an optimiser to minimise; the unpenalised comparison is ",
        "benchmark = TRUE."
      )
    ))
  }

  # The two routes below have no decision layer at all -- there is no
  # shrinkage prior for a submodel to be projected out of -- so they report
  # `full` however the gate read, and the resolution record has to say that
  # rather than the projection it would have picked for a penalised request.
  forced_full <- list(
    requested = decision_resolution$requested,
    resolved_from = if (identical(decision_resolution$requested, "auto")) {
      "auto"
    } else {
      NA_character_
    },
    gate_verdict = decision_resolution$gate_verdict,
    reason = paste0(
      "every candidate is protected on this route, so there is no submodel ",
      "to decide between and the answer is the extended model itself"
    )
  )

  if (identical(workflow, "unpenalized_benchmark")) {
    if (!is.null(penalty)) {
      stop(.fb_refusal_condition(
        reason_code = "selection_penalty_spec_malformed",
        message = paste0(
          "fb_select(): workflow = \"unpenalized_benchmark\" is the ",
          "closed-form generalised-least-squares arm, which has no penalty ",
          "term. Ask for the penalised optimiser by name, or for the ",
          "Bayesian workflow with the penalty and benchmark = TRUE."
        )
      ))
    }
    bench <- .fb_select_benchmark(
      base = base,
      candidates = candidates,
      protected = protected,
      resamples = resamples,
      prediction = prediction,
      seed = seed,
      metric = .fb_select_benchmark_metric(metric),
      prior_family = "none",
      blocks = NULL,
      decision = decision,
      shared_predictive_table = TRUE
    )
    return(.fb_select_assemble(
      base = base,
      axes = list(
        workflow = bench$cv$workflow,
        provider = bench$cv$provider,
        prior = "none",
        decision = "full"
      ),
      cv = bench$cv,
      resamples = resamples,
      candidates = candidates,
      metric = metric,
      calibration = calibration,
      provider_versions = list(),
      structure_gate = structure_gate,
      decision_resolution = forced_full,
      dots = dots
    ))
  }

  if (identical(workflow, "penalized_optimizer")) {
    return(.fb_select_optimizer(
      base = base,
      candidates = candidates,
      penalty = penalty,
      protected = protected,
      protected_labels = protected_labels,
      resamples = resamples,
      prediction = prediction,
      decision = decision,
      seed = seed,
      metric = metric,
      plan = plan,
      fallback_fired = fallback_fired,
      requested_workflow = requested_workflow,
      calibration = calibration,
      structure_gate = structure_gate,
      decision_resolution = list(
        requested = decision_resolution$requested,
        resolved_from = if (
          identical(decision_resolution$requested, "auto")
        ) {
          "auto"
        } else {
          NA_character_
        },
        gate_verdict = decision_resolution$gate_verdict,
        reason = paste0(
          "the penalised optimiser reports the support of its own solution ",
          "at the chosen lambda; there is no posterior for a projection to ",
          "act on"
        )
      )
    ))
  }

  # --- the Bayesian route, no penalty ------------------------------------- #
  if (is.null(penalty)) {
    if (isTRUE(plan)) {
      return(.fb_plan_selection(
        base = base,
        axes = list(
          workflow = "bayes",
          provider = engine_name,
          prior = "none",
          decision = "full"
        ),
        blocks = NULL,
        structure_gate = structure_gate,
        resamples = resamples,
        prediction = prediction,
        fits_per_fold = 2L,
        candidates = candidates,
        decision_resolution = forced_full
      ))
    }
    cv <- .fb_sel_bayes_cv(
      base = base,
      candidates = candidates,
      resamples = resamples,
      prediction = prediction,
      metric = metric,
      engine = engine,
      seed = seed,
      budget = budget
    )
    bench <- NULL
    if (isTRUE(benchmark)) {
      bench <- .fb_select_benchmark(
        base = base,
        candidates = candidates,
        protected = protected,
        resamples = resamples,
        prediction = prediction,
        seed = seed,
        metric = .fb_select_benchmark_metric(metric),
        prior_family = "none",
        blocks = NULL,
        decision = decision,
        # This route always writes a unit-level out-of-fold table, so the
        # comparison the middle equivalence level names is one that
        # happened.
        shared_predictive_table = TRUE
      )
    }
    return(.fb_select_assemble(
      base = base,
      axes = list(
        workflow = cv$workflow,
        provider = cv$provider,
        prior = "none",
        decision = "full"
      ),
      cv = cv,
      resamples = resamples,
      candidates = candidates,
      benchmark = bench,
      metric = metric,
      approximations = character(0),
      calibration = calibration,
      provider_versions = .fb_select_provider_versions(cv$provider),
      structure_gate = structure_gate,
      decision_resolution = forced_full,
      dots = dots
    ))
  }

  # --- the Bayesian route, with a joint shrinkage prior -------------------- #
  resolved <- .fb_select_resolve_penalty(
    base,
    candidates,
    penalty,
    protected_labels,
    calibration = calibration,
    workflow = "bayes"
  )
  prior_family <- paste(
    unique(vapply(resolved$blocks, function(b) b$family, character(1L))),
    collapse = " + "
  )

  # The reference model is fitted by brms whatever fitted the base. That is
  # not a substitution: the parsed model is engine-independent, and the
  # reference is a re-emission of it with an extended fixed part and the
  # base's own priors carried through the brms prior translation.
  #
  # lgm_gate()'s twelfth check therefore does not run on the base. It used
  # to run on the base's IR whenever the base was INLA, so an INLA user
  # asking for a Laplace prior was told the prior was impossible on INLA --
  # on a fit that was never going to INLA. No public call reaches that check
  # now: fb_select() is the only verb that takes a shrinkage prior, and
  # nothing in R/ writes a penalty block onto the IR. It stays as an
  # internal guard.
  reference_engine <- "brms"
  approximations <- resolved$approximations
  if (!identical(engine_name, reference_engine)) {
    approximations <- c(approximations, paste0(
      "reference_fitted_on_brms: the base was fitted by ", engine_name,
      " and the reference model carrying the joint shrinkage prior is ",
      "re-emitted on brms with the base's structure and its own priors ",
      "carried through, because INLA does not fit a joint shrinkage prior: ",
      "a Laplace, horseshoe or R2D2 prior is Gaussian only given a local ",
      "scale per coefficient, and INLA's integration is built for a few ",
      "hyperparameters"
    ))
  }
  # The projection's own k-fold search refits the reference through brms,
  # which reuses the reference fit's data and compiled prior. The benchmark's
  # out-of-fold tables re-resolve both on each fold's analysis rows; the
  # search inside the provider cannot, and the record says so.
  if (identical(decision$method, "projection")) {
    approximations <- c(approximations, paste0(
      "projection_search_folds_reuse_whole_sample_preprocessing: the ",
      "projection's k-fold search refits the reference model through brms, ",
      "which reuses the reference fit's standardised data, compiled prior ",
      "and calibration, so the constants and prior scales inside the ",
      "search's folds are the whole sample's and the elpd-by-size curve ",
      "that chooses the submodel is not strictly out of fold. The ",
      "benchmark's out-of-fold tables learn every constant on each fold's ",
      "analysis rows"
    ))
  }

  if (isTRUE(plan)) {
    return(.fb_plan_selection(
      base = base,
      axes = list(
        workflow = "bayes",
        provider = "brms",
        prior = prior_family,
        decision = decision$method
      ),
      blocks = resolved$blocks,
      structure_gate = structure_gate,
      resamples = resamples,
      prediction = prediction,
      fits_per_fold = if (isTRUE(decision$validate_search)) 1L else 0L,
      candidates = candidates,
      compiled = resolved$compiled,
      approximations = approximations,
      base_engine = engine_name,
      reference_engine = reference_engine,
      decision_resolution = decision_resolution
    ))
  }

  # An explicit projection on a red gate refuses HERE, before the reference
  # model is sampled. It used to refuse inside the projection, after the
  # sampler had run: the same answer, several minutes later, about a
  # structure the gate could read from the parsed model before the first
  # fit. The gate inside the projection stays -- it reads the reference
  # model rather than the base, so it can still red a base that was green.
  if (identical(decision$method, "projection")) {
    .fb_sel_require_projectable(base$fb)
  }
  # The out-of-fold score `benchmark = TRUE` asks for relabels the grouping
  # levels under `new_cluster`; a term that cannot take a new level refuses
  # here, before the reference model is sampled.
  if (isTRUE(benchmark)) {
    .fb_sel_require_relabelable(base$fb, prediction)
  }
  # A formula variable whose value on a row depends on other rows would let a
  # fold's training rows read its assessment rows. The projection's k-fold
  # search and the out-of-fold score `benchmark = TRUE` asks for both refit
  # the formula on each fold, so it refuses here, before the reference model
  # is sampled.
  if (identical(decision$method, "projection") || isTRUE(benchmark)) {
    .fb_sel_require_fold_local_formula(
      list(protected, candidates),
      base$data,
      .fb_sel_folds(resamples, base$data, seed = seed)
    )
  }
  # The out-of-fold score `benchmark = TRUE` asks for compiles the reference
  # prior over the declared candidate block on every fold, so a fold whose
  # training rows leave a candidate unidentified refuses here, before the
  # reference model is sampled (R/selection_candidate_space.R).
  if (isTRUE(benchmark)) {
    .fb_sel_require_no_candidate_drop(.fb_sel_candidate_space(
      protected,
      candidates,
      base$data,
      .fb_sel_folds(resamples, base$data, seed = seed),
      groupings = unique(c(
        .fb_sel_grouping_columns(base$fb)$columns,
        resamples$unit
      )),
      route = "bayes"
    ))
  }

  cost <- .fb_sel_cost_estimate(
    base,
    resamples,
    fits_per_fold = if (isTRUE(decision$validate_search)) 1L else 0L,
    probe_seconds = .fb_select_base_seconds(base)
  )
  # One reference fit on top of the fold refits the projection search makes.
  cost$n_fits <- cost$n_fits + 1L
  cost$seconds <- cost$per_fit_seconds * cost$n_fits
  .fb_sel_bayes_check_budget(cost, budget)

  reference_args <- list(
    base = base,
    candidates = candidates,
    prior = resolved$compiled$prior,
    stanvars = resolved$compiled$stanvars,
    column_map = resolved$column_map,
    engine = engine,
    seed = seed,
    data = resolved$reference_data
  )
  # The sampler's step-size target is a lever, not a default to be argued
  # with: the contract's 0.95 is the floor, and the convergence gate refuses
  # on a single divergent transition, so a caller who hits one needs
  # somewhere to reach.
  adapt_delta <- (.fb_engine_opts(engine) %||% list())$adapt_delta
  if (!is.null(adapt_delta)) {
    reference_args$control <- list(adapt_delta = adapt_delta)
  }
  fit <- do.call(.fb_sel_reference_fit, reference_args)
  adequacy <- .fb_sel_require_adequate(fit, metric = metric, seed = seed)

  kappa <- .fb_sel_kappa(fit, resolved$blocks)

  projection <- NULL
  decision_result <- NULL
  folds <- .fb_sel_folds(resamples, base$data, seed = seed)
  if (identical(decision$method, "projection")) {
    # The search and the predictive evaluation run on CLUSTERED draws
    # rather than on the whole posterior. Left at projpred's own defaults
    # the search projects every draw onto every submodel at every step and
    # inside every fold, which on a 72-row mixed model with four candidates
    # had not finished in twenty minutes -- the projection of a multilevel
    # submodel is a fit per draw. These are the settings this release runs
    # at; they are recorded on the object so a reader can see what the elpd
    # column was computed from.
    projection <- .fb_sel_project(
      fit,
      folds = folds,
      protected_terms = protected_labels,
      decision = decision,
      seed = seed,
      nclusters = .FB_SEL_PROJECTION_NCLUSTERS,
      nclusters_pred = .FB_SEL_PROJECTION_NCLUSTERS_PRED
    )
    projection$settings$nclusters <- .FB_SEL_PROJECTION_NCLUSTERS
    projection$settings$nclusters_pred <- .FB_SEL_PROJECTION_NCLUSTERS_PRED
    decision_result <- projection
  } else if (identical(decision$method, "posterior_threshold")) {
    decision_result <- .fb_sel_posterior_threshold(
      fit,
      blocks = resolved$blocks,
      decision = decision,
      candidate_terms = resolved$candidate_labels,
      protected_terms = protected_labels
    )
  } else {
    decision_result <- .fb_sel_decision_full(
      fit,
      candidate_terms = resolved$candidate_labels,
      protected_terms = protected_labels
    )
  }

  # The primary arm's own out-of-fold score, on the folds the projection
  # used, so a benchmark placed beside it is placed beside something. It is
  # `benchmark = TRUE` only because it costs a reference refit per fold per
  # model; the certificate says which of the two situations a reader is in.
  bench <- NULL
  reference_cv <- NULL
  if (isTRUE(benchmark)) {
    reference_cv <- .fb_select_reference_cv(
      base = base,
      candidates = candidates,
      penalty = penalty,
      protected_labels = protected_labels,
      resamples = resamples,
      prediction = prediction,
      metric = metric,
      engine = engine,
      seed = seed,
      control = reference_args$control,
      # The decision stage searches each fold's submodel terms again on that
      # fold's reference fit, at the size the whole-sample decision chose and
      # with the whole-sample search's own settings.
      decision_search = if (is.null(projection)) {
        NULL
      } else {
        list(
          size = projection$chosen_size,
          forced = projection$forced_terms,
          method = projection$settings$search,
          nterms_max = projection$settings$nterms_max,
          nclusters = .FB_SEL_PROJECTION_NCLUSTERS
        )
      }
    )
    bench <- .fb_select_benchmark(
      base = base,
      candidates = candidates,
      protected = protected,
      resamples = resamples,
      prediction = prediction,
      seed = seed,
      metric = .fb_select_benchmark_metric(metric),
      prior_family = prior_family,
      blocks = resolved$blocks,
      decision = decision,
      shared_predictive_table = !is.null(reference_cv$aggregated)
    )
  }
  # Where the decision stage's submodel terms were searched, read off the
  # out-of-fold run: on each fold's training rows where one scored the
  # submodel, and nowhere where none did.
  if (!is.null(projection)) {
    projection$decision_stage_terms <- reference_cv$derived_terms
    projection$settings$decision_stage_search <- if (
      length(reference_cv$derived_terms) > 0L
    ) {
      "each_fold_training_rows"
    } else {
      "none"
    }
  }

  .fb_select_assemble(
    base = base,
    axes = list(
      workflow = "bayes",
      provider = "brms",
      prior = prior_family,
      decision = decision$method
    ),
    cv = list(
      # The folds the out-of-fold score was computed on where one ran, which
      # is what its score contract is read from.
      folds = reference_cv$folds %||% folds,
      design = reference_cv$design,
      # NULL unless `benchmark = TRUE` asked for the reference model to be
      # scored out of fold; the certificate's two predictive tables name
      # that as the reason where it is absent.
      aggregated = reference_cv$aggregated,
      clusters = reference_cv$clusters,
      rows = reference_cv$rows,
      preprocessing = reference_cv$preprocessing,
      candidate_space = reference_cv$candidate_space,
      metric = metric,
      prediction = prediction,
      seed = seed,
      covariance_estimand = "integrated"
    ),
    resamples = resamples,
    candidates = candidates,
    reference = fit,
    blocks = resolved$blocks,
    column_map = resolved$column_map,
    adequacy = adequacy,
    kappa = kappa,
    projection = projection,
    decision_result = decision_result,
    benchmark = bench,
    metric = metric,
    approximations = approximations,
    calibration = calibration,
    provider_versions = .fb_select_provider_versions("brms", decision),
    reference_engine = reference_engine,
    structure_gate = structure_gate,
    decision_resolution = decision_resolution,
    dots = dots
  )
}


# --- the optimiser workflow ---------------------------------------------- #

# .fb_select_optimizer() --- the penalised optimiser as the primary answer.
#
# Reached two ways, and the object says which: a deliberate
# `workflow = "penalized_optimizer"`, or `fallback = "glmnet"` firing because
# the Bayesian provider is absent. Neither is labelled Bayesian anywhere.
.fb_select_optimizer <- function(base,
                                 candidates,
                                 penalty,
                                 protected,
                                 protected_labels,
                                 resamples,
                                 prediction,
                                 decision,
                                 seed,
                                 metric,
                                 plan,
                                 fallback_fired,
                                 requested_workflow,
                                 calibration = NULL,
                                 structure_gate = NULL,
                                 decision_resolution = NULL) {
  resolved <- .fb_select_resolve_penalty(
    base,
    candidates,
    penalty,
    protected_labels,
    calibration = calibration,
    workflow = "penalized_optimizer"
  )
  prior_family <- paste(
    unique(vapply(resolved$blocks, function(b) b$family, character(1L))),
    collapse = " + "
  )
  optimizer_provider <- if (all(vapply(
    resolved$blocks,
    function(block) block$family %in% .fb_sel_grouped_families(),
    logical(1L)
  ))) {
    "prox"
  } else {
    "glmnet"
  }
  axes <- list(
    workflow = "penalized_optimizer",
    provider = optimizer_provider,
    prior = prior_family,
    decision = "optimizer_support"
  )
  if (isTRUE(plan)) {
    return(.fb_plan_selection(
      base = base,
      axes = axes,
      blocks = resolved$blocks,
      structure_gate = structure_gate,
      resamples = resamples,
      prediction = prediction,
      fits_per_fold = 0L,
      candidates = candidates,
      compiled = resolved$compiled,
      approximations = resolved$approximations,
      decision_resolution = decision_resolution
    ))
  }
  bench <- .fb_select_benchmark(
    base = base,
    candidates = candidates,
    protected = protected,
    resamples = resamples,
    prediction = prediction,
    seed = seed,
    metric = .fb_select_benchmark_metric(metric),
    prior_family = prior_family,
    blocks = resolved$blocks,
    decision = decision,
    shared_predictive_table = TRUE
  )
  if (is.null(bench$arms$glmnet)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_provider_absent",
      message = paste0(
        "fb_select(): the penalised-optimiser workflow produced no ",
        "optimiser arm. ",
        paste(bench$refusals, collapse = " ")
      )
    ))
  }

  .fb_select_assemble(
    base = base,
    axes = axes,
    cv = bench$cv,
    resamples = resamples,
    candidates = candidates,
    blocks = resolved$blocks,
    column_map = resolved$column_map,
    decision_result = list(
      decision = "optimizer_support",
      chosen_terms = bench$arms$glmnet$selection$selected_terms %||%
        character(0),
      chosen_vocabulary = "term",
      chosen_size = length(
        bench$arms$glmnet$selection$selected_terms %||% character(0)
      ),
      label = paste0(
        "the support of the penalised optimiser at its chosen lambda; a ",
        "mode, not a posterior"
      ),
      # Per penalised column: the folds that estimated it, beside the folds
      # whose support held it at the chosen point.
      candidate_availability = bench$arms$glmnet$selection$candidate_availability,
      settings = list(method = "optimizer_support")
    ),
    benchmark = bench,
    metric = metric,
    approximations = resolved$approximations,
    refusals = bench$refusals,
    calibration = calibration,
    provider_versions = .fb_select_provider_versions(optimizer_provider),
    fallback_fired = fallback_fired,
    requested_workflow = requested_workflow,
    structure_gate = structure_gate,
    decision_resolution = decision_resolution,
    dots = list()
  )
}


# .fb_select_assessment_predictions() --- the optimiser decision's scored
# rows, surfaced without rebuilding or predicting on data supplied later.
#
# The optimiser already computed these values inside each fold; earlier
# versions retained them inside the canonical path and discarded them while
# assembling the public result. Keeping that table makes the score auditable
# at row level. It does not create a `predict.fb_select()` method and does not
# license prediction on new rows: the table is restricted to assessment rows
# the recorded fold contract identifies.
.fb_select_empty_assessment_predictions <- function() {
  data.frame(
    stage = character(0), workflow = character(0),
    provider = character(0), model = character(0),
    fold = character(0), repeat_id = integer(0), row = integer(0),
    row_key = character(0), unit = character(0),
    observed_model_scale = numeric(0), predicted_model_scale = numeric(0),
    residual = numeric(0), selected_lambda_index = integer(0),
    selected_s = numeric(0), lambda_canonical = numeric(0),
    prediction = character(0), response = character(0),
    contract_hash = character(0), stringsAsFactors = FALSE
  )
}

.fb_select_assessment_predictions <- function(base,
                                              axes,
                                              benchmark,
                                              resamples,
                                              metric) {
  arm <- benchmark$arms$glmnet
  raw <- arm$assessment_predictions
  selection <- arm$selection
  if (is.null(arm) || !isTRUE(selection$ran) ||
      !is.data.frame(raw) || nrow(raw) == 0L) {
    return(.fb_select_empty_assessment_predictions())
  }
  index <- selection$selected_lambda_index
  raw <- raw[
    raw$arm == "benchmark" & raw$lambda_index == index, , drop = FALSE
  ]
  if (nrow(raw) == 0L) {
    return(.fb_select_empty_assessment_predictions())
  }
  contract <- .fb_select_score_contract(
    folds = benchmark$cv$folds,
    data = base$data,
    response = base$fb$response,
    prediction = benchmark$cv$prediction,
    metric = metric,
    design = benchmark$cv$design,
    learned_on = benchmark$cv$preprocessing$learned_on %||% NA_character_,
    stage = "decision",
    scoring_unit = arm$scoring_unit,
    offsets = base$fb$offsets %||% character(0),
    candidate_space = benchmark$cv$candidate_space
  )
  .fb_select_verify_contract(contract, "optimizer assessment predictions")
  out <- data.frame(
    stage = "decision",
    workflow = axes$workflow,
    provider = arm$provider,
    model = selection$selected_model,
    fold = as.character(raw$fold),
    repeat_id = as.integer(raw$repeat_id),
    row = as.integer(raw$row),
    row_key = paste0("repeat:", raw$repeat_id, "/row:", raw$row),
    unit = as.character(raw$unit),
    observed_model_scale = as.numeric(raw$observed),
    predicted_model_scale = as.numeric(raw$fitted),
    residual = as.numeric(raw$residual),
    selected_lambda_index = as.integer(index),
    selected_s = as.numeric(selection$selected_s),
    lambda_canonical = as.numeric(raw$lambda_canonical),
    prediction = benchmark$cv$prediction,
    response = as.character(base$fb$response),
    contract_hash = contract$hash,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$repeat_id, out$row, method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}


# --- assembly ------------------------------------------------------------ #

.fb_select_assemble <- function(base,
                                axes,
                                cv,
                                resamples,
                                candidates,
                                reference = NULL,
                                blocks = NULL,
                                column_map = NULL,
                                adequacy = NULL,
                                kappa = NULL,
                                projection = NULL,
                                decision_result = NULL,
                                benchmark = NULL,
                                metric = "elpd",
                                approximations = character(0),
                                refusals = character(0),
                                calibration = NULL,
                                provider_versions = list(),
                                fallback_fired = FALSE,
                                requested_workflow = NULL,
                                reference_engine = NA_character_,
                                structure_gate = NULL,
                                decision_resolution = NULL,
                                dots = list()) {
  increment <- NULL
  agg <- cv$aggregated
  if (is.data.frame(agg) &&
      all(c("protected_only", "full") %in% agg$model)) {
    base_row <- agg[agg$model == "protected_only", ]
    full_row <- agg[agg$model == "full", ]
    increment <- list(
      protected_only = base_row$score,
      full = full_row$score,
      # Every metric is carried as a LOSS, so the increment reads the same
      # way whichever one was asked for: positive is better.
      improvement = base_row$score - full_row$score,
      metric = cv$metric %||% metric,
      se_unit = full_row$se_unit,
      n_units = full_row$n_units
    )
  }

  stages <- .fb_select_stage_tables(
    base = base,
    axes = axes,
    cv = cv,
    benchmark = benchmark,
    metric = .fb_select_benchmark_metric(cv$metric %||% metric)
  )

  predictions <- .fb_select_assessment_predictions(
    base = base,
    axes = axes,
    benchmark = benchmark,
    resamples = resamples,
    metric = .fb_select_benchmark_metric(cv$metric %||% metric)
  )

  certificate <- .fb_select_certificate(
    base = base,
    axes = axes,
    cv = cv,
    resamples = resamples,
    candidates = candidates,
    blocks = blocks,
    adequacy = adequacy,
    projection = projection,
    benchmark = benchmark$arms,
    fallback_fired = fallback_fired,
    approximations = approximations,
    refusals = c(refusals, benchmark$refusals %||% character(0)),
    calibration = calibration,
    provider_versions = provider_versions,
    column_map = column_map,
    reference_engine = reference_engine,
    reference_predictive = stages$reference,
    decision_predictive = stages$decision,
    structure_gate = structure_gate,
    decision_resolution = decision_resolution,
    fold_preprocessing = .fb_select_fold_preprocessing(cv),
    # The primary arm's folds where it compiled any, and the benchmark's
    # otherwise: the same rule gives both arms the same space.
    candidate_space = cv$candidate_space %||% benchmark$cv$candidate_space
  )

  structure(
    list(
      workflow = axes$workflow,
      provider = axes$provider,
      prior = axes$prior,
      decision = axes$decision,
      requested_workflow = requested_workflow %||% axes$workflow,
      requested_decision = decision_resolution$requested %||% axes$decision,
      fallback_fired = isTRUE(fallback_fired),
      structure_gate = structure_gate,
      base = base,
      reference = reference,
      column_map = column_map,
      blocks = blocks,
      projection = projection,
      selection = decision_result,
      scores = agg,
      clusters = cv$clusters,
      assessment_predictions = predictions,
      increment = increment,
      diagnostics = list(
        adequacy = adequacy,
        kappa = kappa,
        rank_stability = projection$rank_stability
      ),
      benchmark = benchmark$arms,
      folds = cv$folds,
      certificate = certificate
    ),
    class = "fb_select"
  )
}


# --- small helpers -------------------------------------------------------- #

# The clustered-draw settings the projection search runs at. They are
# constants rather than arguments in this release: they trade a small amount
# of Monte-Carlo error in the elpd column for a search that finishes, and a
# caller who wants the whole posterior projected can call the internal
# search directly. Both travel onto the object and the certificate.
.FB_SEL_PROJECTION_NCLUSTERS <- 20L
.FB_SEL_PROJECTION_NCLUSTERS_PRED <- 100L

# The two arguments the earlier signature carried. `metric` moved into
# fb_decision(); the covariance summary is a property of the benchmark
# workflow and is now read off the folds rather than declared, which is the
# defect the contract's section 9.2 closes. Both are accepted for one release
# so a call written against the earlier signature reports what changed
# instead of failing on an unused argument.
.fb_select_absorb_deprecated <- function(dots) {
  for (nm in c("metric", "covariance")) {
    if (nm %in% names(dots)) {
      .fb_warn_deprecated(c(
        paste0(
          "The `", nm, "` argument of `fb_select()` is deprecated as of ",
          "flexyBayes 0.11.0."
        ),
        if (identical(nm, "metric")) {
          "Use the `metric` argument of `fb_decision()`."
        } else {
          paste0(
            "The covariance summary is a property of the benchmark ",
            "workflow and is read from the fold object that used it, so a ",
            "declared value could disagree with what ran."
          )
        }
      ))
    }
  }
  unknown <- setdiff(names(dots), c("metric", "covariance"))
  if (length(unknown) > 0L) {
    stop(
      "fb_select(): unused argument(s) ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  dots
}

.fb_select_check_decision <- function(decision) {
  if (inherits(decision, "fb_decision")) {
    return(decision)
  }
  if (is.list(decision)) {
    return(do.call(fb_decision, decision))
  }
  stop(.fb_refusal_condition(
    reason_code = "selection_decision_invalid",
    message = paste0(
      "fb_select(): `decision` must be an <fb_decision>; see ?fb_decision."
    )
  ))
}

# The benchmark arms are closed-form and score a loss, and `elpd` is not one
# of the three losses they compute. `nlpd` is the same quantity with the sign
# a loss carries, so it is what the arms are asked for and what the arm's own
# record says they scored.
.fb_select_benchmark_metric <- function(metric) {
  if (identical(metric, "elpd")) "nlpd" else metric
}

.fb_select_provider_versions <- function(provider, decision = NULL) {
  out <- list()
  packages <- switch(
    provider,
    brms = c("brms", "posterior"),
    inla = "INLA",
    glmnet = "glmnet",
    character(0)
  )
  if (!is.null(decision) &&
      identical(decision$method %||% "", "projection")) {
    packages <- c(packages, "projpred")
  }
  for (p in packages) {
    if (.fb_sel_provider_installed(p)) {
      out[[p]] <- as.character(utils::packageVersion(p))
    }
  }
  out
}

# The base fit's own wall time, as the per-fit proxy for a cost estimate that
# must not itself fit anything. It is a proxy and is labelled one: the
# reference model carries a wider fixed part and a shrinkage prior, so it is
# the floor of the eventual cost rather than the estimate of it.
.fb_select_base_seconds <- function(base) {
  seconds <- suppressWarnings(as.numeric(
    base$extras$run_time %||% NA_real_
  ))
  if (length(seconds) != 1L || !is.finite(seconds) || seconds <= 0) {
    return(1)
  }
  seconds
}


# --- methods -------------------------------------------------------------- #

#' @export
print.fb_select <- function(x, digits = 4L, ...) {
  cert <- x$certificate
  if (isTRUE(x$fallback_fired)) {
    cat(
      "<fb_select> SUBSTITUTED: the Bayesian provider was absent, so the ",
      "penalised optimiser ran under its own name.\n",
      sep = ""
    )
  }
  cat(
    "<fb_select> workflow ", x$workflow,
    " / provider ", x$provider,
    " / prior ", x$prior,
    " / decision ", x$decision,
    if (identical(cert$decision_resolved_from, "auto")) {
      " (resolved from \"auto\")"
    } else {
      ""
    },
    if (identical(x$decision, "projection")) {
      paste0(
        ", searched on ", .FB_SEL_PROJECTION_NCLUSTERS,
        " clustered draws and evaluated on ",
        .FB_SEL_PROJECTION_NCLUSTERS_PRED
      )
    } else {
      ""
    },
    "\n",
    sep = ""
  )
  # The gate leads, because it is what decides whether the default decision
  # can run at all and most of what this package fits is red at it.
  .fb_print_structure_gate(x$structure_gate)
  if (identical(cert$decision_resolved_from, "auto")) {
    cat("  decision  : ", .fb_squish(cert$decision_resolution_reason), "\n",
        sep = "")
  }
  cat(
    "  base      : ", cert$base$family, "/", cert$base$link,
    ", ", cert$base$n_observed, " rows, convergence ",
    if (is.na(cert$base$convergence)) "UNVERIFIED"
    else if (cert$base$convergence) "ok" else "FAILED",
    "\n",
    sep = ""
  )
  # Printed only where it says something the row count does not, and only
  # where a prior scale was actually calibrated against it: in a clustered
  # design the scales are stated at the unit, and a reader who saw only the
  # row count would read a different prior.
  if (length(cert$prior_blocks) > 0L &&
        !identical(cert$calibration_n, cert$n_rows)) {
    cat(
      "  prior n   : calibrated at ", cert$calibration_n, " ",
      paste(cert$calibration_unit, collapse = ", "),
      ", not at ", cert$n_rows, " rows\n",
      sep = ""
    )
  }
  cat(
    "  covariance: ", cert$covariance$builder %||% "unknown builder", " / ",
    cert$covariance$estimand %||% "not recorded",
    if (isTRUE(cert$covariance$integrated)) {
      " (integrated by the engine, not plugged in)"
    } else {
      ""
    },
    "\n",
    sep = ""
  )
  cat("  resampling: ", cert$resampling$n_folds, " folds, prediction ",
      cert$resampling$prediction,
      if (is.na(cert$resampling$scoring_unit %||% NA_character_)) {
        ", no out-of-fold score on this route"
      } else {
        paste0(", scored at ", cert$resampling$scoring_unit)
      },
      "\n", sep = "")

  base_engine <- if (inherits(x$base, "flexybayes_inla")) {
    "inla"
  } else if (inherits(x$base, "flexybayes_brms")) {
    "brms"
  } else {
    NA_character_
  }
  if (!is.na(cert$reference_engine %||% NA_character_) &&
        !identical(cert$reference_engine, base_engine)) {
    cat(
      "  reference : fitted by ", cert$reference_engine,
      ", while the base was fitted by ", base_engine,
      " (the parsed model is engine-independent; the base's priors are ",
      "carried through)\n",
      sep = ""
    )
  }
  if (!is.null(x$diagnostics$adequacy)) {
    ad <- x$diagnostics$adequacy
    cat("  reference : adequacy ", ad$verdict, sep = "")
    if (length(ad$reasons) > 0L) {
      cat(" -- ", paste(ad$reasons, collapse = "; "), sep = "")
    }
    cat("\n")
  }

  if (!is.null(x$selection)) {
    cat("\n  decision  : ", x$selection$label %||% x$decision, "\n", sep = "")
    if (!is.null(x$selection$chosen_size)) {
      cat("    size    : ", x$selection$chosen_size, sep = "")
      if (length(x$selection$chosen_terms) > 0L) {
        cat("  (", paste(x$selection$chosen_terms, collapse = ", "), ")",
            sep = "")
      }
      cat("\n")
    }
    if (isTRUE(x$projection$settings$optimistic)) {
      cat(
        "    NOTE    : validate_search = FALSE, so the elpd curve is ",
        "optimistic\n",
        sep = ""
      )
    }
  }

  if (!is.null(x$increment)) {
    cat("\n  out-of-fold ", x$increment$metric, ":\n", sep = "")
    cat("    protected only : ", round(x$increment$protected_only, digits),
        "\n", sep = "")
    cat("    with candidates: ", round(x$increment$full, digits), "\n",
        sep = "")
    cat("    improvement    : ", round(x$increment$improvement, digits),
        "  (", x$increment$n_units, " independent ",
        x$increment$se_unit, ")\n", sep = "")
  }

  if (length(x$benchmark) > 0L) {
    cat("\n  benchmark arms:\n", sep = "")
    for (nm in names(x$benchmark)) {
      arm <- x$benchmark[[nm]]
      cat(
        "    ", format(nm, width = 8), " ", arm$workflow, " / ",
        arm$provider, "  [", arm$equivalence_level$level, "]\n",
        sep = ""
      )
    }
    for (ps in list(cert$reference_predictive, cert$decision_predictive)) {
      if (isTRUE(ps$available)) {
        cat("    ", ps$stage, " stage, out-of-fold ", ps$metric,
            " at the resampling unit:\n", sep = "")
        for (i in seq_len(nrow(ps$arms))) {
          cat("      ", format(ps$arms$arm[[i]], width = 8), " ",
              format(ps$arms$model[[i]], width = 18), " ",
              format(round(ps$arms$score[[i]], digits), nsmall = 2L),
              "  (", ps$arms$n_units[[i]], " ", ps$arms$unit[[i]], ")\n",
              sep = "")
        }
      } else {
        cat("    ", ps$stage, " stage: none -- ", .fb_squish(ps$reason),
            "\n", sep = "")
      }
    }
    mo <- cert$mode_oracle
    cat("    mode oracle: ",
        if (isTRUE(mo$available)) {
          paste(mo$arms, collapse = ", ")
        } else {
          "none"
        },
        " -- ", .fb_squish(mo$reason), "\n", sep = "")
  }

  if (length(cert$approximations) > 0L) {
    cat("\n  approximations:\n", sep = "")
    for (a in cert$approximations) {
      cat("    - ", .fb_squish(a), "\n", sep = "")
    }
  }
  if (length(cert$refusals) > 0L) {
    cat("\n  refusals carried:\n", sep = "")
    for (r in cert$refusals) {
      cat("    - ", .fb_squish(r), "\n", sep = "")
    }
  }
  invisible(x)
}


#' The penalised optimiser's support history, as the accessors read it
#'
#' The `glmnet` arm records what its path did with every declared candidate
#' column at every point of the tuning grid. Both tables are empty where no
#' such arm ran, and where the family has no support to report.
#'
#' @param x An `fb_select` object.
#' @return `list(points, by_fold)`: the membership at each point of each
#'   fold, and one row per fold and column carrying that column's entry, exit
#'   and re-entries.
#' @noRd
#' @keywords internal
.fb_select_support_path <- function(x) {
  support <- x$benchmark$glmnet$support
  if (is.null(support$history) || is.null(support$history_by_fold)) {
    return(.fb_sel_glmnet_support_history(
      supports = NULL,
      penalised_by_fold = list(),
      declared_penalised = character(0),
      alpha = NA_real_,
      selects = FALSE
    ))
  }
  list(points = support$history, by_fold = support$history_by_fold)
}


#' Tidy a selection into one row per candidate coefficient
#'
#' @details
#' Retention is read at the vocabulary the decision actually chose in. A
#' projection, a `full` decision and an optimiser support all choose FORMULA
#' TERMS -- `f`, `x1:f` -- while the design carries contrast columns --
#' `fB`, `x1:fB`. Matching a chosen term against a column name reports every
#' factor and interaction column as dropped however the decision went, so
#' the term recorded on the column map is what the match is made against. A
#' posterior threshold is the one decision that chooses coefficient by
#' coefficient, and its retention is matched on the column.
#'
#' Where a penalised optimiser ran -- as the workflow, or as the `glmnet` arm
#' of `benchmark = TRUE` -- each row also carries what that path did with the
#' column over its lambda grid, on the coordinate the grid was built on:
#' `s = lambda / lambda_max`, with each fold's `lambda_max` read on its own
#' training rows. `folds_available` counts the folds whose training rows
#' identified the column at all, which is a different statement from the
#' folds that never held it: a column a fold could not estimate was never on
#' that fold's path. `path_entry_s` is the largest `s` at which any fold
#' holds the column and `path_exit_s` the smallest `s` at which a fold that
#' held it has let it go, so a column retained across the whole path has an
#' entry and no exit. `path_reentries` is the largest number of returns in
#' any one fold. Those summaries are taken from the attribute
#' `support_path_by_fold`, one row per fold and column, and the membership at
#' every point of every fold is the attribute `support_path`.
#'
#' @param x An `fb_select` object.
#' @param ... Ignored.
#' @returns A `data.frame` with one row per design column: the column, the
#'   formula term it came from, its role, its block and family, whether the
#'   decision retained it, and, where a reference posterior was fitted, its
#'   posterior mean and standard deviation on the scientific scale. Beside
#'   them, the penalised path's `folds_available`, `folds_entered`,
#'   `path_entry_s`, `path_exit_s` and `path_reentries`, all `NA` for a
#'   protected column and on a run with no penalised path, with the path
#'   itself in the attributes `support_path` and `support_path_by_fold`.
#' @export
tidy.fb_select <- function(x, ...) {
  map <- x$column_map
  if (is.null(map)) {
    map <- data.frame(
      column = character(0), term = character(0), role = character(0),
      block = integer(0), family = character(0), centre = numeric(0),
      scale = numeric(0), stringsAsFactors = FALSE
    )
  }
  if (is.null(map$term)) {
    map$term <- map$column
  }
  chosen <- x$selection$chosen_terms %||% character(0)
  # The decision says which vocabulary its chosen set is written in; only
  # the posterior threshold chooses coefficients, and everything else
  # chooses terms. An absent field is a term-level decision, which is what
  # every route that predates the field was.
  vocabulary <- x$selection$chosen_vocabulary %||% "term"
  retained <- if (identical(vocabulary, "coefficient")) {
    map$column %in% chosen
  } else {
    map$term %in% chosen
  }
  n <- nrow(map)
  out <- data.frame(
    term = map$column,
    source_term = map$term,
    role = map$role,
    block = map$block,
    family = map$family,
    retained = retained | map$role == "protected",
    estimate = rep(NA_real_, n),
    std.error = rep(NA_real_, n),
    scale = rep("scientific", n),
    stringsAsFactors = FALSE
  )

  brmsfit <- x$reference$brms
  if (!is.null(brmsfit) && !is.null(x$blocks)) {
    draws <- as.matrix(brmsfit)
    keep <- colnames(draws) %in% c(paste0("b_", map$column), "b_Intercept")
    scaled <- .fb_pen_back_transform(draws[, keep, drop = FALSE], x$blocks)
    for (i in seq_len(nrow(out))) {
      column <- paste0("b_", out$term[[i]])
      if (column %in% colnames(scaled)) {
        v <- as.numeric(scaled[, column])
        out$estimate[[i]] <- mean(v)
        out$std.error[[i]] <- stats::sd(v)
      }
    }
  }
  # Where the projection ran, say which draws a consumer is looking at: the
  # estimates above are the REFERENCE posterior, never the projected one.
  attr(out, "draw_source") <- if (is.null(brmsfit)) {
    "none"
  } else {
    "reference_posterior"
  }

  # The penalised path, summarised per column over the folds. A protected
  # column is on no penalised path and a run with no optimiser has none at
  # all; both read NA rather than zero, which would say the path was walked
  # and the column never entered.
  path <- .fb_select_support_path(x)
  by_fold <- path$by_fold
  out$folds_available <- rep(NA_integer_, n)
  out$folds_entered <- rep(NA_integer_, n)
  out$path_entry_s <- rep(NA_real_, n)
  out$path_exit_s <- rep(NA_real_, n)
  out$path_reentries <- rep(NA_integer_, n)
  for (i in seq_len(n)) {
    rows <- by_fold[by_fold$column == out$term[[i]], , drop = FALSE]
    if (nrow(rows) == 0L) {
      next
    }
    entered <- rows[rows$n_in_support > 0L, , drop = FALSE]
    out$folds_available[[i]] <- sum(rows$available)
    out$folds_entered[[i]] <- nrow(entered)
    if (nrow(entered) > 0L) {
      out$path_entry_s[[i]] <- max(entered$first_entry_s)
      if (any(!is.na(entered$last_exit_s))) {
        out$path_exit_s[[i]] <- min(entered$last_exit_s, na.rm = TRUE)
      }
      out$path_reentries[[i]] <- max(entered$n_reentries)
    }
  }
  attr(out, "support_path") <- path$points
  attr(out, "support_path_by_fold") <- by_fold
  out
}


#' @export
glance.fb_select <- function(x, ...) {
  arm <- x$benchmark$glmnet
  path <- .fb_select_support_path(x)
  by_lambda <- if (is.null(arm$support$by_lambda)) {
    NULL
  } else {
    as.data.frame(arm$support$by_lambda)
  }
  # The grid the support history walks, which is the benchmark arm: the fixed
  # arm is one canonical point rather than a path.
  indices <- if (is.null(by_lambda)) {
    integer(0)
  } else {
    by_lambda$lambda_index[by_lambda$arm == "benchmark"]
  }
  entered <- if (nrow(path$by_fold) == 0L) {
    NA_integer_
  } else {
    sum(tapply(
      path$by_fold$n_in_support, path$by_fold$column,
      function(v) any(v > 0L)
    ))
  }
  data.frame(
    workflow = x$workflow,
    provider = x$provider,
    prior = x$prior,
    decision = x$decision,
    metric = x$increment$metric %||% x$certificate$resampling$metric %||%
      NA_character_,
    protected_only = x$increment$protected_only %||% NA_real_,
    full = x$increment$full %||% NA_real_,
    improvement = x$increment$improvement %||% NA_real_,
    chosen_size = x$selection$chosen_size %||% NA_integer_,
    n_units = x$increment$n_units %||% NA_integer_,
    se_unit = x$increment$se_unit %||% NA_character_,
    n_folds = length(x$folds),
    prediction = x$certificate$resampling$prediction,
    covariance = x$certificate$covariance$estimand %||% NA_character_,
    builder = x$certificate$covariance$builder %||% NA_character_,
    fallback_fired = isTRUE(x$fallback_fired),
    # The penalised optimiser's tuning report and the size of the path its
    # support history walks. A run with no such arm carries them as absent.
    tuning_scale = arm$selection$tuning_scale %||% NA_character_,
    penalty_alpha = arm$support$alpha %||% NA_real_,
    selected_s = arm$selection$selected_s %||% NA_real_,
    path_points = if (length(indices) == 0L) {
      NA_integer_
    } else {
      length(unique(indices))
    },
    path_columns_entered = entered,
    stringsAsFactors = FALSE
  )
}


#' Out-of-fold results from a selection
#'
#' @param x An `fb_select` object.
#' @param data Ignored; the out-of-fold results are the ones the folds
#'   computed and cannot be recomputed on other rows.
#' @param type `"unit_scores"` (the default retained for compatibility) or
#'   `"assessment_predictions"`. The latter returns the optimiser decision's
#'   row-level assessment predictions, with original row positions, repeat and
#'   fold identifiers, scoring units, selected tuning coordinates and the
#'   verified score-contract hash. Observed and predicted values are on the
#'   modelled response scale, after any Gaussian offset has been subtracted.
#' @param ... Ignored.
#' @returns A `data.frame` of per-unit out-of-fold scores or row-level
#'   assessment predictions, according to `type`. A zero-row frame with the
#'   stable schema is returned where that result did not run. These are stored
#'   assessment results, not a prediction method for caller-supplied rows.
#' @export
augment.fb_select <- function(
  x,
  data = NULL,
  type = c("unit_scores", "assessment_predictions"),
  ...
) {
  type <- match.arg(type)
  if (identical(type, "assessment_predictions")) {
    out <- x$assessment_predictions
    if (is.null(out)) {
      out <- .fb_select_empty_assessment_predictions()
    }
    return(as.data.frame(out))
  }
  out <- x$clusters
  if (is.null(out)) {
    return(data.frame(
      model = character(0), unit = character(0), score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  as.data.frame(out)
}


#' Displays for a selection
#'
#' @param x An `fb_select` object.
#' @param type `"elpd_by_size"`, `"rank_stability"`, `"shrinkage"`,
#'   `"benchmark"` or `"entrance_path"`. The last plots the common
#'   out-of-fold loss against `s = lambda / lambda_max` and labels each
#'   candidate at its median first entry: the point on `s` by which at least
#'   half of the folds that could estimate it had admitted it. A candidate
#'   admitted on fewer than half of those folds is not labelled, and
#'   candidates entering at the same point share one label. The ordinate
#'   is the mean over scoring units of each unit's score, not an inclusion
#'   probability. It is available for optimiser workflows and optimiser
#'   benchmark arms whose penalty has a finite null anchor, which excludes
#'   `ridge()` and `exclusive_lasso()`.
#' @param ... Passed to the underlying plotting call.
#' @returns Invisibly, the data the display was drawn from. Called for the
#'   plot it draws.
#' @export
plot.fb_select <- function(x,
                           type = c("elpd_by_size", "rank_stability",
                                    "shrinkage", "benchmark",
                                    "entrance_path"),
                           ...) {
  type <- match.arg(type)
  switch(
    type,
    elpd_by_size = .fb_select_plot_elpd(x, ...),
    rank_stability = .fb_select_plot_stability(x, ...),
    shrinkage = .fb_select_plot_shrinkage(x, ...),
    benchmark = .fb_select_plot_benchmark(x, ...),
    entrance_path = .fb_select_plot_entrance(x, ...)
  )
}

.fb_select_plot_unavailable <- function(what, why) {
  message("plot.fb_select(): ", what, " is not available -- ", why, ".")
  invisible(NULL)
}

.fb_select_plot_elpd <- function(x, ...) {
  tab <- x$projection$elpd_by_size
  if (!is.data.frame(tab) || nrow(tab) == 0L) {
    return(.fb_select_plot_unavailable(
      "the elpd curve",
      "no projection search ran, so there is no elpd by submodel size"
    ))
  }
  graphics::plot(
    tab$size, tab$elpd,
    type = "b", pch = 16L,
    xlab = "submodel size", ylab = "elpd (unit-aggregated)",
    main = "Projection: predictive accuracy by submodel size",
    ...
  )
  if ("se_row_level" %in% names(tab)) {
    ok <- is.finite(tab$se_row_level)
    if (any(ok)) {
      graphics::arrows(
        tab$size[ok], tab$elpd[ok] - tab$se_row_level[ok],
        tab$size[ok], tab$elpd[ok] + tab$se_row_level[ok],
        length = 0.03, angle = 90, code = 3L
      )
    }
  }
  if (!is.null(x$selection$chosen_size)) {
    graphics::abline(v = x$selection$chosen_size, lty = 2L)
  }
  invisible(tab)
}

.fb_select_plot_stability <- function(x, ...) {
  stability <- x$projection$rank_stability
  if (is.null(stability)) {
    return(.fb_select_plot_unavailable(
      "the rank-stability display",
      "no fold-wise ranking was returned by the projection search"
    ))
  }
  m <- as.matrix(stability)
  graphics::image(
    seq_len(ncol(m)), seq_len(nrow(m)), t(m),
    xlab = "rank position", ylab = "term", axes = FALSE,
    main = "Fold-wise ranking agreement",
    ...
  )
  graphics::axis(1L, at = seq_len(ncol(m)), labels = colnames(m))
  graphics::axis(
    2L, at = seq_len(nrow(m)), labels = rownames(m), las = 1L, cex.axis = 0.7
  )
  graphics::box()
  invisible(stability)
}

.fb_select_plot_shrinkage <- function(x, ...) {
  kappa <- x$diagnostics$kappa
  if (is.null(kappa) || !isTRUE(kappa$available)) {
    return(.fb_select_plot_unavailable(
      "the shrinkage display",
      kappa$reason %||% "no shrinkage factor was computed for this prior"
    ))
  }
  ord <- order(kappa$kappa)
  graphics::barplot(
    kappa$kappa[ord],
    horiz = TRUE, las = 1L, xlim = c(0, 1),
    xlab = "shrinkage factor (orthogonal-design approximation)",
    main = "Candidate shrinkage",
    ...
  )
  graphics::abline(v = 0.5, lty = 2L)
  invisible(kappa)
}

.fb_select_plot_benchmark <- function(x, ...) {
  arms <- x$benchmark
  if (length(arms) == 0L) {
    return(.fb_select_plot_unavailable(
      "the benchmark display",
      "no benchmark arm ran; pass benchmark = TRUE"
    ))
  }
  labels <- vapply(
    names(arms),
    function(nm) paste0(nm, "\n[", arms[[nm]]$equivalence_level$level, "]"),
    character(1L)
  )
  scores <- vapply(
    arms,
    function(arm) {
      agg <- arm$scores
      if (is.data.frame(agg) && "full" %in% agg$model) {
        as.numeric(agg$score[agg$model == "full"][[1L]])
      } else {
        NA_real_
      }
    },
    numeric(1L)
  )
  graphics::barplot(
    scores,
    names.arg = labels, las = 1L,
    ylab = "out-of-fold loss on the common folds",
    main = "Benchmark arms, with the equivalence level of each comparison",
    ...
  )
  invisible(data.frame(
    arm = names(arms), score = unname(scores),
    equivalence_level = vapply(
      arms, function(a) a$equivalence_level$level, character(1L)
    ),
    stringsAsFactors = FALSE
  ))
}

# .fb_select_median_entry() --- where each candidate typically enters the
# path, read off the per-fold support history.
#
# The folds counted for a column are the ones whose training rows could
# estimate it. Its median first entry is the point on s by which at least
# half of those folds had admitted it: the ceiling(k / 2)-th largest first
# entry over k folds, which is the median for odd k and the earlier of the
# two middle entries for even k, so the label sits on a grid point. A fold on
# which the column never entered counts against it rather than dropping out,
# so a column admitted on one fold in five is not placed as if that fold were
# typical; a column admitted on fewer than half of its folds is not labelled.
.fb_select_median_entry <- function(history) {
  available <- history[history$available, , drop = FALSE]
  rows <- lapply(unique(available$column), function(column) {
    part <- available[available$column == column, , drop = FALSE]
    entered <- part$n_in_support > 0L & is.finite(part$first_entry_s)
    first <- sort(part$first_entry_s[entered], decreasing = TRUE)
    needed <- ceiling(nrow(part) / 2)
    data.frame(
      variable = column,
      entry_s = if (length(first) >= needed) first[[needed]] else NA_real_,
      n_folds_available = nrow(part),
      n_folds_entered = sum(entered),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    return(data.frame(
      variable = character(0), entry_s = numeric(0),
      n_folds_available = integer(0), n_folds_entered = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  out[is.finite(out$entry_s), , drop = FALSE]
}

.fb_select_plot_entrance <- function(x, max_labels = 15L,
                                     label_cex = 0.7, ...) {
  arm <- x$benchmark$glmnet
  scores <- arm$selection$scores
  history <- arm$support$history_by_fold
  if (!is.data.frame(scores) || nrow(scores) == 0L) {
    return(.fb_select_plot_unavailable(
      "the entrance path",
      "no penalised optimiser path was scored out of fold"
    ))
  }
  # Ridge and the exclusive lasso reach the null only in the limit, so their
  # grid is walked on the penalty itself and carries no s coordinate.
  if (!any(is.finite(scores$s))) {
    return(.fb_select_plot_unavailable(
      "the entrance path",
      paste0(
        "this penalty has no finite null anchor, so its path is indexed by `",
        arm$selection$tuning_scale %||% "the penalty",
        "` and has no s = lambda / lambda_max coordinate"
      )
    ))
  }
  if (!is.data.frame(history) || nrow(history) == 0L) {
    return(.fb_select_plot_unavailable(
      "the entrance path",
      "this penalty has no finite coefficient-support path"
    ))
  }
  entries <- .fb_select_median_entry(history)
  if (nrow(entries) == 0L) {
    return(.fb_select_plot_unavailable(
      "the entrance path",
      paste0(
        "no candidate entered on at least half of the folds that could ",
        "estimate it"
      )
    ))
  }
  entries <- entries[order(entries$entry_s, decreasing = TRUE), , drop = FALSE]
  max_labels <- as.integer(max_labels)
  if (!is.finite(max_labels) || max_labels < 0L) {
    stop("`max_labels` must be a non-negative integer.", call. = FALSE)
  }
  entries <- utils::head(entries, max_labels)
  nearest <- vapply(entries$entry_s, function(s) {
    which.min(abs(scores$s - s))
  }, integer(1L))
  entries$score <- scores$score[nearest]
  rownames(entries) <- NULL
  metric <- arm$selection$metric %||% "loss"
  graphics::plot(
    scores$s, scores$score,
    type = "b", pch = 16L,
    xlim = rev(range(scores$s, finite = TRUE)),
    xlab = expression("relative penalty " * s == lambda / lambda[max] *
                        " (strong to weak)"),
    # The score is the mean over scoring units of each unit's own value,
    # which is not the pooled value over rows.
    ylab = paste0(
      "out-of-fold ", metric, ", mean of per-",
      arm$scoring_unit %||% "unit", " values"
    ),
    main = "Predictive performance along the entrance path",
    ...
  )
  selected_s <- arm$selection$selected_s
  if (is.finite(selected_s)) {
    graphics::abline(v = selected_s, lty = 2L, col = "#0072B2")
  }
  # Candidates entering at the same point share one label; drawn one by one
  # they would overprint each other.
  labels <- data.frame(
    entry_s = numeric(0), score = numeric(0), label = character(0),
    stringsAsFactors = FALSE
  )
  if (nrow(entries) > 0L) {
    labels <- stats::aggregate(
      entries$variable,
      list(entry_s = entries$entry_s, score = entries$score),
      paste,
      collapse = ", "
    )
    names(labels)[[3L]] <- "label"
    labels <- labels[order(labels$entry_s, decreasing = TRUE), , drop = FALSE]
    rownames(labels) <- NULL
    graphics::points(labels$entry_s, labels$score, pch = 21L, bg = "white")
    graphics::text(
      labels$entry_s, labels$score, labels = labels$label,
      pos = 3L, cex = label_cex, xpd = NA
    )
  }
  invisible(list(
    scores = as.data.frame(scores), entries = entries, labels = labels
  ))
}


#' @export
print.fb_select_certificate <- function(x, ...) {
  cat("<fb_select_certificate> ", x$contract, "\n", sep = "")
  cat("  axes       : workflow ", x$workflow, " / provider ", x$provider,
      " / prior ", x$prior, " / decision ", x$decision, "\n", sep = "")
  cat("  covariance : ", x$covariance$builder %||% "unknown builder", " / ",
      x$covariance$estimand %||% "not recorded", "\n", sep = "")
  cat("  prediction : ", x$resampling$prediction, "\n", sep = "")
  cat("  scored at  : ",
      if (is.na(x$resampling$scoring_unit %||% NA_character_)) {
        "no out-of-fold score on this route"
      } else {
        x$resampling$scoring_unit
      }, "\n",
      sep = "")
  if (!is.na(x$structure_verdict %||% NA_character_)) {
    cat("  structure  : ", x$structure_verdict, "\n", sep = "")
  }
  if (identical(x$decision_resolved_from, "auto")) {
    cat("  decision   : resolved from \"auto\" -- ",
        .fb_squish(x$decision_resolution_reason), "\n", sep = "")
  }
  cat("  seed       : ", x$resampling$seed %||% "none", "\n", sep = "")
  cat(
    "  calibration_n: ", x$calibration_n,
    if (identical(x$calibration_basis, "resampling_units")) {
      paste0(
        " (", paste(x$calibration_unit, collapse = ", "),
        "; the prior scales are calibrated against these rather than ",
        "against ", x$n_rows, " rows)"
      )
    } else {
      " (rows; no resampling unit was declared, so a row is a unit)"
    },
    "\n",
    sep = ""
  )
  cat("  fallback   : ", if (isTRUE(x$fallback_fired)) "FIRED" else "none",
      "\n", sep = "")
  if (length(x$prior_blocks) > 0L) {
    cat("  prior blocks:\n", sep = "")
    for (b in x$prior_blocks) {
      cat("    ", b$id, ". ", b$target, " ~ ", b$family,
          " (", b$hyperparameter_mode, ", route ", b$route, ")\n", sep = "")
    }
  }
  if (length(x$hyperparameter_treatment) > 0L) {
    cat("  hyperparameters: ",
        paste(x$hyperparameter_treatment, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$selection_data_reuse)) {
    cat("  data reuse : ", x$selection_data_reuse$level, " -- ",
        .fb_squish(x$selection_data_reuse$reason), "\n", sep = "")
  }
  if (!is.na(x$reference_engine %||% NA_character_)) {
    cat("  reference  : fitted by ", x$reference_engine, "\n", sep = "")
  }
  if (!is.null(x$reference)) {
    cat("  reference  : adequacy ", x$reference$verdict,
        ", ppc tails ", format(x$reference$ppc$p_mean, digits = 3L),
        " / ", format(x$reference$ppc$p_sd, digits = 3L), "\n", sep = "")
  }
  if (!is.null(x$projection)) {
    cat("  projection : ", x$projection$search,
        ", validate_search ", x$projection$validate_search,
        if (isTRUE(x$projection$optimistic)) " (OPTIMISTIC)" else "",
        if (!is.null(x$projection$decision_stage_search)) {
          paste0(", decision-stage search ", x$projection$decision_stage_search)
        } else {
          ""
        },
        "\n", sep = "")
  }
  if (length(x$benchmark) > 0L) {
    cat("  benchmark  :\n", sep = "")
    for (nm in names(x$benchmark)) {
      arm <- x$benchmark[[nm]]
      cat("    ", nm, ": ", arm$equivalence_level,
          if (!is.null(arm$tuning_scale)) {
            # Ridge has no `s`; its point is one lambda on every fold.
            point <- if (is.na(arm$selected_s)) {
              unique(arm$selected_lambda_canonical)
            } else {
              arm$selected_s
            }
            paste0(", tuned on ", arm$tuning_scale, " = ",
                   format(point, digits = 4L))
          } else {
            ""
          },
          "\n", sep = "")
    }
  }
  mo <- x$mode_oracle
  if (!is.null(mo)) {
    cat("  mode oracle: ",
        if (isTRUE(mo$available)) {
          paste0(paste(mo$arms, collapse = ", "), " [", mo$level, "]")
        } else {
          "none"
        },
        "\n", sep = "")
    cat("               ", .fb_squish(mo$reason), "\n", sep = "")
  }
  for (ps in list(x$reference_predictive, x$decision_predictive)) {
    if (is.null(ps)) {
      next
    }
    if (isTRUE(ps$available)) {
      cat("  ", ps$stage, " stage (out-of-fold ", ps$metric,
          ", at the resampling unit, contract ", ps$contract$hash, "):\n",
          sep = "")
      for (i in seq_len(nrow(ps$arms))) {
        cat("    ", format(ps$arms$arm[[i]], width = 8), " ",
            format(ps$arms$model[[i]], width = 18), " ",
            format(round(ps$arms$score[[i]], 4L), nsmall = 4L),
            "  (se ", format(round(ps$arms$se[[i]], 4L), nsmall = 4L),
            ", ", ps$arms$n_units[[i]], " ", ps$arms$unit[[i]], ")\n",
            sep = "")
      }
    } else {
      cat("  ", ps$stage, " stage: none -- ", .fb_squish(ps$reason), "\n",
          sep = "")
    }
  }
  if (!is.null(x$fold_preprocessing)) {
    cat("  fold prep  : constants and prior learned on the ",
        gsub("_", " ", x$fold_preprocessing$learned_on), " of each of ",
        nrow(x$fold_preprocessing$folds), " folds\n", sep = "")
  }
  cs <- x$design$candidate_space
  if (!is.null(cs)) {
    cat("  candidates : ", length(cs$declared), " declared columns; ",
        if (nrow(cs$dropped) == 0L) {
          paste0("each estimable on all ", cs$n_folds, " folds")
        } else {
          paste0(
            nrow(cs$dropped), " fold drops on ",
            length(unique(cs$dropped$fold)), " of ", cs$n_folds,
            " folds, coefficient NA: ",
            paste(unique(cs$dropped$column), collapse = ", ")
          )
        },
        "\n", sep = "")
  }
  cat("  approximations: ",
      if (length(x$approximations) == 0L) "none declared"
      else paste(vapply(x$approximations, .fb_squish, character(1L)),
                 collapse = "; "), "\n", sep = "")
  cat("  refusals   : ",
      if (length(x$refusals) == 0L) "none carried"
      else paste(x$refusals, collapse = "; "), "\n", sep = "")
  invisible(x)
}
