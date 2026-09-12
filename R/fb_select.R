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
# A base is admitted only when every axis passes, and the refusal names the
# measurement rather than the conclusion. `us()` and `fa()` on INLA stay out
# regardless of how a particular fit behaved: the measured genetic
# correlation is 0.142 against ASReml's 0.930, and multiplying the response
# by ten -- an exact transformation under which a correlation is invariant --
# does not reproduce it.
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
  # measurement that put them there.
  types <- vapply(
    base$fb$random_terms %||% list(),
    function(t) t$type %||% NA_character_, character(1L)
  )
  bad <- intersect(types, c("us", "fa", "unstructured", "factor_analytic"))
  if (length(bad) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_not_allowlisted",
      message = paste0(
        "fb_select(): the base carries ", paste(bad, collapse = ", "),
        ", which is not on the selection allowlist. INLA's us() returns a ",
        "genetic correlation of 0.142 where ASReml returns 0.930, and ",
        "rescaling the response by ten does not reproduce it, so a rotation ",
        "built on that covariance would be badly misspecified while looking ",
        "perfectly reproducible."
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
# The Bayesian route with a joint prior is brms and nothing else: INLA never
# approximates a scale-mixture prior, and lgm_gate()'s twelfth check refuses
# such a request there rather than fitting a nearby model. `projpred` joins
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

# .fb_select_resolve_penalty() --- fb_penalty() to a compiled brms prior.
#
# Four steps, in the order the contract fixes them, and each of them is a
# function a sibling unit owns: resolve the targets against the model matrix
# (`.fb_pen_resolve`), learn the scaling constants (`.fb_pen_scaling`),
# decide which of those constants the reference fit can actually apply, and
# compile (`.fb_pen_compile_brms`).
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
.fb_select_resolve_penalty <- function(base,
                                       candidates,
                                       penalty,
                                       protected_labels,
                                       calibration = NULL) {
  data <- base$data
  calibration <- calibration %||% .fb_select_calibration_n(NULL, data)
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
  compiled <- .fb_pen_compile_brms(
    blocks,
    coef_names_all,
    n_train = calibration$n
  )

  # The data the reference fit is emitted on, with the standardisable
  # candidate columns centred and scaled by the constants just stored.
  reference_data <- data
  for (b in blocks) {
    for (k in seq_along(b$coef_names)) {
      column <- candidate_names[b$columns[[k]]]
      if (isTRUE(scalable[[column]]) && !identical(b$scale[[k]], 1)) {
        reference_data[[column]] <-
          (data[[column]] - b$centre[[k]]) / b$scale[[k]]
      }
    }
  }

  scope_all_fixed <- any(vapply(
    blocks,
    function(b) identical(b$hyperparameters$scope %||% "candidates",
                          "all_fixed"),
    logical(1L)
  ))

  # The approximations this resolution commits to. They are settled here,
  # before anything is fitted, so `plan = TRUE` reports the same list the
  # certificate will carry rather than a shorter one.
  approximations <- character(0)
  if (length(left_alone) > 0L) {
    approximations <- c(approximations, paste0(
      "candidate_columns_not_standardised: ",
      paste(unique(left_alone), collapse = ", "),
      " are factor contrasts or interaction products with no data column ",
      "to centre and scale, so their prior scale is on their own scale"
    ))
  }
  if (!identical(penalty$scaling, "none")) {
    approximations <- c(
      approximations,
      paste0(
        "scaling_learned_on_all_rows: the reference model is one fit, so ",
        "the standardising constants are the whole sample's rather than ",
        "each training fold's"
      )
    )
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
  if (isTRUE(scope_all_fixed)) {
    approximations <- c(approximations, "protected_under_global_shrinkage")
  }

  list(
    map = map,
    partition = part,
    blocks = blocks,
    compiled = compiled,
    candidate_labels = candidate_labels,
    candidate_names = candidate_names,
    coef_names_all = coef_names_all,
    reference_data = reference_data,
    calibration = calibration,
    unscaled_columns = unique(left_alone),
    scope_all_fixed = scope_all_fixed,
    approximations = approximations,
    column_map = .fb_select_column_map(part, blocks, protected_labels, map)
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
                                   predictive_side_by_side = NULL,
                                   structure_gate = NULL,
                                   decision_resolution = NULL) {
  conv <- .fb_convergence_verdict(base)
  calibration <- .fb_select_calibration_normalise(calibration, base)
  predictive_side_by_side <- predictive_side_by_side %||%
    .fb_select_predictive_side_by_side(
      NULL, benchmark %||% list(),
      metric = .fb_select_benchmark_metric(cv$metric %||% "elpd")
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
        column_map = column_map
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
            provider_lambda_convention = arm$provider_lambda_convention
          )
        })
      },
      # Two fields where one overloaded label used to sit. `mode_oracle` is
      # the Tier-1 claim and is empty unless the covariance and the residual
      # scale were both held fixed, which is never true of a sampled route;
      # `predictive_side_by_side` is the out-of-fold comparison and is
      # written only when the arms really do share one loss, one unit and
      # one set of folds. Reading a glmnet support as a check on a posterior
      # is the category error the split exists to prevent.
      mode_oracle = .fb_select_mode_oracle(benchmark %||% list()),
      predictive_side_by_side = predictive_side_by_side,
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
# The machinery is `.fb_sel_bayes_cv()`'s, unchanged: the same folds, the same
# relabelling of held-out grouping levels, the same predictive draws and the
# same aggregation. Only the refit differs, and it differs in exactly one
# way -- the candidate block arrives under the compiled joint prior.
.fb_select_reference_cv <- function(base,
                                    candidates,
                                    resolved,
                                    protected_labels,
                                    resamples,
                                    prediction,
                                    metric,
                                    engine,
                                    seed,
                                    control = NULL) {
  refit <- function(fold_data, fixed, fit_seed) {
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
    args <- list(
      base = base,
      candidates = stats::as.formula(
        paste0("~ ", paste(extra, collapse = " + "))
      ),
      prior = resolved$compiled$prior,
      stanvars = resolved$compiled$stanvars,
      column_map = resolved$column_map,
      engine = engine,
      seed = fit_seed,
      data = fold_data
    )
    if (!is.null(control)) {
      args$control <- control
    }
    do.call(.fb_sel_reference_fit, args)
  }

  .fb_sel_bayes_cv(
    base = base,
    candidates = candidates,
    resamples = resamples,
    prediction = prediction,
    metric = metric,
    engine = engine,
    seed = seed,
    # The reference model's data, so the fold refits see the same
    # standardised candidate columns the reference itself was fitted on.
    data = resolved$reference_data,
    refit_fn = refit,
    # The provider is brms whatever fitted the base: a joint shrinkage prior
    # is always emitted there.
    provider = "brms",
    # The base's own recorded wall time, so the cost estimate costs no fit.
    probe_seconds = .fb_select_base_seconds(base)
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


# .fb_select_predictive_side_by_side() --- the arms' out-of-fold loss, on one
# metric, at one unit, on one set of folds, or a named reason there is none.
#
# Every score here is re-aggregated from the arms' own per-cluster tables
# through the same function, so "the same loss" is a property of the code
# rather than of a sentence. The metric is the benchmark metric rather than
# the decision's: the closed-form arms do not compute `elpd`, and a stanza
# that placed an elpd beside an nlpd would be the defect this replaces in a
# new spelling.
#
# The fold labels are compared, not assumed. Two arms built from one
# `fb_resamples()` at one seed do get the same partition, but the stanza's
# whole claim is "on one set of folds", and a claim that rests on an
# argument rather than on a comparison is the class of defect the
# certificate exists to close.
.fb_select_predictive_side_by_side <- function(primary_clusters,
                                               arms,
                                               metric,
                                               primary_folds = NULL,
                                               benchmark_folds = NULL) {
  absent <- function(reason) {
    list(available = FALSE, arms = NULL, metric = metric, reason = reason)
  }
  if (!is.data.frame(primary_clusters) || nrow(primary_clusters) == 0L) {
    return(absent(paste0(
      "the primary arm carries no out-of-fold score table, so there is ",
      "nothing to place beside a benchmark. `benchmark = TRUE` scores the ",
      "reference model on the benchmark's own folds and is what makes this ",
      "comparison available."
    )))
  }
  if (length(arms) == 0L) {
    return(absent(
      "no benchmark arm ran, so the primary arm's table stands alone."
    ))
  }
  if (!is.null(primary_folds) && !is.null(benchmark_folds) &&
        !setequal(primary_folds, benchmark_folds)) {
    return(absent(paste0(
      "the arms were scored on different partitions -- the primary on ",
      length(unique(primary_folds)), " fold(s), the benchmark on ",
      length(unique(benchmark_folds)),
      " -- so their out-of-fold losses are not comparable."
    )))
  }

  score_of <- function(clusters, model) {
    subset <- clusters[clusters$model == model, , drop = FALSE]
    if (nrow(subset) == 0L) {
      return(NULL)
    }
    summary <- tryCatch(
      .fb_sel_aggregate_bayes_scores(subset, metric = metric, unit = "cluster"),
      error = function(e) NULL
    )
    if (is.null(summary)) {
      return(NULL)
    }
    data.frame(
      model = model,
      score = as.numeric(summary$score),
      se = as.numeric(summary$se),
      n_units = as.integer(summary$n_units),
      unit = "cluster",
      stringsAsFactors = FALSE
    )
  }

  primary_model <- if ("full" %in% primary_clusters$model) {
    "full"
  } else {
    utils::tail(unique(primary_clusters$model), 1L)
  }
  primary_row <- score_of(primary_clusters, primary_model)
  if (is.null(primary_row)) {
    return(absent(paste0(
      "the primary arm's per-unit table carries no `", metric, "` column, ",
      "so it cannot be scored on the loss the benchmark arms scored."
    )))
  }
  rows <- list(cbind(arm = "bayes", primary_row, stringsAsFactors = FALSE))

  for (nm in names(arms)) {
    arm <- arms[[nm]]
    agg <- arm$scores
    if (is.data.frame(agg) && "full" %in% agg$model) {
      row <- agg[agg$model == "full", , drop = FALSE]
      rows[[length(rows) + 1L]] <- data.frame(
        arm = nm, model = "full",
        score = as.numeric(row$score[[1L]]),
        se = as.numeric(row$se[[1L]]),
        n_units = as.integer(row$n_units[[1L]]),
        unit = as.character(row$se_unit[[1L]]),
        stringsAsFactors = FALSE
      )
      next
    }
    chosen <- arm$selection
    if (is.list(chosen) && is.data.frame(chosen$scores) &&
          !is.null(chosen$selected_model)) {
      row <- chosen$scores[chosen$scores$model == chosen$selected_model, ,
                           drop = FALSE]
      if (nrow(row) == 1L) {
        rows[[length(rows) + 1L]] <- data.frame(
          arm = nm, model = as.character(chosen$selected_model),
          score = as.numeric(row$score[[1L]]),
          se = as.numeric(row$se[[1L]]),
          n_units = as.integer(row$n_units[[1L]]),
          unit = as.character(chosen$unit %||% NA_character_),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) < 2L) {
    return(absent(paste0(
      "no benchmark arm carries an out-of-fold score at the resampling ",
      "unit, so the primary arm's own table stands alone."
    )))
  }
  table <- do.call(rbind, rows)
  # An arm that aggregated at a different number of units did not score the
  # same loss, whatever the column is called. It is dropped and named rather
  # than averaged in: a two-fold standard error beside a twelve-site one is
  # the defect this stanza replaces, in a new spelling.
  agrees <- table$n_units == table$n_units[[1L]]
  dropped <- table$arm[!agrees]
  table <- table[agrees, , drop = FALSE]
  if (nrow(table) < 2L) {
    return(absent(paste0(
      "the arms aggregated at different numbers of units (",
      paste(paste0(
        do.call(rbind, rows)$arm, ": ", do.call(rbind, rows)$n_units
      ), collapse = "; "),
      "), so no two of them scored one loss at one unit."
    )))
  }
  list(
    available = TRUE,
    excluded = dropped,
    metric = metric,
    arms = table,
    reason = paste0(
      "every arm's score is the out-of-fold `", metric, "` aggregated at ",
      "the resampling unit on one set of folds, under one prediction ",
      "target. Similar predictive risk is what this licenses; identical ",
      "coefficients or supports are not.",
      if (length(dropped) > 0L) {
        paste0(
          " Excluded for aggregating at a different number of units: ",
          paste(dropped, collapse = ", "), "."
        )
      } else {
        ""
      }
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
                                 random_vars,
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
    random = random_vars,
    resamples = resamples,
    seed = seed,
    prediction = prediction,
    metric = metric
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
    if (!.fb_sel_provider_installed("glmnet")) {
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
          equivalence_level = .fb_sel_equivalence_level(
            prior_family = prior_family,
            provider = "glmnet",
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
#'   \item{`provider`}{`"brms"`, `"inla"`, `"glmnet"` or
#'     `"gls_reml_plugin"`.}
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
#' `fallback = "none"` is the default. A Bayesian request whose provider is
#' absent refuses with `selection_provider_absent`, naming the package and
#' the install line. `fallback = "glmnet"` is a deliberate choice; when it
#' fires, `workflow` becomes `penalized_optimizer`, the object prints the
#' substitution on its first line, and the certificate records
#' `fallback_fired = TRUE`. There is no configuration that produces a
#' penalised-likelihood answer under a Bayesian label.
#'
#' @param base A fitted flexyBayes model. Its fixed effects are the protected
#'   block and are never penalised.
#' @param candidates A one-sided formula over columns of the base fit's data.
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
#'   arm and, where a penalty was given, the `glmnet` arm on the same folds
#'   -- and scores the reference model itself out of fold, so the arms have
#'   one loss at one unit to be placed beside. That last part is what the
#'   certificate's `predictive_side_by_side` stanza reports and what the
#'   equivalence level `posterior_predictive_comparable` requires; without
#'   it the certificate says so by name rather than claiming a comparison
#'   nobody made. It costs one reference refit per fold per model.
#' @param resamples An [fb_resamples()] specification.
#' @param prediction `"new_cluster"` (default) or `"known_cluster"`. A
#'   held-out cluster's own effect is not estimable from the training rows,
#'   so scoring it conditionally would reward the model for a random effect
#'   it was handed.
#' @param decision An [fb_decision()]: the method, metric, rule and search.
#' @param engine An [fb_engine()] setting the sampler's chains, iterations
#'   and control for the reference fit.
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
  random_vars <- unique(unlist(lapply(
    base$fb$random_terms %||% list(),
    function(t) if (identical(t$type, "simple")) t$var else c(t$outer, t$inner)
  ))) %||% character(0)

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
      random_vars = random_vars,
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
      random_vars = random_vars,
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
        random_vars = random_vars,
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
    base, candidates, penalty, protected_labels, calibration = calibration
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
  # lgm_gate()'s twelfth check therefore belongs to the engine that WILL fit
  # the reference, not to the engine that fitted the base. It used to run on
  # the base's IR whenever the base was INLA, so an INLA user asking for a
  # Laplace prior was told the prior was impossible on INLA -- on a fit that
  # was never going to INLA. The check still fires where it means something:
  # a genuine INLA fit request carrying a joint block on the IR (fb_inla(),
  # fb_plan() routing at INLA) reaches it unchanged.
  reference_engine <- "brms"
  approximations <- resolved$approximations
  if (!identical(engine_name, reference_engine)) {
    approximations <- c(approximations, paste0(
      "reference_fitted_on_brms: the base was fitted by ", engine_name,
      " and the reference model carrying the joint shrinkage prior is ",
      "re-emitted on brms with the base's structure and its own priors ",
      "carried through, because a scale-mixture prior leaves the latent ",
      "field non-Gaussian given the hyperparameters and the nested Laplace ",
      "approximation is built on that assumption"
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
  adapt_delta <- if (inherits(engine, "fb_engine")) {
    engine$opts$adapt_delta
  } else {
    NULL
  }
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
      resolved = resolved,
      protected_labels = protected_labels,
      resamples = resamples,
      prediction = prediction,
      metric = metric,
      engine = engine,
      seed = seed,
      control = reference_args$control
    )
    bench <- .fb_select_benchmark(
      base = base,
      candidates = candidates,
      protected = protected,
      random_vars = random_vars,
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

  .fb_select_assemble(
    base = base,
    axes = list(
      workflow = "bayes",
      provider = "brms",
      prior = prior_family,
      decision = decision$method
    ),
    cv = list(
      folds = folds,
      # NULL unless `benchmark = TRUE` asked for the reference model to be
      # scored out of fold; the certificate's `predictive_side_by_side`
      # stanza names that as the reason where it is absent.
      aggregated = reference_cv$aggregated,
      clusters = reference_cv$clusters,
      rows = reference_cv$rows,
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
                                 random_vars,
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
    base, candidates, penalty, protected_labels, calibration = calibration
  )
  prior_family <- paste(
    unique(vapply(resolved$blocks, function(b) b$family, character(1L))),
    collapse = " + "
  )
  axes <- list(
    workflow = "penalized_optimizer",
    provider = "glmnet",
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
  .fb_sel_require_glmnet()

  bench <- .fb_select_benchmark(
    base = base,
    candidates = candidates,
    protected = protected,
    random_vars = random_vars,
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
        "fb_select(): the penalised-optimiser workflow produced no `glmnet` ",
        "arm. ",
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
      settings = list(method = "optimizer_support")
    ),
    benchmark = bench,
    metric = metric,
    approximations = resolved$approximations,
    refusals = bench$refusals,
    calibration = calibration,
    provider_versions = .fb_select_provider_versions("glmnet"),
    fallback_fired = fallback_fired,
    requested_workflow = requested_workflow,
    structure_gate = structure_gate,
    decision_resolution = decision_resolution,
    dots = list()
  )
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

  fold_labels <- function(folds) {
    if (length(folds %||% list()) == 0L) {
      return(NULL)
    }
    vapply(folds, function(sp) as.character(sp$label), character(1L))
  }
  side_by_side <- .fb_select_predictive_side_by_side(
    primary_clusters = cv$clusters,
    arms = benchmark$arms %||% list(),
    metric = .fb_select_benchmark_metric(cv$metric %||% metric),
    primary_folds = fold_labels(cv$folds),
    benchmark_folds = fold_labels(benchmark$folds)
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
    predictive_side_by_side = side_by_side,
    structure_gate = structure_gate,
    decision_resolution = decision_resolution
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
      lifecycle::deprecate_warn(
        when = "0.11.0",
        what = paste0("fb_select(", nm, ")"),
        with = if (identical(nm, "metric")) {
          "fb_decision(metric)"
        } else {
          NULL
        },
        details = if (identical(nm, "covariance")) {
          paste0(
            "The covariance summary is a property of the benchmark ",
            "workflow and is read from the fold object that used it, so a ",
            "declared value could disagree with what ran."
          )
        } else {
          NULL
        }
      )
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
    ps <- cert$predictive_side_by_side
    if (isTRUE(ps$available)) {
      cat("    side by side, out-of-fold ", ps$metric,
          " at the resampling unit:\n", sep = "")
      for (i in seq_len(nrow(ps$arms))) {
        cat("      ", format(ps$arms$arm[[i]], width = 8), " ",
            format(round(ps$arms$score[[i]], digits), nsmall = 2L),
            "  (", ps$arms$n_units[[i]], " ", ps$arms$unit[[i]], ")\n",
            sep = "")
      }
    } else {
      cat("    side by side: none -- ", .fb_squish(ps$reason), "\n", sep = "")
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
#' @param x An `fb_select` object.
#' @param ... Ignored.
#' @returns A `data.frame` with one row per design column: the column, the
#'   formula term it came from, its role, its block and family, whether the
#'   decision retained it, and, where a reference posterior was fitted, its
#'   posterior mean and standard deviation on the scientific scale.
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
  out
}


#' @export
glance.fb_select <- function(x, ...) {
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
    stringsAsFactors = FALSE
  )
}


#' Out-of-fold predictions from a selection
#'
#' @param x An `fb_select` object.
#' @param data Ignored; the out-of-fold predictions are the ones the folds
#'   computed and cannot be recomputed on other rows.
#' @param ... Ignored.
#' @returns A `data.frame` of the per-unit out-of-fold scores, or a zero-row
#'   frame with the same columns where no out-of-fold path ran. A projection
#'   assesses inside `projpred` and writes no table here on its own; under
#'   `benchmark = TRUE` the reference model is additionally scored out of
#'   fold and that table is what this returns.
#' @export
augment.fb_select <- function(x, data = NULL, ...) {
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
#' @param type `"elpd_by_size"`, `"rank_stability"`, `"shrinkage"` or
#'   `"benchmark"`.
#' @param ... Passed to the underlying plotting call.
#' @returns Invisibly, the data the display was drawn from. Called for the
#'   plot it draws.
#' @export
plot.fb_select <- function(x,
                           type = c("elpd_by_size", "rank_stability",
                                    "shrinkage", "benchmark"),
                           ...) {
  type <- match.arg(type)
  switch(
    type,
    elpd_by_size = .fb_select_plot_elpd(x, ...),
    rank_stability = .fb_select_plot_stability(x, ...),
    shrinkage = .fb_select_plot_shrinkage(x, ...),
    benchmark = .fb_select_plot_benchmark(x, ...)
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
        "\n", sep = "")
  }
  if (length(x$benchmark) > 0L) {
    cat("  benchmark  :\n", sep = "")
    for (nm in names(x$benchmark)) {
      cat("    ", nm, ": ", x$benchmark[[nm]]$equivalence_level, "\n",
          sep = "")
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
  ps <- x$predictive_side_by_side
  if (!is.null(ps)) {
    if (isTRUE(ps$available)) {
      cat("  side by side (out-of-fold ", ps$metric,
          ", at the resampling unit):\n", sep = "")
      for (i in seq_len(nrow(ps$arms))) {
        cat("    ", format(ps$arms$arm[[i]], width = 8), " ",
            format(round(ps$arms$score[[i]], 4L), nsmall = 4L),
            "  (se ", format(round(ps$arms$se[[i]], 4L), nsmall = 4L),
            ", ", ps$arms$n_units[[i]], " ", ps$arms$unit[[i]], ")\n",
            sep = "")
      }
    } else {
      cat("  side by side: none -- ", .fb_squish(ps$reason), "\n", sep = "")
    }
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
