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
                                       protected_labels) {
  data <- base$data
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

  coef_names_all <- .fb_pen_brms_coef_names(
    colnames(map$X)[map$assign != 0L]
  )
  compiled <- .fb_pen_compile_brms(
    blocks,
    coef_names_all,
    n_train = nrow(data)
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
    unscaled_columns = unique(left_alone),
    scope_all_fixed = scope_all_fixed,
    approximations = approximations,
    column_map = .fb_select_column_map(part, blocks, protected_labels)
  )
}

# .fb_select_column_map() --- one row per design column, with its role, its
# block and the scaling constants that were actually applied to it.
.fb_select_column_map <- function(part, blocks, protected_labels) {
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
  data.frame(
    column = c(protected_names, candidate_names),
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
                                   calibration_n = NA_integer_,
                                   provider_versions = list(),
                                   column_map = NULL) {
  conv <- .fb_convergence_verdict(base)
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
        # The projection route scores at the resampling unit inside the
        # search itself and carries no aggregated table, so the unit is read
        # off the projection's own scoring record.
        scoring_unit = cv$aggregated$se_unit[[1L]] %||%
          projection$scoring$unit %||% NA_character_,
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
      provider_versions = provider_versions,
      # The number of rows a prior scale was calibrated against. It is on the
      # certificate rather than implied because in clustered data the row
      # count overstates the information, and a reader has to be able to see
      # by how much.
      calibration_n = calibration_n,
      fallback_fired = isTRUE(fallback_fired),
      # Named because a reader must be able to see what was not exact.
      approximations = as.character(approximations),
      refusals = as.character(refusals),
      contract = "SPEC_bayesian_selection_v5_2026-09-12"
    ),
    class = "fb_select_certificate"
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
                                 decision = NULL) {
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
      sigma_fixed = FALSE
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
          cluster = NULL
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
            sigma_fixed = FALSE
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
#' them licenses the sentence "the same as glmnet".
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
#'   arm and, where a penalty was given, the `glmnet` arm on the same folds.
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
      decision = decision
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
      calibration_n = nrow(base$data),
      provider_versions = list(),
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
      requested_workflow = requested_workflow
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
        structure_gate = NULL,
        resamples = resamples,
        prediction = prediction,
        fits_per_fold = 2L,
        candidates = candidates
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
        decision = decision
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
      calibration_n = nrow(base$data),
      provider_versions = .fb_select_provider_versions(cv$provider),
      dots = dots
    ))
  }

  # --- the Bayesian route, with a joint shrinkage prior -------------------- #
  resolved <- .fb_select_resolve_penalty(
    base, candidates, penalty, protected_labels
  )
  prior_family <- paste(
    unique(vapply(resolved$blocks, function(b) b$family, character(1L))),
    collapse = " + "
  )

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
      structure_gate = .fb_sel_structure_gate(base$fb),
      resamples = resamples,
      prediction = prediction,
      fits_per_fold = if (isTRUE(decision$validate_search)) 1L else 0L,
      candidates = candidates,
      compiled = resolved$compiled,
      approximations = resolved$approximations
    ))
  }

  # The joint prior is on the IR before the gate reads it, so a request that
  # routes at INLA refuses at the gate rather than at the emitter.
  fb_gate <- base$fb
  fb_gate$penalty_blocks <- resolved$blocks
  if (identical(engine_name, "inla")) {
    gated <- lgm_gate(fb_gate)
    if (is_lgm_refusal(gated)) {
      primary <- gated$failures[[1L]]
      stop(.fb_refusal_condition(
        reason_code = if (
          identical(primary$rule_id, "prior_family_unsupported_inla")
        ) {
          "prior_family_unsupported_inla"
        } else {
          "selection_base_not_allowlisted"
        },
        message = paste0("fb_select(): ", primary$reason)
      ))
    }
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
      prior_family = prior_family,
      blocks = resolved$blocks,
      decision = decision
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
      aggregated = NULL,
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
    approximations = resolved$approximations,
    calibration_n = nrow(base$data),
    provider_versions = .fb_select_provider_versions("brms", decision),
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
                                 requested_workflow) {
  resolved <- .fb_select_resolve_penalty(
    base, candidates, penalty, protected_labels
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
      structure_gate = NULL,
      resamples = resamples,
      prediction = prediction,
      fits_per_fold = 0L,
      candidates = candidates,
      compiled = resolved$compiled,
      approximations = resolved$approximations
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
    decision = decision
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
    calibration_n = nrow(base$data),
    provider_versions = .fb_select_provider_versions("glmnet"),
    fallback_fired = fallback_fired,
    requested_workflow = requested_workflow,
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
                                calibration_n = NA_integer_,
                                provider_versions = list(),
                                fallback_fired = FALSE,
                                requested_workflow = NULL,
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
    calibration_n = calibration_n,
    provider_versions = provider_versions,
    column_map = column_map
  )

  structure(
    list(
      workflow = axes$workflow,
      provider = axes$provider,
      prior = axes$prior,
      decision = axes$decision,
      requested_workflow = requested_workflow %||% axes$workflow,
      fallback_fired = isTRUE(fallback_fired),
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
    " / decision ", x$decision, "\n",
    sep = ""
  )
  cat(
    "  base      : ", cert$base$family, "/", cert$base$link,
    ", ", cert$base$n_observed, " rows, convergence ",
    if (is.na(cert$base$convergence)) "UNVERIFIED"
    else if (cert$base$convergence) "ok" else "FAILED",
    "\n",
    sep = ""
  )
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
      cert$resampling$prediction, ", scored at ",
      cert$resampling$scoring_unit %||% "the fold", "\n", sep = "")

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
#' @param x An `fb_select` object.
#' @param ... Ignored.
#' @returns A `data.frame` with one row per design column: the column, its
#'   role, its block and family, whether the decision retained it, and, where
#'   a reference posterior was fitted, its posterior mean and standard
#'   deviation on the scientific scale.
#' @export
tidy.fb_select <- function(x, ...) {
  map <- x$column_map
  if (is.null(map)) {
    map <- data.frame(
      column = character(0), role = character(0), block = integer(0),
      family = character(0), centre = numeric(0), scale = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  chosen <- x$selection$chosen_terms %||% character(0)
  n <- nrow(map)
  out <- data.frame(
    term = map$column,
    role = map$role,
    block = map$block,
    family = map$family,
    retained = map$column %in% chosen | map$role == "protected",
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
#'   frame with the same columns where no out-of-fold path ran (the
#'   projection assesses inside `projpred` and writes no per-row table here).
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
  cat("  scored at  : ", x$resampling$scoring_unit %||% "the fold", "\n",
      sep = "")
  cat("  seed       : ", x$resampling$seed %||% "none", "\n", sep = "")
  cat("  calibration_n: ", x$calibration_n, "\n", sep = "")
  cat("  fallback   : ", if (isTRUE(x$fallback_fired)) "FIRED" else "none",
      "\n", sep = "")
  if (length(x$prior_blocks) > 0L) {
    cat("  prior blocks:\n", sep = "")
    for (b in x$prior_blocks) {
      cat("    ", b$id, ". ", b$target, " ~ ", b$family,
          " (", b$hyperparameter_mode, ", route ", b$route, ")\n", sep = "")
    }
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
  cat("  approximations: ",
      if (length(x$approximations) == 0L) "none declared"
      else paste(vapply(x$approximations, .fb_squish, character(1L)),
                 collapse = "; "), "\n", sep = "")
  cat("  refusals   : ",
      if (length(x$refusals) == 0L) "none carried"
      else paste(x$refusals, collapse = "; "), "\n", sep = "")
  invisible(x)
}
