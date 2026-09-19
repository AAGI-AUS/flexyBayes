# The reference fit --- the posterior a penalised selection decides from
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md section 5.
#
# The reference model is the base model's random and residual structure,
# its protected fixed effects under their own priors, and the candidate
# block under a compiled joint shrinkage prior. It is a RE-EMISSION of the
# base with an extended fixed formula, not a second fitting path: the same
# IR builder, the same brms emit, the same prior translation. Two things
# reach the emit that the public entry cannot carry -- a compiled prior
# data frame and a `stanvars` object -- and they arrive through the
# emitter's own arguments rather than through a parallel call to brm().
#
# Three gates stand between the fit and any decision, and all three are
# here because a selection read off an unfit reference model is a
# statement about the sampler or about the misfit, not about the
# candidates:
#
#   convergence   the tri-state contract, extended to the NUTS evidence
#   ppc           the observed mean and SD against replicated data
#   kappa         the horseshoe's shrinkage factors, reported against p0
#
# The third is a diagnostic, never a gate: it is exact only for a
# standardised near-orthogonal design, and it says so in its own label.


# --- reading the base ------------------------------------------------- #

# The analysis data a fit was built on. The INLA emit keeps it at `$data`;
# the brms emit does not carry one at all, and the model frame brms itself
# stored is then the only copy of what the fit actually saw.
#
# That copy is not the analysis data. brms keeps the model frame -- the
# variables the base model USED -- so a candidate column the base did not
# mention is absent from it. A reference fit extending a brms base is given
# the analysis data explicitly; the model frame is the fallback, and it is
# checked against the candidates before anything is fitted.
.fb_sel_fit_data <- function(fit) {
  d <- fit$data
  if (is.data.frame(d)) {
    return(d)
  }
  d <- tryCatch(fit$brms$data, error = function(e) NULL)
  if (is.data.frame(d)) {
    return(d)
  }
  NULL
}

# Term labels on the right-hand side of a one-sided candidate formula.
.fb_sel_candidate_labels <- function(candidates) {
  if (!inherits(candidates, "formula")) {
    stop(
      ".fb_sel_reference_fit(): `candidates` must be a one-sided formula.",
      call. = FALSE
    )
  }
  labels <- attr(stats::terms(candidates), "term.labels")
  labels <- labels[nzchar(labels)]
  if (length(labels) == 0L) {
    stop(
      ".fb_sel_reference_fit(): `candidates` names no terms.",
      call. = FALSE
    )
  }
  labels
}

# The protected fixed terms, as written in the base model's formula.
.fb_sel_protected_labels <- function(ir) {
  labels <- vapply(
    ir$fixed_terms %||% list(),
    function(t) t$label %||% t$var %||% NA_character_,
    character(1L)
  )
  labels[!is.na(labels) & nzchar(labels)]
}

# The left-hand side of a grouped slope term, with or without the
# intercept the IR records separately.
.fb_sel_slope_lhs <- function(t) {
  if (isTRUE(t$with_intercept)) t$slope_var else paste0("0 + ", t$slope_var)
}

# The base model's random terms in brms notation --- the spelling projpred
# uses for a group-level term in a ranking and accepts in `search_terms`.
.fb_sel_random_labels <- function(ir) {
  out <- character(0)
  for (t in ir$random_terms %||% list()) {
    type <- t$type %||% ""
    if (identical(type, "simple")) {
      out <- c(out, paste0("(1 | ", t$var, ")"))
    } else if (identical(type, "simple_slope_cor")) {
      lhs <- .fb_sel_slope_lhs(t)
      out <- c(out, paste0("(", lhs, " | ", t$var, ")"))
    } else if (identical(type, "simple_slope_uncor")) {
      lhs <- .fb_sel_slope_lhs(t)
      out <- c(out, paste0("(", lhs, " || ", t$var, ")"))
    } else if (type %in% c("vm", "ped")) {
      out <- c(out, paste0("(1 | gr(", t$var, ", cov = ", t$mat %||% "K", "))"))
    } else {
      # Anything else is a structure the projection gate refuses by name
      # before a search is attempted; it is not given a spelling here.
      out <- c(out, paste0("(", type, " | ", t$var %||% t$outer %||% "?", ")"))
    }
  }
  out
}

# Extend a model formula's right-hand side with further terms.
.fb_sel_extend_fixed <- function(fixed, add) {
  if (!inherits(fixed, "formula") || length(fixed) != 3L) {
    stop(
      ".fb_sel_reference_fit(): the base fit's model formula is not a ",
      "two-sided formula, so it cannot be extended.",
      call. = FALSE
    )
  }
  lhs <- paste(deparse(fixed[[2L]]), collapse = " ")
  rhs <- paste(deparse(fixed[[3L]]), collapse = " ")
  stats::as.formula(
    paste0(lhs, " ~ ", rhs, " + ", paste(add, collapse = " + ")),
    env = environment(fixed)
  )
}


# --- the fit ----------------------------------------------------------- #

# .fb_sel_reference_fit() --- re-emit the base with the candidate block.
#
# `prior` and `stanvars` are what the joint-prior compiler produced for the
# candidate block; they are passed through, never built here. `column_map`
# travels on the result as an attribute because posterior extraction,
# projection and back-transformation all read it and none of them should
# have to rebuild it.
#
# The sampler defaults are the floor the spec sets, not a recommendation:
# four chains, 2000 iterations (1500 kept), `adapt_delta = 0.95`. A
# horseshoe at barley shape produced eight divergences and a tail ESS of 86
# at exactly these settings, which is why the convergence contract, and not
# the defaults, is what decides whether the fit may be used.
.fb_sel_reference_fit <- function(base,
                                  candidates,
                                  prior = NULL,
                                  stanvars = NULL,
                                  column_map = NULL,
                                  engine = NULL,
                                  seed = NULL,
                                  data = NULL,
                                  ...) {
  if (!inherits(base, "flexybayes")) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_not_allowlisted",
      message = paste0(
        "fb_select(): `base` must be a fitted flexyBayes model for the ",
        "reference fit to extend."
      )
    ))
  }
  # A seed is required before anything is fitted. The sampler, the
  # posterior predictive check and the projection search all consume RNG,
  # and a default seed is an unreproducible run that looks reproducible.
  if (
    is.null(seed) ||
      length(seed) != 1L ||
      !is.numeric(seed) ||
      !is.finite(seed)
  ) {
    stop(.fb_refusal_condition(
      reason_code = "selection_seed_required",
      message = paste0(
        "fb_select(): the reference fit samples, and no `seed` was given. ",
        "No seed is invented here: a default seed is an unreproducible run ",
        "that looks reproducible."
      )
    ))
  }

  ir <- .fb_fit_ir(base)
  # The reference model is re-emitted from the IR, which a known variance per
  # row would not reach (R/known_variance.R).
  if (!is.null(.fb_known_variance_term(ir))) {
    .fb_sel_known_variance_refuse_route(
      "the projection's reference fit, whose posterior draws are projected,"
    )
  }
  data <- data %||% .fb_sel_fit_data(base)
  call_info <- base$extras$call_info
  if (is.null(ir) || is.null(data) || is.null(call_info$fixed)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_shape_unreadable",
      message = paste0(
        "fb_select(): this base fit does not carry the parsed model, the ",
        "data and the model formula the reference fit re-emits from."
      )
    ))
  }
  absent <- setdiff(all.vars(candidates), names(data))
  if (length(absent) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_shape_unreadable",
      message = paste0(
        "fb_select(): the candidate column(s) ",
        paste(absent, collapse = ", "),
        " are not in the data this base fit carries. A brms fit keeps the ",
        "model frame of the variables it used, not the analysis data, so a ",
        "candidate the base did not mention has to arrive with the ",
        "analysis data itself."
      )
    ))
  }

  candidate_terms <- .fb_sel_candidate_labels(candidates)
  protected_terms <- .fb_sel_protected_labels(ir)
  extended <- .fb_sel_extend_fixed(call_info$fixed, candidate_terms)

  fb <- .build_ir_polymorphic(
    fixed = extended,
    random = call_info$random,
    residual = call_info$residual,
    data = data,
    family = call_info$family %||% ir$family,
    link = call_info$link %||% ir$link,
    weights = call_info$weights,
    known_matrices = call_info$known_matrices %||% list(),
    prior = NULL,
    prior_fixed_sd = call_info$prior_fixed_sd %||% 100,
    prior_vc_sd = call_info$prior_vc_sd %||% 1,
    syntax = "auto"
  )
  # The base's own resolved prior, carried onto the extended IR rather than
  # re-derived: the reference model differs from the base in its fixed part
  # and nowhere else, and re-deriving would silently substitute this
  # release's default for whatever the base was actually fitted under.
  if (!is.null(ir$priors)) {
    fb$priors <- ir$priors
  }
  if (!is.null(ir$prior_scalars)) {
    fb$prior_scalars <- ir$prior_scalars
  }

  opts <- .fb_engine_opts(engine) %||% list()
  if (inherits(engine, "fb_engine") && !identical(engine$name, "brms")) {
    stop(
      ".fb_sel_reference_fit(): the reference fit is emitted through brms; ",
      "engine '", engine$name, "' cannot carry a joint shrinkage prior.",
      call. = FALSE
    )
  }
  dots <- list(...)
  control <- dots$control %||% list(adapt_delta = 0.95)
  dots$control <- NULL

  fit <- do.call(
    emit_brms,
    c(
      list(
        fb = fb,
        data = data,
        known_matrices = call_info$known_matrices %||% list(),
        weights = call_info$weights,
        n_samples = opts$n_samples %||% 1500L,
        warmup = opts$warmup %||% 500L,
        chains = opts$chains %||% 4L,
        seed = as.integer(seed),
        control = control,
        prior_fixed_sd = call_info$prior_fixed_sd %||% 100,
        prior_vc_sd = call_info$prior_vc_sd %||% 1,
        verbose = FALSE,
        mcmc_verbose = FALSE,
        extra_prior = prior,
        stanvars = stanvars,
        fixed = extended,
        random = call_info$random,
        residual = call_info$residual,
        family = call_info$family %||% ir$family,
        link = call_info$link %||% ir$link,
        data_name = call_info$data_name %||% NA_character_,
        requested_backend = "brms"
      ),
      dots
    )
  )

  # The IR and the analysis data, on the slots the selection substrate
  # reads. The shared dispatch path writes `$extras$fb_terms` and the INLA
  # emit writes `$fb` / `$data`; a reference fit built here never passes
  # through either, so it writes all three itself rather than being a fit
  # that half the module cannot read.
  fit$fb <- fb
  fit$data <- data
  fit$extras$fb_terms <- fb

  attr(fit, "column_map") <- column_map
  attr(fit, "candidate_terms") <- candidate_terms
  attr(fit, "protected_terms") <- protected_terms
  attr(fit, "random_terms") <- .fb_sel_random_labels(fb)
  attr(fit, "reference_seed") <- as.integer(seed)
  fit
}


# --- adequacy ---------------------------------------------------------- #

# The two posterior predictive tail probabilities of section 5.2.
#
# T(y) is the observed mean and the observed standard deviation; the tail
# probability is the posterior predictive probability that the replicated
# statistic falls at or below the observed one. A value near 0 or near 1 is
# a statistic the model cannot reproduce. The response is read from the
# Stan data the fit was actually built with, not from a column looked up by
# name in a data frame that may since have been transformed.
.fb_sel_ppc_tails <- function(fit, ndraws = 500L, seed = NULL) {
  out <- list(p_mean = NA_real_, p_sd = NA_real_, ndraws = NA_integer_)
  brmsfit <- fit$brms
  if (is.null(brmsfit) || !requireNamespace("brms", quietly = TRUE)) {
    return(out)
  }
  y <- tryCatch(as.numeric(brms::standata(brmsfit)$Y), error = function(e) NULL)
  if (is.null(y) || length(y) < 2L || anyNA(y)) {
    return(out)
  }
  draw <- function() {
    brms::posterior_predict(brmsfit, ndraws = as.integer(ndraws))
  }
  yrep <- tryCatch(
    if (is.null(seed)) draw() else .fb_sel_with_seed(seed, draw()),
    error = function(e) NULL
  )
  if (!is.matrix(yrep) || nrow(yrep) < 2L) {
    return(out)
  }
  rep_mean <- rowMeans(yrep)
  rep_sd <- apply(yrep, 1L, stats::sd)
  out$p_mean <- mean(rep_mean <= mean(y))
  out$p_sd <- mean(rep_sd <= stats::sd(y))
  out$ndraws <- nrow(yrep)
  out
}

.FB_SEL_PPC_TAIL_MIN <- 0.005
.FB_SEL_PPC_TAIL_MAX <- 0.995

# .fb_sel_reference_adequacy() --- may this posterior carry a decision?
#
# Tri-state, and the three states are the convergence contract's own: a
# `fail` on either gate is `fail`; evidence that is absent rather than bad
# is `unverified` and is carried onto the certificate rather than being
# read as a pass. Treedepth saturation is reported and does not decide:
# it is evidence the sampler was slow, not that the posterior is wrong.
.fb_sel_reference_adequacy <- function(fit,
                                       metric = c("elpd", "rmse", "nlpd"),
                                       ndraws = 500L,
                                       seed = NULL) {
  metric <- match.arg(metric)
  conv <- .fb_convergence_verdict(fit)
  ppc <- .fb_sel_ppc_tails(fit, ndraws = ndraws, seed = seed)

  reasons <- character(0)
  convergence_ok <- conv$ok
  if (identical(convergence_ok, FALSE)) {
    reasons <- c(reasons, conv$reasons)
  } else if (is.na(convergence_ok)) {
    reasons <- c(reasons, conv$reasons)
  }

  saturated <- conv$evidence$n_max_treedepth
  if (!is.null(saturated) && is.finite(saturated) && saturated > 0L) {
    reasons <- c(reasons, sprintf(
      paste0(
        "%d transition(s) reached the treedepth ceiling of %.0f; ",
        "reported, not decisive"
      ),
      as.integer(saturated),
      conv$evidence$max_treedepth
    ))
  }

  in_range <- function(p) {
    is.finite(p) && p > .FB_SEL_PPC_TAIL_MIN && p < .FB_SEL_PPC_TAIL_MAX
  }
  ppc_ok <- NA
  if (is.finite(ppc$p_mean) && is.finite(ppc$p_sd)) {
    ppc_ok <- in_range(ppc$p_mean) && in_range(ppc$p_sd)
    if (!ppc_ok) {
      if (!in_range(ppc$p_mean)) {
        reasons <- c(reasons, sprintf(
          paste0(
            "the observed mean sits at posterior predictive tail ",
            "probability %.4f"
          ),
          ppc$p_mean
        ))
      }
      if (!in_range(ppc$p_sd)) {
        reasons <- c(reasons, sprintf(
          paste0(
            "the observed standard deviation sits at posterior predictive ",
            "tail probability %.4f"
          ),
          ppc$p_sd
        ))
      }
    }
  } else {
    reasons <- c(reasons, paste0(
      "no replicated data: the posterior predictive check did not run"
    ))
  }

  verdict <- if (identical(convergence_ok, FALSE) || identical(ppc_ok, FALSE)) {
    "fail"
  } else if (is.na(convergence_ok) || is.na(ppc_ok)) {
    "unverified"
  } else {
    "pass"
  }

  list(
    verdict = verdict,
    reasons = as.character(reasons),
    ppc = list(p_mean = ppc$p_mean, p_sd = ppc$p_sd),
    convergence = convergence_ok,
    convergence_reasons = as.character(conv$reasons),
    ppc_ok = ppc_ok,
    evidence = conv$evidence,
    metric = metric
  )
}

# .fb_sel_require_adequate() --- the same measurement, as a gate.
#
# Two refusals, and which one fires says which gate failed: the sampler, or
# the model's description of the data. An `unverified` verdict passes and
# is returned, because absent evidence is not evidence of failure -- it is
# carried on the certificate and printed.
.fb_sel_require_adequate <- function(fit,
                                     metric = c("elpd", "rmse", "nlpd"),
                                     ndraws = 500L,
                                     seed = NULL) {
  adequacy <- .fb_sel_reference_adequacy(
    fit,
    metric = metric,
    ndraws = ndraws,
    seed = seed
  )
  if (identical(adequacy$convergence, FALSE)) {
    stop(.fb_refusal_condition(
      reason_code = "penalty_sampler_diagnostics_failed",
      message = paste0(
        "fb_select(): the reference model did not converge: ",
        paste(adequacy$convergence_reasons, collapse = "; "),
        ". A selection read off this posterior would be a statement about ",
        "the sampler."
      ),
      adequacy = adequacy
    ))
  }
  if (identical(adequacy$ppc_ok, FALSE)) {
    stop(.fb_refusal_condition(
      reason_code = "penalty_reference_model_inadequate",
      message = paste0(
        "fb_select(): the reference model does not describe the data ",
        "(posterior predictive tail probabilities: mean ",
        format(adequacy$ppc$p_mean, digits = 4L), ", SD ",
        format(adequacy$ppc$p_sd, digits = 4L),
        "; both must lie in (", .FB_SEL_PPC_TAIL_MIN, ", ",
        .FB_SEL_PPC_TAIL_MAX, ")). A projection inherits the reference ",
        "model's misfit."
      ),
      adequacy = adequacy
    ))
  }
  invisible(adequacy)
}


# --- shrinkage factors -------------------------------------------------- #

# .fb_sel_kappa() --- the horseshoe's shrinkage factors, and their caveat.
#
#   kappa_j = 1 / (1 + n * sigma^-2 * tau^2 * lambda_tilde_j^2)
#
# brms 2.23.0 does not save `hs_local` in the draws. It saves `sdb_<coef>`,
# and the generated Stan says what that is: `scales_horseshoe()` returns
# `lambda_tilde * tau`, and `sdb = scales[1:Kc]`. So `sdb_j^2` IS
# `tau^2 * lambda_tilde_j^2` and the formula is evaluated on it directly --
# read from the program that was compiled, not from the parameter names a
# paper uses.
#
# `n` is the number of rows the model was fitted on. The label is not
# decoration: the scalar formula is the shrinkage factor of a standardised
# ORTHOGONAL design, and a correlated candidate block does not have one
# shrinkage factor per coefficient at all.
.fb_sel_kappa <- function(fit, blocks) {
  label <- "orthogonal_design_kappa_approx"
  empty <- list(
    kappa = numeric(0),
    label = label,
    effective_nonzero_at_half = NA_integer_,
    p0 = NA_real_,
    n = NA_integer_,
    available = FALSE,
    reason = "no regularized_horseshoe block carries a saved scale parameter"
  )
  brmsfit <- fit$brms
  if (is.null(brmsfit) || !requireNamespace("posterior", quietly = TRUE)) {
    return(empty)
  }

  hs_blocks <- Filter(
    function(b) {
      identical(b$family %||% "", "regularized_horseshoe") ||
        identical(b$family %||% "", "horseshoe")
    },
    blocks %||% list()
  )
  if (length(hs_blocks) == 0L) {
    empty$reason <- "no regularized_horseshoe block in the compiled prior"
    return(empty)
  }
  coef_names <- unique(unlist(
    lapply(hs_blocks, function(b) b$coef_names),
    use.names = FALSE
  ))
  p0 <- suppressWarnings(as.numeric(
    hs_blocks[[1L]]$hyperparameters$p0 %||% NA_real_
  ))

  draws <- posterior::as_draws_matrix(brmsfit)
  available <- colnames(draws)
  wanted <- paste0("sdb_", coef_names)
  present <- wanted %in% available
  if (!any(present) || !("sigma" %in% available)) {
    return(empty)
  }
  n <- nrow(.fb_sel_fit_data(fit) %||% data.frame())
  if (!is.finite(n) || n == 0L) {
    n <- length(tryCatch(
      as.numeric(brms::standata(brmsfit)$Y),
      error = function(e) numeric(0)
    ))
  }

  sigma2 <- as.numeric(draws[, "sigma"])^2
  kappa <- vapply(
    which(present),
    function(j) {
      sdb2 <- as.numeric(draws[, wanted[[j]]])^2
      mean(1 / (1 + n * sdb2 / sigma2))
    },
    numeric(1L)
  )
  names(kappa) <- coef_names[present]

  list(
    kappa = kappa,
    label = label,
    effective_nonzero_at_half = sum(kappa < 0.5),
    p0 = p0,
    n = as.integer(n),
    available = TRUE,
    reason = if (all(present)) {
      NA_character_
    } else {
      paste0(
        "no saved scale parameter for: ",
        paste(coef_names[!present], collapse = ", ")
      )
    }
  )
}
