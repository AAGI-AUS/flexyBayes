# .dispatch_backend() -- shared dispatch helper.
#
# Shared so that both flexybayes() (asreml entry) and fb_brms()
# (brms entry) drive the same backend dispatch. Byte-equivalent
# under snapshot guard (test-emit-greta.md + test-backend-auto.R).
#
# Routes to emit_inla() or emit_greta() per backend = c("greta",
# "inla", "auto"), populates fit$extras$backend_decision with the
# uniform-shape trace, and attaches the fb_terms IR for downstream
# canonical_names() resolution.
#
# The trace shape carries four fields (`preflight_summary`,
# `representation_plan`, `rejected_routes`,
# `routing_policy_version`); the dispatch threads the
# `<fb_preflight>` result into `lgm_gate()` so the
# `memory_feasibility_inla` rule can fire on INLA-specific memory
# infeasibility; and the routing policy is centralised in the
# `.routing_policy_table` registry below.
#
# Internal -- not exported.

# ---------------------------------------------------------------- #
# Routing policy + trace helpers                                    #
# ---------------------------------------------------------------- #

# The routing policy version string. Bump when the policy table
# changes shape (rows added / chosen_backend semantics changed) so
# every fit object produced under a given policy carries the
# version that decided its dispatch.
.ROUTING_POLICY_VERSION <- "stage5a_v1"


# .routing_policy_table -- single source of truth for the routing
# policy. Each row encodes one (user_request, gate_outcome,
# preflight_outcome) tuple and maps to the chosen backend plus a
# canonical reason code from the fixed vocabulary. The table is
# internal-only (no @export); the dispatch logic queries it via
# `.resolve_routing()`.
#
# NA values in the key columns mean "any value" (a wildcard match
# for that input). For example, explicit user requests bypass the
# gate / preflight columns -- a "greta" / "brms" request matches
# regardless of gate / preflight outcomes.
#
# Reason-code vocabulary:
#   explicit_greta               -- user requested greta directly.
#   explicit_brms                -- user requested brms (Stan via brms).
#   explicit_inla_accept         -- user requested inla; gate accepted.
#   explicit_inla_gate_refused   -- user requested inla; gate refused
#                                   (raises rather than fallback).
#   auto_inla_accept             -- auto + gate accept + INLA available.
#   auto_inla_unavailable        -- auto + gate accept + INLA missing
#                                   from .libPaths().
#   auto_gate_refuse_structural  -- auto + gate refuse on a structural
#                                   reason (the structural gate rules
#                                   or the verification rule).
#   auto_gate_refuse_memory      -- auto + gate refuse on the
#                                   memory_feasibility_inla rule.
#   backend_not_activated        -- gretaR slot rejected because the
#                                   backend is a dormant opt-in (enable with
#                                   options(flexyBayes.gretaR_activated = TRUE);
#                                   see R/backend_registry.R for lifecycle).
#
# The table is consulted by `.resolve_routing()`; the actual
# backend invocation lives in `.dispatch_backend()` below. A
# divergence between the table-resolved chosen backend and the
# code-resolved chosen backend surfaces as a test failure in
# tests/testthat/test-backend-routing-trace.R subtest (n).

# ---------------------------------------------------------------- #
# Constructor-noun resolution                                       #
# ---------------------------------------------------------------- #
# Resolve a `backend` argument to its engine-name string. An
# `fb_engine()` object resolves to its `name`; a bare string (or the
# multi-value default vector) passes through unchanged for match.arg().
.resolve_engine_string <- function(backend) {
  if (inherits(backend, "fb_engine")) {
    return(backend$name)
  }
  backend
}

# The tuning options carried by a `backend = fb_engine(...)` argument,
# or NULL for the bare-string form. Consumed by the fitting verbs to
# override the corresponding sampler controls.
.fb_engine_opts <- function(backend) {
  if (inherits(backend, "fb_engine") && length(backend$opts)) {
    backend$opts
  } else {
    NULL
  }
}


.routing_policy_table <- function() {
  data.frame(
    user_request = c(
      "greta",
      "brms",
      "inla",
      "inla",
      "auto",
      "auto",
      "auto",
      "auto",
      "auto"
    ),
    gate_outcome = c(
      NA_character_,
      NA_character_,
      "accept",
      "refuse",
      "accept",
      "accept",
      "refuse_structural",
      "refuse_memory",
      "accept"
    ),
    preflight_outcome = c(
      NA_character_,
      NA_character_,
      NA_character_,
      NA_character_,
      "clear",
      "clear",
      NA_character_,
      NA_character_,
      NA_character_
    ),
    inla_installed = c(NA, NA, NA, NA, TRUE, FALSE, NA, NA, NA),
    chosen_backend = c(
      "greta",
      "brms",
      "inla",
      NA_character_, # explicit-inla refusal raises
      "inla",
      "greta",
      "greta",
      "greta",
      "greta"
    ), # auto + accept + skipped preflight
    reason_code = c(
      "explicit_greta",
      "explicit_brms",
      "explicit_inla_accept",
      "explicit_inla_gate_refused",
      "auto_inla_accept",
      "auto_inla_unavailable",
      "auto_gate_refuse_structural",
      "auto_gate_refuse_memory",
      "auto_inla_accept"
    ),
    stringsAsFactors = FALSE
  )
}


# .resolve_routing() -- given the three input keys and the
# environmental flags (INLA install + gretaR activation), returns
# the canonical chosen_backend + rejected_routes per the
# `.routing_policy_table`. Used to populate the trace at dispatch
# time; the actual backend invocation lives in `.dispatch_backend()`
# below.
#
# Returns list with:
#   $chosen_backend   character(1)
#   $reason_code      character(1)
#   $rejected_routes  list of list(backend, reason_code); empty for
#                     explicit user requests (policy bypassed).
.resolve_routing <- function(
  user_request,
  gate_outcome,
  preflight_outcome,
  inla_installed,
  gretaR_activated
) {
  is_explicit <- user_request %in% c("greta", "brms", "inla")
  table <- .routing_policy_table()

  # Find the matching row -- columns with NA are wildcards.
  match_row <- function(row) {
    (is.na(row$user_request) || row$user_request == user_request) &&
      (is.na(row$gate_outcome) ||
        row$gate_outcome == (gate_outcome %||% NA_character_)) &&
      (is.na(row$preflight_outcome) ||
        row$preflight_outcome == (preflight_outcome %||% NA_character_)) &&
      (is.na(row$inla_installed) ||
        identical(row$inla_installed, inla_installed))
  }
  hits <- vapply(
    seq_len(nrow(table)),
    function(i) match_row(table[i, ]),
    logical(1L)
  )
  if (!any(hits)) {
    # No policy row matches -- fall back to greta with a generic
    # reason code. This branch should never fire in production
    # (the table covers every conceivable tuple); if it does, the
    # test suite catches it via subtest (n).
    return(list(
      chosen_backend = "greta",
      reason_code = "policy_table_no_match_fallback_greta",
      rejected_routes = list()
    ))
  }
  row <- table[which(hits)[1L], ]

  # Rejected routes: empty for explicit user requests (policy
  # bypassed); otherwise enumerate the non-chosen backends with
  # per-backend reason codes derived from the same routing inputs.
  rejected <- list()
  if (!is_explicit) {
    # auto's policy considers {inla, greta, gretaR}; brms never
    # appears in auto's candidate set.
    chosen <- row$chosen_backend
    candidates <- c("inla", "greta", "gretaR")
    for (cand in setdiff(candidates, chosen)) {
      cand_reason <- switch(
        cand,
        "inla" = .inla_rejection_reason(gate_outcome, inla_installed),
        "greta" = "not_chosen_by_policy",
        "gretaR" = if (isTRUE(gretaR_activated)) {
          "not_chosen_by_policy"
        } else {
          "backend_not_activated"
        }
      )
      rejected[[length(rejected) + 1L]] <- list(
        backend = cand,
        reason = cand_reason
      )
    }
  }

  list(
    chosen_backend = row$chosen_backend,
    reason_code = row$reason_code,
    rejected_routes = rejected
  )
}


# Reason code for the INLA rejection in auto's rejected_routes
# list. Derived from the same inputs that drove the policy row
# match; surfaces the specific gate / install signal that
# excluded INLA from the chosen backend on the auto path.
.inla_rejection_reason <- function(gate_outcome, inla_installed) {
  if (identical(gate_outcome, "refuse_memory")) {
    return("memory_infeasibility_inla")
  }
  if (identical(gate_outcome, "refuse_structural")) {
    return("structural_infeasibility_inla")
  }
  if (isTRUE(inla_installed) || is.na(inla_installed)) {
    return("not_chosen_by_policy")
  }
  "backend_not_installed"
}


# .build_routing_decision() -- canonical 8-field trace constructor.
# All six dispatch sites in this file build their trace via this
# helper so the trace shape is uniform. The four base fields
# (backend, path, gate_checks, reason) are preserved verbatim; the
# four extended fields (preflight_summary, representation_plan,
# rejected_routes, routing_policy_version) are populated alongside
# them.
.build_routing_decision <- function(
  backend,
  path,
  gate_checks,
  reason,
  preflight_summary = NULL,
  representation_plan = NULL,
  rejected_routes = list(),
  routing_policy_version = .ROUTING_POLICY_VERSION
) {
  list(
    backend = backend,
    path = path,
    gate_checks = gate_checks,
    reason = reason,
    preflight_summary = preflight_summary,
    representation_plan = representation_plan,
    rejected_routes = rejected_routes,
    routing_policy_version = routing_policy_version
  )
}


# .representation_plan_from_preflight() -- slim per-term plan
# derived from the `<fb_preflight>` per_term_estimate slot. One
# entry per IR term with the representation class + a one-line
# justification. The full preflight per-term
# entries remain available on $preflight_summary; this slim
# representation surfaces the chosen-class summary for trace
# consumers.
.representation_plan_from_preflight <- function(preflight) {
  if (is.null(preflight) || is.null(preflight$per_term_estimate)) {
    return(NULL)
  }
  per_term <- preflight$per_term_estimate
  out <- lapply(names(per_term), function(label) {
    entry <- per_term[[label]]
    list(
      term_id = label,
      representation_class = entry$representation_class %||% NA_character_,
      justification = .representation_justification(entry)
    )
  })
  names(out) <- names(per_term)
  out
}

.representation_justification <- function(entry) {
  if (isTRUE(entry$unknown_representation)) {
    return("representation not characterised by preflight estimator")
  }
  bytes <- entry$design_memory_bytes %||% NA_real_
  if (is.na(bytes)) {
    return("byte estimate unavailable")
  }
  if (isTRUE(entry$aggregated_likelihood_candidate)) {
    sprintf("indexed (%.1f MB); aggregation-eligible", bytes / 1024^2)
  } else {
    sprintf("indexed (%.1f MB)", bytes / 1024^2)
  }
}


# .check_approximate_scheme() -- entry-function guard called from
# `flexybayes()` and `fb_brms()` BEFORE match.arg(backend) so
# approximate-scheme requests refuse with a structured pointer
# rather than match.arg's generic "should be one of" error.
#
# Pattern: any backend name containing "approximate" (case-
# insensitive) or matching ^variational_ is considered an
# approximate-scheme request. Catches the named examples
# (inla_pardiso_approximate, variational_advi) plus future
# `approximate_<scheme>` shapes.
.check_approximate_scheme <- function(backend) {
  if (!is.character(backend) || length(backend) == 0L) {
    return(invisible(NULL))
  }
  # User may pass `c("greta", "inla", "auto")` (the default arg
  # list) -- only the actual single-value selection matters.
  # Inspect every element defensively.
  for (b in backend) {
    if (is.na(b) || !nzchar(b)) {
      next
    }
    looks_approximate <- grepl(
      "approximate",
      b,
      ignore.case = TRUE,
      fixed = FALSE
    ) ||
      grepl("^variational_", b)
    if (!looks_approximate) {
      next
    }
    stop(.fb_refusal_condition(
      reason_code = "approximate_route_not_yet_registered",
      message = paste0(
        "backend = \"",
        b,
        "\": approximate-scheme dispatch is ",
        "not yet registered (reason_code = ",
        "approximate_route_not_yet_registered).\n",
        "The approximation registry does not yet register this ",
        "scheme; the routing layer refuses any backend ",
        "name carrying an exactness label of \"approximate_<scheme>\".\n",
        "The package's identity rests on exactness labelling; ",
        "silently routing an approximate ",
        "request would violate that contract.\n",
        "Re-route via backend = \"greta\" (full MCMC; \"exact\"), ",
        "backend = \"inla\" (Laplace; \"exact\"), or ",
        "backend = \"brms\" (Stan via brms; \"exact\")."
      ),
      family_class = "flexybayes_approximate_route_refusal",
      backend = b
    ))
  }
  invisible(NULL)
}


# ---------------------------------------------------------------- #
# Native-greta dispatch                                             #
# ---------------------------------------------------------------- #
# .dispatch_native_greta() -- route a greta-source IR (a native model
# graph wrapped by fb_from_greta()) to the direct greta::mcmc() fit. A
# native graph is greta-only by construction (greta built it), so:
#   - any backend other than "greta" / "auto" is a structured refusal
#     (the registered native_greta_requires_greta_backend code);
#   - the code-inspection (return_code / review_code) and planning
#     (plan) modes do not apply (there is no flexyBayes-generated code,
#     and no formula triple to plan over) -- they raise a clear,
#     forward-pointing error, preserving the pre-v0.5.0 fb_greta()
#     deferrals.
# On the accepted path it delegates to .fit_native_greta() (R/fb_greta.R).
.dispatch_native_greta <- function(
  fb,
  backend,
  n_samples,
  warmup,
  chains,
  verbose,
  mcmc_verbose,
  return_code,
  review_code,
  plan,
  the_call
) {
  if (!backend %in% c("greta", "auto")) {
    stop(.fb_refusal_condition(
      reason_code = "native_greta_requires_greta_backend",
      message = paste0(
        "A native greta model graph is fit by greta::mcmc() and is ",
        "greta-only by construction; backend = \"",
        backend,
        "\" cannot ",
        "fit it. Use backend = \"greta\" (or the default), or rebuild the ",
        "model in the ASReml / brms formula grammar to reach the ",
        backend,
        " engine."
      ),
      backend = backend
    ))
  }

  if (isTRUE(return_code)) {
    stop(
      "`return_code = TRUE` is not available for a native greta model: ",
      "the graph is already greta-side, so there is no ",
      "flexyBayes-generated code string to return.",
      call. = FALSE
    )
  }

  if (isTRUE(review_code)) {
    stop(
      "`review_code = TRUE` is not available for a native greta model: ",
      "the inspect-then-fit token reflects flexyBayes-generated code, ",
      "and the graph is user-authored. Fit it directly with fb() / ",
      "fb_greta().",
      call. = FALSE
    )
  }

  if (isTRUE(plan)) {
    stop(
      "`plan = TRUE` is not available for a native greta model: the ",
      "planning object summarises a parsed formula triple, which a ",
      "native graph does not carry.",
      call. = FALSE
    )
  }

  .fit_native_greta(
    fb = fb,
    n_samples = n_samples,
    warmup = warmup,
    chains = chains,
    verbose = verbose,
    mcmc_verbose = mcmc_verbose,
    the_call = the_call
  )
}


# ---------------------------------------------------------------- #

.dispatch_backend <- function(
  fb,
  data,
  backend,
  known_matrices,
  weights,
  n_samples,
  warmup,
  chains,
  prior_fixed_sd,
  prior_vc_sd,
  verbose,
  mcmc_verbose,
  return_code,
  the_call,
  fixed,
  random,
  rcov,
  family,
  link,
  data_name,
  aggregate = "auto"
) {
  aggregate <- .normalise_aggregate(aggregate)

  # Factor-dictionary persistence. Build a metadata-only
  # <fb_dataset> descriptor at dispatch time so every exit path can
  # attach it to fit$extras$fb_dataset for downstream
  # predict.flexybayes(newdata, ...) consumption. Stripping $data
  # avoids retaining the full training frame on the fit object;
  # dictionaries + col_types + n_rows are the only fields predict
  # needs.
  fb_dataset_meta_for_fit <- tryCatch(
    .fb_dataset_metadata(.fb_dataset(data)),
    error = function(e) NULL
  )

  # gretaR backend (activated): an explicit opt-in backend driven out of
  # process (greta + gretaR share symbols, cannot co-load). The capability
  # predicate refuses model classes gretaR cannot express (structured
  # covariance) before the worker is launched. See R/emit_gretaR.R +
  # inst/backend_contract.md.
  if (identical(backend, "gretaR")) {
    cap <- .capability_gretaR(fb)
    if (!isTRUE(cap)) {
      stop(.fb_refusal_condition(
        reason_code = cap,
        message = paste0(
          "backend = \"gretaR\" cannot fit this model (",
          cap,
          "). Use backend = \"greta\" / \"inla\" / ",
          "\"brms\"."
        ),
        family_class = "flexybayes_gretaR_refusal"
      ))
    }
    fit <- .backend_emit_fn("gretaR")(
      fb = fb,
      data = data,
      known_matrices = known_matrices,
      weights = weights,
      n_samples = n_samples,
      warmup = warmup,
      chains = chains,
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd,
      verbose = verbose,
      mcmc_verbose = mcmc_verbose,
      return_code = return_code,
      the_call = the_call,
      fixed = fixed,
      random = random,
      rcov = rcov,
      family = family,
      link = link,
      data_name = data_name
    )
    if (!isTRUE(return_code) && !is.null(fit$extras)) {
      fit$extras$fb_terms <- fb
      if (!is.null(fb_dataset_meta_for_fit)) {
        fit$extras$fb_dataset <- fb_dataset_meta_for_fit
      }
    }
    return(fit)
  }

  # Design-memory preflight. Runs only above the 1e5-row threshold
  # so the v0.2 user surface is unchanged at small scales -- the
  # existing test suite never trips this branch. Above the
  # threshold, .fb_preflight() reads only fb_dataset metadata
  # (n_rows, dictionaries, col_types) and refuses upstream with a
  # structured error if any term's design-memory estimate exceeds
  # the active ceiling. The refusal short-circuits dispatch before
  # any backend code runs.
  # The preflight return is captured and threaded into
  # lgm_gate(fb, preflight) below so the memory_feasibility_inla
  # rule can fire on INLA-specific memory infeasibility.
  # .maybe_preflight() still raises on the hard-ceiling refusal
  # (backward compat for test-fb-preflight-dispatch.R); the soft
  # INLA-only ceiling is handled by the gate rule.
  preflight_result <- .maybe_preflight(
    fb = fb,
    data = data,
    the_call = the_call,
    known_matrices = known_matrices
  )

  # Aggregated-gaussian gate. Consumes the
  # `<fb_aggregation_plan>` if the user opted in (aggregate = "auto"
  # silently routes when the plan declares eligibility; aggregate =
  # TRUE refuses with the plan's reason_codes when ineligible;
  # aggregate = FALSE skips this gate entirely). When the gate
  # accepts, builds the `<fb_aggregated>` sufficient-statistics
  # object and dispatches to emit_gaussian_aggregated() on the
  # resolved backend, short-circuiting the per-row paths below.
  if (!isFALSE(aggregate) && !isTRUE(return_code)) {
    agg_fit <- .maybe_aggregate_gaussian(
      fb = fb,
      data = data,
      backend = backend,
      aggregate = aggregate,
      n_samples = n_samples,
      warmup = warmup,
      chains = chains,
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd,
      verbose = verbose,
      mcmc_verbose = mcmc_verbose,
      the_call = the_call,
      fixed = fixed,
      random = random,
      rcov = rcov,
      family = family,
      link = link,
      data_name = data_name
    )
    if (!is.null(agg_fit)) {
      if (!is.null(agg_fit$extras) && !is.null(fb_dataset_meta_for_fit)) {
        agg_fit$extras$fb_dataset <- fb_dataset_meta_for_fit
      }
      return(agg_fit)
    }
  }

  # low_rank_smooth is a greta-backend approximation (v0.4.0). It
  # truncates the dense mgcv smooth
  # basis, which only the greta emit path materialises; INLA
  # represents smooths via rw2 and brms via Stan spline bases, so
  # neither can honour the truncation. An explicit inla / brms request
  # refuses (rather than silently fitting an unrelated smooth); an
  # auto request resolves to greta because the approximation dictates
  # the only engine that can produce it.
  if (length(.collect_approx(fb$random_terms)) > 0L) {
    if (backend %in% c("inla", "brms")) {
      stop(.fb_refusal_condition(
        reason_code = "low_rank_requires_greta",
        message = paste0(
          "A smooth requesting the low_rank_smooth approximation ",
          "requires the greta backend; the '",
          backend,
          "' backend ",
          "represents smooths differently (INLA via rw2, brms via ",
          "Stan spline bases) and cannot honour the rank-K basis ",
          "truncation. Re-fit with backend = \"greta\", or drop the ",
          "representation = ... argument to fit the exact smooth."
        ),
        family_class = "flexybayes_low_rank_requires_greta"
      ))
    }
    if (identical(backend, "auto")) backend <- "greta"
  }

  # Stan passthrough via brms. Opt-in only -- backend = "auto" never
  # routes to Stan because brms's first-call compile latency (30-60s)
  # would silently break the auto-dispatch transparent-fast-path
  # promise. Users wanting Stan triangulation set backend = "brms"
  # explicitly on fb_brms(); flexybayes() (the asreml entry) does
  # not advertise the value because asreml structured-covariance
  # terms (fa, us, ar1) do not translate losslessly to brms.
  if (identical(backend, "brms")) {
    # Capability gate. With flexybayes() (asreml entry) now able to
    # select backend = "brms", an asreml structured-covariance term
    # (vm/ped/fa/us/ar1) can reach this branch -- brms/Stan has no
    # lossless translation for it. Refuse structurally before the emit
    # rather than fail mid-translation. The low_rank case is already
    # handled upstream by the .collect_approx() check, so structured-cov
    # is the only capability failure that reaches here; gating explicitly
    # on that reason_code keeps the raised code registered.
    cap <- .backend_can_fit("brms", fb)
    if (
      !isTRUE(cap$ok) &&
        identical(cap$reason_code, "stan_cannot_represent_structured_cov")
    ) {
      stop(.fb_refusal_condition(
        reason_code = "stan_cannot_represent_structured_cov",
        message = paste0(
          "backend = \"brms\" (Stan) cannot represent this model: it ",
          "contains an asreml structured-covariance term (one of ",
          "vm/ped/fa/us/ar1) that has no lossless brms/Stan translation. ",
          "Re-fit with backend = \"greta\" (full MCMC), or backend = ",
          "\"inla\" when the model is latent-Gaussian feasible."
        ),
        family_class = "flexybayes_stan_cannot_represent_structured_cov"
      ))
    }
    # Emit entry-point resolved through the registry `engine` field
    # rather than a hard-coded symbol.
    fit <- .backend_emit_fn("brms")(
      fb = fb,
      data = data,
      known_matrices = known_matrices,
      weights = weights,
      n_samples = n_samples,
      warmup = warmup,
      chains = chains,
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd,
      verbose = verbose,
      mcmc_verbose = mcmc_verbose,
      return_code = return_code,
      the_call = the_call,
      fixed = fixed,
      random = random,
      rcov = rcov,
      family = family,
      link = link,
      data_name = data_name
    )
    if (return_code) {
      return(fit)
    }
    if (!is.null(fit$extras)) {
      # Explicit-brms request bypasses the policy; rejected_routes
      # is empty.
      fit$extras$backend_decision <- .build_routing_decision(
        backend = "brms",
        path = "explicit_brms",
        gate_checks = NULL,
        reason = paste0(
          "user requested brms explicitly ",
          "(Stan passthrough; lgm_gate() not ",
          "consulted on this path)."
        ),
        preflight_summary = preflight_result,
        representation_plan = .representation_plan_from_preflight(
          preflight_result
        ),
        rejected_routes = list()
      )
      fit$extras$fb_terms <- fb
      if (!is.null(fb_dataset_meta_for_fit)) {
        fit$extras$fb_dataset <- fb_dataset_meta_for_fit
      }
    }
    fit$exactness <- "exact"
    return(fit)
  }

  decision <- NULL
  # Availability is a registry fact, not a raw requireNamespace()
  # here. .available_backend_names() returns the active backends
  # whose required package is installed -- for inla that is the INLA
  # namespace check, so this is behaviour-equivalent while sourcing
  # the fact from the single registry authority.
  inla_installed <- "inla" %in% .available_backend_names()
  gretaR_activated <- isTRUE(getOption("flexyBayes.gretaR_activated", FALSE))
  if (backend %in% c("inla", "auto")) {
    # Gate receives preflight so the memory_feasibility_inla rule
    # can fire on INLA-specific memory infeasibility. preflight is
    # NULL below the 1e5-row threshold; the rule then trivially
    # passes per its v0.3.5 backward-compatible default.
    gated <- lgm_gate(fb, preflight = preflight_result)

    if (is_lgm_refusal(gated)) {
      if (backend == "inla") {
        msg <- paste(utils::capture.output(print(gated)), collapse = "\n")
        stop("backend = \"inla\" refused by lgm_gate():\n", msg, call. = FALSE)
      }
      primary <- gated$failures[[1L]]
      # Distinguish memory vs structural refusal for routing trace.
      is_memory <- identical(primary$rule_id, "memory_feasibility_inla")
      gate_outcome <- if (is_memory) "refuse_memory" else "refuse_structural"
      # One-time per session (auto is now the default, so a per-call
      # note on every non-LGM model would nag). Silenceable outright
      # via the option; the routing trace always carries the decision
      # on backend_decision()$rejected_routes.
      if (
        !isTRUE(getOption(
          "flexyBayes.silence_auto_fallback_note",
          FALSE
        )) &&
          !.emit_state_get("auto_fallback_note")
      ) {
        message(
          "backend = \"auto\": lgm_gate() refused (",
          primary$rule_id,
          ": ",
          primary$reason,
          "); falling back to greta. Pass backend = ",
          "\"greta\" to silence, or backend = \"inla\" ",
          "to force the refusal as an error."
        )
        .emit_state_set("auto_fallback_note", TRUE)
      }
      routing <- .resolve_routing(
        user_request = "auto",
        gate_outcome = gate_outcome,
        preflight_outcome = NA_character_,
        inla_installed = inla_installed,
        gretaR_activated = gretaR_activated
      )
      decision <- .build_routing_decision(
        backend = "greta",
        path = if (is_memory) {
          "auto_lgm_refuse_memory"
        } else {
          "auto_lgm_refuse"
        },
        gate_checks = gated$failures,
        reason = paste0(
          "lgm_gate() refused (",
          primary$rule_id,
          "); auto fell back to greta."
        ),
        preflight_summary = preflight_result,
        representation_plan = .representation_plan_from_preflight(
          preflight_result
        ),
        rejected_routes = routing$rejected_routes
      )
    } else {
      if (!inla_installed) {
        if (backend == "auto") {
          if (
            !isTRUE(getOption(
              "flexyBayes.silence_auto_inla_missing_note",
              FALSE
            )) &&
              !.emit_state_get("auto_inla_missing_note")
          ) {
            message(
              "backend = \"auto\": INLA is not installed; ",
              "routing to greta. Install INLA from ",
              "https://inla.r-inla-download.org for the LGM ",
              "fast path. Silence via ",
              "options(flexyBayes.silence_auto_inla_missing_note ",
              "= TRUE)."
            )
            .emit_state_set("auto_inla_missing_note", TRUE)
          }
          routing <- .resolve_routing(
            user_request = "auto",
            gate_outcome = "accept",
            preflight_outcome = "clear",
            inla_installed = FALSE,
            gretaR_activated = gretaR_activated
          )
          decision <- .build_routing_decision(
            backend = "greta",
            path = "auto_inla_unavailable",
            gate_checks = gated$capabilities,
            reason = paste0(
              "lgm_gate() accepted; INLA not installed; ",
              "routed to greta."
            ),
            preflight_summary = preflight_result,
            representation_plan = .representation_plan_from_preflight(
              preflight_result
            ),
            rejected_routes = routing$rejected_routes
          )
        }
      }

      if (
        backend == "inla" ||
          (backend == "auto" && inla_installed)
      ) {
        # Explicit backend = "inla" surfaces any
        # emit_inla() error directly (the user asked for INLA; the true
        # diagnostic must not be masked, and a flexyBayes-side contract
        # assertion must surface). On the AUTO path -- now the default --
        # an INLA-side numerical / runtime failure on a model the gate
        # accepted structurally falls back to greta with a one-time note,
        # because auto promised a working fit.
        inla_fit <- tryCatch(
          .backend_emit_fn("inla")(
            fb = gated,
            data = data,
            known_matrices = known_matrices,
            verbose = verbose,
            return_code = return_code,
            the_call = the_call,
            fixed = fixed,
            random = random,
            rcov = rcov,
            family = family,
            link = link,
            data_name = data_name
          ),
          error = function(e) {
            if (identical(backend, "auto")) {
              structure(list(cond = e), class = "fb_auto_inla_failure")
            } else {
              stop(e)
            }
          }
        )

        if (inherits(inla_fit, "fb_auto_inla_failure")) {
          # Fallback: numerical INLA failure on the auto path -> greta.
          if (
            !isTRUE(getOption(
              "flexyBayes.silence_auto_fallback_note",
              FALSE
            )) &&
              !.emit_state_get("auto_inla_numerical_fallback_note")
          ) {
            message(
              "backend = \"auto\": INLA failed numerically (",
              conditionMessage(inla_fit$cond),
              "); falling back to greta. Pass backend = \"inla\" ",
              "to surface the error as such, or backend = ",
              "\"greta\" to silence."
            )
            .emit_state_set("auto_inla_numerical_fallback_note", TRUE)
          }
          decision <- .build_routing_decision(
            backend = "greta",
            path = "auto_inla_numerical_fallback",
            gate_checks = gated$capabilities,
            reason = paste0(
              "lgm_gate() accepted but INLA failed ",
              "numerically; auto fell back to greta."
            ),
            preflight_summary = preflight_result,
            representation_plan = .representation_plan_from_preflight(
              preflight_result
            ),
            rejected_routes = list()
          )
          # fall through to the greta emit below (do not return).
        } else {
          fit <- inla_fit
          if (return_code) {
            return(fit)
          }
          if (!is.null(fit$extras)) {
            inla_routing <- if (backend == "auto") {
              .resolve_routing(
                user_request = "auto",
                gate_outcome = "accept",
                preflight_outcome = "clear",
                inla_installed = TRUE,
                gretaR_activated = gretaR_activated
              )
            } else {
              list(rejected_routes = list())
            } # explicit inla bypass
            fit$extras$backend_decision <- .build_routing_decision(
              backend = "inla",
              path = if (backend == "auto") {
                "auto_accept"
              } else {
                "explicit_inla_accept"
              },
              gate_checks = gated$capabilities,
              reason = "lgm_gate() accepted; INLA dispatch.",
              preflight_summary = preflight_result,
              representation_plan = .representation_plan_from_preflight(
                preflight_result
              ),
              rejected_routes = inla_routing$rejected_routes
            )
            fit$extras$fb_terms <- gated
            if (!is.null(fb_dataset_meta_for_fit)) {
              fit$extras$fb_dataset <- fb_dataset_meta_for_fit
            }
          }
          fit$exactness <- "exact"
          return(fit)
        }
      }
    }
  } else {
    # Explicit greta request -- policy bypassed.
    decision <- .build_routing_decision(
      backend = "greta",
      path = "explicit_greta",
      gate_checks = NULL,
      reason = "user requested greta explicitly (no gate run).",
      preflight_summary = preflight_result,
      representation_plan = .representation_plan_from_preflight(
        preflight_result
      ),
      rejected_routes = list()
    )
  }

  if (!requireNamespace("greta", quietly = TRUE)) {
    stop(
      "Package 'greta' is required for the greta backend. ",
      "Install with:\n",
      "  install.packages('greta')\n",
      "  greta::install_greta_deps()",
      call. = FALSE
    )
  }

  # Emit entry-point resolved through the registry `engine` field
  # rather than a hard-coded symbol.
  fit <- .backend_emit_fn("greta")(
    fb = fb,
    data = data,
    known_matrices = known_matrices,
    weights = weights,
    n_samples = n_samples,
    warmup = warmup,
    chains = chains,
    prior_fixed_sd = prior_fixed_sd,
    prior_vc_sd = prior_vc_sd,
    verbose = verbose,
    mcmc_verbose = mcmc_verbose,
    return_code = return_code,
    the_call = the_call,
    fixed = fixed,
    random = random,
    rcov = rcov,
    family = family,
    link = link,
    data_name = data_name
  )

  if (!isTRUE(return_code) && !is.null(fit$extras)) {
    if (!is.null(decision)) {
      fit$extras$backend_decision <- decision
    }
    fit$extras$fb_terms <- fb
    if (!is.null(fb_dataset_meta_for_fit)) {
      fit$extras$fb_dataset <- fb_dataset_meta_for_fit
    }
  }
  # Three-tier exactness vocabulary (v0.4.0). A greta fit carrying
  # a registered approximation scheme
  # is labelled "approximate_<scheme>"; an exact fit keeps "exact"
  # (aggregated paths set "aggregated_exact" upstream). The label is
  # the same string validate_approximation() and the [APPROX:] display
  # badge key off.
  if (!isTRUE(return_code)) {
    approx <- fit$extras$parse_info$approx
    fit$exactness <- if (length(approx) > 0L) {
      paste0("approximate_", approx[[1L]]$scheme)
    } else {
      "exact"
    }
  }
  fit
}


# Normalise the user-facing `aggregate` argument. Accepts:
#   TRUE        -- force aggregation; refuse if plan ineligible.
#   FALSE       -- skip aggregation; always per-row.
#   "auto"      -- aggregate when plan eligible; silently fall through
#                  to per-row otherwise (the v0.3.2 default).
.normalise_aggregate <- function(x) {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) {
    return(x)
  }
  if (is.character(x) && length(x) == 1L && identical(x, "auto")) {
    return("auto")
  }
  stop(
    "`aggregate` must be TRUE, FALSE, or \"auto\"; got: ",
    paste(deparse(x), collapse = " "),
    call. = FALSE
  )
}


# Aggregated-gaussian gate. Returns either a
# fitted <flexybayes_aggregated> object (gate accepted; aggregated
# emit ran) or NULL (gate declined; caller falls through to per-row).
# Refuses (raises a typed condition) when aggregate = TRUE but the
# plan declares ineligibility, or when aggregate = TRUE on a backend
# that does not support the aggregated path (brms, gretaR, stan).
.maybe_aggregate_gaussian <- function(
  fb,
  data,
  backend,
  aggregate,
  n_samples,
  warmup,
  chains,
  prior_fixed_sd,
  prior_vc_sd,
  verbose,
  mcmc_verbose,
  the_call,
  fixed,
  random,
  rcov,
  family,
  link,
  data_name
) {
  # Aggregation only defined for backends with a flexyBayes-side
  # aggregated emit path. v0.3.3 supports greta + inla; brms / gretaR
  # / stan have no aggregated path in scope.
  agg_capable_backend <- backend %in% c("greta", "inla", "auto")
  if (!agg_capable_backend) {
    if (isTRUE(aggregate)) {
      stop(
        "`aggregate = TRUE` is not supported on backend = \"",
        backend,
        "\" (the aggregated path is wired for greta and ",
        "inla only). Pass aggregate = FALSE or switch backend.",
        call. = FALSE
      )
    }
    return(NULL)
  }

  # Run the model-level plan against the fb_dataset metadata. Reads
  # only IR + dictionaries; no design matrix is materialised.
  fb_dataset_meta <- .fb_dataset(data)
  plan <- .fb_aggregation_plan(fb, fb_dataset_meta)

  # Resolve the effective backend for the aggregated emit. v0.3.3
  # supports both INLA and greta on the aggregated path; the greta
  # path uses the call-construct-in-env workaround in
  # .emit_gaussian_aggregated_greta() to bypass greta::model()'s
  # substitute()-deparse trap under do.call().
  eff_backend <- .agg_resolve_backend(backend, fb)
  if (identical(eff_backend, "greta") && !isTRUE(aggregate)) {
    return(NULL)
  }

  # at_units on INLA: the multi-likelihood INLA stack is still
  # deferred at v0.3.3 (queued for a later minor). If we'd otherwise
  # route to INLA AND the rcov uses at_units, force the eff_backend
  # down to greta (which supports the heterogeneous-residual aggregated
  # path) or fall through to per-row when greta isn't an option.
  uses_at_units <- length(fb$rcov_terms) > 0L &&
    any(vapply(
      fb$rcov_terms,
      function(t) identical(t$type %||% "", "at_units"),
      logical(1L)
    ))
  if (identical(eff_backend, "inla") && uses_at_units) {
    if (identical(backend, "inla") && isTRUE(aggregate)) {
      stop(
        "`aggregate = TRUE` with backend = \"inla\" refused: ",
        "heterogeneous residual at_units on INLA requires the ",
        "multi-likelihood INLA stack (deferred to a future release). Pass ",
        "backend = \"greta\" with aggregate = TRUE, or aggregate ",
        "= FALSE for per-row INLA.",
        call. = FALSE
      )
    }
    if (identical(backend, "auto")) {
      eff_backend <- "greta"
    } else {
      return(NULL)
    }
  }

  if (!isTRUE(plan$eligible)) {
    if (isTRUE(aggregate)) {
      cond <- structure(
        class = c("flexybayes_aggregate_refusal", "error", "condition"),
        list(
          message = paste0(
            "`aggregate = TRUE` refused by .fb_aggregation_plan(): ",
            paste(plan$reason_codes, collapse = ", "),
            ". Pass aggregate = FALSE for the per-row path."
          ),
          call = the_call,
          reason_codes = plan$reason_codes,
          plan = plan
        )
      )
      stop(cond)
    }
    return(NULL)
  }

  # Family branch. Gaussian uses the dedicated sufficient-statistics
  # aggregator + gaussian emit; the binomial / poisson families use the
  # family-generic streaming aggregator over the in-memory data (one
  # chunk) + the count emit. Both produce the same <fb_aggregated> shape
  # downstream methods consume.
  if (identical(fb$family, "gaussian")) {
    fb_aggregated <- .fb_aggregate_gaussian(fb, fb_dataset_meta)
    fit <- emit_gaussian_aggregated(
      fb = fb,
      fb_aggregated = fb_aggregated,
      data = data,
      backend = eff_backend,
      n_samples = n_samples,
      warmup = warmup,
      chains = chains,
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd,
      verbose = verbose,
      mcmc_verbose = mcmc_verbose,
      return_code = FALSE,
      the_call = the_call,
      fixed = fixed,
      random = random,
      rcov = rcov,
      family = family,
      link = link,
      data_name = data_name
    )
  } else {
    # Exact aggregation for the count families requires unit-trial
    # binomial (a 0/1 numeric response); a non-Bernoulli binomial needs
    # per-row trial counts the cell sums cannot recover.
    if (identical(fb$family, "binomial")) {
      yv <- data[[fb$response]]
      if (!(is.numeric(yv) && all(yv %in% c(0, 1)))) {
        if (isTRUE(aggregate)) {
          stop(
            "`aggregate = TRUE` refused: a binomial response that is ",
            "not a 0/1 numeric vector needs per-row trial counts the ",
            "cell sums cannot recover. Pass aggregate = FALSE.",
            call. = FALSE
          )
        }
        return(NULL)
      }
    }
    src <- .fb_stream_source(data, chunk_rows = max(nrow(data), 1L))
    fb_aggregated <- .fb_stream_aggregate(
      fb,
      src,
      trials = NULL,
      exposure = NULL,
      verbose = FALSE
    )
    fit <- emit_count_aggregated(
      fb = fb,
      fb_aggregated = fb_aggregated,
      data = data,
      backend = eff_backend,
      n_samples = n_samples,
      warmup = warmup,
      chains = chains,
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd,
      verbose = verbose,
      mcmc_verbose = mcmc_verbose,
      the_call = the_call,
      fixed = fixed,
      random = random,
      rcov = rcov,
      family = family,
      link = link,
      data_name = data_name
    )
  }

  fit$exactness <- "aggregated_exact"
  if (!is.null(fit$extras)) {
    # The aggregated-gaussian path's chosen backend was already
    # resolved by .agg_resolve_backend() above; rejected_routes is
    # empty (the aggregation gate itself bypasses the routing policy
    # table -- this is its own dispatch path). preflight is NULL here
    # because aggregation runs upstream of the >1e5 preflight
    # threshold in practice; if invoked under the metadata-only
    # path the slot stays NULL and reproducibility checks treat
    # it as legitimately absent.
    fit$extras$backend_decision <- .build_routing_decision(
      backend = eff_backend,
      path = if (identical(fb$family, "gaussian")) {
        "aggregated_gaussian"
      } else {
        "aggregated_count"
      },
      gate_checks = NULL,
      reason = sprintf(
        "aggregation plan eligible (N = %d, K = %d, ratio = %.2f:1)",
        plan$N,
        plan$K_est,
        plan$N / plan$K_est
      ),
      preflight_summary = NULL,
      representation_plan = NULL,
      rejected_routes = list()
    )
    fit$extras$fb_terms <- fb
  }
  fit
}


# Choose the effective aggregated-emit backend. greta / inla pass
# through verbatim. auto consults lgm_gate -- if the model is LGM-
# compatible and INLA is installed, prefer INLA (faster deterministic
# on the aggregated path); else greta.
.agg_resolve_backend <- function(backend, fb) {
  if (backend %in% c("greta", "inla")) {
    return(backend)
  }
  if (identical(backend, "auto")) {
    # Let a genuine lgm_gate() fault propagate (charter: never swallow an
    # error into a silent fallback). A refusal is a normal return value
    # (an `lgm_refusal` object), not an error, and steers the choice to
    # greta below; only an actual fault would have been masked here.
    gated <- lgm_gate(fb)
    if (!is_lgm_refusal(gated) && requireNamespace("INLA", quietly = TRUE)) {
      return("inla")
    }
    return("greta")
  }
  # Unknown backend -- caller already filtered, defensive default.
  "greta"
}


# Shared preflight call. Runs IFF the IR carries
# n_rows >= 1e5; below that threshold the v0.2 user surface is
# unchanged. Returns the <fb_preflight> result (or NULL when below
# the threshold); on refusal raises the typed condition via
# .stop_preflight_refusal(). Used by both .dispatch_backend()'s
# entry block and the review_code branches in flexybayes() /
# fb_brms() so the gate fires identically across both paths.
.maybe_preflight <- function(
  fb,
  data,
  the_call,
  known_matrices = NULL,
  threshold = 1e5L
) {
  if (
    !isTRUE(
      !is.null(fb$data_summary$n) &&
        fb$data_summary$n >= threshold
    )
  ) {
    return(NULL)
  }
  # Optional global ceiling override -- lets the user set a stricter
  # cap than the default 60% x RAM (e.g. on a shared host) without
  # passing an argument through every dispatch entry. NULL falls back
  # to the default RAM-probe resolution inside .fb_preflight(). The
  # related option `flexyBayes.preflight_ram_fraction` adjusts the
  # default fraction multiplicatively.
  pf_ceiling_gb <- getOption("flexyBayes.preflight_ceiling_gb", NULL)
  pf_dataset <- .fb_dataset(data)
  pf_result <- .fb_preflight(
    fb,
    pf_dataset,
    memory_ceiling_gb = pf_ceiling_gb,
    known_matrices = known_matrices
  )
  if (inherits(pf_result$refusal, "fb_preflight_refusal")) {
    .stop_preflight_refusal(pf_result$refusal, call = the_call)
  }
  pf_result
}


# Raise a typed preflight refusal as a structured condition. The
# custom class `flexybayes_preflight_refusal` carries the refusal
# object verbatim so downstream tooling can pattern-match on the
# slots without parsing free text.
#
# The headline is branched on refusal$reason_code so that an
# unknown-representation refusal does not falsely suggest raising
# the memory ceiling. The detailed body already prints correctly
# via print.fb_preflight_refusal().
.stop_preflight_refusal <- function(refusal, call) {
  msg <- paste(utils::capture.output(print(refusal)), collapse = "\n")
  headline <- switch(
    refusal$reason_code,
    design_memory_exceeds_ceiling = paste0(
      "flexyBayes preflight refused: the design exceeds the ",
      "active memory ceiling. The dispatch was short-circuited ",
      "before any backend code ran."
    ),
    representation_unknown_for_preflight = paste0(
      "flexyBayes preflight refused: the design representation ",
      "is not characterised by the preflight estimator. The ",
      "dispatch was short-circuited before any backend code ",
      "ran -- raising the memory ceiling will not help."
    ),
    # default: neutral wording carrying the reason code so any future
    # refusal class lands honestly even before its headline is added.
    paste0(
      "flexyBayes preflight refused (reason_code = ",
      refusal$reason_code,
      "). The dispatch was short-circuited ",
      "before any backend code ran."
    )
  )
  stop(.fb_refusal_condition(
    reason_code = refusal$reason_code,
    message = paste0(headline, "\n", msg),
    family_class = "flexybayes_preflight_refusal",
    call = call,
    refusal = refusal,
    binding_term = refusal$binding_term,
    ceiling_bytes = refusal$ceiling_bytes
  ))
}
