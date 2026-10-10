# .dispatch_backend() -- shared dispatch helper.
#
# Shared so that both flexybayes() (asreml entry) and fb_brms()
# (brms entry) drive the same backend dispatch. Byte-equivalent
# under snapshot guard (test-backend-auto.R).
#
# Routes to emit_inla() or emit_brms() per backend = c("inla", "brms",
# "auto"), populates fit$extras$backend_decision with the
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
# gate / preflight columns -- a "brms" request matches regardless of
# gate / preflight outcomes.
#
# Reason-code vocabulary:
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
  } else if (
    is.list(backend) &&
      !inherits(backend, "fb_engine") &&
      length(backend)
  ) {
    backend
  } else {
    NULL
  }
}


.routing_policy_table <- function() {
  data.frame(
    user_request = c(
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
      "clear",
      "clear",
      NA_character_,
      NA_character_,
      NA_character_
    ),
    inla_installed = c(NA, NA, NA, TRUE, FALSE, NA, NA, NA),
    chosen_backend = c(
      "brms",
      "inla",
      NA_character_, # explicit-inla refusal raises
      "inla",
      "pending_fallback",
      "pending_fallback",
      "pending_fallback",
      "pending_fallback"
    ), # auto + accept + skipped preflight
    reason_code = c(
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


# .resolve_routing() -- given the three input keys and the INLA-install
# environmental flag, returns the canonical chosen_backend +
# rejected_routes per the `.routing_policy_table`. Used to populate the
# trace at dispatch time; the actual backend invocation lives in
# `.dispatch_backend()` below.
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
  fb = NULL
) {
  is_explicit <- user_request %in% c("brms", "inla")
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
    # No policy row matches -- refuse rather than fall back to an
    # unrepresented engine. This branch should never fire in production
    # (the table covers every conceivable tuple); if it does, the
    # test suite catches it via subtest (n).
    return(list(
      chosen_backend = NA_character_,
      reason_code = "policy_table_no_match_refuse",
      rejected_routes = list()
    ))
  }
  row <- table[which(hits)[1L], ]
  chosen_backend <- row$chosen_backend
  reason_code <- row$reason_code

  # "pending_fallback" is a sentinel, never a real backend name: every
  # auto row that reaches it (INLA refused / unavailable) re-resolves to
  # the faithful active fallback -- brms when it can represent the model
  # (the multi-stratum designed-experiment case keeps its specific
  # reason for the trace), otherwise no active route and the caller
  # refuses (A3). Explicit requests bypass this (is_explicit); by
  # construction no explicit row ever produces this sentinel.
  if (!is_explicit && identical(chosen_backend, "pending_fallback")) {
    fallback <- if (!is.null(fb)) {
      .resolve_auto_fallback_backend(fb)
    } else {
      NA_character_
    }
    if (identical(fallback, "brms")) {
      chosen_backend <- "brms"
      reason_code <- if (!is.null(fb) && .has_interaction_random_term(fb)) {
        "auto_multistratum_to_brms"
      } else {
        "auto_fallback_to_brms"
      }
    } else {
      chosen_backend <- NA_character_
      reason_code <- "auto_no_active_route"
    }
  }

  # Rejected routes: empty for explicit user requests (policy
  # bypassed); otherwise enumerate the non-chosen backends with
  # per-backend reason codes derived from the same routing inputs.
  # brms is a candidate only via the multi-stratum override above.
  rejected <- list()
  if (!is_explicit) {
    chosen <- chosen_backend
    candidates <- c("inla")
    for (cand in setdiff(candidates, chosen)) {
      cand_reason <- switch(
        cand,
        "inla" = .inla_rejection_reason(gate_outcome, inla_installed)
      )
      rejected[[length(rejected) + 1L]] <- list(
        backend = cand,
        reason = cand_reason
      )
    }
  }

  list(
    chosen_backend = chosen_backend,
    reason_code = reason_code,
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
  # User may pass `c("inla", "brms", "auto")` (the default arg
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
        "Re-route via backend = \"inla\" (Laplace; \"exact\"), or ",
        "backend = \"brms\" (Stan via brms; \"exact\")."
      ),
      family_class = "flexybayes_approximate_route_refusal",
      backend = b
    ))
  }
  invisible(NULL)
}


# .check_known_backend_name() -- entry-function guard called from
# `flexybayes()`, `fb_plan()`, and `flexybayes_stream()` BEFORE
# match.arg(backend) so an unrecognised backend name -- including a
# formerly-registered engine name this package has withdrawn entirely --
# refuses with the structured `unknown_backend` reason naming the active
# engines, rather than match.arg()'s generic "should be one of" error.
#
# Tolerates a multi-value `backend` (the unmodified default arg vector,
# e.g. `c("auto", "inla", "brms")`): only a genuine single-value
# selection is checked, mirroring `.check_approximate_scheme()`'s own
# defensive pattern.
.check_known_backend_name <- function(backend, allowed) {
  if (!is.character(backend) || length(backend) != 1L) {
    return(invisible(NULL))
  }
  if (is.na(backend) || !nzchar(backend) || backend %in% allowed) {
    return(invisible(NULL))
  }
  stop(.fb_refusal_condition(
    reason_code = "unknown_backend",
    message = paste0(
      "backend = \"",
      backend,
      "\" is not a recognised flexyBayes engine. The active engines are ",
      "\"inla\" and \"brms\" (or \"auto\" to let flexyBayes choose)."
    ),
    family_class = "flexybayes_unknown_backend_refusal",
    backend = backend
  ))
}


# ---------------------------------------------------------------- #

# .has_interaction_random_term() -- TRUE when the model carries a nested
# (A:B) or combo (A:B:C) interaction random intercept: the multi-stratum
# designed-experiment shape (gen:env + env:rep + env:rep:block).
.has_interaction_random_term <- function(fb) {
  any(vapply(
    fb$random_terms,
    function(t) isTRUE((t$type %||% "") %in% c("nested", "combo")),
    logical(1L)
  ))
}

# .resolve_auto_fallback_backend() -- when the routing policy would place an
# auto request on the pending-fallback sentinel (INLA refused /
# unavailable), pick the faithful active fallback: brms when it is
# installed AND can represent the model -- otherwise NA, meaning no
# active backend can faithfully fit it and the caller refuses (A3: "no
# faithful route => refuse"). A structured-covariance term (fa / us / ar1)
# has no lossless brms translation, so those models return NA here (a named
# transitional refusal until WP16 lands AR1xAR1 spatial on INLA).
.resolve_auto_fallback_backend <- function(fb) {
  if (
    requireNamespace("brms", quietly = TRUE) &&
      isTRUE(.backend_can_fit("brms", fb)$ok)
  ) {
    return("brms")
  }
  NA_character_
}


# .fb_ir_weights() --- the model's observation-weight vector, read from
# the IR's addition_terms (every ingest adapter records weights there
# under type = "weights"). Shared by the dispatch-level refusal guard
# below and both per-row emits (R/emit_brms.R, R/emit_inla.R), so all
# three read the exact same extraction -- an IR built upstream and a
# direct dispatch-level `weights =` argument cannot disagree about
# whether a weighted fit was requested.
#
# @noRd
# @keywords internal
.fb_ir_weights <- function(fb) {
  for (t in fb$addition_terms %||% list()) {
    if (identical(t$type %||% "", "weights")) {
      return(t$values)
    }
  }
  NULL
}

# .fb_weights_nonconstant() --- TRUE when `w` is a genuine (non-
# constant) weight vector. An all-equal vector is the unweighted fit
# under a different spelling and is let through everywhere a real
# weight vector would be checked -- this keeps `weights = rep(1, n)`
# working, exactly as it did before C6.
#
# @noRd
# @keywords internal
.fb_weights_nonconstant <- function(w) {
  !is.null(w) && is.numeric(w) && length(w) > 0L && !isTRUE(all(w == w[[1L]]))
}

# .refuse_unsupported_weights() --- guard for `weights` on a family /
# link this package does not lower it for.
#
# C6 lowers weights for family = "gaussian" with the identity link
# only: INLA's `scale = w` per-observation precision multiplier and
# brms's sigma-distributional offset (R/emit_brms.R's
# .FB_BRMS_WEIGHTS_OFFSET_COL -- NOT brms's own `weights()` addition
# term, which implements a different, likelihood-power quantity; see
# that constant's definition for the full grounding trail) both give
# Var(y_i) = sigma^2 / w_i on the Gaussian likelihood (the "documented
# ASReml sense" the original all-families refusal named as the one
# unambiguous mapping -- see .populate_refusal_registry_
# field_hardening()'s weights_not_supported comment). Any other family,
# or a non-identity link on gaussian (scale = w's meaning is tied to
# the identity-link Gaussian likelihood, not a transformed one),
# refuses by name naming the family. Reads both the IR's addition
# terms and the dispatch-level argument, so neither an IR built
# upstream nor a direct dispatch call can slip a weighted request past
# the active emitters.
.refuse_unsupported_weights <- function(fb, weights = NULL) {
  w <- .fb_ir_weights(fb) %||% weights
  if (!.fb_weights_nonconstant(w)) {
    return(invisible(NULL))
  }
  fam <- fb$family %||% ""
  link_ok <- is.null(fb$link) || identical(fb$link, "identity")
  if (identical(fam, "gaussian") && link_ok) {
    return(invisible(NULL))
  }
  stop(.fb_refusal_condition(
    reason_code = "weights_requires_gaussian",
    message = paste0(
      "`weights` is lowered for family = \"gaussian\" with the identity ",
      "link only (INLA's `scale = w` per-observation precision ",
      "multiplier; brms's sigma-distributional offset -- both give ",
      "Var(y_i) = sigma^2 / w_i, the same Gaussian likelihood weighting ",
      "on either engine). This model is family = \"", fam, "\"",
      if (!link_ok) paste0(", link = \"", fb$link, "\"") else "",
      ". Drop `weights`, or pre-scale the response if the intended ",
      "semantics is exactly Var(y_i) = sigma^2 / w_i."
    ),
    family_class = "flexybayes_weights_requires_gaussian_refusal",
    family = fam,
    link = fb$link
  ))
}

# .refuse_non_binary_binomial_response() --- the binomial-response guard.
#
# D14 / WP-D. flexybayes() / fb_inla() / fb_brms() (the main,
# non-streaming entry) have no cbind(success, failure) or trials =
# spelling -- that exists only on flexybayes_stream() -- so
# family = "binomial" on this entry requires a numeric 0/1 response.
# Left unguarded, a response outside {0, 1} reached the engine on every
# route and failed there raw: INLA's own subprocess dies ("Binomial
# data ... is void"), which C2's wrapper then reports as
# inla_program_failed with a design-memory diagnosis that has nothing to
# do with the actual cause; brms raises an untyped condition; the
# aggregated route (R/dispatch.R's own .maybe_aggregate_gaussian() call
# site) already had an ad hoc untyped stop() for the identical fact.
# One typed refusal at the dispatch choke point every route passes
# through replaces all three.
.refuse_offset_route <- function(fb, route) {
  offsets <- fb$offsets %||% character(0)
  if (!length(offsets)) {
    return(invisible(NULL))
  }
  label <- paste0("`", offsets, "`", collapse = ", ")
  stop(.fb_refusal_condition(
    reason_code = "fixed_offset_unsupported",
    message = paste0(
      "The fixed formula carries ", label, ", and ", route,
      " does not carry an offset."
    ),
    offsets = offsets,
    route = route
  ))
}

.refuse_non_binary_binomial_response <- function(fb, data) {
  if (!identical(fb$family %||% "", "binomial")) {
    return(invisible(NULL))
  }
  col <- fb$response
  yv <- .fb_response_values(fb, data)
  if (is.null(yv)) {
    return(invisible(NULL)) # nothing to check against
  }
  if (is.numeric(yv) && all(yv %in% c(0, 1) | is.na(yv))) {
    return(invisible(NULL))
  }

  offending <- if (is.numeric(yv)) {
    sort(unique(yv[!is.na(yv) & !(yv %in% c(0, 1))]))
  } else {
    unique(yv)
  }
  values_shown <- paste(utils::head(as.character(offending), 6L),
    collapse = ", "
  )
  more_note <- if (length(offending) > 6L) {
    paste0(", and ", length(offending) - 6L, " more")
  } else {
    ""
  }
  what <- if (is.numeric(yv)) {
    paste0("the value(s) ", values_shown, more_note, " outside {0, 1}")
  } else {
    paste0(
      "class `", class(yv)[[1L]], "`, not a numeric 0/1 vector ",
      "(value(s) ", values_shown, more_note, ")"
    )
  }

  stop(.fb_refusal_condition(
    reason_code = "binomial_response_not_binary",
    message = paste0(
      "family = \"binomial\" requires a numeric 0/1 response on this ",
      "entry point. Column `", col, "` carries ", what, ". There is no ",
      "cbind(success, failure) or trials = spelling on flexybayes() / ",
      "fb_inla() / fb_brms(); pass pre-aggregated trial counts through ",
      "flexybayes_stream(..., trials = ...) instead, or recode `", col,
      "` to a numeric 0/1 vector for a per-row binary outcome."
    ),
    family_class = "flexybayes_binomial_response_not_binary_refusal",
    column = col,
    offending_values = offending
  ))
}

# .refuse_unrepresentable_structures() --- the named-structure guard.
#
# Some structures parse into the IR so the formula catalogue can describe
# them, and no active engine emits them. Left alone they surface as an
# lgm_gate() print on INLA and an untyped stop() on brms, which is two
# different refusal vocabularies for one fact. This guard states the fact
# once, at the dispatch choke point every fitting and code-inspection route
# passes through, with a reason code the caller can catch and an alternative
# the caller can act on.
#
# A residual-side AR1 reaches here when it is not the separable form on a
# Gaussian identity-link response, which R/fb.R carries as a nugget-free
# field: the 1-D and grouped forms, and any other family.
.refuse_unrepresentable_structures <- function(fb) {
  for (term in fb$residual_terms %||% list()) {
    if ((term$type %||% "") %in% c("ar1", "ar1_spatial", "ar1_spatial_at")) {
      stop(.fb_refusal_condition(
        reason_code = "ar1_residual_not_representable",
        message = .ar1_residual_refusal_message(term)
      ))
    }
  }

  all_terms <- c(fb$random_terms %||% list(), fb$residual_terms %||% list())
  for (term in all_terms) {
    ty <- term$type %||% ""
    # fa_gxe is emitted on INLA from 0.10.0.9000: each common factor is one
    # latent field entered per level as a SCALED COPY, and the copy's
    # scaling parameter is the loading. The claim this refusal used to make
    # -- that INLA has no factor-analytic latent model -- is true of INLA's
    # named model list and false of what INLA can express, which is the
    # distinction that matters to a user. The bare `fa(v, k)` spelling (no
    # inner grouping factor) has no emit on either engine and still refuses.
    if (identical(ty, "fa")) {
      stop(.fb_refusal_condition(
        reason_code = "fa_not_representable",
        message = paste0(
          "fa(", term$var, ", k) names a factor-analytic covariance with no ",
          "grouping factor. flexyBayes does not fit a factor-analytic ",
          "covariance in any form. Write diag(", term$var, "):<group> for ",
          "one variance per level, or us(", term$var, "):<group> on brms ",
          "for the full covariance."
        )
      ))
    }
    # INLA can represent fa_gxe; its fits of simulated trials with a known
    # covariance missed the coverage required in 0.18.0.
    #
    # The altitude is dispatch-wide rather than the INLA emit because the
    # limit is not INLA's representation: brms has no factor-analytic lowering
    # at all, so there is no second engine to route to, and a per-engine
    # refusal would send backend = "auto" to brms to collect a second, less
    # informative refusal. One vocabulary, stated once, before any route.
    if (identical(ty, "fa_gxe")) {
      stop(.fb_refusal_condition(
        reason_code = "factor_analytic_not_repeatable",
        message = paste0(
          "fa(", term$outer, ", ", term$k, "):", term$inner,
          " asks for a factor-analytic genotype-by-environment covariance, ",
          "which flexyBayes does not fit on either engine. In a test fixed ",
          "before it ran, INLA fits of simulated trials with a known ",
          "covariance missed the coverage required: ",
          "genetic-correlation intervals covered the true value 134 of 150 ",
          "times for k = 1 and 125 of 150 for k = 2, against 136 required, ",
          "and for k = 2 the genetic correlations were biased upward. ",
          "flexyBayes has no factor-analytic lowering for brms. Write ",
          "diag(", term$outer, "):", term$inner, ", one variance per ",
          "environment, which fits on both engines, or us(", term$outer,
          "):", term$inner, " on brms, which estimates every correlation. ",
          "On INLA, summary(fit)$varcomp marks some variances the data ",
          "barely identify as `weakly identified`."
        )
      ))
    }
    if (identical(ty, "ar2")) {
      stop(.fb_refusal_condition(
        reason_code = "ar2_not_representable",
        message = paste0(
          "ar2(", term$var, ") asks for a second-order autoregressive ",
          "structure. flexyBayes writes first-order autoregressive terms ",
          "only, on INLA. Use ar1(", term$var, ") for the first-order ",
          "process; for a second-order one, INLA's own ar model takes ",
          "higher orders, and ASReml fits ar2()."
        )
      ))
    }
    if (identical(ty, "cor")) {
      stop(.fb_refusal_condition(
        reason_code = "corh_no_equicorrelation_representation",
        message = paste0(
          "cor(", term$var, ") asks for one correlation shared by every ",
          "pair of levels. flexyBayes has no lowering for it on either ",
          "engine. Use us(", term$var, "):g on brms to estimate every ",
          "correlation, or diag(", term$var, "):g for separate variances ",
          "without correlation."
        )
      ))
    }
    # corh(f):g is refused here, on every engine, rather than by the brms
    # emit alone: the INLA gate used to refuse it as an emit gap and point
    # at brms, which then refused it too.
    if (identical(ty, "corh_gxe")) {
      stop(.fb_refusal_condition(
        reason_code = "corh_no_equicorrelation_representation",
        message = .fb_corh_refusal_message(term)
      ))
    }
    if (identical(ty, "str")) {
      stop(.fb_refusal_condition(
        reason_code = "str_not_representable",
        message = paste0(
          "str() declares one covariance structure across a group of ",
          "terms. flexyBayes has no emit path for it on either active ",
          "engine, and the terms inside it would otherwise be fitted ",
          "independently -- a different model. Write the structure ",
          "directly on the term that carries it (us(f):g, diag(f):g, ",
          "ar1(t)), or use ASReml for str()."
        )
      ))
    }
    if (identical(ty, "interaction_generic")) {
      left <- term$left$type %||% "?"
      right <- term$right$type %||% "?"
      multi_trait <- identical(left, "us") && right %in% c("vm", "ped")
      stop(.fb_refusal_condition(
        reason_code = "interaction_not_representable",
        message = paste0(
          "flexyBayes cannot fit the interaction of a \"", left,
          "\" term with a \"", right, "\" term. It matches none of the ",
          "structured interactions flexyBayes recognises: fa(f, k):g, ",
          "us(f):g, corh(f):g, diag(f):g, at(f):units, ar1(r):ar1(c), ",
          "vm(g, K):f and the crossings a:b and a:b:c. ",
          if (multi_trait) {
            paste0(
              "flexyBayes has no lowering for a multi-trait covariance, ",
              "us(trait):vm(gen); fit it with brms directly, as ",
              "(0 + trait | gr(gen, cov = K)), or fit each trait separately."
            )
          } else {
            "Write the structure as one of these."
          }
        )
      ))
    }
  }

  # A term written on both sides of the model is not an unrepresentable
  # structure -- every engine will attempt it -- but it is the same class
  # of fact: a model shape named at the gate rather than discovered
  # downstream. Sharing the guard keeps one call site for both.
  .refuse_term_in_fixed_and_random(fb)

  invisible(NULL)
}

# .refuse_term_in_fixed_and_random() --- the aliasing guard.
#
# A factor written on both the fixed and the random side is confounded
# with itself. The fixed part already spends one parameter per level, so
# the random copy's deviations carry no information the fixed means have
# not already taken: they are identified by their prior, and the variance
# component reported for them reads that prior rather than the data.
#
# It is refused rather than warned about because no practitioner intends
# it and every engine answers it badly in a different way. On INLA the
# marginal solve goes intermittently singular -- measured as an
# empty variance-component table on 4 of 8 small-n runs, and as two
# `Segmentation fault: 11` lines from the inla.run child followed by an
# untyped error at benchmark scale. brms completes, and returns a
# posterior whose two term sets are estimating one quantity between them.
# ASReml accepts the spelling and fits it, so a translated script arrives
# here unchanged and deserves to be told why.
#
# Scope is deliberately the plain aliased pair: a factor MAIN EFFECT on
# the fixed side against an intercept-only random main effect on the same
# variable (`f`, or `id(f)` / `ide(f)`). A random slope over a fixed
# grouping factor -- `y ~ g + (x || g)` -- is a model a user does mean,
# and an interaction such as `random = ~ f:b` is not aliased with the
# main effect at all, so neither fires.
.refuse_term_in_fixed_and_random <- function(fb) {
  fixed_factors <- unlist(lapply(
    fb$fixed_terms %||% list(),
    function(t) if (identical(t$type %||% "", "factor")) t$var else NULL
  ))
  if (length(fixed_factors) == 0L) {
    return(invisible(NULL))
  }

  random_mains <- unlist(lapply(
    fb$random_terms %||% list(),
    function(t) {
      if ((t$type %||% "") %in% c("simple", "ide")) t$var else NULL
    }
  ))
  clash <- intersect(as.character(fixed_factors), as.character(random_mains))
  if (length(clash) == 0L) {
    return(invisible(NULL))
  }

  entry <- .lookup_refusal("term_in_fixed_and_random")
  stop(.fb_refusal_condition(
    reason_code = "term_in_fixed_and_random",
    message = sprintf(entry$message_template, clash[[1L]])
  ))
}


# The refusal text for a residual-side AR1 flexyBayes does not fit. A
# separable one reaches here only off the Gaussian identity link.
.ar1_residual_refusal_message <- function(term) {
  if (
    identical(term$type, "ar1_spatial_at") ||
      (identical(term$type, "ar1_spatial") && isTRUE(term$col_ar1) &&
        is.null(term$at_var))
  ) {
    spelling <- if (identical(term$type, "ar1_spatial_at")) {
      paste0(
        "dsum(~ ar1(", term$row_var, "):ar1(", term$col_var, ") | ",
        term$by_var, ")"
      )
    } else {
      .ar1_term_spelling(term)
    }
    return(paste0(
      "residual = ~ ", spelling, " asks for ASReml's nugget-free ",
      "separable residual, which flexyBayes fits for a Gaussian response ",
      "with the identity link only. For this family, write the field on ",
      "the random side, random = ~ ar1(", term$row_var, "):ar1(",
      term$col_var, "), which fits on INLA."
    ))
  }
  spelling <- .ar1_term_spelling(term)
  asreml_params <- if (identical(term$type, "ar1_spatial")) {
    paste0(
      "two parameters (the correlation along ", term$row_var,
      " and one variance, shared by every level of ", term$col_var, ")"
    )
  } else {
    paste0(
      "two parameters (the correlation along ", term$var, " and one variance)"
    )
  }
  paste0(
    "residual = ~ ", spelling, " asks for ASReml's autoregressive ",
    "residual: one correlated residual with no independent plot error, ",
    asreml_params, ". flexyBayes fits this process as a latent field plus ",
    "the Gaussian observation error, three parameters, which is a ",
    "different model and gives a different correlation when the data ",
    "carry plot-to-plot noise. Write the field on the random side, ",
    "random = ~ ", spelling, ", which fits that model under its own name ",
    "on INLA (one observation per grid node), or use ASReml for the ",
    "residual form (Gilmour, Cullis, & Verbyla, 1997). The separable ",
    "residual ~ ar1(row):ar1(col) is fitted without a nugget."
  )
}


# .note_sampler_args_ignored() --- say so when `seed` / `control` reach an
# engine that does not consume them.
#
# Both arguments are Stan settings threaded from flexybayes() to whichever
# engine dispatch resolves. Only the brms emit passes them on. Accepting a
# seed on a deterministic or non-Stan engine and discarding it silently
# would let a user believe a fit was pinned when nothing was; the message
# is emitted once per fit and is silenceable through the same option as
# the other engine notes.
.note_sampler_args_ignored <- function(engine, seed = NULL, control = NULL) {
  if (is.null(seed) && is.null(control)) {
    return(invisible(NULL))
  }
  if (isTRUE(getOption("flexyBayes.silence_sampler_arg_note", FALSE))) {
    return(invisible(NULL))
  }
  passed <- c(
    if (!is.null(seed)) "`seed`",
    if (!is.null(control)) "`control`"
  )
  tail_note <- if (identical(engine, "INLA")) {
    paste0(
      " INLA's nested Laplace approximation draws no random numbers, so ",
      "`seed` cannot pin a fit. On a well-identified model two runs ",
      "agree because the optimiser does. On a weakly identified variance ",
      "component the landing point can move -- treat that the same way ",
      "as a failed sampler diagnostic, not as something a seed would fix."
    )
  } else {
    ""
  }
  message(
    "flexyBayes: ", paste(passed, collapse = " and "),
    if (length(passed) == 1L) " is" else " are",
    " a Stan sampler setting the ", engine, " backend does not use, so ",
    if (length(passed) == 1L) "it has" else "they have",
    " no effect on this fit.", tail_note
  )
  invisible(NULL)
}


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
  residual,
  family,
  link,
  data_name,
  aggregate = "auto",
  seed = NULL,
  control = NULL,
  na_action = NULL
) {
  aggregate <- .normalise_aggregate(aggregate)

  # Observation weights are parsed into an addition term by every ingest
  # adapter, and consumed by no active emitter -- brms and INLA generate
  # identical code with and without them. Refuse here, at the single
  # choke point every fitting and code-inspection route passes through,
  # rather than return the unweighted posterior under a weighted call.
  # The IR keeps the parsed weights: the request stays on the record,
  # and the parser tests continue to cover it.
  .refuse_unsupported_weights(fb, weights)

  # Structures that parse for the formula catalogue and that no active
  # engine emits -- the residual-side AR1 spelling above all -- are refused
  # by name here, in one vocabulary, rather than as an lgm_gate() print on
  # one engine and an untyped stop() on the other.
  .refuse_unrepresentable_structures(fb)

  # A pooled() prior names its term; one the model does not have is refused
  # here, before either engine is chosen.
  .fb_pooled_blocks(fb)

  # D14 / WP-D: a family = "binomial" response outside {0, 1} on this
  # (non-streaming) entry used to reach the engine unrefused on every
  # route -- an INLA subprocess death, an untyped brms condition, and an
  # ad hoc untyped stop() on the aggregated route (now removed; see
  # below). One typed refusal here, before any route is attempted.
  .refuse_non_binary_binomial_response(fb, data)

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
      known_matrices = known_matrices,
      weights = weights,
      n_samples = n_samples,
      warmup = warmup,
      chains = chains,
      seed = seed,
      control = control,
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd,
      verbose = verbose,
      mcmc_verbose = mcmc_verbose,
      the_call = the_call,
      fixed = fixed,
      random = random,
      residual = residual,
      family = family,
      link = link,
      data_name = data_name,
      na_action = na_action
    )
    if (!is.null(agg_fit)) {
      if (!is.null(agg_fit$extras) && !is.null(fb_dataset_meta_for_fit)) {
        agg_fit$extras$fb_dataset <- fb_dataset_meta_for_fit
      }
      return(agg_fit)
    }
  }

  # low_rank_smooth truncates the dense mgcv smooth basis, and flexyBayes
  # has no route for the truncated basis on either engine, so it refuses
  # unconditionally, including under `auto`.
  if (length(.collect_approx(fb$random_terms)) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "low_rank_smooth_unsupported",
      message = paste0(
        "flexyBayes cannot fit a smooth with the low_rank_smooth ",
        "approximation: it has no route for the approximation on either ",
        "engine. Drop the representation = ... argument and write the ",
        "smooth as spl(x), which fits on INLA."
      ),
      family_class = "flexybayes_low_rank_smooth_unsupported_refusal"
    ))
  }

  # Stan passthrough via brms. Opt-in only -- backend = "auto" never
  # routes to Stan because brms's first-call compile latency (30-60s)
  # would silently break the auto-dispatch transparent-fast-path
  # promise. Users wanting a Stan fit set backend = "brms"
  # explicitly on fb_brms(); flexybayes() (the asreml entry) does
  # not advertise the value because asreml structured-covariance
  # terms (fa, us, ar1) do not translate losslessly to brms.
  if (identical(backend, "brms")) {
    # Capability gate. With flexybayes() (asreml entry) now able to
    # select backend = "brms", a term brms/Stan has no lossless
    # translation for can reach this branch -- an asreml
    # structured-covariance term (vm/ped/fa/us/ar1) or a structured
    # residual (ar1 / separable ar1_spatial). Refuse structurally before
    # the emit rather than fail mid-translation, or worse, emit a model
    # that silently omits the term.
    #
    # This refuses on ANY capability failure and re-raises the
    # predicate's own reason_code, rather than matching one code by
    # name. Enumerating codes here made the gate fail open: a capability
    # the predicate learned to detect later was silently emitted anyway,
    # which is how a requested AR1xAR1 residual reached Stan as an
    # intercept-only iid Gaussian. The low_rank case is handled upstream
    # by the .collect_approx() check.
    cap <- .backend_can_fit("brms", fb)
    if (!isTRUE(cap$ok)) {
      entry <- .lookup_refusal(cap$reason_code)
      stop(.fb_refusal_condition(
        reason_code = cap$reason_code,
        message = if (!is.null(entry)) {
          entry$message_template
        } else {
          paste0(
            "flexyBayes cannot fit this model on brms (",
            cap$reason_code, ")."
          )
        },
        family_class = paste0("flexybayes_", cap$reason_code)
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
      seed = seed,
      control = control,
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd,
      verbose = verbose,
      mcmc_verbose = mcmc_verbose,
      return_code = return_code,
      the_call = the_call,
      fixed = fixed,
      random = random,
      residual = residual,
      family = family,
      link = link,
      data_name = data_name,
      na_action = na_action,
      requested_backend = backend,
      requested_aggregate = aggregate
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
        reason = .fb_pooled_reason(paste0(
          "user requested brms explicitly ",
          "(Stan passthrough; lgm_gate() not ",
          "consulted on this path)."
        ), fb),
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
  if (backend %in% c("inla", "auto")) {
    # Gate receives preflight so the memory_feasibility_inla rule
    # can fire on INLA-specific memory infeasibility. preflight is
    # NULL below the 1e5-row threshold; the rule then trivially
    # passes per its v0.3.5 backward-compatible default.
    gated <- lgm_gate(fb, preflight = preflight_result)

    if (is_lgm_refusal(gated)) {
      if (backend == "inla") {
        # The gate's own vocabulary -- rule ids and a printed check-list --
        # is the useful part and is kept verbatim as the message. What it
        # lacked was a class: an explicit backend = "inla" request that the
        # gate refused raised a bare simpleError, so a caller could catch
        # every other flexyBayes refusal by name and not this one. The
        # check-list object rides on the condition as `lgm_refusal`, and
        # the failing rule supplies a second class so a single rule can be
        # caught on its own.
        .stop_inla_gate_refused(gated)
      }
      primary <- gated$failures[[1L]]
      # Distinguish memory vs structural refusal for routing trace.
      is_memory <- identical(primary$rule_id, "memory_feasibility_inla")
      gate_outcome <- if (is_memory) "refuse_memory" else "refuse_structural"
      # Resolve the fallback engine through the single routing brain
      # (.resolve_routing), passing fb so a multi-stratum designed
      # experiment routes to brms -- the faithful full-HMC path. Passing
      # fb also makes the plan (fb_plan) and the fit agree on the chosen
      # backend.
      routing <- .resolve_routing(
        user_request = "auto",
        gate_outcome = gate_outcome,
        preflight_outcome = NA_character_,
        inla_installed = inla_installed,
        fb = fb
      )
      fallback_backend <- routing$chosen_backend
      is_multistratum <- !is_memory && .has_interaction_random_term(fb)
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
        if (identical(fallback_backend, "brms")) {
          message(
            "backend = \"auto\": the INLA route refused this model (",
            primary$rule_id,
            "), so it is fitted on brms. Pass backend = \"brms\" to ",
            "silence this note."
          )
        } else if (is_multistratum) {
          message(
            "backend = \"auto\": the INLA route refused this model (",
            primary$rule_id,
            "), and brms is not installed or flexyBayes cannot fit the ",
            "model on it, so neither engine fits it. Install brms, or ",
            "change the term the refusal names."
          )
        } else {
          message(
            "backend = \"auto\": the INLA route refused this model (",
            primary$rule_id,
            ": ",
            primary$reason,
            "), and flexyBayes cannot fit it on brms either. Pass ",
            "backend = \"inla\" to see the refusal as an error."
          )
        }
        .emit_state_set("auto_fallback_note", TRUE)
      }
      decision <- .build_routing_decision(
        backend = fallback_backend,
        path = if (is_memory) {
          "auto_lgm_refuse_memory"
        } else if (identical(fallback_backend, "brms")) {
          "auto_multistratum_to_brms"
        } else {
          "auto_lgm_refuse"
        },
        gate_checks = gated$failures,
        reason = paste0(
          "lgm_gate() refused (",
          primary$rule_id,
          "); auto routed to ",
          fallback_backend,
          "."
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
          # The gate accepted, but the engine it accepted for is not on
          # this host. Resolve through the same brain the structural
          # refusal path uses -- brms when it can represent the model,
          # otherwise no active route.
          fallback_backend <- .resolve_auto_fallback_backend(gated)
          if (
            !isTRUE(getOption(
              "flexyBayes.silence_auto_inla_missing_note",
              FALSE
            )) &&
              !.emit_state_get("auto_inla_missing_note")
          ) {
            message(
              "backend = \"auto\": INLA is not installed; ",
              if (identical(fallback_backend, "brms")) {
                "routing to brms. "
              } else {
                paste0(
                  "flexyBayes cannot fit this model on brms either. "
                )
              },
              "Install INLA from ",
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
            inla_installed = FALSE
          )
          decision <- .build_routing_decision(
            backend = fallback_backend,
            path = "auto_inla_unavailable",
            gate_checks = gated$capabilities,
            reason = paste0(
              "lgm_gate() accepted; INLA not installed; auto routed to ",
              if (is.na(fallback_backend)) "no active backend" else
                fallback_backend,
              "."
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
        # accepted structurally (e.g. ar1_spatial_requires_complete_grid)
        # resolves through .resolve_auto_fallback_backend() to brms, or an
        # auto_no_active_route refusal -- this path must never silently
        # substitute an unrepresented engine.
        inla_fit <- tryCatch(
          .backend_emit_fn("inla")(
            fb = gated,
            data = data,
            known_matrices = known_matrices,
            weights = weights,
            n_samples = n_samples,
            warmup = warmup,
            chains = chains,
            seed = seed,
            control = control,
            prior_fixed_sd = prior_fixed_sd,
            prior_vc_sd = prior_vc_sd,
            verbose = verbose,
            return_code = return_code,
            the_call = the_call,
            fixed = fixed,
            random = random,
            residual = residual,
            family = family,
            link = link,
            data_name = data_name,
            na_action = na_action,
            requested_backend = backend,
            requested_aggregate = aggregate
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
          # Fallback: numerical/runtime INLA failure on the auto path.
          # Resolve through the same brain the structural-refusal path
          # uses -- brms if it can represent the model, otherwise no
          # active backend can faithfully fit it and auto refuses (A3),
          # never a silent substitution of an unrepresented engine.
          fallback_backend <- .resolve_auto_fallback_backend(gated)
          if (
            !isTRUE(getOption(
              "flexyBayes.silence_auto_fallback_note",
              FALSE
            )) &&
              !.emit_state_get("auto_inla_numerical_fallback_note")
          ) {
            if (identical(fallback_backend, "brms")) {
              message(
                "backend = \"auto\": INLA failed numerically (",
                conditionMessage(inla_fit$cond),
                "); falling back to brms. Pass backend = \"inla\" ",
                "to surface the error as such, or backend = ",
                "\"brms\" to silence."
              )
            } else {
              message(
                "backend = \"auto\": INLA failed numerically (",
                conditionMessage(inla_fit$cond),
                "), and flexyBayes cannot fit this model on brms either. ",
                "Pass backend = \"inla\" to see the error as such."
              )
            }
            .emit_state_set("auto_inla_numerical_fallback_note", TRUE)
          }
          decision <- .build_routing_decision(
            backend = fallback_backend,
            path = if (identical(fallback_backend, "brms")) {
              "auto_inla_numerical_fallback_to_brms"
            } else {
              "auto_inla_numerical_fallback_no_route"
            },
            gate_checks = gated$capabilities,
            reason = paste0(
              "lgm_gate() accepted but INLA failed numerically; auto ",
              "routed to ",
              if (identical(fallback_backend, "brms")) {
                "brms"
              } else {
                "no active backend"
              },
              "."
            ),
            preflight_summary = preflight_result,
            representation_plan = .representation_plan_from_preflight(
              preflight_result
            ),
            rejected_routes = list(list(
              backend = "inla",
              reason = "inla_numerical_failure"
            ))
          )
          # fall through to the resolved backend emit below (brms, or an
          # auto_no_active_route refusal when brms cannot represent the
          # model either).
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
                inla_installed = TRUE
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
              reason = if (is.null(fit$extras$inla_integration)) {
                "lgm_gate() accepted; INLA dispatch."
              } else {
                paste0(
                  "lgm_gate() accepted; INLA dispatch by empirical Bayes ",
                  "(int.strategy = \"eb\"), past the hyperparameter budget."
                )
              },
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
  }
  # By construction `backend` is "inla", "auto", or "brms" (match.arg);
  # the "brms" branch above always returns before this point, so the
  # `if` above is the only place `decision` is assigned. `decision`
  # stays NULL only if that reachability invariant is ever broken, and
  # the NULL -> NA fallback below still refuses cleanly (auto_no_active_route)
  # rather than dereferencing a NULL `decision`.

  # The auto path resolves the fallback engine into decision$backend (inla,
  # or brms for a model brms can fit); when no active backend can
  # faithfully fit the model it leaves it NA, so refuse (A3), never
  # silently substitute.
  fit_backend <- decision$backend
  if (is.null(fit_backend) || is.na(fit_backend)) {
    stop(.fb_refusal_condition(
      reason_code = "auto_no_active_route",
      message = paste0(
        "backend = \"auto\": flexyBayes cannot fit this model on INLA or ",
        "on brms. Pass backend = \"inla\" or backend = \"brms\" to see ",
        "which term each route refuses, and change that term."
      ),
      family_class = "flexybayes_auto_no_route_refusal",
      backend = "auto"
    ))
  }

  # Sampler settings the call does not carry are omitted here rather than
  # forwarded as NULL. Every route reaching this point is a fallback --
  # the gate refused, INLA is absent, or INLA failed numerically -- and a
  # re-fit of an INLA fit arrives with `n_samples`, `warmup` and `chains`
  # all NULL, because the nested Laplace approximation runs none of them
  # and the INLA record says so. brms builds `iter` by adding two of them
  # and refuses the result outright ("Cannot coerce 'iter' to a single
  # numeric value"), so the user saw a brms argument error where an INLA
  # failure had happened. An absent argument means "use the default"
  # everywhere else in this package, and dropping it here is how the
  # receiving emit gets to apply its own.
  sampler_args <- list(
    n_samples = n_samples,
    warmup = warmup,
    chains = chains
  )
  sampler_args <- sampler_args[
    !vapply(sampler_args, is.null, logical(1L))
  ]

  # Emit entry-point resolved through the registry `engine` field
  # rather than a hard-coded symbol.
  #
  # `quote = TRUE` is load-bearing, not decoration. do.call() otherwise
  # re-evaluates every argument it is handed, and four of these are
  # language objects: `the_call` is a call, so evaluating it would RUN
  # the user's flexybayes() call a second time, and `fixed` / `random` /
  # `residual` are formulas, which would come out re-created against this
  # frame rather than the caller's -- a silent change of the environment
  # model.frame() later resolves the data in.
  fit <- do.call(
    .backend_emit_fn(fit_backend),
    c(
      list(
        fb = fb,
        data = data,
        known_matrices = known_matrices,
        weights = weights,
        seed = seed,
        control = control,
        prior_fixed_sd = prior_fixed_sd,
        prior_vc_sd = prior_vc_sd,
        verbose = verbose,
        mcmc_verbose = mcmc_verbose,
        return_code = return_code,
        the_call = the_call,
        fixed = fixed,
        random = random,
        residual = residual,
        family = family,
        link = link,
        data_name = data_name,
        na_action = na_action,
        requested_backend = backend,
        requested_aggregate = aggregate
      ),
      sampler_args
    ),
    quote = TRUE
  )

  if (!isTRUE(return_code) && !is.null(fit$extras)) {
    if (!is.null(decision)) {
      if (identical(decision$backend, "brms")) {
        decision$reason <- .fb_pooled_reason(decision$reason, fb)
      }
      fit$extras$backend_decision <- decision
    }
    fit$extras$fb_terms <- fb
    if (!is.null(fb_dataset_meta_for_fit)) {
      fit$extras$fb_dataset <- fb_dataset_meta_for_fit
    }
  }
  # Three-tier exactness vocabulary (v0.4.0). A fit carrying
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
  .fb_refuse_aggregate_argument(x)
}


# Raise the typed refusal for an explicit backend = "inla" request that
# lgm_gate() refused. Shared by the per-row dispatch below and by
# .maybe_aggregate_gaussian(): a non-canonical link on the aggregated
# route (reason_codes "non_canonical_link" / "non_identity_link" from
# .fb_aggregation_plan()) is the identical fact lgm_gate()'s
# family_allowlist rule already refuses for gaussian, binomial and
# poisson -- the aggregation canon link (identity / logit / log) IS each
# family's INLA default link -- so the two sites read the same gate the
# same way rather than growing a second, unregistered message for one
# fact.
.stop_inla_gate_refused <- function(gated) {
  msg <- paste(utils::capture.output(print(gated)), collapse = "\n")
  primary_rule <- gated$primary_rule %||% NA_character_
  primary_code <- gated$failures[[1L]]$reason_code
  stop(.fb_refusal_condition(
    reason_code = primary_code %||% "inla_gate_refused",
    message = msg,
    family_class = c(
      if (!is.null(primary_code)) {
        "flexybayes_refusal_inla_gate_refused"
      },
      if (!is.na(primary_rule)) {
        paste0("flexybayes_lgm_", primary_rule)
      }
    ) %||% character(0),
    lgm_refusal = gated,
    rule_id = primary_rule
  ))
}

# Aggregated-gaussian gate. Returns either a
# fitted <flexybayes_aggregated> object (gate accepted; aggregated
# emit ran) or NULL (gate declined; caller falls through to per-row).
# Refuses (raises a typed condition) when aggregate = TRUE but the
# plan declares ineligibility, or when aggregate = TRUE on a backend
# that does not support the aggregated path (brms, Stan).
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
  residual,
  family,
  link,
  data_name,
  known_matrices = NULL,
  weights = NULL,
  seed = NULL,
  control = NULL,
  na_action = NULL
) {
  # A missing response has no sufficient statistic to contribute, and the
  # aggregated emit compresses rows into per-cell means and counts without
  # accounting for that: the cell mean comes back NA and a downstream test
  # fails with "missing value where TRUE/FALSE needed". The failure is
  # reachable from an ordinary call -- a randomised complete block design with
  # a lost plot, fitted with the default aggregate = "auto".
  #
  # Fall through to the per-row path, which represents a missing response as
  # the latent quantity `na_action` says it is. Teaching the compression to
  # skip NA rows would be the wrong repair: it would silently make the
  # aggregated path mean complete-case deletion while the per-row path means
  # augmentation, so the same call would answer two different questions
  # depending on a compression ratio the user never set.
  resp_vals <- .fb_response_values(fb, data)
  if (!is.null(resp_vals) && anyNA(resp_vals)) {
    if (isTRUE(aggregate)) {
      # Carries `flexybayes_aggregate_refusal` beneath the per-code class
      # for the same reason the two refusals further down do: a caller
      # pattern-matching the family catches every aggregated-route
      # refusal rather than only some of them.
      stop(.fb_refusal_condition(
        reason_code = "aggregation_response_incomplete",
        message = paste0(
          "`aggregate = TRUE` refused: the response has missing values, and ",
          "the aggregated path compresses rows into per-cell sufficient ",
          "statistics that are not defined for a missing response. Pass ",
          "aggregate = FALSE (or leave it at \"auto\") for the per-row path, ",
          "which carries a missing response as a latent quantity."
        ),
        family_class = "flexybayes_aggregate_refusal"
      ))
    }
    return(NULL)
  }

  # C6: observation weights are lowered for the per-row Gaussian route
  # only. Folding a per-observation weight into the aggregated route's
  # sufficient statistics (cell sums and counts) would change what
  # those statistics mean -- a weighted cell mean is not the same
  # quantity the aggregated Gaussian emit's closed form assumes -- so
  # weights are not carried through the compression. Same isTRUE(aggregate)
  # -refuses / "auto" falls through to the per-row path -- which DOES
  # honour weights -- shape as the missing-response guard above.
  w <- .fb_ir_weights(fb) %||% weights
  if (.fb_weights_nonconstant(w)) {
    if (isTRUE(aggregate)) {
      stop(.fb_refusal_condition(
        reason_code = "weights_not_aggregatable",
        message = paste0(
          "`aggregate = TRUE` refused: observation weights are lowered ",
          "for the per-row Gaussian route only. The aggregated route's ",
          "closed form compresses rows into per-cell sums and counts, and ",
          "does not fold a per-observation weight into that compression -- ",
          "the weighted cell statistic is a different quantity from the ",
          "one the aggregated Gaussian emit's sufficient statistics ",
          "assume. Pass aggregate = FALSE (or leave it at \"auto\") for ",
          "the per-row weighted fit."
        ),
        family_class = "flexybayes_aggregate_refusal"
      ))
    }
    return(NULL)
  }

  # Aggregation only defined for backends with a flexyBayes-side
  # aggregated emit path. INLA supports it; brms / Stan have no
  # aggregated path in scope.
  agg_capable_backend <- backend %in% c("inla", "auto")
  if (!agg_capable_backend) {
    if (isTRUE(aggregate)) {
      stop(.fb_refusal_condition(
        reason_code = "aggregation_route_unavailable",
        message = paste0(
          "`aggregate = TRUE` is not supported on backend = \"",
          backend,
          "\" (the aggregated path is wired for inla only). Pass ",
          "aggregate = FALSE or switch backend."
        ),
        family_class = "flexybayes_aggregate_refusal",
        backend = backend
      ))
    }
    return(NULL)
  }

  # Run the model-level plan against the fb_dataset metadata. Reads
  # only IR + dictionaries; no design matrix is materialised.
  fb_dataset_meta <- .fb_dataset(data)
  plan <- .fb_aggregation_plan(fb, fb_dataset_meta)

  # Resolve the effective backend for the aggregated emit. INLA is the
  # only backend with an aggregated path.
  eff_backend <- .agg_resolve_backend(backend, fb)
  if (is.na(eff_backend)) {
    # No active aggregated backend (INLA refused / unavailable).
    # aggregate = TRUE cannot be honoured -> refuse; otherwise fall
    # through to the per-row path (main dispatch).
    if (isTRUE(aggregate)) {
      stop(.fb_refusal_condition(
        reason_code = "aggregation_route_unavailable",
        message = paste0(
          "`aggregate = TRUE` refused: no active aggregated backend (INLA ",
          "refused or unavailable). Pass aggregate = FALSE for the ",
          "per-row path."
        ),
        family_class = "flexybayes_aggregate_refusal",
        backend = backend
      ))
    }
    return(NULL)
  }

  # Heterogeneous residual (at_units) on the aggregated path: the aggregated
  # INLA emit has no sectioned residual, and there is no aggregated route on
  # brms. Refuse aggregate = TRUE, else fall through to the per-row path.
  uses_at_units <- length(fb$residual_terms) > 0L &&
    any(vapply(
      fb$residual_terms,
      function(t) identical(t$type %||% "", "at_units"),
      logical(1L)
    ))
  if (uses_at_units) {
    if (isTRUE(aggregate)) {
      # Carries the same condition class as the plan-ineligibility refusal
      # a few lines below, so a caller pattern-matching on
      # `flexybayes_aggregate_refusal` catches every aggregated-route
      # refusal rather than only some of them.
      stop(.fb_refusal_condition(
        reason_code = "aggregation_plan_ineligible",
        message = paste0(
          "`aggregate = TRUE` refused: flexyBayes's aggregated route has ",
          "no sectioned residual. Pass aggregate = FALSE to fit the model ",
          "row by row."
        ),
        family_class = "flexybayes_aggregate_refusal",
        call = the_call,
        reason_codes = "structured_residual_not_aggregatable",
        plan = NULL
      ))
    }
    return(NULL)
  }

  if (!isTRUE(plan$eligible)) {
    if (isTRUE(aggregate)) {
      # A cell count past the integer limit carries its own typed refusal
      # on the plan; raise that one, not the generic ineligibility.
      if (!is.null(plan$overflow_refusal)) {
        stop(plan$overflow_refusal)
      }
      # The family / link scope check runs first inside
      # .fb_aggregation_plan() (before the term-shape loops), so a family
      # or link reason is always reason_codes[1] when either fired. Both
      # facts are already named refusals elsewhere in the package, and
      # this branch used to raise them as an unregistered condition
      # instead of reusing that name:
      #   - non_canonical_link / non_identity_link is exactly the fact
      #     lgm_gate()'s family_allowlist rule refuses for an explicit
      #     backend = "inla" request, for the three aggregatable
      #     families (the aggregation canon link -- identity, logit,
      #     log -- IS each one's INLA default), so it is raised the same
      #     way: inla_gate_refused (.stop_inla_gate_refused()).
      #   - non_aggregatable_family is exactly what aggregation_route_
      #     unavailable already names at the other aggregated-route and
      #     streamed-route refusal sites in this file and in
      #     .check_stream_family() (R/stream_aggregate.R): no aggregated
      #     backend resolves for this family.
      primary_reason <- plan$reason_codes[[1L]]
      if (primary_reason %in% c("non_canonical_link", "non_identity_link")) {
        gated <- lgm_gate(fb)
        if (is_lgm_refusal(gated)) {
          .stop_inla_gate_refused(gated)
        }
        # lgm_gate() agreeing the pair is admitted on INLA cannot happen
        # here: the two reason codes above and the family_allowlist rule
        # read the same family/link pair against the same default-link
        # fact. Fall through to the generic refusal below rather than
        # silently letting an unreachable disagreement through unrefused.
      } else if (identical(primary_reason, "non_aggregatable_family")) {
        stop(.fb_refusal_condition(
          reason_code = "aggregation_route_unavailable",
          message = paste0(
            "`aggregate = TRUE` refused: the aggregated route is exact ",
            "only for gaussian, binomial and poisson, and family = \"",
            fb$family,
            "\" has no aggregated emit.",
            .agg_refusal_note(plan$reason_codes),
            " Pass aggregate = FALSE for the per-row path."
          ),
          family_class = "flexybayes_aggregate_refusal"
        ))
      }
      stop(.fb_refusal_condition(
        reason_code = "aggregation_plan_ineligible",
        message = paste0(
          "`aggregate = TRUE` refused by .fb_aggregation_plan(): ",
          paste(plan$reason_codes, collapse = ", "),
          ".",
          .agg_refusal_note(plan$reason_codes),
          " Pass aggregate = FALSE for the per-row path."
        ),
        family_class = "flexybayes_aggregate_refusal",
        call = the_call,
        reason_codes = plan$reason_codes,
        plan = plan
      ))
    }
    return(NULL)
  }

  # A Gaussian model is not fitted on the aggregated representation. The
  # aggregated Gaussian route composes the within-cell sum of squares with
  # a fixed gamma prior on the residual precision, so a stated or default
  # residual prior would not be the prior fitted while prior_summary()
  # reported it. Until that composition carries the residual prior,
  # "auto" takes the per-row route and aggregate = TRUE refuses. Every
  # refusal above still fires first, because each names a more specific
  # reason. The count families carry no residual prior and aggregate below.
  if (identical(fb$family, "gaussian")) {
    if (isTRUE(aggregate)) {
      stop(.fb_refusal_condition(
        reason_code = "aggregation_residual_prior_not_composed",
        message = paste0(
          "`aggregate = TRUE` refused on a Gaussian model: the aggregated ",
          "Gaussian route cannot yet carry the residual prior, so the fit ",
          "would not run under the prior prior_summary() reports. Pass ",
          "aggregate = FALSE, or leave it at \"auto\", which fits a Gaussian ",
          "model per row. Binomial and Poisson models still aggregate."
        ),
        family_class = "flexybayes_aggregate_refusal"
      ))
    }
    return(NULL)
  }

  # The binomial and Poisson families use the family-generic streaming
  # aggregator over the in-memory data (one chunk) and the count emit.
  # Exact aggregation for the count families requires unit-trial
  # binomial (a 0/1 numeric response); a non-Bernoulli binomial needs
  # per-row trial counts the cell sums cannot recover. This used to be
  # checked here with an untyped stop() -- now dead code, because
  # .refuse_non_binary_binomial_response() (D14 / WP-D) already refused
  # unconditionally at the top of .dispatch_backend(), the only caller
  # of .maybe_aggregate_gaussian(), before this branch is ever reached.
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
    known_matrices = known_matrices,
    weights = weights,
    n_samples = n_samples,
    warmup = warmup,
    chains = chains,
    seed = seed,
    control = control,
    prior_fixed_sd = prior_fixed_sd,
    prior_vc_sd = prior_vc_sd,
    verbose = verbose,
    mcmc_verbose = mcmc_verbose,
    the_call = the_call,
    fixed = fixed,
    random = random,
    residual = residual,
    family = family,
    link = link,
    data_name = data_name,
    na_action = na_action,
    requested_backend = backend,
    requested_aggregate = aggregate
  )

  fit$exactness <- "aggregated_exact"
  if (!is.null(fit$extras)) {
    # The aggregated path's chosen backend was already
    # resolved by .agg_resolve_backend() above; rejected_routes is
    # empty (the aggregation gate itself bypasses the routing policy
    # table -- this is its own dispatch path). preflight is NULL here
    # because aggregation runs upstream of the >1e5 preflight
    # threshold in practice; if invoked under the metadata-only
    # path the slot stays NULL and reproducibility checks treat
    # it as legitimately absent.
    fit$extras$backend_decision <- .build_routing_decision(
      backend = eff_backend,
      path = "aggregated_count",
      gate_checks = NULL,
      # `K` is said twice on one fit and used to mean two things. This is
      # the PLANNER's cell count -- the product of the declared factor
      # level counts, i.e. the complete grid -- and it is what made the
      # route eligible. The print and summary compression lines report
      # the REALISED count, `nrow()` of the cell-key combinations the
      # data actually contains, which is smaller whenever a grid cell
      # carries no rows. Both were labelled plain `K`, so a user reading
      # 135 here and 133 there had nothing to tell them why. The estimate
      # says it is one.
      reason = sprintf(
        "aggregation plan eligible (N = %d, K = %d (estimated), %s)",
        plan$N,
        plan$K_est,
        sprintf("ratio = %.2f:1", plan$N / plan$K_est)
      ),
      preflight_summary = NULL,
      representation_plan = NULL,
      rejected_routes = list()
    )
    fit$extras$fb_terms <- fb
  }
  fit
}


# Choose the effective aggregated-emit backend. Only INLA has an
# aggregated path: an explicit inla passes through; auto consults
# lgm_gate -- INLA if LGM-compatible and installed, else NA (no active
# aggregated route).
.agg_resolve_backend <- function(backend, fb) {
  if (identical(backend, "inla")) {
    return("inla")
  }
  if (identical(backend, "auto")) {
    # Let a genuine lgm_gate() fault propagate (charter: never swallow an
    # error into a silent fallback). A refusal is a normal return value (an
    # `lgm_refusal` object), not an error; INLA is the only active aggregated
    # backend, so a refusal / missing INLA leaves NA (no active aggregated
    # route) and the caller falls through / refuses.
    gated <- lgm_gate(fb)
    if (!is_lgm_refusal(gated) && requireNamespace("INLA", quietly = TRUE)) {
      return("inla")
    }
    return(NA_character_)
  }
  # Unknown backend -- caller already filtered, defensive NA.
  NA_character_
}


# Whether dispatch runs the design-memory preflight on this IR. fb_plan()
# reads the same rule, so a plan refuses only where the fit would.
.preflight_runs_at_dispatch <- function(fb, threshold = 1e5L) {
  isTRUE(!is.null(fb$data_summary$n) && fb$data_summary$n >= threshold)
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
  if (!.preflight_runs_at_dispatch(fb, threshold)) {
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
    # refusal class still reads correctly before its headline is added.
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


# ---------------------------------------------------------------- #
# backend_decision() -- dispatch-trace reader                       #
# ---------------------------------------------------------------- #

#' Backend dispatch trace for a flexyBayes fit
#'
#' Returns the dispatch trace recorded at fit time: which backend was
#' selected, which gate checks ran, and why. Engine-agnostic -- works
#' on an `flexybayes_inla` or `flexybayes_brms` fit alike, reading
#' `fit$extras$backend_decision` as populated by `.build_routing_decision()`
#' during dispatch (see the top of this file).
#'
#' @param fit A `flexybayes` object.
#' @returns A list holding the recorded dispatch trace. The first four
#'   components are present on every fit. The four routing-trace fields
#'   are present on fits from v0.3.6 onwards and `NULL` on earlier ones,
#'   for backward compatibility.
#'
#'   \describe{
#'     \item{`backend`}{Character; one of `"inla"`, `"brms"`.}
#'     \item{`path`}{Character; the dispatch path token.}
#'     \item{`gate_checks`}{List or NULL; the `lgm_gate()` check
#'       trail (failures on refusal; capabilities on accept).}
#'     \item{`reason`}{Character; the dispatch-decision rationale.}
#'     \item{`preflight_summary`}{An `<fb_preflight>` object or
#'       NULL. Populated when `.fb_preflight()` ran (>1e5-row
#'       path); NULL on the small-data fast path.}
#'     \item{`representation_plan`}{Named list of slim per-term
#'       entries `(term_id, representation_class, justification)`
#'       derived from `preflight_summary`; NULL when no preflight.}
#'     \item{`rejected_routes`}{List of `(backend, reason)` pairs
#'       for the backends considered but not chosen. Empty for
#'       explicit user requests (the routing policy is bypassed
#'       when the backend is named directly).}
#'     \item{`routing_policy_version`}{Character; e.g.
#'       `"stage5a_v1"`. The audit-anchor for reproducibility -- a
#'       policy change bumps this string.}
#'   }
#' @examples
#' # backend_decision() reads a slot recorded on an already-fitted
#' # object, so (unlike fb_plan()) an engine has to run first. Wrapped
#' # in \donttest{} and guarded on INLA so this does not fire on a
#' # machine without it.
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
#'   fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
#'                      verbose = FALSE)
#'   bd <- backend_decision(fit)
#'   bd$backend   # "inla"
#'   bd$path      # the dispatch path token
#'   bd$reason    # why this backend was chosen
#' }
#' }
#' @export
backend_decision <- function(fit) {
  if (!inherits(fit, "flexybayes") && !inherits(fit, "flexybayes_inla")) {
    stop(
      "`fit` must be a `flexybayes` (or `flexybayes_inla`) object.",
      call. = FALSE
    )
  }
  decision <- fit$extras$backend_decision
  if (is.null(decision)) {
    # Older fit (predates dispatch-trace recording, v0.3.6): synthesise
    # a minimal legacy trace so downstream tooling never sees NULL.
    # Any fit reaching this branch is either flexybayes_inla or (by
    # elimination, since those are the only two active engine classes)
    # flexybayes_brms.
    return(list(
      backend = if (inherits(fit, "flexybayes_inla")) "inla" else "brms",
      path = "legacy_no_recorded_decision",
      gate_checks = NULL,
      reason = paste0(
        "fit predates dispatch-trace recording; ",
        "backend inferred from class."
      ),
      preflight_summary = NULL,
      representation_plan = NULL,
      rejected_routes = list(),
      routing_policy_version = NA_character_
    ))
  }
  # v0.3.6+: surface the four new fields with
  # NULL defaults so earlier fits (those that recorded the
  # 4-field trace before v0.3.6 landed) read out with the
  # same shape downstream tooling sees on v0.3.6+ fits. Construct
  # a fresh canonical 8-field list with explicit NULL values --
  # `decision[["x"]] <- NULL` removes the slot rather than setting
  # it to NULL, which would silently drop the field on every read.
  rejected <- decision$rejected_routes
  if (is.null(rejected)) {
    rejected <- list()
  }
  version <- decision$routing_policy_version
  if (is.null(version)) {
    version <- NA_character_
  }
  list(
    backend = decision$backend,
    path = decision$path,
    gate_checks = decision$gate_checks,
    reason = decision$reason,
    preflight_summary = decision$preflight_summary,
    representation_plan = decision$representation_plan,
    rejected_routes = rejected,
    routing_policy_version = version
  )
}
