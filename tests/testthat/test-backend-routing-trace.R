# test-backend-routing-trace.R -- Stage 4 (ADR 0024).
#
# Covers the v0.3.6 routing-trace surface end-to-end:
#   - backend_decision(fit) 8-field shape with the four new fields
#     (preflight_summary, representation_plan, rejected_routes,
#     routing_policy_version).
#   - lgm_gate(fb, preflight) 11th rule (memory_feasibility_inla)
#     with the v0.3.5 backward-compatible preflight=NULL default.
#   - .routing_policy_table() invariants -- every (gate_outcome,
#     preflight_outcome, user_request) tuple maps to at most one row.
#   - Approximate-scheme entry-function refusal.
#   - Reproducibility of the trace across two identical calls.

# 0.9.3: a third native engine and its dormant sibling were withdrawn
# entirely (see NEWS.md). Neither is registered any more, so the
# routing-policy table carries no row naming either, and the auto
# path's rejected_routes no longer carries a dormant-sibling entry
# (ADR 0024 (i), below, is removed rather than ported for exactly this
# reason). The fixture default backend and every explicit-request
# subtest move to brms, the remaining full-HMC engine.

suppressPackageStartupMessages({
  library(testthat)
})

old_opts <- options(
  flexyBayes.silence_default_prior_note = TRUE,
  flexyBayes.silence_uniform_inla_approx = TRUE,
  flexyBayes.silence_auto_fallback_note = TRUE,
  flexyBayes.silence_auto_inla_missing_note = TRUE
)
on.exit(options(old_opts), add = TRUE)


# Test fixture -- small brms-backend fit on a 4-level RI model.
# Used by most subtests to keep wall-time bounded.
mk_routing_fit <- function(seed = 2026L, N = 60L, J = 4L, backend = "brms") {
  set.seed(seed)
  dat <- data.frame(
    x = rnorm(N),
    g = factor(sample(letters[seq_len(J)], N, replace = TRUE)),
    y = rnorm(N)
  )
  # Structural fixture (backend_decision() trace shape, not numeric
  # agreement) shared by most subtests in this file. The property
  # under test does not depend on sampler convergence, but a chains=1L
  # / 30-draw budget produces genuine R-hat (not just ESS) warnings on
  # this 4-level random intercept, which .muffle_ess_warnings() does
  # not (and should not) catch -- rather than widen the muffle pattern
  # to hide a real convergence signal, the budget is raised until the
  # fit actually converges (chains=2L, more draws, higher adapt_delta).
  # fb()'s own `seed` argument (distinct from the `set.seed(seed)`
  # call above, which only fixes the *data*) must also be pinned --
  # left at the fb() default of NULL, brms passes Stan `seed = NA`
  # (draw-your-own), so the sampler's convergence was observed to vary
  # run to run at identical data and settings. Pinning it makes the
  # already-verified-clean draw reproducible instead of leaving it to
  # chance.
  fit <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = dat,
    backend = backend,
    n_samples = 500L,
    warmup = 500L,
    chains = 2L,
    seed = seed,
    control = list(adapt_delta = 0.97),
    prior_fixed_sd = 10,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  list(fit = fit, dat = dat)
}


# ---------------------------------------------------------------- #
# ADR §6 (a) -- preflight_summary is non-NULL when preflight ran    #
# ---------------------------------------------------------------- #

test_that("ADR 0024 (a): backend_decision() carries preflight_summary slot (NULL on small N)", {
  testthat::skip_if_not_installed("brms")
  fx <- mk_routing_fit()
  bd <- backend_decision(fx$fit)
  expect_true("preflight_summary" %in% names(bd))
  # Small-N fits skip preflight by design (>1e5 threshold).
  expect_null(bd$preflight_summary)
})


# ---------------------------------------------------------------- #
# ADR §6 (b) -- representation_plan slot present                    #
# ---------------------------------------------------------------- #

test_that("ADR 0024 (b): backend_decision() carries representation_plan slot (NULL on small N)", {
  testthat::skip_if_not_installed("brms")
  fx <- mk_routing_fit()
  bd <- backend_decision(fx$fit)
  expect_true("representation_plan" %in% names(bd))
  expect_null(bd$representation_plan)
})


# ---------------------------------------------------------------- #
# ADR §6 (c) -- rejected_routes lists non-chosen backends           #
# ---------------------------------------------------------------- #

test_that("ADR 0024 (c): backend_decision()$rejected_routes lists non-chosen backends with reason codes", {
  testthat::skip_if_not_installed("INLA")
  # backend = "auto" exercises the rejected_routes population path.
  fx <- mk_routing_fit(backend = "auto")
  bd <- backend_decision(fx$fit)
  expect_true("rejected_routes" %in% names(bd))
  expect_type(bd$rejected_routes, "list")
  # With only two active engines, auto's rejected_routes is typically
  # empty (INLA accepts, or the fallback to brms is a route taken, not
  # rejected) -- the shape check below still holds for whichever rows,
  # if any, a future third engine would add.
  if (length(bd$rejected_routes) > 0L) {
    for (rr in bd$rejected_routes) {
      expect_true(all(c("backend", "reason") %in% names(rr)))
      expect_true(is.character(rr$backend))
      expect_true(is.character(rr$reason))
    }
  }
})


# ---------------------------------------------------------------- #
# ADR §6 (d) -- routing_policy_version == "stage5a_v1"              #
# ---------------------------------------------------------------- #

test_that("ADR 0024 (d): routing_policy_version is 'stage5a_v1' on fresh fits", {
  testthat::skip_if_not_installed("brms")
  fx <- mk_routing_fit()
  bd <- backend_decision(fx$fit)
  expect_identical(bd$routing_policy_version, "stage5a_v1")
})


# ---------------------------------------------------------------- #
# ADR §6 (e) -- memory-infeasibility refusal                        #
# ---------------------------------------------------------------- #

test_that("ADR 0024 (e): memory_feasibility_inla rule fires when INLA 3x ceiling exceeded", {
  # Direct rule-level test -- the dispatch-level memory-infeasibility
  # path requires a >1e5-row preflight; we test the gate rule in
  # isolation with a synthetic preflight result that crosses the
  # INLA 3x ceiling but not the hard ceiling.
  fb <- fb_from_brms(
    y ~ x + (1 | g),
    data = data.frame(
      y = rnorm(50),
      x = rnorm(50),
      g = factor(rep(letters[1:5], 10))
    )
  )
  # Synthetic preflight: indexed total = 500 MB, ceiling = 1 GB.
  # INLA 3x = 1.5 GB > 1 GB ceiling -> rule should refuse.
  synthetic_pf <- structure(
    list(
      per_term_estimate = list(
        x = list(
          design_memory_bytes = 1e7,
          representation_class = "indexed_continuous"
        ),
        g = list(
          design_memory_bytes = 5e8 - 1e7,
          representation_class = "indexed_factor"
        )
      ),
      total_estimate_bytes = 5e8, # 500 MB
      ceiling_bytes = 1e9, # 1 GB
      n_rows = 1e7,
      aggregation_plan = NULL,
      refusal = NULL
    ),
    class = c("fb_preflight", "list")
  )
  gated <- lgm_gate(fb, preflight = synthetic_pf)
  expect_s3_class(gated, "lgm_refusal")
  expect_identical(gated$primary_rule, "memory_feasibility_inla")
})

test_that("ADR 0024 (e+): structural-vs-memory refusal distinguishable from rejected_routes alone", {
  # The .inla_rejection_reason() helper distinguishes memory_infeasibility_inla
  # from structural_infeasibility_inla on the auto path's rejected_routes.
  mem_reason <- flexyBayes:::.inla_rejection_reason("refuse_memory", TRUE)
  str_reason <- flexyBayes:::.inla_rejection_reason("refuse_structural", TRUE)
  expect_identical(mem_reason, "memory_infeasibility_inla")
  expect_identical(str_reason, "structural_infeasibility_inla")
  expect_false(identical(mem_reason, str_reason))
})


# ---------------------------------------------------------------- #
# ADR §6 (f) -- structural-infeasibility refusal                    #
# ---------------------------------------------------------------- #

test_that("ADR 0024 (f): structural LGM-infeasible model surfaces structural reason on auto", {
  # vm() carries random_term_type_inla refusal (per .inla_random_term_type_allowlist).
  # We build a synthetic IR directly to bypass requiring a vm() ingest path.
  fb <- fb_from_brms(
    y ~ x + (1 | g),
    data = data.frame(
      y = rnorm(50),
      x = rnorm(50),
      g = factor(rep(letters[1:5], 10))
    )
  )
  # Inject an out-of-allowlist random-term type to trigger the
  # structural rule.
  fb$random_terms[[1L]]$type <- "vm"
  gated <- lgm_gate(fb)
  expect_s3_class(gated, "lgm_refusal")
  expect_identical(gated$primary_rule, "random_term_type_inla")
  # The auto rejection reason for INLA should be the structural code.
  rj <- flexyBayes:::.inla_rejection_reason("refuse_structural", TRUE)
  expect_identical(rj, "structural_infeasibility_inla")
})


# ---------------------------------------------------------------- #
# ADR §6 (g) -- explicit brms request bypasses policy                #
# ---------------------------------------------------------------- #

test_that("ADR 0024 (g): explicit backend='brms' bypasses policy (rejected_routes empty)", {
  testthat::skip_if_not_installed("brms")
  testthat::skip_on_cran()
  # Use review_code to avoid the brms compile latency; the trace is
  # built on the dispatch path either way.
  set.seed(2026L)
  dat <- data.frame(
    y = rnorm(20),
    x = rnorm(20),
    g = factor(rep(letters[1:4], 5))
  )
  rev <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = dat,
    backend = "brms",
    review_code = TRUE,
    n_samples = 30L,
    warmup = 20L,
    chains = 1L,
    prior_fixed_sd = 10,
    prior_vc_sd = 1,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  # review_code returns a review token rather than a fitted object;
  # the routing decision is on the proceed() path. Skip the
  # backend_decision() check here -- the policy-bypass semantics
  # are verified via the .resolve_routing() unit test below.
  expect_s3_class(rev, "flexybayes_review")
})

test_that("ADR 0024 (g+): .resolve_routing() returns empty rejected_routes for explicit user requests", {
  for (req in c("brms", "inla")) {
    out <- flexyBayes:::.resolve_routing(
      user_request = req,
      gate_outcome = "accept",
      preflight_outcome = "clear",
      inla_installed = TRUE
    )
    expect_length(out$rejected_routes, 0L)
  }
})


# ---------------------------------------------------------------- #
# Multi-stratum designs route to brms (the faithful full-HMC path). #
# INLA collapses these variance components to zero; brms recovers   #
# them (agridat::besag.met).                                        #
# ---------------------------------------------------------------- #

test_that(".has_interaction_random_term() detects nested / combo terms", {
  d <- data.frame(
    y = rnorm(24),
    g = factor(rep(1:4, 6)),
    e = factor(rep(1:3, 8)),
    b = factor(rep(1:2, 12))
  )
  fb_multi <- fb_from_asreml(y ~ 1, random = ~ g + g:e + e:b, data = d)
  fb_simple <- fb_from_asreml(y ~ 1, random = ~ g, data = d)
  expect_true(flexyBayes:::.has_interaction_random_term(fb_multi))
  expect_false(flexyBayes:::.has_interaction_random_term(fb_simple))
})

test_that("auto routes a structurally-refused multi-stratum design to brms", {
  skip_if_not_installed("brms")
  d <- data.frame(
    y = rnorm(24),
    g = factor(rep(1:4, 6)),
    e = factor(rep(1:3, 8)),
    b = factor(rep(1:2, 12))
  )
  fb_multi <- fb_from_asreml(y ~ 1, random = ~ g + g:e + e:b, data = d)
  out <- flexyBayes:::.resolve_routing(
    user_request = "auto",
    gate_outcome = "refuse_structural",
    preflight_outcome = NA_character_,
    inla_installed = TRUE,
    fb = fb_multi
  )
  expect_equal(out$chosen_backend, "brms")
  expect_equal(out$reason_code, "auto_multistratum_to_brms")
})

test_that("auto reroutes a generic structural refusal to brms", {
  skip_if_not_installed("brms")
  d <- data.frame(y = rnorm(24), g = factor(rep(1:4, 6)))
  fb_simple <- fb_from_asreml(y ~ 1, random = ~ g, data = d)
  out <- flexyBayes:::.resolve_routing(
    user_request = "auto",
    gate_outcome = "refuse_structural",
    preflight_outcome = NA_character_,
    inla_installed = TRUE,
    fb = fb_simple
  )
  # brms can fit the simple random-intercept model, so the auto
  # structural-refusal path reroutes to brms via the generic fallback
  # reason (not the interaction multi-stratum reason above, as there
  # is no interaction term here).
  expect_equal(out$chosen_backend, "brms")
  expect_equal(out$reason_code, "auto_fallback_to_brms")
})

test_that(".resolve_routing() without fb refuses (no active route)", {
  out <- flexyBayes:::.resolve_routing(
    user_request = "auto",
    gate_outcome = "refuse_structural",
    preflight_outcome = NA_character_,
    inla_installed = TRUE
  )
  # Without an fb the fallback engine cannot be resolved -- there is no
  # active route -> refuse (NA chosen backend).
  expect_true(is.na(out$chosen_backend))
  expect_equal(out$reason_code, "auto_no_active_route")
})


# ---------------------------------------------------------------- #
# ADR §6 (h) -- approximate-scheme refusal                          #
# ---------------------------------------------------------------- #

test_that("ADR 0024 (h): backend='inla_pardiso_approximate' refuses with approximation-registry pointer", {
  dat <- data.frame(y = rnorm(20), x = rnorm(20))
  err <- expect_error(
    flexybayes(y ~ x, data = dat, backend = "inla_pardiso_approximate"),
    "approximate_route_not_yet_registered"
  )
  expect_match(conditionMessage(err), "approximation registry")
})

test_that("ADR 0024 (h+): backend='variational_advi' also refuses (variational pattern)", {
  dat <- data.frame(y = rnorm(20), x = rnorm(20))
  err <- expect_error(
    flexybayes(y ~ x, data = dat, backend = "variational_advi"),
    "approximate_route_not_yet_registered"
  )
  expect_match(conditionMessage(err), "approximation registry")
})

test_that("ADR 0024 (h++): approximate-scheme refusal works on fb_brms too", {
  dat <- data.frame(y = rnorm(20), x = rnorm(20))
  err <- expect_error(
    fb(y ~ x, data = dat, backend = "inla_pardiso_approximate"),
    "approximate_route_not_yet_registered"
  )
  expect_match(conditionMessage(err), "approximation registry")
})


# ---------------------------------------------------------------- #
# ADR §6 (i) -- [removed] dormant sibling-engine slot in rejected_routes #
# ---------------------------------------------------------------- #
#
# Deleted (not rewritten): the dormant sibling engine this pinned was
# withdrawn entirely at 0.9.3 (see NEWS.md) -- it is not registered at
# all any more, so `rejected_routes` never carries an entry for it.
# There is no successor: with only two active engines and no
# provisioned-but-dormant third slot, auto's rejected_routes is empty
# on the accept path this fixture exercises (see (c) above).


# ---------------------------------------------------------------- #
# ADR §6 (j) -- reproducibility of the trace                        #
# ---------------------------------------------------------------- #

test_that("ADR 0024 (j): two identical calls produce identical backend_decision traces (modulo draws)", {
  testthat::skip_if_not_installed("brms")
  fx1 <- mk_routing_fit(seed = 7L)
  fx2 <- mk_routing_fit(seed = 7L)
  bd1 <- backend_decision(fx1$fit)
  bd2 <- backend_decision(fx2$fit)
  expect_identical(bd1$backend, bd2$backend)
  expect_identical(bd1$path, bd2$path)
  expect_identical(bd1$reason, bd2$reason)
  expect_identical(bd1$rejected_routes, bd2$rejected_routes)
  expect_identical(bd1$routing_policy_version, bd2$routing_policy_version)
  expect_identical(bd1$representation_plan, bd2$representation_plan)
  # preflight_summary may differ if RAM probe is stateful; skip.
})


# ---------------------------------------------------------------- #
# Additional coverage                                               #
# ---------------------------------------------------------------- #

test_that("lgm_gate() backward-compat: preflight=NULL produces 10-rule pass on LGM-feasible IR", {
  fb <- fb_from_brms(
    y ~ x + (1 | g),
    data = data.frame(
      y = rnorm(50),
      x = rnorm(50),
      g = factor(rep(letters[1:5], 10))
    )
  )
  gated <- lgm_gate(fb) # no preflight arg -- v0.3.5 calling convention
  expect_true(is_fb_terms(gated))
  expect_true("lgm_compatible" %in% gated$capabilities)
})

test_that("lgm_gate() with preflight: rule 11 trivially passes when INLA 3x within ceiling", {
  fb <- fb_from_brms(
    y ~ x + (1 | g),
    data = data.frame(
      y = rnorm(50),
      x = rnorm(50),
      g = factor(rep(letters[1:5], 10))
    )
  )
  small_pf <- structure(
    list(
      per_term_estimate = list(),
      total_estimate_bytes = 1e6, # 1 MB -- well within typical ceiling
      ceiling_bytes = 1e10, # 10 GB ceiling
      n_rows = 100L,
      aggregation_plan = NULL,
      refusal = NULL
    ),
    class = c("fb_preflight", "list")
  )
  gated <- lgm_gate(fb, preflight = small_pf)
  expect_true(is_fb_terms(gated))
  # The 11th rule's pass appears in capabilities only as a warning
  # entry (it's a pass without warning); confirm no lgm_refusal.
  expect_false(is_lgm_refusal(gated))
})

test_that("Legacy fit (no Stage 4 fields) accessor surfaces 8-field shape with NULL defaults", {
  testthat::skip_if_not_installed("brms")
  fx <- mk_routing_fit()
  # Simulate a v0.3.5 fit by stripping the new fields.
  legacy <- fx$fit
  legacy$extras$backend_decision <- list(
    backend = "brms",
    path = "explicit_brms",
    gate_checks = NULL,
    reason = "legacy 4-field shape"
  )
  bd <- backend_decision(legacy)
  expect_true(all(
    c(
      "backend",
      "path",
      "gate_checks",
      "reason",
      "preflight_summary",
      "representation_plan",
      "rejected_routes",
      "routing_policy_version"
    ) %in%
      names(bd)
  ))
  expect_null(bd$preflight_summary)
  expect_null(bd$representation_plan)
  expect_length(bd$rejected_routes, 0L)
  expect_identical(bd$routing_policy_version, NA_character_)
})

test_that(".routing_policy_table() invariants: every row carries valid columns", {
  table <- flexyBayes:::.routing_policy_table()
  expect_s3_class(table, "data.frame")
  expect_true(all(
    c(
      "user_request",
      "gate_outcome",
      "preflight_outcome",
      "inla_installed",
      "chosen_backend",
      "reason_code"
    ) %in%
      names(table)
  ))
  expect_true(nrow(table) >= 7L)
  # Every non-NA user_request value is in the expected vocabulary.
  reqs <- na.omit(unique(table$user_request))
  expect_true(all(reqs %in% c("brms", "inla", "auto")))
  # Every non-NA gate_outcome is in the expected vocabulary.
  outs <- na.omit(unique(table$gate_outcome))
  expect_true(all(
    outs %in% c("accept", "refuse", "refuse_structural", "refuse_memory")
  ))
  # Every non-NA reason_code is non-empty.
  expect_true(all(nzchar(na.omit(table$reason_code))))
})
