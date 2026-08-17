# Tests for the refusal-taxonomy scaffold (ADR 0030 C7; v0.3.8 scaffold-
# only posture). The registry is empty at v0.3.8 by design --- the
# ~30 existing refusal sites in dispatch.R / lgm_gate.R /
# structured_cov.R / emit_inla.R are migrated in one batch at v0.4.0
# alongside the user-facing fb_refusals() accessor.
#
# Five acceptance criteria per v038-plan-2026-05-25 section 5.2:
#
#   (a) .refusal_registry exists post-.onLoad() and is locked
#   (b) .register_refusal() with valid args succeeds (synthetic code)
#   (c) .register_refusal() with duplicate reason_code refuses
#   (d) user-side assign() into .refusal_registry raises
#   (e) the registry is empty at v0.3.8 (deliberate scaffold-only)

# ---------------------------------------------------------------- #
# (a) registry exists + is locked                                    #
# ---------------------------------------------------------------- #

test_that(".refusal_registry exists and is locked post-load", {
  reg <- flexyBayes:::.refusal_registry
  expect_true(is.environment(reg))
  expect_true(environmentIsLocked(reg))
})

# ---------------------------------------------------------------- #
# (b) .register_refusal() with valid args (via fresh environment)   #
# ---------------------------------------------------------------- #
#
# The post-load registry is locked, so we cannot test .register_refusal()
# directly against it. Instead we exercise the helper's validation logic
# on a fresh unlocked environment by temporarily swapping. This honours
# the scaffold-only posture (no real entries land in the live registry).

test_that(".register_refusal() validates args before assigning", {
  # Direct validation calls (no environment mutation needed).
  expect_error(
    flexyBayes:::.register_refusal(
      reason_code = NULL,
      description = "x",
      message_template = "x",
      registered_in_adr = "ADR 0030",
      since_version = "0.3.8"
    ),
    regexp = "reason_code"
  )
  expect_error(
    flexyBayes:::.register_refusal(
      reason_code = "x",
      description = c("a", "b"),
      message_template = "x",
      registered_in_adr = "ADR 0030",
      since_version = "0.3.8"
    ),
    regexp = "description"
  )
  expect_error(
    flexyBayes:::.register_refusal(
      reason_code = "x",
      description = "x",
      message_template = "x",
      registered_in_adr = "ADR 0030",
      since_version = ""
    ),
    regexp = "since_version"
  )
})

# ---------------------------------------------------------------- #
# (c) duplicate reason_code refuses (deferred to v0.4.0 migration) #
# ---------------------------------------------------------------- #
#
# Since the v0.3.8 registry is empty + locked, the duplicate-detection
# branch cannot be exercised against the live registry without test-only
# unlock. We assert the behaviour by reading the function body --- the
# error string mentions "already registered".

test_that(".register_refusal() body carries the duplicate-detection branch", {
  body_src <- paste(
    deparse(body(flexyBayes:::.register_refusal)),
    collapse = "\n"
  )
  expect_true(grepl("already registered", body_src, fixed = TRUE))
  expect_true(grepl("environmentIsLocked", body_src, fixed = TRUE))
})

# ---------------------------------------------------------------- #
# (d) user-side assign() raises on locked environment                #
# ---------------------------------------------------------------- #

test_that("user-side assign() into .refusal_registry raises", {
  reg <- flexyBayes:::.refusal_registry
  expect_error(
    assign("synthetic_test_code", list(), envir = reg),
    regexp = "locked|cannot"
  )
})

# ---------------------------------------------------------------- #
# (e) v0.3.10 first-migration entries (ADR 0025 Decisions 3 + 4)    #
# ---------------------------------------------------------------- #
#
# v0.3.8 shipped the registry scaffold empty (zero entries; deliberate
# scaffold-only posture). v0.3.10 lands the first three entries via
# .populate_refusal_registry_v0310(): the two blocks-structural codes
# (block_partition_incomplete + block_not_positive_definite) plus the
# upgraded approximate_route_not_yet_registered. The full ~28-site
# bulk migration follows at v0.4.0 Wave 1 Phase 1C using the same
# .register_refusal() shape.

test_that(".refusal_registry carries the v0.3.10 first-migration entries", {
  reg <- flexyBayes:::.refusal_registry
  entries <- ls(envir = reg, all.names = TRUE)
  expect_true("block_partition_incomplete" %in% entries)
  expect_true("block_not_positive_definite" %in% entries)
  expect_true("approximate_route_not_yet_registered" %in% entries)
})

# ---------------------------------------------------------------- #
# (e2) refusal vocabulary: complete 32-code post-remediation set    #
# ---------------------------------------------------------------- #
#
# The bulk migration (.populate_refusal_registry_v0400()) registers
# the remaining user-facing refusal codes. The routing-decision
# reasons and the internal aggregate-out-of-scope control-flow
# signals are deliberately excluded (see fb_refusals() docs).

test_that(".refusal_registry holds the complete refusal vocabulary", {
  reg <- flexyBayes:::.refusal_registry
  entries <- ls(envir = reg, all.names = TRUE)
  # 31 from the Phase 1C bulk migration + the family-support refusal,
  # the four parse-time spec refusals, the tensor-smooth refusal, and
  # the two entry-point argument guards added in the 2026-05-30 audit
  # remediation (= 39); plus the three fb_cov() carrier-construction
  # refusals added at v0.4.0 Wave 2 Phase 2A (ADR 0030 C3):
  # fb_cov_type_unknown, fb_cov_missing_matrix, cov_arg_not_fb_cov (= 42);
  # plus stan_cannot_represent_structured_cov from ADR 0031 Phase 2b (the
  # brms capability gate now reachable from flexybayes()) (= 43); plus the
  # two grammar-polymorphism guards on the universal entry from ADR 0031
  # Phase 3: grammar_brms_with_asreml_terms,
  # grammar_brms_known_matrices_unsupported (= 45). At v0.5.0 (ADR 0031
  # Phase 3.6) the deferral code grammar_greta_via_fb_deferred is REMOVED
  # (fb() now fits a native greta graph) and two genuine refusals are
  # added: native_greta_requires_greta_backend (a native graph on a
  # non-greta engine) and engine_pin_backend_conflict (a pin handed a
  # conflicting backend) (= 47). At v0.6.0.9000 the gretaR backend is
  # activated: six gretaR refusals are registered --
  # gretaR_not_installed, gretaR_below_version_floor,
  # gretaR_family_unsupported, gretaR_random_group_not_in_data,
  # gretaR_random_term_type_unsupported,
  # gretaR_cannot_represent_structured_cov (= 53). At v0.7.0 the
  # data-aware factor-analytic rank upper bound is added:
  # fa_rank_exceeds_dim (fa(x, k) refused for k >= n_outer), the
  # identifiability complement of the data-free fa_rank_invalid floor
  # (= 54). At 0.9.0 the design-fidelity refusal
  # dsum_structured_inner_unsupported is added (a structured dsum() inner
  # -- ar1/us/... -- is refused rather than silently reduced to a
  # per-region heteroscedastic variance) (= 55).
  # At 0.9.1 the greta / gretaR quarantine adds three refusal codes
  # (backend_quarantined, native_greta_fit_quarantined, auto_no_active_route)
  # (= 58), and WP16 adds ar1_spatial_requires_complete_grid (= 59).
  # The residual-structure fidelity fix adds two: brms has no
  # residual-covariance lowering, so a structured residual reaching it must
  # refuse rather than emit a model without the term
  # (stan_cannot_represent_ar1_residual,
  # stan_cannot_represent_structured_residual) (= 61), and weights, parsed
  # but consumed by no emitter, are refused rather than silently returning
  # the unweighted posterior (weights_not_supported) (= 62). The
  # missing-response augmentation layer adds three: na_action = "fail"
  # with a missing response (missing_response_refused), a missing
  # predictor, which the device does not cover
  # (missing_covariate_not_supported), and a design cell absent from the
  # data whose model variables were never recorded
  # (augment_cell_not_determinable) (= 65). brms drops NA-response rows and
  # its mi() addition term is Gaussian-only, so a missing response on a
  # non-Gaussian family refuses rather than being silently deleted
  # (brms_cannot_augment_nongaussian) (= 66). at(f, level):g conditions the
  # random effect on selected levels, which no active emitter represents and
  # which is a different model from the heterogeneous variance diag(f):g that
  # shares its spelling (at_level_conditioning_unsupported) (= 67).
  # The 0.9.0 fidelity-of-names pass adds nine. The separable AR1 field is
  # respelled: the residual form names ASReml's three-parameter nugget-free
  # process and refuses (ar1_residual_not_representable), and the random-side
  # field is INLA-only (stan_cannot_represent_ar1_field) (= 72). The parser
  # vocabulary closes, so an unrecognised call is refused rather than read as
  # a variable of that name (asreml_function_not_recognised), and the four
  # parse-for-catalogue structures with no emit path are refused by name
  # (ar2_not_representable, str_not_representable, fa_not_representable,
  # interaction_not_representable) (= 77). The brms formula reconstruction's
  # untyped stop becomes a named refusal (brms_cannot_represent_term) (= 78),
  # and a numeric variable inside a random interaction is refused rather than
  # resolved to one of its two readings
  # (numeric_variable_in_random_interaction) (= 79). triangulate() gains a
  # comparability gate, and two fits of different models are refused rather
  # than compared (triangulate_incomparable_fits) (= 80). Aligning the class
  # graph so INLA fits inherit `flexybayes` makes five parent methods
  # reachable on them, and the four that cannot answer from an INLA fit's
  # slots refuse by name rather than erroring obscurely or returning an
  # empty result (fit_lacks_posterior_draws,
  # conditional_loglik_not_available, model_matrix_not_recoverable,
  # update_call_not_reconstructable) (= 84). An explicit backend = "inla"
  # request that the gate refuses stops with a class rather than a bare
  # simpleError, which was the last untyped refusal vocabulary in the
  # package (inla_gate_refused) (= 85). fb_met_summary() abstained with a
  # bare simpleError on both active engines and is now typed
  # (met_summary_not_available) (= 86).
  # A single-file or in-memory streaming source past 2^31 - 1 rows refuses
  # rather than recording NA (row_count_exceeds_integer) (= 87).
  # A `family` argument that is neither a name nor a usable stats family
  # object -- or a family object contradicting an explicit `link =` --
  # refuses rather than reaching base R's scalar `if` and erroring with
  # "the condition has length > 1" (family_argument_not_recognised)
  # (= 88).
  expect_equal(length(entries), 88L)
  expect_true("family_argument_not_recognised" %in% entries)
  expect_true("row_count_exceeds_integer" %in% entries)
  expect_true("met_summary_not_available" %in% entries)
  expect_true("inla_gate_refused" %in% entries)
  expect_true("fit_lacks_posterior_draws" %in% entries)
  expect_true("conditional_loglik_not_available" %in% entries)
  expect_true("model_matrix_not_recoverable" %in% entries)
  expect_true("update_call_not_reconstructable" %in% entries)
  expect_true("triangulate_incomparable_fits" %in% entries)
  expect_true("at_level_conditioning_unsupported" %in% entries)
  expect_true("ar1_residual_not_representable" %in% entries)
  expect_true("stan_cannot_represent_ar1_field" %in% entries)
  expect_true("asreml_function_not_recognised" %in% entries)
  expect_true("brms_cannot_represent_term" %in% entries)
  expect_true("numeric_variable_in_random_interaction" %in% entries)
  expect_true("brms_cannot_augment_nongaussian" %in% entries)
  expect_true("missing_response_refused" %in% entries)
  expect_true("missing_covariate_not_supported" %in% entries)
  expect_true("augment_cell_not_determinable" %in% entries)
  expect_true("stan_cannot_represent_structured_cov" %in% entries)
  expect_true("stan_cannot_represent_ar1_residual" %in% entries)
  expect_true("stan_cannot_represent_structured_residual" %in% entries)
  expect_true("weights_not_supported" %in% entries)
  expect_true("grammar_brms_with_asreml_terms" %in% entries)
  expect_true("native_greta_requires_greta_backend" %in% entries)
  expect_true("engine_pin_backend_conflict" %in% entries)
  expect_true("fa_rank_exceeds_dim" %in% entries)
  expect_false("grammar_greta_via_fb_deferred" %in% entries)

  # spot-check a representative new code from each family
  expect_true("precision_not_symmetric" %in% entries) # structured cov
  expect_true("low_rank_requires_greta" %in% entries) # approximation
  expect_true("residual_type_unsupported_for_aggregation" %in% entries) # aggregate emit
  expect_true("predict_kernel_invalid_include" %in% entries) # prediction
  expect_true("design_memory_exceeds_ceiling" %in% entries) # preflight
  expect_true("unsupported_family" %in% entries) # family gate
  expect_true("fb_cov_type_unknown" %in% entries) # fb_cov carrier
  expect_true("gretaR_not_installed" %in% entries) # gretaR backend
  expect_true("gretaR_cannot_represent_structured_cov" %in% entries) # gretaR gate

  # control-flow / routing reasons must NOT be registered
  expect_false("non_gaussian_family" %in% entries)
  expect_false("smooth_term_not_aggregatable" %in% entries)
  expect_false("policy_table_no_match_fallback_greta" %in% entries)
})

# ---------------------------------------------------------------- #
# Bonus: .lookup_refusal() returns the registered entry shape       #
# ---------------------------------------------------------------- #

test_that(".lookup_refusal() returns NULL for unregistered codes", {
  expect_null(flexyBayes:::.lookup_refusal("explicit_inla_gate_refused"))
  expect_null(flexyBayes:::.lookup_refusal("non_gaussian_family"))
  expect_null(flexyBayes:::.lookup_refusal("any_arbitrary_string"))
  expect_null(flexyBayes:::.lookup_refusal(""))
  expect_null(flexyBayes:::.lookup_refusal(NULL))
})

test_that(".lookup_refusal() returns the registered entry for the v0.3.10 codes", {
  entry <- flexyBayes:::.lookup_refusal("block_partition_incomplete")
  expect_type(entry, "list")
  expect_equal(entry$reason_code, "block_partition_incomplete")
  expect_equal(entry$registered_in_adr, "ADR 0025")
  expect_equal(entry$since_version, "0.3.10")
  expect_equal(entry$plan_field, "representation_plan")

  entry_low <- flexyBayes:::.lookup_refusal(
    "approximate_route_not_yet_registered"
  )
  expect_equal(entry_low$registered_in_adr, "ADR 0025+0027")
  expect_equal(entry_low$since_version, "0.3.10")
})
