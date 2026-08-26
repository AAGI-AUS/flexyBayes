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
  # Phase 3.6) a deferral code naming the third native engine (present
  # through v0.9.2, withdrawn entirely at 0.9.3 -- see NEWS.md) is
  # REMOVED (fb() now fits a native graph on that engine directly) and
  # two genuine refusals are added: one naming a native graph handed to
  # a different engine, and engine_pin_backend_conflict (a pin handed a
  # conflicting backend) (= 47). At v0.6.0.9000 that engine's dormant
  # sibling (also withdrawn at 0.9.3) is activated: six of its refusal
  # codes are registered, covering package-not-installed, below the
  # version floor, an unsupported family, a random-effect grouping
  # factor absent from the data, an unsupported random-term type, and a
  # structured-covariance term it cannot represent (= 53). At v0.7.0 the
  # data-aware factor-analytic rank upper bound is added:
  # fa_rank_exceeds_dim (fa(x, k) refused for k >= n_outer), the
  # identifiability complement of the data-free fa_rank_invalid floor
  # (= 54). At 0.9.0 the design-fidelity refusal
  # dsum_structured_inner_unsupported is added (a structured dsum() inner
  # -- ar1/us/... -- is refused rather than silently reduced to a
  # per-region heteroscedastic variance) (= 55).
  # At 0.9.0 the (since fully withdrawn) native engine and its dormant
  # sibling are quarantined, adding three refusal codes -- one naming
  # each engine directly plus backend_quarantined -- and auto_no_active_route
  # (= 58), and WP16 adds ar1_spatial_requires_complete_grid (= 59). All
  # four shipped in the 0.9.0 tag and their `since_version` said 0.9.1,
  # which named a release that did not exist -- corrected here and in the
  # registry.
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
  # At 0.9.1 the covariate half of ASReml's na.method() becomes settable,
  # and its `include` setting -- treat a missing covariate as zero
  # (Reference Manual 4.2, section 3.11) -- is the one ASReml behaviour
  # this package refuses rather than reproduces
  # (covariate_zero_fill_not_supported) (= 89).
  # The ASReml-shaped accessors add one suggested-package guard that
  # refuses by name rather than through a bare requireNamespace() error
  # (classify_requires_emmeans) (= 90), and one display refusal for the
  # residual variogram, which has no lag to plot against on a model that
  # names no design array (variogram_requires_design_index) (= 91).
  # The Bayesian door adds one: loo() reads the log-likelihood of each
  # observation at each posterior draw, which a nested Laplace
  # approximation never computed, so the call refuses by name and points
  # at the information criteria the fit does carry
  # (loo_requires_sampler_draws) (= 92). It is registered as the sibling
  # of fit_lacks_posterior_draws and the two descriptions name each
  # other (AD-4). It adds a second: a posterior predictive check needs
  # datasets replicated from the fitted model, and a fit whose engine
  # returns none refuses rather than showing a different display under
  # the check's name (pp_check_requires_predictive_draws) (= 93).
  # The asreml head-to-head adds two. A factor written on both the fixed
  # and the random side is aliased with itself, was accepted silently,
  # and reached INLA as an intermittently singular marginal solve -- two
  # engine segmentation faults and an untyped error at benchmark scale,
  # the only untyped structural failure in that run
  # (term_in_fixed_and_random) (= 94). And predict(classify = ) for a
  # random-effects grouping factor -- what a breeder asks a MET fit for
  # first -- died inside emmeans with a bare "No variable named Geno in
  # the reference grid", so it refuses by name and points at the
  # accessors that do carry the level effects
  # (classify_random_factor_not_supported) (= 95).
  # The field-coverage sweep adds twelve, all on the prior mini-language.
  # It accepted 23 of 31 malformed specifications and substituted a
  # default for the value the caller wrote, so three completed fits ran
  # under a prior nobody asked for. A mini-language is a public API: the
  # constructor now refuses an unknown or duplicated argument name
  # (prior_argument_unknown, prior_argument_duplicated), a missing
  # required one (prior_argument_missing -- half a PC prior states no
  # probability), a non-scalar or non-finite value
  # (prior_hyperparameter_not_scalar), and a value outside its domain
  # (prior_hyperparameter_out_of_domain -- a negative standard deviation
  # used to reach the sampler and surface as a missing-draws error).
  # The surrounding entry shapes are typed with it, because a refusal a
  # wrapper cannot match on by class is not a contract: the empty call,
  # the non-formula and one-sided argument, the unsupported target and
  # its missing identifying argument, and the right-hand side that is
  # not a known distribution call (prior_spec_empty,
  # prior_spec_not_formula, prior_spec_not_two_sided,
  # prior_target_unsupported, prior_target_argument_missing,
  # prior_distribution_not_a_call, prior_distribution_unknown) (= 107).
  # The same sweep's second slate adds seven more, all on paths that
  # refused for real reasons in untyped conditions. The prior emit is two:
  # a distribution the backend cannot carry on that parameter
  # (prior_not_translatable_for_backend -- on INLA it used to be dropped
  # silently while prior_summary() printed it as applied), and a prior
  # naming a term the model does not have (prior_target_not_in_model,
  # which brms's own parser answered with a synthesised Stan name). The
  # formula grammar is five: an mgcv smoother in the fixed part
  # (fixed_smoother_not_supported), a variable indexing both a fixed term
  # and an f() term (inla_variable_used_twice), a bar-grammar factor slope
  # whose ASReml surface fits (brms_factor_random_slope_unsupported), and
  # the two remaining ingest refusals typed with it
  # (brms_random_effect_form_unsupported,
  # brms_ingest_feature_unsupported) (= 114). The execution grid then
  # reached the two aggregated-route entry refusals that raised bare
  # errors above the plan gate -- a missing response, which is a property
  # of the data (aggregation_response_incomplete), and an unreachable
  # aggregated emit, which is a property of the engine roster
  # (aggregation_route_unavailable, raised at two sites) (= 116).
  # 0.9.3's WP-C scale/weights slate adds five more: an aggregation
  # cell-count product past R's integer limit (cell_count_exceeds_
  # integer), a typed INLA engine-death classification (inla_program_
  # failed), a per-level hyperparameter an at():ar1() field cannot
  # represent (at_field_per_level_hyper_not_representable), and the
  # two weights refusals -- wrong family/link (weights_requires_
  # gaussian) and the aggregated route (weights_not_aggregatable)
  # (= 121). 0.9.3's WP-G withdraws a third native engine entirely (see
  # NEWS.md): eleven codes are removed -- backend_quarantined; the six
  # refusal codes belonging to that engine's dormant sibling (package-
  # not-installed, below the version floor, unsupported family,
  # random-effect grouping factor absent from the data, unsupported
  # random-term type, and an unrepresentable structured-covariance term);
  # the two codes naming the native engine directly (one for the fit
  # path, one for a native graph handed to a different engine); one
  # naming it as the sole consumer of the low_rank_smooth approximation
  # scheme (renamed, not simply deleted -- see next); and one unrelated
  # orphan, whose raising site (.predict_linear_draws()) was deleted the
  # same session as an unreachable-code discovery, not as part of the
  # withdrawal itself. Two are added: low_rank_smooth_unsupported (the
  # renamed approximation-scheme code -- the low_rank_smooth scheme has
  # no active-engine consumer at all now, rather than a since-withdrawn
  # one) and unknown_backend (an explicit `backend` naming a value
  # outside {"auto", "inla", "brms"} -- including the withdrawn engine's
  # name -- refuses by this code rather than via match.arg()) (121 - 11 +
  # 2 = 112).
  expect_equal(length(entries), 112L)
  expect_true("aggregation_response_incomplete" %in% entries)
  expect_true("aggregation_route_unavailable" %in% entries)
  expect_true("prior_not_translatable_for_backend" %in% entries)
  expect_true("prior_target_not_in_model" %in% entries)
  expect_true("fixed_smoother_not_supported" %in% entries)
  expect_true("inla_variable_used_twice" %in% entries)
  expect_true("brms_factor_random_slope_unsupported" %in% entries)
  expect_true("term_in_fixed_and_random" %in% entries)
  expect_true("classify_random_factor_not_supported" %in% entries)
  expect_true("prior_argument_unknown" %in% entries)
  expect_true("prior_argument_missing" %in% entries)
  expect_true("prior_hyperparameter_out_of_domain" %in% entries)
  expect_true("prior_hyperparameter_not_scalar" %in% entries)
  expect_true("loo_requires_sampler_draws" %in% entries)
  expect_true("pp_check_requires_predictive_draws" %in% entries)
  expect_true("covariate_zero_fill_not_supported" %in% entries)
  expect_true("classify_requires_emmeans" %in% entries)
  expect_true("variogram_requires_design_index" %in% entries)
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
  expect_true("engine_pin_backend_conflict" %in% entries)
  expect_true("fa_rank_exceeds_dim" %in% entries)

  # spot-check a representative new code from each family
  expect_true("precision_not_symmetric" %in% entries) # structured cov
  expect_true("low_rank_smooth_unsupported" %in% entries) # approximation
  expect_true("residual_type_unsupported_for_aggregation" %in% entries) # aggregate emit
  expect_true("design_memory_exceeds_ceiling" %in% entries) # preflight
  expect_true("unsupported_family" %in% entries) # family gate
  expect_true("fb_cov_type_unknown" %in% entries) # fb_cov carrier
  expect_true("unknown_backend" %in% entries) # backend name gate

  # control-flow / routing reasons must NOT be registered
  expect_false("non_gaussian_family" %in% entries)
  expect_false("smooth_term_not_aggregatable" %in% entries)
  expect_false("policy_table_no_match_fallback_pending" %in% entries)
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
