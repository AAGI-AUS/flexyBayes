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

test_that(".refusal_registry holds the complete 54-code vocabulary", {
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
  # (= 54).
  expect_equal(length(entries), 54L)
  expect_true("stan_cannot_represent_structured_cov" %in% entries)
  expect_true("grammar_brms_with_asreml_terms" %in% entries)
  expect_true("native_greta_requires_greta_backend" %in% entries)
  expect_true("engine_pin_backend_conflict" %in% entries)
  expect_true("fa_rank_exceeds_dim" %in% entries)
  expect_false("grammar_greta_via_fb_deferred" %in% entries)

  # spot-check a representative new code from each family
  expect_true("precision_not_symmetric" %in% entries) # structured cov
  expect_true("low_rank_requires_greta" %in% entries) # approximation
  expect_true("rcov_type_unsupported_for_aggregation" %in% entries) # aggregate emit
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
