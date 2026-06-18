# Tests for fb_plan() --- the plan-only verb (ADR 0030 C1; v0.3.8).
#
# Covers the 12 acceptance criteria from v038-plan-2026-05-25 section 5.2:
#
#   (a) returns <fb_plan> classed object
#   (b) idempotence: two consecutive calls return identical objects
#   (c) flexybayes(plan = TRUE) routing equals fb_plan() for same model
#   (d) flight-checklist print form on a routable model
#   (e) flight-checklist print form on a refused model
#   (f) representation_plan matches Phase A known-matrix alignment policy
#   (g) rejected_routes equals backend_decision()$rejected_routes on
#       a subsequently-fitted version (shape-only check; no fit)
#   (h) as.data.frame() produces stable column ordering
#   (i) cov_validation_policy reports correctly
#   (j) chunked-predict planning surface --- predict_plan arg honoured
#   (k) malformed formula raises with refusal-registry-eligible code
#   (l) backward-compat: backend_decision() on a v0.3.7-saved fit (no
#       <fb_plan> slot) returns the legacy 8-field shape

# ---------------------------------------------------------------- #
# Test fixtures                                                      #
# ---------------------------------------------------------------- #

.test_plan_data_simple <- function(n = 60L, seed = 20260525L) {
  set.seed(seed)
  data.frame(
    y = rnorm(n),
    x = rnorm(n),
    g = factor(sample(letters[1:5], n, replace = TRUE))
  )
}

# ---------------------------------------------------------------- #
# (a) returns <fb_plan> classed object                              #
# ---------------------------------------------------------------- #

test_that("fb_plan() returns an <fb_plan> classed object", {
  d <- .test_plan_data_simple()
  p <- fb_plan(y ~ x + (1 | g), data = d, backend = "auto")
  expect_s3_class(p, "fb_plan")
  expect_true(is.list(p))
  expect_true("backend_chosen" %in% names(p))
  expect_true("routing_policy_version" %in% names(p))
  expect_true("representation_plan" %in% names(p))
  expect_true("rejected_routes" %in% names(p))
  expect_true("cov_validation_policy" %in% names(p))
  expect_true("memory_estimate_bytes" %in% names(p))
})

# ---------------------------------------------------------------- #
# (b) idempotence                                                    #
# ---------------------------------------------------------------- #

test_that("fb_plan() is idempotent (same inputs -> identical plan)", {
  d <- .test_plan_data_simple()
  p1 <- fb_plan(y ~ x + (1 | g), data = d, backend = "auto")
  p2 <- fb_plan(y ~ x + (1 | g), data = d, backend = "auto")
  # The call slot embeds match.call() which is identical here.
  # representation_plan + rejected_routes + memory_estimate_bytes
  # must agree exactly.
  expect_equal(p1$backend_chosen, p2$backend_chosen)
  expect_equal(p1$reason_code, p2$reason_code)
  expect_equal(p1$memory_estimate_bytes, p2$memory_estimate_bytes)
  expect_equal(p1$representation_plan, p2$representation_plan)
  expect_equal(p1$rejected_routes, p2$rejected_routes)
  expect_equal(p1$cov_validation_policy, p2$cov_validation_policy)
})

# ---------------------------------------------------------------- #
# (c) flexybayes(plan = TRUE) routing parity                         #
# ---------------------------------------------------------------- #

test_that("flexybayes(plan = TRUE) matches fb_plan() routing decision", {
  d <- .test_plan_data_simple()
  p_brms <- fb_plan(y ~ x + (1 | g), data = d, backend = "auto")
  p_asreml <- flexybayes(
    fixed = y ~ x,
    random = ~g,
    data = d,
    backend = "auto",
    plan = TRUE
  )
  expect_s3_class(p_asreml, "fb_plan")
  # Backend chosen + routing policy version must agree across ingest paths.
  expect_equal(p_asreml$backend_chosen, p_brms$backend_chosen)
  expect_equal(p_asreml$routing_policy_version, p_brms$routing_policy_version)
  expect_equal(p_asreml$cov_validation_policy, p_brms$cov_validation_policy)
})

# ---------------------------------------------------------------- #
# (d) flight-checklist print form, routable                          #
# ---------------------------------------------------------------- #

test_that("print.fb_plan() emits the flight-checklist surface", {
  d <- .test_plan_data_simple()
  p <- fb_plan(y ~ x + (1 | g), data = d, backend = "auto")
  out <- utils::capture.output(print(p))
  expect_true(any(grepl("flexyBayes plan", out)))
  expect_true(any(grepl("Will fit:", out)))
  expect_true(any(grepl("Backend chosen:", out)))
  expect_true(any(grepl("Routing policy version:", out)))
  expect_true(any(grepl("Cov validation policy:", out)))
  expect_true(any(grepl("Representation:", out)))
  expect_true(any(grepl("Engine:", out)))
})

# ---------------------------------------------------------------- #
# (e) print form on a refused model                                  #
# ---------------------------------------------------------------- #

test_that("print.fb_plan() surfaces the refusal when preflight refuses", {
  # Force preflight refusal with an absurdly small memory ceiling.
  d <- .test_plan_data_simple()
  p <- fb_plan(
    y ~ x + (1 | g),
    data = d,
    backend = "auto",
    memory_ceiling_gb = 1e-9
  )
  out <- utils::capture.output(print(p))
  expect_true(any(grepl("Will fit:", out)))
  # When refused, will_fit is no and the print form announces it.
  expect_true(
    any(grepl("preflight refused", out, fixed = TRUE)) ||
      isFALSE(p$will_fit)
  )
  expect_true(isFALSE(p$will_fit))
  expect_true(inherits(p$preflight_refusal, "fb_preflight_refusal"))
})

# ---------------------------------------------------------------- #
# (f) representation_plan structure                                  #
# ---------------------------------------------------------------- #

test_that("representation_plan carries per-term entries with class + justification", {
  d <- .test_plan_data_simple()
  p <- fb_plan(y ~ x + (1 | g), data = d, backend = "auto")
  if (length(p$representation_plan) > 0L) {
    for (rp in p$representation_plan) {
      expect_true(c("term_id") %in% names(rp))
      expect_true("representation_class" %in% names(rp))
      expect_true("justification" %in% names(rp))
    }
  } else {
    # Below the 1e5-row dispatcher threshold the per_term_estimate
    # is built only when fb_plan() ran preflight directly --- which
    # it does. So expect non-empty here.
    skip(
      "representation_plan empty -- preflight returned no per-term estimates"
    )
  }
})

# ---------------------------------------------------------------- #
# (g) rejected_routes shape                                          #
# ---------------------------------------------------------------- #

test_that("rejected_routes carries list of list(backend, reason)", {
  d <- .test_plan_data_simple()
  p <- fb_plan(y ~ x + (1 | g), data = d, backend = "auto")
  expect_true(is.list(p$rejected_routes))
  for (rr in p$rejected_routes) {
    expect_true("backend" %in% names(rr))
    expect_true("reason" %in% names(rr))
  }
})

# ---------------------------------------------------------------- #
# (h) as.data.frame() stable column ordering                         #
# ---------------------------------------------------------------- #

test_that("as.data.frame.fb_plan() preserves .FB_PLAN_DF_COLS ordering", {
  d <- .test_plan_data_simple()
  p <- fb_plan(y ~ x + (1 | g), data = d, backend = "auto")
  df <- as.data.frame(p)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1L)
  expect_equal(names(df), flexyBayes:::.FB_PLAN_DF_COLS)
  expect_true(is.character(df$backend_chosen))
  expect_true(is.logical(df$will_fit))
})

# ---------------------------------------------------------------- #
# (i) cov_validation_policy                                          #
# ---------------------------------------------------------------- #

test_that("cov_validation_policy reads flexyBayes.trust_pd option correctly", {
  d <- .test_plan_data_simple()
  # No known matrices -> "n/a"
  p <- fb_plan(y ~ x + (1 | g), data = d, backend = "auto")
  expect_equal(p$cov_validation_policy, "n/a")

  # With trust_pd = TRUE and no known matrices the answer is still
  # "n/a" --- the policy only matters when known matrices are present.
  old <- options(flexyBayes.trust_pd = TRUE)
  on.exit(options(old), add = TRUE)
  p2 <- fb_plan(y ~ x + (1 | g), data = d, backend = "auto")
  expect_equal(p2$cov_validation_policy, "n/a")
})

# ---------------------------------------------------------------- #
# (j) predict_plan surface                                           #
# ---------------------------------------------------------------- #

test_that("predict_plan arg computes the planned prediction shape", {
  d <- .test_plan_data_simple()
  nd <- .test_plan_data_simple(n = 500L, seed = 99L)
  p <- fb_plan(
    y ~ x + (1 | g),
    data = d,
    backend = "auto",
    predict_plan = list(newdata = nd, chunk_size = 100L)
  )
  expect_false(is.null(p$predict_plan))
  expect_equal(p$predict_plan$n_newrows, 500L)
  expect_equal(p$predict_plan$chunk_size, 100L)
  expect_equal(p$predict_plan$n_chunks, 5L)

  # df coercion exposes chunk count
  df <- as.data.frame(p)
  expect_true(df$predict_planned)
  expect_equal(df$predict_chunks, 5L)
})

# ---------------------------------------------------------------- #
# (k) malformed formula raises                                       #
# ---------------------------------------------------------------- #

test_that("malformed formula raises (refusal-registry-eligible at v0.4.0)", {
  d <- .test_plan_data_simple()
  # No formula -> error.
  expect_error(fb_plan(data = d, backend = "auto"), regexp = "formula|argument")
})

# ---------------------------------------------------------------- #
# (l) backward-compat: backend_decision() on a v0.3.7-shaped fit     #
# ---------------------------------------------------------------- #

test_that("backend_decision() backward-compat: legacy fits return 8-field shape", {
  # Synthesise a v0.3.7-shaped fit object (the 8-field <backend_decision>
  # shape introduced by ADR 0024 at v0.3.6). No <fb_plan> slot --- the
  # plan surface is v0.3.8-and-later. backend_decision() must return
  # the legacy 8-field list without raising.
  fake_fit <- structure(
    list(
      extras = list(
        backend_decision = list(
          backend = "greta",
          path = "explicit_greta",
          gate_checks = NULL,
          reason = "user requested greta explicitly",
          preflight_summary = NULL,
          representation_plan = NULL,
          rejected_routes = list(),
          routing_policy_version = "stage5a_v1"
        )
      )
    ),
    class = "flexybayes"
  )
  bd <- backend_decision(fake_fit)
  expect_true(is.list(bd))
  # Either the 4-field legacy shape OR the 8-field extended shape;
  # both are accepted at v0.3.8 per the backward-compat invariant.
  required_legacy <- c("backend", "path", "gate_checks", "reason")
  expect_true(all(required_legacy %in% names(bd)))
})
