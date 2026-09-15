# =============================================================================
# fb_select(): the projection's metric is validated before anything runs.
#
# The defect this closes. `fb_decision(method = "projection", metric = "rmse")`
# was constructible, and the plain stop() that rejected it sat inside
# `.fb_sel_project()`, after the reference model had been sampled. The same
# held for `method = "auto"` resolving to the projection. The pair is now
# validated on the RESOLVED method, before the provider check, the budget and
# the first fit, and `plan = TRUE` returns the same verdict.
#
# No sampler is involved. The base is a parsed model with its data and no
# posterior, which is everything the path up to the check reads; the
# reference fit and the budget are replaced by functions that fail loudly, so
# a check that fired late would surface as their error rather than as the
# refusal.
# =============================================================================

fbdm_base <- function() {
  set.seed(11L)
  d <- data.frame(site = factor(rep(sprintf("S%02d", seq_len(8L)), each = 5L)))
  d$p1 <- stats::rnorm(nrow(d))
  d$x1 <- stats::rnorm(nrow(d))
  d$x2 <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$p1 + d$x1 +
    stats::rnorm(8L, sd = 0.7)[as.integer(d$site)] +
    stats::rnorm(nrow(d), sd = 0.5)
  fb <- fb_from_brms(y ~ p1 + (1 | site), data = d, family = "gaussian")
  fb$link <- fb$link %||% "identity"
  structure(
    list(fb = fb, data = d),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
}

fbdm_call <- function(base, decision, plan = FALSE) {
  fb_select(
    base, ~ x1 + x2,
    penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    resamples = fb_resamples(unit = ~site, v = 2L),
    decision = decision,
    plan = plan,
    seed = 1L
  )
}

fbdm_tripwires <- function(env = parent.frame()) {
  local_mocked_bindings(
    .fb_sel_reference_fit = function(...) {
      stop("a reference fit started before the metric was checked")
    },
    .fb_sel_cost_estimate = function(...) {
      stop("the budget was computed before the metric was checked")
    },
    .env = env
  )
}


test_that("an explicit projection on rmse refuses before the budget and any fit", {
  base <- fbdm_base()
  fbdm_tripwires()

  t0 <- proc.time()
  err <- tryCatch(
    fbdm_call(base, fb_decision(method = "projection", metric = "rmse")),
    error = function(e) e
  )
  elapsed <- unname((proc.time() - t0)["elapsed"])

  expect_s3_class(err, "flexybayes_refusal_projection_metric_unsupported")
  expect_identical(err$metric, "rmse")
  expect_identical(err$decision_method, "projection")
  expect_match(conditionMessage(err), "Nothing was fitted", fixed = TRUE)
  expect_match(conditionMessage(err), "metric = \"elpd\"", fixed = TRUE)
  expect_lt(elapsed, 1)
})


test_that("auto resolving to the projection is held to the same metric", {
  base <- fbdm_base()
  # The precondition the test needs: this model's gate is green, so `auto`
  # really does resolve to the projection. On a red gate it resolves to
  # `full`, where rmse is a valid out-of-fold metric and nothing refuses.
  expect_identical(
    .fb_select_gate_verdict(.fb_sel_structure_gate(base$fb)),
    "supported"
  )
  fbdm_tripwires()

  err <- tryCatch(
    fbdm_call(base, fb_decision(metric = "nlpd")),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_projection_metric_unsupported")
  expect_identical(err$requested_method, "auto")
  expect_identical(err$decision_method, "projection")
  expect_match(conditionMessage(err), "method = \"auto\"", fixed = TRUE)
})


test_that("plan = TRUE returns the verdict the run would", {
  base <- fbdm_base()
  fbdm_tripwires()

  for (decision in list(
    fb_decision(method = "projection", metric = "rmse"),
    fb_decision(metric = "rmse")
  )) {
    t0 <- proc.time()
    expect_error(
      fbdm_call(base, decision, plan = TRUE),
      class = "flexybayes_refusal_projection_metric_unsupported"
    )
    expect_lt(unname((proc.time() - t0)["elapsed"]), 1)
  }
})


test_that("the check refuses only the pair the projection cannot choose on", {
  # The pairs that pass the metric check reach the provider check, which
  # needs brms for the Bayesian workflow.
  skip_if_not_installed("brms")
  base <- fbdm_base()
  fbdm_tripwires()

  # A decision that does not choose on elpd keeps rmse as its out-of-fold
  # metric; so does the route with every candidate protected, which has no
  # decision layer at all.
  full_plan <- fbdm_call(
    base, fb_decision(method = "full", metric = "rmse"), plan = TRUE
  )
  expect_s3_class(full_plan, "fb_select_plan")
  expect_identical(full_plan$selection$decision, "full")

  null_plan <- fb_select(
    base, ~ x1 + x2,
    resamples = fb_resamples(unit = ~site, v = 2L),
    decision = fb_decision(method = "projection", metric = "rmse"),
    plan = TRUE, seed = 1L
  )
  expect_s3_class(null_plan, "fb_select_plan")

  # And the projection on elpd plans.
  skip_if_not_installed("brms")
  elpd_plan <- fbdm_call(
    base, fb_decision(method = "projection", metric = "elpd"), plan = TRUE
  )
  expect_s3_class(elpd_plan, "fb_select_plan")
  expect_identical(elpd_plan$selection$decision, "projection")
})


test_that("the validator itself refuses by class, and passes what it should", {
  expect_error(
    .fb_sel_require_projection_metric("projection", "rmse"),
    class = "flexybayes_refusal_projection_metric_unsupported"
  )
  expect_error(
    .fb_sel_require_projection_metric("projection", "nlpd"),
    class = "flexybayes_refusal_projection_metric_unsupported"
  )
  expect_true(.fb_sel_require_projection_metric("projection", "elpd"))
  expect_true(.fb_sel_require_projection_metric("full", "rmse"))
  expect_true(.fb_sel_require_projection_metric("posterior_threshold", "nlpd"))
})
