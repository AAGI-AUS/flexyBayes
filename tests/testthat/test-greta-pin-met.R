# PR3: an explicit `backend = "greta"` request is the opt-in path by which
# flexyBayes fits the MET structures INLA (and therefore the auto default)
# refuses -- crossed interaction random effects and a heteroscedastic
# per-environment residual. The auto default keeps its honest structural
# refusal; only an explicit greta pin plans as will_fit. These are
# plan-level tests (no greta install required).

test_that("an explicit greta pin plans will_fit for interaction + heteroscedastic models", {
  set.seed(1)
  d <- expand.grid(g = factor(1:20), e = factor(1:6), rep = 1:4)
  d$y <- rnorm(nrow(d))

  # g:e crossed random effect + a per-e heteroscedastic residual: INLA
  # cannot represent either, but the greta codegen gathers both.
  p_greta <- flexybayes(
    y ~ e, random = ~ g + g:e, residual = ~ dsum(~ units | e),
    data = d, plan = TRUE, backend = "greta"
  )
  expect_true(p_greta$will_fit)

  # the same model on the auto default keeps the honest INLA refusal
  p_auto <- flexybayes(
    y ~ e, random = ~ g + g:e, residual = ~ dsum(~ units | e),
    data = d, plan = TRUE, backend = "auto"
  )
  expect_false(p_auto$will_fit)
  expect_identical(p_auto$gate_outcome, "refuse_structural")

  # a simple random intercept is unaffected on both routes
  expect_true(flexybayes(y ~ e, random = ~ g, data = d, plan = TRUE)$will_fit)
})

test_that(".preflight_random_term sizes nested and combo interaction terms", {
  d <- expand.grid(a = factor(1:5), b = factor(1:4), cc = factor(1:3))
  ds <- .fb_dataset(d)

  # combo A:B:C -- an index gather into a k-length latent vector
  combo <- .preflight_random_term(
    list(type = "combo", vars = c("a", "b", "cc")),
    n_rows = nrow(d), fb_dataset = ds, family_ok = TRUE
  )
  expect_false(isTRUE(combo$unknown_representation))
  expect_true(is.finite(combo$design_memory_bytes))

  # nested A:B -- inner / outer fields
  nested <- .preflight_random_term(
    list(type = "nested", inner = "a", outer = "b"),
    n_rows = nrow(d), fb_dataset = ds, family_ok = TRUE
  )
  expect_false(isTRUE(nested$unknown_representation))
  expect_true(is.finite(nested$design_memory_bytes))
})
