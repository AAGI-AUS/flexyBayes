# PR3, as amended after the greta quarantine: crossed interaction random
# effects and a heteroscedastic per-environment residual are the MET
# structures INLA cannot represent, and an explicit `backend = "greta"`
# request was once the only route that planned as will_fit. It is not any
# more -- the router resolves the same model to brms under `auto`, and
# brms fits it (tests/testthat/test-met-combined.R runs it). The plan
# surface must say so: reading the LGM gate's verdict on INLA as the
# fit's verdict told a user this model would not run and then ran it.
# These are plan-level tests (no greta install required).

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

  # The same model on the auto default: the INLA gate still refuses on
  # structure, and the route resolved past it to brms, which fits.
  skip_if_not_installed("brms")
  p_auto <- flexybayes(
    y ~ e, random = ~ g + g:e, residual = ~ dsum(~ units | e),
    data = d, plan = TRUE, backend = "auto"
  )
  expect_identical(p_auto$gate_outcome, "refuse_structural")
  expect_identical(p_auto$backend_chosen, "brms")
  expect_identical(p_auto$path, "auto_multistratum_to_brms")
  expect_true(p_auto$will_fit)

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
