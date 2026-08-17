# =============================================================================
# Residual-structure fidelity: a requested residual covariance must appear in
# the emitted model, or the call must refuse.
#
# These tests exist because a green suite of 2571 tests did not catch a
# requested AR1(row) x AR1(col) residual being emitted as an intercept-only
# iid Gaussian Stan program. Every test that touched the path asserted on the
# RETURNED OBJECT (its class, its recorded backend) rather than on the model
# that was actually generated -- and the returned object looked perfect.
#
# The discipline here is therefore: assert on the emitted code string, and on
# refusal reason codes. A wrapper's class is not evidence about its contents.
#
# brms has no residual-structure lowering at all: .fb_to_brms_formula()
# reconstructs the fixed and random blocks only. So for brms the fidelity
# requirement can only be met by refusing.
#
# 0.9.0 respelling: the separable AR1 field is written on the RANDOM side,
# because that is what the INLA emit builds -- a latent autoregressive field
# plus the Gaussian observation nugget. The residual spelling names ASReml's
# three-parameter nugget-free process, which neither engine fits, and now
# refuses on every route. Both halves are asserted here.
# =============================================================================

.mk_resid_ir <- function(residual_terms) {
  flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    link = "identity",
    fixed_terms = list(),
    random_terms = list(),
    residual_terms = residual_terms,
    source = "asreml"
  )
}

# A complete 5 x 4 grid, one observation per (row, col) node.
.complete_grid <- function(n_row = 5L, n_col = 4L, seed = 11L) {
  set.seed(seed)
  g <- expand.grid(row = factor(seq_len(n_row)), col = factor(seq_len(n_col)))
  g$y <- stats::rnorm(nrow(g))
  g$x <- stats::rnorm(nrow(g))
  g
}

# ---------------------------------------------------------------- #
# 1. The capability predicate sees residual terms                   #
# ---------------------------------------------------------------- #

test_that(".capability_brms() inspects residual_terms, not just random_terms", {
  # The iid `units` residual is the one form brms carries -- as the family
  # scale parameter sigma.
  expect_true(
    isTRUE(flexyBayes:::.capability_brms(
      .mk_resid_ir(list(list(type = "units")))
    ))
  )
  expect_true(isTRUE(flexyBayes:::.capability_brms(.mk_resid_ir(list()))))

  # Structured residuals have no lowering and must be reported as such.
  expect_identical(
    flexyBayes:::.capability_brms(
      .mk_resid_ir(list(list(type = "ar1", var = "row")))
    ),
    "stan_cannot_represent_ar1_residual"
  )
  expect_identical(
    flexyBayes:::.capability_brms(
      .mk_resid_ir(list(list(
        type = "ar1_spatial", row_var = "row", col_var = "col"
      )))
    ),
    "stan_cannot_represent_ar1_residual"
  )
})

test_that("an unrecognised residual type is refused rather than dropped", {
  # The allowlist is deliberately positive: a residual type the parser
  # learns later must default to a refusal, not to silent omission.
  expect_identical(
    flexyBayes:::.capability_brms(
      .mk_resid_ir(list(list(type = "some_future_residual_form")))
    ),
    "stan_cannot_represent_structured_residual"
  )
})

# ---------------------------------------------------------------- #
# 2. Every route refuses rather than emitting a different model     #
# ---------------------------------------------------------------- #

test_that("explicit brms refuses a random-side AR1 field", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .complete_grid()
  expect_error(
    suppressMessages(flexybayes(
      y ~ 1, random = ~ ar1(row):ar1(col), data = g,
      backend = "brms", return_code = TRUE
    )),
    class = "flexybayes_refusal_stan_cannot_represent_ar1_field"
  )
})

test_that("the residual AR1 spelling refuses on every route", {
  # ASReml's separable residual is a three-parameter nugget-free process.
  # The INLA emit builds a latent field plus the observation nugget -- four
  # parameters, a different model -- so the spelling refuses rather than
  # being fitted under a name it does not have, whichever engine is asked.
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .complete_grid()
  for (engine in c("inla", "brms", "auto")) {
    expect_error(
      suppressMessages(flexybayes(
        y ~ 1, residual = ~ ar1(row):ar1(col), data = g,
        backend = engine, return_code = TRUE
      )),
      class = "flexybayes_refusal_ar1_residual_not_representable"
    )
  }
  # The 1-D spelling refuses on the same grounds.
  expect_error(
    suppressMessages(flexybayes(
      y ~ 1, residual = ~ ar1(row), data = g, backend = "inla",
      return_code = TRUE
    )),
    class = "flexybayes_refusal_ar1_residual_not_representable"
  )
})

test_that("auto + return_code refuses when brms cannot represent the model", {
  # The code-inspection modes resolve backend = "auto" to brms, since brms
  # is the only active code-producing engine. That resolution must be gated
  # on capability: returning brms code for a model brms cannot express
  # answers a different question from the one asked.
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .complete_grid()
  expect_error(
    suppressMessages(flexybayes(
      y ~ 1, random = ~ ar1(row):ar1(col), data = g,
      backend = "auto", return_code = TRUE
    )),
    class = "flexybayes_refusal_auto_no_active_route"
  )
})

test_that("auto + review_code refuses when brms cannot represent the model", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .complete_grid()
  expect_error(
    suppressMessages(flexybayes(
      y ~ 1, random = ~ ar1(row):ar1(col), data = g,
      backend = "auto", review_code = TRUE
    )),
    class = "flexybayes_refusal_auto_no_active_route"
  )
})

# ---------------------------------------------------------------- #
# 3. Emit-level positive control: the structure is really there     #
# ---------------------------------------------------------------- #

test_that("INLA emits the separable AR1xAR1 field for a complete grid", {
  # The positive half of the contract. Refusing everything would pass every
  # test above; this asserts the structure survives into the emitted model.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .complete_grid()
  code <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(row):ar1(col), data = g,
    backend = "inla", return_code = TRUE
  ))
  f <- paste(deparse(code$formula), collapse = " ")
  expect_match(f, "model\\s*=\\s*\"ar1\"")
  expect_match(f, "group\\s*=")
  expect_match(f, "control\\.group")
})

test_that("brms code for an unstructured model still carries its terms", {
  # Guards the other direction: the capability gate must not have closed
  # routes that were correct. A plain fixed + random model still emits.
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .complete_grid()
  g$blk <- factor(rep(1:4, length.out = nrow(g)))
  code <- suppressMessages(flexybayes(
    y ~ x, random = ~ blk, data = g, backend = "auto", return_code = TRUE
  ))
  txt <- paste(as.character(code), collapse = "\n")
  expect_match(txt, "Intercept")
  # brms names grouping-factor data blocks J_1 / N_1 / Z_1.
  expect_match(txt, "J_1|N_1|Z_1")
})
