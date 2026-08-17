# =============================================================================
# Observation weights: refused, not ignored.
#
# `weights` were parsed into an addition term by every ingest adapter and
# consumed by no active emitter. The brms Stan program was byte-identical
# with and without them; the INLA formula, family, hyperparameters and
# family control were identical too; and the standard-R compatibility view
# reported weights = 1 regardless. A weighted call therefore returned the
# UNWEIGHTED posterior, correctly formatted, with nothing to indicate it.
#
# The word covers several distinct models -- inverse-residual-variance (the
# ASReml sense the documentation promised), frequency, likelihood-power,
# trials, exposure. They are not interchangeable: for a Gaussian with
# unknown sigma, raising each likelihood contribution to w_i is not the
# normalised model Var(y_i) = sigma^2 / w_i. So the argument stays refused
# until each mapping is implemented and checked against an analytic oracle.
#
# The IR still records the request -- the refusal is a dispatch-level gate,
# not a parser change, so the ingest adapters keep their coverage.
# =============================================================================

.wdat <- function(n = 8L, seed = 3L) {
  set.seed(seed)
  data.frame(y = stats::rnorm(n), x = stats::rnorm(n))
}

test_that("a non-constant weight vector is refused on every active backend", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wdat()
  w <- as.numeric(seq_len(nrow(d)))

  for (be in c("auto", "inla", "brms")) {
    expect_error(
      suppressMessages(flexybayes(
        y ~ x, data = d, weights = w, backend = be, return_code = TRUE
      )),
      class = "flexybayes_refusal_weights_not_supported",
      label = paste0("backend = ", be)
    )
  }
})

test_that("the refusal fires before the fit, not after", {
  # return_code = TRUE never reaches a sampler, so a refusal here proves
  # the gate sits ahead of emission rather than being a post-hoc check on
  # a fitted object.
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wdat()
  expect_error(
    suppressMessages(flexybayes(
      y ~ x, data = d, weights = as.numeric(seq_len(nrow(d))),
      backend = "brms", return_code = TRUE
    )),
    class = "flexybayes_refusal_weights_not_supported"
  )
})

test_that("constant weights are the unweighted model and pass through", {
  # rep(1, n) -- and any all-equal vector -- specifies the same model as
  # no weights at all, so refusing it would be gratuitous.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wdat()
  expect_no_error(
    suppressMessages(flexybayes(
      y ~ x, data = d, weights = rep(1, nrow(d)),
      backend = "inla", return_code = TRUE
    ))
  )
})

test_that("NULL weights are unaffected", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .wdat()
  expect_no_error(
    suppressMessages(flexybayes(
      y ~ x, data = d, backend = "inla", return_code = TRUE
    ))
  )
})

test_that("the ingest adapters still record the requested weights", {
  # The refusal is a dispatch gate. The representation must keep the
  # request, so that implementing the semantics later is a matter of
  # consuming a term that is already there.
  d <- .wdat(5L)
  w <- c(1, 2, 3, 4, 5)
  ir <- flexyBayes:::fb_from_asreml(y ~ x, data = d, weights = w)
  types <- vapply(ir$addition_terms, function(t) t$type, character(1L))
  expect_true("weights" %in% types)
  expect_equal(ir$addition_terms[[match("weights", types)]]$values, w)
})

test_that("the guard reads the IR, not only the dispatch argument", {
  # An IR carrying weights must be refused even when the dispatch-level
  # `weights` argument is NULL -- otherwise a caller that built the IR
  # upstream would slip past the gate.
  d <- .wdat(5L)
  ir <- flexyBayes:::fb_from_asreml(
    y ~ x, data = d, weights = c(1, 2, 3, 4, 5)
  )
  expect_error(
    flexyBayes:::.refuse_unsupported_weights(ir, weights = NULL),
    class = "flexybayes_refusal_weights_not_supported"
  )
})
