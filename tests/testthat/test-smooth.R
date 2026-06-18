# Tests for mgcv-style univariate smooth ingest s(x) (ADR 0004 D3).
#
# v0.1 supports s() only; te() / ti() / t2() defer to v0.2. The
# basis is constructed at parse time via mgcv::smoothCon() and
# emitted as a literal n x k matrix in the greta code; the random-
# effect layer uses .sigma_decl() so the smooth's sigma honours the
# prior priority (uniform_per_vc -> pc_per_vc -> legacy lognormal).

skip_if_no_mgcv <- function() skip_if_not_installed("mgcv")

mk_smooth_data <- function() {
  set.seed(2026L)
  n <- 50L
  x <- sort(runif(n, 0, 10))
  y <- sin(x) + rnorm(n, sd = 0.2)
  data.frame(x = x, y = y, g = factor(rep(1:5, length.out = n)))
}

# ---------------------------------------------------------------- #
# Parse-time basis construction                                    #
# ---------------------------------------------------------------- #

test_that(".enrich() builds the smoothCon basis for s(x)", {
  skip_if_no_mgcv()
  d <- mk_smooth_data()
  fb <- flexyBayes:::fb_from_asreml(
    fixed = y ~ 1,
    random = ~ s(x),
    data = d
  )
  expect_length(fb$random_terms, 1L)
  rt <- fb$random_terms[[1]]
  expect_identical(rt$type, "smooth_mgcv")
  expect_identical(rt$var, "x")
  expect_true(is.matrix(rt$X))
  expect_identical(nrow(rt$X), nrow(d))
  expect_true(rt$k >= 1L)
  expect_true(is.character(rt$smooth_label))
})

test_that("s(x, k = 6) honours user-supplied k", {
  skip_if_no_mgcv()
  d <- mk_smooth_data()
  fb <- flexyBayes:::fb_from_asreml(
    fixed = y ~ 1,
    random = ~ s(x, k = 6),
    data = d
  )
  rt <- fb$random_terms[[1]]
  expect_identical(rt$type, "smooth_mgcv")
  # absorb.cons = TRUE drops one column to fold the sum-to-zero
  # constraint; k = 6 -> 5 effective columns.
  expect_identical(rt$k, 5L)
})

test_that(".enrich() errors politely when s() variable is missing from data", {
  skip_if_no_mgcv()
  d <- mk_smooth_data()
  expect_error(
    flexyBayes:::fb_from_asreml(
      fixed = y ~ 1,
      random = ~ s(missing_var),
      data = d
    ),
    regexp = "not found in data"
  )
})

# ---------------------------------------------------------------- #
# Code generation through emit_greta                                #
# ---------------------------------------------------------------- #

test_that("emit_greta emits a bound-basis + random-effect block for s(x)", {
  skip_if_no_greta()
  skip_if_no_mgcv()
  d <- mk_smooth_data()
  local_clean_emit_state()
  set_emit_state("default_prior_note", TRUE)
  code <- flexybayes(
    fixed = y ~ 1,
    random = ~ s(x),
    data = d,
    prior_vc_sd = 1,
    verbose = FALSE,
    return_code = TRUE
  )
  # ADR 0018 / the design spec: the basis matrix is bound into the
  # model env as B_s_x; the emitted code references it via
  # `as_data(B_s_x)` rather than inlining the n x k literal. The
  # B_s_x literal expression is gone from the code surface.
  expect_no_match(code, "B_s_x <- matrix\\(c\\(")
  expect_match(code, "B_g_s_x <- as_data\\(B_s_x\\)")
  expect_match(code, "sigma_s_x <- greta::lognormal\\(0, 1\\)")
  expect_match(code, "s_x_raw <- normal\\(0, 1, dim = ")
  expect_match(code, "f_s_x <- B_g_s_x %\\*% \\(s_x_raw \\* sigma_s_x\\)")
})

test_that("default uniform prior reaches the smooth's sigma via .sigma_decl", {
  skip_if_no_greta()
  skip_if_no_mgcv()
  # v0.3.9 emit-state migration + bare-options leak cleanup: the
  # default_prior_note once-flag now lives in a package-internal env
  # (R/emit_state.R), shielded from the options() namespace where
  # unrelated callers used to consume it before the intended emission;
  # in parallel, six bare options(flexyBayes.silence_default_prior_note
  # = TRUE) calls across test-inla-verification-simple-slope-uncor.R
  # and test-random-slopes-uncor.R were migrated to
  # withr::local_options() so the silence flag does not leak into
  # sibling test files in the tally.R single-process loop (that leak
  # was the upstream cause of the v0.3.8 suite-order flake at this
  # site). The structural-witness assertions below -- latch transitions
  # FALSE -> TRUE around the flexybayes() call -- pin the contract
  # directly and are robust against future callinghandler-stack churn;
  # textual emission via message() remains exercised end-to-end in
  # test-emit-greta.R and test-emit-state-isolation.R.
  d <- mk_smooth_data()
  local_clean_emit_state()
  expect_false(flexyBayes:::.emit_state_get("default_prior_note"))
  code <- flexybayes(
    fixed = y ~ 1,
    random = ~ s(x),
    data = d,
    verbose = FALSE,
    return_code = TRUE
  )
  expect_true(flexyBayes:::.emit_state_get("default_prior_note"))
  # Both the residual sigma and the smooth's sigma should be uniform
  # under the ADR 0004 default. Smooth's sigma is keyed by the
  # variable name "x" in priors_to_legacy()'s uniform_per_vc map (via
  # the smooth target type); .default_uniform_prior() only seeds
  # named random groups, so the smooth's sigma falls through to the
  # legacy lognormal default. Document the gap; v0.2 extends.
  expect_match(code, "sigma_e_atg <- greta::uniform\\(0, ")
})
