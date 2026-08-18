# =============================================================================
# A term on both the fixed and the random side: refused, not attempted.
#
# `yield ~ Variety` with `random = ~ Block + Variety` writes Variety twice.
# The fixed part already estimates one mean per level, so the random copy's
# deviations are identified by their prior alone -- the variance component
# reported for them reads the prior, not the data.
#
# It is refused rather than warned about because every engine answers it
# badly and each does so differently. On INLA the marginal solve goes
# intermittently singular: an empty variance-component table on 4 of 8
# small-n runs, and at benchmark scale two `Segmentation fault: 11` lines
# from the inla.run child followed by an untyped error carrying no reason
# code -- the only untyped structural failure the ASReml head-to-head found.
# brms completes and returns a posterior in which the two term sets are
# estimating one quantity between them. ASReml accepts the spelling and
# fits it, so a script translated from ASReml arrives here unchanged.
#
# Scope is the plain aliased pair on purpose: a factor main effect on the
# fixed side against an intercept-only random main effect on the same
# variable. A random slope over a fixed grouping factor and a random
# interaction involving a fixed main effect are both models a user means,
# and neither fires -- pinned below.
# =============================================================================

.aliased_trial <- function(n_var = 6L, n_blk = 4L, seed = 20260817L) {
  set.seed(seed)
  d <- expand.grid(
    Variety = factor(sprintf("V%02d", seq_len(n_var))),
    Block = factor(sprintf("B%d", seq_len(n_blk)))
  )
  d$Row <- factor(rep(seq_len(n_blk), length.out = nrow(d)))
  d$yield <- 5 +
    rep(stats::rnorm(n_var, 0, 1), times = n_blk) +
    rep(stats::rnorm(n_blk, 0, 0.8), each = n_var) +
    stats::rnorm(nrow(d), 0, 0.5)
  d
}

# ---------------------------------------------------------------- #
# 1. The refusal fires, on every route into the package             #
# ---------------------------------------------------------------- #

test_that("the aliased model is refused on every backend, before any fit", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .aliased_trial()

  # return_code = TRUE never reaches a sampler, so a refusal on that route
  # proves the gate sits ahead of emission rather than reading a fit.
  for (be in c("auto", "inla", "brms")) {
    expect_error(
      suppressMessages(flexybayes(
        yield ~ Variety,
        random = ~ Block + Variety,
        data = d,
        backend = be,
        return_code = TRUE
      )),
      class = "flexybayes_refusal_term_in_fixed_and_random",
      label = paste0("backend = ", be)
    )
  }
})

test_that("the refusal carries the reason code and names the term", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .aliased_trial()

  cnd <- tryCatch(
    suppressMessages(flexybayes(
      yield ~ Variety,
      random = ~ Block + Variety,
      data = d,
      return_code = TRUE
    )),
    error = function(e) e
  )

  expect_s3_class(cnd, "flexybayes_refusal_term_in_fixed_and_random")
  expect_s3_class(cnd, "flexybayes_refusal")
  expect_identical(cnd$reason_code, "term_in_fixed_and_random")

  msg <- conditionMessage(cnd)
  # Names the offending term, states the aliasing, and gives BOTH remedies.
  expect_match(msg, "Variety", fixed = TRUE)
  expect_match(msg, "both the fixed and the random side", fixed = TRUE)
  expect_match(msg, "drop `Variety` from `random`", fixed = TRUE)
  expect_match(msg, "drop it from the fixed part", fixed = TRUE)
  expect_match(msg, "ASReml accepts this spelling", fixed = TRUE)
})

test_that("id() and ide() spell the same aliased model and are refused", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .aliased_trial()

  for (spelling in c("id", "ide")) {
    rand <- stats::as.formula(paste0("~ Block + ", spelling, "(Variety)"))
    expect_error(
      suppressMessages(flexybayes(
        yield ~ Variety,
        random = rand,
        data = d,
        return_code = TRUE
      )),
      class = "flexybayes_refusal_term_in_fixed_and_random",
      label = spelling
    )
  }
})

test_that("update() into the aliased model refuses instead of segfaulting", {
  # The v0.9.1 execution spec's own exit-gate line. `yield ~ Variety` is
  # already fixed, so `update(fit, random = ~ Block + Variety)` asks for the
  # aliased model -- and used to reach INLA, which failed untyped.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .aliased_trial()

  fit <- suppressMessages(flexybayes(
    yield ~ Variety,
    random = ~Block,
    data = d,
    backend = "inla"
  ))

  expect_error(
    suppressMessages(update(fit, random = ~ Block + Variety)),
    class = "flexybayes_refusal_term_in_fixed_and_random"
  )
})

test_that("planning refuses it too, on both grammars", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .aliased_trial()

  # ASReml grammar, through flexybayes(plan = TRUE).
  expect_error(
    suppressMessages(flexybayes(
      yield ~ Variety,
      random = ~ Block + Variety,
      data = d,
      plan = TRUE
    )),
    class = "flexybayes_refusal_term_in_fixed_and_random"
  )

  # brms grammar, through fb_plan().
  expect_error(
    fb_plan(yield ~ Variety + (1 | Variety), data = d),
    class = "flexybayes_refusal_term_in_fixed_and_random"
  )
})

# ---------------------------------------------------------------- #
# 2. It does not fire on models a user means                        #
# ---------------------------------------------------------------- #

test_that("a well-posed randomised block model is untouched", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .aliased_trial()

  expect_no_error(suppressMessages(flexybayes(
    yield ~ Variety,
    random = ~Block,
    data = d,
    backend = "brms",
    return_code = TRUE
  )))
})

test_that("a random interaction involving a fixed main effect is untouched", {
  # `random = ~ Variety:Block` is not aliased with the Variety main effect
  # -- it is the interaction stratum, and refusing it would be wrong.
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .aliased_trial()

  expect_no_error(suppressMessages(flexybayes(
    yield ~ Variety,
    random = ~ Block + Variety:Block,
    data = d,
    backend = "brms",
    return_code = TRUE
  )))
})

test_that("a random effect on a factor absent from the fixed part passes", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .aliased_trial()

  expect_no_error(fb_plan(yield ~ Block + (1 | Variety), data = d))
})

# ---------------------------------------------------------------- #
# 3. Registry and documentation surfaces                            #
# ---------------------------------------------------------------- #

test_that("the code is registered and reachable from fb_refusals()", {
  entry <- flexyBayes:::.lookup_refusal("term_in_fixed_and_random")
  expect_type(entry, "list")
  expect_identical(entry$since_version, "0.9.1")
  expect_identical(entry$plan_field, "rejected_routes")

  tbl <- fb_refusals(reason_code = "term_in_fixed_and_random")
  expect_s3_class(tbl, "fb_refusals_table")
  expect_identical(nrow(tbl), 1L)
  expect_match(tbl$description, "aliases", fixed = TRUE)
})
