# test-space-level-factor.R -- C4 / FS-26: the INLA emit dies untyped on
# a non-syntactic factor level.
#
# R/emit_inla.R died with `object 'NAPPLIED2low N' not found` -- INLA
# expands a factor fixed term to per-level names built as
# paste0(var, level) and evaluates them internally; a level containing a
# space breaks that. The fix legalises every factor the model uses on a
# copy of `data` before anything reaches INLA, records the map on the
# fit, and every reader of the fit's own labels (coef(), the tidier
# summary() reads, the random-effect table ranef() reads, the emmeans /
# marginaleffects design-matrix reconciliation predict(classify =) and
# confint() go through) restores the user's own labels from that map.
#
# Note on the live-crash reproduction: the raw INLA-program failure
# mode at small / unbalanced N is not deterministic in this environment
# -- repeated draws of an unbalanced n = 60 fixture with the identical
# space-containing level surfaced the exact FS-26 text
# ("object 'NAPPLIED2low N' not found") on some runs and an unrelated
# GSL numerical / subprocess-segfault failure (caught by C2's
# `inla_program_failed`) on others, evidently from the INLA subprocess's
# own multi-threaded solver. A live crash is therefore not a reliable
# test oracle for THIS specific defect. The mechanism-level test below
# (mocking `.inla_call()` to intercept the data before it ever reaches
# the engine, per the same pattern C2 uses) pins the actual claim
# deterministically: no factor level flexyBayes hands to INLA carries a
# space, or anything else make.names() would touch. The live-fit tests
# after it use a larger, balanced, signal-bearing design (spot-checked
# numerically stable across five seeds and three repeats) to prove the
# fix does not preclude a normal, working fit and that all four
# accessors round-trip correctly end to end.
#
# `aggregate = FALSE` is explicit throughout: the identical space-level
# design is also aggregation-eligible at small N, and the aggregated
# INLA route (R/emit_gaussian_aggregated.R) is a SEPARATE code path this
# slice does not touch -- see the WP-C report's Handoffs. These tests
# exercise the per-row emit_inla() fix this file owns.

skip_if_no_inla <- function() skip_if_not_installed("INLA")
skip_if_no_brms <- function() skip_if_not_installed("brms")

# Two crossed random-intercept groups (4 and 3 levels) with 8 reps per
# cell and a real fixed + random signal (not pure noise) -- a smaller /
# noisier version of this shape hit INLA's own hyperparameter grid-
# search numerical instability (unrelated to the level-name defect this
# file exercises) on several seeds during development; this size and
# signal was spot-checked stable across five seeds and three repeats.
.space_level_data <- function(seed = 14L) {
  set.seed(seed)
  trial_u <- stats::rnorm(4L, 0, 0.4)
  ntrial_u <- stats::rnorm(3L, 0, 0.3)
  d <- expand.grid(
    NAPPLIED2 = factor(c("nil N", "low N", "high N")),
    TRIAL = factor(seq_len(4L)),
    N_TRIAL = factor(seq_len(3L)),
    rep = seq_len(8L)
  )
  napplied_eff <- c("high N" = 1.2, "low N" = 0.6, "nil N" = 0)
  d$GRAIN_YIELD_THA <- 3 +
    napplied_eff[as.character(d$NAPPLIED2)] +
    trial_u[as.integer(d$TRIAL)] +
    ntrial_u[as.integer(d$N_TRIAL)] +
    stats::rnorm(nrow(d), 0, 0.4)
  d
}


# ---------------------------------------------------------------- #
# Mechanism-level proof: INLA never sees the space (deterministic,   #
# no live engine call, no numerical variability)                    #
# ---------------------------------------------------------------- #

test_that(".inla_call() receives data with every non-syntactic level legalised, never the space", {
  skip_if_no_inla()
  d <- .space_level_data()
  fb <- fb_from_asreml(
    fixed = GRAIN_YIELD_THA ~ NAPPLIED2,
    random = ~ TRIAL + N_TRIAL,
    data = d
  )
  gated <- flexyBayes:::lgm_gate(fb)

  captured <- NULL
  testthat::local_mocked_bindings(
    .inla_call = function(args) {
      captured <<- args
      stop("mocked -- intercepted before the real engine call")
    },
    .package = "flexyBayes"
  )
  tryCatch(
    flexyBayes:::emit_inla(gated, data = d, verbose = FALSE),
    error = function(e) NULL
  )

  expect_false(is.null(captured))
  expect_true(is.factor(captured$data$NAPPLIED2))
  expect_identical(
    sort(levels(captured$data$NAPPLIED2)),
    sort(make.names(levels(d$NAPPLIED2), unique = TRUE))
  )
  # The direct FS-26 claim: no level INLA is handed carries a space --
  # checked the same way the implementation itself decides legality
  # (a letter-led probe prefix, since INLA's own coefficient name
  # concatenates a real, letter-led variable name onto the level; see
  # .inla_legalise_factor_levels()'s banner). NOT bare
  # make.names(lv) == lv: that test is a false positive on a purely
  # numeric level like TRIAL's "1"/"2"/... -- concatenated onto a
  # variable name ("TRIAL1") it is already a valid identifier, so it
  # must NOT be legalised, which the second half of this test asserts
  # directly (a regression guard for exactly this false positive,
  # which once silently diverged the per-row and streamed-aggregated
  # INLA emits' coefficient names -- test-aggregation-equivalence-
  # backend.R is where that surfaced).
  for (v in c("NAPPLIED2", "TRIAL", "N_TRIAL")) {
    lv <- levels(captured$data[[v]])
    probed <- paste0("Q9zZ", lv)
    expect_identical(probed, make.names(probed, unique = TRUE))
  }
  expect_identical(levels(captured$data$TRIAL), levels(d$TRIAL))
  expect_identical(levels(captured$data$N_TRIAL), levels(d$N_TRIAL))
})


# ---------------------------------------------------------------- #
# The model fits, and legalises exactly the model's own factors     #
# ---------------------------------------------------------------- #

test_that("a space-containing factor level fits cleanly on INLA", {
  skip_if_no_inla()
  d <- .space_level_data()
  fit <- flexybayes(
    GRAIN_YIELD_THA ~ NAPPLIED2,
    random = ~ TRIAL + N_TRIAL,
    data = d,
    backend = "inla",
    aggregate = FALSE,
    verbose = FALSE
  )
  expect_s3_class(fit, "flexybayes_inla")
  expect_true(fit$num_check$pass)
})

test_that("the level map legalises exactly the factors the model touches", {
  skip_if_no_inla()
  d <- .space_level_data()
  fit <- flexybayes(
    GRAIN_YIELD_THA ~ NAPPLIED2,
    random = ~ TRIAL + N_TRIAL,
    data = d,
    backend = "inla",
    aggregate = FALSE,
    verbose = FALSE
  )
  expect_setequal(names(fit$level_labels), c("NAPPLIED2", "TRIAL", "N_TRIAL"))
  # NAPPLIED2's space-containing levels are the only ones this fixture
  # actually needs legalised; TRIAL / N_TRIAL's purely-numeric levels
  # ("1", "2", ...) are already valid once concatenated onto the
  # variable name ("TRIAL1"), so their map entries are the identity --
  # present (every touched factor gets one, unconditionally), but
  # unchanged. Asserting the identity map here is itself the regression
  # guard for the false positive test-aggregation-equivalence-backend.R
  # caught: an earlier version of this legalisation mapped "1" -> "X1"
  # unconditionally, which the two INLA emits (this per-row one, and
  # the untouched streamed/aggregated one) then disagreed on.
  expect_identical(
    unname(fit$level_labels$NAPPLIED2[c("high.N", "low.N", "nil.N")]),
    c("high N", "low N", "nil N")
  )
  expect_identical(unname(fit$level_labels$TRIAL[["1"]]), "1")
  expect_identical(fit$level_labels$TRIAL, stats::setNames(
    levels(.space_level_data()$TRIAL), levels(.space_level_data()$TRIAL)
  ))
})


# ---------------------------------------------------------------- #
# All four accessors round-trip the user's own labels                #
# ---------------------------------------------------------------- #

test_that("summary() shows the user's own labels", {
  skip_if_no_inla()
  d <- .space_level_data()
  fit <- flexybayes(
    GRAIN_YIELD_THA ~ NAPPLIED2,
    random = ~ TRIAL + N_TRIAL,
    data = d,
    backend = "inla",
    aggregate = FALSE,
    verbose = FALSE
  )
  s <- summary(fit)
  expect_setequal(
    s$fixed$term,
    c("(Intercept)", "NAPPLIED2low N", "NAPPLIED2nil N")
  )
  out <- capture.output(print(s))
  expect_true(any(grepl("NAPPLIED2low N", out, fixed = TRUE)))
  expect_false(any(grepl("NAPPLIED2low.N", out, fixed = TRUE)))
})

test_that("coef() shows the user's own labels", {
  skip_if_no_inla()
  d <- .space_level_data()
  fit <- flexybayes(
    GRAIN_YIELD_THA ~ NAPPLIED2,
    random = ~ TRIAL + N_TRIAL,
    data = d,
    backend = "inla",
    aggregate = FALSE,
    verbose = FALSE
  )
  cf <- coef(fit)
  expect_setequal(
    names(cf),
    c("(Intercept)", "NAPPLIED2low N", "NAPPLIED2nil N")
  )
})

test_that("ranef() shows the user's own labels", {
  skip_if_no_inla()
  d <- .space_level_data()
  fit <- flexybayes(
    GRAIN_YIELD_THA ~ NAPPLIED2,
    random = ~ TRIAL + N_TRIAL,
    data = d,
    backend = "inla",
    aggregate = FALSE,
    verbose = FALSE
  )
  re <- ranef(fit)
  expect_setequal(re$TRIAL$level, c("1", "2", "3", "4"))
  expect_setequal(re$N_TRIAL$level, c("1", "2", "3"))
})

test_that("predict(classify = ) shows the user's own labels", {
  skip_if_no_inla()
  skip_if_not_installed("emmeans")
  d <- .space_level_data()
  fit <- flexybayes(
    GRAIN_YIELD_THA ~ NAPPLIED2,
    random = ~ TRIAL + N_TRIAL,
    data = d,
    backend = "inla",
    aggregate = FALSE,
    verbose = FALSE
  )
  pc <- predict(fit, classify = "NAPPLIED2")
  expect_s3_class(pc, "fb_predict_classify")
  expect_setequal(as.character(pc$NAPPLIED2), c("nil N", "low N", "high N"))
  expect_identical(nrow(pc), 3L)
})


# ---------------------------------------------------------------- #
# backend = "auto" now routes to INLA where it silently fell to      #
# brms before the fix (FS-26's "Class" note)                        #
# ---------------------------------------------------------------- #

test_that("backend = 'auto' now routes to INLA on a space-level design (was silently brms)", {
  skip_if_no_inla()
  d <- .space_level_data()
  p <- flexybayes(
    GRAIN_YIELD_THA ~ NAPPLIED2,
    random = ~ TRIAL + N_TRIAL,
    data = d,
    backend = "auto",
    aggregate = FALSE,
    plan = TRUE
  )
  expect_identical(p$backend_chosen, "inla")
  expect_identical(p$gate_outcome, "accept")

  fit <- flexybayes(
    GRAIN_YIELD_THA ~ NAPPLIED2,
    random = ~ TRIAL + N_TRIAL,
    data = d,
    backend = "auto",
    aggregate = FALSE,
    verbose = FALSE
  )
  expect_s3_class(fit, "flexybayes_inla")
})


# ---------------------------------------------------------------- #
# brms control: the identical model fits regardless of the space    #
# (confirms the defect was INLA-path-specific, per FS-26's own       #
# control experiment: "The same data on backend = brms fits either  #
# way, in 23.4 s")                                                   #
# ---------------------------------------------------------------- #
#
# "Confirm brms output labels are also the user's" (C4 item text) is
# a verification step, not a claim both engines format labels
# identically. Grounded here rather than assumed: brms/Stan's OWN
# fixed-effect parameter-naming convention (independent of anything
# this slice touches) strips the space rather than substituting a
# dot -- "NAPPLIED2low N" becomes "NAPPLIED2lowN", not
# "NAPPLIED2low.N" (INLA/make.names()'s convention) and not
# "NAPPLIED2low N" (a literal round-trip). What IS the user's own,
# unmangled, on brms is anything that passes through as a bare DATA
# VALUE rather than a synthesised parameter name -- the random-effect
# group levels (TRIAL, N_TRIAL) ranef() reads, verified below.

test_that("brms fits the identical space-level model without legalisation, and its random-effect levels are the user's own unmangled values", {
  skip_if_no_brms()
  d <- .space_level_data()
  # A label round-trip control (no numeric assertion below), so the
  # tiny deliberate sampler budget's ESS warning is irrelevant to what
  # this test checks -- muffled rather than fixed by a bigger budget
  # (see helper-ess-warnings.R).
  fit <- .muffle_ess_warnings(flexybayes(
    GRAIN_YIELD_THA ~ NAPPLIED2,
    random = ~ TRIAL + N_TRIAL,
    data = d,
    backend = "brms",
    n_samples = 200L,
    warmup = 200L,
    chains = 1L,
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_brms")
  # No level_labels map on a brms fit -- there was nothing to legalise.
  expect_null(fit$level_labels)
  re <- ranef(fit)
  expect_setequal(re$TRIAL$level, c("1", "2", "3", "4"))
  expect_setequal(re$N_TRIAL$level, c("1", "2", "3"))
})
