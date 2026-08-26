# test-plan-tells-truth.R -- C3 / S6 + FS-22: the plan surface tells the
# truth.
#
# Two independent defects, one file:
#
# S6 -- print.fb_plan() and print.fb_aggregation_plan() reported an
# eligibility VERDICT ("not eligible (compression_unproductive)") with no
# numbers behind it. A user could not see the rows-per-cell that produced
# the verdict, or the threshold the verdict was tested against, without
# reading the source (SCALE_STRATEGY_2026-08-22's own contrast: "1.06
# rows per cell against the shipped demo's 1,200 cells").
#
# FS-22 -- fb_plan() is a brms-formula-only entry point
# (fb_from_brms(formula = formula, ...)) that captured `random =` /
# `residual =` into its unused `...` and silently dropped them. A caller
# who reached for fb_plan() the way they would reach for flexybayes()
# ASReml-style got a plan for a FIXED-EFFECTS-ONLY model: no random rows,
# no residual row, and a backend prediction (INLA, accept) that live
# dispatch on the real (fuller) model does not honour.
#
# Reproduced against the pre-fix source (see the induced-defect section
# of the report) exactly as filed in
# review/phase_reports_091/field_sweep/FIELD_SWEEP_MISSES.md #FS-22.

skip_if_no_inla <- function() skip_if_not_installed("INLA")
skip_if_no_brms <- function() skip_if_not_installed("brms")


# ---------------------------------------------------------------- #
# S6 -- rows-per-cell + threshold-stated verdict                    #
# ---------------------------------------------------------------- #

test_that("print.fb_aggregation_plan() states N, K, rows/cell and the threshold -- eligible case", {
  set.seed(1L)
  n <- 1000L
  df <- data.frame(
    y = rnorm(n),
    f = factor(sample(letters[1:5], n, replace = TRUE)),
    g = factor(sample.int(20L, n, replace = TRUE))
  )
  fb <- fb_from_brms(y ~ f + (1 | g), data = df)
  plan <- .fb_aggregation_plan(fb, .fb_dataset(df))
  out <- capture.output(print(plan))

  # K = 100, N = 1000 -> rows/cell = 10.00, well above the 2.0 threshold.
  expect_true(any(grepl("rows/cell = 10.00", out, fixed = TRUE)))
  expect_true(any(grepl("threshold: rows/cell >= 2.0", out, fixed = TRUE)))
  expect_true(any(grepl("aggregation will pay", out, fixed = TRUE)))
})

test_that("print.fb_aggregation_plan() states the numbers even when NOT eligible (compression_unproductive)", {
  # Hand-built IR, K = 80 on N = 100 -> rows/cell = 1.25, below threshold.
  # Before the fix, the "not eligible" branch printed only the reason
  # code -- no N, K, or rows/cell -- so a user could not see WHY.
  fb <- structure(
    list(
      response = "y",
      family = "gaussian",
      link = "identity",
      intercept = TRUE,
      fixed_terms = list(),
      random_terms = list(list(type = "simple", var = "g", var_n = 80L)),
      residual_terms = list(),
      addition_terms = list(),
      priors = list(),
      data_summary = list(n = 100L),
      capabilities = character(),
      source = "test"
    ),
    class = c("fb_terms", "list")
  )
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 100L,
    col_types = list(y = "double", g = "factor"),
    dictionaries = list(g = as.character(seq_len(80L)))
  )
  plan <- .fb_aggregation_plan(fb, ds)
  expect_false(plan$eligible)
  expect_identical(plan$reason_codes, "compression_unproductive")

  out <- capture.output(print(plan))
  expect_true(any(grepl("N = 100; K = 80; rows/cell = 1.25", out, fixed = TRUE)))
  expect_true(any(grepl("threshold: rows/cell >= 2.0", out, fixed = TRUE)))
  expect_true(any(grepl("will NOT pay", out, fixed = TRUE)))
})

test_that("print.fb_plan()'s Aggregation line carries the same N/K/rows-per-cell/verdict", {
  skip_if_no_inla()
  set.seed(2L)
  n <- 1000L
  df <- data.frame(
    y = rnorm(n),
    f = factor(sample(letters[1:5], n, replace = TRUE)),
    g = factor(sample.int(20L, n, replace = TRUE))
  )
  p <- fb_plan(y ~ f + (1 | g), data = df, backend = "auto")
  out <- capture.output(print(p))
  expect_true(any(grepl("rows/cell = 10.00", out, fixed = TRUE)))
  expect_true(any(grepl("threshold: rows/cell >= 2.0", out, fixed = TRUE)))
  expect_true(any(grepl("aggregation will pay", out, fixed = TRUE)))
})

test_that("print.fb_plan()'s Aggregation line shows numbers on an unproductive-but-known-K model", {
  skip_if_no_inla()
  set.seed(3L)
  # 90 levels on 100 rows -> K/N close to 1, compression_unproductive,
  # but K IS resolvable, so the numbers should print.
  df <- data.frame(
    y = rnorm(100L),
    f = factor(sample(paste0("f", seq_len(90L)), 100L, replace = TRUE))
  )
  p <- fb_plan(y ~ f, data = df, backend = "auto")
  expect_false(isTRUE(p$aggregation$eligible))
  expect_false(is.na(p$aggregation$K))
  out <- capture.output(print(p))
  expect_true(any(grepl("N = 100; K = ", out, fixed = TRUE)))
  expect_true(any(grepl("will NOT pay", out, fixed = TRUE)))
})


# ---------------------------------------------------------------- #
# FS-22 -- fb_plan() must run the same structural checks live      #
# dispatch runs, on both spellings                                  #
# ---------------------------------------------------------------- #

.fs22_data <- function() {
  set.seed(20260822L)
  loc <- factor(paste0("L", 1:3))
  yearf <- factor(paste0("Y", 1:2))
  gen <- factor(paste0("G", 1:5))
  rep <- factor(1:2)
  d <- expand.grid(loc = loc, yearf = yearf, gen = gen, rep = rep)
  d$env <- interaction(d$loc, d$yearf, drop = TRUE)
  d$yield <- rnorm(nrow(d))
  d
}

test_that("FS-22 reproduction: fb_plan() names brms under backend = 'auto', matching flexybayes(plan = TRUE)", {
  skip_if_no_brms()
  d <- .fs22_data()

  p_fb_plan <- fb_plan(
    yield ~ loc + yearf + loc:yearf,
    random = ~ gen + rep:loc:yearf + gen:yearf + gen:loc + gen:loc:yearf,
    residual = ~ dsum(~ units | env),
    data = d,
    backend = "auto"
  )
  p_flexybayes_plan <- suppressMessages(flexybayes(
    fixed = yield ~ loc + yearf + loc:yearf,
    random = ~ gen + rep:loc:yearf + gen:yearf + gen:loc + gen:loc:yearf,
    residual = ~ dsum(~ units | env),
    data = d,
    backend = "auto",
    plan = TRUE,
    verbose = FALSE
  ))

  # "Both spellings" of the entry point (fb_plan() vs. flexybayes(plan =
  # TRUE)) predict the identical route.
  expect_identical(p_fb_plan$backend_chosen, p_flexybayes_plan$backend_chosen)
  expect_identical(p_fb_plan$gate_outcome, p_flexybayes_plan$gate_outcome)
  expect_identical(p_fb_plan$will_fit, p_flexybayes_plan$will_fit)

  # And the FS-22 headline: brms, not INLA/aggregated_inla, with the
  # gate's real verdict reported (not silently accepted on a
  # fixed-effects-only misreading of the model).
  expect_identical(p_fb_plan$backend_chosen, "brms")
  expect_false(identical(p_fb_plan$path, "aggregated_inla"))
  expect_identical(p_fb_plan$gate_outcome, "refuse_structural")

  # The representation table enumerates random AND residual terms, not
  # only the three fixed terms (FS-22's second observed symptom).
  ids <- vapply(p_fb_plan$representation_plan, function(rp) rp$term_id, character(1L))
  expect_true(any(grepl("^\\(", ids))) # at least one random-term row, "(... | ...)"
  expect_true(any(grepl("dsum(~ units | env)", ids, fixed = TRUE)))
})

test_that("FS-22 reproduction: fb_plan() matches the live refusal under backend = 'inla' (both spellings of backend)", {
  skip_if_no_inla()
  d <- .fs22_data()

  p_inla <- fb_plan(
    yield ~ loc + yearf + loc:yearf,
    random = ~ gen + rep:loc:yearf + gen:yearf + gen:loc + gen:loc:yearf,
    residual = ~ dsum(~ units | env),
    data = d,
    backend = "inla"
  )
  expect_false(p_inla$will_fit)
  expect_true(is.na(p_inla$backend_chosen))

  live_err <- tryCatch(
    suppressMessages(flexybayes(
      fixed = yield ~ loc + yearf + loc:yearf,
      random = ~ gen + rep:loc:yearf + gen:yearf + gen:loc + gen:loc:yearf,
      residual = ~ dsum(~ units | env),
      data = d,
      backend = "inla",
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(live_err, "flexybayes_refusal_inla_gate_refused")
})

test_that("FS-22 reproduction: fb_plan() predicts the live auto fit's backend (plan == live)", {
  skip_if_no_brms()
  skip_on_cran()
  d <- .fs22_data()

  p <- fb_plan(
    yield ~ loc + yearf + loc:yearf,
    random = ~ gen + rep:loc:yearf + gen:yearf + gen:loc + gen:loc:yearf,
    residual = ~ dsum(~ units | env),
    data = d,
    backend = "auto"
  )
  # This assertion is routing parity (which backend the fit landed on),
  # not a numeric agreement -- Bulk/Tail-ESS is irrelevant to it, so the
  # tiny deliberate n_samples = warmup = 200 / chains = 1 budget's ESS
  # warning is muffled rather than fixed by raising the budget (see
  # helper-ess-warnings.R).
  fit <- .muffle_ess_warnings(suppressMessages(flexybayes(
    fixed = yield ~ loc + yearf + loc:yearf,
    random = ~ gen + rep:loc:yearf + gen:yearf + gen:loc + gen:loc:yearf,
    residual = ~ dsum(~ units | env),
    data = d,
    backend = "auto",
    n_samples = 200L,
    warmup = 200L,
    chains = 1L,
    verbose = FALSE
  )))
  expect_identical(fit$extras$backend_decision$backend, p$backend_chosen)
  expect_identical(fit$extras$backend_decision$backend, "brms")
})


# ---------------------------------------------------------------- #
# Parity roster -- plan$backend == fit$backend, or plan refuses      #
# where fit refuses, across a small spread of gate outcomes          #
# ---------------------------------------------------------------- #

.plan_parity_roster <- list(
  list(
    label = "simple random intercept -> INLA accept",
    data = function() {
      set.seed(11L)
      n <- 150L
      data.frame(
        y = rnorm(n),
        x = rnorm(n),
        g = factor(sample(letters[1:6], n, replace = TRUE))
      )
    },
    formula = y ~ x + (1 | g),
    expect_backend = "inla"
  ),
  list(
    label = "factor + factor -> aggregation-eligible aggregated INLA",
    data = function() {
      set.seed(12L)
      n <- 1000L
      data.frame(
        y = rnorm(n),
        f = factor(sample(paste0("f", seq_len(5L)), n, replace = TRUE)),
        g = factor(sample(paste0("g", seq_len(20L)), n, replace = TRUE))
      )
    },
    formula = y ~ f + (1 | g),
    expect_backend = "inla"
  ),
  list(
    label = "colon-interaction random term -> INLA structural refusal -> brms",
    data = function() {
      set.seed(13L)
      n <- 120L
      data.frame(
        y = rnorm(n),
        x = rnorm(n),
        a = factor(sample(letters[1:3], n, replace = TRUE)),
        b = factor(sample(LETTERS[1:3], n, replace = TRUE)),
        cc = factor(sample(1:3, n, replace = TRUE))
      )
    },
    formula = y ~ x + (1 | a:b:cc),
    expect_backend = "brms"
  )
)

test_that("parity roster: plan$backend_chosen equals the live fit's backend across a spread of gate outcomes", {
  skip_if_no_inla()
  skip_if_no_brms()
  skip_on_cran()

  for (spec in .plan_parity_roster) {
    d <- spec$data()
    p <- fb_plan(spec$formula, data = d, backend = "auto")
    expect_identical(
      p$backend_chosen,
      spec$expect_backend,
      info = paste0(spec$label, " (plan)")
    )
    # Routing parity again (see the FS-22 reproduction above) -- the
    # tiny sampler budget's ESS warning is muffled, not fixed by a
    # bigger budget, for the same reason.
    fit <- .muffle_ess_warnings(suppressMessages(flexybayes(
      spec$formula,
      data = d,
      backend = "auto",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE
    )))
    expect_identical(
      fit$extras$backend_decision$backend,
      p$backend_chosen,
      info = paste0(spec$label, " (live parity)")
    )
  }
})
