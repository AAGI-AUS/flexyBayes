# test-plan-tells-truth.R: the plan surface tells the truth.
#
# Two independent defects, one file:
#
# First, print.fb_plan() and print.fb_aggregation_plan() reported an
# eligibility VERDICT ("not eligible (compression_unproductive)") with no
# numbers behind it. A user could not see the rows-per-cell that produced
# the verdict, or the threshold the verdict was tested against, without
# reading the source (for example "1.06 rows per cell against the shipped demo's
# 1,200 cells").
#
# Second, fb_plan() is a brms-formula-only entry point
# (fb_from_brms(formula = formula, ...)) that captured `random =` /
# `residual =` into its unused `...` and silently dropped them. A caller
# who reached for fb_plan() the way they would reach for flexybayes()
# ASReml-style got a plan for a FIXED-EFFECTS-ONLY model: no random rows,
# no residual row, and a backend prediction (INLA, accept) that live
# dispatch on the real (fuller) model does not honour.
#
# Both were reproduced against the pre-fix source.

skip_if_no_inla <- function() skip_if_not_installed("INLA")
skip_if_no_brms <- function() skip_if_not_installed("brms")


# ---------------------------------------------------------------- #
# rows-per-cell + threshold-stated verdict                        #
# ---------------------------------------------------------------- #

test_that("print.fb_aggregation_plan() states N, K, rows/cell and the threshold, eligible case", {
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
  # code, no N, K, or rows/cell, so a user could not see WHY.
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
# fb_plan() must run the same structural checks live      #
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

test_that("the plan and the gate agree under backend = 'auto'", {
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

  # The reported defect was that the plan announced an INLA route while the
  # gate refused the model, and the live call went to brms. The five
  # strata and the sectioned residual all reach INLA from 0.10.0.9000, so
  # the gate accepts and INLA is the correct answer, the plan and the
  # gate agree, which is what the defect was about. What must still hold is
  # that the plan does not name the AGGREGATED route: a sectioned
  # residual is not aggregatable whatever the gate says.
  expect_identical(p_fb_plan$backend_chosen, "inla")
  expect_false(identical(p_fb_plan$path, "aggregated_inla"))
  expect_identical(p_fb_plan$gate_outcome, "accept")

  # The representation table enumerates random AND residual terms, not
  # only the three fixed terms (the second observed symptom).
  ids <- vapply(p_fb_plan$representation_plan, function(rp) rp$term_id, character(1L))
  expect_true(any(grepl("^\\(", ids))) # at least one random-term row, "(... | ...)"
  expect_true(any(grepl("dsum(~ units | env)", ids, fixed = TRUE)))
})

test_that("fb_plan() matches the live fit under backend = 'inla'", {
  skip_if_no_inla()
  d <- .fs22_data()

  p_inla <- fb_plan(
    yield ~ loc + yearf + loc:yearf,
    random = ~ gen + rep:loc:yearf + gen:yearf + gen:loc + gen:loc:yearf,
    residual = ~ dsum(~ units | env),
    data = d,
    backend = "inla"
  )
  expect_true(p_inla$will_fit)
  expect_identical(p_inla$backend_chosen, "inla")
})

test_that("plan and live agree when INLA genuinely refuses", {
  # The parity contract has to be shown on a refusal as well as on a fit,
  # or it only proves the easy half. The five-stratum MET model used to
  # be the refusal case and is now the fit case, so the refusal case
  # moves to a family INLA's likelihood roster has no counterpart for.
  skip_if_no_inla()
  d <- .fs22_data()
  d$yield <- pmax(d$yield, 0.01)

  p <- fb_plan(
    yield ~ loc,
    random = ~gen,
    data = d,
    family = "hurdle_gamma",
    backend = "inla"
  )
  expect_false(p$will_fit)

  live_err <- tryCatch(
    suppressMessages(flexybayes(
      fixed = yield ~ loc,
      random = ~gen,
      data = d,
      family = "hurdle_gamma",
      backend = "inla",
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(live_err, "error")
})

test_that("fb_plan() predicts the live auto fit's backend (plan == live)", {
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
  # not a numeric agreement. Bulk/Tail-ESS is irrelevant to it, so the
  # tiny deliberate n_samples = warmup = 200 / chains = 1 budget's ESS
  # warning is muffled, not fixed by raising the budget (see
  # helper-ess-warnings.R). The yield is pure noise, so a variance
  # component is weakly identified when the fit lands on INLA.
  withr::local_options(flexyBayes.silence_weak_identification_warning = TRUE)
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
  # The parity assertion is the first one: what the plan named is where
  # the fit landed. The second pins which engine that is today.
  expect_identical(fit$extras$backend_decision$backend, p$backend_chosen)
  expect_identical(fit$extras$backend_decision$backend, "inla")
})


# ---------------------------------------------------------------- #
# Parity roster, plan$backend == fit$backend, or plan refuses      #
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
    # This entry used to read "INLA structural refusal -> brms". The
    # three-way interaction random term is an iid field over the
    # interaction index on INLA from 0.10.0.9000, so it stays on INLA.
    label = "colon-interaction random term -> INLA accept",
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
    expect_backend = "inla"
  ),
  list(
    # The roster needs a route that still leaves INLA, or it stops
    # testing the branch it was built for. A hurdle gamma has no
    # counterpart in INLA's likelihood roster, so the family gate refuses
    # and `auto` routes to brms.
    label = "hurdle gamma -> INLA family refusal -> brms",
    data = function() {
      set.seed(14L)
      n <- 120L
      data.frame(
        y = pmax(rnorm(n, 2, 1), 0),
        x = rnorm(n),
        g = factor(sample(letters[1:6], n, replace = TRUE))
      )
    },
    formula = y ~ x + (1 | g),
    family = "hurdle_gamma",
    expect_backend = "brms"
  )
)

test_that("parity roster: plan$backend_chosen equals the live fit's backend across a spread of gate outcomes", {
  skip_if_no_inla()
  skip_if_no_brms()
  skip_on_cran()
  skip_if_not_long_tier()

  for (spec in .plan_parity_roster) {
    d <- spec$data()
    args <- list(spec$formula, data = d, backend = "auto")
    if (!is.null(spec$family)) args$family <- spec$family
    p <- do.call(fb_plan, args)
    expect_identical(
      p$backend_chosen,
      spec$expect_backend,
      info = paste0(spec$label, " (plan)")
    )
    # Routing parity again (see the plan tests above), the
    # tiny sampler budget's ESS warning is muffled, not fixed by a
    # bigger budget, for the same reason.
    fit_args <- list(
      spec$formula,
      data = d,
      backend = "auto",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE
    )
    if (!is.null(spec$family)) fit_args$family <- spec$family
    fit <- .muffle_ess_warnings(suppressMessages(
      do.call(flexybayes, fit_args)
    ))
    expect_identical(
      fit$extras$backend_decision$backend,
      p$backend_chosen,
      info = paste0(spec$label, " (live parity)")
    )
  }
})
