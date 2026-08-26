# ---------------------------------------------------------------- #
# Family support -- .resolve_family() is the single authoritative   #
# family gate every user entry passes through (asreml via fb.R,     #
# brms via fb_from_brms.R). It admits only the families flexyBayes   #
# can emit; any other family -- including those INLA's roster        #
# recognises but flexyBayes cannot emit (survival / time-to-event)   #
# -- is refused up front with a structured, registry-backed refusal, #
# never silently fitted. (Audit 2026-05-30, finding                  #
# D-gap-survival-silent: corrected -- the refusal is clear and       #
# up-front, not a silent malformed model.)                           #
# ---------------------------------------------------------------- #

test_that(".resolve_family() accepts the eight supported families", {
  for (fam in c(
    "gaussian",
    "binomial",
    "binary",
    "poisson",
    "negative_binomial",
    "negbinom",
    "gamma",
    "beta"
  )) {
    fl <- flexyBayes:::.resolve_family(fam, NULL)
    expect_type(fl, "list")
    expect_equal(fl$family, tolower(fam))
  }
})

test_that(".resolve_family() refuses an unsupported family with a structured refusal", {
  err <- expect_error(
    flexyBayes:::.resolve_family("weibullsurv", NULL),
    class = "flexybayes_refusal_unsupported_family"
  )
  expect_match(conditionMessage(err), "Unsupported family") # legacy phrase retained
  expect_match(conditionMessage(err), "survival") # roadmap note
  expect_match(conditionMessage(err), "fb_refusals") # discovery pointer
  # retained family class so any legacy class-based handler keeps working
  expect_s3_class(err, "flexybayes_unsupported_family")
})

test_that("all five survival families are refused at the family gate", {
  for (fam in c(
    "exponentialsurv",
    "weibullsurv",
    "loggaussiansurv",
    "lognormalsurv",
    "coxph"
  )) {
    expect_error(
      flexyBayes:::.resolve_family(fam, NULL),
      class = "flexybayes_refusal_unsupported_family"
    )
  }
})

# Integration -- the refusal fires through the real flexybayes()
# dispatch path, before any backend code runs (so neither brms nor
# INLA need be installed): .resolve_family() is reached inside
# fb_from_asreml() at the top of flexybayes().

test_that("flexybayes() refuses a survival family end-to-end", {
  dat <- data.frame(y = abs(rnorm(20)) + 0.1, x = rnorm(20))
  err <- expect_error(
    flexybayes(
      y ~ x,
      data = dat,
      family = "weibullsurv",
      backend = "inla",
      verbose = FALSE
    ),
    class = "flexybayes_refusal_unsupported_family"
  )
  expect_match(conditionMessage(err), "Unsupported family")
})

# Discovery surface lists the refusal, free of internal-history tokens.

test_that("fb_refusals() lists the unsupported-family refusal cleanly", {
  refs <- fb_refusals()
  expect_true("unsupported_family" %in% refs$reason_code)
  row <- refs[refs$reason_code == "unsupported_family", ]
  expect_false(any(grepl(
    "ADR|Stage|Wave|Phase",
    unlist(row, use.names = FALSE)
  )))
})


# ---------------------------------------------------------------- #
# The allowlist against the engines' live rosters                   #
# (field-sweep FS-4 / field finding C2)                             #
# ---------------------------------------------------------------- #
#
# The entry allowlist is narrower than either engine's family roster,
# deliberately -- a family is admitted when this package has an emit for
# it and a test exercises it. What it must not do is refuse a family the
# engine behind it carries with no documented boundary and no named
# alternative. `hurdle_gamma` was that case: brms-native, refused at the
# door, absent from every documented boundary in the tree.
#
# The rosters are read from the installed engines rather than recalled.
# The first draft of this work asserted that brms carries Tweedie
# natively; a live `brms::brmsfamily()` call said otherwise, and getting
# it backwards would have changed what the fix is.

test_that("hurdle_gamma is admitted, and brms declares it natively", {
  skip_if_not_installed("brms")
  fl <- flexyBayes:::.resolve_family("hurdle_gamma", NULL)
  expect_equal(fl$family, "hurdle_gamma")
  expect_equal(fl$link, "log")
  bf <- brms::brmsfamily("hurdle_gamma")
  expect_setequal(bf$dpars, c("mu", "shape", "hu"))
  # The flexyBayes emit hands brms that same family object.
  mapped <- flexyBayes:::.fb_family_to_brms("hurdle_gamma", "log")
  expect_equal(mapped$family, "hurdle_gamma")
})

test_that("hurdle_gamma carries no residual sigma, so no sigma row is emitted", {
  # The family-traits table (R/family_traits.R) is the one source; brms
  # declares mu / shape / hu, none of which is a standard deviation.
  expect_false(flexyBayes:::.fb_family_has_brms_sigma("hurdle_gamma"))
})

test_that("hurdle_gamma fits end-to-end on brms at a reachability budget", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(20260821L)
  n <- 60L
  g <- factor(rep(seq_len(10L), each = 6L))
  x <- stats::rnorm(n)
  eta <- 0.5 + 0.4 * x + stats::rnorm(10L, 0, 0.7)[as.integer(g)]
  d <- data.frame(
    y = ifelse(
      stats::rbinom(n, 1L, 0.7) == 1L,
      stats::rgamma(n, shape = 2, rate = 2 / exp(eta)),
      0
    ),
    x = x,
    g = g
  )
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ x,
    random = ~g,
    data = d,
    family = "hurdle_gamma",
    backend = "brms",
    chains = 1L,
    n_samples = 300L,
    warmup = 150L,
    seed = 1L,
    verbose = FALSE
  )))
  expect_s3_class(fit, "flexybayes")
  expect_equal(fit$brms$family$family, "hurdle_gamma")
  ps <- as.data.frame(brms::prior_summary(fit$brms))
  # The default prior reached the variance component, and no sigma row was
  # sent to a family that has none.
  expect_true(any(ps$source == "user" & ps$class == "sd"))
  expect_false("sigma" %in% ps$class)
})

test_that("hurdle_gamma refuses on INLA at the family gate, not in the engine", {
  skip_if_not_installed("INLA")
  d <- data.frame(
    y = c(rep(0, 10L), abs(stats::rnorm(50L)) + 0.1),
    x = stats::rnorm(60L),
    g = factor(rep(seq_len(10L), each = 6L))
  )
  err <- expect_error(
    suppressMessages(flexybayes(
      y ~ x,
      random = ~g,
      data = d,
      family = "hurdle_gamma",
      backend = "inla",
      verbose = FALSE
    )),
    class = "flexybayes_refusal_inla_gate_refused"
  )
  expect_match(conditionMessage(err), "family_allowlist")
  # And INLA's live roster really has no counterpart.
  expect_false("hurdle_gamma" %in% names(INLA::inla.models()$likelihood))
})

test_that("the three remaining C2 families are documented boundaries", {
  skip_if_not_installed("brms")
  # brms's own constructor is the oracle for the first column.
  for (fam in c("zero_inflated_gamma", "tweedie", "compound_poisson")) {
    expect_error(brms::brmsfamily(fam), label = fam)
  }
  for (fam in c("zero_inflated_gamma", "tweedie", "compound_poisson")) {
    err <- expect_error(
      flexyBayes:::.resolve_family(fam, NULL),
      class = "flexybayes_refusal_unsupported_family"
    )
    msg <- conditionMessage(err)
    expect_match(msg, "Boundary note", label = fam)
    # Every one of the three names an implemented alternative.
    expect_match(msg, "hurdle_gamma", label = fam)
  }
})

test_that("the tweedie boundary is attributed to flexyBayes, not to INLA", {
  skip_if_not_installed("INLA")
  # INLA does carry the likelihood, so the message must not imply the
  # engine lacks it -- that is the error the boundary note exists to stop.
  expect_true("tweedie" %in% names(INLA::inla.models()$likelihood))
  err <- expect_error(
    flexyBayes:::.resolve_family("tweedie", NULL),
    class = "flexybayes_refusal_unsupported_family"
  )
  expect_match(
    conditionMessage(err),
    "flexyBayes boundary rather than an engine one"
  )
})

test_that("the documented boundaries are recorded in KNOWN_ISSUES", {
  path <- system.file("KNOWN_ISSUES.md", package = "flexyBayes")
  skip_if(!nzchar(path))
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  for (fam in c(
    "hurdle_gamma",
    "zero_inflated_gamma",
    "tweedie",
    "compound_poisson"
  )) {
    expect_match(txt, fam, fixed = TRUE, label = fam)
  }
})
