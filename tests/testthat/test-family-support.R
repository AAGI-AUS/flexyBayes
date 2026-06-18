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
# dispatch path, before any backend code runs (so neither greta nor
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
