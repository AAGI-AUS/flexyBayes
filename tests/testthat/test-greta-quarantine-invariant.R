# =============================================================================
# 0.9.3: the native engine formerly named here was withdrawn entirely, not
# quarantined -- no code path, export, registry row, or `Suggests` entry
# remains (see NEWS.md). The invariant this file pins therefore changed
# shape: it used to be "every fitting entry point refuses a quarantined
# name"; it is now "the name is unknown to every fitting entry point, and
# the engine it used to name has no trace left in R/ at all".
#
# These tests are written to fail when a NEW fitting route is added without
# going through the shared unknown-backend check, rather than to enumerate
# the routes known today. The source scan is the part that does that work;
# the behavioural tests pin the routes already known.
# =============================================================================

test_that("naming the withdrawn engine (or its dormant sibling) raises unknown_backend at every known fitting entry point", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(1L)
  d <- data.frame(y = stats::rnorm(20L), x = stats::rnorm(20L))

  # The mixed-model surface. Both names raise the same reason code and both
  # active engines are named in the message.
  for (nm in c("greta", "gretaR")) {
    err <- tryCatch(
      suppressMessages(flexybayes(y ~ x, data = d, backend = nm)),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_unknown_backend_refusal")
    expect_identical(err$reason_code, "unknown_backend", label = nm)
    expect_match(conditionMessage(err), '"inla"', fixed = TRUE, label = nm)
    expect_match(conditionMessage(err), '"brms"', fixed = TRUE, label = nm)
  }

  # The specialist estimator that, before 0.9.3, called the native engine
  # directly rather than going through the backend registry: fb_dirichlet()
  # no longer offers a method by that name at all, so the request fails
  # `match.arg()` before any fitting code runs.
  m <- matrix(stats::rgamma(60L, 2), 20L, 3L)
  m <- m / rowSums(m)
  err <- tryCatch(fb_dirichlet(m, method = "greta"), error = function(e) e)
  expect_s3_class(err, "error")
  expect_false(inherits(err, "flexybayes_refusal"))
  expect_match(conditionMessage(err), "ml", fixed = TRUE)
})

test_that("the maximum-likelihood Dirichlet fit is unaffected", {
  set.seed(2L)
  m <- matrix(stats::rgamma(90L, 2), 30L, 3L)
  m <- m / rowSums(m)
  fit <- fb_dirichlet(m)
  expect_s3_class(fit, "fb_dirichlet_fit")
  expect_identical(fit$method, "ml")
  expect_equal(sum(fit$mean_composition), 1, tolerance = 1e-8)
})

test_that("no R source calls the withdrawn engine's package at all", {
  # 0.9.3 strengthens the old exempt-list source scan (which allowed the
  # retained emitters to call greta::) into an unconditional invariant:
  # nothing under R/ may call greta:: any more, because nothing that ever
  # did still exists. This is the source-level twin of the goal-grep
  # acceptance criterion (WP-G spec) -- it runs from the installed
  # package's perspective rather than the source tree's.
  skip_on_cran()
  # Sources are present under load_all() and during development, and
  # absent from an installed package, where R/ has been byte-compiled
  # away. Skip rather than reach for a root-finder: this scan exists to
  # catch a route added during development, which is exactly when the
  # sources are there.
  r_dir <- system.file("..", "R", package = "flexyBayes")
  skip_if_not(
    length(r_dir) == 1L && nzchar(r_dir) && dir.exists(r_dir),
    "package sources not available (installed package)"
  )

  files <- list.files(r_dir, pattern = "\\.R$")
  offenders <- character(0)
  for (f in files) {
    lines <- readLines(file.path(r_dir, f), warn = FALSE)
    code <- lines[!grepl("^\\s*#", lines)]
    if (any(grepl("greta::", code, fixed = TRUE))) {
      offenders <- c(offenders, f)
    }
  }
  expect_identical(
    offenders, character(0),
    info = paste0(
      "These files call greta:: after the 0.9.3 withdrawal. The engine ",
      "was removed entirely, not quarantined -- there is no longer an ",
      "exempt set of retained emitters to call it from."
    )
  )
})
