# =============================================================================
# The greta quarantine is an invariant, not a property of one code path.
#
# "greta is quarantined" was true of the backend registry and false of
# fb_dirichlet(), which called greta directly and initialised TensorFlow while
# `backend = "greta"` refused two lines away. A quarantine that holds at the
# entry points someone remembered to check is not a quarantine.
#
# These tests are therefore written to fail when a NEW fitting route is added
# without the guard, rather than to enumerate the routes known today. The
# source scan is the part that does that work; the behavioural tests pin the
# routes already known.
#
# What is deliberately NOT quarantined: reading a greta model that was fitted
# elsewhere. fb_from_greta() imports an existing object's draws, and
# fb_log_posterior() operates on one. Neither fits anything, and dropping them
# would strand saved work. Re-entry for the fitting side is repair plus
# conformance, so the emit sources stay in the tree.
# =============================================================================

test_that("greta and gretaR are quarantined in the registry", {
  for (b in c("greta", "gretaR")) {
    expect_true(flexyBayes:::.backend_is_quarantined(b), label = b)
  }
  for (b in c("inla", "brms")) {
    expect_false(flexyBayes:::.backend_is_quarantined(b), label = b)
  }
})

test_that("every known fitting entry point refuses greta", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(1L)
  d <- data.frame(y = stats::rnorm(20L), x = stats::rnorm(20L))

  # The mixed-model surface.
  expect_error(
    suppressMessages(flexybayes(y ~ x, data = d, backend = "greta")),
    class = "flexybayes_refusal_backend_quarantined"
  )
  expect_error(
    suppressMessages(flexybayes(y ~ x, data = d, backend = "gretaR")),
    class = "flexybayes_refusal_backend_quarantined"
  )

  # The specialist estimator that used to bypass the registry entirely.
  m <- matrix(stats::rgamma(60L, 2), 20L, 3L)
  m <- m / rowSums(m)
  expect_error(
    fb_dirichlet(m, method = "greta"),
    class = "flexybayes_refusal_backend_quarantined"
  )
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

test_that("no NEW fitting route reaches greta without the quarantine guard", {
  # Source scan, so that adding a fitting route which calls greta fails here
  # rather than being found by an audit months later. Files that hold the
  # retained emitters and the import path are exempt by name: they are
  # unreachable while the registry says quarantined, and that unreachability
  # is what the behavioural tests above pin.
  skip_on_cran()
  # Sources are present under load_all() and during development, and absent
  # from an installed package, where R/ has been byte-compiled away. Skip
  # rather than reach for a root-finder: the behavioural tests above cover
  # the routes that exist, and this scan exists to catch a route added
  # during development, which is exactly when the sources are there.
  r_dir <- system.file("..", "R", package = "flexyBayes")
  skip_if_not(
    length(r_dir) == 1L && nzchar(r_dir) && dir.exists(r_dir),
    "package sources not available (installed package)"
  )

  exempt <- c(
    "emit_greta.R", "emit_gretaR.R", "codegen.R",
    "fb_from_greta.R", "fb_greta.R", "gretaR_slot.R",
    # Operates on an already-fitted greta object; imports, does not fit.
    "fb_log_posterior.R",
    # Holds the quarantine machinery and its message text.
    "backend_registry.R", "refusal_taxonomy.R", "dispatch.R",
    "backend_status.R", "flexybayes.R", "emit_gaussian_aggregated.R",
    "fb_dirichlet.R"
  )
  files <- setdiff(basename(list.files(r_dir, pattern = "\\.R$")), exempt)
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
      "These files call greta:: outside the retained set. If a new fitting ",
      "route was added, guard it with .backend_is_quarantined(\"greta\"); ",
      "if it only reads an existing fit, add it to `exempt` with a reason."
    )
  )
})
