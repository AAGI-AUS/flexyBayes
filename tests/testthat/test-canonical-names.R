# Tests for the canonical parameter-name registry (ADR 0005).
#
# Contract:
#   - canonical_names(fit) returns a list(map, transform, source, unmapped)
#     resolving the backend-native parameter names to canonical
#     (brms-style) names; transforms carry value-side conversions
#     such as INLA's sqrt(1/prec) precision-to-SD mapping.
#   - The brms mapper translates raw Stan parameter names: b_Intercept ->
#     (Intercept); b_<level> -> <level>; sd_<group>__Intercept ->
#     sd_<group>; sigma stays sigma (already canonical).
#   - The INLA mapper translates summary.fixed (identity) +
#     summary.hyperpar (Precision for ... -> sd_<group> / sigma)
#     with the precision-to-SD transform attached.
#   - The draws accessors (posterior::as_draws_df()) resolve both
#     engines' names through the registry, with no user-supplied map.
#   - register_canonical_mapper() round-trips a custom mapper.
#
# A prior contract clause -- the withdrawn native engine's pin
# accepting a user-supplied canonical_names argument that won over a
# verbatim fallback -- is removed at 0.9.3, not ported (see NEWS.md):
# it was specific to hand-built native model graphs, which needed a
# user-supplied map because there was no formula to parse names from.
# That capability has no successor -- fb_from_brms() / fb_from_asreml()
# both parse a real formula and resolve canonical names from the
# registry automatically, so no ingest path needs (or accepts) a
# user-supplied canonical_names argument today.


mk_canon_data <- function() {
  set.seed(20260522L)
  n <- 30L
  data.frame(
    yield = rnorm(n, 100, 10),
    env = factor(rep(1:3, length.out = n)),
    geno = factor(rep(1:5, length.out = n))
  )
}

# ---------------------------------------------------------------- #
# (a) [removed] canonical mapping for the withdrawn native engine   #
# ---------------------------------------------------------------- #
#
# Deleted (not rewritten): the withdrawn native engine's raw-name
# mapping (mu_atg, sigma_e_atg, tau_<tag>[i,1], ...) has no successor
# -- there is no active engine with that naming convention to map from
# (see NEWS.md, 0.9.3). The equivalent registry coverage for the two
# active engines is subtests (b) (INLA) below and the brms mapping
# exercised in (c)/(d).


# ---------------------------------------------------------------- #
# (b) INLA canonical mapping + sqrt(1/prec) transform               #
# ---------------------------------------------------------------- #

test_that("canonical_names() on an INLA fit translates Precision-for and attaches sqrt(1/prec)", {
  testthat::skip_if_not_installed("INLA")
  d <- mk_canon_data()
  # aggregate = FALSE: this test exercises the per-row INLA fit's
  # hyperparameter naming; the Stage 3A aggregated path (ADR 0022)
  # uses a different fit-object shape (no `$inla$summary.hyperpar`
  # of the per-row form). Aggregated-path canonical-name coverage
  # ships with Phase C tests.
  fit <- suppressMessages(flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    backend = "inla",
    aggregate = FALSE,
    verbose = FALSE
  ))
  cn <- canonical_names(fit)

  expect_identical(cn$source, "registry")

  # Fixed-effect rownames are identity.
  expect_identical(cn$map[["(Intercept)"]], "(Intercept)")
  expect_identical(cn$map[["env2"]], "env2")

  # Hyperparameter rename + transform.
  expect_identical(
    cn$map[["Precision for the Gaussian observations"]],
    "sigma"
  )
  expect_identical(cn$map[["Precision for geno"]], "sd_geno")

  # The transform is keyed by the *native* hyperpar name, not the
  # canonical target. .fb_canonical_draws() applies transforms to a fit's
  # draws while they still carry native names (before renaming), so a
  # canonical key would silently no-op and leave INLA precision draws
  # un-converted. Guard the native keying explicitly.
  prec_g <- "Precision for the Gaussian observations"
  expect_true(is.function(cn$transform[[prec_g]]))
  expect_true(is.function(cn$transform[["Precision for geno"]]))
  expect_null(cn$transform[["sigma"]]) # NOT keyed by canonical name
  expect_null(cn$transform[["sd_geno"]])

  # sqrt(1/4) = 0.5 -- the transform is the precision-to-SD form.
  expect_equal(cn$transform[[prec_g]](4), 0.5)
  expect_equal(cn$transform[["Precision for geno"]](100), 0.1)
})


# ---------------------------------------------------------------- #
# (c) brms draws resolve to canonical names with no user map        #
# ---------------------------------------------------------------- #
# A brms fit's draws should resolve to the canonical parameters via
# the registry without a user-supplied map. Cross-engine (brms-vs-INLA)
# resolution lives behind the skip_if_not_installed("INLA") guard in (d).

test_that("as_draws_df() resolves a brms fit's parameters via the registry", {
  # Fits on brms: a Stan compile, so it is out of the
  # devtools::check() path and runs under devtools::test().
  skip_on_cran()
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not_installed("posterior")
  d <- mk_canon_data()
  fit_a <- suppressMessages(flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 1L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  vars <- posterior::variables(posterior::as_draws_df(fit_a))

  # The draws speak canonical names after the rename: the four
  # canonical parameters present on a Gaussian random-intercept model.
  expect_true("(Intercept)" %in% vars)
  expect_true("sd_geno" %in% vars)
  expect_true("sigma" %in% vars)
  # At least one env-level coefficient survives the rename.
  expect_true(any(grepl("^env", vars)))
})


# ---------------------------------------------------------------- #
# (d) Cross-engine brms / INLA draws via the registry                #
# ---------------------------------------------------------------- #

test_that("as_draws_df() resolves brms and INLA to common names on a Gaussian random-intercept model", {
  # Fits on brms: a Stan compile, so it is out of the
  # devtools::check() path and runs under devtools::test().
  skip_on_cran()
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not_installed("INLA")
  testthat::skip_if_not_installed("posterior")
  skip_if_not_long_tier()
  d <- mk_canon_data()
  fit_b <- suppressMessages(flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    backend = "brms",
    n_samples = 500L,
    warmup = 500L,
    chains = 1L,
    seed = 20260523L,
    control = list(adapt_delta = 0.97),
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  fit_i <- suppressMessages(flexybayes(
    yield ~ env,
    random = ~geno,
    data = d,
    backend = "inla",
    verbose = FALSE
  ))

  common <- intersect(
    posterior::variables(posterior::as_draws_df(fit_b)),
    posterior::variables(posterior::as_draws_df(fit_i))
  )
  # Canonical names that survive in both fits' renames.
  expect_true("(Intercept)" %in% common)
  expect_true("sd_geno" %in% common)
  expect_true("sigma" %in% common)
})


# ---------------------------------------------------------------- #
# (e) register_canonical_mapper() round-trip                         #
# ---------------------------------------------------------------- #

test_that("register_canonical_mapper() round-trips a custom mapper", {
  # Save / restore the registry slot to avoid contaminating the
  # package-level state for other tests.
  prev <- if (
    exists(
      "test_backend",
      envir = flexyBayes:::.canonical_mappers,
      inherits = FALSE
    )
  ) {
    get(
      "test_backend",
      envir = flexyBayes:::.canonical_mappers,
      inherits = FALSE
    )
  } else {
    NULL
  }
  on.exit({
    if (is.null(prev)) {
      rm("test_backend", envir = flexyBayes:::.canonical_mappers)
    } else {
      assign("test_backend", prev, envir = flexyBayes:::.canonical_mappers)
    }
  })

  mapper <- function(fit, fb_terms) {
    list(map = c(x = "X_canon"), transform = list())
  }
  flexyBayes:::register_canonical_mapper("test_backend", mapper)
  out <- get(
    "test_backend",
    envir = flexyBayes:::.canonical_mappers,
    inherits = FALSE
  )
  expect_identical(out, mapper)

  # Invalid inputs raise clean errors.
  expect_error(
    flexyBayes:::register_canonical_mapper(123, mapper),
    "backend.*non-empty"
  )
  expect_error(
    flexyBayes:::register_canonical_mapper("x", "not a function"),
    "must be a function"
  )
})


# ---------------------------------------------------------------- #
# (f) [removed] withdrawn engine's pin: user-supplied              #
# ---------------------------------------------------------------- #
#
# Deleted (not rewritten): this test called the withdrawn engine's
# pin and its native-model-graph ingest adapter directly, and the
# setup code called that withdrawn package's own API to build the toy
# model graph -- exactly the kind of call site the 0.9.3 withdrawal
# removes (see NEWS.md). No successor: `fb_from_brms()` and
# `fb_from_asreml()` both parse a real formula and resolve canonical
# names from the registry automatically, so neither accepts (or
# needs) a user-supplied `canonical_names` argument -- that capability
# was specific to importing a hand-built native model graph, which no
# longer exists as an ingest path.
