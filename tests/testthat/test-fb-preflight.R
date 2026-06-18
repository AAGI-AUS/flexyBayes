# Tests for the internal Stage 2 MVP `.fb_preflight()` design-memory
# preflight layer (ADR 0021 / v0.3.0). Covers:
#   - the accept path on small in-scope designs
#   - the refusal path with the binding-term identity
#   - the ceiling override (memory_ceiling_gb) accept-after-refuse
#   - aggregated_likelihood_candidate gating: in-scope per term-class,
#     in-scope per family, closed by random-slope, closed by smooth,
#     closed by non-gaussian family
#   - per-term estimate accuracy against object.size() of the actually-
#     allocated representation (the 5% ADR acceptance tolerance)
#   - the metadata-only path produces estimates without touching $data
#   - the refusal-object print surfaces both binding_term and ceiling
#   - the system-RAM probe's fallback to 8 GiB when probes absent
#   - default ceiling resolution yields a positive byte count

# Hand-build a minimal <fb_terms> IR for a given term spec; convenient
# helper for the agg_candidate / smooth / refusal tests where we want
# to bypass the brms walker's smooth refusal.
.test_make_fb_ir <- function(
  family = "gaussian",
  link = "identity",
  n_rows = 1000L,
  fixed_terms = list(),
  random_terms = list()
) {
  structure(
    list(
      response = "y",
      family = family,
      link = link,
      intercept = TRUE,
      fixed_terms = fixed_terms,
      random_terms = random_terms,
      rcov_terms = list(),
      addition_terms = list(),
      priors = list(),
      data_summary = list(n = n_rows),
      capabilities = character(),
      source = "test"
    ),
    class = c("fb_terms", "list")
  )
}


# ---------------------------------------------------------------- #
# Accept / refusal core                                             #
# ---------------------------------------------------------------- #

test_that(".fb_preflight() accepts a small in-scope design", {
  df <- data.frame(
    y = rnorm(50),
    x = rnorm(50),
    g = factor(rep(letters[1:5], 10))
  )
  fb <- fb_from_brms(y ~ x + (1 | g), data = df)
  pf <- .fb_preflight(fb, .fb_dataset(df))

  expect_s3_class(pf, "fb_preflight")
  expect_null(pf$refusal)
  expect_true(pf$total_estimate_bytes > 0)
  expect_true(pf$ceiling_bytes > 0)
  expect_identical(pf$n_rows, 50)
  # Two terms (fixed numeric x + random intercept g)
  expect_identical(length(pf$per_term_estimate), 2L)
  expect_true("x" %in% names(pf$per_term_estimate))
  expect_true("(1 | g)" %in% names(pf$per_term_estimate))
  # v0.3.10 ADR 0024 amendment (a): per-term INLA memory estimator
  # surfaces a stable-shape breakdown plus a backward-compatible
  # scalar coercion via as.numeric().
  expect_s3_class(pf$memory_estimate, "fb_memory_estimate")
  expect_identical(
    names(pf$memory_estimate),
    c("total", "breakdown", "overhead_factor", "raw_total")
  )
  expect_identical(
    names(pf$memory_estimate$breakdown),
    c("term_label", "representation", "bytes", "share")
  )
  expect_identical(
    as.numeric(pf$memory_estimate),
    pf$memory_estimate$total
  )
})

test_that(".fb_preflight() per-term memory_estimate accumulates the indexed-random + fixed contributions", {
  df <- data.frame(
    y = rnorm(60),
    x = rnorm(60),
    g = factor(rep(letters[1:5], 12))
  )
  fb <- fb_from_brms(y ~ x + (1 | g), data = df)
  pf <- .fb_preflight(fb, .fb_dataset(df))
  bd <- pf$memory_estimate$breakdown
  # Random intercept term + a single "(fixed effects)" row.
  expect_true(nrow(bd) >= 2L)
  expect_true("(fixed effects)" %in% bd$term_label)
  expect_true("indexed_random_intercept" %in% bd$representation)
  expect_true(all(bd$bytes >= 0))
  expect_true(
    abs(sum(bd$share, na.rm = TRUE) - 1) < 1e-6 ||
      all(is.na(bd$share))
  )
  expect_true(pf$memory_estimate$overhead_factor == 2)
})

test_that(".fb_preflight() refuses when design exceeds tight ceiling", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1e7,
    col_types = list(y = "double", x = "double", g = "factor"),
    dictionaries = list(g = as.character(seq_len(1e5)))
  )
  fb <- .test_make_fb_ir(
    n_rows = 1e7,
    fixed_terms = list(list(type = "smooth", var = "x", k = 10L)),
    random_terms = list(list(type = "simple", var = "g", var_n = 1e5))
  )
  pf <- .fb_preflight(fb, ds, memory_ceiling_gb = 0.5)

  expect_s3_class(pf$refusal, "fb_preflight_refusal")
  expect_identical(pf$refusal$reason_code, "design_memory_exceeds_ceiling")
  # Smooth at 1e7 x 10 dense block (~800 MB) is the binding term
  expect_identical(pf$refusal$binding_term, "s(x)")
  expect_true(pf$refusal$binding_bytes > pf$refusal$ceiling_bytes / 2)
})

test_that(".fb_preflight() ceiling override accepts a previously-refused design", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1e7,
    col_types = list(y = "double", x = "double", g = "factor"),
    dictionaries = list(g = as.character(seq_len(1e5)))
  )
  fb <- .test_make_fb_ir(
    n_rows = 1e7,
    fixed_terms = list(list(type = "smooth", var = "x", k = 10L)),
    random_terms = list(list(type = "simple", var = "g", var_n = 1e5))
  )
  expect_s3_class(
    .fb_preflight(fb, ds, memory_ceiling_gb = 0.5)$refusal,
    "fb_preflight_refusal"
  )
  expect_null(.fb_preflight(fb, ds, memory_ceiling_gb = 16)$refusal)
})

test_that(".fb_preflight() rejects malformed memory_ceiling_gb", {
  ds <- .fb_dataset(data.frame(y = 1:3, x = 1:3))
  fb <- fb_from_brms(y ~ x, data = data.frame(y = rnorm(3), x = rnorm(3)))
  expect_error(
    .fb_preflight(fb, ds, memory_ceiling_gb = -1),
    regexp = "positive"
  )
  expect_error(
    .fb_preflight(fb, ds, memory_ceiling_gb = "16"),
    regexp = "positive numeric"
  )
  expect_error(
    .fb_preflight(fb, ds, memory_ceiling_gb = c(1, 2)),
    regexp = "positive numeric"
  )
})

test_that(".fb_preflight() requires fb_terms IR + fb_dataset wrapper", {
  fb <- fb_from_brms(y ~ x, data = data.frame(y = rnorm(3), x = rnorm(3)))
  ds <- .fb_dataset(data.frame(y = 1:3, x = 1:3))
  expect_error(.fb_preflight(list(), ds), regexp = "fb_terms")
  expect_error(.fb_preflight(fb, list()), regexp = "fb_dataset")
})


# ---------------------------------------------------------------- #
# aggregated_likelihood_candidate flag                              #
# ---------------------------------------------------------------- #

test_that("agg_candidate: gaussian + fixed numeric + random intercept = TRUE", {
  df <- data.frame(y = rnorm(100), x = rnorm(100), g = factor(rep(1:10, 10)))
  fb <- fb_from_brms(y ~ x + (1 | g), data = df)
  pf <- .fb_preflight(fb, .fb_dataset(df))

  expect_true(pf$per_term_estimate$x$aggregated_likelihood_candidate)
  expect_true(pf$per_term_estimate[["(1 | g)"]]$aggregated_likelihood_candidate)
})

test_that("agg_candidate: gaussian + fixed factor = TRUE", {
  df <- data.frame(y = rnorm(100), f = factor(rep(letters[1:4], 25)))
  fb <- fb_from_brms(y ~ f, data = df)
  pf <- .fb_preflight(fb, .fb_dataset(df))

  expect_true(pf$per_term_estimate$f$aggregated_likelihood_candidate)
  expect_identical(
    pf$per_term_estimate$f$representation_class,
    "indexed_fixed_factor"
  )
})

test_that("agg_candidate: (x || g) random slope = FALSE", {
  df <- data.frame(y = rnorm(100), x = rnorm(100), g = factor(rep(1:10, 10)))
  fb <- fb_from_brms(y ~ (x || g), data = df)
  pf <- .fb_preflight(fb, .fb_dataset(df))

  # The slope term is the one with " || " in its label
  slope_idx <- grepl("||", names(pf$per_term_estimate), fixed = TRUE)
  expect_true(any(slope_idx))
  slope_entry <- pf$per_term_estimate[[which(slope_idx)[[1L]]]]
  expect_false(slope_entry$aggregated_likelihood_candidate)
  expect_identical(slope_entry$representation_class, "dense_baseline")
})

test_that("agg_candidate: smooth s(x) = FALSE; representation_class = dense_baseline", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", x = "double"),
    dictionaries = list()
  )
  fb <- .test_make_fb_ir(
    n_rows = 1000L,
    fixed_terms = list(list(type = "smooth", var = "x", k = 10L))
  )
  pf <- .fb_preflight(fb, ds)

  expect_false(pf$per_term_estimate[["s(x)"]]$aggregated_likelihood_candidate)
  expect_identical(
    pf$per_term_estimate[["s(x)"]]$representation_class,
    "dense_baseline"
  )
})

test_that("agg_candidate: non-gaussian family closes the flag on every term", {
  df <- data.frame(y = rpois(100, 2), x = rnorm(100), g = factor(rep(1:10, 10)))
  fb <- fb_from_brms(y ~ x + (1 | g), data = df, family = "poisson")
  pf <- .fb_preflight(fb, .fb_dataset(df))

  expect_false(pf$per_term_estimate$x$aggregated_likelihood_candidate)
  expect_false(
    pf$per_term_estimate[["(1 | g)"]]$aggregated_likelihood_candidate
  )
})


# ---------------------------------------------------------------- #
# Estimate accuracy vs object.size() (5% acceptance tolerance)      #
# ---------------------------------------------------------------- #

test_that("estimate accuracy: fixed numeric within 5% of object.size()", {
  N <- 5000L
  df <- data.frame(y = rnorm(N), x = rnorm(N))
  fb <- fb_from_brms(y ~ x, data = df)
  pf <- .fb_preflight(fb, .fb_dataset(df))

  predicted <- pf$per_term_estimate$x$design_memory_bytes - 256
  actual <- as.numeric(object.size(df$x)) +
    as.numeric(object.size(numeric(1L)))
  expect_lte(abs(predicted - actual) / actual, 0.05)
})

test_that("estimate accuracy: fixed factor within 5% of object.size()", {
  N <- 5000L
  k <- 8L
  df <- data.frame(y = rnorm(N), f = factor(sample.int(k, N, replace = TRUE)))
  fb <- fb_from_brms(y ~ f, data = df)
  pf <- .fb_preflight(fb, .fb_dataset(df))

  predicted <- pf$per_term_estimate$f$design_memory_bytes - 256
  # Integer index (4*N) + (k-1) beta-dummy vector
  actual <- as.numeric(object.size(as.integer(df$f))) +
    as.numeric(object.size(numeric(k - 1L)))
  expect_lte(abs(predicted - actual) / actual, 0.05)
})

test_that("estimate accuracy: random intercept within 5% of object.size()", {
  N <- 5000L
  k_g <- 200L
  g <- factor(sample.int(k_g, N, replace = TRUE))
  df <- data.frame(y = rnorm(N), g = g)
  fb <- fb_from_brms(y ~ (1 | g), data = df)
  pf <- .fb_preflight(fb, .fb_dataset(df))

  predicted <- pf$per_term_estimate[["(1 | g)"]]$design_memory_bytes - 256
  actual <- as.numeric(object.size(as.integer(g))) +
    as.numeric(object.size(numeric(k_g)))
  expect_lte(abs(predicted - actual) / actual, 0.05)
})

test_that("estimate accuracy: smooth s(x) basis block within 5%", {
  N <- 5000L
  k <- 10L
  ds <- .fb_dataset(
    data = NULL,
    n_rows = N,
    col_types = list(y = "double", x = "double"),
    dictionaries = list()
  )
  fb <- .test_make_fb_ir(
    n_rows = N,
    fixed_terms = list(list(type = "smooth", var = "x", k = k))
  )
  pf <- .fb_preflight(fb, ds)

  predicted <- pf$per_term_estimate[["s(x)"]]$design_memory_bytes - 256
  # Dense N x k double block + k basis coefs
  actual <- as.numeric(object.size(matrix(0, N, k))) +
    as.numeric(object.size(numeric(k)))
  expect_lte(abs(predicted - actual) / actual, 0.05)
})


# ---------------------------------------------------------------- #
# Metadata-only path: O(terms), independent of n_rows               #
# ---------------------------------------------------------------- #

test_that(".fb_preflight() runs on metadata-only fb_dataset without touching $data", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1e6,
    col_types = list(y = "double", x = "double", g = "factor"),
    dictionaries = list(g = as.character(seq_len(1000L)))
  )
  expect_true(.fb_dataset_is_metadata(ds))
  expect_null(ds$data)

  fb <- .test_make_fb_ir(
    n_rows = 1e6,
    fixed_terms = list(list(type = "numeric", var = "x")),
    random_terms = list(list(type = "simple", var = "g", var_n = 1000L))
  )
  pf <- .fb_preflight(fb, ds, memory_ceiling_gb = 16)

  expect_null(pf$refusal)
  expect_identical(length(pf$per_term_estimate), 2L)
  # The random-intercept estimate uses the dictionary's 1000 levels
  expect_true(pf$per_term_estimate[["(1 | g)"]]$design_memory_bytes > 4 * 1e6)
})


# ---------------------------------------------------------------- #
# Ceiling resolution + RAM probe fallback                           #
# ---------------------------------------------------------------- #

test_that("default ceiling resolves to a positive byte count", {
  df <- data.frame(y = rnorm(50), x = rnorm(50))
  fb <- fb_from_brms(y ~ x, data = df)
  pf <- .fb_preflight(fb, .fb_dataset(df))

  expect_true(pf$ceiling_bytes > 0)
  # Host-robust: the default ceiling is the documented formula -- 60% of
  # probed RAM, or an 8 GiB fallback when the probe returns NA -- so
  # assert against the actual probe rather than assuming the host has
  # >= 8 GiB (CI runners may have less). Audit P2.3 (2026-05-23) lowered
  # the default fraction from 0.8 to 0.6.
  ram <- flexyBayes:::.fb_probe_system_ram()
  frac <- getOption("flexyBayes.preflight_ram_fraction", 0.6)
  expected <- if (is.na(ram)) 8 * 1024^3 else as.numeric(frac) * ram
  expect_equal(pf$ceiling_bytes, expected, tolerance = 1)
})

test_that("ceiling fallback returns 8 GiB when probes return NA", {
  # Stub the RAM probe to simulate the no-probe environment. Save
  # the original via getFromNamespace() and restore via on.exit so
  # the test does not leak across to the rest of the suite.
  ns <- asNamespace("flexyBayes")
  original <- get(".fb_probe_system_ram", envir = ns, inherits = FALSE)
  stub_fn <- function(...) NA_real_
  assignInNamespace(".fb_probe_system_ram", stub_fn, ns = "flexyBayes")
  on.exit(
    assignInNamespace(".fb_probe_system_ram", original, ns = "flexyBayes"),
    add = TRUE
  )

  df <- data.frame(y = rnorm(50), x = rnorm(50))
  fb <- fb_from_brms(y ~ x, data = df)
  pf <- .fb_preflight(fb, .fb_dataset(df))
  expect_identical(pf$ceiling_bytes, 8 * 1024^3)
})


# ---------------------------------------------------------------- #
# Print contents + class invariants                                 #
# ---------------------------------------------------------------- #

test_that("print.fb_preflight_refusal() surfaces binding_term + ceiling", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1e7,
    col_types = list(y = "double", x = "double", g = "factor"),
    dictionaries = list(g = as.character(seq_len(1e5)))
  )
  fb <- .test_make_fb_ir(
    n_rows = 1e7,
    fixed_terms = list(list(type = "smooth", var = "x", k = 10L)),
    random_terms = list(list(type = "simple", var = "g", var_n = 1e5))
  )
  pf <- .fb_preflight(fb, ds, memory_ceiling_gb = 0.5)

  out <- capture.output(print(pf$refusal))
  expect_true(any(grepl("binding_term:", out, fixed = TRUE)))
  expect_true(any(grepl("s(x)", out, fixed = TRUE)))
  expect_true(any(grepl("ceiling_bytes:", out, fixed = TRUE)))
  expect_true(any(grepl("override:", out, fixed = TRUE)))
  expect_true(any(grepl("memory_ceiling_gb", out, fixed = TRUE)))
})

test_that("class invariants: <fb_preflight> + <fb_preflight_refusal>", {
  df <- data.frame(y = rnorm(50), x = rnorm(50))
  fb <- fb_from_brms(y ~ x, data = df)
  pf <- .fb_preflight(fb, .fb_dataset(df))

  expect_s3_class(pf, "fb_preflight")
  # Inner class is "list" so length() / names() work as expected
  expect_true("list" %in% class(pf))

  # Refusal class only attached when refusing
  ds_big <- .fb_dataset(
    data = NULL,
    n_rows = 1e8,
    col_types = list(y = "double", g = "factor"),
    dictionaries = list(g = as.character(seq_len(1e7)))
  )
  fb_big <- .test_make_fb_ir(
    n_rows = 1e8,
    random_terms = list(list(type = "simple", var = "g", var_n = 1e7))
  )
  pf_big <- .fb_preflight(fb_big, ds_big, memory_ceiling_gb = 0.1)
  expect_s3_class(pf_big$refusal, "fb_preflight_refusal")
})

test_that(".fb_preflight is internal -- no exported binding", {
  ns <- asNamespace("flexyBayes")
  exp <- getNamespaceExports(ns)
  expect_false(".fb_preflight" %in% exp)
  expect_false("fb_preflight" %in% exp)
  expect_false(".fb_preflight_refusal" %in% exp)
  expect_true(exists(".fb_preflight", envir = ns, inherits = FALSE))
})


# ---------------------------------------------------------------- #
# fb_from_brms(data = NULL, carry_n_rows) + dataset round-trip      #
# ---------------------------------------------------------------- #

test_that("fb_from_brms(data = NULL, carry_n_rows): metadata-only IR shape", {
  fb <- fb_from_brms(y ~ x + (1 | g), data = NULL, carry_n_rows = 1e7)

  expect_s3_class(fb, "fb_terms")
  expect_identical(fb$source, "brms_metadata_only")
  expect_identical(fb$data_summary$n, 10000000L)
  # Random-term level cache stripped -- preflight will read from
  # <fb_dataset>$dictionaries.
  expect_true(is.na(fb$random_terms[[1L]]$var_n))
  expect_null(fb$random_terms[[1L]]$var_levels)
})

test_that("fb_from_brms(data = NULL) without carry_n_rows refuses", {
  expect_error(fb_from_brms(y ~ x, data = NULL), regexp = "carry_n_rows")
  expect_error(
    fb_from_brms(y ~ x, data = NULL, carry_n_rows = -10),
    regexp = "positive numeric"
  )
})

test_that("fb_from_brms(data, carry_n_rows): override path preserves real var_n", {
  df <- data.frame(y = rnorm(10), x = rnorm(10), g = factor(letters[1:10]))
  fb <- fb_from_brms(y ~ x + (1 | g), data = df, carry_n_rows = 5e6)

  expect_identical(fb$source, "brms")
  expect_identical(fb$data_summary$n, 5000000L)
  # Real var_n from the placeholder data is preserved -- it is the
  # caller's responsibility to ensure dictionaries on the <fb_dataset>
  # match the actual large-scale workflow when they override n.
  expect_identical(fb$random_terms[[1L]]$var_n, 10L)
})

test_that("carry_n_rows round-trip: metadata-only IR + dataset -> preflight runs", {
  fb <- fb_from_brms(y ~ x + (1 | g), data = NULL, carry_n_rows = 1e7)
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1e7,
    col_types = list(y = "double", x = "double", g = "factor"),
    dictionaries = list(g = as.character(seq_len(1e5)))
  )
  pf <- .fb_preflight(fb, ds, memory_ceiling_gb = 16)

  expect_null(pf$refusal)
  expect_identical(pf$n_rows, 1e7)
  # (1 | g) estimate uses the dictionary's 1e5 level count plus the
  # 4 * 1e7 integer index column. Sanity-check the order of magnitude.
  re_bytes <- pf$per_term_estimate[["(1 | g)"]]$design_memory_bytes
  expect_gte(re_bytes, 4e7)
  expect_lte(re_bytes, 4.2e7)
})


# ---------------------------------------------------------------- #
# Phase B hardening (2026-05-23) -- audit-driven new estimators     #
# ---------------------------------------------------------------- #

test_that("smooth_mgcv parser type estimates N x k dense basis", {
  # Audit P1.2: parser emits type = "smooth_mgcv"; previously fell
  # through to random_other (12N + 16k) for a 3.3x under-estimate at
  # N = 5000, k = 5. Expected estimate now matches the audit number.
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 5000L,
    col_types = list(y = "double", x = "double"),
    dictionaries = list()
  )
  fb <- .test_make_fb_ir(
    n_rows = 5000L,
    random_terms = list(list(type = "smooth_mgcv", var = "x", k = 5L))
  )
  pf <- .fb_preflight(fb, ds)

  entry <- pf$per_term_estimate[["s(x)"]]
  expect_identical(entry$term_kind, "random_smooth")
  expect_identical(entry$representation_class, "dense_baseline")
  # 8 * 5000 * 5 + 8 * 5 + 256 = 200,296.
  expect_equal(entry$design_memory_bytes, 8 * 5000 * 5 + 8 * 5 + 256)
  expect_null(pf$refusal)
})

test_that("spline (spl) parser type estimates N x df dense basis at df=8", {
  # Audit P1.2: parser emits type = "spline"; codegen hardcodes
  # df = 8. Previously fell through to random_other (12N + 16k) for
  # a 5.3x under-estimate at N = 5000.
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 5000L,
    col_types = list(y = "double", x = "double"),
    dictionaries = list()
  )
  fb <- .test_make_fb_ir(
    n_rows = 5000L,
    random_terms = list(list(type = "spline", var = "x"))
  )
  pf <- .fb_preflight(fb, ds)

  entry <- pf$per_term_estimate[["spl(x)"]]
  expect_identical(entry$term_kind, "random_spline")
  # 8 * 5000 * 8 + 8 * 8 + 256 = 320,320.
  expect_equal(entry$design_memory_bytes, 8 * 5000 * 8 + 8 * 8 + 256)
})

test_that("factor:factor interaction sizes dense treatment-coded block", {
  # Audit P1.3: high-cardinality factor:factor is the named scary
  # case. Previously fell through to a single-column estimate
  # (8N + 8) -- the dense block is (K-1)x larger.
  df <- data.frame(
    y = rnorm(1000),
    f1 = factor(rep(letters[1:10], 100)),
    f2 = factor(rep(LETTERS[1:5], each = 200))
  )
  ds <- .fb_dataset(df)
  fb <- .test_make_fb_ir(
    n_rows = 1000L,
    fixed_terms = list(list(
      type = "factor_interaction",
      vars = c("f1", "f2"),
      label = "f1:f2"
    ))
  )
  pf <- .fb_preflight(fb, ds)

  entry <- pf$per_term_estimate[["f1:f2"]]
  expect_identical(entry$term_kind, "fixed_factor_interaction")
  # Reduced-rank treatment coding: model.matrix(~ f1 * f2) interaction
  # block alone has prod(L_i - 1) = 9 * 4 = 36 columns. The estimate
  # matches the common additive-with-interaction shape exactly.
  expect_equal(entry$design_memory_bytes, 8 * 1000 * 36 + 8 * 36 + 256)
})

test_that("factor:factor refuses when cardinality unresolvable", {
  # Metadata-only dataset with no dictionaries -- the factor:factor
  # estimator cannot compute K and flags unknown_representation.
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", f1 = "factor", f2 = "factor"),
    dictionaries = list()
  )
  fb <- .test_make_fb_ir(
    n_rows = 1000L,
    fixed_terms = list(list(
      type = "factor_interaction",
      vars = c("f1", "f2"),
      label = "f1:f2"
    ))
  )
  pf <- .fb_preflight(fb, ds)

  expect_false(is.null(pf$refusal))
  expect_identical(
    pf$refusal$reason_code,
    "representation_unknown_for_preflight"
  )
  expect_identical(pf$refusal$binding_term, "f1:f2")
  expect_true(is.na(pf$refusal$binding_bytes))
  expect_true(is.na(pf$refusal$total_bytes))
})

test_that("factor_numeric_interaction (ADR 0019) sizes indexed slope block", {
  df <- data.frame(
    y = rnorm(1000),
    x = rnorm(1000),
    f = factor(rep(letters[1:5], 200))
  )
  ds <- .fb_dataset(df)
  # Slot names match parse_formula.R:67 ($factor, $continuous, $n_levels).
  fb <- .test_make_fb_ir(
    n_rows = 1000L,
    fixed_terms = list(list(
      type = "factor_numeric_interaction",
      factor = "f",
      continuous = "x",
      vars = c("f", "x"),
      n_levels = 5L,
      label = "f:x"
    ))
  )
  pf <- .fb_preflight(fb, ds)

  entry <- pf$per_term_estimate[["f:x"]]
  expect_identical(entry$term_kind, "fixed_factor_numeric_interaction")
  # 4 * 1000 + 8 * 1000 + 8 * (5 - 1) + 256 = 12,288.
  expect_equal(entry$design_memory_bytes, 4 * 1000 + 8 * 1000 + 8 * 4 + 256)
})

test_that("polynomial term sizes degree-d dense basis", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", x = "double"),
    dictionaries = list()
  )
  fb <- .test_make_fb_ir(
    n_rows = 1000L,
    fixed_terms = list(list(
      type = "polynomial",
      var = "x",
      degree = 4L,
      label = "pol(x,4)"
    ))
  )
  pf <- .fb_preflight(fb, ds)

  entry <- pf$per_term_estimate[["pol(x,4)"]]
  expect_identical(entry$term_kind, "fixed_polynomial")
  # 8 * 1000 * 4 + 8 * 4 + 256 = 32,288.
  expect_equal(entry$design_memory_bytes, 8 * 1000 * 4 + 8 * 4 + 256)
})

test_that("expression term (I() / arithmetic) sizes one column", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", x = "double"),
    dictionaries = list()
  )
  fb <- .test_make_fb_ir(
    n_rows = 1000L,
    fixed_terms = list(list(type = "expression", label = "I(x^2)"))
  )
  pf <- .fb_preflight(fb, ds)

  entry <- pf$per_term_estimate[["I(x^2)"]]
  expect_identical(entry$term_kind, "fixed_expression")
  # 8 * 1000 + 8 + 256 = 8,264. NOT flagged unknown.
  expect_equal(entry$design_memory_bytes, 8 * 1000 + 8 + 256)
  expect_null(pf$refusal)
})

test_that("interaction_generic random term triggers representation_unknown", {
  # parse_formula.R:214 fallback for random terms whose interaction
  # shape didn't match any structured-cov pattern. The audit P1.3
  # framing: unknown shape != "estimate a small column".
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", g = "factor"),
    dictionaries = list(g = letters[1:5])
  )
  fb <- .test_make_fb_ir(
    n_rows = 1000L,
    random_terms = list(list(
      type = "interaction_generic",
      left = list(type = "vm", var = "g", mat = "M"),
      right = list(type = "polynomial", var = "x")
    ))
  )
  pf <- .fb_preflight(fb, ds)

  expect_false(is.null(pf$refusal))
  expect_identical(
    pf$refusal$reason_code,
    "representation_unknown_for_preflight"
  )
})

test_that("unknown fixed type triggers representation_unknown refusal", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", x = "double"),
    dictionaries = list()
  )
  fb <- .test_make_fb_ir(
    n_rows = 1000L,
    fixed_terms = list(list(
      type = "exotic_unsupported_type",
      var = "x",
      label = "exotic(x)"
    ))
  )
  pf <- .fb_preflight(fb, ds)

  expect_false(is.null(pf$refusal))
  expect_identical(
    pf$refusal$reason_code,
    "representation_unknown_for_preflight"
  )
  expect_identical(pf$refusal$binding_term, "exotic(x)")
})


# ---------------------------------------------------------------- #
# flexyBayes.preflight_ram_fraction option                          #
# ---------------------------------------------------------------- #

test_that("preflight_ram_fraction scales the default ceiling", {
  # Stub the probe so we always exercise the fraction code path
  # (8 GiB fallback ignores the fraction by design).
  ns <- asNamespace("flexyBayes")
  original <- get(".fb_probe_system_ram", envir = ns, inherits = FALSE)
  on.exit(assignInNamespace(
    ".fb_probe_system_ram",
    original,
    ns = "flexyBayes"
  ))
  assignInNamespace(
    ".fb_probe_system_ram",
    function(...) 16 * 1024^3,
    ns = "flexyBayes"
  )

  withr::local_options(flexyBayes.preflight_ram_fraction = 0.3)
  c_low <- .fb_resolve_ceiling(NULL)

  withr::local_options(flexyBayes.preflight_ram_fraction = 0.9)
  c_high <- .fb_resolve_ceiling(NULL)

  expect_equal(c_low, 0.3 * 16 * 1024^3)
  expect_equal(c_high, 0.9 * 16 * 1024^3)
  expect_equal(c_high / c_low, 3, tolerance = 1e-9)
})

test_that("preflight_ram_fraction default is 0.6 (audit P2.3)", {
  ns <- asNamespace("flexyBayes")
  original <- get(".fb_probe_system_ram", envir = ns, inherits = FALSE)
  on.exit(assignInNamespace(
    ".fb_probe_system_ram",
    original,
    ns = "flexyBayes"
  ))
  assignInNamespace(
    ".fb_probe_system_ram",
    function(...) 16 * 1024^3,
    ns = "flexyBayes"
  )

  withr::local_options(flexyBayes.preflight_ram_fraction = NULL)
  ceiling_default <- .fb_resolve_ceiling(NULL)
  expect_equal(ceiling_default, 0.6 * 16 * 1024^3)
})

test_that("preflight_ram_fraction rejects out-of-range values", {
  ns <- asNamespace("flexyBayes")
  original <- get(".fb_probe_system_ram", envir = ns, inherits = FALSE)
  on.exit(assignInNamespace(
    ".fb_probe_system_ram",
    original,
    ns = "flexyBayes"
  ))
  assignInNamespace(
    ".fb_probe_system_ram",
    function(...) 8 * 1024^3,
    ns = "flexyBayes"
  )

  withr::local_options(flexyBayes.preflight_ram_fraction = 0)
  expect_error(.fb_resolve_ceiling(NULL), regexp = "preflight_ram_fraction")

  withr::local_options(flexyBayes.preflight_ram_fraction = 1.5)
  expect_error(.fb_resolve_ceiling(NULL), regexp = "preflight_ram_fraction")

  withr::local_options(flexyBayes.preflight_ram_fraction = -0.1)
  expect_error(.fb_resolve_ceiling(NULL), regexp = "preflight_ram_fraction")
})


# ---------------------------------------------------------------- #
# print.fb_preflight_refusal -- representation_unknown shape         #
# ---------------------------------------------------------------- #

test_that("print.fb_preflight_refusal: NA bytes render cleanly", {
  ref <- structure(
    list(
      reason_code = "representation_unknown_for_preflight",
      binding_term = "exotic(x)",
      binding_bytes = NA_real_,
      total_bytes = NA_real_,
      ceiling_bytes = 8 * 1024^3,
      n_rows = 1000L
    ),
    class = c("fb_preflight_refusal", "list")
  )
  out <- capture.output(print(ref))
  expect_true(any(grepl(
    "representation_unknown_for_preflight",
    out,
    fixed = TRUE
  )))
  expect_true(any(grepl("NA (not characterised)", out, fixed = TRUE)))
  expect_true(any(grepl("remedy:", out, fixed = TRUE)))
  # No override-suggest line on unknown_for_preflight.
  expect_false(any(grepl("override:", out, fixed = TRUE)))
})


# ---------------------------------------------------------------- #
# Audit P1.3 (2026-05-24): .stop_preflight_refusal() must branch   #
# its headline on refusal$reason_code so that an unknown-          #
# representation refusal does NOT falsely advertise a memory-      #
# ceiling excess (raising the ceiling cannot help) and does NOT    #
# suggest a ceiling override. The detailed body via                #
# print.fb_preflight_refusal() is exercised separately above.      #
# ---------------------------------------------------------------- #

test_that("dispatch refusal headline branches on reason_code", {
  ref_mem <- structure(
    list(
      reason_code = "design_memory_exceeds_ceiling",
      binding_term = "(1 | g)",
      binding_bytes = 2 * 1024^3,
      total_bytes = 2 * 1024^3,
      ceiling_bytes = 1 * 1024^3,
      n_rows = 1e6L
    ),
    class = c("fb_preflight_refusal", "list")
  )
  ref_unk <- structure(
    list(
      reason_code = "representation_unknown_for_preflight",
      binding_term = "exotic(x)",
      binding_bytes = NA_real_,
      total_bytes = NA_real_,
      ceiling_bytes = 8 * 1024^3,
      n_rows = 1e6L
    ),
    class = c("fb_preflight_refusal", "list")
  )

  err_mem <- tryCatch(
    flexyBayes:::.stop_preflight_refusal(ref_mem, call = NULL),
    flexybayes_preflight_refusal = function(c) c
  )
  err_unk <- tryCatch(
    flexyBayes:::.stop_preflight_refusal(ref_unk, call = NULL),
    flexybayes_preflight_refusal = function(c) c
  )

  expect_s3_class(err_mem, "flexybayes_preflight_refusal")
  expect_s3_class(err_unk, "flexybayes_preflight_refusal")

  # design_memory_exceeds_ceiling: headline mentions the ceiling.
  expect_true(grepl(
    "exceeds the active memory ceiling",
    conditionMessage(err_mem),
    fixed = TRUE
  ))

  # representation_unknown_for_preflight: headline must NOT mention
  # ceiling excess, and must NOT suggest raising the ceiling.
  expect_false(grepl(
    "exceeds the active memory ceiling",
    conditionMessage(err_unk),
    fixed = TRUE
  ))
  expect_true(grepl(
    "representation is not characterised",
    conditionMessage(err_unk),
    fixed = TRUE
  ))
  expect_true(grepl(
    "raising the memory ceiling will not help",
    conditionMessage(err_unk),
    fixed = TRUE
  ))
})


# ---------------------------------------------------------------- #
# Audit P1.4 (2026-05-24): structured-covariance preflight split   #
# -- known-indexed forms (vm, ped, ar1, ar1_spatial) keep the      #
# tested indexed heuristic; every other structured-cov form        #
# (at, us, fa, at_simple, us_gxe, fa_gxe, vm_gxe, ...) routes to   #
# representation_unknown_for_preflight rather than silently        #
# under-estimating at large N. The .preflight_term_level_count()   #
# fallback that previously returned `1L` was removed; NA level     #
# count escalates every consumer branch to unknown_representation. #
# ---------------------------------------------------------------- #

test_that("vm() with resolvable k uses the indexed_structured_known estimator", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", g = "factor"),
    dictionaries = list(g = letters[1:10])
  )
  fb <- .test_make_fb_ir(
    n_rows = 1000L,
    random_terms = list(list(type = "vm", var = "g", var_n = 10L, mat = "G"))
  )
  pf <- .fb_preflight(fb, ds)

  expect_null(pf$refusal)
  entry <- pf$per_term_estimate[[1]]
  expect_identical(entry$representation_class, "indexed_structured_known")
  expect_false(isTRUE(entry$unknown_representation))
  expect_false(is.na(entry$design_memory_bytes))
})

test_that("ped() with resolvable k uses the indexed_structured_known estimator", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 500L,
    col_types = list(y = "double", anim = "factor"),
    dictionaries = list(anim = as.character(1:50))
  )
  fb <- .test_make_fb_ir(
    n_rows = 500L,
    random_terms = list(list(type = "ped", var = "anim", var_n = 50L))
  )
  pf <- .fb_preflight(fb, ds)

  expect_null(pf$refusal)
  entry <- pf$per_term_estimate[[1]]
  expect_identical(entry$representation_class, "indexed_structured_known")
})

test_that("vm() with NA k escalates to representation_unknown", {
  # Metadata-only dataset without a dictionary for the grouping
  # factor -- the level count is unresolvable.
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", g = "factor"),
    dictionaries = list()
  )
  fb <- .test_make_fb_ir(
    n_rows = 1000L,
    random_terms = list(list(type = "vm", var = "g"))
  )
  pf <- .fb_preflight(fb, ds)

  expect_false(is.null(pf$refusal))
  expect_identical(
    pf$refusal$reason_code,
    "representation_unknown_for_preflight"
  )
})

test_that("at() routes to representation_unknown regardless of k", {
  # at() is a known dense-matrix structured-cov form that the
  # indexed heuristic under-estimates. Per P1.4 it routes to
  # representation_unknown_for_preflight pending a term-specific
  # estimator.
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", env = "factor"),
    dictionaries = list(env = letters[1:10])
  )
  fb <- .test_make_fb_ir(
    n_rows = 1000L,
    random_terms = list(list(type = "at", var = "env", var_n = 10L))
  )
  pf <- .fb_preflight(fb, ds)

  expect_false(is.null(pf$refusal))
  expect_identical(
    pf$refusal$reason_code,
    "representation_unknown_for_preflight"
  )
})

test_that("us() routes to representation_unknown regardless of k", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", env = "factor"),
    dictionaries = list(env = letters[1:10])
  )
  fb <- .test_make_fb_ir(
    n_rows = 1000L,
    random_terms = list(list(type = "us", var = "env", var_n = 10L))
  )
  pf <- .fb_preflight(fb, ds)

  expect_false(is.null(pf$refusal))
  expect_identical(
    pf$refusal$reason_code,
    "representation_unknown_for_preflight"
  )
})

test_that("fa() routes to representation_unknown regardless of k", {
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 2000L,
    col_types = list(y = "double", geno = "factor"),
    dictionaries = list(geno = as.character(1:100))
  )
  fb <- .test_make_fb_ir(
    n_rows = 2000L,
    random_terms = list(list(type = "fa", var = "geno", var_n = 100L, k = 2L))
  )
  pf <- .fb_preflight(fb, ds)

  expect_false(is.null(pf$refusal))
  expect_identical(
    pf$refusal$reason_code,
    "representation_unknown_for_preflight"
  )
})

test_that("random intercept with NA k escalates to representation_unknown", {
  # Audit P1.4 (2026-05-24): the previous `1L` fallback masked
  # this case; the random-intercept branch now treats NA k as
  # the signal to refuse via representation_unknown.
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 5000L,
    col_types = list(y = "double", g = "factor"),
    dictionaries = list() # no dictionary; level count unknown
  )
  fb <- .test_make_fb_ir(
    n_rows = 5000L,
    random_terms = list(list(type = "simple", var = "g"))
  )
  pf <- .fb_preflight(fb, ds)

  expect_false(is.null(pf$refusal))
  expect_identical(
    pf$refusal$reason_code,
    "representation_unknown_for_preflight"
  )
})

test_that("factor fixed term with NA k escalates to representation_unknown", {
  # Mirror of the random-intercept regression for the fixed-factor
  # branch -- the dummy-block size depends on k - 1.
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 5000L,
    col_types = list(y = "double", f = "factor"),
    dictionaries = list()
  )
  fb <- .test_make_fb_ir(
    n_rows = 5000L,
    fixed_terms = list(list(type = "factor", var = "f"))
  )
  pf <- .fb_preflight(fb, ds)

  expect_false(is.null(pf$refusal))
  expect_identical(
    pf$refusal$reason_code,
    "representation_unknown_for_preflight"
  )
})
