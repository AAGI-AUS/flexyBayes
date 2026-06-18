# Grammar polymorphism on the universal entry (ADR 0031 Phase 3).
#
# fb() / flexybayes() detect the input grammar from the call shape and
# route to the matching ingest adapter, producing the backend-agnostic
# fb_terms IR. These tests pin the detection rule, the syntax override,
# the structured refusals, and the exported ingest adapters.

mk_df <- function(n = 30) {
  set.seed(11)
  data.frame(
    y = rnorm(n),
    x = rnorm(n),
    g = factor(rep(letters[1:5], length.out = n)),
    env = factor(rep(c("a", "b"), length.out = n))
  )
}

# ---------------------------------------------------------------- #
# Detection rule                                                    #
# ---------------------------------------------------------------- #

test_that("a bar-grouped formula is detected as brms grammar", {
  expect_identical(flexyBayes:::.detect_grammar(y ~ x + (1 | g)), "brms")
  expect_identical(flexyBayes:::.detect_grammar(y ~ x + (1 || g)), "brms")
})

test_that("a bar-free formula defaults to asreml grammar", {
  expect_identical(flexyBayes:::.detect_grammar(y ~ x), "asreml")
  expect_identical(flexyBayes:::.detect_grammar(y ~ x + env), "asreml")
})

test_that("a greta_model object is detected as greta grammar", {
  fake <- structure(list(), class = "greta_model")
  expect_identical(flexyBayes:::.detect_grammar(fake), "greta")
})

test_that("an explicit syntax argument overrides shape detection", {
  expect_identical(flexyBayes:::.detect_grammar(y ~ x, syntax = "brms"), "brms")
  expect_identical(
    flexyBayes:::.detect_grammar(y ~ x + (1 | g), syntax = "asreml"),
    "asreml"
  )
})

# ---------------------------------------------------------------- #
# IR routing through the universal entry                            #
# ---------------------------------------------------------------- #

test_that("flexybayes() builds an asreml-source IR from the ASReml form", {
  df <- mk_df()
  ir <- flexyBayes:::.build_ir_polymorphic(
    y ~ x,
    ~g,
    NULL,
    df,
    "gaussian",
    NULL,
    NULL,
    list(),
    NULL,
    100,
    1,
    "auto"
  )
  expect_identical(ir$source, "asreml")
})

test_that("flexybayes() builds a brms-source IR from a bar-grouped formula", {
  df <- mk_df()
  ir <- flexyBayes:::.build_ir_polymorphic(
    y ~ x + (1 | g),
    NULL,
    NULL,
    df,
    "gaussian",
    NULL,
    NULL,
    list(),
    NULL,
    100,
    1,
    "auto"
  )
  expect_identical(ir$source, "brms")
})

test_that("the asreml IR is byte-identical to a direct fb_from_asreml() call", {
  df <- mk_df()
  via_entry <- flexyBayes:::.build_ir_polymorphic(
    y ~ x,
    ~g,
    NULL,
    df,
    "gaussian",
    NULL,
    NULL,
    list(),
    NULL,
    100,
    1,
    "auto"
  )
  direct <- fb_from_asreml(y ~ x, random = ~g, data = df)
  expect_identical(via_entry, direct)
})

# ---------------------------------------------------------------- #
# Structured refusals                                               #
# ---------------------------------------------------------------- #

test_that("a brms formula combined with ASReml random/rcov refuses", {
  df <- mk_df()
  expect_error(
    flexyBayes:::.build_ir_polymorphic(
      y ~ x + (1 | g),
      ~g,
      NULL,
      df,
      "gaussian",
      NULL,
      NULL,
      list(),
      NULL,
      100,
      1,
      "auto"
    ),
    class = "flexybayes_refusal_grammar_brms_with_asreml_terms"
  )
})

test_that("brms grammar with known_matrices refuses", {
  df <- mk_df()
  expect_error(
    flexyBayes:::.build_ir_polymorphic(
      y ~ x + (1 | g),
      NULL,
      NULL,
      df,
      "gaussian",
      NULL,
      NULL,
      list(G = diag(2)),
      NULL,
      100,
      1,
      "auto"
    ),
    class = "flexybayes_refusal_grammar_brms_known_matrices_unsupported"
  )
})

test_that("a native greta_model with ASReml random/rcov refuses", {
  # v0.5.0: the universal entry now FITS a native greta_model (no longer
  # deferred). But a native graph encodes its full structure itself, so
  # combining it with the ASReml `random` / `rcov` slots is a category
  # error and refuses.
  df <- mk_df()
  fake <- structure(list(), class = "greta_model")
  expect_error(
    flexyBayes:::.build_ir_polymorphic(
      fake,
      ~g,
      NULL,
      df,
      "gaussian",
      NULL,
      NULL,
      list(),
      NULL,
      100,
      1,
      "auto"
    ),
    "cannot be combined with `random`"
  )
})

test_that("a native greta_model on the universal entry builds a greta-source IR", {
  skip_if_greta_backend_unusable()
  # A real (fittable) graph lowers to a greta-source IR via fb_from_greta().
  set.seed(3)
  yy <- greta::as_data(rnorm(20))
  xx <- greta::as_data(rnorm(20))
  b0 <- greta::normal(0, 10)
  b1 <- greta::normal(0, 10)
  s <- greta::uniform(0, 5)
  greta::distribution(yy) <- greta::normal(b0 + b1 * xx, s)
  m <- greta::model(b0, b1, s)
  ir <- suppressMessages(
    flexyBayes:::.build_ir_polymorphic(
      m,
      NULL,
      NULL,
      NULL,
      "gaussian",
      NULL,
      NULL,
      list(),
      NULL,
      100,
      1,
      "auto"
    )
  )
  expect_s3_class(ir, "fb_terms")
  expect_identical(ir$source, "greta")
  expect_true(inherits(ir$greta_meta$model, "greta_model"))
})

test_that("a prebuilt greta-source IR passes through the universal entry", {
  skip_if_greta_backend_unusable()
  set.seed(4)
  yy <- greta::as_data(rnorm(20))
  xx <- greta::as_data(rnorm(20))
  b0 <- greta::normal(0, 10)
  b1 <- greta::normal(0, 10)
  s <- greta::uniform(0, 5)
  greta::distribution(yy) <- greta::normal(b0 + b1 * xx, s)
  m <- greta::model(b0, b1, s)
  ir <- suppressMessages(
    fb_from_greta(m, canonical_names = c(b0 = "(Intercept)", b1 = "x"))
  )
  out <- flexyBayes:::.build_ir_polymorphic(
    ir,
    NULL,
    NULL,
    NULL,
    "gaussian",
    NULL,
    NULL,
    list(),
    NULL,
    100,
    1,
    "auto"
  )
  expect_identical(out, ir)
  # A prebuilt asreml IR is NOT accepted (its emit-display path needs the
  # formula triple); only greta-source IRs pass through.
  asreml_ir <- fb_from_asreml(y ~ x, random = ~g, data = mk_df())
  expect_error(
    flexyBayes:::.build_ir_polymorphic(
      asreml_ir,
      NULL,
      NULL,
      mk_df(),
      "gaussian",
      NULL,
      NULL,
      list(),
      NULL,
      100,
      1,
      "auto"
    ),
    "prebuilt IR only for the"
  )
})

# ---------------------------------------------------------------- #
# Exported ingest adapters                                          #
# ---------------------------------------------------------------- #

test_that("the three ingest adapters are exported and return fb_terms", {
  exports <- getNamespaceExports("flexyBayes")
  for (fn in c("fb_from_asreml", "fb_from_brms", "fb_from_greta")) {
    expect_true(fn %in% exports, info = fn)
  }
  df <- mk_df()
  expect_s3_class(fb_from_asreml(y ~ x, random = ~g, data = df), "fb_terms")
  expect_s3_class(
    suppressMessages(fb_from_brms(y ~ x + (1 | g), data = df)),
    "fb_terms"
  )
})

# ---------------------------------------------------------------- #
# End-to-end: brms grammar through fb() reaches a non-greta backend #
# ---------------------------------------------------------------- #

test_that("fb() fits a bar-grouped formula and (auto) reaches INLA", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  df <- mk_df(40)
  fit <- suppressMessages(
    fb(
      y ~ x + (1 | g),
      data = df,
      backend = "auto",
      verbose = FALSE,
      mcmc_verbose = FALSE
    )
  )
  expect_s3_class(fit, "flexybayes_inla")
  expect_identical(fit$extras$backend_decision$backend, "inla")
})
