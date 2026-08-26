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

test_that("a brms formula combined with ASReml random/residual refuses", {
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

# ---------------------------------------------------------------- #
# Exported ingest adapters                                          #
# ---------------------------------------------------------------- #

test_that("the two ingest adapters are exported and return fb_terms", {
  exports <- getNamespaceExports("flexyBayes")
  for (fn in c("fb_from_asreml", "fb_from_brms")) {
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
# End-to-end: brms grammar through fb() reaches an active backend    #
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
