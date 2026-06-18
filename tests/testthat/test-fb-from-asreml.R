# Tests for fb_from_asreml() — phase 0.B of deliverable 0.
#
# Verifies that the existing asreml-format call shapes accepted by
# flexybayes() produce a valid fb_terms object via the new IR
# (intermediate representation) wrapper, with the parsed information
# preserved exactly. Round-trip transparency is the contract.
#
# Internal helpers accessed via `flexyBayes:::` (matches the pattern
# used in test-parse.R).

# ---------------------------------------------------------------- #
# Test fixture                                                     #
# ---------------------------------------------------------------- #

mk_asreml_data <- function() {
  set.seed(123)
  n <- 40L
  data.frame(
    yield = rnorm(n),
    env = factor(rep(1:4, length.out = n)),
    geno = factor(rep(1:5, length.out = n)),
    row = factor(rep(1:5, length.out = n)),
    col = factor(rep(1:8, length.out = n)),
    x = rnorm(n)
  )
}

# ---------------------------------------------------------------- #
# Minimal cases                                                    #
# ---------------------------------------------------------------- #

test_that("fb_from_asreml() handles a minimal gaussian fixed-only model", {
  d <- mk_asreml_data()
  fb <- flexyBayes:::fb_from_asreml(yield ~ x, data = d)

  expect_s3_class(fb, "fb_terms")
  expect_identical(fb$response, "yield")
  expect_identical(fb$family, "gaussian")
  expect_identical(fb$link, "identity")
  expect_true(fb$intercept)
  expect_length(fb$fixed_terms, 1L)
  expect_identical(fb$fixed_terms[[1]]$type, "continuous")
  expect_identical(fb$fixed_terms[[1]]$var, "x")
  expect_length(fb$random_terms, 0L)
  expect_length(fb$rcov_terms, 1L)
  expect_identical(fb$rcov_terms[[1]]$type, "units")
  expect_length(fb$addition_terms, 0L)
  expect_identical(fb$source, "asreml")
  expect_identical(fb$data_summary$n, nrow(d))
})

test_that("fb_from_asreml() default rcov is ~ units (matches flexybayes)", {
  d <- mk_asreml_data()
  fb <- flexyBayes:::fb_from_asreml(yield ~ x, data = d)
  expect_identical(fb$rcov_terms[[1]]$type, "units")
})

test_that("fb_from_asreml() handles intercept-only model", {
  d <- mk_asreml_data()
  fb <- flexyBayes:::fb_from_asreml(yield ~ 1, data = d)
  expect_true(fb$intercept)
  expect_length(fb$fixed_terms, 0L)
})

test_that("fb_from_asreml() handles no-intercept factor model", {
  d <- mk_asreml_data()
  fb <- flexyBayes:::fb_from_asreml(yield ~ 0 + env, data = d)
  expect_false(fb$intercept)
  expect_length(fb$fixed_terms, 1L)
  expect_identical(fb$fixed_terms[[1]]$type, "factor")
  expect_identical(fb$fixed_terms[[1]]$n_levels, 4L)
})

# ---------------------------------------------------------------- #
# Random effects                                                   #
# ---------------------------------------------------------------- #

test_that("fb_from_asreml() handles simple random effect", {
  d <- mk_asreml_data()
  fb <- flexyBayes:::fb_from_asreml(yield ~ env, random = ~geno, data = d)
  expect_length(fb$random_terms, 1L)
  expect_identical(fb$random_terms[[1]]$type, "simple")
  expect_identical(fb$random_terms[[1]]$var, "geno")
  expect_identical(fb$random_terms[[1]]$var_n, 5L)
})

test_that("fb_from_asreml() handles structured GxE (fa_gxe)", {
  d <- mk_asreml_data()
  fb <- flexyBayes:::fb_from_asreml(
    yield ~ env,
    random = ~ fa(env, 2):id(geno),
    data = d
  )
  expect_length(fb$random_terms, 1L)
  expect_identical(fb$random_terms[[1]]$type, "fa_gxe")
  expect_identical(fb$random_terms[[1]]$k, 2L)
})

test_that("fb_from_asreml() handles us_gxe", {
  d <- mk_asreml_data()
  fb <- flexyBayes:::fb_from_asreml(
    yield ~ env,
    random = ~ us(env):id(geno),
    data = d
  )
  expect_identical(fb$random_terms[[1]]$type, "us_gxe")
})

test_that("fb_from_asreml() handles ar1 spatial", {
  d <- mk_asreml_data()
  fb <- flexyBayes:::fb_from_asreml(
    yield ~ env,
    random = ~ ar1(row):id(col),
    data = d
  )
  expect_identical(fb$random_terms[[1]]$type, "ar1_spatial")
})

# ---------------------------------------------------------------- #
# Heterogeneous residual                                           #
# ---------------------------------------------------------------- #

test_that("fb_from_asreml() handles heterogeneous residual (at_units)", {
  d <- mk_asreml_data()
  fb <- flexyBayes:::fb_from_asreml(
    yield ~ env,
    random = ~geno,
    rcov = ~ at(env):units,
    data = d
  )
  expect_length(fb$rcov_terms, 1L)
  expect_identical(fb$rcov_terms[[1]]$type, "at_units")
})

# ---------------------------------------------------------------- #
# Family + link                                                    #
# ---------------------------------------------------------------- #

test_that("fb_from_asreml() handles binomial with default logit link", {
  d <- mk_asreml_data()
  d$y_bin <- as.integer(d$yield > 0)
  fb <- flexyBayes:::fb_from_asreml(y_bin ~ x, data = d, family = "binomial")
  expect_identical(fb$family, "binomial")
  expect_identical(fb$link, "logit")
})

test_that("fb_from_asreml() handles binomial with custom probit link", {
  d <- mk_asreml_data()
  d$y_bin <- as.integer(d$yield > 0)
  fb <- flexyBayes:::fb_from_asreml(
    y_bin ~ x,
    data = d,
    family = "binomial",
    link = "probit"
  )
  expect_identical(fb$family, "binomial")
  expect_identical(fb$link, "probit")
})

test_that("fb_from_asreml() handles poisson", {
  d <- mk_asreml_data()
  d$y_count <- pmax(0L, as.integer(round(d$yield + 5)))
  fb <- flexyBayes:::fb_from_asreml(y_count ~ x, data = d, family = "poisson")
  expect_identical(fb$family, "poisson")
  expect_identical(fb$link, "log")
})

test_that("fb_from_asreml() rejects unsupported family", {
  d <- mk_asreml_data()
  expect_error(
    flexyBayes:::fb_from_asreml(yield ~ x, data = d, family = "weibull"),
    "Unsupported family"
  )
})

# ---------------------------------------------------------------- #
# Weights, known_matrices, priors                                  #
# ---------------------------------------------------------------- #

test_that("fb_from_asreml() captures `weights` as an addition_term", {
  d <- mk_asreml_data()
  w <- runif(nrow(d))
  fb <- flexyBayes:::fb_from_asreml(yield ~ x, data = d, weights = w)
  expect_length(fb$addition_terms, 1L)
  expect_identical(fb$addition_terms[[1]]$type, "weights")
  expect_identical(fb$addition_terms[[1]]$values, w)
})

test_that("fb_from_asreml() rejects mismatched weights length", {
  d <- mk_asreml_data()
  expect_error(
    flexyBayes:::fb_from_asreml(yield ~ x, data = d, weights = 1:5),
    "must be a numeric vector of length nrow"
  )
})

test_that("fb_from_asreml() records known_matrices names in data_summary", {
  d <- mk_asreml_data()
  G <- diag(5)
  rownames(G) <- colnames(G) <- as.character(1:5)
  fb <- flexyBayes:::fb_from_asreml(
    yield ~ env,
    random = ~ vm(geno, Gmat),
    data = d,
    known_matrices = list(Gmat = G)
  )
  expect_identical(fb$random_terms[[1]]$type, "vm")
  expect_identical(fb$random_terms[[1]]$mat, "Gmat")
  expect_identical(fb$data_summary$known_matrices, "Gmat")
})

test_that("fb_from_asreml() captures legacy prior scalars", {
  d <- mk_asreml_data()
  fb <- flexyBayes:::fb_from_asreml(
    yield ~ x,
    data = d,
    prior_fixed_sd = 5,
    prior_vc_sd = 0.5
  )
  expect_true(isTRUE(fb$priors$legacy))
  expect_identical(fb$priors$fixed_sd, 5)
  expect_identical(fb$priors$vc_sd, 0.5)
})

# ---------------------------------------------------------------- #
# Composite round-trip + validator                                 #
# ---------------------------------------------------------------- #

test_that("fb_from_asreml() result passes validate_fb_terms() on a complex model", {
  d <- mk_asreml_data()
  fb <- flexyBayes:::fb_from_asreml(
    yield ~ env + x,
    random = ~ geno + fa(env, 2):id(geno),
    rcov = ~ at(env):units,
    data = d
  )
  expect_silent(flexyBayes:::validate_fb_terms(fb))
  expect_length(fb$fixed_terms, 2L)
  expect_length(fb$random_terms, 2L)
  expect_identical(fb$rcov_terms[[1]]$type, "at_units")
  expect_identical(fb$source, "asreml")
})

test_that("fb_from_asreml() preserves parse_formula.R term descriptors verbatim", {
  d <- mk_asreml_data()
  fb <- flexyBayes:::fb_from_asreml(
    yield ~ env,
    random = ~geno,
    data = d
  )

  # Compare against direct parse_formula calls — they should match
  # exactly (round-trip transparency).
  fixed_direct <- flexyBayes:::.parse_fixed(yield ~ env, d)
  random_direct <- flexyBayes:::.parse_formula(~geno, d)

  expect_identical(fb$fixed_terms, fixed_direct$terms)
  expect_identical(fb$intercept, fixed_direct$intercept)
  expect_identical(fb$response, fixed_direct$response)
  expect_identical(fb$random_terms, random_direct)
})
