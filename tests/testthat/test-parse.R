# Tests for formula parsing

test_that("parse_fixed handles intercept-only", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  info <- flexyBayes:::.parse_fixed(y ~ 1, dat)
  expect_equal(info$response, "y")
  expect_true(info$intercept)
  expect_length(info$terms, 0)
})

test_that("parse_fixed handles no-intercept with factor", {
  dat <- data.frame(y = rnorm(20), env = factor(rep(1:4, 5)))
  info <- flexyBayes:::.parse_fixed(y ~ 0 + env, dat)
  expect_false(info$intercept)
  expect_length(info$terms, 1)
  expect_equal(info$terms[[1]]$type, "factor")
  expect_equal(info$terms[[1]]$var, "env")
  expect_equal(info$terms[[1]]$n_levels, 4)
})

test_that("parse_fixed handles continuous covariate", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  info <- flexyBayes:::.parse_fixed(y ~ x, dat)
  expect_true(info$intercept)
  expect_length(info$terms, 1)
  expect_equal(info$terms[[1]]$type, "continuous")
  expect_equal(info$terms[[1]]$var, "x")
})

test_that("parse_fixed handles factor interaction", {
  dat <- data.frame(
    y = rnorm(20),
    A = factor(rep(1:2, 10)),
    B = factor(rep(1:5, 4))
  )
  info <- flexyBayes:::.parse_fixed(y ~ A:B, dat)
  expect_length(info$terms, 1)
  expect_equal(info$terms[[1]]$type, "factor_interaction")
})

test_that("parse_fixed handles I() expressions", {
  dat <- data.frame(y = rnorm(10), x = rnorm(10))
  info <- flexyBayes:::.parse_fixed(y ~ I(x^2), dat)
  expect_length(info$terms, 1)
  expect_equal(info$terms[[1]]$type, "expression")
})

test_that("parse_fixed errors on missing response", {
  dat <- data.frame(x = rnorm(10))
  expect_error(
    flexyBayes:::.parse_fixed(y ~ x, dat),
    "not found"
  )
})

test_that("parse_formula handles simple random", {
  dat <- data.frame(geno = factor(rep(1:5, 4)))
  terms <- flexyBayes:::.parse_formula(~geno, dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "simple")
  expect_equal(terms[[1]]$var, "geno")
  expect_equal(terms[[1]]$var_n, 5)
})

test_that("parse_formula handles crossed random", {
  dat <- data.frame(
    geno = factor(rep(1:5, 4)),
    env = factor(rep(1:4, each = 5))
  )
  terms <- flexyBayes:::.parse_formula(~ geno + env, dat)
  expect_length(terms, 2)
  expect_equal(terms[[1]]$type, "simple")
  expect_equal(terms[[2]]$type, "simple")
})

test_that("parse_formula handles nested (colon)", {
  dat <- data.frame(block = factor(rep(1:3, 4)), rep = factor(rep(1:2, 6)))
  terms <- flexyBayes:::.parse_formula(~ block:rep, dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "nested")
})

test_that("parse_formula handles vm()", {
  dat <- data.frame(geno = factor(1:5))
  terms <- flexyBayes:::.parse_formula(~ vm(geno, Gmat), dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "vm")
  expect_equal(terms[[1]]$var, "geno")
  expect_equal(terms[[1]]$mat, "Gmat")
})

test_that("parse_formula handles at(env):geno", {
  dat <- data.frame(
    geno = factor(rep(1:5, 4)),
    env = factor(rep(1:4, each = 5))
  )
  terms <- flexyBayes:::.parse_formula(~ at(env):geno, dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "at_simple")
})

test_that("parse_formula handles us(env):id(geno)", {
  dat <- data.frame(
    geno = factor(rep(1:5, 4)),
    env = factor(rep(1:4, each = 5))
  )
  terms <- flexyBayes:::.parse_formula(~ us(env):id(geno), dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "us_gxe")
})

test_that("parse_formula handles fa(env,2):id(geno)", {
  dat <- data.frame(
    geno = factor(rep(1:5, 4)),
    env = factor(rep(1:4, each = 5))
  )
  terms <- flexyBayes:::.parse_formula(~ fa(env, 2):id(geno), dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "fa_gxe")
  expect_equal(terms[[1]]$k, 2L)
})

test_that("parse_formula errors on fa(env, 0)", {
  dat <- data.frame(
    geno = factor(rep(1:5, 4)),
    env = factor(rep(1:4, each = 5))
  )
  expect_error(
    flexyBayes:::.parse_formula(~ fa(env, 0):id(geno), dat),
    "k >= 1"
  )
})

test_that("parse_formula handles ar1(row):id(col)", {
  dat <- data.frame(row = factor(rep(1:5, 3)), col = factor(rep(1:3, each = 5)))
  terms <- flexyBayes:::.parse_formula(~ ar1(row):id(col), dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "ar1_spatial")
})

test_that("parse_formula handles spl(x)", {
  dat <- data.frame(x = rnorm(20))
  terms <- flexyBayes:::.parse_formula(~ spl(x), dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "spline")
})

test_that("parse_formula handles at(env):units (rcov)", {
  dat <- data.frame(env = factor(rep(1:4, each = 5)))
  terms <- flexyBayes:::.parse_formula(~ at(env):units, dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "at_units")
})

test_that("parse_formula handles units", {
  terms <- flexyBayes:::.parse_formula(~units, data.frame(x = 1))
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "units")
})

test_that("parse_formula handles dsum(~units|env)", {
  dat <- data.frame(env = factor(rep(1:3, 5)))
  terms <- flexyBayes:::.parse_formula(~ dsum(~ units | env), dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "at_units")
  expect_equal(terms[[1]]$var, "env")
})

test_that("parse_formula handles three-way combo", {
  dat <- data.frame(a = factor(1:2), b = factor(1:2), c = factor(1:2))
  terms <- flexyBayes:::.parse_formula(~ a:b:c, dat)
  expect_length(terms, 1)
  expect_equal(terms[[1]]$type, "combo")
  expect_length(terms[[1]]$vars, 3)
})

test_that("resolve_family works", {
  expect_equal(flexyBayes:::.resolve_family("gaussian", NULL)$link, "identity")
  expect_equal(flexyBayes:::.resolve_family("binomial", NULL)$link, "logit")
  expect_equal(
    flexyBayes:::.resolve_family("binomial", "probit")$link,
    "probit"
  )
  expect_equal(flexyBayes:::.resolve_family("poisson", NULL)$link, "log")
  expect_error(flexyBayes:::.resolve_family("foo", NULL), "Unsupported")
})

# ---------------------------------------------------------------- #
# Parse-time refusals route through the structured registry         #
#                                                                   #
# The early formula / spec validation in .parse_fixed() and the     #
# term walker raises registered `flexybayes_refusal_*` conditions   #
# rather than bare stop(), so downstream tooling can pattern-match   #
# on the class. The user-facing message text is preserved.          #
# ---------------------------------------------------------------- #

test_that("the parse-time refusal codes are registered + user-visible", {
  codes <- c(
    "formula_not_two_sided",
    "response_not_in_data",
    "fa_rank_invalid",
    "smooth_variable_not_in_data"
  )
  tab <- fb_refusals()
  expect_true(all(codes %in% tab$reason_code))
  expect_true(all(tab$since_version[tab$reason_code %in% codes] == "0.4.0"))
})

test_that("a one-sided formula raises formula_not_two_sided", {
  d <- data.frame(y = rnorm(5), x = rnorm(5))
  err <- tryCatch(flexyBayes:::.parse_fixed(~x, d), error = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_formula_not_two_sided")
  expect_s3_class(err, "flexybayes_refusal")
  expect_match(conditionMessage(err), "must be two-sided")
})

test_that("a response absent from data raises response_not_in_data", {
  d <- data.frame(y = rnorm(5), x = rnorm(5))
  err <- tryCatch(flexyBayes:::.parse_fixed(z ~ x, d), error = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_response_not_in_data")
  expect_s3_class(err, "flexybayes_refusal")
  expect_match(conditionMessage(err), "not found in data")
})

test_that("fa(x, k) with k < 1 raises fa_rank_invalid", {
  d <- data.frame(
    y = rnorm(20),
    env = factor(rep(letters[1:4], 5)),
    geno = factor(rep(1:5, each = 4))
  )
  err <- tryCatch(
    fb_from_asreml(fixed = y ~ 1, random = ~ fa(env, 0):geno, data = d),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_fa_rank_invalid")
  expect_s3_class(err, "flexybayes_refusal")
  expect_match(conditionMessage(err), "fa\\(\\) requires k >= 1")
})

test_that("fa(x, k) with k >= n_outer raises fa_rank_exceeds_dim", {
  # Four environments => an identifiable factor-analytic structure needs
  # k < 4. k = 4 (= n_outer) is an over-parameterised reparameterisation
  # of the unstructured form; k = 5 (> n_outer) leaves empty loading
  # columns. Both are refused by the data-aware preflight in .enrich().
  d <- data.frame(
    y = rnorm(20),
    env = factor(rep(letters[1:4], 5)),
    geno = factor(rep(1:5, each = 4))
  )
  for (kk in c(4L, 5L)) {
    f <- stats::as.formula(paste0("~ fa(env, ", kk, "):geno"))
    err <- tryCatch(
      flexyBayes:::.parse_formula(f, d),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_refusal_fa_rank_exceeds_dim")
    expect_s3_class(err, "flexybayes_refusal")
    expect_match(conditionMessage(err), "fa\\(\\) requires k <")
    expect_match(conditionMessage(err), "Use k <= 3")
  }

  # The boundary just below n_outer (k = 3 with 4 levels) is accepted and
  # enriched with the data-derived n_outer.
  terms <- flexyBayes:::.parse_formula(~ fa(env, 3):geno, d)
  fa_term <- Find(function(tm) identical(tm$type, "fa_gxe"), terms)
  expect_false(is.null(fa_term))
  expect_identical(fa_term$k, 3L)
  expect_identical(fa_term$n_outer, 4L)

  # since_version metadata is correctly recorded as the new code's release.
  tab <- fb_refusals()
  expect_identical(
    tab$since_version[tab$reason_code == "fa_rank_exceeds_dim"],
    "0.7.0"
  )
})

test_that("s() on a variable absent from data raises smooth_variable_not_in_data", {
  testthat::skip_if_not_installed("mgcv")
  d <- data.frame(y = rnorm(20), x = rnorm(20))
  err <- tryCatch(
    fb_from_asreml(fixed = y ~ 1, random = ~ s(missing_var), data = d),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_smooth_variable_not_in_data")
  expect_s3_class(err, "flexybayes_refusal")
  expect_match(conditionMessage(err), "not found in data")
})

test_that("tensor smooths te()/ti()/t2() are refused, not silently accepted", {
  d <- data.frame(y = rnorm(20), x = rnorm(20), z = rnorm(20))
  for (f in c("~ te(x, z)", "~ ti(x, z)", "~ t2(x, z)")) {
    err <- tryCatch(
      fb_from_asreml(fixed = y ~ 1, random = as.formula(f), data = d),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_refusal_tensor_smooth_unsupported")
    expect_s3_class(err, "flexybayes_refusal")
    expect_match(conditionMessage(err), "univariate penalised splines")
  }
})
