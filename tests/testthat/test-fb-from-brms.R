# fb_from_brms(): the minimum subset (response + linear fixed effects +
# random intercepts (1 | g) / (1 || g)), the family + link
# resolution path, the weights-as-addition_term capture, and the
# fail-fast refusals for unsupported brms specials.

# ---------------------------------------------------------------- #
# Test fixture                                                     #
# ---------------------------------------------------------------- #

mk_brms_data <- function() {
  set.seed(2026L)
  n <- 40L
  data.frame(
    y = rnorm(n),
    y01 = as.integer(rbinom(n, 1, 0.5)),
    x = rnorm(n),
    g = factor(rep(1:5, length.out = n)),
    site = factor(rep(1:4, length.out = n)),
    lon = rnorm(n),
    lat = rnorm(n)
  )
}

# ---------------------------------------------------------------- #
# Minimal cases                                                    #
# ---------------------------------------------------------------- #

test_that("fb_from_brms() reads minimal gaussian fixed-only model as a `continuous` term", {
  d <- mk_brms_data()
  fb <- flexyBayes:::fb_from_brms(y ~ x, data = d)
  expect_s3_class(fb, "fb_terms")
  expect_identical(fb$response, "y")
  expect_identical(fb$family, "gaussian")
  expect_identical(fb$link, "identity")
  expect_identical(fb$source, "brms")
  expect_true(fb$intercept)
  expect_length(fb$fixed_terms, 1L)
  expect_identical(fb$fixed_terms[[1]]$type, "continuous")
  expect_identical(fb$fixed_terms[[1]]$var, "x")
  expect_length(fb$random_terms, 0L)
  expect_length(fb$residual_terms, 0L)
  expect_length(fb$addition_terms, 0L)
})

test_that("fb_from_brms() reads an intercept-only model as one fixed term", {
  d <- mk_brms_data()
  fb <- flexyBayes:::fb_from_brms(y ~ 1, data = d)
  expect_true(fb$intercept)
  expect_length(fb$fixed_terms, 0L)
})

test_that("fb_from_brms() reads factor fixed effect as a `factor` term", {
  d <- mk_brms_data()
  fb <- flexyBayes:::fb_from_brms(y ~ g, data = d)
  expect_length(fb$fixed_terms, 1L)
  expect_identical(fb$fixed_terms[[1]]$type, "factor")
  expect_identical(fb$fixed_terms[[1]]$n_levels, 5L)
})

# ---------------------------------------------------------------- #
# Random intercepts                                                #
# ---------------------------------------------------------------- #

test_that("fb_from_brms() reads `(1 | g)` random intercept as a `simple` term", {
  d <- mk_brms_data()
  fb <- flexyBayes:::fb_from_brms(y ~ x + (1 | g), data = d)
  expect_length(fb$random_terms, 1L)
  expect_identical(fb$random_terms[[1]]$type, "simple")
  expect_identical(fb$random_terms[[1]]$var, "g")
  expect_identical(fb$random_terms[[1]]$var_n, 5L)
})

test_that("fb_from_brms() reads `(1 as a `simple` term", {
  d <- mk_brms_data()
  fb <- flexyBayes:::fb_from_brms(y ~ x + (1 || g), data = d)
  expect_length(fb$random_terms, 1L)
  expect_identical(fb$random_terms[[1]]$type, "simple")
  expect_identical(fb$random_terms[[1]]$var, "g")
})

test_that("fb_from_brms() reads two bar terms as two random effects", {
  d <- mk_brms_data()
  fb <- flexyBayes:::fb_from_brms(y ~ x + (1 | g) + (1 | site), data = d)
  expect_length(fb$random_terms, 2L)
  vars <- vapply(fb$random_terms, function(r) r$var, character(1))
  expect_setequal(vars, c("g", "site"))
})

# ---------------------------------------------------------------- #
# Family + link                                                    #
# ---------------------------------------------------------------- #

test_that("fb_from_brms() gives binomial the logit link by default", {
  d <- mk_brms_data()
  fb <- flexyBayes:::fb_from_brms(y01 ~ x, data = d, family = "binomial")
  expect_identical(fb$family, "binomial")
  expect_identical(fb$link, "logit")
})

test_that("fb_from_brms() accepts a `family()` object", {
  d <- mk_brms_data()
  fb <- flexyBayes:::fb_from_brms(
    y01 ~ x,
    data = d,
    family = stats::binomial(link = "probit")
  )
  expect_identical(fb$family, "binomial")
  expect_identical(fb$link, "probit")
})

test_that("fb_from_brms() rejects an unsupported character family", {
  d <- mk_brms_data()
  expect_error(
    flexyBayes:::fb_from_brms(y ~ x, data = d, family = "weibull5p_hyp"),
    "Unsupported family"
  )
})

# ---------------------------------------------------------------- #
# Weights                                                          #
# ---------------------------------------------------------------- #

test_that("fb_from_brms() captures `weights` as an addition_term", {
  d <- mk_brms_data()
  w <- runif(nrow(d))
  fb <- flexyBayes:::fb_from_brms(y ~ x, data = d, weights = w)
  expect_length(fb$addition_terms, 1L)
  expect_identical(fb$addition_terms[[1]]$type, "weights")
  expect_identical(fb$addition_terms[[1]]$values, w)
})

test_that("fb_from_brms() rejects mismatched weights length", {
  d <- mk_brms_data()
  expect_error(
    flexyBayes:::fb_from_brms(y ~ x, data = d, weights = 1:5),
    "must be a numeric vector of length nrow"
  )
})

# ---------------------------------------------------------------- #
# Fail-fast refusals for unsupported v0.1 features                 #
# ---------------------------------------------------------------- #

test_that("fb_from_brms() carries correlated random slopes (x | g)", {
  d <- mk_brms_data()
  # Both bar forms are ordinary terms now: (x || g) is two independent
  # variances and (x | g) is the same two with the correlation freed.
  fb <- flexyBayes:::fb_from_brms(y ~ x + (x | g), data = d)
  expect_identical(fb$random_terms[[1]]$type, "simple_slope_cor")
  expect_identical(fb$random_terms[[1]]$slope_var, "x")
  expect_identical(fb$random_terms[[1]]$var, "g")
})

test_that("fb_from_brms() refuses smoothers s(), naming spl()", {
  # The same refusal the ASReml grammar raises, naming the flexyBayes
  # spelling.
  d <- mk_brms_data()
  err <- expect_error(
    flexyBayes:::fb_from_brms(y ~ s(x), data = d),
    class = "flexybayes_refusal_fixed_smoother_not_supported"
  )
  expect_match(conditionMessage(err), "spl(x)", fixed = TRUE)
})

test_that("fb_from_brms() refuses gp() for v0.1", {
  d <- mk_brms_data()
  expect_error(
    flexyBayes:::fb_from_brms(y ~ gp(lon, lat), data = d),
    class = "flexybayes_refusal_brms_ingest_feature_unsupported"
  )
})

test_that("fb_from_brms() refuses ar() for v0.1", {
  d <- mk_brms_data()
  expect_error(
    flexyBayes:::fb_from_brms(y ~ x + ar(p = 1, gr = g), data = d),
    class = "flexybayes_refusal_brms_ingest_feature_unsupported"
  )
})

test_that("fb_from_brms() refuses cens() for v0.1", {
  d <- mk_brms_data()
  expect_error(
    flexyBayes:::fb_from_brms(y ~ x + cens(g), data = d),
    class = "flexybayes_refusal_brms_ingest_feature_unsupported"
  )
})

# ---------------------------------------------------------------- #
# Validator passes the produced object                             #
# ---------------------------------------------------------------- #

test_that("fb_from_brms() result passes validate_fb_terms()", {
  d <- mk_brms_data()
  fb <- flexyBayes:::fb_from_brms(y ~ x + (1 | g), data = d)
  expect_silent(flexyBayes:::validate_fb_terms(fb))
})

test_that("fb_from_brms() rejects malformed input", {
  d <- mk_brms_data()
  expect_error(
    flexyBayes:::fb_from_brms(~x, data = d),
    "must be two-sided"
  )
  expect_error(
    flexyBayes:::fb_from_brms(y ~ x, data = d, family = list("not_family")),
    "`family` must be"
  )
})

test_that("fb_from_brms() rejects response not in data", {
  d <- mk_brms_data()
  expect_error(
    flexyBayes:::fb_from_brms(missing_y ~ x, data = d),
    class = "flexybayes_refusal_response_not_in_data"
  )
})
