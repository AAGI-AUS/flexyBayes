# Deprecation conditions.
#
# flexyBayes signals deprecation with base R conditions. A deprecated export,
# argument or keyword raises a warning of class `flexybayes_deprecated`. The
# warning also carries `lifecycle_warning_deprecated`, the class the package
# raised before 0.13.0, and base R's `deprecatedWarning`. A defunct argument
# raises a `defunctError`. Every call warns: there is no once-per-session
# latch. Nothing here needs an engine.

withr::local_options(flexyBayes.deprecation_warnings = TRUE)

.dep_data <- function() {
  withr::with_seed(1L, data.frame(
    y = stats::rnorm(60L),
    x = stats::rnorm(60L),
    g = factor(rep(1:6, 10L)),
    env = factor(rep(c("a", "b"), 30L))
  ))
}

.dep_cov <- function() {
  K <- diag(3)
  dimnames(K) <- list(letters[1:3], letters[1:3])
  K
}

# A deferred-execution object built without brms: the code slot is a fixed
# string and nothing is proceeded into.
.dep_review <- function(d) {
  flexyBayes:::.new_flexybayes_review(
    code = "data {\n}",
    backend = "stan_via_brms",
    ir = flexyBayes:::.fb_from_asreml_impl(y ~ x, data = d),
    prior = NULL,
    data_name = "d",
    call = quote(flexybayes(y ~ x, data = d, review_code = TRUE)),
    seed = NULL,
    proceed_args = list()
  )
}

.dep_count <- function(expr) {
  n <- 0L
  withCallingHandlers(
    expr,
    flexybayes_deprecated = function(w) {
      n <<- n + 1L
      invokeRestart("muffleWarning")
    }
  )
  n
}

test_that("the class predicates warn with the documented class and still answer", {
  K <- .dep_cov()
  expect_warning(
    out <- is_fb_cov(fb_cov(K, type = "dense")),
    "inherits(x, \"fb_cov\")", fixed = TRUE, class = "flexybayes_deprecated"
  )
  expect_true(out)
  expect_warning(out <- is_fb_cov(K), class = "flexybayes_deprecated")
  expect_false(out)
  expect_warning(
    out <- is_fb_approx(fb_approx("low_rank_smooth")),
    "inherits(x, \"fb_approx\")", fixed = TRUE, class = "flexybayes_deprecated"
  )
  expect_true(out)
  expect_warning(
    out <- is_fb_engine(suppressWarnings(fb_engine("inla"))),
    "inherits(x, \"fb_engine\")", fixed = TRUE, class = "flexybayes_deprecated"
  )
  expect_true(out)
})

test_that("fb() warns, names flexybayes(), and returns what flexybayes() returns", {
  d <- .dep_data()
  expect_warning(
    via_fb <- fb(y ~ x + (1 | g), data = d, plan = TRUE, verbose = FALSE),
    "`flexybayes()`", fixed = TRUE, class = "flexybayes_deprecated"
  )
  direct <- flexybayes(y ~ x + (1 | g), data = d, plan = TRUE, verbose = FALSE)
  expect_s3_class(via_fb, "fb_plan")
  expect_identical(as.data.frame(via_fb), as.data.frame(direct))
  expect_identical(formals(fb), formals(flexybayes))
})

test_that("the 0.13.0 compatibility surfaces warn and keep their result", {
  d <- .dep_data()
  expect_warning(
    ir_a <- fb_from_asreml(y ~ x, random = ~ g, data = d),
    "fb_plan()", fixed = TRUE, class = "flexybayes_deprecated"
  )
  expect_s3_class(ir_a, "fb_terms")
  expect_warning(
    ir_b <- fb_from_brms(y ~ x + (1 | g), data = d),
    "fb_plan()", fixed = TRUE, class = "flexybayes_deprecated"
  )
  expect_s3_class(ir_b, "fb_terms")
  expect_warning(
    eng <- fb_engine("brms", chains = 2L),
    "directly to `flexybayes()`", fixed = TRUE,
    class = "flexybayes_deprecated"
  )
  expect_identical(eng$opts$chains, 2L)
  expect_warning(
    expect_error(
      fb_brms(y ~ x, data = d, backend = "inla"),
      class = "flexybayes_refusal_engine_pin_backend_conflict"
    ),
    "backend = \"brms\"", fixed = TRUE,
    class = "flexybayes_deprecated"
  )
})

test_that("fb_inla() warns before its own checks and names flexybayes(backend = \"inla\")", {
  d <- .dep_data()
  expect_warning(
    expect_error(
      fb_inla(y ~ x + (1 | g), data = d, backend = "brms"),
      class = "flexybayes_refusal_engine_pin_backend_conflict"
    ),
    "backend = \"inla\"", fixed = TRUE, class = "flexybayes_deprecated"
  )
  skip_if_not_installed("INLA")
  expect_warning(
    via_pin <- fb_inla(y ~ x + (1 | g), data = d, plan = TRUE, verbose = FALSE),
    class = "flexybayes_deprecated"
  )
  direct <- flexybayes(
    y ~ x + (1 | g), data = d, backend = "inla", plan = TRUE, verbose = FALSE
  )
  expect_identical(as.data.frame(via_pin), as.data.frame(direct))
})

test_that("cat_code() and proceed() warn with the documented class and still work", {
  d <- .dep_data()
  rev <- .dep_review(d)
  expect_warning(
    emitted <- utils::capture.output(returned <- cat_code(rev)),
    "return_code = TRUE", fixed = TRUE, class = "flexybayes_deprecated"
  )
  expect_identical(returned, rev$code)
  expect_identical(emitted, c("data {", "}"))

  rev$proceeded <- TRUE
  rev$fit <- "the cached fit"
  withr::local_options(flexyBayes.silence_review_cached_note = TRUE)
  expect_warning(
    out <- proceed(rev),
    "review_code = FALSE", fixed = TRUE, class = "flexybayes_deprecated"
  )
  expect_identical(out, "the cached fit")
})

test_that("a deprecation warning is raised on every call, not once per session", {
  K <- .dep_cov()
  expect_identical(.dep_count({
    is_fb_cov(K)
    is_fb_cov(K)
    is_fb_cov(K)
  }), 3L)
  dat <- data.frame(geno = factor(1:5))
  expect_identical(.dep_count({
    flexyBayes:::.parse_formula(~ vm(geno, chol = L), dat)
    flexyBayes:::.parse_formula(~ vm(geno, chol = L), dat)
  }), 2L)
})

test_that("the warning and the defunct error carry the documented class vectors", {
  w <- tryCatch(is_fb_cov(.dep_cov()), warning = identity)
  expect_identical(
    class(w),
    c("flexybayes_deprecated", "lifecycle_warning_deprecated",
      "deprecatedWarning", "warning", "condition")
  )
  e <- tryCatch(
    suppressWarnings(fb_from_asreml(
      y ~ env, rcov = ~ units, data = .dep_data()
    )),
    error = identity
  )
  expect_identical(
    class(e),
    c("lifecycle_error_deprecated", "defunctError", "error", "condition")
  )
  expect_match(conditionMessage(e), "`residual`", fixed = TRUE)
})

test_that("the defunct rcov argument refuses with a defunctError on every entry", {
  d <- .dep_data()
  expect_error(
    flexybayes(y ~ env, rcov = ~ units, data = d, backend = "inla"),
    "residual", class = "defunctError"
  )
  expect_error(
    suppressWarnings(fb_from_asreml(y ~ env, rcov = ~ units, data = d)),
    "residual", class = "defunctError"
  )
  expect_warning(
    expect_error(
      fb(y ~ env, rcov = ~ units, data = d, backend = "inla"),
      "residual", class = "defunctError"
    ),
    class = "flexybayes_deprecated"
  )
  expect_warning(
    expect_error(
      fb_inla(y ~ env, rcov = ~ units, data = d),
      "residual", class = "defunctError"
    ),
    class = "flexybayes_deprecated"
  )
})


test_that("the replacements named in the warnings raise no deprecation warning", {
  d <- .dep_data()
  expect_no_warning(
    flexybayes(y ~ x + (1 | g), data = d, plan = TRUE, verbose = FALSE),
    class = "flexybayes_deprecated"
  )
  expect_no_warning(
    inherits(fb_cov(.dep_cov(), type = "dense"), "fb_cov"),
    class = "flexybayes_deprecated"
  )
  expect_no_warning(
    utils::capture.output(print(.dep_review(d))),
    class = "flexybayes_deprecated"
  )
})
