# Class guard DC-17: a refusal raised as a bare error, without a
# condition class.
#
# A corpus of models the package cannot fit is sent through every entry
# route: flexybayes() on INLA, flexybayes() returning brms code (no
# sampling), flexybayes(aggregate = TRUE), fb_plan() and
# flexybayes_stream(fit = FALSE). Every error a route raises must carry
# the class `flexybayes_refusal` and a reason code. A route that returns
# instead of stopping is not judged here.
#
# The cases come from the class's register rows: an unsupported family
# (DEF-017), an mgcv smoother in the fixed part (DEF-076), a factor slope
# in bar grammar (DEF-079), the aggregated route's entry refusals (DEF-080,
# DEF-159), a response outside its family's domain (DEF-163) and the
# stream route's model-shape refusals (DEF-164).
#
# `.g17_open` lists the case and route pairs that still stop untyped, each
# with its register row. The guard fails on any untyped error not listed,
# and on a listed pair that is now typed, so the list is kept exact.

.g17_data <- function() {
  set.seed(17)
  d <- data.frame(
    g = factor(rep(seq_len(6L), each = 5L)),
    trt = factor(rep(c("a", "b", "c"), 10L)),
    x = factor(rep(c("e1", "e2", "e3"), 10L)),
    z = stats::rnorm(30L)
  )
  d$y <- 1 + stats::rnorm(30L)
  d$yneg <- d$y - 10
  d$n <- stats::rpois(30L, 3)
  d
}

.g17_corpus <- list(
  unknown_family = list(fixed = y ~ x, random = ~g, family = "no_such_family"),
  survival_family = list(fixed = y ~ x, random = ~g, family = "weibull"),
  smoother_fixed = list(fixed = y ~ s(z), random = ~g),
  factor_slope_bar = list(fixed = y ~ x + (trt | g)),
  lognormal_nonpositive = list(fixed = yneg ~ x, random = ~g, family = "lognormal"),
  missing_column = list(fixed = y ~ nosuch, random = ~g),
  missing_random_column = list(fixed = y ~ x, random = ~nosuch),
  fa_term = list(fixed = y ~ x, random = ~ fa(trt, 2):g),
  bad_backend = list(fixed = y ~ x, random = ~g, backend = "stan_direct"),
  bad_aggregate = list(fixed = y ~ x, random = ~g, aggregate = "foo"),
  bad_link = list(fixed = n ~ x, random = ~g, family = "poisson", link = "probit"),
  bad_prior = list(fixed = y ~ x, random = ~g, prior = "sd ~ banana(1)"),
  weights_wrong_length = list(fixed = y ~ x, random = ~g, weights = c(1, 2)),
  data_not_frame = list(fixed = y ~ x, random = ~g, data = "d")
)

.g17_routes <- function(d) {
  fit_route <- function(extra) {
    function(a) {
      do.call(flexybayes, utils::modifyList(c(list(data = d, verbose = FALSE), extra), a))
    }
  }
  list(
    inla = fit_route(list(backend = "inla")),
    brms_code = fit_route(list(backend = "brms", return_code = TRUE)),
    aggregate = fit_route(list(backend = "inla", aggregate = TRUE)),
    plan = function(a) {
      a$formula <- a$fixed
      a$fixed <- NULL
      do.call(fb_plan, utils::modifyList(list(data = d), a))
    },
    stream = function(a) {
      a$source <- if (is.null(a$data)) d else a$data
      a[c("data", "aggregate", "weights", "link")] <- NULL
      do.call(flexybayes_stream, utils::modifyList(list(verbose = FALSE, fit = FALSE), a))
    }
  )
}

# Case and route pairs that stop untyped today, with their register rows.
.g17_open <- c(
  "missing_column inla" = "DEF-212",
  "missing_column brms_code" = "DEF-212",
  "missing_column aggregate" = "DEF-212",
  "missing_random_column inla" = "DEF-212",
  "missing_random_column brms_code" = "DEF-212",
  "missing_random_column aggregate" = "DEF-216",
  "bad_aggregate inla" = "DEF-213",
  "bad_aggregate brms_code" = "DEF-213",
  "bad_aggregate aggregate" = "DEF-213",
  "bad_aggregate plan" = "DEF-213",
  "weights_wrong_length inla" = "DEF-214",
  "weights_wrong_length brms_code" = "DEF-214",
  "weights_wrong_length aggregate" = "DEF-214",
  "weights_wrong_length plan" = "DEF-214",
  "bad_prior inla" = "DEF-215",
  "bad_prior brms_code" = "DEF-215",
  "bad_prior aggregate" = "DEF-216",
  "factor_slope_bar stream" = "DEF-216",
  "fa_term stream" = "DEF-216",
  "data_not_frame stream" = "DEF-217"
)

test_that("DC-17 guard: every route refuses an unfittable model with a typed refusal", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  d <- .g17_data()
  routes <- .g17_routes(d)
  untyped <- character(0)
  refused <- 0L
  for (nm in names(.g17_corpus)) {
    for (r in names(routes)) {
      out <- tryCatch(
        suppressWarnings(suppressMessages(
          utils::capture.output(res <- routes[[r]](.g17_corpus[[nm]]))
        )),
        error = function(e) e
      )
      if (!inherits(out, "error")) next
      refused <- refused + 1L
      typed <- inherits(out, "flexybayes_refusal") &&
        is.character(out$reason_code) && nzchar(out$reason_code)
      if (!typed) {
        untyped[[paste(nm, r)]] <- paste0(
          class(out)[1L], ": ", substr(conditionMessage(out), 1L, 80L)
        )
      }
    }
  }
  expect_gt(refused, 40L)
  new_untyped <- untyped[setdiff(names(untyped), names(.g17_open))]
  now_typed <- setdiff(names(.g17_open), names(untyped))
  expect_identical(
    new_untyped,
    stats::setNames(character(0), character(0)),
    label = "untyped errors not in the open list"
  )
  expect_identical(
    now_typed,
    character(0),
    label = "open-list pairs that are now typed (remove them from .g17_open)"
  )
})
