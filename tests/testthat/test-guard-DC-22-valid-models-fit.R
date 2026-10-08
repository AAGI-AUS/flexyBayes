# Class guard DC-22: a valid model is refused or stops with an error.
#
# Each model below is fitted by lme4 (or glm) on the same data, so it is a
# model the peers accept; flexybayes() on INLA must fit it too. The cases
# come from the class's register rows: a stats family object (DEF-059), a
# transformed response (DEF-103), a stated prior on a three-way random
# interaction (DEF-136, DEF-198) and a model with no random terms
# (DEF-173), with neighbours of each: an offset, a polynomial term, a
# factor-by-numeric interaction, a nested intercept, weights and the bar
# grammar.
#
# `.g22_open` lists the cases refused today, with their register rows; the
# list must stay exact.

.g22_data <- function() {
  set.seed(22)
  d <- expand.grid(a = factor(1:3), b = factor(1:2), cc = factor(1:2), rep = 1:5)
  d$g <- factor(rep(1:6, length.out = nrow(d)))
  d$x <- stats::rnorm(nrow(d))
  d$t <- stats::runif(nrow(d), 1, 3)
  d$y <- exp(1 + 0.3 * d$x + stats::rnorm(6L, 0, 0.3)[d$g] +
    stats::rnorm(nrow(d), 0, 0.2))
  d$n <- stats::rpois(nrow(d), d$t * 2)
  d$k <- stats::rbinom(nrow(d), 1L, 0.4)
  d
}

# Each case: the flexybayes() arguments and the peer call on the same data.
.g22_corpus <- function(d) {
  list(
    stats_binomial = list(
      fb = list(k ~ x, random = ~g, family = stats::binomial()),
      peer = quote(lme4::glmer(k ~ x + (1 | g), data = d, family = stats::binomial()))
    ),
    stats_poisson = list(
      fb = list(n ~ x, random = ~g, family = stats::poisson()),
      peer = quote(lme4::glmer(n ~ x + (1 | g), data = d, family = stats::poisson()))
    ),
    transformed_response = list(
      fb = list(log(y) ~ x, random = ~g),
      peer = quote(lme4::lmer(log(y) ~ x + (1 | g), data = d))
    ),
    three_way_stated_prior = list(
      fb = list(
        log(y) ~ x, random = ~ a:b:cc,
        prior = fb_prior(sd(group = "a:b:cc") ~ half_normal(scale = 1))
      ),
      peer = quote(lme4::lmer(log(y) ~ x + (1 | a:b:cc), data = d))
    ),
    no_random_terms = list(
      fb = list(log(y) ~ x + a),
      peer = quote(stats::lm(log(y) ~ x + a, data = d))
    ),
    offset = list(
      fb = list(n ~ x + offset(log(t)), random = ~g, family = "poisson"),
      peer = quote(lme4::glmer(
        n ~ x + offset(log(t)) + (1 | g), data = d, family = stats::poisson()
      ))
    ),
    polynomial = list(
      fb = list(log(y) ~ x + I(x^2), random = ~g),
      peer = quote(lme4::lmer(log(y) ~ x + I(x^2) + (1 | g), data = d))
    ),
    factor_by_numeric = list(
      fb = list(log(y) ~ a * x, random = ~g),
      peer = quote(lme4::lmer(log(y) ~ a * x + (1 | g), data = d))
    ),
    nested_intercept = list(
      fb = list(log(y) ~ x, random = ~ g + g:a),
      peer = quote(lme4::lmer(log(y) ~ x + (1 | g) + (1 | g:a), data = d))
    ),
    weights = list(
      fb = list(log(y) ~ x, random = ~g, weights = d$t),
      peer = quote(lme4::lmer(log(y) ~ x + (1 | g), data = d, weights = t))
    ),
    bar_grammar = list(
      fb = list(y ~ x + (1 | g)),
      peer = quote(lme4::lmer(y ~ x + (1 | g), data = d))
    ),
    bar_grammar_transformed = list(
      fb = list(log(y) ~ x + (1 | g)),
      peer = quote(lme4::lmer(log(y) ~ x + (1 | g), data = d))
    )
  )
}

.g22_open <- c(bar_grammar_transformed = "DEF-215")

test_that("DC-22 guard: every model the peers fit, flexybayes() fits on INLA", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  d <- .g22_data()
  corpus <- .g22_corpus(d)
  refused <- character(0)
  for (nm in names(corpus)) {
    peer <- tryCatch(
      suppressWarnings(suppressMessages(eval(corpus[[nm]]$peer))),
      error = function(e) e
    )
    expect_false(inherits(peer, "error"), label = paste(nm, "peer fit"))
    fit <- tryCatch(
      suppressWarnings(suppressMessages(do.call(
        flexybayes, c(corpus[[nm]]$fb, list(data = d, backend = "inla", verbose = FALSE))
      ))),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      refused[[nm]] <- substr(conditionMessage(fit), 1L, 100L)
    } else {
      expect_s3_class(fit, "flexybayes")
    }
  }
  expect_identical(
    setdiff(names(refused), names(.g22_open)),
    character(0),
    label = "valid models refused"
  )
  expect_identical(
    setdiff(names(.g22_open), names(refused)),
    character(0),
    label = "open-list cases that now fit (remove them from .g22_open)"
  )
})
