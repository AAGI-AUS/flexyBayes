# =============================================================================
# What the prints say, and what they must never say.
#
# print() always opened on `MCMC : n chain(s) x n samples`, on every
# engine. On INLA that is false -- a nested Laplace approximation runs no
# chains and discards no warmup -- and to an ASReml user it reads as
# "this will be slow and random", which is the opposite of what the
# engine is. The three prints now share one header, and the sampler lines
# appear only where a sampler ran.
#
# The summary banner is the other half. The variance-component table
# looks like an ASReml variance-component table and is not one: every
# number in it is a posterior summary under a prior. The banner says so
# above the table, once, in one line.
# =============================================================================

suppressPackageStartupMessages(library(testthat))

.pb_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    # The fixture below is a 12-cell grid of pure noise with two plots
    # lost, so its autoregressive field genuinely does not identify and
    # the fit says so. That is correct behaviour and it is asserted in
    # test-spatial-field-collapse.R; this file's subject is the shape of
    # the printed banner, so the warning is silenced rather than allowed
    # to stand in for a failure here.
    flexyBayes.silence_spatial_collapse_warning = TRUE,
    .local_envir = parent.frame()
  )
}

.pb_cache <- new.env(parent = emptyenv())

# A 12-cell field with two lost plots: the smallest fixture that makes
# the design count and the observation count different numbers.
.pb_holed_grid <- function() {
  set.seed(20260817L)
  g <- expand.grid(
    row = factor(seq_len(4L)), col = factor(seq_len(3L)),
    KEEP.OUT.ATTRS = FALSE
  )
  g$y <- stats::rnorm(nrow(g))
  g$y[c(2L, 7L)] <- NA
  g
}

.pb_field_fit <- function() {
  if (is.null(.pb_cache$field)) {
    .pb_cache$field <- suppressMessages(flexybayes(
      y ~ 1,
      random = ~ ar1(row):ar1(col),
      data = .pb_holed_grid(),
      backend = "inla",
      verbose = FALSE
    ))
  }
  .pb_cache$field
}

.pb_brms_fit <- function() {
  if (is.null(.pb_cache$brms)) {
    set.seed(3L)
    d <- data.frame(g = factor(rep(letters[1:5], 12L)))
    d$y <- stats::rnorm(60L) + as.integer(d$g) * 0.2
    .pb_cache$brms <- suppressMessages(suppressWarnings(fb(
      y ~ 1 + (1 | g),
      data = d,
      backend = "brms",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )))
  }
  .pb_cache$brms
}


# ---------------------------------------------------------------- #
# 1. The banner above the variance components                       #
# ---------------------------------------------------------------- #

.PB_BANNER <- paste0(
  "Estimate = posterior mean (not REML). ",
  "Prior column is the prior that was used."
)

test_that("the banner is one line, exactly, above the table", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pb_silence()
  out <- utils::capture.output(summary(.pb_field_fit()))

  hits <- which(out == .PB_BANNER)
  expect_length(hits, 1L)

  # It sits between the section rule and the first table row.
  heading <- grep("Variance components", out, fixed = TRUE)
  expect_length(heading, 1L)
  expect_equal(hits[[1L]], heading[[1L]] + 1L)

  # The interval is named for what it is, in the heading right above.
  expect_match(out[[heading[[1L]]]], "credible interval", fixed = TRUE)
  expect_match(
    paste(out, collapse = "\n"), "posterior mean", fixed = TRUE
  )
})

test_that("a field fit names its fourth parameter after the table", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pb_silence()
  out <- utils::capture.output(summary(.pb_field_fit()))
  joined <- paste(out, collapse = "\n")

  sentence <- paste0(
    "The field and the nugget are separate parameters: this is not ",
    "ASReml's nugget-free residual."
  )
  expect_true(any(out == sentence))
  expect_match(joined, "not ASReml", fixed = TRUE)
  expect_match(joined, "nugget-free", fixed = TRUE)

  # After the table it explains, not before it.
  expect_gt(
    which(out == sentence)[[1L]],
    grep("Variance components", out, fixed = TRUE)[[1L]]
  )
})

test_that("a fit with no field does not claim one", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .pb_silence()
  joined <- paste(
    utils::capture.output(summary(.pb_brms_fit())),
    collapse = "\n"
  )
  expect_false(grepl("nugget", joined, fixed = TRUE))
  # The banner is not conditional on the engine.
  expect_match(joined, "not REML", fixed = TRUE)
})


# ---------------------------------------------------------------- #
# 2. Engine truthfulness in print() (D-p2)                          #
# ---------------------------------------------------------------- #

test_that("an INLA print says Laplace and never MCMC or warmup", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pb_silence()
  for (out in list(
    utils::capture.output(print(.pb_field_fit())),
    utils::capture.output(summary(.pb_field_fit()))
  )) {
    joined <- paste(out, collapse = "\n")
    expect_match(joined, "Laplace", fixed = TRUE)
    expect_false(grepl("MCMC", joined, fixed = TRUE))
    expect_false(grepl("warmup", joined, fixed = TRUE))
    expect_false(grepl("chain", joined, fixed = TRUE))
  }
})

test_that("a brms print keeps its sampler settings and diagnostics", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .pb_silence()
  printed <- paste(
    utils::capture.output(print(.pb_brms_fit())),
    collapse = "\n"
  )
  expect_match(printed, "chain", fixed = TRUE)
  expect_match(printed, "warmup", fixed = TRUE)
  expect_match(printed, "Stan", fixed = TRUE)

  summarised <- paste(
    utils::capture.output(summary(.pb_brms_fit())),
    collapse = "\n"
  )
  expect_match(summarised, "R-hat", fixed = TRUE)
  expect_match(summarised, "ESS", fixed = TRUE)
})

test_that("the header derives the model line from the representation", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pb_silence()
  out <- utils::capture.output(print(.pb_field_fit()))
  model <- grep("  Model    :", out, fixed = TRUE, value = TRUE)
  expect_length(model, 1L)
  expect_match(model, "ar1(row):ar1(col) field + nugget", fixed = TRUE)
  expect_match(model, "4 parameters", fixed = TRUE)
  # Not read back off the emitted INLA formula, which spells the same
  # term as f(row_id, model = "ar1", group = col_id, ...).
  expect_false(grepl("row_id", model, fixed = TRUE))

  # The two models share a spelling and differ by a parameter, and the
  # print says so where the model is described.
  joined <- paste(out, collapse = "\n")
  expect_match(joined, "not ASReml residual", fixed = TRUE)
  expect_match(joined, "3 parameters", fixed = TRUE)
})

test_that("a fit with no field makes no claim about ASReml's residual", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .pb_silence()
  joined <- paste(
    utils::capture.output(print(.pb_brms_fit())),
    collapse = "\n"
  )
  expect_false(grepl("ASReml", joined, fixed = TRUE))
})


# ---------------------------------------------------------------- #
# 3. The two counts                                                 #
# ---------------------------------------------------------------- #

test_that("the N line separates design rows from observed responses", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pb_silence()
  out <- utils::capture.output(print(.pb_field_fit()))
  n_line <- grep("  N        :", out, fixed = TRUE, value = TRUE)
  expect_length(n_line, 1L)
  expect_match(n_line, "12 design rows, 10 observed responses", fixed = TRUE)

  na_line <- grep("  na_action:", out, fixed = TRUE, value = TRUE)
  expect_length(na_line, 1L)
  expect_match(na_line, "augment", fixed = TRUE)
  expect_match(na_line, "2 missing response(s)", fixed = TRUE)
})

test_that("nobs() defaults to the design count and offers the other", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pb_silence()
  fit <- .pb_field_fit()

  expect_identical(stats::nobs(fit), 12L)
  expect_identical(stats::nobs(fit, type = "design"), 12L)
  expect_identical(stats::nobs(fit, type = "observed"), 10L)
  # The acceptance the work package states: observed is design less the
  # responses that were not observed.
  expect_identical(
    stats::nobs(fit, type = "observed"),
    stats::nobs(fit) - fit$extras$na_action$n_missing_response
  )
  expect_error(stats::nobs(fit, type = "sown"), "should be one of")
})

test_that("on complete data the two counts agree", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pb_silence()
  g <- .pb_holed_grid()
  g$y[is.na(g$y)] <- 0
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(row):ar1(col), data = g,
    backend = "inla", verbose = FALSE
  ))
  expect_identical(stats::nobs(fit), stats::nobs(fit, type = "observed"))
})


# ---------------------------------------------------------------- #
# 4. The formula, and the prior, at reading precision               #
# ---------------------------------------------------------------- #

.pb_iid_fit <- function() {
  if (is.null(.pb_cache$iid)) {
    set.seed(41L)
    d <- data.frame(g = factor(rep(letters[1:5], 12L)))
    d$y <- stats::rnorm(60L) + as.integer(d$g) * 0.2
    .pb_cache$iid <- suppressMessages(flexybayes(
      y ~ 1, random = ~g, data = d, backend = "inla", aggregate = FALSE,
      verbose = FALSE, mcmc_verbose = FALSE
    ))
  }
  .pb_cache$iid
}

test_that("the printed INLA formula names its priors instead of coding", {
  # A uniform-on-SD prior reaches INLA as a C expression: about a hundred
  # characters of theta, -1.0e10 and return(...) spliced into the
  # formula. Printing it buries the model in the encoding of one prior.
  skip_if_not_installed("INLA")
  skip_on_cran()
  .pb_silence()
  fit <- .pb_iid_fit()

  # The premise: this fit's emitted formula really does carry the blob.
  full <- paste(deparse(fit$extras$formula), collapse = " ")
  expect_match(full, "hyper = list(", fixed = TRUE)
  expect_match(full, "-1.0e10", fixed = TRUE)

  txt <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_match(txt, "INLA formula", fixed = TRUE)
  expect_false(grepl("-1.0e10", txt, fixed = TRUE))
  expect_false(grepl("return(", txt, fixed = TRUE))
  expect_false(grepl("expression:", txt, fixed = TRUE))
  expect_false(grepl("theta", txt, fixed = TRUE))
  # The model itself is still legible, and the prior is named.
  expect_match(txt, "f(g, model = \"iid\"", fixed = TRUE)
  expect_match(txt, "<prior: uniform-SD(0, ", fixed = TRUE)

  # A field print still shows its own formula, blob-free.
  field_txt <- paste(
    utils::capture.output(print(.pb_field_fit())), collapse = "\n"
  )
  expect_match(field_txt, "f(row_id", fixed = TRUE)
  expect_false(grepl("-1.0e10", field_txt, fixed = TRUE))
})

test_that("a compressed prior block is tagged by what it encodes", {
  compress <- flexyBayes:::.fb_compress_hyper_blobs
  uniform0 <- paste0(
    "f(g, model = \"iid\", hyper = list(prec = list(prior = ",
    "\"expression: U=5.596432198; lb=-2*log(U); ",
    "ld=-log(U)-log(2)-theta/2; return( theta<lb ? -1.0e10 : ld );\")))"
  )
  out <- compress(uniform0)
  expect_match(out, "<prior: uniform-SD(0, 5.596)>", fixed = TRUE)
  expect_false(grepl("theta", out, fixed = TRUE))
  # Everything outside the block is untouched.
  expect_match(out, "f(g, model = \"iid\", ", fixed = TRUE)

  bounded <- sub(
    "expression: U=5.596432198;",
    "expression: L=0.25; U=5.596432198;",
    uniform0,
    fixed = TRUE
  )
  expect_match(
    compress(bounded), "<prior: uniform-SD(0.25, 5.596)>",
    fixed = TRUE
  )

  halfnorm <- sub(
    "expression: U=5.596432198; lb=-2*log(U);",
    "expression: s=1.234567; sig=exp(-theta/2);",
    uniform0,
    fixed = TRUE
  )
  expect_match(
    compress(halfnorm), "<prior: half-normal-SD(1.235)>",
    fixed = TRUE
  )

  named <- paste0(
    "f(g, model = \"iid\", hyper = list(prec = list(prior = \"pc.prec\", ",
    "param = c(1.0000001, 0.01))))"
  )
  expect_match(compress(named), "<prior: pc.prec(1, 0.01)>", fixed = TRUE)

  # Two blocks in one formula both compress.
  two <- paste(uniform0, "+", uniform0)
  expect_identical(
    lengths(regmatches(compress(two), gregexpr("<prior:", compress(two)))),
    2L
  )

  # A formula with no prior block is returned unchanged.
  plain <- "y ~ 1 + f(g, model = \"iid\")"
  expect_identical(compress(plain), plain)
})

test_that("the prior cell rounds for display and not in the object", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .pb_silence()
  fit <- .pb_field_fit()

  out <- utils::capture.output(s <- summary(fit))
  txt <- paste(out, collapse = "\n")

  cells <- s$varcomp$prior
  cells <- cells[!is.na(cells) & grepl("[0-9]", cells)]
  skip_if(length(cells) == 0L, "no numeric prior cell on this fixture")

  # The stored cell is the resolved prior verbatim -- the basis the
  # projection identity is checked on -- and carries more digits than the
  # printed one.
  expect_true(any(nchar(cells) > nchar(
    flexyBayes:::.fb_round_prior_cell(cells)
  )))
  expect_false(any(grepl("[0-9]{8}", flexyBayes:::.fb_round_prior_cell(cells))))
  # And the printed table shows the short form.
  for (short in flexyBayes:::.fb_round_prior_cell(cells)) {
    expect_match(txt, short, fixed = TRUE)
  }
})

test_that("rounding the prior cell leaves the family name alone", {
  round_cell <- flexyBayes:::.fb_round_prior_cell
  expect_identical(
    round_cell("uniform(lower=0, upper=4.08242179176)"),
    "uniform(lower=0, upper=4.082)"
  )
  expect_identical(round_cell("engine default"), "engine default")
  expect_identical(round_cell(character(0)), character(0))
  expect_identical(
    round_cell("half_normal(scale=1.23456789)"),
    "half_normal(scale=1.235)"
  )
})
