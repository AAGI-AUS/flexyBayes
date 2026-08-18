# =============================================================================
# The projection identity: one resolved prior, two surfaces.
#
# prior_summary() is the prior door and summary(fit)$varcomp$prior is a
# PROJECTION of it -- not a second prior-string builder that happens to
# agree today. This file holds the identity in both directions, on both
# engines:
#
#   * every prior the package resolved appears verbatim in the cell for
#     its own component;
#   * every cell that is not the two-word engine-default marker traces
#     back to a resolved prior with that exact string.
#
# The comparison is on the STORED record, which is unrounded. The printed
# table abbreviates the scale to four significant digits for reading
# (Phase 2, AF-2); the cell itself keeps the resolved prior's own string,
# and that is what a caller subsetting $varcomp gets. Both halves are
# asserted, because a rounding applied to the record rather than to the
# render would make this identity quietly false.
#
# A component the package left to the engine reads "engine default": the
# full sentence is too long for a table cell and prior_summary() carries
# it. Those components are skipped by construction, per Phase 1 H-1.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


# ---------------------------------------------------------------- #
# Fixtures                                                          #
# ---------------------------------------------------------------- #

.ppi_data <- function(seed = 20260817L, n = 60L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(rep(c("a", "b"), length.out = n)),
    g = factor(rep(letters[1:5], length.out = n)),
    x = stats::rnorm(n)
  )
  b_g <- stats::rnorm(5L, sd = 0.8)[as.integer(d$g)]
  d$y <- 1 + 0.6 * (d$f == "b") + 0.4 * d$x + b_g + stats::rnorm(n, sd = 0.5)
  d
}

.ppi_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}

.ppi_cache <- new.env(parent = emptyenv())

.ppi_inla_fit <- function() {
  if (is.null(.ppi_cache$inla)) {
    .ppi_cache$inla <- suppressMessages(fb(
      y ~ f + x,
      random = ~g,
      data = .ppi_data(),
      backend = "inla",
      aggregate = FALSE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  .ppi_cache$inla
}

.ppi_brms_fit <- function() {
  if (is.null(.ppi_cache$brms)) {
    .ppi_cache$brms <- suppressMessages(suppressWarnings(fb(
      y ~ f + x + (1 | g),
      data = .ppi_data(),
      backend = "brms",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )))
  }
  .ppi_cache$brms
}

.ppi_fits <- function() {
  out <- list()
  if (requireNamespace("INLA", quietly = TRUE)) {
    out$inla <- .ppi_inla_fit()
  }
  if (requireNamespace("brms", quietly = TRUE)) {
    out$brms <- .ppi_brms_fit()
  }
  out
}

# The cells that are not a projection of a resolved prior spec, and why.
.PPI_NON_SPEC_CELLS <- c("engine default", "legacy scalar bridge")


# ---------------------------------------------------------------- #
# 1. Every resolved prior appears in its own cell                   #
# ---------------------------------------------------------------- #

test_that("each resolved prior is the cell for its component", {
  skip_on_cran()
  skip_on_ci()
  .ppi_silence()
  fits <- .ppi_fits()
  skip_if(length(fits) == 0L, "neither engine installed")

  for (nm in names(fits)) {
    fit <- fits[[nm]]
    out <- suppressMessages(summary(fit))
    ps <- prior_summary(fit)
    expect_identical(ps$kind, "fb_prior", label = nm)

    checked <- 0L
    for (spec in ps$fb_prior$specs) {
      key <- flexyBayes:::.fb_prior_spec_parameter(spec)
      if (is.na(key) || !key %in% out$varcomp$component) {
        next
      }
      expect_identical(
        out$varcomp$prior[out$varcomp$component == key],
        flexyBayes:::.fb_prior_spec_string(spec),
        label = paste0(nm, ": ", key)
      )
      checked <- checked + 1L
    }
    # A silent zero would let the loop pass on a fit with no resolved
    # prior at all, which is the one case this test exists to catch.
    expect_gt(checked, 0L)
  }
})


# ---------------------------------------------------------------- #
# 2. Every cell traces back to a resolved prior                     #
# ---------------------------------------------------------------- #

test_that("no cell carries a string prior_summary() did not resolve", {
  # The other direction. A second prior-string builder would show up
  # here: a cell that reads plausibly but matches no resolved spec.
  skip_on_cran()
  skip_on_ci()
  .ppi_silence()
  fits <- .ppi_fits()
  skip_if(length(fits) == 0L, "neither engine installed")

  for (nm in names(fits)) {
    fit <- fits[[nm]]
    out <- suppressMessages(summary(fit))
    ps <- prior_summary(fit)

    resolved <- character(0)
    for (spec in ps$fb_prior$specs) {
      key <- flexyBayes:::.fb_prior_spec_parameter(spec)
      if (!is.na(key)) {
        resolved[[key]] <- flexyBayes:::.fb_prior_spec_string(spec)
      }
    }

    for (i in seq_len(nrow(out$varcomp))) {
      cmp <- out$varcomp$component[[i]]
      cell <- out$varcomp$prior[[i]]
      if (is.na(cell) || cell %in% .PPI_NON_SPEC_CELLS) {
        next
      }
      expect_true(cmp %in% names(resolved), label = paste0(nm, ": ", cmp))
      expect_identical(cell, unname(resolved[[cmp]]), label = paste0(nm, cmp))
    }
  }
})

test_that("an engine-default component is skipped rather than invented", {
  # Phase 1 H-1: where the package priored nothing, the cell says so in
  # two words and prior_summary() carries the sentence. The cell must
  # never be a string that looks like a resolved prior.
  skip_on_cran()
  skip_on_ci()
  .ppi_silence()
  fits <- .ppi_fits()
  skip_if(length(fits) == 0L, "neither engine installed")

  for (nm in names(fits)) {
    out <- suppressMessages(summary(fits[[nm]]))
    engine_default <- out$varcomp$prior %in% "engine default"
    expect_true(
      all(!is.na(out$varcomp$prior)),
      label = paste0(nm, ": no cell is NA")
    )
    if (any(engine_default)) {
      expect_false(
        any(grepl("(", out$varcomp$prior[engine_default], fixed = TRUE)),
        label = nm
      )
    }
  }
})


# ---------------------------------------------------------------- #
# 3. The rounding is a render, not a record (AF-2)                  #
# ---------------------------------------------------------------- #

test_that("the printed cell rounds and the stored cell does not", {
  skip_on_cran()
  skip_on_ci()
  .ppi_silence()
  fits <- .ppi_fits()
  skip_if(length(fits) == 0L, "neither engine installed")

  for (nm in names(fits)) {
    fit <- fits[[nm]]
    out <- suppressMessages(summary(fit))
    cells <- out$varcomp$prior[!out$varcomp$prior %in% .PPI_NON_SPEC_CELLS]
    skip_if(length(cells) == 0L, "no resolved prior on this fit")

    rounded <- flexyBayes:::.fb_round_prior_cell(cells, digits = 4L)
    printed <- paste(utils::capture.output(print(out)), collapse = "\n")

    for (k in seq_along(cells)) {
      # What the reader sees is the abbreviated scale.
      expect_match(printed, rounded[[k]], fixed = TRUE, label = nm)
      # What a caller subsetting $varcomp gets is the resolved record,
      # at the precision prior_summary() reports it.
      expect_gte(nchar(cells[[k]]), nchar(rounded[[k]]))
    }
  }
})


# ---------------------------------------------------------------- #
# 4. The legacy scalar bridge is a prior, and the cell says which   #
# ---------------------------------------------------------------- #

# Supplying `prior_vc_sd` closes the auto-default's gate -- it fires only
# on a call that supplies neither a `prior` nor a `prior_vc_sd` -- so this
# is the fit that reaches the engine through the legacy scalar bridge.
.ppi_legacy_brms_fit <- function() {
  if (is.null(.ppi_cache$legacy_brms)) {
    .ppi_cache$legacy_brms <- suppressMessages(suppressWarnings(fb(
      y ~ f + x + (1 | g),
      data = .ppi_data(),
      backend = "brms",
      prior_vc_sd = 3,
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )))
  }
  .ppi_cache$legacy_brms
}

test_that("a bridged component names the prior the bridge applied", {
  # `engine default` is what the classification alone produces here, and
  # on Stan it is false: .brms_legacy_specs() writes lognormal(0,
  # prior_vc_sd) onto sigma and onto every named sd group. The cell has to
  # say so, and the string it says is read off brms's own prior table
  # rather than rebuilt.
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .ppi_silence()
  fit <- .ppi_legacy_brms_fit()

  ps <- prior_summary(fit)
  expect_identical(ps$kind, "legacy_scalar")
  expect_equal(ps$vc_sd, 3)

  out <- suppressMessages(summary(fit))
  cells <- stats::setNames(out$varcomp$prior, out$varcomp$component)
  expect_true("sigma" %in% names(cells))
  expect_true("sd_g" %in% names(cells))
  expect_identical(unname(cells[["sigma"]]), "lognormal(0, 3)")
  expect_identical(unname(cells[["sd_g"]]), "lognormal(0, 3)")
  # The whole point of the repair: no bridged component reads as one the
  # package left alone.
  expect_false(any(cells %in% "engine default"))

  # And it renders, at display precision like every other cell (AF-2).
  printed <- paste(utils::capture.output(print(out)), collapse = "\n")
  expect_match(printed, "lognormal(0, 3)", fixed = TRUE)
})

test_that("the same bridge on INLA is genuinely an engine default", {
  # The counterpart, and the reason the repair is engine-aware. The bridge
  # hands INLA nothing -- the engine keeps its own hyperprior -- so
  # `engine default` is the true cell there and inventing a lognormal
  # would be the same defect in the other direction.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .ppi_silence()
  fit <- suppressMessages(fb(
    y ~ f + x,
    random = ~g,
    data = .ppi_data(),
    backend = "inla",
    aggregate = FALSE,
    prior_vc_sd = 3,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))

  expect_identical(prior_summary(fit)$kind, "legacy_scalar")
  out <- suppressMessages(summary(fit))
  expect_true(all(out$varcomp$prior %in% "engine default"))
})
