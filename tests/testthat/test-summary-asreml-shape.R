# =============================================================================
# One summary object, both engines.
#
# `summary(fit)$varcomp` is the first thing an ASReml user subsets. Before
# 0.9.1 it was NULL on every fit, and the two engines returned two
# incomparable objects: INLA a bare four-slot list of its own tables, brms
# its own `brmssummary`. This file fixes the shared object's shape, the
# column contract, and the three properties that make the numbers in it
# mean what the printed banner says they mean:
#
#   * INLA's variance components reach the standard-deviation scale
#     through the precision MARGINAL, never by transforming a tabulated
#     posterior mean -- the transform is nonlinear, so the two answers
#     differ, and the difference grows exactly where the component is
#     interesting. The bounds use the quantile identity, which a
#     strictly monotone transform carries exactly.
#   * the `prior` column is a projection of prior_summary(), not a second
#     prior-string builder that could drift from it.
#   * `converge` is whatever the engine reports. A Laplace approximation
#     has no R-hat, and one is not invented for it.
#
# The stored `$extras$variance_comps` table is deliberately NOT renamed:
# its five column names are a broom contract (R/tidiers.R) read by three
# internal callers and six test files.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


# ---------------------------------------------------------------- #
# Fixtures                                                          #
# ---------------------------------------------------------------- #

.sas_data <- function(seed = 20260817L, n = 60L) {
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

.sas_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}

.sas_cache <- new.env(parent = emptyenv())

.sas_inla_fit <- function() {
  if (is.null(.sas_cache$inla)) {
    .sas_cache$inla <- suppressMessages(fb(
      y ~ f + x,
      random = ~g,
      data = .sas_data(),
      backend = "inla",
      aggregate = FALSE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  .sas_cache$inla
}

.sas_brms_fit <- function() {
  if (is.null(.sas_cache$brms)) {
    .sas_cache$brms <- suppressMessages(suppressWarnings(fb(
      y ~ f + x + (1 | g),
      data = .sas_data(),
      backend = "brms",
      n_samples = 200L,
      warmup = 200L,
      chains = 1L,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )))
  }
  .sas_cache$brms
}

.SAS_SLOTS <- c(
  "fixed", "varcomp", "random", "missing", "converge",
  "n_design", "n_observed", "na_action", "model", "engine", "call"
)

.SAS_VARCOMP_COLS <- c(
  "component", "estimate", "std.error", "conf.low", "conf.high",
  "prior", "note"
)

# The frequentist vocabulary the table must never grow. Each of these
# names implies an AI-REML iteration and a Wald test, neither of which
# happened.
.SAS_FORBIDDEN <- c("z.ratio", "bound", "%ch", "Pr(>")


# ---------------------------------------------------------------- #
# 1. Shape, on both engines                                         #
# ---------------------------------------------------------------- #

test_that("summary() returns the eleven-slot object on an INLA fit", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_inla_fit()
  captured <- utils::capture.output(res <- withVisible(summary(fit)))
  s <- res$value

  expect_false(res$visible)
  expect_true(length(captured) > 0L)
  expect_identical(class(s), c("summary.flexybayes", "list"))
  expect_identical(names(s), .SAS_SLOTS)
  expect_identical(s$engine, "inla")
})

test_that("summary() returns the eleven-slot object on a brms fit", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_brms_fit()
  captured <- utils::capture.output(res <- withVisible(summary(fit)))
  s <- res$value

  expect_false(res$visible)
  expect_true(length(captured) > 0L)
  expect_identical(class(s), c("summary.flexybayes", "list"))
  expect_identical(names(s), .SAS_SLOTS)
  expect_identical(s$engine, "brms")
})

test_that("$varcomp carries the frozen columns and rows on both engines", {
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fits <- list()
  if (requireNamespace("INLA", quietly = TRUE)) {
    fits$inla <- .sas_inla_fit()
  }
  if (requireNamespace("brms", quietly = TRUE)) {
    fits$brms <- .sas_brms_fit()
  }
  skip_if(length(fits) == 0L, "neither engine installed")

  for (nm in names(fits)) {
    s <- suppressMessages(
      utils::capture.output(out <- summary(fits[[nm]]))
    )
    expect_true(is.data.frame(out$varcomp), label = nm)
    expect_true(
      all(.SAS_VARCOMP_COLS %in% names(out$varcomp)),
      label = nm
    )
    # A model with a residual scale and a random term has at least two
    # rows to report.
    expect_gte(nrow(out$varcomp), 1L)
    expect_true("sigma" %in% out$varcomp$component, label = nm)
    expect_true("sd_g" %in% out$varcomp$component, label = nm)
  }
})

test_that("$missing is an empty typed frame and never NULL", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  s <- suppressMessages(
    utils::capture.output(out <- summary(.sas_inla_fit()))
  )
  expect_true(is.data.frame(out$missing))
  expect_identical(
    names(out$missing),
    c("row", "estimate", "std.error", "conf.low", "conf.high")
  )
})

test_that("the counts and the model line come off the fit, not a formula", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_inla_fit()
  s <- suppressMessages(utils::capture.output(out <- summary(fit)))

  expect_identical(out$n_design, 60L)
  expect_identical(out$n_observed, 60L)
  expect_identical(out$n_design, stats::nobs(fit))
  # Derived from the model representation: the INLA emit's own formula
  # says f(g, model = "iid", hyper = ...), which is the program and not
  # the model.
  expect_match(out$model, "G: g iid; R: units", fixed = TRUE)
  expect_false(grepl("hyper", out$model, fixed = TRUE))
})


# ---------------------------------------------------------------- #
# 2. Nothing in frequentist dress                                   #
# ---------------------------------------------------------------- #

test_that("no Wald vocabulary reaches the object or the print", {
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fits <- list()
  if (requireNamespace("INLA", quietly = TRUE)) {
    fits$inla <- .sas_inla_fit()
  }
  if (requireNamespace("brms", quietly = TRUE)) {
    fits$brms <- .sas_brms_fit()
  }
  skip_if(length(fits) == 0L, "neither engine installed")

  for (nm in names(fits)) {
    printed <- utils::capture.output(out <- summary(fits[[nm]]))
    joined <- paste(printed, collapse = "\n")
    for (tok in .SAS_FORBIDDEN) {
      expect_false(
        grepl(tok, joined, fixed = TRUE),
        label = paste0(nm, ": printed output carries '", tok, "'")
      )
      expect_false(
        any(grepl(tok, names(out$varcomp), fixed = TRUE)),
        label = paste0(nm, ": $varcomp column named '", tok, "'")
      )
      expect_false(
        any(grepl(tok, names(out$fixed), fixed = TRUE)),
        label = paste0(nm, ": $fixed column named '", tok, "'")
      )
    }
  }
})


# ---------------------------------------------------------------- #
# 3. The INLA components are marginal-transformed (A-p1, A-p4)      #
# ---------------------------------------------------------------- #

test_that("INLA SD components equal the marginal transform, not 1/sqrt(mean)", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_inla_fit()
  s <- suppressMessages(utils::capture.output(out <- summary(fit)))

  hp <- fit$inla$summary.hyperpar
  marg <- fit$inla$marginals.hyperpar
  nm <- "Precision for the Gaussian observations"
  expect_true(nm %in% rownames(hp))
  expect_true(nm %in% names(marg))

  from_marginal <- INLA::inla.emarginal(function(x) 1 / sqrt(x), marg[[nm]])
  reported <- out$varcomp$estimate[out$varcomp$component == "sigma"]
  expect_equal(reported, from_marginal, tolerance = 1e-6)

  # And it is genuinely a different number from the forbidden route.
  # 1 / sqrt() is nonlinear, so transforming the posterior mean precision
  # does not give the posterior mean standard deviation. If these two
  # ever agreed to 1e-6 the test above would prove nothing.
  wrong_route <- 1 / sqrt(hp[nm, "mean"])
  expect_gt(abs(reported - wrong_route), 1e-6)

  # The bounds come from the precision marginal's own quantiles, with
  # the probabilities reversed: 1 / sqrt() is strictly decreasing, so
  # the SD's lower bound is carried by the precision's upper quantile.
  # A monotone transform carries quantiles exactly, which is why the
  # bounds may be transformed and the mean above may not.
  prec_q <- INLA::inla.qmarginal(c(0.975, 0.025), marg[[nm]])
  row <- out$varcomp[out$varcomp$component == "sigma", ]
  expect_equal(row$conf.low, 1 / sqrt(prec_q[[1L]]), tolerance = 1e-10)
  expect_equal(row$conf.high, 1 / sqrt(prec_q[[2L]]), tolerance = 1e-10)
  expect_lt(row$conf.low, row$estimate)
  expect_lt(row$estimate, row$conf.high)
})

test_that("a wide-support precision marginal still yields a table", {
  # Regression guard. Reading the bounds off a numerically transformed
  # density (INLA::inla.tmarginal) fails on this fit: the transform is
  # evaluated at inla.qmarginal((1:2048)/2049, m), whose extreme-left
  # values fall below the marginal's own support on a precision spanning
  # five orders of magnitude, so 1 / sqrt() returns NaN and the call
  # errors. The fit here is an ordinary random intercept -- nothing
  # exotic -- and before the quantile identity replaced the transformed
  # density it took the whole fit down with it, not just the table.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE
  )
  set.seed(202L)
  d <- data.frame(
    y = stats::rnorm(60L),
    x = stats::rnorm(60L),
    g = factor(rep(seq_len(6L), length.out = 60L))
  )
  fit <- suppressMessages(fb(
    y ~ x + (1 | g),
    data = d, backend = "inla", aggregate = FALSE,
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_inla")

  vc <- fit$extras$variance_comps
  expect_true(is.data.frame(vc))
  expect_true(all(c("sigma", "sd_g") %in% vc$component))
  expect_true(all(is.finite(vc$estimate)))
  expect_true(all(is.finite(vc$q2.5)))
  expect_true(all(is.finite(vc$q97.5)))
  expect_true(all(vc$q2.5 > 0))
  expect_true(all(vc$q2.5 <= vc$q97.5))
})

test_that("a correlation hyperparameter is not put through the SD transform", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  # The response is pure noise, so the field genuinely does not identify
  # and the fit says so. That is correct behaviour and not this test's
  # subject, which is whether a correlation is left off the SD transform.
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_spatial_collapse_warning = TRUE
  )
  set.seed(4L)
  g <- expand.grid(
    row = factor(seq_len(8L)), col = factor(seq_len(6L)),
    KEEP.OUT.ATTRS = FALSE
  )
  g$y <- stats::rnorm(nrow(g))
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(row):ar1(col), data = g,
    backend = "inla", verbose = FALSE
  ))
  s <- suppressMessages(utils::capture.output(out <- summary(fit)))

  expect_true(all(c("rho_row", "rho_col") %in% out$varcomp$component))
  rho <- out$varcomp[out$varcomp$component == "rho_row", ]
  hp <- fit$inla$summary.hyperpar
  expect_equal(rho$estimate, hp["Rho for row_id", "mean"], tolerance = 1e-10)
  expect_equal(rho$std.error, hp["Rho for row_id", "sd"], tolerance = 1e-10)
  # A correlation lives on [-1, 1]; a botched SD transform would not.
  expect_gte(rho$conf.low, -1)
  expect_lte(rho$conf.high, 1)
})


# ---------------------------------------------------------------- #
# 4. The prior column is a projection (A-p2)                        #
# ---------------------------------------------------------------- #

test_that("the prior cell is the resolved prior, not a second rendering", {
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fits <- list()
  if (requireNamespace("INLA", quietly = TRUE)) {
    fits$inla <- .sas_inla_fit()
  }
  if (requireNamespace("brms", quietly = TRUE)) {
    fits$brms <- .sas_brms_fit()
  }
  skip_if(length(fits) == 0L, "neither engine installed")

  for (nm in names(fits)) {
    fit <- fits[[nm]]
    printed <- utils::capture.output(out <- summary(fit))
    ps <- prior_summary(fit)
    expect_identical(ps$kind, "fb_prior", label = nm)

    # Every spec prior_summary() resolves appears verbatim in the cell
    # for its own component. Same record, same formatter, one source.
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
    }
    expect_false(any(is.na(out$varcomp$prior)), label = nm)
  }
})


# ---------------------------------------------------------------- #
# 5. Convergence stays engine-native (A-p3)                         #
# ---------------------------------------------------------------- #

test_that("an INLA fit reports its own diagnostics and no R-hat", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  printed <- utils::capture.output(out <- summary(.sas_inla_fit()))
  cv <- out$converge

  expect_true(is.list(cv))
  expect_match(cv$engine, "Laplace")
  expect_true(all(
    c("mode_status", "mlik", "kld_max", "numerical_confirm") %in% names(cv)
  ))
  expect_equal(cv$mode_status, 0)
  expect_true(cv$numerical_confirm)
  # No fabricated sampler diagnostics: a deterministic approximation has
  # no chains to compare and no draws to count.
  expect_false(any(c("max_rhat", "min_ess_bulk") %in% names(cv)))
  expect_false(any(grepl("R-hat", printed, fixed = TRUE)))
})

test_that("a brms fit reports R-hat, both ESS flavours and divergences", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  printed <- utils::capture.output(out <- summary(.sas_brms_fit()))
  cv <- out$converge

  expect_true(all(
    c("max_rhat", "min_ess_bulk", "min_ess_tail", "n_divergent") %in%
      names(cv)
  ))
  expect_true(is.finite(cv$max_rhat))
  expect_true(is.finite(cv$min_ess_bulk))
  expect_true(is.finite(cv$min_ess_tail))
  expect_false(any(c("mode_status", "kld_max") %in% names(cv)))
})


# ---------------------------------------------------------------- #
# 6. The boundary-collapse display flag (E-p3)                      #
# ---------------------------------------------------------------- #
#
# One pinned heuristic: flagged when the component's 97.5% quantile on
# the SD scale sits below 1% of the posterior-median residual SD. The
# fixtures below are constructed rather than fitted, because the point is
# the rule and not any engine's arithmetic.

.sas_shell <- function(cls, vc, medians) {
  attr(vc, "posterior_median") <- medians
  structure(
    list(extras = list(variance_comps = vc)),
    class = c(cls, "flexybayes", "list")
  )
}

.sas_collapse_table <- function() {
  data.frame(
    component = c("sigma", "sd_g", "sd_block"),
    estimate = c(1.0, 0.5, 0.0011),
    sd = c(0.1, 0.05, 0.0007),
    q2.5 = c(0.8, 0.4, 0.0001),
    q97.5 = c(1.2, 0.6, 0.0060),
    stringsAsFactors = FALSE
  )
}

test_that("a component piled against zero is flagged, on either engine", {
  medians <- c(sigma = 1.0, sd_g = 0.5, sd_block = 0.001)
  for (cls in c("flexybayes_inla", "flexybayes_brms")) {
    shell <- .sas_shell(cls, .sas_collapse_table(), medians)
    out <- flexyBayes:::.fb_summary_varcomp(shell)
    expect_identical(out$note, c("", "", "collapsed"), label = cls)
  }
})

test_that("the flag is a comparison against the residual scale", {
  # Same component, ten times the residual SD: 0.006 is now above 1% of
  # 0.1, so the same row is not flagged. The rule reads the model, not a
  # fixed threshold on the number.
  medians <- c(sigma = 0.1, sd_g = 0.5, sd_block = 0.001)
  vc <- .sas_collapse_table()
  vc$estimate[1L] <- 0.1
  vc$q97.5[1L] <- 0.12
  shell <- .sas_shell("flexybayes_inla", vc, medians)
  expect_identical(
    flexyBayes:::.fb_summary_varcomp(shell)$note,
    c("", "", "")
  )
})

test_that("with no residual scale to compare against, nothing is flagged", {
  # A sectioned residual leaves no scalar sigma, so there is no reference
  # and every cell stays blank rather than being flagged against a
  # number the model does not have.
  vc <- .sas_collapse_table()[-1L, ]
  shell <- .sas_shell("flexybayes_brms", vc, c(sd_g = 0.5, sd_block = 0.001))
  expect_identical(flexyBayes:::.fb_summary_varcomp(shell)$note, c("", ""))
})

test_that("a correlation is never flagged as collapsed", {
  vc <- data.frame(
    component = c("sigma", "rho_row"),
    estimate = c(1.0, 0.002),
    sd = c(0.1, 0.4),
    q2.5 = c(0.8, -0.99),
    q97.5 = c(1.2, 0.005),
    stringsAsFactors = FALSE
  )
  shell <- .sas_shell(
    "flexybayes_inla", vc, c(sigma = 1.0, rho_row = 0.002)
  )
  expect_identical(flexyBayes:::.fb_summary_varcomp(shell)$note, c("", ""))
})


# ---------------------------------------------------------------- #
# 7. What the new object must not disturb                           #
# ---------------------------------------------------------------- #

test_that("$extras$variance_comps keeps its five broom column names", {
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fits <- list()
  if (requireNamespace("INLA", quietly = TRUE)) {
    fits$inla <- .sas_inla_fit()
  }
  if (requireNamespace("brms", quietly = TRUE)) {
    fits$brms <- .sas_brms_fit()
  }
  skip_if(length(fits) == 0L, "neither engine installed")

  for (nm in names(fits)) {
    vc <- fits[[nm]]$extras$variance_comps
    expect_identical(
      names(vc),
      c("component", "estimate", "sd", "q2.5", "q97.5"),
      label = nm
    )
  }
})

test_that("$extras$summary stays populated for the seven internal readers", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_inla_fit()
  expect_false(is.null(fit$extras$summary))
  expect_true(is.data.frame(fit$extras$summary$fixed))
  expect_gte(nrow(fit$extras$summary$fixed), 1L)
})

test_that("tidy() column names are unchanged by the new summary", {
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  cols <- c("term", "estimate", "std.error", "conf.low", "conf.high")

  if (requireNamespace("INLA", quietly = TRUE)) {
    expect_identical(names(generics::tidy(.sas_inla_fit())), cols)
  }
  if (requireNamespace("brms", quietly = TRUE)) {
    fit <- .sas_brms_fit()
    expect_identical(names(generics::tidy(fit)), cols)
    expect_identical(
      names(generics::tidy(fit, effects = "random")),
      cols
    )
  }
})


# ---------------------------------------------------------------- #
# 5. The conditional twelfth slot                                   #
# ---------------------------------------------------------------- #
#
# The 0.9.0 INLA summary returned the AR1 field's parameters as
# `$spatial_field`, and the unified object dropped it. The removal was
# silent -- the caller got NULL and whatever arithmetic followed -- and
# it broke the parameter-recovery table in the spatio-temporal vignette,
# which no test reached. The slot is engine-native and conditional: on a
# fit with a field it is there, on a fit without one it is absent rather
# than NULL-valued, so the eleven names above stay a fixed contract.

.sas_field_data <- function(seed = 20260818L, n_row = 8L, n_col = 8L) {
  set.seed(seed)
  d <- expand.grid(col = factor(seq_len(n_col)), row = factor(seq_len(n_row)))
  ar1 <- function(k, rho) rho^abs(outer(seq_len(k), seq_len(k), "-"))
  chol_field <- t(chol(
    1.5^2 * kronecker(ar1(n_row, 0.7), ar1(n_col, 0.5))
  ))
  d$yield <- 20 + as.numeric(chol_field %*% stats::rnorm(n_row * n_col)) +
    stats::rnorm(n_row * n_col, 0, 0.6)
  d
}

test_that("summary() carries $spatial_field on a fit with an AR1 field", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- suppressMessages(flexybayes(
    yield ~ 1,
    random = ~ ar1(row):ar1(col),
    data = .sas_field_data(),
    backend = "inla",
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  invisible(utils::capture.output(s <- summary(fit)))

  expect_identical(names(s), c(.SAS_SLOTS, "spatial_field"))
  fld <- s$spatial_field
  expect_true(is.data.frame(fld))
  expect_identical(names(fld), c("parameter", "median", "lower", "upper"))
  # Two correlations, the field SD and the nugget SD.
  expect_identical(nrow(fld), 4L)
  expect_true(all(vapply(
    fld[, c("median", "lower", "upper")], is.numeric, logical(1L)
  )))
  expect_true(all(fld$lower <= fld$upper))
})

test_that("a fit without a field has no $spatial_field slot at all", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  invisible(utils::capture.output(s <- summary(.sas_inla_fit())))
  expect_identical(names(s), .SAS_SLOTS)
  expect_false("spatial_field" %in% names(s))
})


# ---------------------------------------------------------------- #
# 6. A table that cannot be built warns; the posterior survives      #
# ---------------------------------------------------------------- #
#
# A summary table is a reading of the fit, not part of it, so a table
# that cannot be built must not destroy a posterior the engine computed
# successfully. emit_inla() wraps the variance-component build in a
# tryCatch, warns naming the cause, and returns the fit with an empty
# table.
#
# The failure is FORCED here rather than provoked, and that choice is
# load-bearing. The model that reaches it in the wild -- a factor on both
# the fixed and the random side, so confounded with itself -- does not
# reach it reliably: INLA's solve on such a model is not bit-reproducible
# between runs, and on ten runs with identical inputs the table came back
# empty three times and built fine seven times. Pinning the contract on
# that would be pinning a coin toss. Pinning it on a forced failure holds
# the behaviour that is actually being claimed: warn loudly, name the
# cause, keep the posterior, and say where the numbers still are.

test_that("an unbuildable variance-component table warns and keeps the fit", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("testthat", "3.2.0")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()

  testthat::local_mocked_bindings(
    .inla_variance_comps = function(...) {
      stop("zero non-NA points", call. = FALSE)
    },
    .package = "flexyBayes"
  )

  fit <- NULL
  expect_warning(
    fit <- suppressMessages(flexybayes(
      y ~ x,
      random = ~g,
      data = .sas_data(),
      backend = "inla",
      aggregate = FALSE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )),
    "variance-component table could not be built"
  )

  # The posterior is untouched by a summary that could not be rendered.
  expect_s3_class(fit, "flexybayes_inla")
  expect_false(is.null(fit$inla))
  expect_gt(nrow(fit$inla$summary.hyperpar), 0L)
  # And the warning told the truth about what is empty, and about where
  # the numbers still are.
  expect_null(fit$extras$variance_comps)
  expect_true(is.data.frame(fit$inla$summary.hyperpar))
})


# ---------------------------------------------------------------- #
# 8. The aggregated representation returns the same object          #
# ---------------------------------------------------------------- #
#
# Aggregation is the DEFAULT route for an ordinary Gaussian call with a
# random term, so `summary(fit)$varcomp` being NULL on an aggregated fit
# was not a corner case: it was the commonest fit the package produces.
# `summary.flexybayes_aggregated()` returned its own list of the
# aggregated posterior's raw pieces -- beta_means, sigma_means,
# tau_means -- and carried no variance-component table at all.
#
# The aggregated emits record their components as posterior MEANS only,
# so the table is rebuilt from the engine's own posterior: on the INLA
# route through the same precision-marginal transform the per-row route
# uses. That is what §3 above pins for a per-row fit and what the second
# test here pins for an aggregated one -- the two routes must agree on
# what "the posterior mean of the standard deviation" is.

.FB_RANDOM_COLS_EXPECTED <- c(
  "group", "level", "estimate", "std.error", "conf.low", "conf.high"
)

.sas_agg_data <- function(seed = 20260817L, n = 1000L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(sample(paste0("f", seq_len(5L)), n, replace = TRUE)),
    g = factor(sample(paste0("g", seq_len(20L)), n, replace = TRUE))
  )
  b_g <- stats::rnorm(20L, sd = 0.6)[as.integer(d$g)]
  d$y <- 1 + 0.8 * (d$f == "f2") + b_g + stats::rnorm(n, sd = 0.9)
  d
}

.sas_agg_fit <- function() {
  if (is.null(.sas_cache$agg)) {
    .sas_cache$agg <- suppressMessages(flexybayes(
      y ~ f,
      random = ~g,
      data = .sas_agg_data(),
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  .sas_cache$agg
}

test_that("an aggregated fit returns the same eleven-slot object", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_agg_fit()
  # The premise, asserted rather than assumed.
  expect_s3_class(fit, "flexybayes_aggregated")
  expect_identical(fit$exactness, "aggregated_exact")

  captured <- utils::capture.output(res <- withVisible(summary(fit)))
  s <- res$value

  expect_false(res$visible)
  expect_true(length(captured) > 0L)
  expect_identical(class(s), c("summary.flexybayes", "list"))
  expect_true(all(.SAS_SLOTS %in% names(s)))
  expect_identical(s$engine, "inla")
})

test_that("$varcomp on an aggregated fit is the frozen table, not NULL", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_agg_fit()
  invisible(utils::capture.output(s <- summary(fit)))

  expect_false(is.null(s$varcomp))
  expect_true(is.data.frame(s$varcomp))
  expect_identical(names(s$varcomp), .SAS_VARCOMP_COLS)
  expect_true("sigma" %in% s$varcomp$component)
  expect_true("sd_g" %in% s$varcomp$component)
  # Every cell of the interval columns is a number, because the INLA
  # route has the marginals to compute one from.
  expect_true(all(is.finite(s$varcomp$estimate)))
  expect_true(all(is.finite(s$varcomp$conf.low)))
  expect_true(all(is.finite(s$varcomp$conf.high)))
  # The prior column is the same projection the per-row route gets.
  expect_false(any(is.na(s$varcomp$prior)))
})

test_that("aggregated SD components use the marginal transform", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_agg_fit()
  invisible(utils::capture.output(s <- summary(fit)))

  hp <- fit$inla$summary.hyperpar
  row <- grep("^Precision for g$", rownames(hp))
  skip_if(length(row) != 1L, "no random-effect precision on this fit")

  naive <- sqrt(1 / hp[row, "mean"])
  reported <- s$varcomp$estimate[s$varcomp$component == "sd_g"]
  marginal <- INLA::inla.emarginal(
    function(x) 1 / sqrt(x),
    fit$inla$marginals.hyperpar[["Precision for g"]]
  )
  expect_equal(reported, marginal, tolerance = 1e-8)
  # And the naive transform of the tabulated mean is a different number,
  # so the assertion above is not vacuous.
  expect_false(isTRUE(all.equal(reported, naive, tolerance = 1e-6)))
})

test_that("the aggregated object keeps every other frozen promise", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_agg_fit()
  printed <- utils::capture.output(s <- summary(fit))
  joined <- paste(printed, collapse = "\n")

  # No Wald vocabulary, on the object or the print.
  for (tok in .SAS_FORBIDDEN) {
    expect_false(grepl(tok, joined, fixed = TRUE), label = tok)
    expect_false(any(grepl(tok, names(s$varcomp), fixed = TRUE)), label = tok)
  }
  # $missing is the typed zero-row frame, never NULL.
  expect_true(is.data.frame(s$missing))
  expect_identical(
    names(s$missing),
    c("row", "estimate", "std.error", "conf.low", "conf.high")
  )
  expect_identical(nrow(s$missing), 0L)
  # The counts are the PRE-aggregation design rows: aggregation is a
  # computational route, not a smaller dataset.
  expect_identical(s$n_design, 1000L)
  expect_identical(s$n_observed, 1000L)
  expect_identical(s$n_design, stats::nobs(fit))
  # $random is the same six-column table per grouping factor.
  expect_true(is.list(s$random))
  expect_true("g" %in% names(s$random))
  expect_identical(names(s$random$g), .FB_RANDOM_COLS_EXPECTED)
})

test_that("the aggregated summary names its representation and compression", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_agg_fit()
  printed <- utils::capture.output(s <- summary(fit))

  # The banner names the representation where a per-row fit names the
  # engine, and the compression line is on the header.
  expect_true(any(grepl(
    "[flexyBayes / aggregated-gaussian]", printed, fixed = TRUE
  )))
  expect_true(any(grepl("aggregation:.*->.*cells", printed)))
  expect_true(any(grepl("ratio.*:1", printed)))
  # And the model slot carries the same story, so a caller reading the
  # object rather than the print still learns the fit was aggregated.
  expect_match(s$model, "G: g iid", fixed = TRUE)
  expect_match(s$model, "aggregated: N = 1 000 rows ->", fixed = TRUE)
})

test_that("no numerical confirm is reported where none ran", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_agg_fit()
  printed <- utils::capture.output(s <- summary(fit))

  # The aggregated route records no `num_check`. isTRUE(NULL) is FALSE,
  # so the slot used to report a failed check on a fit that never had
  # one.
  expect_null(fit$num_check)
  expect_true(is.na(s$converge$numerical_confirm))
  expect_false(any(grepl("Numerical confirm", printed, fixed = TRUE)))
  # The diagnostics that DID happen are still reported.
  expect_identical(
    s$converge$engine,
    "INLA nested Laplace approximation"
  )
  expect_false(is.null(s$converge$mode_status))
})

test_that("the ASReml-hands accessors answer on an aggregated fit", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_agg_fit()

  fixed <- coef(fit)
  expect_true(is.numeric(fixed))
  expect_true(length(fixed) > 0L)
  expect_true("(Intercept)" %in% names(fixed))

  rand <- coef(fit, what = "random")
  expect_true(is.list(rand))
  expect_identical(names(rand$g), .FB_RANDOM_COLS_EXPECTED)
  expect_identical(rand, ranef(fit))

  miss <- coef(fit, what = "missing")
  expect_true(is.data.frame(miss))
  expect_identical(nrow(miss), 0L)

  expect_identical(stats::nobs(fit, type = "design"), 1000L)
  expect_identical(stats::nobs(fit, type = "observed"), 1000L)

  cls <- predict(fit, classify = "f")
  expect_s3_class(cls, "fb_predict_classify")
  expect_true(all(
    c("estimate", "std.error", "conf.low", "conf.high") %in% names(cls)
  ))
})

test_that("an aggregated count fit carries the table too", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  d <- .sas_agg_data()
  d$y <- stats::rpois(nrow(d), lambda = 3)
  fit <- suppressMessages(flexybayes(
    y ~ f,
    random = ~g,
    data = d,
    family = "poisson",
    backend = "inla",
    aggregate = TRUE,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")

  printed <- utils::capture.output(s <- summary(fit))
  expect_identical(class(s), c("summary.flexybayes", "list"))
  expect_true(all(.SAS_SLOTS %in% names(s)))
  expect_identical(names(s$varcomp), .SAS_VARCOMP_COLS)
  # A Poisson likelihood has no residual scale, so `sigma` has no row and
  # the random-effect standard deviation does.
  expect_true("sd_g" %in% s$varcomp$component)
  expect_false("sigma" %in% s$varcomp$component)
  expect_true(any(grepl(
    "[flexyBayes / aggregated-poisson]", printed, fixed = TRUE
  )))
})

test_that("a route that recorded only means says so rather than inventing", {
  # No live fixture reaches this branch: the aggregated INLA route always
  # carries marginals, and the only route that did not is a since-withdrawn
  # engine (see NEWS.md). The contract is still the one that matters -- a
  # point estimate with no interval is reported as a point estimate with
  # no interval.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .sas_silence()
  fit <- .sas_agg_fit()

  stripped <- fit
  stripped$inla <- NULL
  stripped$extras$variance_comps <- list(sigma = 0.9, tau = 0.6)
  stripped$extras$summary$sigma_means <- 0.9
  stripped$extras$summary$tau_means <- 0.6

  vc <- flexyBayes:::.fb_summary_varcomp(stripped)
  expect_identical(names(vc), .SAS_VARCOMP_COLS)
  expect_identical(vc$component, c("sigma", "sd_g"))
  expect_identical(vc$estimate, c(0.9, 0.6))
  expect_true(all(is.na(vc$std.error)))
  expect_true(all(is.na(vc$conf.low)))
  expect_true(all(is.na(vc$conf.high)))
  expect_identical(vc$note, rep("no interval recorded", 2L))
})
