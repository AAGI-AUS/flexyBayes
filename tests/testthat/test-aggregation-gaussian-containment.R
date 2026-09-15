# =============================================================================
# Gaussian models are not fitted on the aggregated representation.
#
# The aggregated Gaussian route composes the within-cell sum of squares with a
# fixed gamma prior on the residual precision. A stated residual prior, or the
# package's own bounded default, therefore never reached the fit, while
# prior_summary() reported it: two aggregated fits with residual bounds 0.05
# and 50 returned the same residual SD. Until that composition exists,
# Gaussian models are fitted per row: aggregate = "auto" takes the per-row
# route, and aggregate = TRUE or a streamed Gaussian fit refuses by name.
# Binomial and Poisson aggregation is unaffected, because the count families
# carry no residual prior.
#
# The marginal likelihood the engine reports on the aggregated representation
# is not the marginal likelihood of the rows: it omits constants that depend
# on the data and, under standardisation, scales by the number of cells
# instead of the number of rows. It is not reported on any aggregated fit.
# =============================================================================

# A design the aggregation plan accepts: factor-only fixed and random terms,
# 1000 rows compressing to 100 cells.
.containment_data <- function(n = 1000L, seed = 20260915L) {
  set.seed(seed)
  d <- data.frame(
    f = factor(sample(paste0("f", seq_len(5L)), n, replace = TRUE)),
    g = factor(sample(paste0("g", seq_len(20L)), n, replace = TRUE))
  )
  d$y <- 2 + stats::rnorm(n)
  d$count <- stats::rpois(n, 3)
  d
}

# ---------------------------------------------------------------- #
# 1. aggregate = TRUE on a Gaussian model refuses by name            #
# ---------------------------------------------------------------- #

test_that("aggregate = TRUE on a Gaussian model refuses by name", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .containment_data()

  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ f,
      random = ~g,
      data = d,
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_aggregation_residual_prior_not_composed")
  expect_s3_class(err, "flexybayes_aggregate_refusal")
  expect_match(conditionMessage(err), "per row", fixed = TRUE)
  expect_match(conditionMessage(err), "aggregate = FALSE", fixed = TRUE)
})

# ---------------------------------------------------------------- #
# 2. The plan never names the aggregated route for a Gaussian model  #
# ---------------------------------------------------------------- #

test_that("the plan of an eligible Gaussian model names the per-row route", {
  skip_if_not_installed("INLA")
  d <- .containment_data()
  p <- flexybayes(y ~ f, random = ~g, data = d, backend = "auto", plan = TRUE)

  expect_true(isTRUE(p$aggregation$eligible))
  expect_false(identical(p$path, "aggregated_inla"))
})

test_that("the plan of an eligible Poisson model still names the aggregated route", {
  skip_if_not_installed("INLA")
  d <- .containment_data()
  p <- flexybayes(
    count ~ f,
    random = ~g,
    data = d,
    family = "poisson",
    backend = "auto",
    plan = TRUE
  )

  expect_identical(p$path, "aggregated_inla")
  expect_identical(p$reason_code, "count_aggregated_eligible")
})

# ---------------------------------------------------------------- #
# 3. A streamed Gaussian fit refuses before any chunk is read        #
# ---------------------------------------------------------------- #

test_that("flexybayes_stream() refuses a Gaussian fit by name", {
  d <- .containment_data()
  expect_error(
    flexybayes_stream(y ~ f, random = ~g, source = d, family = "gaussian", verbose = FALSE),
    class = "flexybayes_refusal_aggregation_residual_prior_not_composed"
  )
})

test_that("a chunk source is never read when a Gaussian stream refuses", {
  d <- .containment_data()
  chunks_read <- 0L
  src <- function(i) {
    chunks_read <<- chunks_read + 1L
    if (i == 1L) d else NULL
  }
  expect_error(
    flexybayes_stream(y ~ f, random = ~g, source = src, family = "gaussian", verbose = FALSE),
    class = "flexybayes_refusal_aggregation_residual_prior_not_composed"
  )
  expect_identical(chunks_read, 0L)
})

test_that("flexybayes_stream(fit = FALSE) still returns Gaussian sufficient statistics", {
  d <- .containment_data()
  agg <- flexybayes_stream(
    y ~ f,
    random = ~g,
    source = d,
    family = "gaussian",
    fit = FALSE,
    verbose = FALSE
  )
  expect_s3_class(agg, "fb_aggregated")
  expect_identical(as.integer(agg$N), nrow(d))
})

# ---------------------------------------------------------------- #
# 4. Every call of the aggregated Gaussian emitter is guarded        #
# ---------------------------------------------------------------- #

# The aggregated Gaussian emitter is kept for an exact composition of the
# residual prior, and nothing calls it. A caller would fit a Gaussian model on
# the aggregated representation without that prior, so the call sites are
# counted rather than assumed.
test_that("no code path calls the aggregated Gaussian emitter", {
  dir_r <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(dir_r), "package sources not present")
  file_vec <- list.files(dir_r, pattern = "[.][Rr]$", full.names = TRUE)
  callers <- character(0)
  for (i1 in seq_along(file_vec)) {
    lines_local <- readLines(file_vec[i1], warn = FALSE)
    code_local <- sub("#.*$", "", lines_local)
    hit <- grep("emit_gaussian_aggregated(", code_local, fixed = TRUE)
    hit <- hit[!grepl("emit_gaussian_aggregated <- function", code_local[hit], fixed = TRUE)]
    hit <- hit[!grepl("\"[^\"]*emit_gaussian_aggregated\\(", code_local[hit])]
    if (length(hit)) {
      callers <- c(callers, rep(basename(file_vec[i1]), length(hit)))
    }
  }
  expect_identical(callers, character(0))
})

# ---------------------------------------------------------------- #
# 5. Both reason codes are registered                                #
# ---------------------------------------------------------------- #

test_that("the containment and evidence codes are registered", {
  codes <- fb_refusals()$reason_code
  expect_true("aggregation_residual_prior_not_composed" %in% codes)
  expect_true("aggregated_evidence_unavailable" %in% codes)
})

# ---------------------------------------------------------------- #
# 6. Live: the default route and the evidence line                   #
# ---------------------------------------------------------------- #

test_that("aggregate = \"auto\" fits an eligible Gaussian model per row", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .containment_data()

  fit <- suppressMessages(flexybayes(
    y ~ f,
    random = ~g,
    data = d,
    backend = "inla",
    verbose = FALSE
  ))
  expect_false(inherits(fit, "flexybayes_aggregated"))
  expect_false(identical(fit$extras$backend_decision$path, "aggregated_gaussian"))
  expect_true(is.finite(summary(fit)$converge$mlik))
})

test_that("an aggregated Poisson fit reports no marginal likelihood, and says why", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- .containment_data()

  fit <- suppressMessages(flexybayes(
    count ~ f,
    random = ~g,
    data = d,
    family = "poisson",
    backend = "inla",
    aggregate = TRUE,
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes_aggregated")
  s <- summary(fit)
  expect_true(is.na(s$converge$mlik))
  expect_identical(s$converge$evidence_unavailable, "aggregated_evidence_unavailable")

  printed <- paste(utils::capture.output(print(s)), collapse = "\n")
  expect_match(printed, "aggregated_evidence_unavailable", fixed = TRUE)
  expect_false(grepl("Marginal log-likelihood: -?[0-9]", printed))
})
