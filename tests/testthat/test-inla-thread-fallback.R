# A host may report no core count. INLA turns that into `NA:1`, emits a parser
# failure and silently chooses another plan. The package supplies `1:1` only
# for that case. Valid defaults and explicit caller settings are left alone.

test_that("an unknown host core count has a valid conservative INLA fallback", {
  expect_identical(.fb_inla_threads_if_cores_unknown(NA_integer_), "1:1")
  expect_identical(.fb_inla_threads_if_cores_unknown(integer()), "1:1")
  expect_identical(.fb_inla_threads_if_cores_unknown(0L), "1:1")
  expect_identical(.fb_inla_threads_if_cores_unknown(-1L), "1:1")
  expect_null(.fb_inla_threads_if_cores_unknown(1L))
  expect_null(.fb_inla_threads_if_cores_unknown(8L))
})

test_that("this host never hands INLA the malformed NA thread spelling", {
  value <- .fb_inla_threads_if_cores_unknown()
  if (is.null(value)) {
    expect_true(parallel::detectCores(logical = TRUE) >= 1L)
  } else {
    expect_identical(value, "1:1")
  }
  expect_false(identical(value, "NA:1"))
})

# The thread plan of every INLA fit. At any plan other than "1:1", identical
# calls on a model with a flat variance direction return different variance
# components. The package therefore runs INLA at "1:1" unless told otherwise,
# and reads INLA's own global option no more.

test_that("the thread plan is 1:1 by default, set by option, and an explicit value wins", {
  withr::local_options(flexyBayes.inla_threads = NULL)
  expect_identical(.fb_inla_thread_plan(), "1:1")
  withr::local_options(flexyBayes.inla_threads = "2:1")
  expect_identical(.fb_inla_thread_plan(), "2:1")
  expect_identical(.fb_inla_thread_plan("3:2"), "3:2")
  withr::local_options(flexyBayes.inla_threads = "inla")
  expect_identical(.fb_inla_thread_plan(), .fb_inla_threads_if_cores_unknown())
})

test_that("a malformed thread option stops with its value named", {
  for (bad in list("four", 4, "0:1", "4:", c("1:1", "2:1"), NA_character_)) {
    withr::local_options(flexyBayes.inla_threads = bad)
    expect_error(.fb_inla_thread_plan(), "flexyBayes.inla_threads must be")
  }
})

.th_data <- function(seed = 20260926L, n = 120L) {
  set.seed(seed)
  d <- data.frame(
    g = factor(rep(sprintf("g%02d", 1:10), length.out = n)),
    h = factor(rep(c("A", "B", "C"), length.out = n))
  )
  d$y <- 1 +
    c(0, 0.4, -0.3)[as.integer(d$h)] +
    stats::rnorm(10L, sd = 0.8)[as.integer(d$g)] +
    stats::rnorm(n, sd = 0.6)
  d$cnt <- stats::rpois(
    n,
    exp(0.5 + stats::rnorm(10L, sd = 0.4)[as.integer(d$g)])
  )
  d
}

# The aggregated Gaussian emitter is refused on every user route while it
# cannot carry the residual prior, so no fit reaches it here.
.th_fit <- function(d, family = "gaussian", route = "per_row") {
  response <- if (identical(family, "gaussian")) "y" else "cnt"
  fixed <- stats::as.formula(paste(response, "~ h"))
  suppressMessages(suppressWarnings(switch(
    route,
    per_row = flexybayes(
      fixed,
      random = ~g,
      data = d,
      family = family,
      backend = "inla",
      aggregate = FALSE,
      verbose = FALSE
    ),
    aggregated = flexybayes(
      fixed,
      random = ~g,
      data = d,
      family = family,
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    ),
    stream = flexybayes_stream(
      fixed,
      random = ~g,
      source = d,
      family = family,
      backend = "inla",
      chunk_rows = 50000,
      verbose = FALSE
    )
  )))
}

test_that("every INLA route runs at the package plan, whatever INLA's own option says", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = NULL
  )
  old <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "4:1")
  withr::defer(INLA::inla.setOption(num.threads = old))
  d <- .th_data()
  for (route in list(
    c("gaussian", "per_row"),
    c("poisson", "per_row"),
    c("poisson", "aggregated"),
    c("poisson", "stream")
  )) {
    fit <- .th_fit(d, route[1L], route[2L])
    expect_identical(
      fit$extras$reproducibility$threads$num_threads,
      "1:1",
      label = paste(route, collapse = " ")
    )
  }
  withr::local_options(flexyBayes.inla_threads = "2:1")
  expect_identical(.th_fit(d)$extras$reproducibility$threads$num_threads, "2:1")
  expect_identical(INLA::inla.getOption("num.threads"), "4:1")
})
