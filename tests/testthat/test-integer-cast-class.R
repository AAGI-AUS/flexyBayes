# Integer casts of a product of counts.
#
# `as.integer()` past 2^31 - 1 returns NA with only a coercion warning.
# A cell count of prod(Ks) passes that limit on real designs: about 2.4e9
# cells at 54,208 genotypes of `agridat::barrero.maize`.
#
# Two tests: (a) a grep gate over R/ so the class cannot re-enter one
# member at a time; (b) the aggregation planner refuses by name, built
# from level counts (not a data frame that size).

# ---------------------------------------------------------------- #
# (a) grep gate                                                     #
# ---------------------------------------------------------------- #

test_that("grep gate: no un-allowlisted as.integer(prod/nrow/length( cast", {
  root <- testthat::test_path("..", "..")
  r_dir <- file.path(root, "R")
  skip_if_not(
    dir.exists(r_dir),
    "package source tree not available (installed tree, not source checkout)"
  )

  patterns <- c("as.integer(prod(", "as.integer(nrow(", "as.integer(length(")

  # Sites exempt from the gate, each named with a one-line reason.
  allowlist <- character(0)

  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  hits <- list()
  for (f in r_files) {
    lines <- readLines(f, warn = FALSE)
    for (i in seq_along(lines)) {
      line <- lines[[i]]
      if (grepl("^\\s*#", line)) next # full-line comment, not code
      for (p in patterns) {
        if (grepl(p, line, fixed = TRUE)) {
          hits[[length(hits) + 1L]] <- list(
            file = basename(f),
            line = i,
            text = trimws(line)
          )
        }
      }
    }
  }

  offenders <- Filter(function(h) !(h$file %in% names(allowlist)), hits)
  expect_identical(
    length(offenders),
    0L,
    label = paste0(
      "un-allowlisted cast(s): ",
      paste(
        vapply(
          offenders,
          function(h) paste0(h$file, ":", h$line, " `", h$text, "`"),
          character(1L)
        ),
        collapse = "; "
      )
    )
  )

  # An allowlist entry must still name a cast.
  for (nm in names(allowlist)) {
    expect_true(
      any(vapply(hits, function(h) identical(h$file, nm), logical(1L))),
      label = paste0(
        "allowlist entry '", nm, "' names no surviving cast: remove it"
      )
    )
  }
})


# ---------------------------------------------------------------- #
# (b) planner refuses on cell-count overflow                        #
# ---------------------------------------------------------------- #

# Minimal hand-built `<fb_terms>` IR. A copy of `.test_make_plan_ir()`
# in test-aggregate-plan.R, so this file also runs on its own.
.icc_make_plan_ir <- function(fixed_terms = list(), n = 1000L) {
  structure(
    list(
      response = "y",
      family = "gaussian",
      link = "identity",
      intercept = TRUE,
      fixed_terms = fixed_terms,
      random_terms = list(),
      residual_terms = list(),
      addition_terms = list(),
      priors = list(),
      data_summary = list(n = n),
      capabilities = character(),
      source = "test"
    ),
    class = c("fb_terms", "list")
  )
}

test_that("aggregation planner refuses by name past the integer cell-count limit", {
  # Two factors at 50,000 levels each: product = 2.5e9, past
  # 2^31 - 1 = 2,147,483,647. Built from dictionaries of labels, not a
  # 2.5e9-row frame: .fb_dataset_levels() reads
  # length(dictionaries[[var]]).
  fb <- .icc_make_plan_ir(
    fixed_terms = list(
      list(type = "factor", var = "f1"),
      list(type = "factor", var = "f2")
    ),
    n = 1000L
  )
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", f1 = "factor", f2 = "factor"),
    dictionaries = list(
      f1 = as.character(seq_len(50000L)),
      f2 = as.character(seq_len(50000L))
    )
  )

  # The planner does not refuse: it records the reason, so the default
  # aggregate = "auto" takes the per-row route as 0.9.2 did (found by the
  # 27-Aug barrero side-by-side script, which 0.9.3 had refused outright).
  plan <- .fb_aggregation_plan(fb, ds)
  expect_false(plan$eligible)
  expect_identical(plan$reason_codes, "cell_count_exceeds_integer")
  expect_true(is.na(plan$K_est))
  expect_true(plan$K_est_double > .Machine$integer.max)
  err <- plan$overflow_refusal
  expect_s3_class(err, "flexybayes_refusal_cell_count_exceeds_integer")
  expect_identical(err$reason_code, "cell_count_exceeds_integer")
  expect_match(err$message, "2,500,000,000|2500000000", perl = TRUE)
})

test_that("aggregate = TRUE on an overflowing cell key refuses by name before any fit, and auto does not", {
  skip_if_not_installed("INLA")
  # The overflow plan reaches the dispatch gate through a mocked
  # planner. A real 1,000-row frame can stop first on
  # `compression_unproductive`.
  fb_over <- .icc_make_plan_ir(
    fixed_terms = list(
      list(type = "factor", var = "f1"),
      list(type = "factor", var = "f2")
    ),
    n = 1000L
  )
  ds_over <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", f1 = "factor", f2 = "factor"),
    dictionaries = list(
      f1 = as.character(seq_len(50000L)),
      f2 = as.character(seq_len(50000L))
    )
  )
  plan_over <- .fb_aggregation_plan(fb_over, ds_over)
  expect_identical(plan_over$reason_codes, "cell_count_exceeds_integer")
  testthat::local_mocked_bindings(
    .fb_aggregation_plan = function(...) plan_over
  )
  d <- data.frame(
    y = rnorm(24L),
    f1 = factor(rep(c("a", "b"), 12L)),
    f2 = factor(rep(c("p", "q", "r"), 8L))
  )
  expect_error(
    suppressMessages(flexybayes(y ~ f1, random = ~ f2, data = d,
      backend = "inla", aggregate = TRUE)),
    class = "flexybayes_refusal_cell_count_exceeds_integer"
  )
  # auto: the plan surface carries the reason and the route is per row.
  p <- fb_plan(y ~ f1, random = ~ f2, data = d, backend = "inla",
    aggregate = "auto")
  expect_true("cell_count_exceeds_integer" %in% p$aggregation$reason_codes)
})

test_that("aggregation planner does not refuse below the integer cell-count limit", {
  # Same shape, under the limit.
  fb <- .icc_make_plan_ir(
    fixed_terms = list(
      list(type = "factor", var = "f1"),
      list(type = "factor", var = "f2")
    ),
    n = 1000L
  )
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(y = "double", f1 = "factor", f2 = "factor"),
    dictionaries = list(
      f1 = as.character(seq_len(5L)),
      f2 = as.character(seq_len(20L))
    )
  )
  plan <- .fb_aggregation_plan(fb, ds)
  expect_identical(plan$K_est, 100L)
})


# ---------------------------------------------------------------- #
# (c) .preflight_factor_interaction_dummies()                        #
# ---------------------------------------------------------------- #

test_that(paste(
  ".preflight_factor_interaction_dummies() carries the product as a",
  "double and returns NA_integer_ cleanly past the limit"
), {
  # Three factors at 3,000 levels each: (3000 - 1)^3 ~= 2.7e10, past
  # 2^31 - 1. Built from dictionaries, as in (b).
  ds <- .fb_dataset(
    data = NULL,
    n_rows = 1000L,
    col_types = list(
      y = "double", f1 = "factor", f2 = "factor", f3 = "factor"
    ),
    dictionaries = list(
      f1 = as.character(seq_len(3000L)),
      f2 = as.character(seq_len(3000L)),
      f3 = as.character(seq_len(3000L))
    )
  )
  vars <- c("f1", "f2", "f3")

  # The package function returns NA_integer_ with no coercion warning.
  expect_no_warning(
    result <- .preflight_factor_interaction_dummies(vars, ds)
  )
  expect_identical(result, NA_integer_)

  # Induced defect: the same body without the `is.finite()` /
  # `> .Machine$integer.max` guard. It also returns NA, but with a
  # coercion warning.
  .icc_unfixed_dummies <- function(vars, fb_dataset) {
    if (!length(vars)) {
      return(NA_integer_)
    }
    ks <- vapply(
      vars,
      function(v) .fb_dataset_levels(fb_dataset, as.character(v)),
      numeric(1L)
    )
    if (anyNA(ks)) {
      return(NA_integer_)
    }
    ks_minus_one <- pmax(ks - 1, 0)
    as.integer(prod(ks_minus_one)) # unguarded cast
  }
  expect_warning(
    bad <- .icc_unfixed_dummies(vars, ds),
    "NAs introduced by coercion to integer range"
  )
  expect_identical(bad, NA_integer_)
})
