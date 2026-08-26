# test-integer-cast-class.R -- the integer-cast class (C1, FS-24).
#
# R/aggregate_plan.R:242 `K_est <- as.integer(prod(Ks))` overflowed to NA
# past 2^31 - 1 cells (a design-growth ladder on `agridat::barrero.maize`
# reached ~2.4e9 at 54,208 genotypes). `as.integer()` on a value past the
# limit returns NA with only a base-R coercion warning, and
# `compression_est` then silently divided by that NA -- the aggregation
# plan lost both its cell count and its compression estimate rather than
# refusing. The identical defect in the ROW-count casts was already found
# and fixed (`row_count_exceeds_integer`); the fix was applied to one
# member of the class rather than to the class.
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

  # Sites verified safe and exempt from the double-cast fix, each with a
  # one-line reason kept HERE (not only in the source comment) so the
  # gate and its justification cannot drift apart. Every other site the
  # audit found (R/aggregate_plan.R:242, R/fb_preflight.R:781, and a
  # third in the native-model ingest adapter withdrawn entirely at
  # 0.9.3 -- see NEWS.md) was fixed to carry the product as a double
  # before casting, and no longer matches these patterns (the withdrawn
  # file is simply gone, so the grep gate below no longer scans it at
  # all).
  #
  # The one exemption this allowlist used to carry (fb_log_posterior.R's
  # n_elem sizing of the coord_names / lower / upper vectors) was
  # removed at 0.9.3 along with the log-density producer for the
  # withdrawn native engine that owned it -- the surviving
  # fb_log_posterior() abstains unconditionally and carries no cast in
  # this class. No exemption is currently needed; the allowlist is
  # empty rather than stale.
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

  # An allowlist entry that names no surviving cast is stale
  # documentation for a fix already made -- catch drift in the other
  # direction too.
  for (nm in names(allowlist)) {
    expect_true(
      any(vapply(hits, function(h) identical(h$file, nm), logical(1L))),
      label = paste0(
        "allowlist entry '", nm, "' names no surviving cast -- remove it"
      )
    )
  }
})


# ---------------------------------------------------------------- #
# (b) planner refuses on cell-count overflow                        #
# ---------------------------------------------------------------- #

# Minimal hand-built `<fb_terms>` IR -- mirrors
# test-aggregate-plan.R's `.test_make_plan_ir()`, kept local (test files
# run individually via testthat::test_file() per the WP-C spec, so a
# helper defined in a sibling test file is not guaranteed loaded).
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
  # 2^31 - 1 = 2,147,483,647. Built entirely from metadata-only
  # dictionaries (character vectors of labels), never a 2.5e9-row data
  # frame -- .fb_dataset_levels() reads length(dictionaries[[var]]),
  # not the underlying data.
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

  expect_error(
    .fb_aggregation_plan(fb, ds),
    class = "flexybayes_refusal_cell_count_exceeds_integer"
  )

  err <- tryCatch(
    .fb_aggregation_plan(fb, ds),
    flexybayes_refusal_cell_count_exceeds_integer = function(e) e
  )
  expect_identical(err$reason_code, "cell_count_exceeds_integer")
  expect_true(err$K_est > .Machine$integer.max)
  expect_match(err$message, "2,500,000,000|2500000000", perl = TRUE)
})

test_that("aggregation planner does not refuse below the integer cell-count limit", {
  # Same shape, comfortably under the limit -- confirms the guard is a
  # ceiling check, not a blanket refusal on every two-factor model.
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
# (c) the other class member fixed without a dedicated behavioural  #
# test (F6): .preflight_factor_interaction_dummies(). (A second      #
# member, the native-model ingest adapter's model_dim reducer, was   #
# withdrawn entirely with that adapter at 0.9.3 -- see NEWS.md -- so  #
# its regression test is deleted rather than kept unreachable.)     #
# ---------------------------------------------------------------- #

test_that(paste(
  ".preflight_factor_interaction_dummies() carries the product as a",
  "double and returns NA_integer_ cleanly past the limit"
), {
  # Three factors at 3,000 levels each: (3000 - 1)^3 ~= 2.7e10, past
  # 2^31 - 1. Metadata-only dictionaries, never a 2.7e10-row frame --
  # same construction as (b) above.
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

  # Passes on the fix: the guarded, real implementation returns
  # NA_integer_ with no coercion warning.
  expect_no_warning(
    result <- .preflight_factor_interaction_dummies(vars, ds)
  )
  expect_identical(result, NA_integer_)

  # Induced defect: the pre-fix cast, `as.integer(prod(ks_minus_one))`
  # with no `is.finite()` / `> .Machine$integer.max` guard, restored
  # here as a scratch copy of the function body (R/fb_preflight.R:781
  # before this fix). Same overflow class C1 repairs elsewhere
  # (R/aggregate_plan.R:242) -- the un-guarded cast is not silent in
  # the sense of returning a wrong value (it still lands on NA), but it
  # raises an uncaught base-R coercion warning on every ordinary,
  # in-range-cardinality preflight call that happens to cross the
  # limit, which the WP-C release gate (0-warning suite footer) treats
  # as a real defect.
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
    as.integer(prod(ks_minus_one)) # the pre-fix cast, restored here only
  }
  expect_warning(
    bad <- .icc_unfixed_dummies(vars, ds),
    "NAs introduced by coercion to integer range"
  )
  expect_identical(bad, NA_integer_)
})
