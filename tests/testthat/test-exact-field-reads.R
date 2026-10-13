# Reading list fields by exact name.
#
# `$` on a list falls back to a unique prefix match when no element carries
# the exact name. A multi-way interaction
# term has `vars` and no `var`, so `term$var` returned the interaction's
# columns. The aggregation plan names its cell count `K_est`, so
# `agg_plan$K` returned it. A convention whose response scale had not been
# filled returned its `response_scale_rule` string for `response_scale`.
# Each reader names the field it means. The tests run the readers with
# partial matching reported and collect every match raised in this package's
# own code. Matches raised inside an engine are not collected.

# The namespace a warning was raised in: the deepest frame below the warning
# signal whose function is not base R's.
.pm_warning_origin <- function() {
  calls <- sys.calls()
  signal_at <- which(vapply(
    calls,
    function(cl) identical(cl[[1L]], quote(.signalSimpleWarning)),
    logical(1L)
  ))
  if (!length(signal_at)) {
    return(NA_character_)
  }
  for (k in rev(seq_len(max(signal_at) - 1L))) {
    env <- environment(sys.function(k))
    if (is.null(env)) {
      next
    }
    origin <- environmentName(topenv(env))
    if (!identical(origin, "base")) {
      return(origin)
    }
  }
  NA_character_
}

# Evaluates `expr` with `$` partial matching reported. Returns the value, or
# the error the call raised, and the partial-match messages raised in this
# package's code.
.pm_collect <- function(expr) {
  matches <- character(0)
  value <- withr::with_options(
    list(
      warnPartialMatchDollar = TRUE,
      flexyBayes.silence_default_prior_note = TRUE
    ),
    withCallingHandlers(
      tryCatch(suppressMessages(expr), error = function(e) e),
      warning = function(w) {
        if (
          grepl("partial match of", conditionMessage(w), fixed = TRUE) &&
            identical(.pm_warning_origin(), "flexyBayes")
        ) {
          matches <<- c(matches, conditionMessage(w))
        }
        invokeRestart("muffleWarning")
      }
    )
  )
  list(value = value, matches = matches)
}


test_that("a multi-way random term is labelled from its vars", {
  term <- list(type = "combo", vars = c("a", "b", "cc"))
  res <- .pm_collect(.preflight_term_label(term, "combo", kind = "random"))

  expect_identical(res$matches, character(0))
  expect_identical(res$value, "(combo | a:b:cc)")
})

test_that("a multi-way fixed term with no label is named from its vars", {
  term <- list(type = "interaction", vars = c("x", "z"))
  res <- .pm_collect(.preflight_term_label(term, "interaction", kind = "fixed"))

  expect_identical(res$matches, character(0))
  expect_identical(res$value, "x:z")
})

test_that("a multi-way random term's level count is its observed combinations", {
  d <- data.frame(
    a = factor(c(1, 1, 2, 2, 3, 3)),
    b = factor(c(1, 2, 1, 2, 1, 1)),
    cc = factor(c(1, 1, 1, 2, 2, 2))
  )
  term <- list(type = "combo", vars = c("a", "b", "cc"))
  res <- .pm_collect(.preflight_term_level_count(term, .fb_dataset(d)))

  expect_identical(res$matches, character(0))
  # Six rows, the last two the same combination.
  expect_identical(res$value, nrow(unique(d)))
})

test_that("a three-way interaction random term plans with no partial match", {
  set.seed(1)
  d <- expand.grid(a = factor(1:4), b = factor(1:3), cc = factor(1:2), rep = 1:3)
  d$y <- rnorm(nrow(d))
  res <- .pm_collect(flexybayes(y ~ 1, random = ~ a:b:cc, data = d, plan = TRUE))

  expect_identical(res$matches, character(0))
  expect_s3_class(res$value, "fb_plan")
})

test_that("the plan reads the aggregation cell count by its name", {
  set.seed(2)
  d <- data.frame(
    f = factor(rep(1:3, each = 100)),
    g = factor(rep(1:10, times = 30))
  )
  d$count <- stats::rpois(nrow(d), 3)
  res <- .pm_collect(flexybayes(
    count ~ f,
    random = ~g,
    data = d,
    family = "poisson",
    backend = "auto",
    plan = TRUE
  ))
  n_cells <- nrow(unique(d[c("f", "g")]))

  expect_identical(res$matches, character(0))
  expect_equal(res$value$aggregation$K, n_cells)
  expect_equal(res$value$aggregation$ratio, nrow(d) / n_cells)
})

