# =============================================================================
# fb_select() held-out relabelling is verified by identity.
#
# Under prediction = "new_cluster" every assessment row's grouping value is
# replaced by the held-out label of its own level. A check on counts (the
# number of levels, or rows per level) accepts a relabelling that permutes
# those identities: the held-out rows of one site carrying another site's
# label, or a training site's own label. The oracle is the one-to-one level
# map, written out literally here and compared element by element.
#
# No block fits a model: the check runs while a fold's engine inputs are
# built, before any engine is called.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

# Three sites of four rows and a three-level block factor crossed with them.
fbri_data <- function() {
  data.frame(
    site = factor(rep(c("s1", "s2", "s3"), each = 4L)),
    blk = factor(rep(c("b1", "b2", "b3"), times = 4L)),
    x = seq(0.1, 1.2, by = 0.1),
    y = c(1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6)
  )
}

# One fold that splits every site: the first row of each site is held out,
# and those rows carry blocks b1, b2 and b3.
fbri_split <- function() {
  list(
    fold_id = 1L,
    repeat_id = 1L,
    label = "1",
    assess_rows = c(1L, 5L, 9L),
    train_rows = c(2:4, 6:8, 10:12)
  )
}

# `data` with column `v` rewritten on `rows` to `values`.
fbri_rewrite <- function(data, v, rows, values) {
  column <- as.character(data[[v]])
  column[rows] <- values
  data[[v]] <- factor(column)
  data
}

# The two counts a count-only check compares on the assessment rows: the
# number of distinct levels, and the sorted rows per level.
fbri_counts <- function(column, rows) {
  values <- as.character(column)[rows]
  list(
    n_levels = length(unique(values)),
    per_level = sort(as.integer(table(values)))
  )
}

fbri_refusal <- "flexybayes_refusal_selection_new_cluster_not_representable"


# --- a correct relabelling passes ---------------------------------------- #

test_that("a correct held-out relabelling carries each row's own level and passes", {
  d <- fbri_data()
  split <- fbri_split()
  inputs <- .fb_sel_bayes_fold_inputs(
    d, split, "y", ~x, NULL, random = c("site", "blk")
  )
  e <- inputs$engine_data
  a <- split$assess_rows

  expect_identical(
    as.character(e$site[a]),
    c("s1::fb_heldout", "s2::fb_heldout", "s3::fb_heldout")
  )
  expect_identical(
    as.character(e$blk[a]),
    c("b1::fb_heldout", "b2::fb_heldout", "b3::fb_heldout")
  )
  expect_identical(
    as.character(e$site[split$train_rows]),
    as.character(d$site[split$train_rows])
  )
  expect_identical(
    as.character(e$blk[split$train_rows]),
    as.character(d$blk[split$train_rows])
  )
  expect_null(.fb_sel_check_heldout_labels(d, e, split, c("site", "blk")))

  # known_cluster relabels nothing.
  known <- .fb_sel_bayes_fold_inputs(
    d, split, "y", ~x, NULL,
    random = c("site", "blk"), prediction = "known_cluster"
  )
  expect_identical(known$engine_data$site, d$site)
  expect_identical(known$engine_data$blk, d$blk)
})


# --- a count-preserving permutation of identities refuses ---------------- #

test_that("a count-preserving permutation of held-out identities refuses by name", {
  d <- fbri_data()
  split <- fbri_split()
  a <- split$assess_rows

  # s1 -> s2, s2 -> s3, s3 -> s1 on the assessment rows: once with the
  # held-out suffix, and once onto the training sites' own labels.
  cases <- list(
    held_out_labels = c("s2::fb_heldout", "s3::fb_heldout", "s1::fb_heldout"),
    training_labels = c("s2", "s3", "s1")
  )
  correct <- c("s1::fb_heldout", "s2::fb_heldout", "s3::fb_heldout")
  for (nm in names(cases)) {
    e <- fbri_rewrite(d, "site", a, cases[[nm]])
    # Counts agree with the correct relabelling and with the original
    # levels, so a check on them sees nothing.
    expect_identical(
      fbri_counts(e$site, a),
      fbri_counts(fbri_rewrite(d, "site", a, correct)$site, a),
      info = nm
    )
    expect_identical(fbri_counts(e$site, a), fbri_counts(d$site, a), info = nm)

    err <- expect_error(
      .fb_sel_check_heldout_labels(d, e, split, "site"),
      class = fbri_refusal
    )
    expect_match(conditionMessage(err), "`site`", fixed = TRUE)
  }
})

test_that("a permuted relabelling refuses where the fold's engine inputs are built", {
  d <- fbri_data()
  split <- fbri_split()
  reached <- 0L
  testthat::local_mocked_bindings(
    .fb_sel_bayes_heldout_labels = function(data, assess_rows, random_vars) {
      reached <<- reached + 1L
      fbri_rewrite(
        data, "site", assess_rows,
        c("s2::fb_heldout", "s3::fb_heldout", "s1::fb_heldout")
      )
    }
  )
  expect_error(
    .fb_sel_bayes_fold_inputs(d, split, "y", ~x, NULL, random = "site"),
    class = fbri_refusal
  )
  expect_identical(reached, 1L)
})

test_that("a permutation in any one relabelled column refuses and names that column", {
  d <- fbri_data()
  split <- fbri_split()
  a <- split$assess_rows
  e <- fbri_rewrite(
    d, "site", a, c("s1::fb_heldout", "s2::fb_heldout", "s3::fb_heldout")
  )
  e <- fbri_rewrite(
    e, "blk", a, c("b2::fb_heldout", "b3::fb_heldout", "b1::fb_heldout")
  )
  expect_identical(fbri_counts(e$blk, a), fbri_counts(d$blk, a))
  err <- expect_error(
    .fb_sel_check_heldout_labels(d, e, split, c("site", "blk")),
    class = fbri_refusal
  )
  expect_match(conditionMessage(err), "`blk`", fixed = TRUE)
  expect_no_match(conditionMessage(err), "`site`", fixed = TRUE)

  # A training row moved is caught as well.
  moved <- fbri_rewrite(
    d, "site", a, c("s1::fb_heldout", "s2::fb_heldout", "s3::fb_heldout")
  )
  moved <- fbri_rewrite(moved, "site", c(2L, 6L), c("s2", "s1"))
  expect_identical(
    fbri_counts(moved$site, split$train_rows),
    fbri_counts(d$site, split$train_rows)
  )
  expect_error(
    .fb_sel_check_heldout_labels(d, moved, split, "site"),
    class = fbri_refusal
  )
})


# --- a held-out label a training row already carries refuses ------------ #

test_that("a training level spelled like a held-out label refuses before any fit", {
  d <- fbri_data()
  d$site <- factor(c(rep("s1", 4L), rep("s2", 4L), rep("s1::fb_heldout", 4L)))
  split <- list(
    fold_id = 1L, repeat_id = 1L, label = "1",
    assess_rows = 1:2, train_rows = 3:12
  )
  err <- expect_error(
    .fb_sel_bayes_fold_inputs(d, split, "y", ~x, NULL, random = "site"),
    class = fbri_refusal
  )
  expect_match(conditionMessage(err), "s1::fb_heldout", fixed = TRUE)

  known <- .fb_sel_bayes_fold_inputs(
    d, split, "y", ~x, NULL, random = "site", prediction = "known_cluster"
  )
  expect_identical(known$engine_data$site, d$site)
})
