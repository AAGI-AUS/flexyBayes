# ranef() names the group the user wrote, not the index INLA was given.
#
# emit_inla() builds an index column per structured random term --
# `rep_block_id` for `rep:block`, `gen_by_env_id` for an unstructured
# covariance -- and INLA names each summary.random entry after the column it
# was indexed by. The variance components were mapped back to user labels
# when that machinery was built; the random-effect table was not, so
# ranef(fit) came back with a grouping factor called `rep_block_id`. Same
# defect class as the label leak closed for summary()$varcomp, reached
# through a different accessor.

test_that("a nested term keeps its own name in ranef()", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  data(john.alpha, package = "agridat")
  fit <- suppressMessages(flexybayes(
    yield ~ rep, random = ~ rep:block + gen,
    data = john.alpha, backend = "inla", verbose = FALSE
  ))
  r <- ranef(fit)
  expect_true("rep:block" %in% names(r))
  expect_true("gen" %in% names(r))
  expect_false(any(grepl("_id$", names(r))))
  # And inside the table, not only on the list name.
  expect_identical(unique(r[["rep:block"]]$group), "rep:block")
})

test_that("a plain grouping factor is left alone", {
  # The map must not rename a group whose index IS the variable the user
  # named, which is the common case and the one a blanket rewrite would
  # break.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  data(john.alpha, package = "agridat")
  fit <- suppressMessages(flexybayes(
    yield ~ rep, random = ~gen, data = john.alpha,
    backend = "inla", verbose = FALSE
  ))
  expect_identical(names(ranef(fit)), "gen")
})

test_that("the label map covers every structured term type it claims", {
  # Unit-level, so it runs without INLA: build the IR for each structured
  # term and check the generated index name maps to something a reader
  # would recognise.
  d <- expand.grid(rep = factor(1:2), gen = factor(1:6), env = factor(1:3))
  d$y <- stats::rnorm(nrow(d))
  d$x <- stats::rnorm(nrow(d))

  cases <- list(
    list(random = ~ rep:gen, index = "rep_gen_id", label = "rep:gen"),
    list(random = ~ diag(env):gen, index = "gen_at_env_1",
         label = "gen (env = 1)"),
    list(random = ~ us(env):gen, index = "gen_by_env_id", label = "gen:env")
  )
  for (cs in cases) {
    fb <- fb_from_asreml(y ~ 1, random = cs$random, data = d)
    map <- .fb_inla_group_labels(fb)
    expect_true(cs$index %in% names(map), info = cs$index)
    expect_identical(unname(map[[cs$index]]), cs$label)
  }
})

test_that("the correlated slope's two index blocks take distinct labels", {
  # They used to share one. The iid2d field is laid out intercepts then
  # slopes, so the two index names address two different blocks of
  # deviations, and mapping both to `x | g` returned a list with two
  # elements of the same name -- `ranef(fit)[["x | g"]]` reached only the
  # first, and no caller could tell which block that was.
  d <- data.frame(y = stats::rnorm(40), x = stats::rnorm(40),
                  g = factor(rep(letters[1:8], each = 5)))
  fb <- fb_from_brms(y ~ x + (x | g), data = d)
  map <- .fb_inla_group_labels(fb)
  expect_identical(unname(map[["g_cor_x_id"]]), "(Intercept) | g")
  expect_identical(unname(map[["g_cor_x_slope_id"]]), "x | g")
})


# --- the levels inside the table, not only the name on it -------------- #
#
# Finding S01's defect D2. `.fb_random_inla()` restored a level label only
# where INLA's `ID` column came back CHARACTER, which is what a plain `~ g`
# produces. Every structured term is indexed by an integer column
# emit_inla() built, so `ID` is 1..n and the guard never fired: `rep:block`
# on john.alpha reported levels "1" to "18" against ASReml's
# `rep_R1:block_B1 ...`, and `(Days | Subject)` reported "1" to "36" with
# the subject identifiers sitting unused on the fit. Two of the four BLUP
# comparisons in the S01 programme could not be made at all, and aligning by
# row order was refused because `interaction()` orders `rep` fastest while
# ASReml orders it slowest -- a positional match would have been silently
# wrong on exactly the row where it looked easiest.

test_that("an interaction group reports the interaction's own levels", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  e <- new.env()
  utils::data("john.alpha", package = "agridat", envir = e)
  d <- e$john.alpha

  fit <- suppressMessages(flexybayes(
    yield ~ rep, random = ~ rep:block + gen,
    data = d, backend = "inla", verbose = FALSE
  ))
  got <- ranef(fit)[["rep:block"]]$level

  # The expected labels are written out here rather than taken from
  # `interaction()`, which is the call the emitter itself makes: an oracle
  # that shares the code under test proves self-consistency. The index is
  # lexicographic with the FIRST variable slowest, so the order is
  # `order(rep, block)` over the observed cells.
  cells <- unique(d[c("rep", "block")])
  cells <- cells[order(
    as.integer(factor(cells$rep, levels = levels(factor(d$rep)))),
    as.integer(factor(cells$block, levels = levels(factor(d$block))))
  ), ]
  expected <- paste0(cells$rep, ":", cells$block)

  expect_identical(got, expected)
  expect_length(got, 18L)
  # And nothing that looks like a bare index survived.
  expect_false(any(grepl("^[0-9]+$", got)))
})

test_that("a correlated slope reports two named blocks over the group's levels", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  e <- new.env()
  utils::data("sleepstudy", package = "lme4", envir = e)
  d <- e$sleepstudy

  fit <- suppressWarnings(suppressMessages(flexybayes(
    Reaction ~ Days + (Days | Subject),
    data = d, backend = "inla", verbose = FALSE
  )))
  r <- ranef(fit)

  # Two elements, and a caller can tell which is which. Both used to be
  # called "Days | Subject", so `r[["Days | Subject"]]` reached one of the
  # two blocks and there was no way to say which.
  expect_identical(anyDuplicated(names(r)), 0L)
  expect_setequal(names(r), c("(Intercept) | Subject", "Days | Subject"))

  subjects <- levels(factor(d$Subject))
  expect_length(subjects, 18L)
  for (nm in names(r)) {
    expect_identical(r[[nm]]$level, subjects, info = nm)
    expect_identical(unique(r[[nm]]$group), nm)
  }

  # The two blocks are different numbers, not one block reported twice --
  # the copy INLA emits for the slope carries the whole field, so reading
  # the same rows out of both would return the intercepts under both names.
  expect_false(isTRUE(all.equal(
    r[["(Intercept) | Subject"]]$estimate,
    r[["Days | Subject"]]$estimate
  )))
})
