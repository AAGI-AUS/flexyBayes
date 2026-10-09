# A child started with callr's `supervise = TRUE` makes processx start one
# supervisor process for the session, which holds two fifo connections open.
# R CMD check stops an example that leaves connections open, so a supervisor
# that flexyBayes started is ended when its last child is collected, and one
# the session already had is left running.

.fic_supervisor_fifos <- function() {
  sum(grepl("supervisor_std", showConnections(all = TRUE)[, "description"],
            fixed = TRUE))
}

.fic_skip <- function() {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("processx")
  if (.Platform$OS.type != "unix") {
    skip("the processx supervisor talks over fifos on unix only")
  }
}

test_that("an isolated call leaves no connection open", {
  .fic_skip()
  processx::supervisor_kill()
  n_before <- nrow(showConnections(all = TRUE))
  value <- .fb_isolated_call(function() 1 + 1, engine = "inla",
                             entry_point = "flexybayes()")
  expect_identical(value, 2)
  expect_identical(nrow(showConnections(all = TRUE)), n_before)
  expect_identical(.fic_supervisor_fifos(), 0L)
})

test_that("jobs on two workers leave no connection open", {
  .fic_skip()
  processx::supervisor_kill()
  n_before <- nrow(showConnections(all = TRUE))
  out <- fb_isolated_lapply(1:3, function(i) i * 2, plan = list(workers = 2L))
  expect_identical(out, list(2, 4, 6))
  expect_identical(nrow(showConnections(all = TRUE)), n_before)
})

test_that("a child that fails leaves no connection open", {
  .fic_skip()
  processx::supervisor_kill()
  n_before <- nrow(showConnections(all = TRUE))
  expect_error(
    .fb_isolated_call(function() stop("planted"), engine = "inla",
                      entry_point = "flexybayes()"),
    "planted"
  )
  expect_identical(nrow(showConnections(all = TRUE)), n_before)
})

test_that("a supervisor the session already had is left running", {
  .fic_skip()
  processx::supervisor_kill()
  own <- processx::process$new("sleep", "30", supervise = TRUE)
  on.exit({
    own$kill()
    processx::supervisor_kill()
  }, add = TRUE)
  expect_identical(.fic_supervisor_fifos(), 2L)
  .fb_isolated_call(function() 1, engine = "inla",
                    entry_point = "flexybayes()")
  expect_identical(.fic_supervisor_fifos(), 2L)
})
