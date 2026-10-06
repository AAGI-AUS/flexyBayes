# fb_isolated_lapply(): jobs in separate R processes, for packages that
# extend flexyBayes.

test_that("the plan is checked field by field before any process starts", {
  expect_error(
    fb_isolated_lapply(1:2, identity, list(workers = 0L)),
    "`plan$workers` must be a whole number",
    fixed = TRUE
  )
  expect_error(
    fb_isolated_lapply(1:2, identity, list(threads_per_worker = 1.5)),
    "`plan$threads_per_worker` must be NA or a whole number",
    fixed = TRUE
  )
  expect_error(
    fb_isolated_lapply(1:2, identity, list(warm_first = NA)),
    "`plan$warm_first` must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    fb_isolated_lapply(1:2, identity, list(engine = "")),
    "`plan$engine` must be a single string",
    fixed = TRUE
  )
  expect_error(
    fb_isolated_lapply(1:2, identity, list(worker = 2L)),
    "`plan` has unknown field `worker`",
    fixed = TRUE
  )
  expect_error(fb_isolated_lapply(1:2, "identity"), "`fn` must be a function")

  plan <- .fb_isolation_plan(list())
  expect_identical(plan$workers, 1L)
  expect_identical(plan$threads_per_worker, NA_integer_)
  expect_false(plan$warm_first)
  expect_identical(plan$engine, "none")
  expect_identical(plan$entry_point, "fb_isolated_lapply()")
})

test_that("without callr the call refuses by name and names the entry point", {
  testthat::local_mocked_bindings(.fb_isolation_callr_available = function() {
    FALSE
  })
  refused <- tryCatch(
    fb_isolated_lapply(1:2, identity, list(entry_point = "mypkg_verb()")),
    error = function(e) e
  )
  expect_s3_class(refused, "flexybayes_refusal_isolation_requires_callr")
  expect_identical(refused$entry_point, "mypkg_verb()")
  expect_match(conditionMessage(refused), "^mypkg_verb\\(\\): ")
})

test_that("an empty job list returns an empty list without starting a process", {
  expect_identical(fb_isolated_lapply(list(), identity), list())
})

test_that("values come back in job order from processes other than this one", {
  skip_on_cran()
  skip_if_not_installed("callr")
  me <- Sys.getpid()
  out <- fb_isolated_lapply(1:3, function(j) {
    Sys.sleep(c(0.6, 0.3, 0)[[j]])
    list(j = j, pid = Sys.getpid(), omp = Sys.getenv("OMP_NUM_THREADS"))
  }, list(workers = 3L, threads_per_worker = 2L))
  expect_identical(vapply(out, function(z) z$j, integer(1)), 1:3)
  expect_false(any(vapply(out, function(z) z$pid, integer(1)) == me))
  expect_identical(vapply(out, function(z) z$omp, ""), rep("2", 3L))
})

test_that("the lowest failing job's error is raised, keeping its class", {
  skip_on_cran()
  skip_if_not_installed("callr")
  refused <- tryCatch(
    fb_isolated_lapply(1:3, function(j) {
      if (j >= 2L) {
        fb_refuse("isolated_process_died", paste0("job ", j, " refusing"))
      }
      j
    }, list(workers = 1L)),
    error = function(e) e
  )
  expect_s3_class(refused, "flexybayes_refusal_isolated_process_died")
  expect_match(conditionMessage(refused), "job 2 refusing")
})

test_that("a function from another namespace runs in a child that loads it", {
  skip_on_cran()
  skip_if_not_installed("callr")
  out <- fb_isolated_lapply(list(c(1, 5, 2), c(4, 4, 9)), stats::median)
  expect_identical(out, list(2, 4))
})

test_that("one worker carries the random-number state from job to job", {
  skip_on_cran()
  skip_if_not_installed("callr")
  set.seed(11L)
  isolated <- fb_isolated_lapply(1:2, function(j) stats::runif(1L))
  set.seed(11L)
  serial <- lapply(1:2, function(j) stats::runif(1L))
  expect_identical(isolated, serial)
})
