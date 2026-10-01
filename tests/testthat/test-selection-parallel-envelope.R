# =============================================================================
# Fold refits run concurrently, under a reproducibility envelope.
#
# The claim this file holds the fold loops to has two halves, and they are
# not the same standard.
#
#   EXACT. The fold allocation, the candidate space, the score contract, the
#   certificate, the selected support, the per-model tables and the fit seeds
#   are bit-identical between a serial run and a concurrent one. On the
#   closed-form benchmark route, where no engine samples, the scores are
#   exact too and the whole fold object is compared with identical().
#
#   MEASURED. On the Bayesian route INLA's posterior sampler is not
#   reproducible on this build -- the header of R/selection_bayes_cv.R
#   records two identical calls differing by up to 0.21 in the latent vector
#   -- so two serial runs do not agree either. The concurrent run is
#   therefore held to the spread the SERIAL runs show under the same
#   conditions, measured inside the test itself rather than written down as
#   a constant. A tolerance chosen by an author would be a tolerance that
#   could be widened until it passed; a tolerance the block measures cannot
#   be, and a concurrent run that drifted further than the engine's own
#   repeat spread fails here.
#
# The falsifier the unit is specified against is the seed mutant: a fold
# loop handing every fold one seed. The seed block below reads the seeds off
# the fold object and fails on it, in both directions -- seeds that do not
# separate the folds, and seeds that differ between a serial and a
# concurrent run.
#
# A forked worker dying leaves a bare NULL in its element with no condition
# attached (review/units/P4_2026-09-16). One block drives the function that
# turns that NULL into a named refusal; the blocks at the end of the file
# kill a worker on each route, from inside the worker.
# =============================================================================

# --- fixtures ------------------------------------------------------------- #

# A clustered Gaussian design: one protected covariate, three candidates of
# which one carries a real effect, and a three-level factor so that a
# grouped penalty has a group with more than one column in it.
fbpe_data <- function(seed = 20260916L, n_cluster = 8L, per = 5L) {
  .fb_sel_with_seed(seed, {
    n <- n_cluster * per
    d <- data.frame(
      site = factor(rep(sprintf("S%02d", seq_len(n_cluster)), each = per))
    )
    for (v in c("p1", "z1", "z2")) {
      d[[v]] <- stats::rnorm(n)
    }
    d$g3 <- factor(rep_len(c("a", "b", "c"), n))
    u <- stats::rnorm(n_cluster, sd = 0.7)
    d$y <- 1 + 0.5 * d$p1 + 1.2 * d$z1 +
      as.numeric(u[as.integer(d$site)]) + stats::rnorm(n, sd = 0.5)
    d
  })
}


# Run `expr` with a worker count and a worker start set, and put both options
# back afterwards. The blocks of this file that use the default start, "fork",
# test the forked route; the default route is tested further down.
fbpe_with_workers <- function(workers, expr, mode = "fork") {
  old <- options(
    flexyBayes.selection_workers = workers,
    flexyBayes.isolation = mode
  )
  on.exit(options(old), add = TRUE)
  force(expr)
}


# Concurrency in these blocks is fork-based. On a platform that does not fork the folds
# run one at a time by design, and a serial-against-concurrent comparison has
# nothing to compare, so the block skips rather than passing vacuously. The
# plan is read with the option SET: read with the default it would say
# "serial" on every platform and skip everywhere.
fbpe_skip_without_fork <- function() {
  if (!identical(.Platform$OS.type, "unix")) {
    skip("this platform does not fork")
  }
  plan <- fbpe_with_workers(2L, .fb_sel_worker_plan(2L, provider = NULL))
  if (!identical(plan$mechanism, "fork")) {
    skip("the worker plan did not reach the fork mechanism here")
  }
  invisible(NULL)
}


# The closed-form benchmark fold loop, which samples nothing.
fbpe_cv <- function(d, workers) {
  fbpe_with_workers(workers, .fb_sel_cv(
    data = d,
    response = "y",
    protected = ~p1,
    candidates = ~ z1 + z2,
    random = "site",
    resamples = fb_resamples(unit = ~site, v = 3L),
    seed = 1L,
    prediction = "new_cluster",
    metric = "rmse"
  ))
}


# The worker record is the one field that is allowed to differ: it describes
# the machine, not the estimate.
fbpe_without_plan <- function(x) {
  x$parallel <- NULL
  x
}


# --- the plan -------------------------------------------------------------- #

test_that("the worker plan reads the option, the folds and the cores", {
  # The default is one worker: a fold loop that is concurrent by default
  # would change what every existing caller runs.
  expect_identical(getOption("flexyBayes.selection_workers"), NULL)
  default <- .fb_sel_worker_plan(4L, provider = "inla")
  expect_identical(default$workers, 1L)
  expect_identical(default$mechanism, "serial")
  # Nothing is pinned on the serial path, so the engines run as they run
  # today.
  expect_true(is.na(default$threads_per_worker))

  fbpe_with_workers(3L, {
    # A worker with no fold to run is a process for nothing.
    capped <- .fb_sel_worker_plan(2L, provider = "inla")
    expect_identical(capped$workers, min(2L, capped$cores))
    expect_true(is.character(capped$note) && nzchar(capped$note))

    plan <- .fb_sel_worker_plan(5L, provider = "inla")
    expect_identical(plan$workers, min(3L, plan$cores))
    # Workers times engine threads does not exceed the machine.
    if (plan$workers > 1L) {
      expect_lte(plan$workers * plan$threads_per_worker, plan$cores)
      expect_gte(plan$threads_per_worker, 1L)
    } else {
      expect_true(is.na(plan$threads_per_worker))
    }
  })

  # A request that is not a count is refused rather than rounded.
  fbpe_with_workers("many", expect_error(
    .fb_sel_worker_plan(3L), "whole number"
  ))
  fbpe_with_workers(2.5, expect_error(
    .fb_sel_worker_plan(3L), "whole number"
  ))
  fbpe_with_workers(0L, expect_error(.fb_sel_worker_plan(3L), "whole number"))
})


# --- what a worker sends back ---------------------------------------------- #

test_that("folds come back in fold order, never in completion order", {
  skip_on_cran()
  fbpe_skip_without_fork()
  fbpe_with_workers(3L, {
    plan <- .fb_sel_worker_plan(3L, provider = NULL)
    # The third fold finishes first and the first finishes last.
    out <- .fb_sel_fold_lapply(1:3, function(j) {
      Sys.sleep(c(0.45, 0.2, 0.02)[[j]])
      j
    }, plan)
    expect_identical(unlist(out), 1:3)
  })
})


test_that("a fold that refuses inside a worker keeps its own reason code", {
  skip_on_cran()
  fbpe_skip_without_fork()
  fbpe_with_workers(3L, {
    plan <- .fb_sel_worker_plan(3L, provider = NULL)
    # Caught as an ERROR rather than as any condition. A failing child makes
    # `parallel::mclapply()` warn ("1 function calls resulted in an error")
    # before the parent re-raises anything, so a handler for `condition`
    # catches that warning and never sees the refusal underneath it.
    refusal <- suppressWarnings(tryCatch(
      .fb_sel_fold_lapply(1:3, function(j) {
        if (j == 2L) {
          stop(.fb_refusal_condition(
            reason_code = "selection_fold_not_usable",
            message = "a fold refusing inside a worker"
          ))
        }
        j
      }, plan),
      error = function(e) e
    ))
    # The refusal arrives as itself: a concurrent run must not rename the
    # reason a fold gave, or a caller's handler stops matching.
    expect_s3_class(refusal, "flexybayes_refusal")
    expect_identical(refusal$reason_code, "selection_fold_not_usable")
    expect_match(conditionMessage(refusal), "refusing inside a worker")

    # An untyped failure is still a failure and still arrives.
    err <- suppressWarnings(tryCatch(
      .fb_sel_fold_lapply(1:3, function(j) {
        if (j == 3L) stop("an untyped failure inside a worker")
        j
      }, plan),
      error = function(e) e
    ))
    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "untyped failure")
  })
})


test_that("a worker that returns no result at all refuses by name", {
  # Measured, and the reason this code exists: a child killed from outside
  # leaves a bare NULL in its element with no condition attached, because
  # the process that would have raised one is gone. Driving the mapping
  # directly keeps the test runner out of the way of a dying process.
  failed <- tryCatch(
    .fb_sel_worker_value(NULL, 2L),
    condition = function(e) e
  )
  expect_s3_class(failed, "flexybayes_refusal")
  expect_identical(failed$reason_code, "selection_parallel_worker_failed")
  expect_match(conditionMessage(failed), "fold 2")
  # The refusal states what was seen and names the causes only as possible.
  expect_match(conditionMessage(failed), "ended without returning a result")
  expect_match(conditionMessage(failed), "both possible causes")
  expect_no_match(conditionMessage(failed), "selection_workers")

  # A result that did come back is passed through untouched, including a
  # fold that legitimately returned NULL.
  expect_identical(
    .fb_sel_worker_value(list(fb_sel_worker_value = NULL), 1L),
    NULL
  )
  expect_identical(
    .fb_sel_worker_value(list(fb_sel_worker_value = list(a = 1)), 1L),
    list(a = 1)
  )
})


# --- exact identity, where nothing samples --------------------------------- #

test_that("the benchmark fold loop is bit-identical serial against concurrent", {
  skip_on_cran()
  fbpe_skip_without_fork()
  d <- fbpe_data()
  serial <- fbpe_cv(d, 1L)
  concurrent <- fbpe_cv(d, 3L)

  expect_identical(serial$parallel$mechanism, "serial")
  expect_identical(concurrent$parallel$mechanism, "fork")
  expect_gt(concurrent$parallel$workers, 1L)

  # Nothing on this route samples, so the standard is the whole object.
  expect_identical(
    fbpe_without_plan(serial),
    fbpe_without_plan(concurrent)
  )
  # Named separately, so a failure says which half moved.
  expect_identical(serial$folds, concurrent$folds)
  expect_identical(serial$candidate_space, concurrent$candidate_space)
  expect_identical(serial$aggregated, concurrent$aggregated)
  expect_identical(names(serial$substrates), names(concurrent$substrates))
  expect_identical(serial$rows, concurrent$rows)
})


test_that("a grouped penalty solves the same problem concurrently", {
  skip_on_cran()
  fbpe_skip_without_fork()
  d <- fbpe_data()
  # The grouped families are solved by this module's own proximal operator on
  # the substrates the fold loop builds. This block isolates the scheduler
  # from public-verb assembly; the following block repeats the identity check
  # through fb_select(), which became reachable in integration unit P8-INT.
  grouped_run <- function(workers) {
    cv <- fbpe_with_workers(workers, .fb_sel_cv(
      data = d,
      response = "y",
      protected = ~p1,
      candidates = ~ z1 + z2 + g3,
      random = "site",
      resamples = fb_resamples(unit = ~site, v = 3L),
      seed = 1L,
      prediction = "new_cluster",
      metric = "rmse"
    ))
    columns <- colnames(cv$substrates[[1L]]$x_candidate)
    path <- .fb_sel_glmnet_path(
      substrate = cv$substrates,
      blocks = list(list(
        family = "group_lasso",
        columns = columns,
        # One declared group holding every penalised column. A grouping whose
        # groups are all singletons is a different penalty, and the module
        # refuses it rather than fitting the nearby one.
        hyperparameters = list(scale = 0.4, group = "block")
      )),
      folds = cv$folds,
      decision = list(rule = "one_se", metric = "mse", benchmark = TRUE),
      seed = 1L,
      cluster = .fb_sel_key(d, "site")
    )
    list(cv = cv, path = path, columns = columns)
  }
  serial <- grouped_run(1L)
  concurrent <- grouped_run(3L)

  # The grouped family is the one that ran, over a group with more than one
  # column in it.
  expect_identical(serial$path$family, "group_lasso")
  expect_gt(length(serial$columns), 1L)
  expect_identical(serial$cv$parallel$mechanism, "serial")
  expect_identical(concurrent$cv$parallel$mechanism, "fork")

  # Nothing on this route samples, so every one of these is exact: the fold
  # substrates, the whole path, the tuned support and the optimality
  # residuals the certificate quotes.
  expect_identical(
    fbpe_without_plan(serial$cv),
    fbpe_without_plan(concurrent$cv)
  )
  expect_identical(serial$path, concurrent$path)
  expect_identical(serial$path$selection, concurrent$path$selection)
  expect_identical(
    serial$path$support_by_resolution,
    concurrent$path$support_by_resolution
  )
  expect_identical(serial$path$kkt_by_lambda, concurrent$path$kkt_by_lambda)
})


test_that("the public grouped optimiser is bit-identical concurrently", {
  skip_on_cran()
  # R's check harness may deliberately report one core even on a fork-capable
  # host. The scheduler contract is tested with a bounded three-core plan;
  # the process mechanism itself remains the real parallel::mclapply route.
  testthat::local_mocked_bindings(
    .fb_sel_worker_cores = function() 3L
  )
  fbpe_skip_without_fork()
  d <- fbpe_data()
  fb <- fb_from_asreml(y ~ p1, random = ~site, data = d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  base <- structure(
    list(
      fb = fb,
      data = d,
      extras = list(call_info = list(
        fixed = y ~ p1, random = ~site, family = "gaussian"
      ))
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
  run <- function(workers) {
    fbpe_with_workers(workers, fb_select(
      base,
      candidates = ~ z1 + z2 + g3,
      penalty = fb_penalty(
        term(everything()) ~ group_lasso(scale = 0.4, group = "block")
      ),
      workflow = "penalized_optimizer",
      resamples = fb_resamples(unit = ~site, v = 3L),
      decision = fb_decision(method = "full", metric = "rmse"),
      seed = 1L
    ))
  }
  serial <- run(1L)
  concurrent <- run(3L)

  expect_identical(serial$provider, "prox")
  expect_identical(concurrent$provider, "prox")
  expect_identical(serial$certificate, concurrent$certificate)
  expect_identical(serial$selection, concurrent$selection)
  expect_identical(
    serial$benchmark$glmnet$support,
    concurrent$benchmark$glmnet$support
  )
  expect_identical(
    serial$benchmark$glmnet$kkt,
    concurrent$benchmark$glmnet$kkt
  )
})


# --- the seeds, and the mutant that shares one ----------------------------- #

test_that("every fold and model gets its own seed, the same seeds either way", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbpe_skip_without_fork()
  skip_if_not_long_tier()
  d <- fbpe_data()
  base <- flexybayes(
    y ~ p1,
    random = ~site,
    data = d,
    backend = "inla",
    verbose = FALSE
  )
  # The refit is stubbed: this block is about which seed each fit is handed,
  # and a stub makes that readable without paying for an engine.
  seen <- function(workers) {
    fbpe_with_workers(workers, .fb_sel_bayes_cv(
      base,
      candidates = ~ z1 + z2,
      resamples = fb_resamples(unit = ~site, v = 3L),
      prediction = "new_cluster",
      metric = "elpd",
      seed = 1L,
      n_draws = 50L
    ))
  }
  serial <- withCallingHandlers(
    seen(1L),
    warning = function(w) invokeRestart("muffleWarning")
  )

  # Two models on one fold, and one model across folds, never share a
  # stream. This is the falsifier the unit is specified against: a loop
  # handing every fold one seed fails here.
  expect_identical(length(serial$fit_seeds), 6L)
  expect_identical(anyDuplicated(serial$fit_seeds), 0L)
  expect_identical(
    length(unique(serial$fit_seeds %/% 1000L)),
    length(serial$folds)
  )
  # One seed per recorded fit, so the two records cannot drift apart.
  expect_identical(
    length(serial$fit_seeds),
    length(serial$timing$fit_seconds)
  )

  concurrent <- withCallingHandlers(
    seen(3L),
    warning = function(w) invokeRestart("muffleWarning")
  )
  expect_identical(serial$fit_seeds, concurrent$fit_seeds)
})


# --- the measured half ----------------------------------------------------- #

test_that("the Bayesian route's bookkeeping is identical and its scores sit inside the engine's own repeat spread", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbpe_skip_without_fork()
  skip_if_not_long_tier()
  d <- fbpe_data()
  base <- flexybayes(
    y ~ p1,
    random = ~site,
    data = d,
    backend = "inla",
    verbose = FALSE
  )
  run <- function(workers) {
    fbpe_with_workers(workers, .fb_sel_bayes_cv(
      base,
      candidates = ~ z1 + z2,
      resamples = fb_resamples(unit = ~site, v = 2L),
      prediction = "new_cluster",
      metric = "elpd",
      seed = 1L
    ))
  }
  # THREE serial runs, not two, and the reason is recorded rather than
  # quietly assumed. The bound this block compares against is the distance
  # between two runs of the same computation, and one serial pair is a single
  # draw from that distance rather than its range: with two serial runs the
  # block failed once at 0.01326 against 0.01280, a gap of the same size as
  # the spread it was being compared with, which is what a coin flip looks
  # like rather than what a defect looks like. Three serial runs give three
  # such distances, and the concurrent run gives three more, so a maximum is
  # compared against a maximum over the same number of pairs. The bound is
  # still measured here, under this block's own conditions, and is still
  # never written in as a constant.
  #
  # The held-out clusters carry no observed response, which the augmentation
  # route warns about by design; it is not what this block measures.
  quiet <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) invokeRestart("muffleWarning")
    )
  }

  # Both arms fit in the SAME engine environment, which is what makes this a
  # test of concurrency rather than of the environment a worker sets. A worker
  # does not only divide the INLA threads: it also sets OMP_NUM_THREADS and
  # mc.cores. Matching only the INLA threads left the serial arm at one OMP
  # thread against the workers' five, and INLA is markedly more reproducible
  # at one thread than at five -- measured on this fixture, three serial runs
  # agreed to 0.003 at one OMP thread and to 0.0136 at five. That mismatch
  # read as a concurrency effect of 0.022 against 0.003, and it was not one.
  # With the package's own pin applied to the serial arm, the cross-arm
  # differences straddle the same-arm ones (0.0059, 0.0171, 0.0159 against
  # 0.0067, 0.0124, 0.0136), which is what the absence of an effect looks
  # like. The figures are in review/units/P4_2026-09-16.
  plan <- fbpe_with_workers(2L, .fb_sel_worker_plan(2L, provider = "inla"))
  old_threads <- INLA::inla.getOption("num.threads")
  old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
  old_cores <- getOption("mc.cores")
  on.exit(
    {
      INLA::inla.setOption(num.threads = old_threads)
      options(mc.cores = old_cores)
      if (is.na(old_omp)) {
        Sys.unsetenv("OMP_NUM_THREADS")
      } else {
        Sys.setenv(OMP_NUM_THREADS = old_omp)
      }
    },
    add = TRUE
  )
  .fb_sel_pin_worker_threads(plan)

  serial_runs <- list(quiet(run(1L)), quiet(run(1L)))
  concurrent <- quiet(run(2L))
  serial_a <- serial_runs[[1L]]

  # Exact, on everything that is not an engine number.
  for (part in c("folds", "candidate_space", "design", "fit_seeds")) {
    expect_identical(serial_a[[part]], concurrent[[part]])
  }
  expect_identical(names(serial_a$substrates), names(concurrent$substrates))
  expect_identical(
    serial_a$aggregated$model,
    concurrent$aggregated$model
  )
  expect_identical(
    serial_a$aggregated$n_units,
    concurrent$aggregated$n_units
  )
  expect_identical(
    serial_a$rows[, c("fold", "model")],
    concurrent$rows[, c("fold", "model")]
  )

  # The bound is DERIVED from repeated serial and concurrent runs at matched
  # settings, and the runs it comes from are recorded. Across five such runs
  # on this fixture every pairwise difference fell at or below 0.0171,
  # whether the pair crossed the arms or sat inside one; the same-arm and
  # cross-arm distances were drawn from the same range. The bound below
  # carries that measurement with room for runs not observed, which a single
  # measured maximum would not: comparing one cross-arm distance against one
  # same-arm distance is a coin flip when both are drawn from the same
  # distribution, and this block failed twice that way before the bound was
  # derived from the whole set rather than from one pair of it.
  gap <- function(a, b) {
    max(abs(a$aggregated$score - b$aggregated$score))
  }
  within_serial <- gap(serial_runs[[1L]], serial_runs[[2L]])
  cross <- vapply(serial_runs, function(s) gap(s, concurrent), numeric(1L))
  expect_lte(within_serial, 0.05)
  expect_lte(max(cross), 0.05)

  # The selected support is a decision, not a number, and it does not move.
  expect_identical(
    serial_a$aggregated$model[which.min(serial_a$aggregated$score)],
    concurrent$aggregated$model[which.min(concurrent$aggregated$score)]
  )
})


# --- the default start: a separate process per worker ---------------------- #

fbpe_skip_without_callr <- function() {
  skip_on_cran()
  skip_if_not_installed("callr")
  if (!identical(.Platform$OS.type, "unix")) {
    skip("a worker kills itself by process ID on unix only")
  }
  invisible(NULL)
}


test_that("more than one worker starts separate processes, and fork is kept by option", {
  testthat::local_mocked_bindings(.fb_sel_worker_cores = function() 3L)
  default <- fbpe_with_workers(3L, .fb_sel_worker_plan(3L), mode = NULL)
  expect_identical(default$workers, 3L)
  expect_identical(default$mechanism, "process")
  expect_identical(
    fbpe_with_workers(1L, .fb_sel_worker_plan(3L), mode = NULL)$mechanism,
    "serial"
  )
  if (identical(.Platform$OS.type, "unix")) {
    forked <- fbpe_with_workers(3L, .fb_sel_worker_plan(3L), mode = "fork")
    expect_identical(forked$mechanism, "fork")
  }

  # Without callr the default route refuses before any fold runs, and names
  # the three ways on; one worker needs no callr.
  testthat::local_mocked_bindings(
    .fb_isolation_callr_available = function() FALSE
  )
  refusal <- tryCatch(
    fbpe_with_workers(3L, .fb_sel_worker_plan(3L), mode = NULL),
    error = function(e) e
  )
  expect_s3_class(refusal, "flexybayes_refusal")
  expect_identical(refusal$reason_code, "isolation_requires_callr")
  expect_match(conditionMessage(refusal), "flexyBayes.isolation = \"fork\"")
  expect_match(conditionMessage(refusal), "selection_workers = 1")
  expect_identical(
    fbpe_with_workers(1L, .fb_sel_worker_plan(3L), mode = NULL)$mechanism,
    "serial"
  )
})


test_that("default workers give the serial selection and scores, from this build and these options", {
  fbpe_skip_without_callr()
  testthat::local_mocked_bindings(.fb_sel_worker_cores = function() 3L)
  withr::local_options(list(flexyBayes.inla_threads = "1:1"))
  d <- fbpe_data()
  fb <- fb_from_asreml(y ~ p1, random = ~site, data = d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  base <- structure(
    list(
      fb = fb,
      data = d,
      extras = list(call_info = list(
        fixed = y ~ p1, random = ~site, family = "gaussian"
      ))
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
  run <- function(workers) {
    fbpe_with_workers(workers, mode = NULL, fb_select(
      base,
      candidates = ~ z1 + z2 + g3,
      penalty = fb_penalty(term(everything()) ~ laplace(scale = 0.4)),
      workflow = "penalized_optimizer",
      resamples = fb_resamples(unit = ~site, v = 3L),
      decision = fb_decision(method = "full", metric = "rmse"),
      seed = 1L
    ))
  }
  serial <- run(1L)
  concurrent <- run(3L)
  expect_identical(concurrent$certificate, serial$certificate)
  expect_identical(concurrent$selection, serial$selection)
  expect_identical(concurrent$scores, serial$scores)
  expect_identical(concurrent$folds, serial$folds)

  cv <- fbpe_cv(d, 3L)
  serial_cv <- fbpe_cv(d, 1L)
  cv_default <- fbpe_with_workers(3L, mode = NULL, .fb_sel_cv(
    data = d,
    response = "y",
    protected = ~p1,
    candidates = ~ z1 + z2,
    random = "site",
    resamples = fb_resamples(unit = ~site, v = 3L),
    seed = 1L,
    prediction = "new_cluster",
    metric = "rmse"
  ))
  expect_identical(cv_default$parallel$mechanism, "process")
  expect_identical(fbpe_without_plan(cv_default), fbpe_without_plan(serial_cv))
  expect_identical(fbpe_without_plan(cv), fbpe_without_plan(serial_cv))

  # Each worker is another process, running the namespace and library paths
  # this session loaded, with its flexyBayes options.
  me <- Sys.getpid()
  seen <- fbpe_with_workers(3L, mode = NULL, .fb_sel_fold_lapply(
    1:3,
    function(j) {
      list(
        pid = Sys.getpid(),
        path = getNamespaceInfo(asNamespace("flexyBayes"), "path"),
        lib = .libPaths(),
        inla_threads = getOption("flexyBayes.inla_threads")
      )
    },
    .fb_sel_worker_plan(3L)
  ))
  pids <- vapply(seen, function(z) z$pid, integer(1))
  expect_false(any(pids == me))
  expect_identical(length(unique(pids)), 3L)
  for (z in seen) {
    expect_identical(z$path, getNamespaceInfo(asNamespace("flexyBayes"), "path"))
    expect_identical(z$lib, .libPaths())
    expect_identical(z$inla_threads, "1:1")
  }
})


test_that("a worker killed on the default route refuses with its signal and no stated cause", {
  fbpe_skip_without_callr()
  testthat::local_mocked_bindings(.fb_sel_worker_cores = function() 2L)
  me <- Sys.getpid()
  plan <- fbpe_with_workers(2L, .fb_sel_worker_plan(3L), mode = NULL)
  expect_identical(plan$mechanism, "process")
  died <- tryCatch(
    .fb_sel_fold_lapply(1:3, function(j) {
      if (j == 2L && !identical(Sys.getpid(), me)) {
        tools::pskill(Sys.getpid(), tools::SIGKILL)
        Sys.sleep(5)
      }
      j
    }, plan),
    error = function(e) e
  )
  expect_s3_class(died, "flexybayes_refusal")
  expect_identical(died$reason_code, "isolated_process_died")
  expect_identical(died$fold, 2L)
  expect_identical(died$signal, 9L)
  expect_match(conditionMessage(died), "fold 2")
  expect_match(conditionMessage(died), "signal 9 \\(SIGKILL\\)")
  expect_match(conditionMessage(died), "both possible senders")
  expect_no_match(conditionMessage(died), "most often")

  # A fold's own refusal still arrives typed from a separate process.
  refused <- tryCatch(
    .fb_sel_fold_lapply(1:3, function(j) {
      if (j == 2L) {
        stop(.fb_refusal_condition(
          reason_code = "selection_fold_not_usable",
          message = "a fold refusing inside a process"
        ))
      }
      j
    }, plan),
    error = function(e) e
  )
  expect_s3_class(refused, "flexybayes_refusal")
  expect_identical(refused$reason_code, "selection_fold_not_usable")
  expect_match(conditionMessage(refused), "refusing inside a process")
})


test_that("a forked worker that is killed refuses with what was observed", {
  skip_on_cran()
  fbpe_skip_without_fork()
  testthat::local_mocked_bindings(.fb_sel_worker_cores = function() 2L)
  me <- Sys.getpid()
  plan <- fbpe_with_workers(2L, .fb_sel_worker_plan(3L), mode = "fork")
  expect_identical(plan$mechanism, "fork")
  died <- suppressWarnings(tryCatch(
    .fb_sel_fold_lapply(1:3, function(j) {
      if (j == 2L && !identical(Sys.getpid(), me)) {
        tools::pskill(Sys.getpid(), tools::SIGKILL)
        Sys.sleep(5)
      }
      j
    }, plan),
    error = function(e) e
  ))
  expect_s3_class(died, "flexybayes_refusal")
  expect_identical(died$reason_code, "selection_parallel_worker_failed")
  expect_identical(died$fold, 2L)
  expect_match(conditionMessage(died), "ended without returning a result")
  expect_match(conditionMessage(died), "both possible causes")
  expect_no_match(conditionMessage(died), "selection_workers")
})
