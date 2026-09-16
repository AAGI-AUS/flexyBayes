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
# What is NOT done here, and where it is done instead. No block kills a
# worker. A worker dying leaves a bare NULL in its element with no condition
# attached, which was measured twice by two different means and recorded in
# review/units/P4_2026-09-16; the block below drives the function that turns
# that NULL into a named refusal, rather than killing a process inside the
# test runner.
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


# Run `expr` with a worker count set, and put the option back afterwards.
fbpe_with_workers <- function(workers, expr) {
  old <- getOption("flexyBayes.selection_workers")
  on.exit(options(flexyBayes.selection_workers = old), add = TRUE)
  options(flexyBayes.selection_workers = workers)
  force(expr)
}


# Concurrency here is fork-based. On a platform that does not fork the folds
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
    expect_identical(capped$workers, 2L)
    expect_true(is.character(capped$note) && nzchar(capped$note))

    plan <- .fb_sel_worker_plan(5L, provider = "inla")
    expect_identical(plan$workers, min(3L, plan$cores))
    # Workers times engine threads does not exceed the machine.
    expect_lte(plan$workers * plan$threads_per_worker, plan$cores)
    expect_gte(plan$threads_per_worker, 1L)
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
  # The refusal points at the lever that fixes it, which is memory and not
  # the core count.
  expect_match(conditionMessage(failed), "selection_workers")

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
  # the substrates the fold loop builds, so the grouped check is run where
  # the concurrency is. It does not go through fb_select() because in this
  # tree it cannot: the optimiser route resolves its penalty through the brms
  # compiler, which refuses every grouped family by name
  # (selection_penalty_family_optimizer_only). That is measured and recorded
  # in review/units/P4_2026-09-16, and belongs to the unit that owns the
  # grouped families rather than to this one.
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


# --- the seeds, and the mutant that shares one ----------------------------- #

test_that("every fold and model gets its own seed, the same seeds either way", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbpe_skip_without_fork()
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
  serial_runs <- list(quiet(run(1L)), quiet(run(1L)), quiet(run(1L)))
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

  # Measured, not chosen: the engine's own repeat spread is the bound, and it
  # is computed here from the serial runs rather than written in.
  gap <- function(a, b) {
    max(abs(a$aggregated$score - b$aggregated$score))
  }
  serial_pairs <- c(
    gap(serial_runs[[1L]], serial_runs[[2L]]),
    gap(serial_runs[[1L]], serial_runs[[3L]]),
    gap(serial_runs[[2L]], serial_runs[[3L]])
  )
  concurrent_pairs <- vapply(
    serial_runs,
    function(s) gap(s, concurrent),
    numeric(1L)
  )
  expect_lte(
    max(concurrent_pairs),
    max(max(serial_pairs), .Machine$double.eps)
  )

  # The selected support is a decision, not a number, and it does not move.
  expect_identical(
    serial_a$aggregated$model[which.min(serial_a$aggregated$score)],
    concurrent$aggregated$model[which.min(concurrent$aggregated$score)]
  )
})
