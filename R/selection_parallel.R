# Concurrency for the fold loops, under a reproducibility envelope.
#
# `fb_select()` refits the base model inside every fold and scores the
# held-out rows on posterior predictive draws. Measured on a package fixture
# at the pilot's shape -- 10,890 rows, 82 candidates, 3 folds -- one call
# costs 103 s against a 2.9 s base fit, and posterior sampling is about 82
# per cent of each fit (review/units/P1-FIX_2026-09-16). The folds do not
# read each other, so that wall time is a scheduling property rather than a
# property of the estimator, and this file runs the folds at the same time.
#
# What must not change, and why it does not.
#
#   The answer. Every site in a fold that draws is seeded explicitly: the
#   fit seed is `base_seed + 1000 * fold + model`, both posterior samplers
#   receive it, and the brms predictor runs inside `.fb_sel_with_seed()`.
#   Nothing in a fold reads the random stream it inherits, so a fold
#   computes the same thing whichever process runs it and in whatever order.
#   That is the property the envelope rests on, and the tests assert it
#   rather than assume it.
#
#   The bookkeeping. A worker returns its fold's results and the parent
#   assembles them in FOLD order, never in completion order, so the scored
#   rows, the per-model tables, the substrates, the fit seconds and the fit
#   seeds are built in the order a serial run builds them.
#
#   The certificate. The worker count is a property of the machine that ran
#   the call, not of the estimate, so it is recorded beside the timing on the
#   fold object and is deliberately absent from the certificate. A
#   certificate naming its worker count could not be compared between a
#   serial and a parallel run, which is the comparison this file exists to
#   survive.
#
# Nested parallelism. A brms fold would run its chains on `mc.cores`, so f
# workers times t engine threads is f * t busy processes on a machine with c
# cores. Each worker is therefore pinned to `cores %/% workers` OpenMP threads
# and one chain process before it fits anything, and the plan records what was
# set. INLA folds run at option flexyBayes.inla_threads, "1:1" by default.
#
# Memory, and why the worker count is not read from the core count. On the
# same fixture, aggregate memory is driven by the rows rather than by the
# fold count: about 4 GiB per fold at 10,890 rows, whether 3 folds ran or 6.
# Concurrency multiplies that, so on a 32 GiB machine the useful worker count
# is set by memory well below the core count. The option therefore defaults
# to one worker and a caller raises it against a measured figure.
#
# How workers start. More than one worker runs each fold in a separate R
# process (R/process_isolation.R). A forked worker inherits the session's
# loaded engine libraries after fork, and a forked worker of a session with
# INLA, OpenMP and TBB loaded has died with SIGABRT in R's garbage collector
# with memory to spare. Option flexyBayes.isolation = "fork" keeps forking.
#
# A worker that dies. A forked child that ends without returning leaves no
# condition to re-raise and no exit status. Every child therefore returns
# its value inside a marked wrapper, and a result without that marker
# refuses by name rather than leaving one fold quietly missing from the
# aggregation. A separate process reports its exit status or signal.

# --- the plan -------------------------------------------------------------- #

#' The worker count the caller asked for
#'
#' @return A single positive integer.
#' @noRd
#' @keywords internal
.fb_sel_worker_request <- function() {
  value <- getOption("flexyBayes.selection_workers", 1L)
  if (
    !is.numeric(value) ||
      length(value) != 1L ||
      !is.finite(value) ||
      value < 1 ||
      value != round(value)
  ) {
    stop(
      "option flexyBayes.selection_workers must be a single whole number of ",
      "fold workers, one or more. It is: ",
      .fb_prior_deparse(value),
      ".",
      call. = FALSE
    )
  }
  as.integer(value)
}


#' The cores this machine reports
#'
#' @return A single positive integer; 1 where the count cannot be read.
#' @noRd
#' @keywords internal
.fb_sel_worker_cores <- function() {
  cores <- suppressWarnings(as.integer(parallel::detectCores()))
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    return(1L)
  }
  cores
}


#' How a fold loop will be run
#'
#' The requested worker count is capped by the number of folds and by the
#' cores. More than one worker runs each fold in a separate R process, unless
#' option `flexyBayes.isolation` is `"fork"`, which forks them and drops to
#' one worker where the platform cannot fork. The engine threads each
#' worker may use are divided out of the cores at the same time, so the
#' product of the two does not exceed the machine.
#'
#' @param n_jobs The number of folds to run.
#' @param provider The engine that will fit inside a fold, or `NULL` where
#'   no engine does (the closed-form benchmark route).
#' @return A record of the plan: what was asked for, what will run, the
#'   mechanism, the engine threads per worker, and a note where the request
#'   was cut.
#' @noRd
#' @keywords internal
.fb_sel_worker_plan <- function(n_jobs, provider = NULL) {
  requested <- .fb_sel_worker_request()
  cores <- .fb_sel_worker_cores()
  n_jobs <- max(1L, as.integer(n_jobs))
  forkable <- identical(.Platform$OS.type, "unix")
  mode <- .fb_isolation_mode()
  # Under options(flexyBayes.isolation = "process") every fold runs in a
  # child R process of its own (R/process_isolation.R), which needs no fork.
  isolated <- .fb_isolation_active("fb_select()")

  workers <- max(1L, min(requested, n_jobs, cores))
  if (workers > 1L && identical(mode, "none")) {
    .fb_sel_require_callr(workers)
    isolated <- TRUE
  }
  note <- NA_character_
  if (workers > 1L && !forkable && !isolated) {
    note <- paste0(
      "this platform does not fork, so the ",
      n_jobs,
      " folds ran one at a time"
    )
    workers <- 1L
  } else if (requested > workers) {
    note <- paste0(
      "the request of ",
      requested,
      " was cut to ",
      workers,
      ": ",
      n_jobs,
      " folds on ",
      cores,
      " cores"
    )
  }

  list(
    requested = requested,
    workers = workers,
    mechanism = if (isolated) {
      "process"
    } else if (workers > 1L) {
      "fork"
    } else {
      "serial"
    },
    cores = cores,
    n_jobs = n_jobs,
    # Left NA on the serial path, where nothing is pinned and the engines
    # run exactly as they run today.
    threads_per_worker = if (workers > 1L) {
      max(1L, cores %/% workers)
    } else {
      NA_integer_
    },
    provider = provider %||% NA_character_,
    note = note
  )
}


#' Refuse concurrent fold workers when callr is missing
#'
#' @param workers The worker count the plan arrived at.
#' @return `NULL`, invisibly; refuses when callr cannot be loaded.
#' @noRd
#' @keywords internal
.fb_sel_require_callr <- function(workers) {
  if (.fb_isolation_callr_available()) {
    return(invisible(NULL))
  }
  stop(.fb_refusal_condition(
    reason_code = "isolation_requires_callr",
    message = paste0(
      "fb_select(): ",
      workers,
      " fold workers run as separate R processes started through callr, ",
      "and callr is not installed. Install it with install.packages(",
      "\"callr\"), set options(flexyBayes.selection_workers = 1) to run the ",
      "folds one at a time, or set options(flexyBayes.isolation = \"fork\") ",
      "to fork the workers from this session."
    ),
    entry_point = "fb_select()"
  ))
}


#' Pin one worker's engine threads
#'
#' Called inside the worker, before it fits. The settings are process-local:
#' a forked child's options and environment are its own, so the parent's
#' session is left as the caller set it.
#'
#' @param plan A `.fb_sel_worker_plan()` record.
#' @return `NULL`, invisibly.
#' @noRd
#' @keywords internal
.fb_sel_pin_worker_threads <- function(plan) {
  threads <- plan$threads_per_worker
  if (length(threads) != 1L || is.na(threads)) {
    return(invisible(NULL))
  }
  Sys.setenv(OMP_NUM_THREADS = as.character(threads))
  # One chain process per worker: brms reads `cores` from this option, and a
  # worker fanning out its chains would multiply the workers again.
  options(mc.cores = 1L)
  # INLA's threads are not divided here: every fit takes option
  # flexyBayes.inla_threads ("1:1" by default), which the worker inherits,
  # because any other plan makes a fold's variance components depend on the run.
  invisible(NULL)
}


# --- running the folds ----------------------------------------------------- #

#' One worker's result, or the reason there is none
#'
#' @param result One element of the `parallel::mclapply()` result.
#' @param index The fold's position, for the message.
#' @return The value the fold function returned.
#' @noRd
#' @keywords internal
.fb_sel_worker_value <- function(result, index) {
  if (is.list(result) && identical(names(result), "fb_sel_worker_value")) {
    return(result[["fb_sel_worker_value"]])
  }
  # A fold that refused, or failed, inside a live child: its own condition
  # travels back and is raised here, so a typed refusal stays typed and
  # names the same reason it would name on the serial path.
  condition <- attr(result, "condition")
  if (inherits(condition, "condition")) {
    stop(condition)
  }
  stop(.fb_refusal_condition(
    reason_code = "selection_parallel_worker_failed",
    message = paste0(
      "fb_select(): the forked worker running fold ",
      index,
      " ended without returning a result and left no error; a forked ",
      "worker's exit status and signal are not reported. The operating ",
      "system ending it to reclaim memory and a crash in compiled code ",
      "inside it are both possible causes. Nothing is aggregated from a ",
      "partial set of folds. With options(flexyBayes.isolation) unset, ",
      "each fold runs in a separate R process and a fold that ends this ",
      "way is reported with its exit status or signal."
    ),
    fold = index
  ))
}


#' Run a fold loop, serially or in workers
#'
#' @param jobs The fold positions to run.
#' @param fn A function of one fold position, returning that fold's results.
#' @param plan A `.fb_sel_worker_plan()` record.
#' @param warm_first `TRUE` runs the first fold in the parent before the
#'   rest fork. brms compiles a Stan program per model shape and caches it
#'   by the hash of that program; children compiling the same shape at the
#'   same time write the same cache entry, which is the concurrent-load
#'   failure recorded as N29. The first fold compiles every shape the loop
#'   needs -- the models differ by their fixed terms, not by the fold -- so
#'   the children that follow read a warm cache and compile nothing.
#' @return The folds' results, in fold order.
#' @noRd
#' @keywords internal
.fb_sel_fold_lapply <- function(jobs, fn, plan, warm_first = FALSE) {
  if (identical(plan$mechanism, "process")) {
    return(.fb_sel_isolated_lapply(jobs, fn, plan, warm_first = warm_first))
  }
  if (identical(plan$mechanism, "serial") || length(jobs) <= 1L) {
    return(lapply(jobs, fn))
  }
  if (isTRUE(warm_first)) {
    return(c(
      list(fn(jobs[[1L]])),
      .fb_sel_fork_lapply(jobs[-1L], fn, plan)
    ))
  }
  .fb_sel_fork_lapply(jobs, fn, plan)
}


#' The forked half of `.fb_sel_fold_lapply()`
#'
#' `mc.preschedule = FALSE` gives each fold its own child, so one fold's
#' memory returns to the machine as soon as that fold finishes.
#' `mc.set.seed = FALSE` leaves the child's random stream as the parent's:
#' every draw inside a fold is seeded from the fold's own seed, so a
#' per-child stream would be a source of variation where the design has none.
#'
#' @inheritParams .fb_sel_fold_lapply
#' @return The folds' results, in fold order.
#' @noRd
#' @keywords internal
.fb_sel_fork_lapply <- function(jobs, fn, plan) {
  wrapped <- function(job) {
    .fb_sel_pin_worker_threads(plan)
    list(fb_sel_worker_value = fn(job))
  }
  results <- parallel::mclapply(
    jobs,
    wrapped,
    mc.cores = plan$workers,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE,
    mc.allow.recursive = FALSE
  )
  if (length(results) != length(jobs)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_parallel_worker_failed",
      message = paste0(
        "fb_select(): ",
        length(jobs),
        " folds were sent to the workers and ",
        length(results),
        " results came back. Nothing is aggregated from a partial set of ",
        "folds."
      ),
      fold = NA_integer_
    ))
  }
  lapply(
    seq_along(results),
    function(i) .fb_sel_worker_value(results[[i]], jobs[[i]])
  )
}
