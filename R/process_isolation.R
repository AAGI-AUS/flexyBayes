# Engine fits and fold loops in a separate R process.
#
# `options(flexyBayes.isolation = "process")` runs the engine fit behind
# flexybayes() (and so behind fb(), fb_inla() and fb_brms()), the fit behind
# proceed(), and each fold of fb_select()'s fold loops in a separate R process
# started through callr. The default, "none", runs everything in the session,
# exactly as before.
#
# How a session dies, path by path, and what isolation changes. Read from the
# code and measured on this machine (review/units/U-ISOLATE_2026-09-19,
# probe 01), because the design follows from it:
#
#   * INLA fits in its own program, `inla.run`, started with system(). A
#     signal that ends that program (a segmentation fault, or the operating
#     system reclaiming memory) comes back to R as a non-zero exit status;
#     INLA retries once, and emit_inla() turns a second failure into the
#     typed refusal `inla_program_failed`. The R session survives this path
#     with or without isolation. The INLA R package has no compiled code.
#   * brms samples in CmdStan's own program under cmdstanr, so a sampler
#     crash there is also an error in R. Under rstan, which brms uses when
#     `options(brms.backend)` is unset, the sampler runs inside the R
#     process, and a crash in it ends the session.
#   * R itself: when the operating system kills the R process under memory
#     pressure, or compiled code running inside it (a BLAS call on a dense
#     matrix, a Stan sampler under rstan) crashes, the session ends with no
#     condition to catch. The selection layer's dense covariance and the
#     proximal solvers run inside R, so this is the path a wide selection
#     call takes.
#   * fb_select()'s forked fold workers (R/selection_parallel.R): a worker
#     killed from outside leaves a bare NULL, which is refused by name, and
#     the session survives.
#
# Isolation moves the third path into a child process, where it becomes a
# child that ended with a signal or an exit status. That is reported as the
# typed refusal `isolated_process_died`, naming the engine, the entry point
# and the signal, and the session goes on. It does not prevent the death, and
# it does not retry: a retry of the same call is attempted only where a
# measurement shows that retrying an identical call ever succeeds, and the
# measurement in the unit's evidence found none.
#
# What the child is given, so that it fits what the session would have fitted.
#
#   The code. The namespace this session loaded, found from the namespace
#   itself: the source tree under pkgload, otherwise the library the
#   namespace was loaded from. A child that loaded whichever build is
#   installed first could run different code under the same name.
#   The inputs. The same function and arguments, serialised by the parent
#   after the namespace is loaded in the child, so that closures and
#   formulas carry their environments across.
#   The settings. The options that change a fit or its printed output
#   (flexyBayes, brms and cmdstanr options, `mc.cores`, `contrasts`,
#   `na.action`, and the printing options), INLA's options when INLA is
#   loaded, the CmdStan path when cmdstanr is loaded, and the environment
#   variables as they stand in the session (the child reads no startup
#   files, so a thread count set with Sys.setenv() reaches it).
#   The state. The random-number state goes in and the child's comes back,
#   so a draw inside an isolated call is the draw an in-process call makes;
#   the once-per-session notes go in and come back the same way.
#
# What comes back. The value, and every message and warning in the order
# they were raised, re-signalled in the session with its class; printed
# output in its place among them; an error raised in the child re-raised as
# itself, so a typed refusal keeps its reason code.

# --- the option ---------------------------------------------------------- #

#' The isolation mode the caller asked for
#'
#' @return `"none"` or `"process"`.
#' @noRd
#' @keywords internal
.fb_isolation_mode <- function() {
  value <- getOption("flexyBayes.isolation", "none")
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !value %in% c("none", "process")
  ) {
    stop(
      "option flexyBayes.isolation must be \"none\" or \"process\". It is: ",
      .fb_prior_deparse(value),
      ".",
      call. = FALSE
    )
  }
  value
}


#' Whether callr can be loaded
#'
#' A function of its own so that its absence can be simulated in a test.
#'
#' @return `TRUE` or `FALSE`.
#' @noRd
#' @keywords internal
.fb_isolation_callr_available <- function() {
  requireNamespace("callr", quietly = TRUE)
}


#' Whether this entry point runs its fit in a separate process
#'
#' Reads the option and, when isolation is asked for, refuses by name if
#' callr is not installed, before anything is fitted.
#'
#' @param entry_point The public verb asking, for the refusal.
#' @return `TRUE` when the fit runs in a separate process.
#' @noRd
#' @keywords internal
.fb_isolation_active <- function(entry_point) {
  if (!identical(.fb_isolation_mode(), "process")) {
    return(FALSE)
  }
  if (!.fb_isolation_callr_available()) {
    stop(.fb_refusal_condition(
      reason_code = "isolation_requires_callr",
      message = paste0(
        entry_point,
        ": options(flexyBayes.isolation = \"process\") runs each fit in a ",
        "separate R process through the callr package, and callr is not ",
        "installed. Install it with install.packages(\"callr\"), or set ",
        "options(flexyBayes.isolation = \"none\") to fit in this session."
      ),
      entry_point = entry_point
    ))
  }
  TRUE
}


# --- what the child is given --------------------------------------------- #

#' How the child loads the namespace this session loaded
#'
#' @return A list: `dev` (loaded from a source tree by pkgload), `path` (the
#'   tree, or the installed package directory) and `lib` (its library).
#' @noRd
#' @keywords internal
.fb_isolation_load_spec <- function() {
  ns <- asNamespace("flexyBayes")
  path <- getNamespaceInfo(ns, "path")
  list(
    dev = exists(".__DEVTOOLS__", envir = ns, inherits = FALSE),
    path = path,
    lib = dirname(path)
  )
}


#' The session settings and state a child starts from
#'
#' @return A list of options, INLA options, the CmdStan path and output
#'   directory, the random-number state and the once-per-session flags.
#' @noRd
#' @keywords internal
.fb_isolation_state <- function() {
  opts <- options()
  keep <- grepl(
    paste0(
      "^(flexyBayes[.]|brms[.]|cmdstanr|mc[.]cores$|contrasts$|",
      "na[.]action$|digits$|scipen$|OutDec$|width$)"
    ),
    names(opts)
  )
  cmdstan_path <- NULL
  if (isNamespaceLoaded("cmdstanr")) {
    cmdstan_path <- tryCatch(
      getExportedValue("cmdstanr", "cmdstan_path")(),
      error = function(e) NULL
    )
  }
  list(
    options = opts[keep],
    inla = if (isNamespaceLoaded("INLA")) INLA::inla.getOption() else NULL,
    cmdstan_path = cmdstan_path,
    # A Stan program compiled in the session is found again by the child: the
    # child writes its programs where the session writes them.
    stan_file_dir = getOption("cmdstanr_write_stan_file_dir", tempdir()),
    seed = if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      NULL
    },
    emit_state = as.list(.flexybayes_emit_state)
  )
}


#' Apply a session's settings and state inside the child
#'
#' @param state A `.fb_isolation_state()` record.
#' @return `NULL`, invisibly.
#' @noRd
#' @keywords internal
.fb_isolation_apply_state <- function(state) {
  options(state$options)
  options(cmdstanr_write_stan_file_dir = state$stan_file_dir)
  if (!is.null(state$inla)) {
    inla_opts <- state$inla[!vapply(state$inla, is.null, logical(1))]
    do.call(INLA::inla.setOption, inla_opts)
  }
  if (!is.null(state$cmdstan_path)) {
    suppressMessages(
      getExportedValue("cmdstanr", "set_cmdstan_path")(state$cmdstan_path)
    )
  }
  .fb_isolation_set_seed(state$seed)
  .fb_isolation_set_emit_state(state$emit_state)
  invisible(NULL)
}


#' Set, or remove, the random-number state in the global environment
#'
#' @param seed A `.Random.seed` vector, or `NULL` for none.
#' @return `NULL`, invisibly.
#' @noRd
#' @keywords internal
.fb_isolation_set_seed <- function(seed) {
  if (is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", seed, envir = globalenv())
  }
  invisible(NULL)
}


#' Set the once-per-session flags
#'
#' @param flags A named list of logical flags.
#' @return `NULL`, invisibly.
#' @noRd
#' @keywords internal
.fb_isolation_set_emit_state <- function(flags) {
  for (key in names(flags)) {
    .emit_state_set(key, isTRUE(flags[[key]]))
  }
  invisible(NULL)
}


# --- the child ------------------------------------------------------------ #

#' The child's entry point
#'
#' callr runs this function in the child with its environment set to the
#' child's global environment, so it names nothing from this namespace until
#' it has loaded it: the namespace this session loaded, found through
#' `load`. Only then are the payload's closures read, so that their
#' environments attach to that namespace.
#'
#' @param payload_path The file holding the call.
#' @param load A `.fb_isolation_load_spec()` record.
#' @return The child's record, from `.fb_isolation_child_run()`.
#' @noRd
#' @keywords internal
.fb_isolation_child_entry <- function(payload_path, load) {
  if (isTRUE(load$dev)) {
    load_all <- getExportedValue("pkgload", "load_all")
    suppressPackageStartupMessages(suppressMessages(load_all(
      load$path,
      quiet = TRUE,
      helpers = FALSE,
      attach_testthat = FALSE,
      export_all = FALSE
    )))
  } else {
    loadNamespace("flexyBayes", lib.loc = load$lib)
  }
  run <- get(".fb_isolation_child_run", envir = asNamespace("flexyBayes"))
  run(payload_path)
}


#' Run the call inside the child and record everything that comes back
#'
#' @param payload_path The file holding the call: `fn_path` (the function),
#'   `args`, `quote`, `state` and `pin`.
#' @return A list: `ok`, `value` or `condition`, the `events` in the order
#'   they happened, and the child's random-number state and flags.
#' @noRd
#' @keywords internal
.fb_isolation_child_run <- function(payload_path) {
  payload <- readRDS(payload_path)
  fn <- readRDS(payload$fn_path)
  .fb_isolation_apply_state(payload$state)
  # The child is the isolated process; it does not start another.
  options(flexyBayes.isolation = "none")
  if (!is.null(payload$pin)) {
    .fb_sel_pin_worker_threads(payload$pin)
  }

  events <- list()
  printed <- character(0)
  out <- textConnection("printed", open = "w", local = TRUE)
  sink(out)
  closed <- FALSE
  n_seen <- 0L
  # Complete lines printed so far; once the connection is closed, `printed`
  # also holds a last line that ended without a newline.
  flush_printed <- function() {
    lines <- if (closed) printed else textConnectionValue(out)
    if (length(lines) > n_seen) {
      events[[length(events) + 1L]] <<- list(
        type = "output",
        text = lines[seq.int(n_seen + 1L, length(lines))]
      )
      n_seen <<- length(lines)
    }
  }
  result <- tryCatch(
    withCallingHandlers(
      list(
        ok = TRUE,
        value = do.call(fn, payload$args, quote = isTRUE(payload$quote))
      ),
      message = function(m) {
        flush_printed()
        events[[length(events) + 1L]] <<- list(type = "message", condition = m)
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        flush_printed()
        events[[length(events) + 1L]] <<- list(type = "warning", condition = w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) list(ok = FALSE, condition = e)
  )
  sink()
  close(out)
  closed <- TRUE
  flush_printed()

  result$events <- events
  has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  result$seed <- if (has_seed) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
  result$emit_state <- as.list(.flexybayes_emit_state)
  result$pid <- Sys.getpid()
  result
}


# --- the parent ----------------------------------------------------------- #

#' Start one child
#'
#' @param fn_path A file holding the function the child calls.
#' @param args Its arguments.
#' @param pin A `.fb_sel_worker_plan()` record whose engine threads the child
#'   pins, or `NULL`.
#' @param quote Passed to `do.call()` in the child. `TRUE` hands a call object
#'   among `args` over as a value, as a direct call hands over its argument;
#'   `FALSE` evaluates it, as `do.call()` does by default.
#' @return A handle: the process, its PID and the files it uses.
#' @noRd
#' @keywords internal
.fb_isolation_start <- function(fn_path, args = list(), pin = NULL,
                                quote = TRUE) {
  payload_path <- tempfile("fb_isolation_", fileext = ".rds")
  saveRDS(
    list(
      fn_path = fn_path,
      args = args,
      quote = quote,
      state = .fb_isolation_state(),
      pin = pin
    ),
    payload_path,
    compress = FALSE
  )
  out_path <- tempfile("fb_isolation_", fileext = ".out")
  err_path <- tempfile("fb_isolation_", fileext = ".err")
  proc <- callr::r_bg(
    .fb_isolation_child_entry,
    args = list(payload_path = payload_path, load = .fb_isolation_load_spec()),
    stdout = out_path,
    stderr = err_path,
    # No startup files: the child's environment is the session's as it
    # stands, and no profile runs code the session did not.
    cmdargs = c("--slave", "--no-save", "--no-restore", "--no-environ"),
    system_profile = FALSE,
    user_profile = FALSE,
    # The child ends if this session ends, and whatever it started can be
    # found again by the marker every descendant carries.
    supervise = TRUE,
    cleanup_tree = TRUE
  )
  list(
    proc = proc,
    pid = proc$get_pid(),
    files = c(payload_path, out_path, err_path),
    out_path = out_path,
    err_path = err_path
  )
}


#' The readable end of a text file, for a message
#'
#' @param path A file.
#' @param n The number of characters kept from its end.
#' @return A string, possibly empty.
#' @noRd
#' @keywords internal
.fb_isolation_tail <- function(path, n = 600L) {
  txt <- tryCatch(
    paste(readLines(path, warn = FALSE), collapse = "\n"),
    error = function(e) ""
  )
  if (nchar(txt) > n) {
    txt <- substr(txt, nchar(txt) - n + 1L, nchar(txt))
  }
  trimws(txt)
}


#' The refusal for a child that ended without returning
#'
#' @param status The exit status; a negative value is the signal that ended
#'   the process.
#' @param pid The child's process ID.
#' @param engine The engine the child was running, or `"none"`.
#' @param entry_point The public verb.
#' @param fold The fold the child was running, or `NA`.
#' @param stderr_tail What the child wrote to its error stream last.
#' @return A refusal condition.
#' @noRd
#' @keywords internal
.fb_isolation_died_condition <- function(status,
                                         pid,
                                         engine,
                                         entry_point,
                                         fold = NA_integer_,
                                         stderr_tail = "") {
  status <- as.integer(status)
  signal <- if (!is.na(status) && status < 0L) -status else NA_integer_
  signal_names <- c(
    "9" = " (SIGKILL)",
    "11" = " (SIGSEGV)",
    "6" = " (SIGABRT)",
    "15" = " (SIGTERM)"
  )
  named <- !is.na(signal) && as.character(signal) %in% names(signal_names)
  signal_name <- if (named) {
    signal_names[[as.character(signal)]]
  } else {
    ""
  }
  doing_labels <- c(
    inla = "fitting with the INLA engine",
    brms = "fitting with the brms engine",
    auto = "fitting with the engine backend = \"auto\" chose",
    none = "running the call"
  )
  doing <- if (engine %in% names(doing_labels)) {
    doing_labels[[engine]]
  } else {
    paste0("fitting with ", engine)
  }
  how <- if (!is.na(signal)) {
    paste0("it was ended by signal ", signal, signal_name)
  } else {
    paste0("it exited with status ", status)
  }
  cause <- if (identical(signal, 9L)) {
    paste0(
      " A process ended by SIGKILL is most often one the operating system ",
      "stopped to reclaim memory."
    )
  } else if (identical(signal, 11L) || identical(signal, 6L)) {
    " That signal is a crash in compiled code running inside the process."
  } else {
    ""
  }
  where <- if (is.na(fold)) {
    ""
  } else {
    paste0(" while it was running fold ", fold)
  }
  stop_text <- if (nzchar(stderr_tail)) {
    paste0(" Its last output: ", stderr_tail)
  } else {
    ""
  }
  .fb_refusal_condition(
    reason_code = "isolated_process_died",
    message = paste0(
      entry_point,
      ": the separate R process ",
      doing,
      where,
      " ended before it returned: ",
      how,
      ". Nothing was returned from it, and this session is unaffected.",
      cause,
      stop_text
    ),
    engine = engine,
    entry_point = entry_point,
    fold = as.integer(fold),
    exit_status = status,
    signal = signal,
    pid = as.integer(pid)
  )
}


#' Wait for one child and read what it left
#'
#' @param handle A `.fb_isolation_start()` handle.
#' @inheritParams .fb_isolation_died_condition
#' @return The child's record from `.fb_isolation_child_run()`, with the text
#'   it wrote outside R's own output (an engine's program writing to the
#'   terminal) as `raw_out` and `raw_err`. A child that died raises the
#'   refusal instead.
#' @noRd
#' @keywords internal
.fb_isolation_harvest <- function(handle, engine, entry_point,
                                  fold = NA_integer_) {
  on.exit(unlink(handle$files), add = TRUE)
  handle$proc$wait()
  status <- handle$proc$get_exit_status()
  if (is.na(status) || status != 0L) {
    # The processes the child started would otherwise outlive it and keep
    # their memory. They carry the child's marker, so only they are ended.
    tryCatch(handle$proc$kill_tree(), error = function(e) NULL)
    if (!is.na(status) && status > 0L) {
      # A child that failed before the call ran (it could not load the
      # namespace, say) leaves its own error; that is not a death.
      started <- tryCatch(handle$proc$get_result(), error = function(e) e)
      parent <- if (inherits(started, "error")) started$parent else NULL
      if (inherits(parent, "error")) {
        stop(
          entry_point,
          ": flexyBayes could not start the separate R process for this ",
          "fit: ",
          conditionMessage(parent),
          call. = FALSE
        )
      }
    }
    stop(.fb_isolation_died_condition(
      status,
      handle$pid,
      engine,
      entry_point,
      fold = fold,
      stderr_tail = .fb_isolation_tail(handle$err_path)
    ))
  }
  record <- handle$proc$get_result()
  record$raw_out <- tryCatch(
    readLines(handle$out_path, warn = FALSE),
    error = function(e) character(0)
  )
  record$raw_err <- tryCatch(
    readLines(handle$err_path, warn = FALSE),
    error = function(e) character(0)
  )
  record
}


#' Deliver a child's record in the session
#'
#' @param record A `.fb_isolation_harvest()` record.
#' @param adopt_seed `TRUE` to move the session's random-number stream to
#'   where the child's ended, as an in-process call would have moved it.
#' @return The child's value; a condition raised in the child is re-raised.
#' @noRd
#' @keywords internal
.fb_isolation_deliver <- function(record, adopt_seed = TRUE) {
  if (isTRUE(adopt_seed)) {
    .fb_isolation_set_seed(record$seed)
  }
  flags <- record$emit_state
  for (key in names(flags)) {
    if (isTRUE(flags[[key]])) {
      .emit_state_set(key, TRUE)
    }
  }
  for (event in record$events) {
    if (identical(event$type, "output")) {
      cat(paste0(event$text, "\n"), sep = "")
    } else if (identical(event$type, "message")) {
      message(event$condition)
    } else if (identical(event$type, "warning")) {
      warning(event$condition)
    }
  }
  if (length(record$raw_out)) {
    cat(paste0(record$raw_out, "\n"), sep = "")
  }
  if (length(record$raw_err)) {
    cat(paste0(record$raw_err, "\n"), sep = "", file = stderr())
  }
  if (!isTRUE(record$ok)) {
    stop(record$condition)
  }
  record$value
}


#' End a child that is still running, with whatever it started
#'
#' The processes are found by the marker callr gives this child's tree, so
#' nothing else is touched, and the child is waited for so that it does not
#' linger as a finished process nobody has collected.
#'
#' @param handle A `.fb_isolation_start()` handle.
#' @return `NULL`, invisibly.
#' @noRd
#' @keywords internal
.fb_isolation_stop <- function(handle) {
  if (handle$proc$is_alive()) {
    tryCatch(handle$proc$kill_tree(), error = function(e) NULL)
    tryCatch(handle$proc$wait(2000), error = function(e) NULL)
  }
  unlink(handle$files)
  invisible(NULL)
}


#' Call a function in a separate R process
#'
#' @param fn The function the child calls.
#' @param args Its arguments.
#' @param engine The engine the call fits with, for a refusal.
#' @param entry_point The public verb, for a refusal.
#' @param with_pid `TRUE` returns `list(value, pid)` rather than the value.
#' @inheritParams .fb_isolation_start
#' @return The value `fn` returned in the child.
#' @noRd
#' @keywords internal
.fb_isolated_call <- function(fn, args = list(), engine, entry_point,
                              with_pid = FALSE, quote = TRUE) {
  fn_path <- tempfile("fb_isolation_fn_", fileext = ".rds")
  on.exit(unlink(fn_path), add = TRUE)
  saveRDS(fn, fn_path, compress = FALSE)
  handle <- .fb_isolation_start(fn_path, args, quote = quote)
  # An interrupt in the session ends the child and what it started.
  on.exit(.fb_isolation_stop(handle), add = TRUE)
  record <- .fb_isolation_harvest(handle, engine, entry_point)
  value <- .fb_isolation_deliver(record)
  if (isTRUE(with_pid)) {
    return(list(value = value, pid = handle$pid))
  }
  value
}


# --- the fold loop -------------------------------------------------------- #

#' Run a fold loop with each fold in its own R process
#'
#' Up to `plan$workers` children run at once. What comes back is delivered
#' in FOLD order, never in completion order, so the printed output, the
#' messages and warnings, and the first failure raised are the ones a serial
#' run gives: a fold that fails stops new folds from starting, the folds
#' before it are waited for, and the lowest failing fold's condition is the
#' one raised. One fold at a time, each child starts from the random-number
#' state the previous one left, as a serial loop does; concurrent children
#' start from the session's state and leave it where it was, as forked
#' workers do.
#'
#' @inheritParams .fb_sel_fold_lapply
#' @return The folds' results, in fold order.
#' @noRd
#' @keywords internal
.fb_sel_isolated_lapply <- function(jobs, fn, plan, warm_first = FALSE) {
  n <- length(jobs)
  engine <- if (is.null(plan$provider) || is.na(plan$provider)) {
    "none"
  } else {
    plan$provider
  }
  fn_path <- tempfile("fb_isolation_fn_", fileext = ".rds")
  saveRDS(fn, fn_path, compress = FALSE)
  running <- list()
  on.exit(
    {
      for (h in running) {
        .fb_isolation_stop(h)
      }
      unlink(fn_path)
    },
    add = TRUE
  )
  sequential <- plan$workers <= 1L
  pin <- if (sequential) NULL else plan
  records <- vector("list", n)
  failures <- vector("list", n)
  done <- rep(FALSE, n)
  next_i <- 1L
  first_failure <- Inf

  harvest <- function(i) {
    h <- running[[as.character(i)]]
    running[[as.character(i)]] <<- NULL
    outcome <- tryCatch(
      list(record = .fb_isolation_harvest(
        h,
        engine,
        "fb_select()",
        fold = jobs[[i]]
      )),
      error = function(e) list(condition = e)
    )
    done[[i]] <<- TRUE
    if (is.null(outcome$condition)) {
      records[[i]] <<- outcome$record
      if (!isTRUE(outcome$record$ok)) {
        first_failure <<- min(first_failure, i)
      }
      if (sequential) {
        # One fold at a time: the next child starts where this one left off.
        .fb_isolation_set_seed(outcome$record$seed)
      }
    } else {
      failures[[i]] <<- outcome$condition
      first_failure <<- min(first_failure, i)
    }
  }

  while (!all(done[seq_len(min(n, first_failure))])) {
    limit <- if (isTRUE(warm_first) && !done[[1L]]) 1L else plan$workers
    while (
      next_i <= n && next_i < first_failure && length(running) < limit
    ) {
      running[[as.character(next_i)]] <- .fb_isolation_start(
        fn_path,
        list(jobs[[next_i]]),
        pin = pin
      )
      next_i <- next_i + 1L
    }
    finished <- names(running)[!vapply(
      running,
      function(h) h$proc$is_alive(),
      logical(1)
    )]
    if (length(finished) == 0L) {
      Sys.sleep(0.05)
      next
    }
    for (key in finished) {
      harvest(as.integer(key))
    }
  }

  # Deliver in fold order, up to and including the first failure.
  last <- min(n, first_failure)
  out <- vector("list", n)
  for (i in seq_len(last)) {
    if (!is.null(failures[[i]])) {
      stop(failures[[i]])
    }
    out[[i]] <- .fb_isolation_deliver(records[[i]], adopt_seed = FALSE)
  }
  out
}
