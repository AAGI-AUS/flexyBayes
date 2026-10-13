# =============================================================================
# Engine fits in a separate R process.
#
# `options(flexyBayes.isolation = "process")` runs the engine fit behind
# flexybayes() (and so behind fb(), fb_inla() and fb_brms()) and the fit
# behind proceed() in a separate R process started through callr. Four
# claims are tested here. test-isolated-lapply.R covers loops of jobs.
#
#   1. A process killed by a signal inside a fit ends in a typed refusal that
#      names the engine, the entry point and the signal, and the session goes
#      on: the next fit, isolated or not, returns.
#   2. The child fits what the session would have fitted: the namespace this
#      session loaded (not whichever build is installed), the same options,
#      engine settings, environment, random-number state and once-per-session
#      notes. Printed output, messages and warnings come back in the order
#      they were raised, and a refusal raised in the child arrives as itself.
#   3. INLA does not return bit-identical numbers for two identical calls in
#      one session on this build: over ten in-process fits of the fixture
#      below, 36 of 45 pairs differed, by at most 2.18e-5 relative to the
#      scale of the quantity (measured 2026-09-19). The engine's numbers
#      are therefore held to bounds derived from the measured in-process
#      spread. What flexyBayes builds
#      around the engine is compared exactly, and the fields that differ only
#      because an object crossed a process boundary (its environment is a
#      copy) are named, not skipped.
#   4. brms with a seed returns identical draws, and with no seed the draw
#      comes from the session's random-number stream in both places.
#
# How a process is killed here. A watcher process started by this file looks
# only at the children of this R session, skips itself, and sends SIGKILL to
# one PID: the first child that is itself running a program (the engine), so
# the signal lands inside the fit, or, where no engine runs, the first child.
# It never selects a process by name.
# =============================================================================

# --- fixtures ------------------------------------------------------------- #

fbiso_data <- function(seed = 20260919L, n_group = 12L, per = 8L) {
  .fb_with_seed(seed, {
    n <- n_group * per
    d <- data.frame(
      g = factor(rep(sprintf("G%02d", seq_len(n_group)), each = per))
    )
    d$x <- stats::rnorm(n)
    d$z1 <- stats::rnorm(n)
    d$z2 <- stats::rnorm(n)
    u <- stats::rnorm(n_group, sd = 0.7)
    d$y <- 1 + 0.5 * d$x + 0.8 * d$z1 + u[as.integer(d$g)] +
      stats::rnorm(n, sd = 0.5)
    d
  })
}


# Run `expr` with the isolation option set, and put it back afterwards.
fbiso_with <- function(mode, expr) {
  old <- options(flexyBayes.isolation = mode)
  on.exit(options(old), add = TRUE)
  force(expr)
}


fbiso_skip <- function() {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("ps")
  if (!identical(.Platform$OS.type, "unix")) {
    skip("a signal is sent to one process ID on unix only")
  }
  invisible(NULL)
}


# The watcher. It returns the PID it signalled, or NA if it found none before
# its deadline. Only a child started after the watcher itself, and still
# running, is a candidate: an earlier block's finished process can still be
# listed as a child until it is reaped. `busy = TRUE` waits for a child that
# has a child of its own, which for an isolated fit is the engine running
# inside it.
fbiso_watcher <- function(busy = TRUE, timeout = 120) {
  callr::r_bg(
    function(target, busy, timeout) {
      me <- Sys.getpid()
      since <- ps::ps_create_time(ps::ps_handle(me))
      parent <- ps::ps_handle(target)
      deadline <- Sys.time() + timeout
      while (Sys.time() < deadline) {
        kids <- tryCatch(ps::ps_children(parent), error = function(e) list())
        for (k in kids) {
          pid <- ps::ps_pid(k)
          if (identical(as.integer(pid), as.integer(me))) {
            next
          }
          fresh <- tryCatch(
            ps::ps_create_time(k) >= since &&
              !identical(ps::ps_status(k), "zombie"),
            error = function(e) FALSE
          )
          if (!fresh) {
            next
          }
          running <- !busy || length(tryCatch(
            ps::ps_children(k),
            error = function(e) list()
          )) > 0L
          if (running) {
            tools::pskill(pid, tools::SIGKILL)
            return(as.integer(pid))
          }
        }
        Sys.sleep(0.02)
      }
      NA_integer_
    },
    args = list(target = Sys.getpid(), busy = busy, timeout = timeout),
    supervise = TRUE
  )
}


# What the watcher did: the PID it signalled, or NA. A watcher still running
# is stopped through its own handle.
fbiso_watcher_result <- function(watcher) {
  watcher$wait(5000)
  if (watcher$is_alive()) {
    watcher$kill()
    return(NA_integer_)
  }
  watcher$get_result()
}


# Every leaf at which two objects differ, as a path, with a label saying what
# kind of difference it is. A formula or a function whose text is identical
# and whose environment is not differs only in where it lives: "environment".
fbiso_diff <- function(a, b, path = "") {
  if (identical(a, b)) {
    return(character(0))
  }
  if (is.environment(a) || is.environment(b)) {
    return(paste0(path, " <environment>"))
  }
  if (
    (is.function(a) && is.function(b)) ||
      (inherits(a, "formula") && inherits(b, "formula"))
  ) {
    same_text <- identical(deparse(a), deparse(b))
    return(paste0(path, if (same_text) " <environment>" else " <value>"))
  }
  if (
    is.list(a) && is.list(b) && length(a) == length(b) &&
      identical(names(a), names(b))
  ) {
    nm <- names(a) %||% as.character(seq_along(a))
    nm[!nzchar(nm)] <- as.character(which(!nzchar(nm)))
    out <- unlist(lapply(seq_along(a), function(i) {
      fbiso_diff(a[[i]], b[[i]], paste0(path, "$", nm[[i]]))
    }))
    if (length(out) == 0L) {
      return(paste0(path, " <attributes>"))
    }
    return(out)
  }
  paste0(path, " <value>")
}


# The value at a path `fbiso_diff()` reported.
fbiso_get <- function(x, path) {
  keys <- strsplit(sub(" <[a-z]+>$", "", path), "$", fixed = TRUE)[[1L]]
  for (key in keys[nzchar(keys)]) {
    x <- x[[key]]
  }
  x
}


# The engine's numbers, and their scale-relative distance.
fbiso_inla_numbers <- function(fit) {
  list(
    fixed_mean = fit$inla$summary.fixed$mean,
    fixed_sd = fit$inla$summary.fixed$sd,
    hyper_mean = fit$inla$summary.hyperpar$mean,
    random_mean = fit$inla$summary.random$g$mean,
    fitted_mean = fit$inla$summary.fitted.values$mean,
    mlik = as.numeric(fit$inla$mlik)
  )
}

fbiso_distance <- function(a, b) {
  x <- fbiso_inla_numbers(a)
  y <- fbiso_inla_numbers(b)
  vapply(names(x), function(nm) {
    max(abs(x[[nm]] - y[[nm]])) / max(abs(c(x[[nm]], y[[nm]])))
  }, numeric(1))
}

# From ten in-process fits (45 pairs, largest 2.18e-5): about ten times
# the largest measured distance, fixed before any isolated fit was compared.
FBISO_INLA_BOUND <- 2e-4

# The numbers flexyBayes derives from INLA's and keeps under $extras (the
# summary tables, the variance components, the manifest's evidence). Their
# in-process spread over ten fits: largest 3.45e-4, on the approximation's
# own kld diagnostic, and at most 7.07e-5 on every other field. The bound is
# about ten times the largest.
FBISO_DERIVED_BOUND <- 3.5e-3


# --- the option ------------------------------------------------------------ #

test_that("isolation is off unless asked for, and the option takes three values", {
  expect_null(getOption("flexyBayes.isolation"))
  expect_identical(.fb_isolation_mode(), "none")
  fbiso_with("none", expect_identical(.fb_isolation_mode(), "none"))
  fbiso_with("process", expect_identical(.fb_isolation_mode(), "process"))
  fbiso_with("fork", expect_identical(.fb_isolation_mode(), "fork"))
  # "fork" changes fold workers only: a single fit stays in this session.
  fbiso_with("fork", expect_false(.fb_isolation_active("flexybayes()")))
  fbiso_with("thread", expect_error(.fb_isolation_mode(), "\"process\""))
  fbiso_with(TRUE, expect_error(.fb_isolation_mode(), "\"process\""))
  fbiso_with(c("none", "process"), expect_error(.fb_isolation_mode()))
})


test_that("asking for isolation without callr refuses by name, before any fit", {
  testthat::local_mocked_bindings(
    .fb_isolation_callr_available = function() FALSE
  )
  d <- fbiso_data()
  refusal <- fbiso_with("process", tryCatch(
    flexybayes(y ~ x, random = ~g, data = d, backend = "inla", verbose = FALSE),
    error = function(e) e
  ))
  expect_s3_class(refusal, "flexybayes_refusal")
  expect_identical(refusal$reason_code, "isolation_requires_callr")
  expect_match(conditionMessage(refusal), "callr")
  expect_match(conditionMessage(refusal), "flexyBayes.isolation")
})


# --- what the child is ----------------------------------------------------- #

test_that("the child runs this session's namespace, options, engine settings and state", {
  fbiso_skip()
  skip_if_not_installed("INLA")
  old_threads <- INLA::inla.getOption("num.threads")
  on.exit(INLA::inla.setOption(num.threads = old_threads), add = TRUE)
  INLA::inla.setOption(num.threads = "1:1")
  withr::local_envvar(OMP_NUM_THREADS = "3")
  withr::local_options(list(flexyBayes.silence_convergence_warning = TRUE))
  .reset_emit_state_for_test()
  on.exit(.reset_emit_state_for_test(), add = TRUE)

  set.seed(42L)
  expected_draw <- stats::rnorm(1L)
  expected_state <- get(".Random.seed", envir = globalenv())
  set.seed(42L)

  info <- fbiso_with("process", .fb_isolated_call(
    function() {
      ns <- asNamespace("flexyBayes")
      was <- ns$.emit_state_get("auto_fallback_note")
      ns$.emit_state_set("auto_fallback_note", TRUE)
      list(
        pid = Sys.getpid(),
        path = getNamespaceInfo(ns, "path"),
        isolation = getOption("flexyBayes.isolation"),
        option = getOption("flexyBayes.silence_convergence_warning"),
        inla_threads = INLA::inla.getOption("num.threads"),
        omp = Sys.getenv("OMP_NUM_THREADS"),
        note_was = was,
        draw = stats::rnorm(1L)
      )
    },
    engine = "none",
    entry_point = "a test"
  ))

  # Another process, running the namespace this session loaded: under
  # pkgload that is the source tree, not an installed build.
  expect_false(identical(info$pid, Sys.getpid()))
  expect_identical(info$path, getNamespaceInfo(asNamespace("flexyBayes"), "path"))
  # No second level of isolation inside the child.
  expect_identical(info$isolation, "none")
  expect_true(info$option)
  expect_identical(info$inla_threads, "1:1")
  expect_identical(info$omp, "3")
  # The session's random-number state went in, and the child's came back, so
  # a draw inside an isolated call is the draw an in-process call makes and
  # the session's stream moves on as it would have.
  expect_identical(info$draw, expected_draw)
  expect_identical(get(".Random.seed", envir = globalenv()), expected_state)
  # A once-per-session note is once per session, not once per process.
  expect_false(info$note_was)
  expect_true(.emit_state_get("auto_fallback_note"))
})


test_that("output, messages and warnings come back in order, and a refusal arrives as itself", {
  fbiso_skip()
  events <- character(0)
  printed <- utils::capture.output(value <- withCallingHandlers(
    fbiso_with("process", .fb_isolated_call(
      function() {
        cat("line one\n")
        message("a message")
        warning(structure(
          class = c("flexybayes_deprecated", "warning", "condition"),
          list(message = "a classed warning", call = NULL)
        ))
        cat("line two\n")
        7
      },
      engine = "none",
      entry_point = "a test"
    )),
    message = function(m) {
      events <<- c(events, paste0("message: ", trimws(conditionMessage(m))))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      events <<- c(
        events,
        paste0("warning: ", conditionMessage(w), " [", class(w)[[1L]], "]")
      )
      invokeRestart("muffleWarning")
    }
  ))
  expect_identical(value, 7)
  expect_identical(printed, c("line one", "line two"))
  expect_identical(
    events,
    c("message: a message", "warning: a classed warning [flexybayes_deprecated]")
  )

  refusal <- tryCatch(
    fbiso_with("process", .fb_isolated_call(
      function() {
        stop(asNamespace("flexyBayes")$.fb_refusal_condition(
          reason_code = "isolated_process_died",
          message = "a refusal inside the child"
        ))
      },
      engine = "none",
      entry_point = "a test"
    )),
    error = function(e) e
  )
  expect_s3_class(refusal, "flexybayes_refusal")
  expect_identical(refusal$reason_code, "isolated_process_died")
  expect_match(conditionMessage(refusal), "inside the child")
})


# --- a process that dies ---------------------------------------------------- #

test_that("a child killed by a signal inside an INLA fit refuses by name, and the session goes on", {
  fbiso_skip()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()
  withr::local_options(list(flexyBayes.silence_default_prior_note = TRUE))
  # Large enough that the engine runs long enough to be caught running.
  d <- fbiso_data(n_group = 40L, per = 25L)
  fit_once <- function() {
    flexybayes(y ~ x, random = ~g, data = d, backend = "inla", verbose = FALSE)
  }

  watcher <- fbiso_watcher(busy = TRUE)
  refusal <- fbiso_with("process", tryCatch(fit_once(), error = function(e) e))
  killed <- fbiso_watcher_result(watcher)

  expect_false(is.na(killed))
  expect_s3_class(refusal, "flexybayes_refusal")
  expect_identical(refusal$reason_code, "isolated_process_died")
  expect_identical(refusal$engine, "inla")
  expect_identical(refusal$entry_point, "flexybayes()")
  expect_identical(refusal$signal, 9L)
  expect_identical(refusal$exit_status, -9L)
  expect_identical(refusal$pid, killed)
  expect_match(conditionMessage(refusal), "signal 9")
  expect_match(conditionMessage(refusal), "INLA")
  expect_match(conditionMessage(refusal), "flexybayes()", fixed = TRUE)

  # The session goes on: the next isolated fit returns, and so does the next
  # fit in this process.
  again <- fbiso_with("process", fit_once())
  expect_s3_class(again, "flexybayes_inla")
  expect_true(isTRUE(again$num_check$pass))
  here <- fit_once()
  expect_s3_class(here, "flexybayes_inla")
})


# --- the same answer ------------------------------------------------------- #

test_that("an isolated INLA fit matches in-process fits, exactly outside the engine and within its spread inside it", {
  fbiso_skip()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()
  old_threads <- INLA::inla.getOption("num.threads")
  on.exit(INLA::inla.setOption(num.threads = old_threads), add = TRUE)
  INLA::inla.setOption(num.threads = "1:1")
  withr::local_options(list(flexyBayes.silence_default_prior_note = TRUE))
  d <- fbiso_data()
  fit_once <- function() {
    flexybayes(y ~ x, random = ~g, data = d, backend = "inla", verbose = FALSE)
  }
  here_a <- fit_once()
  here_b <- fit_once()
  isolated <- fbiso_with("process", fit_once())
  expect_s3_class(isolated, "flexybayes_inla")

  # The one field isolation adds records where the fit ran. It is set aside
  # before the comparison below.
  expect_identical(isolated$extras$isolation$mode, "process")
  expect_false(identical(isolated$extras$isolation$pid, Sys.getpid()))
  expect_null(here_a$extras$isolation)
  isolated$extras$isolation <- NULL

  # Outside the engine's own output ($inla), a field of the isolated fit may
  # differ from an in-process fit's in three ways only, each named: a formula
  # whose text is identical and whose environment is a copy (two in-process
  # calls from a function differ the same way), the run time, or a number
  # derived from the engine's numbers, held to the same bound.
  outside <- function(x) x[!startsWith(x, "$inla$")]
  for (here in list(here_a, here_b)) {
    cross <- outside(fbiso_diff(here, isolated))
    boundary <- cross[endsWith(cross, "<environment>")]
    expect_true(all(
      sub(" <environment>$", "", boundary) %in%
        c("$extras$call_info$fixed", "$extras$call_info$random", "$extras$formula")
    ))
    derived <- setdiff(cross, c(boundary, "$extras$run_time <value>"))
    for (path in derived) {
      a <- fbiso_get(here, path)
      b <- fbiso_get(isolated, path)
      expect_true(is.numeric(a) && is.numeric(b), info = path)
      expect_lte(
        max(abs(a - b)) / max(abs(c(a, b))),
        FBISO_DERIVED_BOUND,
        label = path
      )
    }
  }

  # Inside the engine's output, the numbers sit within the bound derived from
  # the engine's own in-process repeat spread.
  expect_true(all(fbiso_distance(here_a, isolated) <= FBISO_INLA_BOUND))
  expect_true(all(fbiso_distance(here_b, isolated) <= FBISO_INLA_BOUND))
  # Everything flexyBayes built around the engine is exact.
  for (part in c("fb", "data", "level_labels", "exactness")) {
    expect_identical(isolated[[part]], here_a[[part]], info = part)
  }
})


test_that("an isolated brms fit returns the in-process draws, with a seed and from the session's stream", {
  fbiso_skip()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  skip_if_not_long_tier()
  withr::local_options(list(
    brms.backend = "cmdstanr",
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE
  ))
  d <- fbiso_data()
  fit_once <- function(seed) {
    suppressWarnings(flexybayes(
      y ~ x + (1 | g),
      data = d,
      backend = "brms",
      chains = 1L,
      n_samples = 150L,
      warmup = 150L,
      seed = seed,
      verbose = FALSE,
      mcmc_verbose = FALSE
    ))
  }
  draws <- function(fit) posterior::as_draws_matrix(fit$brms)

  here <- fit_once(11L)
  isolated <- fbiso_with("process", fit_once(11L))
  expect_s3_class(isolated, "flexybayes_brms")
  expect_identical(isolated$extras$isolation$mode, "process")
  expect_false(identical(isolated$extras$isolation$pid, Sys.getpid()))
  expect_identical(draws(isolated), draws(here))

  # No seed: Stan's seed is drawn from the session's stream, in the session
  # or in the child, and the stream moves on the same way in both.
  set.seed(8L)
  here_stream <- fit_once(NULL)
  state_here <- get(".Random.seed", envir = globalenv())
  set.seed(8L)
  isolated_stream <- fbiso_with("process", fit_once(NULL))
  state_isolated <- get(".Random.seed", envir = globalenv())
  expect_identical(draws(isolated_stream), draws(here_stream))
  expect_identical(state_isolated, state_here)

  # proceed() fits a reviewed brms program, and the option reaches it too.
  # A review object keeps its fit once proceeded, so each arm has its own.
  review <- function() {
    suppressWarnings(flexybayes(
      y ~ x + (1 | g),
      data = d,
      backend = "brms",
      chains = 1L,
      n_samples = 150L,
      warmup = 150L,
      seed = 11L,
      verbose = FALSE,
      mcmc_verbose = FALSE,
      review_code = TRUE
    ))
  }
  proceeded_here <- suppressWarnings(proceed(review()))
  proceeded_isolated <- fbiso_with(
    "process",
    suppressWarnings(proceed(review()))
  )
  expect_null(proceeded_here$extras$isolation)
  expect_identical(proceeded_isolated$extras$isolation$mode, "process")
  expect_identical(draws(proceeded_isolated), draws(proceeded_here))
})


