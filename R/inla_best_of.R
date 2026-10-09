# inla_best_of.R -- every INLA fit is run k times and the fit with the highest
# marginal likelihood is kept.
#
# Identical one-thread INLA calls in fresh R processes can reach different
# posterior modes, while repeated fits within one R session agree
# (inst/KNOWN_ISSUES.md, item 26). The k fits are therefore the fit in this
# session and k - 1 fits in fresh R processes. The fresh-process fits start
# before the session fit and run beside it, up to option
# `flexyBayes.inla_best_of_workers` fits at once, fewer when free memory
# does not hold that many. The spread of the k fits is recorded at
# `$extras$reproducibility$best_of`.

# Relative difference of a standard deviation across the k fits above which
# print() and summary() say so.
.FB_BEST_OF_SD_NOTE <- 0.05

# The memory of one INLA fit, in bytes: .FB_BEST_OF_BASE_BYTES plus
# .FB_BEST_OF_ESTIMATE_FACTOR times the design-memory preflight's INLA
# estimate. Set from two peak resident sets measured with /usr/bin/time -l:
# 0.96 GB of R plus up to 0.5 GB of the INLA program on the six-trial barley
# subset (preflight estimate 8 MB), and 8.3 GB on the 48-trial series
# (112 MB). The smallest per-row fit in inst/validation/benchmark_scaling.csv
# peaked at 0.99 GB. The factor is set so both measured peaks fall below the
# estimate (about 1.7 GB and 10 GB): a low estimate would start fits into swap.
.FB_BEST_OF_BASE_BYTES <- 1024^3
.FB_BEST_OF_ESTIMATE_FACTOR <- 80

#' The number of INLA fits each model is run with
#'
#' Called before the first fit.
#'
#' @return A whole number, one or more: option `flexyBayes.inla_best_of`,
#'   3 when unset.
#' @noRd
#' @keywords internal
.fb_inla_best_of_k <- function() {
  value <- getOption("flexyBayes.inla_best_of", 3L)
  ok <- is.numeric(value) &&
    length(value) == 1L &&
    is.finite(value) &&
    value >= 1 &&
    value == round(value)
  if (!ok) {
    stop(.fb_refusal_condition(
      reason_code = "inla_best_of_invalid",
      message = paste0(
        "option flexyBayes.inla_best_of must be a whole number, one or more ",
        "(the number of INLA fits each model is run with). It is: ",
        .fb_prior_deparse(value),
        "."
      )
    ))
  }
  k <- as.integer(value)
  if (k > 1L && !.fb_isolation_callr_available()) {
    stop(.fb_refusal_condition(
      reason_code = "isolation_requires_callr",
      message = paste0(
        "flexybayes(): option flexyBayes.inla_best_of = ",
        k,
        " runs ",
        k - 1L,
        " of the INLA fits in separate R processes through callr, and ",
        "callr is not installed. Install it with install.packages(\"callr\"), ",
        "or set options(flexyBayes.inla_best_of = 1) to fit once."
      ),
      entry_point = "flexybayes()"
    ))
  }
  k
}


#' Start one INLA fit in a fresh R process
#'
#' Tests replace this and the three functions below it.
#'
#' @param args The arguments of `INLA::inla()`.
#' @return The handle of the running process.
#' @noRd
#' @keywords internal
.fb_inla_isolated_start <- function(args) {
  fn_path <- tempfile("fb_isolation_fn_", fileext = ".rds")
  saveRDS(.fb_inla_child_fit, fn_path, compress = FALSE)
  handle <- .fb_isolation_start(fn_path, list(args = args))
  handle$files <- c(handle$files, fn_path)
  handle
}

#' Whether a fit started by .fb_inla_isolated_start() is still running
#' @noRd
.fb_inla_isolated_running <- function(handle) {
  handle$proc$is_alive()
}

#' Wait for a fit and return what its process sent back, undelivered
#' @noRd
.fb_inla_isolated_collect <- function(handle) {
  .fb_isolation_harvest(handle, "inla", "flexybayes()")
}

#' End a fit's process and what it started
#' @noRd
.fb_inla_isolated_stop <- function(handle) {
  .fb_isolation_stop(handle)
}

# Defined at the top level so that only `args`, and no frame of the caller,
# is sent to the separate process.
.fb_inla_child_fit <- function(args) {
  .inla_call(args)
}


#' The log marginal likelihood of an INLA fit, `-Inf` when it has none
#' @noRd
.fb_best_of_mlik <- function(fit) {
  m <- fit$mlik
  value <- if (is.null(m) || !length(m)) NA_real_ else as.numeric(m[1L, 1L])
  if (is.finite(value)) value else -Inf
}


#' The posterior mean of each standard deviation of an INLA fit
#'
#' @return A named numeric vector, one entry per `Precision for` row of the
#'   hyperparameter table, `NA` where the fit kept no marginal.
#' @noRd
.fb_best_of_sds <- function(fit) {
  hp <- fit$summary.hyperpar
  nm <- if (is.null(hp)) character(0) else rownames(hp)
  nm <- nm[startsWith(nm, "Precision for ")]
  marginals <- fit$marginals.hyperpar %||% list()
  out <- vapply(
    nm,
    function(x) {
      m <- marginals[[x]]
      if (is.null(m)) {
        return(NA_real_)
      }
      INLA::inla.emarginal(function(p) 1 / sqrt(p), m)
    },
    numeric(1)
  )
  names(out) <- sub("^Precision for ", "", nm)
  out
}


#' The number of the k fits that run at once
#'
#' Option `flexyBayes.inla_best_of_workers` (default `k`, at most `k`),
#' lowered to the number of fits the free memory holds. A lowered count is
#' said in one message.
#'
#' @param k The number of fits.
#' @param fb,data The model and data the fits are made from.
#' @return A whole number from 1 to `k`; 1 runs the fits one after another.
#' @noRd
#' @keywords internal
.fb_inla_best_of_workers <- function(k, fb, data) {
  value <- getOption("flexyBayes.inla_best_of_workers", k)
  ok <- is.numeric(value) &&
    length(value) == 1L &&
    is.finite(value) &&
    value >= 1 &&
    value == round(value)
  if (!ok) {
    stop(
      "option flexyBayes.inla_best_of_workers must be a whole number, one ",
      "or more (the number of INLA fits run at once). It is: ",
      .fb_prior_deparse(value),
      ".",
      call. = FALSE
    )
  }
  workers <- min(as.integer(value), k)
  if (workers == 1L) {
    return(1L)
  }
  fit_bytes <- .fb_inla_fit_bytes(fb, data)
  free_bytes <- .fb_free_memory_bytes()
  gb <- function(x) paste0(format(round(x / 1024^3, 1), nsmall = 1), " GB")
  why <- if (is.na(free_bytes)) {
    "the free memory could not be read"
  } else if (is.na(fit_bytes)) {
    paste0(
      "the memory of one fit could not be estimated (",
      attr(fit_bytes, "reason") %||% "no estimate",
      ")"
    )
  } else {
    paste0(
      gb(free_bytes),
      " of memory is free and one fit needs about ",
      gb(fit_bytes)
    )
  }
  fits <- if (is.na(free_bytes) || is.na(fit_bytes)) {
    1L
  } else {
    max(1L, as.integer(floor(free_bytes / fit_bytes)))
  }
  if (fits >= workers) {
    return(workers)
  }
  message(
    "flexyBayes: ",
    if (fits == 1L) {
      paste0("the ", k, " INLA fits run one after another: ")
    } else {
      paste0(fits, " of the ", k, " INLA fits run at once: ")
    },
    why,
    " (option flexyBayes.inla_best_of_workers)."
  )
  fits
}


#' The memory one INLA fit of this model needs, in bytes
#'
#' @return A number, or `NA` with attribute `reason` when the design-memory
#'   preflight cannot size the model.
#' @noRd
#' @keywords internal
.fb_inla_fit_bytes <- function(fb, data) {
  estimate <- tryCatch(
    as.numeric(.fb_preflight(fb, .fb_dataset(data))$memory_estimate),
    error = function(e) structure(NA_real_, reason = conditionMessage(e))
  )
  if (length(estimate) != 1L || !is.finite(estimate)) {
    return(structure(
      NA_real_,
      reason = attr(estimate, "reason") %||% "no estimate"
    ))
  }
  .FB_BEST_OF_BASE_BYTES + .FB_BEST_OF_ESTIMATE_FACTOR * estimate
}


#' Memory a new process can be given now, in bytes
#'
#' Free and inactive pages from `vm_stat` on macOS, `MemAvailable` from
#' `/proc/meminfo` on Linux, `NA` elsewhere.
#' @noRd
#' @keywords internal
.fb_free_memory_bytes <- function() {
  sysname <- Sys.info()[["sysname"]]
  if (identical(sysname, "Darwin")) {
    if (!nzchar(Sys.which("vm_stat")[[1L]])) {
      return(NA_real_)
    }
    out <- tryCatch(
      suppressWarnings(system2("vm_stat", stdout = TRUE, stderr = FALSE)),
      error = function(e) character(0)
    )
    page <- regmatches(out, regexpr("page size of [0-9]+", out))
    pages <- function(label) {
      line <- grep(paste0("^", label, ":"), out, value = TRUE)
      if (length(line) != 1L) {
        return(NA_real_)
      }
      as.numeric(gsub("[^0-9]", "", sub("^[^:]*:", "", line)))
    }
    if (length(page) != 1L) {
      return(NA_real_)
    }
    page <- as.numeric(sub("page size of ", "", page))
    return(page * (pages("Pages free") + pages("Pages inactive")))
  }
  if (identical(sysname, "Linux")) {
    line <- tryCatch(
      grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE),
      error = function(e) character(0)
    )
    if (length(line) != 1L) {
      return(NA_real_)
    }
    return(as.numeric(gsub("[^0-9]", "", line)) * 1024)
  }
  NA_real_
}


#' Start the fresh-process fits that run beside the session fit
#'
#' Called just before the session fit, with the arguments it is made with.
#' The caller ends the pool with .fb_inla_best_of_stop() on exit, so a
#' failed session fit leaves no process behind.
#'
#' @param args The arguments of `INLA::inla()`.
#' @param k The number of fits, from `.fb_inla_best_of_k()`.
#' @param fb,data The model and data, for the memory guard.
#' @return `NULL` when `k` is 1, otherwise the pool .fb_inla_best_of() reads.
#' @noRd
#' @keywords internal
.fb_inla_best_of_launch <- function(args, k, fb, data) {
  if (k == 1L) {
    return(NULL)
  }
  pool <- .fb_best_of_pool(args, k, .fb_inla_best_of_workers(k, fb, data))
  # The session fit is one of the fits running at once.
  .fb_best_of_fill(pool, pool$workers - 1L)
  pool$session_start <- proc.time()[["elapsed"]]
  pool
}

.fb_best_of_pool <- function(args, k, workers) {
  pool <- new.env(parent = emptyenv())
  pool$args <- args
  pool$k <- k
  pool$workers <- workers
  pool$handles <- vector("list", k)
  pool$records <- vector("list", k)
  pool$failures <- vector("list", k)
  pool$started <- rep(NA_real_, k)
  pool$wall_s <- rep(NA_real_, k)
  pool$next_j <- 2L
  pool$first_failure <- Inf
  pool$session_start <- NA_real_
  pool
}

# Start fits until `limit` run, none past the first failure.
.fb_best_of_fill <- function(pool, limit) {
  running <- sum(!vapply(pool$handles, is.null, logical(1)))
  while (
    pool$next_j <= pool$k && pool$next_j < pool$first_failure &&
      running < limit
  ) {
    j <- pool$next_j
    pool$started[[j]] <- proc.time()[["elapsed"]]
    pool$handles[[j]] <- .fb_inla_isolated_start(pool$args)
    pool$next_j <- j + 1L
    running <- running + 1L
  }
  invisible(pool)
}

# Collect fit j: its record, or the condition it failed with.
.fb_best_of_collect <- function(pool, j) {
  handle <- pool$handles[[j]]
  pool$handles[j] <- list(NULL)
  record <- tryCatch(
    .fb_inla_isolated_collect(handle),
    error = function(e) e
  )
  pool$wall_s[[j]] <- proc.time()[["elapsed"]] - pool$started[[j]]
  if (inherits(record, "error")) {
    pool$failures[[j]] <- record
  } else {
    pool$records[[j]] <- record
  }
  if (inherits(record, "error") || !isTRUE(record$ok)) {
    pool$first_failure <- min(pool$first_failure, j)
  }
  invisible(pool)
}

#' End every fit of the pool still running
#' @noRd
.fb_inla_best_of_stop <- function(pool) {
  if (is.null(pool)) {
    return(invisible(NULL))
  }
  for (j in seq_len(pool$k)) {
    if (!is.null(pool$handles[[j]])) {
      .fb_inla_isolated_stop(pool$handles[[j]])
      pool$handles[j] <- list(NULL)
    }
  }
  invisible(NULL)
}

# The value of fit j, its messages and warnings re-signalled; a failure
# stops with the fit's number.
.fb_best_of_value <- function(pool, j, adopt_seed) {
  tryCatch(
    {
      if (!is.null(pool$failures[[j]])) {
        stop(pool$failures[[j]])
      }
      .fb_isolation_deliver(pool$records[[j]], adopt_seed = adopt_seed)
    },
    error = function(e) {
      e$message <- paste0(
        "INLA fit ",
        j,
        " of ",
        pool$k,
        ", run in a separate R process, failed: ",
        conditionMessage(e)
      )
      stop(e)
    }
  )
}


#' Collect the remaining k - 1 fits and keep the best of the k
#'
#' @param fit The fit made in this session.
#' @param args The arguments it was made with.
#' @param k The number of fits, from `.fb_inla_best_of_k()`, read before
#'   the first fit.
#' @param pool From `.fb_inla_best_of_launch()`; `NULL` runs the k - 1 fits
#'   one after another from here.
#' @param fa_term The model's fa(f, k):g term, or `NULL`. With one, the
#'   record also holds `fa_gate`, the agreement of the k fits on the
#'   genetic covariance (`.fb_fa_gate()`).
#' @return A list: `fit`, the fit with the highest marginal likelihood (the
#'   first of equals), and `record`, `NULL` when `k` is 1, otherwise
#'   `k`, `mlik`, `mode_status`, `chosen`, `mlik_range`, `sd` (one row per
#'   fit), `sd_rel_diff`, the largest difference of each standard
#'   deviation across the fits relative to its smallest value, `wall_s`,
#'   the seconds each fit took, `workers`, the number of fits run at once,
#'   and `parallel`, whether that number was above 1.
#' @noRd
#' @keywords internal
.fb_inla_best_of <- function(fit, args, k, pool = NULL, fa_term = NULL) {
  if (k == 1L) {
    return(list(fit = fit, record = NULL))
  }
  if (is.null(pool)) {
    pool <- .fb_best_of_pool(args, k, 1L)
  }
  on.exit(.fb_inla_best_of_stop(pool), add = TRUE)
  pool$wall_s[[1L]] <- proc.time()[["elapsed"]] - pool$session_start
  sequential <- pool$workers == 1L
  fits <- vector("list", k)
  fits[[1L]] <- fit
  # Fits up to the first failure are waited for; one run at a time is
  # delivered as it finishes, so each process starts from the random-number
  # state the one before it left.
  repeat {
    .fb_best_of_fill(pool, pool$workers)
    open <- which(!vapply(pool$handles, is.null, logical(1)))
    open <- open[open < pool$first_failure]
    if (!length(open)) {
      break
    }
    done <- open[!vapply(
      pool$handles[open],
      .fb_inla_isolated_running,
      logical(1)
    )]
    if (!length(done)) {
      Sys.sleep(0.05)
      next
    }
    for (j in done) {
      .fb_best_of_collect(pool, j)
      if (sequential) {
        fits[[j]] <- .fb_best_of_value(pool, j, adopt_seed = TRUE)
      }
    }
  }
  if (!sequential) {
    for (j in seq.int(2L, min(k, pool$first_failure))) {
      fits[[j]] <- .fb_best_of_value(pool, j, adopt_seed = FALSE)
    }
  }
  mlik <- vapply(fits, .fb_best_of_mlik, numeric(1))
  status <- vapply(
    fits,
    function(f) as.integer(f$mode$mode.status %||% NA_integer_),
    integer(1)
  )
  sds <- lapply(fits, .fb_best_of_sds)
  chosen <- which.max(mlik)
  sd_mat <- do.call(rbind, sds)
  sd_rel_diff <- if (is.null(sd_mat) || !ncol(sd_mat)) {
    numeric(0)
  } else {
    apply(sd_mat, 2L, function(s) (max(s) - min(s)) / min(s))
  }
  mlik[!is.finite(mlik)] <- NA_real_
  list(
    fit = fits[[chosen]],
    record = list(
      k = k,
      mlik = mlik,
      mode_status = status,
      chosen = chosen,
      mlik_range = if (all(is.na(mlik))) {
        NA_real_
      } else {
        diff(range(mlik, na.rm = TRUE))
      },
      sd = sd_mat,
      sd_rel_diff = sd_rel_diff,
      wall_s = pool$wall_s,
      workers = pool$workers,
      parallel = pool$workers > 1L,
      fa_gate = if (is.null(fa_term)) NULL else .fb_fa_gate(fits, fa_term)
    )
  )
}


#' Move the recorded marginal likelihoods to the original response scale
#'
#' @param record A `.fb_inla_best_of()` record, or `NULL`.
#' @param shift What the back-transform subtracts from the fit's marginal
#'   likelihood.
#' @noRd
.fb_best_of_shift <- function(record, shift) {
  if (!is.null(record)) {
    record$mlik <- record$mlik - shift
  }
  record
}


#' The line print() and summary() show when the k fits disagree
#'
#' @param record A `.fb_inla_best_of()` record.
#' @return A string, or `NULL` when no standard deviation differs by more
#'   than 5% across the fits.
#' @noRd
.fb_best_of_line <- function(record) {
  d <- record$sd_rel_diff
  if (!length(d) || !any(d > .FB_BEST_OF_SD_NOTE, na.rm = TRUE)) {
    return(NULL)
  }
  worst <- which.max(d)
  paste0(
    "fit ",
    record$chosen,
    " of ",
    record$k,
    " kept (highest marginal likelihood); the SD for ",
    names(d)[[worst]],
    " differed by ",
    format(round(100 * d[[worst]], 1), nsmall = 1),
    "% across the ",
    record$k,
    " fits"
  )
}


# Thrown by .fb_inla_best_of_k() and triggered in
# tests/testthat/test-inla-best-of.R.
.populate_refusal_registry_inla_best_of <- function() {
  .register_refusal(
    reason_code = "inla_best_of_invalid",
    description = paste0(
      "Option flexyBayes.inla_best_of, the number of INLA fits each model ",
      "is run with, is not a whole number of one or more. The request is ",
      "refused before anything is fitted."
    ),
    message_template = paste0(
      "option flexyBayes.inla_best_of must be a whole number, one or more. ",
      "It is: %s."
    ),
    registered_in_adr = "PLAN_0.17.0_repeatability_fa_speed_09Oct2026 (R1)",
    plan_field = NA_character_,
    since_version = "0.17.0"
  )
}
