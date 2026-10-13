# reproducibility_record.R, the versions, thread settings and seed a fit ran
# with.
#
# Every fit records, at `fit$extras$reproducibility`, what a second run has to
# match before its numbers can be expected to agree with the first. Matching it
# is not enough on INLA: a refit in a fresh R process can reach another mode
# (inst/KNOWN_ISSUES.md).
#
#   r_version            R.version.string
#   flexybayes_version   this package's version
#   engine               "inla" or "brms"
#   backend              brms only: the Stan interface that ran, "cmdstanr" or
#                        "rstan"
#   engine_versions      INLA, or brms with cmdstanr and CmdStan, or with rstan
#   threads              INLA: num_threads. brms: the chains run, the cores
#                        given to brms::brm(), the threads per chain. Both:
#                        OMP_NUM_THREADS as read when the fit was built, NA
#                        when unset.
#   seed                 supplied: the seed the emitter received, NA for none.
#                        as_run: the seed Stan ran with, which the backend
#                        drew when none was supplied, NA on INLA, which takes
#                        no seed.
#   best_of              INLA, when the fit was run more than once
#                        (option flexyBayes.inla_best_of): each fit's log
#                        marginal likelihood and mode status, the fit kept,
#                        the range of the marginal likelihoods, and each
#                        standard deviation per fit with its largest
#                        difference across the fits, each fit's seconds and
#                        whether the fits ran at once (R/inla_best_of.R).
#
# A value is read from the engine's own return wherever the engine reports it:
# INLA keeps the num.threads it ran with at `$.args$num.threads`, and CmdStan
# and rstan keep each chain's seed in the stanfit's `stan_args`. The option a
# setting defaults to is not read in its place, because an emitter's `...`
# can give the engine a different value. brms does not keep `cores` on the
# fit, so that one is the value brm() was given: its `cores` argument, or the
# `mc.cores` option that argument defaults to.

# --- building the record --------------------------------------------------- #

.fb_repro_version <- function(package) {
  as.character(utils::packageVersion(package))
}

.fb_repro_omp <- function() {
  Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
}

.fb_repro_seed_supplied <- function(seed) {
  if (is.null(seed)) NA_integer_ else as.integer(seed)
}

#' The reproducibility record of an INLA fit
#'
#' @param inla_fit The object `INLA::inla()` returned.
#' @param seed The `seed` the emitter received, or `NULL`.
#' @param best_of The `.fb_inla_best_of()` record, or `NULL` for one fit.
#' @returns The record, a named list (see the file header).
#' @noRd
#' @keywords internal
.fb_repro_record_inla <- function(inla_fit, seed = NULL, best_of = NULL) {
  num_threads <- inla_fit[[".args"]][["num.threads"]]
  rec <- list(
    r_version = R.version.string,
    flexybayes_version = .fb_repro_version("flexyBayes"),
    engine = "inla",
    engine_versions = c(INLA = .fb_repro_version("INLA")),
    threads = list(
      num_threads = if (is.null(num_threads)) {
        NA_character_
      } else {
        as.character(num_threads)
      },
      omp_num_threads = .fb_repro_omp()
    ),
    seed = list(
      supplied = .fb_repro_seed_supplied(seed),
      as_run = NA_integer_
    )
  )
  rec$best_of <- best_of
  rec
}

#' The reproducibility record of a brms fit
#'
#' @param brmsfit The object `brms::brm()` returned.
#' @param seed The `seed` the emitter received, or `NULL`.
#' @param cores The `cores` passed to `brms::brm()` through the emitter's
#'   `...`, or `NULL` when none was.
#' @returns The record, a named list (see the file header).
#' @noRd
#' @keywords internal
.fb_repro_record_brms <- function(brmsfit, seed = NULL, cores = NULL) {
  stanfit <- brmsfit[["fit"]]
  stan_args <- if (inherits(stanfit, "stanfit")) stanfit@stan_args else list()

  backend <- brmsfit[["backend"]] %||% NA_character_
  versions <- c(brms = .fb_repro_version("brms"))
  if (identical(backend, "cmdstanr")) {
    # The CmdStan version is the one the executable wrote into its output,
    # so it names the build that ran. cmdstanr::cmdstan_version() names the
    # installation configured now, and cmdstanr is not a declared dependency.
    written <- (stan_args[[1L]] %||% list())[c(
      "stan_version_major",
      "stan_version_minor",
      "stan_version_patch"
    )]
    versions <- c(
      versions,
      cmdstanr = .fb_repro_version("cmdstanr"),
      CmdStan = if (all(lengths(written) == 1L)) {
        paste(unlist(written), collapse = ".")
      } else {
        NA_character_
      }
    )
  } else if (identical(backend, "rstan")) {
    versions <- c(versions, rstan = .fb_repro_version("rstan"))
  }

  seeds <- unique(vapply(
    stan_args,
    function(a) as.character(a[["seed"]] %||% NA_character_),
    character(1L)
  ))
  # brms threads within a chain only when asked to. Without that a chain runs
  # on one thread, which is the `num_threads` CmdStan writes.
  per_chain <- brmsfit[["threads"]][["threads"]]

  list(
    r_version = R.version.string,
    flexybayes_version = .fb_repro_version("flexyBayes"),
    engine = "brms",
    backend = backend,
    engine_versions = versions,
    threads = list(
      chains = length(stan_args),
      cores = as.integer(cores %||% getOption("mc.cores", 1L)),
      threads_per_chain = if (is.null(per_chain)) 1L else as.integer(per_chain),
      omp_num_threads = .fb_repro_omp()
    ),
    seed = list(
      supplied = .fb_repro_seed_supplied(seed),
      as_run = if (length(seeds)) as.integer(seeds) else NA_integer_
    )
  )
}


# --- reading it ------------------------------------------------------------ #

#' The reproducibility record of a fit
#'
#' @param fit A flexyBayes fit.
#' @returns The record. Stops, naming each missing field, when the fit carries
#'   no record or a record short of any field its engine requires.
#' @noRd
#' @keywords internal
.fb_reproducibility <- function(fit) {
  rec <- if (is.list(fit) && is.list(fit[["extras"]])) {
    fit[["extras"]][["reproducibility"]]
  } else {
    NULL
  }
  if (is.null(rec)) {
    stop(
      "flexyBayes: this fit carries no reproducibility record at ",
      "`$extras$reproducibility`.",
      call. = FALSE
    )
  }
  lacking <- .fb_repro_missing(rec)
  if (length(lacking)) {
    stop(
      "flexyBayes: the reproducibility record of this fit lacks ",
      paste0("`", lacking, "`", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  rec
}

# The fields a record lacks, as paths into it. Which fields are required
# depends on the engine, and on a brms fit on the backend.
.fb_repro_missing <- function(rec) {
  engine <- rec[["engine"]]
  top <- c(
    "r_version",
    "flexybayes_version",
    "engine",
    "engine_versions",
    "threads",
    "seed"
  )
  if (identical(engine, "brms")) {
    top <- c(top, "backend")
  }
  lacking <- setdiff(top, names(rec))
  if (!identical(engine, "inla") && !identical(engine, "brms")) {
    return(union(lacking, "engine"))
  }

  inner <- if (identical(engine, "inla")) {
    list(
      engine_versions = "INLA",
      threads = c("num_threads", "omp_num_threads")
    )
  } else {
    list(
      engine_versions = c(
        "brms",
        switch(
          rec[["backend"]] %||% "",
          cmdstanr = c("cmdstanr", "CmdStan"),
          rstan = "rstan",
          character(0)
        )
      ),
      threads = c("chains", "cores", "threads_per_chain", "omp_num_threads")
    )
  }
  inner$seed <- c("supplied", "as_run")

  for (slot in names(inner)) {
    absent <- setdiff(inner[[slot]], names(rec[[slot]]))
    if (length(absent)) {
      lacking <- c(lacking, paste0(slot, "$", absent))
    }
  }
  lacking
}

# .fb_repro_lines(): the record as the lines print() and summary() show,
# named by their label. A fit without a record says so in one line.
#
# @noRd
# @keywords internal
.fb_repro_lines <- function(fit) {
  if (is.null(fit[["extras"]][["reproducibility"]])) {
    return(c(Versions = "not recorded on this fit"))
  }
  rec <- .fb_reproducibility(fit)

  engine_versions <- rec$engine_versions
  engine_versions[is.na(engine_versions)] <- "not recorded"
  versions <- paste(
    c(
      rec$r_version,
      paste("flexyBayes", rec$flexybayes_version),
      paste(names(engine_versions), engine_versions)
    ),
    collapse = "; "
  )
  omp <- rec$threads$omp_num_threads
  omp <- paste("OMP_NUM_THREADS", if (is.na(omp)) "unset" else omp)
  supplied <- rec$seed$supplied

  if (identical(rec$engine, "inla")) {
    threads <- paste0(
      "INLA num.threads ",
      .fb_repro_value(rec$threads$num_threads),
      "; ",
      omp
    )
    seed <- paste0(
      if (is.na(supplied)) {
        "none supplied"
      } else {
        paste(supplied, "supplied, not used")
      },
      "; INLA takes no seed"
    )
  } else {
    threads <- paste0(
      rec$threads$chains,
      " chain(s), ",
      rec$threads$cores,
      " core(s), ",
      rec$threads$threads_per_chain,
      " thread(s) per chain; ",
      omp
    )
    seed <- paste0(
      "Stan ran with ",
      .fb_repro_value(rec$seed$as_run),
      if (is.na(supplied)) " (none supplied)" else " (supplied)"
    )
  }
  c(
    Versions = versions,
    Threads = threads,
    Seed = seed,
    "Best of" = .fb_best_of_line(rec$best_of)
  )
}

.fb_repro_value <- function(x) {
  x[is.na(x)] <- "not recorded"
  if (length(x) == 0L) "not recorded" else paste(x, collapse = ", ")
}
