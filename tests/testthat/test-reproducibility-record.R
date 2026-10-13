# The reproducibility record.
#
# Every fit an entry point returns records, at `fit$extras$reproducibility`,
# the R and flexyBayes versions, the engine versions, the thread settings and
# the seed it ran with, and print() and summary() show them. The expected
# values are read here from R, the engines and the environment, not from the
# record: `R.version.string`, `utils::packageVersion()`,
# CmdStan's own version, the `num.threads` INLA actually received and the
# `num_threads` CmdStan writes. A seed Stan drew for itself is checked by
# refitting under the recorded value and comparing the draws.

.vr_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_sampler_arg_note = TRUE,
    .local_envir = parent.frame()
  )
}

.vr_version <- function(package) {
  as.character(utils::packageVersion(package))
}

.vr_omp <- function() {
  Sys.getenv("OMP_NUM_THREADS", unset = NA_character_)
}

# The CmdStan installation cmdstanr reports. cmdstanr is not a declared
# dependency, so it is reached by name at run time, behind
# skip_if_not_installed("cmdstanr") in every block that calls this.
.vr_cmdstan_version <- function() {
  as.character(getExportedValue("cmdstanr", "cmdstan_version")())
}

# The Gaussian and count design of test-prior-support-audit.R.
.vr_inla_data <- function(seed = 20260916L, n = 80L) {
  set.seed(seed)
  d <- data.frame(
    g = factor(rep(sprintf("g%02d", 1:8), length.out = n)),
    h = factor(rep(c("A", "B", "C", "D"), length.out = n)),
    x = stats::rnorm(n)
  )
  d$y <- 2 +
    0.5 * d$x +
    stats::rnorm(8L, sd = 1)[as.integer(d$g)] +
    stats::rnorm(n, sd = 0.7)
  d$cnt <- stats::rpois(n, exp(1 + stats::rnorm(8L, sd = 0.5)[as.integer(d$g)]))
  d
}

# A small clustered Gaussian design for brms.
.vr_brms_data <- function(seed = 20260917L, n = 36L) {
  set.seed(seed)
  d <- data.frame(
    site = factor(rep(sprintf("S%02d", 1:6), length.out = n)),
    p1 = stats::rnorm(n),
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n)
  )
  d$y <- 1 +
    0.5 * d$p1 +
    stats::rnorm(6L, sd = 0.7)[as.integer(d$site)] +
    stats::rnorm(n, sd = 0.5)
  d
}

# One chain of 100 draws: the record is under test here, not the posterior, so
# the sampler's small-sample convergence warnings are not.
.vr_brms_fit <- function(d, seed = NULL) {
  suppressWarnings(flexybayes(
    y ~ p1,
    random = ~site,
    data = d,
    backend = "brms",
    n_samples = 100L,
    warmup = 100L,
    chains = 1L,
    seed = seed,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
}

.vr_skip_brms <- function() {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not(
    identical(getOption("brms.backend"), "cmdstanr"),
    "the brms backend of record is not cmdstanr"
  )
}

# The fields a record of each engine carries (a brms record on cmdstanr).
.vr_fields <- list(
  inla = list(
    top = c(
      "r_version",
      "flexybayes_version",
      "engine",
      "engine_versions",
      "threads",
      "seed"
    ),
    engine_versions = "INLA",
    threads = c("num_threads", "omp_num_threads"),
    seed = c("supplied", "as_run")
  ),
  brms = list(
    top = c(
      "r_version",
      "flexybayes_version",
      "engine",
      "backend",
      "engine_versions",
      "threads",
      "seed"
    ),
    engine_versions = c("brms", "cmdstanr", "CmdStan"),
    threads = c("chains", "cores", "threads_per_chain", "omp_num_threads"),
    seed = c("supplied", "as_run")
  )
)

.vr_show <- function(x) {
  if (is.null(x)) "absent" else paste(format(x), collapse = ", ")
}

.vr_compare <- function(label, field, got, want) {
  if (identical(got, want)) {
    return(character(0))
  }
  sprintf(
    "%s: %s is %s, expected %s",
    label,
    field,
    .vr_show(got),
    .vr_show(want)
  )
}

# Every way the record of `fit` differs from what was read here, as text, so a
# failure names the entry point and the field.
.vr_record_problems <- function(
  fit,
  engine,
  seed = NULL,
  chains = 1L,
  label = "fit"
) {
  rec <- fit[["extras"]][["reproducibility"]]
  if (!is.list(rec)) {
    return(sprintf("%s: no record at $extras$reproducibility", label))
  }
  want <- .vr_fields[[engine]]
  out <- c(
    .vr_compare(label, "fields", sort(names(rec)), sort(want$top)),
    .vr_compare(
      label,
      "engine_versions names",
      sort(names(rec$engine_versions)),
      sort(want$engine_versions)
    ),
    .vr_compare(
      label,
      "threads names",
      sort(names(rec$threads)),
      sort(want$threads)
    ),
    .vr_compare(label, "seed names", sort(names(rec$seed)), sort(want$seed)),
    .vr_compare(label, "r_version", rec$r_version, R.version.string),
    .vr_compare(
      label,
      "flexybayes_version",
      rec$flexybayes_version,
      .vr_version("flexyBayes")
    ),
    .vr_compare(label, "engine", rec$engine, engine),
    .vr_compare(
      label,
      "omp_num_threads",
      rec$threads$omp_num_threads,
      .vr_omp()
    ),
    .vr_compare(
      label,
      "seed supplied",
      rec$seed$supplied,
      if (is.null(seed)) NA_integer_ else as.integer(seed)
    )
  )
  if (identical(engine, "inla")) {
    out <- c(
      out,
      .vr_compare(
        label,
        "INLA version",
        rec$engine_versions[["INLA"]],
        .vr_version("INLA")
      ),
      .vr_compare(
        label,
        "num_threads",
        rec$threads$num_threads,
        as.character(fit$inla$.args$num.threads)
      ),
      .vr_compare(label, "seed as_run", rec$seed$as_run, NA_integer_)
    )
  } else {
    stanfit <- fit$brms$fit
    out <- c(
      out,
      .vr_compare(
        label,
        "backend",
        rec$backend,
        getOption("brms.backend", "rstan")
      ),
      .vr_compare(
        label,
        "brms version",
        rec$engine_versions[["brms"]],
        .vr_version("brms")
      ),
      .vr_compare(
        label,
        "cmdstanr version",
        rec$engine_versions[["cmdstanr"]],
        .vr_version("cmdstanr")
      ),
      .vr_compare(
        label,
        "CmdStan version",
        rec$engine_versions[["CmdStan"]],
        .vr_cmdstan_version()
      ),
      .vr_compare(label, "chains", rec$threads$chains, as.integer(chains)),
      .vr_compare(
        label,
        "cores",
        rec$threads$cores,
        as.integer(getOption("mc.cores", 1L))
      ),
      .vr_compare(
        label,
        "threads_per_chain",
        rec$threads$threads_per_chain,
        as.integer(stanfit@stan_args[[1L]]$num_threads)
      ),
      if (is.null(seed)) {
        # A drawn seed has no value to compare with here. The refit in the
        # brms record test shows it is the seed that ran.
        if (
          is.integer(rec$seed$as_run) &&
            length(rec$seed$as_run) == 1L &&
            !is.na(rec$seed$as_run)
        ) {
          character(0)
        } else {
          sprintf(
            "%s: seed as_run is %s, expected one drawn seed",
            label,
            .vr_show(rec$seed$as_run)
          )
        }
      } else {
        .vr_compare(label, "seed as_run", rec$seed$as_run, as.integer(seed))
      }
    )
  }
  out
}

# What print() (or print(summary())) must show of the record.
.vr_print_problems <- function(fit, engine, label = "fit", summary = FALSE) {
  rec <- fit[["extras"]][["reproducibility"]]
  printed <- utils::capture.output(
    if (summary) print(summary(fit)) else print(fit)
  )
  omp <- .vr_omp()
  want <- c(
    R.version.string,
    paste("flexyBayes", .vr_version("flexyBayes")),
    paste("OMP_NUM_THREADS", if (is.na(omp)) "unset" else omp)
  )
  if (identical(engine, "inla")) {
    want <- c(
      want,
      paste("INLA", .vr_version("INLA")),
      paste("INLA num.threads", rec$threads$num_threads),
      "INLA takes no seed"
    )
  } else {
    want <- c(
      want,
      paste("brms", .vr_version("brms")),
      paste("CmdStan", .vr_cmdstan_version()),
      paste("Stan ran with", .vr_show(rec$seed$as_run))
    )
  }
  shown <- vapply(
    want,
    function(w) any(grepl(w, printed, fixed = TRUE)),
    logical(1L)
  )
  if (all(shown)) {
    return(character(0))
  }
  sprintf(
    "%s: %s lacks \"%s\"",
    label,
    if (summary) "summary()" else "print()",
    want[!shown]
  )
}

# Each field a record must carry, as a path into it.
.vr_paths <- function(engine) {
  want <- .vr_fields[[engine]]
  c(
    want$top,
    paste0("engine_versions$", want$engine_versions),
    paste0("threads$", want$threads),
    paste0("seed$", want$seed)
  )
}

.vr_drop <- function(fit, path) {
  parts <- strsplit(path, "$", fixed = TRUE)[[1L]]
  rec <- fit$extras$reproducibility
  if (length(parts) == 1L) {
    rec[[parts]] <- NULL
  } else {
    inner <- rec[[parts[[1L]]]]
    rec[[parts[[1L]]]] <- inner[setdiff(names(inner), parts[[2L]])]
  }
  fit$extras$reproducibility <- rec
  fit
}

# The accessor returns a whole record and refuses every record short of one
# field, naming it.
.vr_expect_contract <- function(fit, engine) {
  expect_identical(.fb_reproducibility(fit), fit$extras$reproducibility)
  for (path in .vr_paths(engine)) {
    expect_error(
      .fb_reproducibility(.vr_drop(fit, path)),
      paste0("`", path, "`"),
      fixed = TRUE
    )
  }
  bare <- fit
  bare$extras$reproducibility <- NULL
  expect_error(
    .fb_reproducibility(bare),
    "no reproducibility record",
    fixed = TRUE
  )
}


# ---------------------------------------------------------------- #
# 1. The record on each engine                                       #
# ---------------------------------------------------------------- #

test_that("an INLA fit records and prints its versions, threads and seed", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .vr_silence()
  d <- .vr_inla_data()

  unseeded <- flexybayes(
    y ~ x,
    random = ~g,
    data = d,
    backend = "inla",
    verbose = FALSE
  )
  seeded <- flexybayes(
    y ~ x,
    random = ~g,
    data = d,
    backend = "inla",
    seed = 11L,
    verbose = FALSE
  )

  expect_identical(
    c(
      .vr_record_problems(unseeded, "inla", label = "no seed"),
      .vr_record_problems(seeded, "inla", seed = 11L, label = "seed 11"),
      .vr_print_problems(unseeded, "inla"),
      .vr_print_problems(unseeded, "inla", summary = TRUE)
    ),
    character(0)
  )
})

test_that("an INLA record missing any field fails the accessor", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .vr_silence()
  fit <- flexybayes(
    y ~ x,
    random = ~g,
    data = .vr_inla_data(),
    backend = "inla",
    verbose = FALSE
  )

  .vr_expect_contract(fit, "inla")

  bare <- fit
  bare$extras$reproducibility <- NULL
  expect_true(any(grepl(
    "not recorded",
    utils::capture.output(print(bare)),
    fixed = TRUE
  )))
})

test_that("a brms fit records its versions, threads and a repeatable seed", {
  .vr_skip_brms()
  skip_if_not_long_tier()
  .vr_silence()
  d <- .vr_brms_data()

  unseeded <- .vr_brms_fit(d)
  seeded <- .vr_brms_fit(d, seed = 1001L)

  expect_identical(
    c(
      .vr_record_problems(unseeded, "brms", label = "no seed"),
      .vr_record_problems(seeded, "brms", seed = 1001L, label = "seed 1001"),
      .vr_print_problems(unseeded, "brms"),
      .vr_print_problems(seeded, "brms", summary = TRUE)
    ),
    character(0)
  )

  # The seed recorded for a run given none is the one Stan ran with: a refit
  # under it returns the same draws.
  as_run <- unseeded$extras$reproducibility$seed$as_run
  expect_false(is.null(as_run) || anyNA(as_run))
  repeated <- .vr_brms_fit(d, seed = as_run)
  expect_identical(as.matrix(repeated$brms), as.matrix(unseeded$brms))

  .vr_expect_contract(seeded, "brms")
})


# ---------------------------------------------------------------- #
# 2. Every entry point                                               #
# ---------------------------------------------------------------- #

test_that("every INLA fit an entry point returns carries the record", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_long_tier()
  .vr_silence()
  d <- .vr_inla_data()

  counts <- new.env(parent = emptyenv())
  counts$engine <- 0L
  real_call <- .inla_call
  testthat::local_mocked_bindings(
    .inla_call = function(args) {
      counts$engine <- counts$engine + 1L
      real_call(args)
    },
    # The restart probe calls INLA for a fit it may discard. It is switched
    # off so that every engine call is a fit that is returned.
    .fb_inla_restart_probe = function(...) NULL
  )

  per_row <- flexybayes(
    y ~ x,
    random = ~g,
    data = d,
    backend = "inla",
    verbose = FALSE
  )
  alias <- fb(
    y ~ x,
    random = ~g,
    data = d,
    backend = "inla",
    seed = 3L,
    verbose = FALSE
  )
  pinned <- fb_inla(y ~ x, random = ~g, data = d, verbose = FALSE)
  aggregated <- flexybayes(
    cnt ~ h,
    random = ~g,
    data = d,
    family = "poisson",
    backend = "inla",
    aggregate = TRUE,
    verbose = FALSE
  )
  streamed <- flexybayes_stream(
    fixed = cnt ~ h,
    random = ~g,
    source = d,
    family = "poisson",
    verbose = FALSE
  )
  # y carries no h effect, so the added sd_h is weakly identified.
  refit <- withr::with_options(
    list(flexyBayes.silence_weak_identification_warning = TRUE),
    stats::update(per_row, random = ~ g + h)
  )
  d_fold <- d
  d_fold$y[seq(1L, nrow(d), by = 4L)] <- NA
  fold <- fb_refit(
    per_row,
    data = d_fold,
    fixed = y ~ x,
    backend = "inla",
    na_action = "augment",
    verbose = FALSE
  )
  # No entry point reaches the aggregated Gaussian emitter
  # (test-aggregation-gaussian-containment.R). Its build site is called here
  # directly, on the IR of a per-row fit of the same model.
  factor_row <- flexybayes(
    y ~ h,
    random = ~g,
    data = d,
    backend = "inla",
    aggregate = FALSE,
    verbose = FALSE
  )
  gaussian_cells <- flexybayes_stream(
    y ~ h,
    random = ~g,
    source = d,
    family = "gaussian",
    fit = FALSE,
    verbose = FALSE
  )
  aggregated_gaussian <- emit_gaussian_aggregated(
    fb = factor_row$fb,
    fb_aggregated = gaussian_cells,
    data = d,
    backend = "inla",
    verbose = FALSE
  )

  fits <- list(
    flexybayes = per_row,
    fb = alias,
    fb_inla = pinned,
    aggregated = aggregated,
    stream = streamed,
    update = refit,
    fold_refit = fold,
    per_row_factor = factor_row,
    aggregated_gaussian = aggregated_gaussian
  )
  seeds <- list(fb = 3L)
  problems <- lapply(names(fits), function(nm) {
    c(
      .vr_record_problems(fits[[nm]], "inla", seed = seeds[[nm]], label = nm),
      .vr_print_problems(fits[[nm]], "inla", label = nm)
    )
  })
  expect_identical(unlist(problems), character(0))
  expect_equal(counts$engine, length(fits))
  expect_equal(sum(lengths(problems) == 0L), counts$engine)
})

test_that("every brms fit an entry point returns carries the record", {
  .vr_skip_brms()
  skip_if_not_long_tier()
  .vr_silence()
  d <- .vr_brms_data()

  counts <- new.env(parent = emptyenv())
  counts$engine <- 0L
  real_brm <- brms::brm
  testthat::local_mocked_bindings(
    brm = function(...) {
      counts$engine <- counts$engine + 1L
      real_brm(...)
    },
    .package = "brms"
  )
  quiet <- function(expr) suppressWarnings(expr)

  per_row <- .vr_brms_fit(d, seed = 1L)
  alias <- quiet(fb(
    y ~ p1,
    random = ~site,
    data = d,
    backend = "brms",
    n_samples = 100L,
    warmup = 100L,
    chains = 1L,
    seed = 2L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  pinned <- quiet(fb_brms(
    y ~ p1 + (1 | site),
    data = d,
    n_samples = 100L,
    warmup = 100L,
    chains = 1L,
    seed = 3L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  refit <- quiet(stats::update(per_row, seed = 4L))
  review <- flexybayes(
    y ~ p1,
    random = ~site,
    data = d,
    backend = "brms",
    n_samples = 100L,
    warmup = 100L,
    chains = 1L,
    review_code = TRUE,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )
  proceeded <- quiet(proceed(review))
  fold <- quiet(fb_refit(per_row, data = d, fixed = y ~ p1, backend = "brms",
                        seed = 5L, verbose = FALSE, mcmc_verbose = FALSE))

  fits <- list(
    flexybayes = per_row,
    fb = alias,
    fb_brms = pinned,
    update = refit,
    proceed = proceeded,
    fold_refit = fold
  )
  seeds <- list(
    flexybayes = 1L,
    fb = 2L,
    fb_brms = 3L,
    update = 4L,
    fold_refit = 5L
  )
  problems <- lapply(names(fits), function(nm) {
    c(
      .vr_record_problems(fits[[nm]], "brms", seed = seeds[[nm]], label = nm),
      .vr_print_problems(fits[[nm]], "brms", label = nm)
    )
  })
  expect_identical(unlist(problems), character(0))
  expect_equal(counts$engine, length(fits))
  expect_equal(sum(lengths(problems) == 0L), counts$engine)
})


# ---------------------------------------------------------------- #
# 3. The value the engine ran with, not the option                   #
# ---------------------------------------------------------------- #

# No entry point passes a thread setting to an engine. The emitters' `...` do.
# Each fit below runs on one thread or core while the option says two, so a
# record read from the option before the fit is wrong.

test_that("INLA's record holds the num.threads it ran with, not the option", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .vr_silence()
  base <- flexybayes(
    y ~ x,
    random = ~g,
    data = .vr_inla_data(),
    backend = "inla",
    aggregate = FALSE,
    verbose = FALSE
  )

  old <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "2:1")
  withr::defer(INLA::inla.setOption(num.threads = old))
  fit <- emit_inla(
    fb = base$fb,
    data = base$data,
    verbose = FALSE,
    num.threads = "1:1"
  )

  expect_identical(INLA::inla.getOption("num.threads"), "2:1")
  expect_identical(fit$inla$.args$num.threads, "1:1")
  expect_identical(fit$extras$reproducibility$threads$num_threads, "1:1")
})

test_that("a brms record holds the cores brm() was given, not the option", {
  .vr_skip_brms()
  skip_if_not_long_tier()
  .vr_silence()
  base <- .vr_brms_fit(.vr_brms_data(), seed = 7L)

  given <- NULL
  real_brm <- brms::brm
  testthat::local_mocked_bindings(
    brm = function(...) {
      given <<- list(...)$cores
      real_brm(...)
    },
    .package = "brms"
  )
  withr::local_options(mc.cores = 2L)
  fit <- suppressWarnings(emit_brms(
    fb = base$fb,
    data = base$data,
    n_samples = 100L,
    warmup = 100L,
    chains = 1L,
    seed = 7L,
    verbose = FALSE,
    mcmc_verbose = FALSE,
    cores = 1L
  ))

  expect_identical(given, 1L)
  expect_identical(fit$extras$reproducibility$threads$cores, 1L)
})
