# Tests for the compiled model manifest (R/model_manifest.R).
#
# Every fit records which stated or defaulted prior became which engine slot,
# what the engine reports it received, the route the fit took, and whether its
# marginal likelihood is available. The manifest lives at
# `fit$extras$manifest` and is read through `.fb_manifest()`.
#
# The eight assertions T1-T8 are the manifest specification's acceptance
# tests. The engine read-back is the engine's own record: `result$all.hyper`
# on INLA, the object INLA::inla.priors.used() prints, and
# brms::prior_summary() on brms. A slot the package did not emit is classified
# against the default the engine declares for it, so the read-back is not a
# copy of what was sent.

.man_data <- function() {
  set.seed(20260915L)
  d <- data.frame(
    a = factor(rep(seq_len(3L), each = 40L)),
    b = factor(rep(rep(seq_len(2L), each = 20L), 3L)),
    cc = factor(rep(seq_len(4L), 30L)),
    env = factor(rep(c("E1", "E2"), 60L)),
    g = factor(rep(sprintf("g%02d", seq_len(12L)), each = 10L)),
    x = stats::rnorm(120L)
  )
  d$y <- 10 + stats::rnorm(12L)[as.integer(d$g)] + 0.5 * d$x +
    stats::rnorm(120L, sd = c(1, 2)[as.integer(d$env)])
  d$count <- stats::rpois(
    120L,
    exp(1 + stats::rnorm(12L, sd = 0.3)[as.integer(d$g)])
  )
  d
}

# The namespace a warning was raised in: the deepest frame below the warning
# signal whose function is not base R's. `$` partial matching is reported
# from inside the engines too (INLA 26.8.7 matches `gcpo$type` to `type.cv`
# on every fit), and those are not this package's to fix.
.man_warning_origin <- function() {
  calls <- sys.calls()
  signal_at <- which(vapply(
    calls,
    function(cl) identical(cl[[1L]], quote(.signalSimpleWarning)),
    logical(1L)
  ))
  if (!length(signal_at)) {
    return(NA_character_)
  }
  for (k in rev(seq_len(max(signal_at) - 1L))) {
    env <- environment(sys.function(k))
    if (is.null(env)) {
      next
    }
    origin <- environmentName(topenv(env))
    if (!identical(origin, "base")) {
      return(origin)
    }
  }
  NA_character_
}

# Evaluates `expr` with `$` partial matching reported, and turns a partial
# match raised in this package's code into an error. Every other
# warning is muffled: the fits below are small and the numerical-confirm gate
# may flag them, which is not what these tests are about.
.man_strict <- function(expr) {
  withr::local_options(
    warnPartialMatchDollar = TRUE,
    flexyBayes.silence_default_prior_note = TRUE
  )
  withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      if (
        grepl("partial match of", conditionMessage(w), fixed = TRUE) &&
          identical(.man_warning_origin(), "flexyBayes")
      ) {
        stop(
          "partial match on the emitter path: ",
          conditionMessage(w),
          call. = FALSE
        )
      }
      invokeRestart("muffleWarning")
    }
  )
}

.man_of <- function(fit) {
  manifest <- flexyBayes:::.fb_manifest(fit)
  expect_false(is.null(manifest))
  manifest
}

.man_package_slots <- function(manifest) {
  manifest$slots[manifest$slots$source == "package", , drop = FALSE]
}


# --- T1 ---------------------------------------------------------------------

test_that("T1: no prior input goes unused without a reason, on every route", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .man_data()

  per_row <- .man_strict(flexybayes(
    fixed = y ~ x, random = ~ g, data = d, backend = "inla",
    prior = fb_prior(
      sigma ~ uniform(lower = 0, upper = 30),
      sd(group = "g") ~ uniform(lower = 0, upper = 20),
      b("x") ~ normal(0, 5)
    ),
    verbose = FALSE
  ))
  aggregated <- .man_strict(flexybayes(
    fixed = count ~ 1, random = ~ g, family = "poisson", data = d,
    backend = "inla", aggregate = TRUE, verbose = FALSE
  ))
  streamed <- .man_strict(flexybayes_stream(
    count ~ 1, random = ~ g, source = d, family = "poisson",
    backend = "inla", verbose = FALSE
  ))

  for (fit in list(per_row, aggregated, streamed)) {
    inputs <- .man_of(fit)$inputs
    expect_gt(nrow(inputs), 0L)
    unused <- inputs[inputs$n_slots == 0L, , drop = FALSE]
    expect_true(all(!is.na(unused$reason_code) & nzchar(unused$reason_code)))
    expect_true(all(unused$cardinality == "not_consumed"))
  }

  # Each stated prior on the per-row fit reaches exactly one slot.
  inputs <- .man_of(per_row)$inputs
  expect_setequal(inputs$input_id[inputs$source == "explicit"],
                  c("sigma:sigma", "sd:g", "b:x"))
  expect_true(all(inputs$n_slots[inputs$source == "explicit"] == 1L))

  # The package default states a residual prior on every family, and a
  # Poisson likelihood has no residual scale to carry it.
  for (fit in list(aggregated, streamed)) {
    inputs <- .man_of(fit)$inputs
    sigma <- inputs[inputs$input_id == "sigma:sigma", , drop = FALSE]
    expect_equal(nrow(sigma), 1L)
    expect_identical(sigma$cardinality, "not_consumed")
    expect_identical(sigma$reason_code, "prior_family_has_no_residual_scale")
  }

  # prior_vc_sd written beside an fb_prior() reaches no slot on INLA: the
  # standard deviation the prior leaves unstated takes the package default.
  superseded <- .man_strict(flexybayes(
    fixed = y ~ x, random = ~ g + a, data = d, backend = "inla",
    prior = fb_prior(sd(group = "g") ~ uniform(lower = 0, upper = 20)),
    prior_vc_sd = 2, verbose = FALSE
  ))
  inputs <- .man_of(superseded)$inputs
  scalar <- inputs[inputs$input_id == "scalar:prior_vc_sd", , drop = FALSE]
  expect_equal(nrow(scalar), 1L)
  expect_identical(scalar$cardinality, "not_consumed")
  expect_identical(scalar$reason_code, "prior_vc_sd_superseded_by_prior")
  filled <- inputs[inputs$input_id == "sd:a", , drop = FALSE]
  expect_identical(filled$source, "filled_default")
  expect_equal(filled$n_slots, 1L)
})


# --- T2 and T3 --------------------------------------------------------------

.man_inla_fits <- function(d) {
  list(
    sectioned = flexybayes(
      fixed = y ~ x, random = ~ g, residual = ~ dsum(~ units | env),
      data = d, backend = "inla", verbose = FALSE
    ),
    heterogeneous = flexybayes(
      fixed = y ~ 1, random = ~ diag(env):g, data = d, backend = "inla",
      prior = fb_prior(sd(group = "g") ~ uniform(lower = 0, upper = 15)),
      verbose = FALSE
    ),
    slope = flexybayes(y ~ x + (x || g), data = d, backend = "inla",
                       verbose = FALSE),
    interaction = flexybayes(
      fixed = y ~ 1, random = ~ a:b:cc, data = d, backend = "inla",
      prior = fb_prior(sd(group = "a:b:cc") ~ uniform(lower = 0, upper = 20)),
      verbose = FALSE
    ),
    legacy = flexybayes(
      fixed = y ~ x, random = ~ g, data = d, backend = "inla",
      prior_vc_sd = 2, prior_fixed_sd = 50, verbose = FALSE
    ),
    counts = flexybayes(
      fixed = count ~ 1, random = ~ g, family = "poisson", data = d,
      backend = "inla", aggregate = TRUE, verbose = FALSE
    )
  )
}

test_that("T2: every slot the package emitted names the input it came from", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .man_data()
  fits <- .man_strict(.man_inla_fits(d))

  for (nm in names(fits)) {
    manifest <- .man_of(fits[[nm]])
    emitted <- .man_package_slots(manifest)
    expect_gt(nrow(emitted), 0L)
    expect_true(all(!is.na(emitted$input_id)), label = nm)
    expect_true(all(emitted$input_id %in% manifest$inputs$input_id), label = nm)
    expect_false(anyDuplicated(manifest$slots$slot_id) > 0L, label = nm)
  }
})

test_that("T3: the emitted slots are the ones INLA reports, and nothing else", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .man_data()
  fits <- .man_strict(.man_inla_fits(d))

  for (nm in names(fits)) {
    manifest <- .man_of(fits[[nm]])
    # Every hyperparameter INLA counts in its mode is in the read-back.
    expect_true(isTRUE(manifest$readback_complete), label = nm)
    # Nothing the engine reports is outside what was emitted or declared as
    # the engine's own default.
    expect_false(any(manifest$slots$source == "unexpected"), label = nm)
    # Everything emitted arrived as sent.
    emitted <- .man_package_slots(manifest)
    expect_true(all(emitted$readback_agrees %in% TRUE), label = nm)
  }

  # The residual precision the interaction fit left unstated carries INLA's
  # own default, and the manifest names it as that.
  slots <- .man_of(fits$interaction)$slots
  family_slot <- slots[slots$section == "family", , drop = FALSE]
  expect_equal(nrow(family_slot), 1L)
  expect_identical(family_slot$source, "engine_default")
  expect_match(family_slot$readback, "loggamma", fixed = TRUE)
})

test_that("T3: brms reports every emitted prior row back as the user's", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  d <- .man_data()
  # The parsed model is engine-agnostic; an INLA fit is the quickest way to
  # hold one with the stated prior attached.
  fb_ir <- .man_strict(flexybayes(
    fixed = y ~ x, random = ~ g, data = d, backend = "inla",
    prior = fb_prior(
      sigma ~ uniform(lower = 0, upper = 30),
      sd(group = "g") ~ half_normal(scale = 3),
      b("x") ~ normal(0, 5)
    ),
    verbose = FALSE
  ))$extras$fb_terms
  rows <- flexyBayes:::.priors_to_brms(fb_ir$priors, fb_ir)
  unsampled <- brms::brm(
    y ~ x + (1 | g), data = d, prior = rows, empty = TRUE, silent = 2L
  )
  manifest <- flexyBayes:::.fb_manifest_brms(
    fb = fb_ir, emitted_prior = rows, brmsfit = unsampled
  )
  emitted <- .man_package_slots(manifest)
  expect_equal(nrow(emitted), 3L)
  expect_true(all(emitted$readback_agrees %in% TRUE))
  expect_setequal(emitted$input_id, c("sigma:sigma", "sd:g", "b:x"))
  expect_false(any(manifest$slots$source == "unexpected"))
  expect_identical(manifest$evidence$status, "unavailable")
  expect_identical(manifest$evidence$reason_code, "evidence_not_computed")
})


# --- T4 and T5 --------------------------------------------------------------

test_that("T4: reported priors, route and evidence trace to the manifest", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .man_data()
  per_row <- .man_strict(flexybayes(
    fixed = y ~ x, random = ~ g, data = d, backend = "inla",
    prior = fb_prior(sigma ~ uniform(lower = 0, upper = 30)),
    verbose = FALSE
  ))
  counts <- .man_strict(flexybayes(
    fixed = count ~ 1, random = ~ g, family = "poisson", data = d,
    backend = "inla", aggregate = TRUE, verbose = FALSE
  ))

  for (fit in list(per_row, counts)) {
    manifest <- .man_of(fit)
    ps <- prior_summary(fit)
    # The summary's slot table is the manifest's.
    expect_identical(ps$slots$slot_id, manifest$slots$slot_id)
    expect_identical(ps$slots$emitted, manifest$slots$emitted)
    # Every emitted slot is printed with its identifier, and every unused
    # input with its reason.
    printed <- paste(utils::capture.output(print(ps)), collapse = "\n")
    for (id in .man_package_slots(manifest)$slot_id) {
      expect_match(printed, id, fixed = TRUE)
    }
    unused <- manifest$inputs[manifest$inputs$n_slots == 0L, , drop = FALSE]
    for (code in unused$reason_code) {
      expect_match(printed, code, fixed = TRUE)
    }
  }

  # Evidence: the value summary() reports is the manifest's, and on the
  # aggregated route so is the reason it is absent.
  conv <- summary(per_row)$converge
  expect_identical(.man_of(per_row)$evidence$status, "available")
  expect_equal(conv$mlik, .man_of(per_row)$evidence$value)
  conv <- summary(counts)$converge
  expect_true(is.na(conv$mlik))
  expect_identical(
    conv$evidence_unavailable,
    .man_of(counts)$evidence$reason_code
  )

  # Route: the manifest's route agrees with the dispatch record.
  expect_identical(.man_of(per_row)$route$route, "per_row")
  expect_false(identical(backend_decision(per_row)$path, "aggregated_count"))
  expect_identical(.man_of(counts)$route$route, "aggregated_count")
  expect_identical(backend_decision(counts)$path, "aggregated_count")

  # The aggregated fit's prior label is read from the input sources.
  expect_identical(
    counts$extras$aggregation_meta$prior_parametrization,
    "package_default"
  )
})

test_that("T5: stated priors round-trip; one-to-many inputs list each slot", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .man_data()
  sectioned <- .man_strict(flexybayes(
    fixed = y ~ x, random = ~ g, residual = ~ dsum(~ units | env),
    data = d, backend = "inla",
    prior = fb_prior(sigma ~ uniform(lower = 0, upper = 30)),
    verbose = FALSE
  ))
  manifest <- .man_of(sectioned)
  sigma <- manifest$inputs[manifest$inputs$input_id == "sigma:sigma", ]
  expect_identical(sigma$cardinality, "one_to_many")
  expect_equal(sigma$n_slots, nlevels(d$env))
  expect_match(sigma$expansion, "2", fixed = TRUE)
  ps <- prior_summary(sectioned)
  sigma_slots <- ps$slots[ps$slots$input_id %in% "sigma:sigma", , drop = FALSE]
  expect_equal(nrow(sigma_slots), 2L)
  # The bound as INLA received it: the stated 30 on the standardised response.
  u_readback <- as.numeric(
    sub("^.*U=([-+0-9.eE]+);.*$", "\\1", sigma_slots$readback)
  )
  expect_equal(u_readback, rep(30 / stats::sd(d$y), 2L), tolerance = 1e-10)
  # INLA applies the uniform on each section's standard deviation as stated,
  # so the summary does not describe the log-sigma retarget brms performs.
  expect_null(ps$residual_lowered_to)
  printed <- utils::capture.output(print(ps))
  expect_false(any(grepl("retargeted", printed, fixed = TRUE)))

  heterogeneous <- .man_strict(flexybayes(
    fixed = y ~ 1, random = ~ diag(env):g, data = d, backend = "inla",
    prior = fb_prior(sd(group = "g") ~ uniform(lower = 0, upper = 15)),
    verbose = FALSE
  ))
  inputs <- .man_of(heterogeneous)$inputs
  sd_g <- inputs[inputs$input_id == "sd:g", ]
  expect_identical(sd_g$cardinality, "one_to_many")
  expect_equal(sd_g$n_slots, nlevels(d$env))
})


# --- T6 ---------------------------------------------------------------------

test_that("T6: on a three-way interaction only its own key takes a prior", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .man_data()
  for (key in c("a", "b", "cc", "a:b", "zz")) {
    # fb_prior() reads its formulas unevaluated, so the key is substituted
    # into the call; `sd(group = key)` would state a group named "key".
    prior <- eval(bquote(
      fb_prior(sd(group = .(key)) ~ uniform(lower = 0, upper = 20))
    ))
    expect_identical(prior$specs[[1L]]$target$group, key)
    err <- tryCatch(
      .man_strict(flexybayes(
        fixed = y ~ 1, random = ~ a:b:cc, data = d, backend = "inla",
        prior = prior, return_code = TRUE
      )),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_prior_target_absent")
    expect_identical(err$reason_code, "prior_target_not_in_model")
  }

  fit <- .man_strict(flexybayes(
    fixed = y ~ 1, random = ~ a:b:cc, data = d, backend = "inla",
    prior = fb_prior(sd(group = "a:b:cc") ~ uniform(lower = 0, upper = 20)),
    verbose = FALSE
  ))
  manifest <- .man_of(fit)
  slot <- manifest$slots[
    manifest$slots$input_id %in% "sd:a:b:cc", ,
    drop = FALSE
  ]
  expect_equal(nrow(slot), 1L)
  expect_match(slot$engine_slot, "a_b_cc_id", fixed = TRUE)
  expect_true(isTRUE(slot$readback_agrees))

  # The aggregated count route checks the same thing against its own slots.
  err <- tryCatch(
    .man_strict(flexybayes(
      fixed = count ~ 1, random = ~ g, family = "poisson", data = d,
      backend = "inla", aggregate = TRUE,
      prior = fb_prior(sd(group = "zz") ~ uniform(lower = 0, upper = 3)),
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_prior_target_absent")
  expect_identical(err$reason_code, "prior_target_not_in_model")
})


# --- T7 ---------------------------------------------------------------------

test_that("T7: every emitter call site returns a manifest or returns code", {
  dir_r <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(dir_r), "package sources are not available")
  emitters <- c(
    "emit_inla", "emit_brms", "emit_count_aggregated",
    "emit_gaussian_aggregated", ".backend_emit_fn"
  )
  sites <- list()
  for (path in list.files(dir_r, pattern = "\\.R$", full.names = TRUE)) {
    pd <- utils::getParseData(parse(path, keep.source = TRUE))
    hits <- pd[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text %in% emitters, ]
    for (i in seq_len(nrow(hits))) {
      call_expr <- pd$parent[pd$id == hits$parent[[i]]]
      sites[[length(sites) + 1L]] <- data.frame(
        file = basename(path),
        emitter = hits$text[[i]],
        text = utils::getParseText(pd, call_expr)
      )
    }
  }
  sites <- do.call(rbind, sites)
  counts <- table(paste(sites$file, sites$emitter))
  expect_equal(
    as.list(counts),
    list(
      "dispatch.R .backend_emit_fn" = 3L,
      "dispatch.R emit_count_aggregated" = 1L,
      "flexybayes.R emit_brms" = 1L,
      "stream_aggregate.R emit_count_aggregated" = 1L
    ),
    ignore_attr = TRUE
  )
  expect_equal(nrow(sites), 6L)
  # No R/ code reaches the aggregated Gaussian emit.
  expect_false("emit_gaussian_aggregated" %in% sites$emitter)
  # The one direct emit_brms() call returns Stan code.
  review <- sites[sites$emitter == "emit_brms", ]
  expect_match(review$text, "return_code = TRUE", fixed = TRUE)
  # The resolved emitters are the two that attach a manifest.
  expect_identical(
    flexyBayes:::.backend_emit_fn("inla"),
    flexyBayes:::emit_inla
  )
  expect_identical(
    flexyBayes:::.backend_emit_fn("brms"),
    flexyBayes:::emit_brms
  )
  # Each emitter function attaches the manifest itself.
  for (fn in c("emit_inla", "emit_brms", ".agg_count_build_fit")) {
    body_txt <- paste(deparse(get(fn, envir = asNamespace("flexyBayes"))),
                      collapse = " ")
    expect_match(body_txt, "manifest", fixed = TRUE, label = fn)
  }
})

test_that("T7: a live brms fit carries a manifest that agrees with brms", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  d <- .man_data()
  fit <- .man_strict(flexybayes(
    fixed = y ~ x, random = ~ g, data = d, backend = "brms",
    prior = fb_prior(sigma ~ uniform(lower = 0, upper = 30)),
    n_samples = 100L, warmup = 100L, chains = 1L, seed = 1L,
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  manifest <- .man_of(fit)
  expect_identical(manifest$engine, "brms")
  emitted <- .man_package_slots(manifest)
  expect_gt(nrow(emitted), 0L)
  expect_true(all(emitted$readback_agrees %in% TRUE))
  expect_false(any(manifest$slots$source == "unexpected"))
})


# --- T8 ---------------------------------------------------------------------

test_that("T8: the emitter paths make no partial match of a list name", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  d <- .man_data()
  expect_no_error(.man_strict(flexybayes(
    fixed = y ~ 1, random = ~ a:b:cc, data = d, backend = "inla",
    prior = fb_prior(sd(group = "a:b:cc") ~ uniform(lower = 0, upper = 20)),
    return_code = TRUE
  )))
  expect_no_error(.man_strict(.man_inla_fits(d)))
  expect_no_error(.man_strict(flexybayes_stream(
    count ~ 1, random = ~ g, source = d, family = "poisson",
    backend = "inla", verbose = FALSE
  )))
})


# --- registry ---------------------------------------------------------------

test_that("every translation-registry row states its prior cardinality", {
  reg <- flexyBayes:::.fb_translation_registry()
  expect_true("cardinality" %in% names(reg))
  head_word <- sub(":.*$", "", reg$cardinality)
  expect_true(all(head_word %in% flexyBayes:::.FB_TRANSLATION_CARDINALITIES))
})
