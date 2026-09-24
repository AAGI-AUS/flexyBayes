# =============================================================================
# What a refusal says about the engines, checked through the public verbs.
#
# A refusal names what cannot be fitted, why, and what to do instead. When
# the reason is flexyBayes's own -- no lowering written, a result that did
# not repeat -- the message says so, and does not describe it as something
# INLA or brms cannot do. Before 0.14.0 several messages did: "no active
# engine represents" a multi-trait covariance that brms writes directly,
# "INLA's f() has no equicorrelation model" when INLA has an exchangeable
# group model, "INLA takes the family's default link" when INLA takes a
# probit link, and "brms's mi() term is admitted only for the gaussian and
# lognormal families" when brms takes a gamma response too. The engine facts
# the messages still state are checked against the engines in
# test-engine-claims.R.
#
# Every call below goes through flexybayes(), fb() or fb_select() with
# return_code = TRUE or plan = TRUE, so nothing is sampled.
# =============================================================================

.rw_data <- function() {
  set.seed(21L)
  d <- expand.grid(
    rep = factor(seq_len(2L)),
    gen = factor(seq_len(6L)),
    env = factor(seq_len(3L))
  )
  d$trait <- factor(rep(seq_len(2L), length.out = nrow(d)))
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rnorm(nrow(d), 10, 2)
  d
}

.rw_grid <- function() {
  set.seed(22L)
  g <- expand.grid(row = factor(seq_len(6L)), col = factor(seq_len(5L)))
  g$y <- stats::rnorm(nrow(g))
  g
}

.rw_kinship <- function() {
  k <- diag(6)
  dimnames(k) <- list(as.character(seq_len(6L)), as.character(seq_len(6L)))
  k
}

# The condition a public call raises, or NULL when it does not refuse.
.rw_refusal <- function(..., backend, known_matrices = list(),
                        return_code = TRUE) {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE
  )
  res <- tryCatch(
    suppressMessages(flexybayes(
      ...,
      known_matrices = known_matrices,
      backend = backend,
      return_code = return_code,
      verbose = FALSE
    )),
    error = function(e) e
  )
  if (inherits(res, "error")) res else NULL
}

# Phrases no refusal uses any more: each described flexyBayes's own gap as a
# limitation of an engine.
.rw_retired <- c(
  "no active engine", "no active backend", "cannot represent",
  "neither active engine", "INLA cannot"
)

.rw_expect_plain <- function(msg, label) {
  for (phrase in .rw_retired) {
    expect_false(
      grepl(phrase, msg, fixed = TRUE),
      label = paste0(label, " does not say '", phrase, "'")
    )
  }
}

# --- structures flexyBayes does not fit on either engine ----------------- #

test_that("a multi-trait covariance refusal points at brms rather than denying it", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  d <- .rw_data()
  for (engine in c("inla", "brms")) {
    err <- .rw_refusal(
      y ~ env,
      random = ~ us(trait):vm(gen),
      data = d,
      known_matrices = list(K = .rw_kinship()),
      backend = engine
    )
    expect_s3_class(err, "flexybayes_refusal_interaction_not_representable")
    msg <- conditionMessage(err)
    expect_match(msg, "flexyBayes has no lowering for a multi-trait",
      fixed = TRUE)
    expect_match(msg, "(0 + trait | gr(gen, cov = K))", fixed = TRUE)
    .rw_expect_plain(msg, engine)
  }
})

test_that("corh() refuses by name on both engines and names the equivalent it fits", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  d <- .rw_data()
  for (engine in c("inla", "brms")) {
    err <- .rw_refusal(y ~ env, random = ~ corh(env):gen, data = d,
      backend = engine)
    expect_s3_class(
      err, "flexybayes_refusal_corh_no_equicorrelation_representation"
    )
    msg <- conditionMessage(err)
    expect_match(msg, "flexyBayes has no lowering for it on either engine",
      fixed = TRUE)
    expect_match(msg, "random = ~ gen + env:gen", fixed = TRUE)
    # INLA used to answer with a gate refusal that pointed at brms, which
    # then refused the same model.
    expect_false(grepl("backend = \"brms\"", msg, fixed = TRUE))
    .rw_expect_plain(msg, engine)
  }
})

test_that("cor() refuses without a claim about INLA's models", {
  skip_if_not_installed("INLA")
  d <- .rw_data()
  err <- .rw_refusal(y ~ 1, random = ~ cor(env), data = d, backend = "inla")
  expect_s3_class(
    err, "flexybayes_refusal_corh_no_equicorrelation_representation"
  )
  msg <- conditionMessage(err)
  expect_match(msg, "flexyBayes has no lowering for it", fixed = TRUE)
  expect_false(grepl("INLA's f()", msg, fixed = TRUE))
  .rw_expect_plain(msg, "cor()")
})

test_that("a factor-analytic refusal says the brms gap is flexyBayes's", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  d <- .rw_data()
  for (engine in c("inla", "brms")) {
    err <- .rw_refusal(y ~ env, random = ~ fa(env, 1):gen, data = d,
      backend = engine)
    expect_s3_class(err, "flexybayes_refusal_factor_analytic_not_repeatable")
    msg <- conditionMessage(err)
    expect_match(msg, "flexyBayes has no factor-analytic lowering for brms",
      fixed = TRUE)
    expect_false(grepl("brms has no", msg, fixed = TRUE))
    expect_false(grepl("0.10.0.9000", msg, fixed = TRUE))
    expect_match(msg, "diag(env):gen", fixed = TRUE)
    .rw_expect_plain(msg, engine)
  }
  bare <- .rw_refusal(y ~ 1, random = ~ fa(env, 1), data = d,
    backend = "inla")
  expect_s3_class(bare, "flexybayes_refusal_fa_not_representable")
  # fa(f, k):g is refused too, so the bare refusal no longer points at it.
  expect_false(grepl("INLA emits", conditionMessage(bare), fixed = TRUE))
})

test_that("ar2() says flexyBayes writes first-order terms only", {
  skip_if_not_installed("INLA")
  g <- .rw_grid()
  err <- .rw_refusal(y ~ 1, random = ~ ar2(row), data = g, backend = "inla")
  expect_s3_class(err, "flexybayes_refusal_ar2_not_representable")
  msg <- conditionMessage(err)
  expect_match(msg, "flexyBayes writes first-order autoregressive terms only",
    fixed = TRUE)
  expect_match(msg, "INLA's own ar model takes higher orders", fixed = TRUE)
  .rw_expect_plain(msg, "ar2()")
})

test_that("the residual AR1 refusal describes flexyBayes's field, not INLA", {
  skip_if_not_installed("INLA")
  g <- .rw_grid()
  err <- .rw_refusal(y ~ 1, residual = ~ ar1(row):ar1(col), data = g,
    backend = "inla")
  expect_s3_class(err, "flexybayes_refusal_ar1_residual_not_representable")
  msg <- conditionMessage(err)
  expect_match(msg, "flexyBayes fits an autoregressive process", fixed = TRUE)
  expect_false(grepl("INLA represents", msg, fixed = TRUE))
  expect_match(msg, "random = ~ ar1(row):ar1(col)", fixed = TRUE)
})

test_that("a low-rank smooth refusal points at a smooth that fits", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("mgcv")
  set.seed(23L)
  d <- data.frame(x = stats::runif(40L), y = stats::rnorm(40L))
  err <- .rw_refusal(
    y ~ 1,
    random = ~ s(x, representation = list(scheme = "low_rank_smooth",
      rank = 4L)),
    data = d,
    backend = "inla",
    return_code = FALSE
  )
  expect_identical(err$reason_code, "low_rank_smooth_unsupported")
  msg <- conditionMessage(err)
  expect_match(msg, "write the smooth as spl(x), which fits on INLA",
    fixed = TRUE)
  # The exact s() smooth it used to point at is refused on both engines.
  expect_false(grepl("fit the exact smooth", msg, fixed = TRUE))
  .rw_expect_plain(msg, "low-rank smooth")
})

# --- the brms route ------------------------------------------------------ #

test_that("a factor slope refusal states what fits on INLA", {
  skip_if_not_installed("INLA")
  d <- .rw_data()
  err <- .rw_refusal(y ~ env + (env | gen), data = d, backend = "inla")
  expect_s3_class(err, "flexybayes_refusal_brms_factor_random_slope_unsupported")
  msg <- conditionMessage(err)
  expect_match(msg, "On INLA, diag(env):gen fits and", fixed = TRUE)
  expect_false(grepl("f() machinery", msg, fixed = TRUE))
})

test_that("an mgcv smooth on brms says brms itself fits s()", {
  skip_if_not_installed("brms")
  d <- .rw_data()
  err <- .rw_refusal(y ~ 1, random = ~ s(x), data = d, backend = "brms")
  expect_s3_class(err, "flexybayes_refusal_brms_cannot_represent_term")
  msg <- conditionMessage(err)
  expect_match(msg, "flexyBayes cannot fit the random term type", fixed = TRUE)
  expect_match(msg, "brms fits s() in its own formula", fixed = TRUE)
  .rw_expect_plain(msg, "s(x) on brms")
})

test_that("an AR1 field on brms is flexyBayes's gap", {
  skip_if_not_installed("brms")
  g <- .rw_grid()
  err <- .rw_refusal(y ~ 1, random = ~ ar1(row):ar1(col), data = g,
    backend = "brms")
  expect_s3_class(err, "flexybayes_refusal_stan_cannot_represent_ar1_field")
  msg <- conditionMessage(err)
  expect_match(msg, "It has no brms lowering for the field", fixed = TRUE)
  expect_false(grepl("brms has no", msg, fixed = TRUE))
})

test_that("a blocks carrier is refused on each route in flexyBayes's words", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  d <- .rw_data()
  blocks <- list(diag(3), diag(3))
  spell <- ~ vm(gen, cov = fb_cov(B, type = "blocks"))

  brms_err <- .rw_refusal(y ~ env, random = spell, data = d,
    known_matrices = list(B = blocks), backend = "brms")
  expect_s3_class(brms_err, "flexybayes_refusal_stan_cannot_represent_structured_cov")
  brms_msg <- conditionMessage(brms_err)
  expect_match(brms_msg, "flexyBayes cannot fit this model on brms", fixed = TRUE)
  # INLA refuses a blocks carrier too, so the message no longer sends the
  # user there.
  expect_false(grepl("backend = \"inla\"", brms_msg, fixed = TRUE))
  .rw_expect_plain(brms_msg, "brms")

  inla_err <- .rw_refusal(y ~ env, random = spell, data = d,
    known_matrices = list(B = blocks), backend = "inla")
  expect_s3_class(inla_err, "flexybayes_refusal_inla_gate_refused")
  expect_match(conditionMessage(inla_err),
    "flexyBayes's INLA route would fit each block", fixed = TRUE)

  # auto: neither route fits it. The one-time routing note and the refusal
  # both say so without a claim about the engines.
  flexyBayes:::.reset_emit_state_for_test()
  msgs <- character(0)
  auto_err <- withr::with_options(
    list(
      flexyBayes.silence_default_prior_note = TRUE,
      flexyBayes.silence_auto_fallback_note = FALSE
    ),
    withCallingHandlers(
      tryCatch(
        flexybayes(y ~ env, random = spell, data = d,
          known_matrices = list(B = blocks), backend = "auto",
          verbose = FALSE),
        error = function(e) e
      ),
      message = function(m) {
        msgs[[length(msgs) + 1L]] <<- conditionMessage(m)
        invokeRestart("muffleMessage")
      }
    )
  )
  expect_s3_class(auto_err, "flexybayes_refusal_auto_no_active_route")
  auto_msg <- conditionMessage(auto_err)
  expect_match(auto_msg, "flexyBayes cannot fit this model on INLA or on brms",
    fixed = TRUE)
  .rw_expect_plain(auto_msg, "auto")
  note <- msgs[grepl("backend = \"auto\"", msgs, fixed = TRUE)]
  expect_length(note, 1L)
  expect_match(note, "flexyBayes cannot fit it on brms either", fixed = TRUE)
  .rw_expect_plain(note, "auto note")
})

test_that("a missing count response on brms names the families flexyBayes checked", {
  skip_if_not_installed("brms")
  d <- .rw_data()
  d$cnt <- stats::rpois(nrow(d), 3)
  d$cnt[[1L]] <- NA
  err <- .rw_refusal(cnt ~ env, random = ~gen, data = d, family = "poisson",
    backend = "brms")
  expect_s3_class(err, "flexybayes_refusal_brms_cannot_augment_nongaussian")
  msg <- conditionMessage(err)
  expect_match(msg, "gaussian and lognormal families, the two checked",
    fixed = TRUE)
  # brms's mi() takes a gamma response as well (test-engine-claims.R).
  expect_false(grepl("admitted only for", msg, fixed = TRUE))
})

test_that("an untranslated brms prior target no longer cites a refusal that is gone", {
  skip_if_not_installed("brms")
  d <- .rw_data()
  err <- .rw_refusal(y ~ x + (x | gen), data = d, backend = "brms",
    prior = fb_prior(cor(group = "gen") ~ lkj(2)))
  expect_s3_class(err, "flexybayes_refusal_prior_not_translatable_for_backend")
  msg <- conditionMessage(err)
  expect_match(msg, "brms keeps its own LKJ prior with eta = 1", fixed = TRUE)
  # Correlated random slopes fit on brms now.
  expect_false(grepl("refused at ingest", msg, fixed = TRUE))
  expect_match(msg, "flexyBayes translates these priors for backend",
    fixed = TRUE)
})

# --- the INLA gate ------------------------------------------------------- #

test_that("a prior INLA is not given says flexyBayes does not translate it", {
  skip_if_not_installed("INLA")
  d <- .rw_data()
  err <- .rw_refusal(y ~ x, random = ~gen, data = d, backend = "inla",
    prior = fb_prior(b("x") ~ student_t(3, 0, 1)))
  expect_s3_class(err, "flexybayes_refusal_prior_not_translatable_for_backend")
  msg <- conditionMessage(err)
  expect_match(msg, "flexyBayes does not translate 1 of the supplied prior",
    fixed = TRUE)
  expect_false(grepl("cannot carry", msg, fixed = TRUE))
  expect_false(grepl("0.9.2", msg, fixed = TRUE))
})

test_that("the gate's hurdle-gamma refusal reads INLA's own likelihood list", {
  skip_if_not_installed("INLA")
  d <- .rw_data()
  set.seed(24L)
  d$hg <- ifelse(stats::rbinom(nrow(d), 1L, 0.7) == 1L,
    stats::rgamma(nrow(d), shape = 2, rate = 0.5), 0)
  err <- .rw_refusal(hg ~ env, random = ~gen, data = d,
    family = "hurdle_gamma", backend = "inla")
  expect_s3_class(err, "flexybayes_refusal_inla_gate_refused")
  msg <- conditionMessage(err)
  expect_match(msg, "flexyBayes cannot fit this model on INLA", fixed = TRUE)
  expect_match(msg, "is not among INLA's likelihoods", fixed = TRUE)
  expect_false(grepl("Laplace machinery", msg, fixed = TRUE))
  # brms fits a hurdle gamma, so the refusal points there.
  expect_match(msg, "Fit it with backend = \"brms\"", fixed = TRUE)
})

test_that("the gate's link refusal describes flexyBayes's INLA route", {
  skip_if_not_installed("INLA")
  d <- .rw_data()
  d$bin <- stats::rbinom(nrow(d), 1L, 0.5)
  err <- .rw_refusal(bin ~ env, random = ~gen, data = d,
    family = stats::binomial(link = "probit"), backend = "inla")
  expect_s3_class(err, "flexybayes_refusal_inla_gate_refused")
  msg <- conditionMessage(err)
  expect_match(msg, "flexyBayes's INLA route passes the family's default link",
    fixed = TRUE)
  # INLA itself takes a probit link (test-engine-claims.R).
  expect_false(grepl("INLA takes the family's default link", msg, fixed = TRUE))
})

test_that("the gate does not point at brms for a model brms refuses too", {
  skip_if_not_installed("INLA")
  d <- .rw_data()
  for (spell in list(~ at(env, "1"):gen, ~ s(x))) {
    err <- .rw_refusal(y ~ env, random = spell, data = d, backend = "inla")
    expect_s3_class(err, "flexybayes_refusal_inla_gate_refused")
    msg <- conditionMessage(err)
    label <- deparse(spell)
    expect_match(msg, "flexyBayes has no INLA lowering for this term",
      fixed = TRUE, label = label)
    expect_match(msg, "flexyBayes cannot fit it on brms either", fixed = TRUE,
      label = label)
    expect_false(grepl("Fit it with backend = \"brms\"", msg, fixed = TRUE),
      label = label)
    expect_false(grepl("flexyBayes gap, not an INLA one", msg, fixed = TRUE),
      label = label)
    .rw_expect_plain(msg, label)
  }
  residual <- .rw_refusal(y ~ env, random = ~gen,
    residual = ~ at(env, "1"):units, data = d, backend = "inla")
  expect_s3_class(residual, "flexybayes_refusal_inla_gate_refused")
  expect_match(conditionMessage(residual),
    "flexyBayes has no INLA lowering for this residual term", fixed = TRUE)
})

test_that("the gate's us() refusals make no claim about INLA's models", {
  skip_if_not_installed("INLA")
  set.seed(25L)
  d <- data.frame(
    y = stats::rnorm(60L),
    env = factor(rep(letters[1:6], each = 10L)),
    gen = factor(rep(seq_len(10L), times = 6L))
  )
  err <- .rw_refusal(y ~ 1, random = ~ us(env):gen, data = d,
    backend = "inla")
  expect_s3_class(err, "flexybayes_refusal_correlated_covariance_requires_brms")
  msg <- conditionMessage(err)
  # INLA's iidkd carries more than five levels; the cap is flexyBayes's.
  expect_match(msg, "flexyBayes caps it at 5 levels on INLA", fixed = TRUE)
  expect_false(grepl("iid2d..iid5d", msg, fixed = TRUE))
  # fa() is refused on both engines, so it is no longer the suggestion.
  expect_false(grepl("Use fa(", msg, fixed = TRUE))
  expect_match(msg, "is above the 15 flexyBayes allows on INLA", fixed = TRUE)
  expect_false(grepl("CCD", msg, fixed = TRUE))
})

# --- the aggregated route, the prior notice, the selection certificate --- #

test_that("an aggregated sectioned residual is refused as flexyBayes's gap", {
  skip_if_not_installed("INLA")
  set.seed(13L)
  d <- expand.grid(
    rep = factor(seq_len(20L)),
    gen = factor(seq_len(6L)),
    env = factor(seq_len(3L))
  )
  d$y <- stats::rnorm(nrow(d), 10, 2)
  err <- .rw_refusal(y ~ env, random = ~gen, residual = ~ at(env):units,
    data = d, backend = "inla", aggregate = TRUE, return_code = FALSE)
  expect_s3_class(err, "flexybayes_aggregate_refusal")
  msg <- conditionMessage(err)
  expect_match(msg, "flexyBayes's aggregated route has no sectioned residual",
    fixed = TRUE)
  expect_false(grepl("multi-likelihood stack", msg, fixed = TRUE))
})

test_that("the default-prior notice is two sentences with the scale and the off switch", {
  flexyBayes:::.reset_emit_state_for_test()
  d <- .rw_data()
  msgs <- character(0)
  withr::with_options(
    list(flexyBayes.silence_default_prior_note = FALSE),
    withCallingHandlers(
      invisible(flexybayes(y ~ env, random = ~gen, data = d, plan = TRUE,
        verbose = FALSE)),
      message = function(m) {
        msgs[[length(msgs) + 1L]] <<- conditionMessage(m)
        invokeRestart("muffleMessage")
      }
    )
  )
  note <- msgs[grepl("half_normal(scale = ", msgs, fixed = TRUE)]
  expect_length(note, 1L)
  # The scale is 5 sd(y) on an identity link.
  expect_match(note, paste0("half_normal(scale = ",
    flexyBayes:::.fb_signif_string(5 * stats::sd(d$y)), ")"), fixed = TRUE)
  expect_match(note, "?fb_prior", fixed = TRUE)
  expect_match(note, "options(flexyBayes.silence_default_prior_note = TRUE)",
    fixed = TRUE)
  # Two sentences: two full stops that end a sentence. The notice used to
  # be a six-line paragraph.
  body <- sub("^flexyBayes: ", "", trimws(note))
  expect_identical(lengths(regmatches(body, gregexpr("\\.( |$)", body))), 2L)
  expect_lte(nchar(note), 320L)
})

test_that("the selection certificate says the shrinkage prior is fitted on brms by flexyBayes's choice", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  d <- fb_sel_bayes_data()
  base <- suppressMessages(flexybayes(
    y ~ p1, random = ~site, data = d, backend = "inla", verbose = FALSE
  ))
  plan <- fb_select(
    base, ~ x1 + x2 + z1 + z2,
    penalty = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    resamples = fb_resamples(unit = ~site, v = 3L),
    decision = fb_decision(method = "full"),
    plan = TRUE, seed = 1L
  )
  approx <- plan$selection$approximations
  hit <- approx[grepl("reference_fitted_on_brms", approx, fixed = TRUE)]
  expect_length(hit, 1L)
  expect_match(hit, "flexyBayes fits a joint shrinkage prior on brms only",
    fixed = TRUE)
  expect_false(grepl("INLA does not fit", hit, fixed = TRUE))
})
