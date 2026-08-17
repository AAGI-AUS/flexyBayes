# =============================================================================
# Round trip: every structure a user can write that no active engine fits
# comes back as a NAMED refusal.
#
# The August 2026 external audit found two ways out of the formula surface
# that were not named refusals. An unrecognised call -- ar2(row),
# corgh(site):gen -- was read as a plain variable whose name happened to be
# its own source text, so a model nobody could fit was silently replaced by
# a model nobody asked for. And a term type with no branch in the brms
# formula reconstruction hit a bare stop(), outside the refusal vocabulary
# and uncatchable by class.
#
# The contract asserted here is one sentence: from the formula surface, no
# untyped stop is reachable. Each case below is written the way a user would
# write it, dispatched through the real gate and emit path, and required to
# raise a condition carrying its reason code.
# =============================================================================

.rt_data <- function(seed = 21L) {
  set.seed(seed)
  d <- expand.grid(
    rep = factor(seq_len(2L)),
    gen = factor(seq_len(6L)),
    env = factor(seq_len(3L))
  )
  d$site <- factor(rep(seq_len(3L), length.out = nrow(d)))
  d$trait <- factor(rep(seq_len(2L), length.out = nrow(d)))
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rnorm(nrow(d), 10, 2)
  d
}

.rt_grid <- function(n_row = 6L, n_col = 5L, seed = 22L) {
  set.seed(seed)
  g <- expand.grid(
    row = factor(seq_len(n_row)),
    col = factor(seq_len(n_col))
  )
  g$y <- stats::rnorm(nrow(g))
  g
}

# Walk the whole gate + emit path without sampling.
.rt_emit <- function(..., backend = "auto", known_matrices = list()) {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE
  )
  tryCatch(
    suppressMessages(flexybayes(
      ...,
      known_matrices = known_matrices,
      backend = backend,
      return_code = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
}

# --- the closed parser vocabulary ------------------------------------- #

test_that("an unrecognised call is refused, not read as a variable", {
  d <- .rt_data()
  bad <- .rt_emit(y ~ 1, random = ~ foo(x), data = d)
  expect_s3_class(bad, "flexybayes_refusal_asreml_function_not_recognised")
  expect_match(conditionMessage(bad), "foo(", fixed = TRUE)
  # The pre-0.9.0 behaviour: the call became a variable of that name. If the
  # message ever stops naming the token, that reading is back.
  expect_match(conditionMessage(bad), "closed set", fixed = TRUE)
})

test_that("a near-miss spelling is pointed at the implemented one", {
  d <- .rt_data()
  bad <- .rt_emit(y ~ 1, random = ~ corgh(site):gen, data = d)
  expect_s3_class(bad, "flexybayes_refusal_asreml_function_not_recognised")
  # corgh is one edit from corh and two from cor; the suggestion must reach
  # a spelling the parser reads, so the next refusal names an alternative.
  expect_match(conditionMessage(bad), "nearest spelling")
  expect_match(conditionMessage(bad), "corh()", fixed = TRUE)
})

test_that("every name in the vocabulary is recognised by the walker", {
  # The vocabulary vector and the switch() it documents must not drift: a
  # name listed as recognised that the switch does not answer to would be
  # refused as unknown, and the refusal message would then be lying about
  # its own contents.
  d <- .rt_data()
  d$t <- factor(rep(seq_len(6L), length.out = nrow(d)))
  calls <- list(
    ar1 = quote(ar1(t)), ar2 = quote(ar2(t)), at = quote(at(env)),
    cor = quote(cor(env)), corh = quote(corh(env)), diag = quote(diag(env)),
    dsum = quote(dsum(~units | env)), fa = quote(fa(env, 2)),
    id = quote(id(gen)), ide = quote(ide(gen)), idh = quote(idh(env)),
    lin = quote(lin(x)), ped = quote(ped(gen)), pol = quote(pol(x)),
    s = quote(s(x)), spl = quote(spl(x)), str = quote(str(~gen, ~us(2))),
    t2 = quote(t2(x, y)), te = quote(te(x, y)), ti = quote(ti(x, y)),
    us = quote(us(env)), vm = quote(vm(gen))
  )
  expect_setequal(names(calls), flexyBayes:::.asreml_function_vocabulary())
  for (nm in names(calls)) {
    out <- tryCatch(flexyBayes:::.walk(calls[[nm]]), error = function(e) e)
    expect_false(
      inherits(out, "flexybayes_refusal_asreml_function_not_recognised"),
      label = paste0(nm, "() is in the vocabulary")
    )
  }
})

# --- parsed for the catalogue, refused by name ------------------------ #

test_that("ar2() is refused and points at ar1()", {
  g <- .rt_grid()
  bad <- .rt_emit(y ~ 1, random = ~ ar2(row), data = g)
  expect_s3_class(bad, "flexybayes_refusal_ar2_not_representable")
  expect_match(conditionMessage(bad), "ar1(row)", fixed = TRUE)
})

test_that("cor() is refused with the equicorrelation reason", {
  d <- .rt_data()
  bad <- .rt_emit(y ~ 1, random = ~ cor(env), data = d)
  expect_s3_class(
    bad, "flexybayes_refusal_corh_no_equicorrelation_representation"
  )
})

test_that("str() is refused and names what to write instead", {
  d <- .rt_data()
  bad <- .rt_emit(y ~ 1, random = ~ str(~gen, ~us(2):id(gen)), data = d)
  expect_s3_class(bad, "flexybayes_refusal_str_not_representable")
  expect_match(conditionMessage(bad), "us(f):g", fixed = TRUE)
})

test_that("a factor-analytic covariance is refused on both engines", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  d <- .rt_data()
  for (engine in c("inla", "brms")) {
    bad <- .rt_emit(y ~ env, random = ~ fa(env, 1):gen, data = d,
                    backend = engine)
    expect_s3_class(bad, "flexybayes_refusal_fa_not_representable")
    expect_match(conditionMessage(bad), "us(env):gen", fixed = TRUE)
  }
  # It still parses, which is what makes fb_terms() a formula catalogue.
  fb <- fb_from_asreml(y ~ env, random = ~ fa(env, 1):gen, data = d)
  expect_identical(fb$random_terms[[1L]]$type, "fa_gxe")
})

test_that("an unmatched structured interaction is refused by name", {
  # us(trait):vm(gen) -- the multi-trait idiom -- matches no interaction
  # pattern, so it lands on the catalogue's generic node. That node has no
  # emit anywhere and must not reach one.
  d <- .rt_data()
  k <- diag(6)
  dimnames(k) <- list(as.character(seq_len(6L)), as.character(seq_len(6L)))
  bad <- .rt_emit(
    y ~ env, random = ~ us(trait):vm(gen), data = d,
    known_matrices = list(K = k)
  )
  expect_s3_class(bad, "flexybayes_refusal_interaction_not_representable")
})

# --- the brms formula reconstruction ---------------------------------- #

test_that("a term brms cannot lower is refused by name, not dropped", {
  skip_if_not_installed("brms")
  d <- .rt_data()
  bad <- .rt_emit(y ~ 1, random = ~ spl(x), data = d, backend = "brms")
  expect_s3_class(bad, "flexybayes_refusal_brms_cannot_represent_term")
  expect_match(conditionMessage(bad), "\"spline\"", fixed = TRUE)
  expect_match(conditionMessage(bad), "backend = \"inla\"", fixed = TRUE)
})

# --- the AR1 respelling ----------------------------------------------- #

test_that("the residual AR1 spelling refuses and the random side fits", {
  skip_if_not_installed("INLA")
  g <- .rt_grid()
  bad <- .rt_emit(y ~ 1, residual = ~ ar1(row):ar1(col), data = g,
                  backend = "inla")
  expect_s3_class(bad, "flexybayes_refusal_ar1_residual_not_representable")
  ok <- .rt_emit(y ~ 1, random = ~ ar1(row):ar1(col), data = g,
                 backend = "inla")
  expect_false(inherits(ok, "error"))
})

test_that("the random-side field is refused on brms", {
  skip_if_not_installed("brms")
  g <- .rt_grid()
  bad <- .rt_emit(y ~ 1, random = ~ ar1(row):ar1(col), data = g,
                  backend = "brms")
  expect_s3_class(bad, "flexybayes_refusal_stan_cannot_represent_ar1_field")
})

# --- a numeric variable inside a random interaction ------------------- #

test_that("Subject:Days with numeric Days names (Days || Subject)", {
  skip_if_not_installed("lme4")
  data(sleepstudy, package = "lme4")
  bad <- .rt_emit(
    Reaction ~ Days, random = ~ Subject + Subject:Days, data = sleepstudy
  )
  expect_s3_class(
    bad, "flexybayes_refusal_numeric_variable_in_random_interaction"
  )
  msg <- conditionMessage(bad)
  expect_match(msg, "(Days || Subject)", fixed = TRUE)
  expect_match(msg, "factor(data$Days)", fixed = TRUE)
  # The count in the message is the size of the model the crossing would
  # actually build: 18 subjects x 10 days.
  expect_match(msg, "180", fixed = TRUE)
})

test_that("a factor-by-factor nested term still fits", {
  # The guard must be about the variable's type, not about the spelling. A
  # legitimate crossing of two factors is the multi-stratum designed
  # experiment brms is the faithful backend for.
  skip_if_not_installed("brms")
  d <- .rt_data()
  ok <- .rt_emit(y ~ env, random = ~ gen + gen:env, data = d,
                 backend = "brms")
  expect_false(inherits(ok, "error"))
  # And a three-way combination of factors likewise.
  ok3 <- .rt_emit(y ~ env, random = ~ gen:env:rep, data = d,
                  backend = "brms")
  expect_false(inherits(ok3, "error"))
})

test_that("a numeric inside a three-way random combination refuses too", {
  d <- .rt_data()
  bad <- .rt_emit(y ~ 1, random = ~ gen:env:x, data = d)
  expect_s3_class(
    bad, "flexybayes_refusal_numeric_variable_in_random_interaction"
  )
})
