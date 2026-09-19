# =============================================================================
# The reference fit: the posterior a penalised selection decides from.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md section 5.
#
# Four things are under test and they are separable on purpose.
#
#   the re-emission   the candidate block reaches the Stan program, on both
#                     of its routes (a class-wide row, and one row per
#                     coefficient), and the compiled prior REPLACES the
#                     default row on the class it names rather than joining
#                     it -- brms reads two rows for one parameter as a
#                     specification error
#   the seed          required before anything samples
#   the two gates     each shown red on a deliberately broken fit: an
#                     under-warmed sampler for the convergence gate, and a
#                     reference model whose residual scale is pinned at a
#                     twentieth of the truth for the posterior predictive
#                     gate. A gate nobody has seen fail is not a gate
#   kappa             the horseshoe's shrinkage factors, computed from the
#                     scale parameter brms actually saves, and carrying the
#                     label that says where the scalar formula holds
# =============================================================================


# --- the fits these tests share ---------------------------------------- #

fb_sel_ref_bits <- function(n_samples = 2500L, warmup = 1000L, ...) {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  d <- fb_sel_bayes_data()
  cand <- c("x1", "x2", "z1", "z2")
  comp <- .fb_sel_test_blocks(cand, family = "regularized_horseshoe", p0 = 1)
  base <- fb_sel_bayes_memo("base", fb_sel_bayes_base(d))
  list(d = d, cand = cand, comp = comp, base = base)
}

fb_sel_reference <- function(bits) {
  # A distinct key: test-selection-projection.R memoises its own reference
  # fit at a smaller sampler budget into the same shared cache, and the
  # adequacy gates below are assertions about THIS fit's diagnostics.
  fb_sel_bayes_memo(
    "reference_adequacy",
    .fb_sel_reference_fit(
      bits$base,
      ~ x1 + x2 + z1 + z2,
      prior = bits$comp$prior,
      stanvars = bits$comp$stanvars,
      column_map = .fb_sel_test_column_map(bits$cand),
      data = bits$d,
      engine = fb_engine("brms", chains = 2L, n_samples = 2500L, warmup = 1000L),
      seed = 20260912L
    )
  )
}


# --- the re-emission ---------------------------------------------------- #

test_that("the reference fit is the base re-emitted with the candidate block", {
  skip_if_not_long_tier()
  bits <- fb_sel_ref_bits()
  ref <- fb_sel_reference(bits)

  expect_s3_class(ref, "flexybayes_brms")
  # The slots the selection substrate reads. A brms fit does not carry them
  # by default, and half the module cannot read a fit that lacks them.
  expect_false(is.null(ref$fb))
  expect_true(is.data.frame(ref$data))
  expect_false(is.null(ref$extras$fb_terms))

  expect_identical(attr(ref, "candidate_terms"), bits$cand)
  expect_identical(attr(ref, "protected_terms"), "p1")
  expect_identical(attr(ref, "random_terms"), "(1 | site)")
  expect_true(is.data.frame(attr(ref, "column_map")))
  expect_identical(attr(ref, "reference_seed"), 20260912L)

  # Every candidate carries a coefficient in the posterior, and the
  # protected one is still there beside them.
  vars <- posterior::variables(ref$brms)
  expect_true(all(paste0("b_", bits$cand) %in% vars))
  expect_true("b_p1" %in% vars)

  # The compiled block replaced the default `b` row rather than joining it.
  prior <- as.data.frame(ref$extras$brms_prior)
  b_rows <- prior[prior$class == "b" & !nzchar(as.character(prior$coef)), ]
  expect_identical(nrow(b_rows), 1L)
  expect_match(b_rows$prior[[1L]], "horseshoe")
})


test_that("the coefficient-row route puts one prior on each candidate", {
  # The other compilation route, checked on the generated program rather
  # than on a fit: `return_code = TRUE` authors the Stan code without
  # compiling or sampling, so this costs a formula walk.
  skip_if_not_installed("brms")
  d <- fb_sel_bayes_data()
  comp <- .fb_sel_test_blocks(c("x1", "x2", "z1", "z2"), family = "ridge", scale = 0.4)
  expect_identical(nrow(as.data.frame(comp$prior)), 4L)

  fb <- .build_ir_polymorphic(
    fixed = y ~ p1 + x1 + x2 + z1 + z2 + (1 | site),
    random = NULL, residual = NULL, data = d,
    family = "gaussian", link = NULL, weights = NULL,
    known_matrices = list(), prior = NULL,
    prior_fixed_sd = 100, prior_vc_sd = 1, syntax = "auto"
  )
  code <- emit_brms(
    fb = fb, data = d, chains = 1L, n_samples = 100L, warmup = 100L,
    seed = 1L, verbose = FALSE, mcmc_verbose = FALSE, return_code = TRUE,
    extra_prior = comp$prior
  )
  # One normal(0, 0.4) contribution per candidate coefficient, and the
  # protected coefficient keeps the model's own default.
  hits <- gregexpr("normal_lpdf(b[", code, fixed = TRUE)[[1L]]
  expect_true(length(hits) >= 4L)
  expect_match(code, "0.4", fixed = TRUE)
})


test_that("a stanvars block reaches the generated program", {
  skip_if_not_installed("brms")
  d <- fb_sel_bayes_data()
  fb <- .build_ir_polymorphic(
    fixed = y ~ p1 + x1 + (1 | site),
    random = NULL, residual = NULL, data = d,
    family = "gaussian", link = NULL, weights = NULL,
    known_matrices = list(), prior = NULL,
    prior_fixed_sd = 100, prior_vc_sd = 1, syntax = "auto"
  )
  code <- emit_brms(
    fb = fb, data = d, chains = 1L, n_samples = 100L, warmup = 100L,
    seed = 1L, verbose = FALSE, mcmc_verbose = FALSE, return_code = TRUE,
    stanvars = brms::stanvar(
      scode = "  // flexybayes joint block", block = "model", position = "end"
    ),
    extra_prior = brms::set_prior(
      "target += normal_lpdf(b[1] | 0, 0.25)", check = FALSE
    )
  )
  # The verbatim row has no class, no coef and no group, so the
  # reachability oracle has no triple to find; it must reach the program
  # regardless, which is why the merge happens after that check.
  expect_match(code, "flexybayes joint block", fixed = TRUE)
  expect_match(code, "normal_lpdf(b[1] | 0, 0.25)", fixed = TRUE)
})


# --- the seed, and a candidate the data does not carry ------------------ #

test_that("the reference fit refuses without a seed", {
  skip_if_not_installed("brms")
  base <- structure(
    list(extras = list(call_info = list(fixed = y ~ p1))),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
  expect_error(
    .fb_sel_reference_fit(base, ~x1, seed = NULL),
    class = "flexybayes_refusal_selection_seed_required"
  )
  # Not a seed either: a vector of them, or a missing value.
  expect_error(
    .fb_sel_reference_fit(base, ~x1, seed = c(1L, 2L)),
    class = "flexybayes_refusal_selection_seed_required"
  )
})


test_that("a candidate absent from the fit's data refuses by name", {
  # A brms fit now carries the analysis data at `$data`, as an INLA fit
  # always did, so a base fitted on a frame that HOLDS the candidates finds
  # them there and the refusal is not reachable from it. The refusal is
  # still the right answer when the data genuinely lacks the column, and
  # this is the base that lacks it: the same fit with its data narrowed to
  # the variables the base model used, which is the shape a fit reloaded
  # from a brms model frame has.
  bits <- fb_sel_ref_bits()
  narrowed <- bits$base
  narrowed$data <- bits$d[, c("y", "p1", "site"), drop = FALSE]

  expect_error(
    .fb_sel_reference_fit(narrowed, ~ x1 + x2, seed = 1L),
    class = "flexybayes_refusal_selection_base_shape_unreadable"
  )
  err <- tryCatch(
    .fb_sel_reference_fit(narrowed, ~ x1 + x2, seed = 1L),
    condition = function(e) e
  )
  expect_match(conditionMessage(err), "x1, x2")

  # And the negative control the emit fix earns: the unnarrowed base does
  # NOT refuse, because the candidates are in the data it carries.
  expect_true(all(c("x1", "x2") %in% names(bits$base$data)))
})


# --- the adequacy gates -------------------------------------------------- #

test_that("a converged reference model passes both adequacy gates", {
  bits <- fb_sel_ref_bits()
  ref <- fb_sel_reference(bits)
  adequacy <- .fb_sel_reference_adequacy(ref, metric = "elpd", ndraws = 300L, seed = 7L)

  expect_true(adequacy$verdict %in% c("pass", "unverified"))
  expect_false(identical(adequacy$verdict, "fail"))
  expect_true(adequacy$ppc$p_mean > 0.005 && adequacy$ppc$p_mean < 0.995)
  expect_true(adequacy$ppc$p_sd > 0.005 && adequacy$ppc$p_sd < 0.995)
  expect_identical(adequacy$metric, "elpd")
  # The NUTS evidence the summary slot does not carry, now on the verdict.
  expect_true(all(
    c("min_ess_tail", "min_ebfmi", "n_max_treedepth") %in%
      names(adequacy$evidence)
  ))
  expect_true(is.finite(adequacy$evidence$min_ebfmi))
  expect_silent(.fb_sel_require_adequate(ref, ndraws = 100L, seed = 7L))
})


test_that("the convergence gate is red on a deliberately under-warmed fit", {
  # The injected defect: ten warm-up iterations. Nothing else changes.
  bits <- fb_sel_ref_bits()
  broken <- suppressWarnings(.fb_sel_reference_fit(
    bits$base,
    ~ x1 + x2 + z1 + z2,
    prior = bits$comp$prior,
    data = bits$d,
    engine = fb_engine("brms", chains = 2L, n_samples = 40L, warmup = 10L),
    seed = 20260912L
  ))
  adequacy <- .fb_sel_reference_adequacy(broken, ndraws = 100L, seed = 7L)
  expect_identical(adequacy$verdict, "fail")
  expect_false(isTRUE(adequacy$convergence))
  expect_gt(length(adequacy$reasons), 0L)

  expect_error(
    .fb_sel_require_adequate(broken, ndraws = 100L, seed = 7L),
    class = "flexybayes_refusal_penalty_sampler_diagnostics_failed"
  )
})


test_that("the posterior predictive gate is red when the model cannot reproduce the spread", {
  skip_if_not_long_tier()
  # The injected defect: the residual and group scales are pinned at a
  # twentieth of the truth through the same prior pass-through the joint
  # block uses. The model then cannot generate data as spread out as the
  # data it was fitted to, and the observed SD lands in the far tail.
  bits <- fb_sel_ref_bits()
  pinned <- rbind(
    as.data.frame(brms::set_prior("constant(0.05)", class = "sigma")),
    as.data.frame(brms::set_prior("constant(0.05)", class = "sd", group = "site"))
  )
  pinned <- structure(
    pinned,
    class = class(brms::set_prior("constant(0.05)", class = "sigma"))
  )
  misfit <- suppressWarnings(.fb_sel_reference_fit(
    bits$base,
    ~ x1 + x2 + z1 + z2,
    prior = pinned,
    data = bits$d,
    engine = fb_engine("brms", chains = 2L, n_samples = 1000L, warmup = 500L),
    seed = 20260912L
  ))
  adequacy <- .fb_sel_reference_adequacy(misfit, ndraws = 400L, seed = 7L)
  expect_false(isTRUE(adequacy$ppc_ok))
  expect_identical(adequacy$verdict, "fail")
  expect_match(
    paste(adequacy$reasons, collapse = " | "),
    "standard deviation"
  )
  err <- tryCatch(
    .fb_sel_require_adequate(misfit, ndraws = 400L, seed = 7L),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal")
  # Which gate failed is part of the answer: the sampler explored this
  # posterior perfectly well, and the model is still the wrong model.
  expect_true(inherits(
    err,
    "flexybayes_refusal_penalty_reference_model_inadequate"
  ) || inherits(err, "flexybayes_refusal_penalty_sampler_diagnostics_failed"))
})


# --- shrinkage factors --------------------------------------------------- #

test_that("kappa separates the planted signals from the noise and says where it holds", {
  bits <- fb_sel_ref_bits()
  ref <- fb_sel_reference(bits)
  kappa <- .fb_sel_kappa(ref, bits$comp$blocks)

  expect_true(kappa$available)
  expect_identical(kappa$label, "orthogonal_design_kappa_approx")
  expect_identical(names(kappa$kappa), bits$cand)
  expect_true(all(kappa$kappa > 0 & kappa$kappa < 1))
  expect_identical(kappa$n, nrow(bits$d))
  expect_true(is.numeric(kappa$effective_nonzero_at_half))

  # The substantive check: the two planted signals are shrunk far less than
  # the two noise columns. Without it the test passes for any function
  # returning numbers in the unit interval.
  expect_lt(max(kappa$kappa[c("x1", "x2")]), min(kappa$kappa[c("z1", "z2")]))
})


test_that("kappa reports its absence rather than inventing a number", {
  bits <- fb_sel_ref_bits()
  ref <- fb_sel_reference(bits)
  # A ridge block saves no scale parameter, so there is nothing to compute
  # the shrinkage factor from and the result says so.
  ridge <- .fb_sel_test_blocks(bits$cand, family = "ridge")
  kappa <- .fb_sel_kappa(ref, ridge$blocks)
  expect_false(kappa$available)
  expect_identical(kappa$label, "orthogonal_design_kappa_approx")
  expect_true(is.character(kappa$reason))
})
