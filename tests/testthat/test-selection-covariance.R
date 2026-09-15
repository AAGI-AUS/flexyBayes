# =============================================================================
# Selection covariance: oracles 1 and 2 of the correctness ledger.
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md sections 4.3 and 12.
#
# Everything fb_select() reports is computed on data whitened by the operator
# these tests gate. If V is assembled wrongly, every increment, every entry
# rank and every variance share downstream is wrong in a way no amount of
# resampling will reveal, because the error is common to every fold.
#
#   Oracle 1  the builder's implied V equals a dense reference assembled the
#             obvious way, elementwise.
#   Oracle 2  W V W' = I.
#
# Both oracles ship with the falsifier that makes them red. A gate nobody has
# seen fail is not evidence -- audit A25 and A27 are two gates that passed for
# every input because nobody had ever made them fail on purpose, and this file
# exists partly so that class does not recur here.
#
# The dense reference is deliberately naive. Its job is to be obviously right
# rather than fast, and it refuses above a size cap instead of exhausting
# memory: the deposit's V at 73,341 observations is about 43 GB.
# =============================================================================


# --- shared fixture -------------------------------------------------- #

# john.alpha is the nested alpha design the ASReml oracle uses for scenario 02
# and it is the shape the nert barley arm has: replicate, and block within
# replicate. One fit, reused, because each costs a few seconds.
fb_sel_fixture <- function() {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  d <- agridat::john.alpha
  d$gen <- factor(d$gen)
  d$rep <- factor(d$rep)
  d$block <- factor(d$block)
  flexybayes(
    fixed = yield ~ gen, random = ~ rep + rep:block,
    data = d, backend = "inla", verbose = FALSE,
    prior = .fb_test_john_alpha_prior(d)
  )
}


# --- the estimand -------------------------------------------------- #

test_that("plugin_mean_covariance is E[sigma^2], not the square of E[sigma]", {
  fit <- fb_sel_fixture()
  vc <- summary(fit)$varcomp

  plug <- .fb_sel_variance_components(fit, "plugin_mean_covariance")
  point <- .fb_sel_variance_components(fit, "point_estimate_covariance")

  # E[sigma^2] = mean^2 + sd^2 by definition, so the plugin mean is the point
  # estimate plus the posterior variance of the component. It is never smaller.
  expect_true(all(plug$terms >= point$terms))
  expect_gte(plug$residual, point$residual)

  # And the gap is the posterior variance, exactly -- not approximately.
  sd_rows <- startsWith(vc$component, "sd_")
  expect_equal(
    unname(plug$terms - point$terms),
    unname(vc$std.error[sd_rows]^2),
    tolerance = 1e-12
  )

  # The gap is material where a component is weakly identified. sd_rep here
  # rests on three levels, and squaring its mean would understate the variance
  # by more than half. This is the case D1 exists for.
  rep_ratio <- plug$terms[["rep"]] / point$terms[["rep"]]
  expect_gt(rep_ratio, 1.5)
})


# --- block structure ------------------------------------------------- #

test_that("the block structure of a nested design is found, and named", {
  fit <- fb_sel_fixture()
  bs <- .fb_sel_block_structure(fit)

  expect_true(bs$applicable)
  expect_identical(bs$factor, "rep")
  expect_identical(bs$n_blocks, 3L)
  expect_identical(sum(lengths(bs$blocks)), nrow(fit$data))
  # Blocks partition the rows: no row in two blocks, none missing.
  expect_setequal(unlist(bs$blocks), seq_len(nrow(fit$data)))
})


test_that("a crossed design is refused rather than given a coarse partition", {
  # Two crossed factors with a third spanning both: there is no grouping
  # factor that contains every random term, which is the deposit's shape
  # (site crossed with year, variety across everything).
  d <- expand.grid(site = factor(1:4), yr = factor(1:3))
  d$variety <- factor(rep(1:2, length.out = nrow(d)))
  d$y <- rnorm(nrow(d))

  mock <- list(
    data = d,
    fb = list(
      random_terms = list(
        list(type = "simple", var = "site"),
        list(type = "simple", var = "yr")
      ),
      residual_terms = list(list(type = "units"))
    )
  )

  bs <- .fb_sel_block_structure(mock)
  expect_false(bs$applicable)
  expect_match(bs$reason, "crossed")
})


test_that("a non-units residual term is refused by the block builder", {
  d <- data.frame(g = factor(rep(1:3, each = 4)), y = rnorm(12))
  mock <- list(
    data = d,
    fb = list(
      random_terms = list(list(type = "simple", var = "g")),
      residual_terms = list(list(type = "ar1"))
    )
  )
  bs <- .fb_sel_block_structure(mock)
  expect_false(bs$applicable)
  expect_match(bs$reason, "independent units")
})


# --- ORACLE 1: V is assembled correctly ------------------------------ #

test_that("oracle 1: the builder's implied V equals the dense reference", {
  fit <- fb_sel_fixture()

  V_ref <- .fb_sel_V_dense(fit, "plugin_mean_covariance")
  w <- .fb_sel_whitener(fit, "plugin_mean_covariance")

  # Reassemble V from what the builder actually holds: V_b = L_b L_b'.
  V_built <- matrix(0, nrow = w$n, ncol = w$n)
  for (f in w$factors) {
    V_built[f$rows, f$rows] <- tcrossprod(f$L)
  }

  expect_equal(V_built, V_ref, tolerance = 1e-10)
})


test_that("oracle 1 falsifier: a perturbed component makes the check red", {
  fit <- fb_sel_fixture()
  V_ref <- .fb_sel_V_dense(fit, "plugin_mean_covariance")

  # A one per cent error in one variance component. If the comparison cannot
  # see this, it is not comparing anything.
  vc <- .fb_sel_variance_components(fit, "plugin_mean_covariance")
  Z <- .fb_sel_random_designs(fit)
  V_wrong <- diag(vc$residual, nrow = nrow(fit$data))
  bump <- 1.01
  for (nm in names(Z)) {
    mult <- if (nm == "rep:block") bump else 1
    V_wrong <- V_wrong + mult * vc$terms[[nm]] * tcrossprod(Z[[nm]])
  }

  expect_false(isTRUE(all.equal(V_wrong, V_ref, tolerance = 1e-10)))
})


test_that("oracle 1 falsifier: the two estimands give different V", {
  # If the estimand label were decorative, these would coincide. They must not.
  fit <- fb_sel_fixture()
  V_plug <- .fb_sel_V_dense(fit, "plugin_mean_covariance")
  V_point <- .fb_sel_V_dense(fit, "point_estimate_covariance")
  expect_false(isTRUE(all.equal(V_plug, V_point, tolerance = 1e-8)))
})


# --- ORACLE 2: the whitening identity -------------------------------- #

test_that("oracle 2: W V W' = I", {
  fit <- fb_sel_fixture()
  w <- .fb_sel_whitener(fit, "plugin_mean_covariance")
  V <- .fb_sel_V_dense(fit, "plugin_mean_covariance")

  WV <- .fb_sel_whiten(w, V)
  WVWt <- t(.fb_sel_whiten(w, t(WV)))

  expect_equal(WVWt, diag(1, w$n), tolerance = 1e-9)
})


test_that("oracle 2 falsifier: a residual-only whitener fails the identity", {
  # The most plausible wrong whitener is one that forgets the random effects
  # and divides by sigma. It must not pass.
  fit <- fb_sel_fixture()
  V <- .fb_sel_V_dense(fit, "plugin_mean_covariance")
  vc <- .fb_sel_variance_components(fit, "plugin_mean_covariance")

  W_wrong <- diag(1 / sqrt(vc$residual), nrow = nrow(V))
  WVWt_wrong <- W_wrong %*% V %*% t(W_wrong)

  expect_false(isTRUE(all.equal(WVWt_wrong, diag(1, nrow(V)), tolerance = 1e-9)))
})


test_that("whitening is a forward solve: no inverse is ever formed", {
  fit <- fb_sel_fixture()
  w <- .fb_sel_whitener(fit, "plugin_mean_covariance")

  # Structural, not aspirational: the object holds lower-triangular factors
  # and nothing of size n by n.
  for (f in w$factors) {
    expect_true(all(f$L[upper.tri(f$L)] == 0))
    expect_identical(nrow(f$L), length(f$rows))
  }
  expect_false(any(vapply(
    w$factors, function(f) nrow(f$L) == w$n && w$n_blocks > 1L, logical(1)
  )))
})


test_that("whitening a vector matches whitening it as a one-column matrix", {
  fit <- fb_sel_fixture()
  w <- .fb_sel_whitener(fit, "plugin_mean_covariance")
  y <- fit$data$yield
  expect_equal(
    .fb_sel_whiten(w, y),
    .fb_sel_whiten(w, matrix(y, ncol = 1L)),
    tolerance = 1e-12
  )
})


test_that("a whitener applied to the wrong number of rows refuses", {
  fit <- fb_sel_fixture()
  w <- .fb_sel_whitener(fit, "plugin_mean_covariance")
  expect_error(
    .fb_sel_whiten(w, matrix(0, nrow = w$n + 1L, ncol = 1L)),
    "was built for"
  )
})


# --- refusals --------------------------------------------------------- #

test_that("the dense reference refuses rather than exhaust memory", {
  fit <- fb_sel_fixture()
  expect_error(
    .fb_sel_V_dense(fit, "plugin_mean_covariance", max_n = 10L),
    "dense covariance reference"
  )
})


test_that("Z and theta are matched by name, not by position", {
  fit <- fb_sel_fixture()
  Z <- .fb_sel_random_designs(fit)
  vc <- .fb_sel_variance_components(fit, "plugin_mean_covariance")

  # The names the two sides use must agree, or a variance is silently
  # attached to the wrong structure and every downstream number is wrong
  # while every test that only checks dimensions still passes.
  expect_setequal(names(Z), names(vc$terms))
  expect_identical(sort(names(Z)), c("rep", "rep:block"))
})


test_that("the crossed base refuses at the whitener, with the named reason", {
  d <- expand.grid(site = factor(1:4), yr = factor(1:3))
  d$y <- rnorm(nrow(d))
  mock <- list(
    data = d,
    fb = list(
      random_terms = list(
        list(type = "simple", var = "site"),
        list(type = "simple", var = "yr")
      ),
      residual_terms = list(list(type = "units"))
    )
  )
  # The refusal must precede any use of the fit's summary, so a mock reaches
  # it. It must also name the route that does apply, or the message tells the
  # user they are stuck when they are not.
  expect_error(.fb_sel_whitener(mock), "random-term reduction")
})


test_that("every refusal this module can throw is registered", {
  # The package raises a contract violation when an unregistered code is
  # thrown, which means an unregistered code surfaces as the wrong error and
  # a bare expect_error() passes on it. This test is the direct check.
  codes <- c(
    "selection_covariance_not_buildable",
    "selection_covariance_not_positive_definite",
    "selection_dense_reference_too_large"
  )
  registered <- fb_refusals()$reason_code
  expect_true(all(codes %in% registered))
})
