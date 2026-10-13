# Nested / interaction random effects on INLA.
#
# These refused with a message asserting that INLA's f() machinery cannot
# represent the term class without an SPDE or kronecker expansion. It can:
# the term is an iid field over the interaction index, which is what brms's
# `(1 | a:b)` means too. agridat::john.alpha, the alpha-design in the report
# that raised this (issue #14, fit2), fits in seconds.

test_that("the gate accepts a nested random term for INLA", {
  tab <- data.frame(
    y = rnorm(24),
    a = factor(rep(letters[1:2], each = 12L)),
    b = factor(rep(letters[3:4], times = 12L)),
    g = factor(rep(paste0("g", 1:4), each = 6L))
  )
  fb <- fb_from_asreml(
    fixed = y ~ a,
    random = ~ a:b + g,
    data = tab,
    family = "gaussian"
  )
  # lgm_gate() returns the IR, marking a pass by setting capabilities.
  res <- lgm_gate(fb)
  expect_true("lgm_compatible" %in% res$capabilities)
})

test_that("the emitted INLA formula carries the interaction index", {
  tab <- data.frame(
    y = rnorm(24),
    a = factor(rep(letters[1:2], each = 12L)),
    b = factor(rep(letters[3:4], times = 12L)),
    g = factor(rep(paste0("g", 1:4), each = 6L))
  )
  fb <- fb_from_asreml(
    fixed = y ~ a,
    random = ~ a:b + g,
    data = tab,
    family = "gaussian"
  )
  form <- .build_inla_formula(fb)
  txt <- paste(deparse(form), collapse = " ")
  expect_match(txt, "f(a_b_id, model = \"iid\"", fixed = TRUE)
})

test_that("a three-way interaction emits one index over all three factors", {
  fb <- list(random_terms = list(list(type = "combo", vars = c("x", "y", "z"))))
  expect_identical(.fb_nested_vars(fb$random_terms[[1L]]), c("x", "y", "z"))
  expect_identical(.fb_inla_interaction_index(c("x", "y", "z")), "x_y_z_id")
})

test_that("a generated index does not silently overwrite a user column", {
  # The failure this prevents is silent: the field would be indexed by
  # whatever the user's column holds, and nothing downstream would say so.
  # emit_inla() checks that INLA is installed before it reads the data.
  skip_if_not_installed("INLA")
  tab <- data.frame(
    y = rnorm(24),
    a = factor(rep(letters[1:2], each = 12L)),
    b = factor(rep(letters[3:4], times = 12L)),
    a_b_id = seq_len(24L)
  )
  fb <- fb_from_asreml(
    fixed = y ~ 1,
    random = ~ a:b,
    data = tab,
    family = "gaussian"
  )
  expect_error(emit_inla(fb, tab, verbose = FALSE), "already has one")
})

test_that("the interaction index does not leak into the summary label", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_if_not_installed("agridat")
  fit <- flexybayes(
    yield ~ rep,
    random = ~ rep:block + gen,
    data = agridat::john.alpha,
    backend = "inla"
  )
  comps <- summary(fit)$varcomp$component
  expect_true("sd_rep:block" %in% comps)
  expect_false(any(grepl("_id", comps, fixed = TRUE)))
})

test_that("INLA recovers the same variance components as lme4 REML", {
  # The independent oracle. lme4 fits `(1 | rep:block) + (1 | gen)` by REML
  # on the same data, and never sees flexyBayes's parser or emit. Posterior
  # means of variance components sit a little above REML point estimates
  # because the posteriors are right-skewed, so the tolerance is on the
  # ratio, not the difference.
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_if_not_installed("agridat")
  skip_on_cran()

  fit <- flexybayes(
    yield ~ rep,
    random = ~ rep:block + gen,
    data = agridat::john.alpha,
    backend = "inla"
  )
  v <- summary(fit)$varcomp
  got <- stats::setNames(v$estimate, v$component)

  m <- lme4::lmer(
    yield ~ rep + (1 | rep:block) + (1 | gen),
    data = agridat::john.alpha,
    REML = TRUE
  )
  vc <- as.data.frame(lme4::VarCorr(m))
  ref <- stats::setNames(vc$sdcor, vc$grp)

  expect_equal(unname(got[["sigma"]]), unname(ref[["Residual"]]),
               tolerance = 0.15)
  expect_equal(unname(got[["sd_gen"]]), unname(ref[["gen"]]),
               tolerance = 0.15)
  expect_equal(unname(got[["sd_rep:block"]]), unname(ref[["rep:block"]]),
               tolerance = 0.20)
})
