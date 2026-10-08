# Heterogeneous variance by factor level on INLA: diag(f):g, idh(f):g, at(f):g.
#
# One free standard deviation per level of the outer factor and no covariance
# between levels. INLA writes this as k independent iid fields, each on an
# index that is NA outside its level. This refused on INLA until 0.10.0.9000
# with no reason recorded beyond the allowlist; it is one of ASReml's two
# multi-environment-trial workhorses.

.hv_fixture <- function(seed = 21L, n_gen = 40L,
                        sd_by_env = c(0.4, 1.0, 2.0), n_rep = 4L,
                        sd_resid = 0.30) {
  set.seed(seed)
  n_env <- length(sd_by_env)
  tab <- expand.grid(
    gen = factor(seq_len(n_gen)),
    env = factor(seq_len(n_env)),
    rep = seq_len(n_rep)
  )
  mat_g <- vapply(sd_by_env, function(s) rnorm(n_gen, 0, s), numeric(n_gen))
  tab$y <- 5 + mat_g[cbind(as.integer(tab$gen), as.integer(tab$env))] +
    rnorm(nrow(tab), 0, sd_resid)
  tab
}

test_that("the gate accepts diag(f):g for INLA", {
  tab <- .hv_fixture(n_gen = 8L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ diag(env):gen,
    data = tab,
    family = "gaussian"
  )
  res <- lgm_gate(fb)
  expect_true("lgm_compatible" %in% res$capabilities)
})

test_that("the emitted formula carries one iid field per level", {
  tab <- .hv_fixture(n_gen = 8L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ diag(env):gen,
    data = tab,
    family = "gaussian"
  )
  # deparse() wraps long formulas mid-call, so collapse and normalise the
  # whitespace before matching -- otherwise the third f() splits and the
  # match fails on formatting rather than on content.
  txt <- gsub("\\s+", " ", paste(deparse(.build_inla_formula(fb)),
                                 collapse = " "))
  for (lev in 1:3) {
    expect_match(txt, paste0("f(gen_at_env_", lev, ", model = \"iid\""),
                 fixed = TRUE)
  }
})

test_that("INLA recovers per-level variances and agrees with lme4 REML", {
  # lme4 fits `(0 + env || gen)`, the same k-variances-no-covariance model,
  # by REML and never sees this package's parser or emit.
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)

  sd_by_env <- c(0.4, 1.0, 2.0)
  tab <- .hv_fixture(sd_by_env = sd_by_env)
  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~ diag(env):gen,
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  v <- summary(fit)$varcomp
  got <- stats::setNames(v$estimate, v$component)

  # Recovery of the generating values.
  for (i in seq_along(sd_by_env)) {
    nm <- paste0("sd_gen__env", i)
    expect_true(nm %in% names(got))
    expect_equal(unname(got[[nm]]), sd_by_env[[i]], tolerance = 0.25)
  }
  expect_equal(unname(got[["sigma"]]), 0.30, tolerance = 0.15)

  # Agreement with the independent REML fit.
  m <- lme4::lmer(y ~ env + (0 + env || gen), data = tab, REML = TRUE)
  vc <- as.data.frame(lme4::VarCorr(m))
  ref <- vc$sdcor[vc$grp == "gen" & vc$var1 %in% paste0("env", 1:3) &
                    is.na(vc$var2)]
  ref <- ref[seq_along(sd_by_env)]
  got_sd <- vapply(seq_along(sd_by_env),
                   function(i) unname(got[[paste0("sd_gen__env", i)]]),
                   numeric(1))
  expect_equal(got_sd, ref, tolerance = 0.20)
})

test_that("the per-level index does not leak into the summary label", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  tab <- .hv_fixture(n_gen = 15L, n_rep = 3L)
  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~ diag(env):gen,
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  comps <- summary(fit)$varcomp$component
  # brms names the same component sd_<inner>__<outer><level>; the two
  # engines must agree, and the generated index must not surface.
  expect_true("sd_gen__env1" %in% comps)
  expect_false(any(grepl("_at_", comps, fixed = TRUE)))
})

test_that("at(f, level) passes the gate with one field for the named level", {
  tab <- .hv_fixture(n_gen = 8L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ at(env, "1"):gen,
    data = tab,
    family = "gaussian"
  )
  res <- lgm_gate(fb)
  expect_true("lgm_compatible" %in% res$capabilities)
  expect_identical(fb$random_terms[[1L]]$outer_levels, "1")
})
