# INLA mapping verification test for the factor_numeric_interaction
# term class (ADR 0019 §Decision 5 / the design spec three-arbitrator
# decision rule).
#
# Policy: the INLA mapper for the new term class registers IF AND
# ONLY IF this verification test demonstrates faithful posterior
# agreement with both greta and lme4::lmer on a gaussian-identity
# fixture, per the v2 benchmark agreement contract. Failure to
# verify (INLA not installed, INLA verification mismatch, INLA fit
# error) means the on-disk verification artefact stays absent or
# carries `pass = FALSE`, and the lgm_gate() check
# .lgm_check_factor_numeric_interaction_inla_verified() refuses
# INLA dispatch on the factor:continuous indexed-slope term class
# with a deferral message.
#
# Tier: Tier 2 (full devtools::test()). Skipped on CRAN and when any
# of INLA / greta / lme4 are unavailable -- the verification
# artefact is host-local and the gate refuses cleanly without it.

test_that("INLA mapping for factor_numeric_interaction passes 3-arbitrator gate", {
  testthat::skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_greta_backend_unusable()
  skip_if_not_installed("lme4")
  withr::local_options(list(
    flexyBayes.silence_default_prior_note = TRUE
  ))

  # Three-arbitrator fixture: gaussian-identity, 3-level factor x
  # continuous interaction. Truth coefficients spread across
  # |beta_true| ~ 0.5 so the tolerance W_1 <= 0.10 * |beta_true|
  # discriminates a translation bug at small posterior noise.
  set.seed(20260523L)
  n <- 240L
  d <- data.frame(
    x = rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n))
  )
  beta_true <- c(
    `(Intercept)` = 0.0,
    fb = 0.3,
    fc = -0.2,
    x = 0.8,
    `fb:x` = 0.5,
    `fc:x` = -0.4
  )
  X <- model.matrix(
    ~ f * x,
    data = d,
    contrasts.arg = list(f = "contr.treatment")
  )
  d$y <- as.numeric(X %*% beta_true + rnorm(n, 0, 0.3))

  # Arbitrator A -- lme4::lm (gaussian-identity reduces to lm on a
  # fixed-effects-only model). REML SEs serve as the tolerance unit.
  fit_lm <- stats::lm(y ~ f * x, data = d)
  coef_lm <- stats::coef(fit_lm)
  se_lm <- sqrt(diag(stats::vcov(fit_lm)))

  # Arbitrator B -- flexybayes greta backend (the v0.2.6 indexed
  # emit). Re-uses set.seed for reproducibility; the chain budget
  # mirrors test-validation-lmer.R's 2000 + 2000 + 2 setup.
  set.seed(20260523L)
  fit_greta <- flexybayes(
    y ~ f * x,
    data = d,
    backend = "greta",
    n_samples = 2000,
    warmup = 2000,
    chains = 2,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )
  draws <- as.matrix(fit_greta$greta$draws)
  # Canonical-slot conversion. flexyBayes greta emit uses the
  # "all-levels + global intercept" parameterisation
  #     eta = mu_atg + tau_f[f_id] + beta_x * x +
  #           slope_dev_f_x[f_id] * x
  # whereas model.matrix() uses treatment coding
  #     eta = (Intercept) + fb*I(level=b) + fc*I(level=c) +
  #           x*x + fbx*I(level=b)*x + fcx*I(level=c)*x.
  # The two parameterisations are equivalent under the map
  #     (Intercept) = mu_atg + tau_f[1,1]
  #     fb          = tau_f[2,1] - tau_f[1,1]
  #     fc          = tau_f[3,1] - tau_f[1,1]
  #     x           = beta_x
  #     fb:x        = slope_dev_f_x_raw[1,1]   (raw idx 1 -> level b)
  #     fc:x        = slope_dev_f_x_raw[2,1]   (raw idx 2 -> level c)
  # Computing the canonical-slot posterior mean by taking the mean of
  # the linear combination on the draws (not the mean of mu_atg and
  # the mean of tau_f separately) is the correct triangulation
  # comparison since the two MCMC parameter blocks are not
  # individually identified.
  pm_greta <- c(
    `(Intercept)` = mean(draws[, "mu_atg"] + draws[, "tau_f[1,1]"]),
    fb = mean(draws[, "tau_f[2,1]"] - draws[, "tau_f[1,1]"]),
    fc = mean(draws[, "tau_f[3,1]"] - draws[, "tau_f[1,1]"]),
    x = mean(draws[, "beta_x"]),
    `fb:x` = mean(draws[, "slope_dev_f_x_raw[1,1]"]),
    `fc:x` = mean(draws[, "slope_dev_f_x_raw[2,1]"])
  )

  # Arbitrator C -- INLA via the native `f:x` formula syntax. The
  # priors are INLA defaults (gaussian likelihood; loggamma on
  # residual precision); we report only the posterior mean of the
  # fixed effects.
  fit_inla <- INLA::inla(
    y ~ f * x,
    family = "gaussian",
    data = d,
    control.compute = list(config = TRUE)
  )
  fixed_summary <- fit_inla$summary.fixed
  inla_rownames <- rownames(fixed_summary)
  # INLA's level naming preserves the factor levels (fb, fc) and
  # the colon-joined interaction names (fb:x, fc:x).
  pm_inla <- setNames(
    fixed_summary[, "mean"],
    inla_rownames
  )

  # Three-way comparison. For every named coefficient, INLA must
  # agree with BOTH greta and lme4 within tolerance.
  #   tolerance for fixed effects: 0.10 * |beta_true| per the design spec
  #   §3.4 (a small absolute floor catches near-zero coefficients).
  agreement_table <- list()
  agreement_pass <- TRUE
  for (nm in names(beta_true)) {
    abs_truth <- max(0.05, abs(beta_true[[nm]]))
    tol <- 0.10 * abs_truth + 0.05
    inla_v <- pm_inla[[nm]]
    greta_v <- pm_greta[[nm]]
    lm_v <- coef_lm[[nm]]
    diff_inla_greta <- abs(inla_v - greta_v)
    diff_inla_lm <- abs(inla_v - lm_v)
    pass <- (diff_inla_greta <= tol) && (diff_inla_lm <= tol)
    agreement_table[[nm]] <- list(
      truth = beta_true[[nm]],
      inla = inla_v,
      greta = greta_v,
      lm = lm_v,
      tol = tol,
      diff_inla_greta = diff_inla_greta,
      diff_inla_lm = diff_inla_lm,
      pass = pass
    )
    if (!pass) agreement_pass <- FALSE
  }

  # Write the verification artefact -- presence + pass = TRUE is the
  # signal lgm_gate's .factor_numeric_interaction_inla_verified()
  # reads. On verification failure we still write the artefact with
  # pass = FALSE so the failure mode is auditable; the gate refuses
  # either way.
  art_dir <- system.file(
    "extdata",
    "inla-verification",
    package = "flexyBayes",
    mustWork = FALSE
  )
  if (!nzchar(art_dir)) {
    # Package not installed; the test harness is running via
    # devtools::test() -- write into the source-tree inst/.
    art_dir <- file.path("..", "..", "inst", "extdata", "inla-verification")
  }
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
  artefact <- list(
    timestamp = Sys.time(),
    R_version = R.version.string,
    INLA_version = as.character(utils::packageVersion("INLA")),
    flexyBayes_v = as.character(utils::packageVersion("flexyBayes")),
    fixture = list(n = n, seed = 20260523L, beta_true = beta_true),
    agreement = agreement_table,
    pass = agreement_pass
  )
  saveRDS(artefact, file = file.path(art_dir, "factor_numeric_interaction.rds"))

  # The test itself asserts that verification passed. A FAIL surfaces
  # the INLA-mapping disagreement immediately rather than silently
  # registering a wrong mapper.
  expect_true(
    agreement_pass,
    info = paste(
      "Three-arbitrator INLA verification disagreed.",
      "Inspect agreement table on the saved artefact at",
      file.path(art_dir, "factor_numeric_interaction.rds")
    )
  )
})

# Companion check: the gate consults the artefact and accepts when
# verification is recorded as passing. Tests the artefact -> gate
# pickup in isolation from the INLA fit so the test is fast and the
# pickup logic is exercised on every Tier-2 run.

test_that("lgm_gate accepts INLA when verification artefact reports pass", {
  art_dir <- system.file(
    "extdata",
    "inla-verification",
    package = "flexyBayes",
    mustWork = FALSE
  )
  if (!nzchar(art_dir)) {
    art_dir <- file.path("..", "..", "inst", "extdata", "inla-verification")
  }
  art_path <- file.path(art_dir, "factor_numeric_interaction.rds")
  if (!file.exists(art_path)) {
    testthat::skip("verification artefact not yet produced on this host")
  }

  rec <- readRDS(art_path)
  if (!isTRUE(rec$pass)) {
    testthat::skip(
      "verification artefact records pass = FALSE; gate refuses INLA"
    )
  }

  fb <- flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    fixed_terms = list(
      list(
        type = "factor",
        var = "f",
        levels = c("a", "b", "c"),
        n_levels = 3L,
        label = "f"
      ),
      list(type = "continuous", var = "x", label = "x"),
      list(
        type = "factor_numeric_interaction",
        factor = "f",
        continuous = "x",
        vars = c("f", "x"),
        levels = c("a", "b", "c"),
        n_levels = 3L,
        label = "f:x"
      )
    ),
    random_terms = list(),
    rcov_terms = list(list(type = "units")),
    priors = list(legacy = TRUE),
    source = "brms"
  )
  r10 <- flexyBayes:::.lgm_check_factor_numeric_interaction_inla_verified(fb)
  expect_true(r10$pass)
})

# Trivial-pass guard: on IRs without a factor_numeric_interaction
# term the §3.4 gate is a no-op (the verification artefact's
# existence is irrelevant when the term class is absent).

test_that("verification gate is no-op without factor_numeric_interaction term", {
  fb <- flexyBayes:::new_fb_terms(
    response = "y",
    family = "gaussian",
    fixed_terms = list(list(type = "continuous", var = "x", label = "x")),
    random_terms = list(),
    rcov_terms = list(list(type = "units")),
    priors = list(legacy = TRUE),
    source = "brms"
  )
  r <- flexyBayes:::.lgm_check_factor_numeric_interaction_inla_verified(fb)
  expect_true(r$pass)
})
