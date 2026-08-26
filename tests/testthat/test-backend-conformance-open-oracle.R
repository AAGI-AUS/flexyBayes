# =============================================================================
# test-backend-conformance-open-oracle.R
#
# The re-pointed backend-conformance battery for the two active engines
# (brms + INLA), the v6-I13 gate that predates the eventual 0.9.3 withdrawal
# of the other two once-registered engines (see NEWS.md). It differs from
# the legacy test-backend-conformance.R in the
# reference tier: each engine is tested as a CANDIDATE against an INDEPENDENT
# OPEN ORACLE (base `lm` for a Gaussian GLM, `lme4` REML for a Gaussian LMM),
# not against a peer engine -- so fidelity is proven against something neither
# engine authored (the fidelity-first / independent-oracle principle). This
# is the check the peer-referenced legacy battery could not make for the
# structured cells, and it is why an engine's covariance defect can no longer
# pass silently.
#
# Candidate parameters are extracted through the canonical accessor
# `.fb_canonical_draws()`, so INLA's raw `inla.posterior.sample()` names
# ("(Intercept):1", "Precision for ...", "Predictor:k") resolve to the same
# canonical tokens ((Intercept), the fixed terms, sigma, sd_<group>) as every
# other backend, with variance components on the SD scale.
#
# Heavy: brms compiles + samples a Stan model; INLA fits a Laplace model.
# Gated `skip_on_cran()` + `skip_on_ci()`, and skipped when an engine (or the
# lme4 oracle) is unavailable. An engine that IS available and whose fit then
# errors is a failure, not a skip -- see `.ooc_fit_ok()`. A backend joins by
# adding one descriptor to `.ooc_backends()`; the open oracles live in the
# per-cell fixtures.
# =============================================================================

# ---- tolerances -------------------------------------------------------------
# Fixed effects + residual SD recover the oracle tightly; a variance-component
# SD legitimately differs more (REML point vs Bayesian posterior mean under a
# weakly-informative prior), so its band is wider by design.
.ooc_tol <- list(fixed = 0.06, sigma = 0.12, re_sd = 0.25)

# ---- fixtures + their open oracles ------------------------------------------
.ooc_glm <- function() {
  set.seed(424242L)
  n <- 150L
  x <- stats::rnorm(n)
  list(
    data = data.frame(y = 1 + 2 * x + stats::rnorm(n, sd = 0.5), x = x),
    formula = y ~ x,
    random = NULL,
    family = "gaussian"
  )
}

.ooc_lmm <- function() {
  set.seed(20260724L)
  n_e <- 4L; n_g <- 25L; n_r <- 8L
  env <- factor(rep(seq_len(n_e), each = n_g * n_r))
  gen <- factor(rep(rep(seq_len(n_g), each = n_r), times = n_e))
  e_e <- stats::rnorm(n_e, 10, 2)
  e_g <- stats::rnorm(n_g, 0, 1.2)
  list(
    data = data.frame(
      yield = e_e[as.integer(env)] + e_g[as.integer(gen)] +
        stats::rnorm(length(env), 0, 0.8),
      environment = env,
      genotype = gen
    ),
    formula = yield ~ environment,
    random = ~ genotype,
    family = "gaussian"
  )
}

.ooc_vm <- function() {
  set.seed(20260724L)
  n_g <- 40L; n_rep <- 6L; n_marker <- 200L
  zm <- matrix(stats::rbinom(n_g * n_marker, 2L, 0.3), n_g, n_marker)
  zc <- scale(zm, scale = FALSE)
  g_mat <- tcrossprod(zc) / n_marker + diag(n_g) * 1e-3
  dimnames(g_mat) <- list(paste0("g", 1:n_g), paste0("g", 1:n_g))
  u <- as.numeric(t(chol(g_mat)) %*% stats::rnorm(n_g))   # sd_g truth = 1
  geno <- factor(rep(paste0("g", 1:n_g), times = n_rep),
                 levels = paste0("g", 1:n_g))
  list(
    data = data.frame(
      y = 5 + u[as.integer(geno)] + stats::rnorm(length(geno), 0, 0.7),
      geno = geno
    ),
    G = g_mat
  )
}

# ---- canonical extraction ---------------------------------------------------
.ooc_canon <- function(fit) flexyBayes:::.fb_canonical_draws(fit)
.ooc_mean <- function(co, nm) if (nm %in% names(co)) mean(co[[nm]]) else NA_real_
.ooc_re_sd <- function(co, key = "gen") {
  nm <- grep(key, names(co), ignore.case = TRUE, value = TRUE)
  nm <- nm[grepl("sd|sigma", nm, ignore.case = TRUE)]
  if (length(nm)) mean(co[[nm[1]]]) else NA_real_
}

# ---- candidate descriptors (a backend joins by adding one entry) ------------
.ooc_backends <- function() {
  list(
    inla = list(
      name = "inla",
      skip_reason = function() {
        if (nzchar(system.file(package = "INLA"))) NULL else "INLA not installed"
      },
      capability = flexyBayes:::.capability_inla
    ),
    brms = list(
      name = "brms",
      skip_reason = function() {
        if (nzchar(system.file(package = "brms"))) NULL else "brms not installed"
      },
      capability = flexyBayes:::.capability_brms
    )
  )
}

# Sampling controls keep the MCMC candidate well enough mixed to clear the
# OO.5 ESS floor (a 25-level variance component mixes slowly at defaults);
# INLA ignores them.
.ooc_sampling <- list(n_samples = 2500L, warmup = 1000L, chains = 4L)

# Fail, do not skip, when an installed engine's fit errors.
#
# Each battery pass already skips on `desc$skip_reason()`, so by the time
# a fit is attempted the engine is present. Skipping on the error as well
# conflated "this engine is not on the host" with "this engine is on the
# host and the fit broke" -- and the second is the regression the battery
# exists to catch (2026-08-16 adversarial review, P2-4). Returns TRUE
# when the caller may proceed.
.ooc_fit_ok <- function(fit, engine, cell_id) {
  if (!inherits(fit, "error")) {
    return(TRUE)
  }
  testthat::fail(paste0(
    "[", engine, "] ", cell_id, ": the engine is installed and the fit ",
    "errored, which is a regression rather than an unavailable engine -- ",
    conditionMessage(fit)
  ))
  FALSE
}

.ooc_fit <- function(engine, cell) {
  args <- c(
    list(
      cell$formula,
      data = cell$data,
      family = cell$family,
      backend = engine,
      verbose = FALSE
    ),
    .ooc_sampling
  )
  if (!is.null(cell$random)) args$random <- cell$random
  suppressMessages(do.call(flexybayes, args))
}

# vm/GBLUP needs the carrier matched per engine: INLA's generic0 path takes the
# PRECISION (Q = G^-1); brms takes the dense covariance (its gr(cov = K)).
.ooc_fit_vm <- function(engine, cell) {
  if (engine == "inla") {
    q_prec <- solve(cell$G)
    suppressMessages(flexybayes(
      y ~ 1, random = ~ vm(geno, cov = fb_cov(Qprec, type = "precision")),
      data = cell$data, known_matrices = list(Qprec = q_prec),
      backend = "inla", verbose = FALSE
    ))
  } else {
    suppressMessages(do.call(flexybayes, c(
      list(
        y ~ 1, random = ~ vm(geno, Gmat), data = cell$data,
        known_matrices = list(Gmat = cell$G), backend = engine, verbose = FALSE
      ),
      .ooc_sampling
    )))
  }
}

# ---- the battery (one parameterised pass per candidate) ---------------------
for (.ooc_key in names(.ooc_backends())) {
  local({
    desc <- .ooc_backends()[[.ooc_key]]
    reason <- desc$skip_reason()

    test_that(sprintf("[%s] OO.1 GLM recovers the open oracle (base lm)", desc$name), {
      skip_on_cran(); skip_on_ci()
      if (!is.null(reason)) skip(reason)
      cell <- .ooc_glm()
      orc <- stats::lm(cell$formula, cell$data)
      of <- stats::coef(orc); os <- summary(orc)$sigma
      fit <- tryCatch(.ooc_fit(desc$name, cell), error = function(e) e)
      if (!.ooc_fit_ok(fit, desc$name, "OO.1")) return(invisible(NULL))
      co <- .ooc_canon(fit)
      # Canonical tokens present (the INLA-canonicalisation contract) ...
      expect_true(all(c("(Intercept)", "x", "sigma") %in% names(co)))
      # ... and the candidate recovers the INDEPENDENT lm oracle.
      expect_equal(.ooc_mean(co, "(Intercept)"), of[["(Intercept)"]],
                   tolerance = .ooc_tol$fixed)
      expect_equal(.ooc_mean(co, "x"), of[["x"]], tolerance = .ooc_tol$fixed)
      expect_equal(.ooc_mean(co, "sigma"), os, tolerance = .ooc_tol$sigma)
    })

    test_that(sprintf("[%s] OO.2 LMM recovers the open oracle (lme4 REML)", desc$name), {
      skip_on_cran(); skip_on_ci()
      if (!is.null(reason)) skip(reason)
      skip_if_not_installed("lme4")
      cell <- .ooc_lmm()
      lm4 <- lme4::lmer(yield ~ environment + (1 | genotype),
                        data = cell$data, REML = TRUE)
      vc <- as.data.frame(lme4::VarCorr(lm4))
      o_int <- lme4::fixef(lm4)[["(Intercept)"]]
      o_sdg <- vc$sdcor[vc$grp == "genotype"]
      o_sdr <- vc$sdcor[vc$grp == "Residual"]
      fit <- tryCatch(.ooc_fit(desc$name, cell), error = function(e) e)
      if (!.ooc_fit_ok(fit, desc$name, "OO.2")) return(invisible(NULL))
      co <- .ooc_canon(fit)
      expect_true(all(c("(Intercept)", "sigma") %in% names(co)))
      expect_equal(.ooc_mean(co, "(Intercept)"), o_int, tolerance = .ooc_tol$fixed)
      expect_equal(.ooc_mean(co, "sigma"), o_sdr, tolerance = .ooc_tol$sigma)
      gsd <- .ooc_re_sd(co, "gen")
      expect_false(is.na(gsd))                      # canonical sd_<group> resolved
      expect_equal(gsd, o_sdg, tolerance = .ooc_tol$re_sd)
    })

    test_that(sprintf("[%s] OO.4 structured-cov vm/GBLUP recovers sommer REML", desc$name), {
      skip_on_cran(); skip_on_ci()
      if (!is.null(reason)) skip(reason)
      skip_if_not_installed("sommer")
      cell <- .ooc_vm()
      som <- sommer::mmer(y ~ 1, random = ~ sommer::vsr(geno, Gu = cell$G),
                          rcov = ~ units, data = cell$data, verbose = FALSE)
      vc <- summary(som)$varcomp
      o_sdg <- sqrt(vc$VarComp[grepl("geno", rownames(vc), ignore.case = TRUE)][1])
      o_sde <- sqrt(vc$VarComp[grepl("units|resid", rownames(vc), ignore.case = TRUE)][1])
      fit <- tryCatch(.ooc_fit_vm(desc$name, cell), error = function(e) e)
      if (!.ooc_fit_ok(fit, desc$name, "OO.4")) return(invisible(NULL))
      co <- .ooc_canon(fit)
      gsd <- .ooc_re_sd(co, "geno")
      # The known-covariance carrier is faithful iff the recovered variance
      # components match the INDEPENDENT sommer GBLUP REML (an iid-instead-of-K
      # or Cholesky-vs-covariance defect would miss sommer's sd_g).
      expect_false(is.na(gsd))                       # canonical sd_<vm-group> resolved
      expect_equal(.ooc_mean(co, "sigma"), o_sde, tolerance = .ooc_tol$sigma)
      expect_equal(gsd, o_sdg, tolerance = .ooc_tol$re_sd)
    })

    test_that(sprintf("[%s] OO.5 MCMC diagnostics meet the floor (Laplace exempt)", desc$name), {
      skip_on_cran(); skip_on_ci()
      if (!is.null(reason)) skip(reason)
      skip_if_not_installed("posterior")
      cell <- .ooc_lmm()                              # the slowest-mixing Gaussian cell
      fit <- tryCatch(.ooc_fit(desc$name, cell), error = function(e) e)
      if (!.ooc_fit_ok(fit, desc$name, "OO.5")) return(invisible(NULL))
      d <- flexyBayes:::.fb_mcmc_diagnostics(fit)
      if (!isTRUE(d$applicable)) {
        # INLA / Laplace is deterministic -- R-hat / bulk-ESS are undefined on
        # samples that are not Markov chains, so the gate is exempt here.
        expect_false(isTRUE(d$applicable))
      } else {
        # MCMC candidate: enforce a convergence floor before its posterior
        # means are trusted by the recovery clauses (OO.1/OO.2/OO.4).
        expect_lte(d$max_rhat, 1.05)
        expect_gte(d$min_ess_bulk, 100)
      }
    })

    test_that(sprintf("[%s] OO.3 capability predicate is a total decision", desc$name), {
      skip_on_cran(); skip_on_ci()
      if (!is.null(reason)) skip(reason)
      if (is.null(desc$capability)) skip("descriptor declares no capability predicate")
      d <- .ooc_glm()$data
      fb_glm <- fb_from_asreml(fixed = y ~ x, data = d, family = "gaussian")
      expect_true(isTRUE(desc$capability(fb_glm)))    # admits a plain GLM
      d$g <- factor(rep(seq_len(5L), length.out = nrow(d)))
      fb_vm <- suppressWarnings(tryCatch(
        fb_from_asreml(fixed = y ~ x, random = ~ vm(g), data = d, family = "gaussian"),
        error = function(e) NULL
      ))
      skip_if(is.null(fb_vm), "structured-cov IR not constructible in this build")
      cap <- desc$capability(fb_vm)
      # Total decision: either admits (TRUE) or refuses with a structured
      # reason string -- never NULL / NA / an error.
      expect_true(isTRUE(cap) || (is.character(cap) && nzchar(cap)))
    })
  })
}
