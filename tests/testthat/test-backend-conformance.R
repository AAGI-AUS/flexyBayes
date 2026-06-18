# =============================================================================
# test-backend-conformance.R -- the backend-contract conformance battery.
#
# The backend contract is the interface any inference backend must honour to be
# a first-class flexyBayes engine. This
# file is the *executable* form of that contract's "Conformance + drift
# detection" section: it asserts a backend honours the required clauses --
#
#   1. Construction  -- an IR-targetable, re-entrant model builder.
#   2. Output        -- a `posterior::draws_array`.
#   3. Canonical names + posterior recovery -- canonical parameter tokens, and a
#      known posterior recovered to tolerance.
#   4. Agreement     -- the fit agrees with a reference engine within the
#      SBC-calibrated `triangulate()` threshold (the live drift signal).
#   5. Capability declaration -- the predicate admits what the backend can fit
#      and refuses its catch-up boundary with a structured reason.
#
# The battery is BACKEND-AGNOSTIC by construction. Each backend under test is one
# descriptor in `.conformance_backends()`; bringing koine (the synthesised fourth
# opinion) onto the engine axis is adding one descriptor here, not a new test
# file. That reusability is the backend-contract design (depend on the
# contract, never a producer's internals).
#
# Heavy: each backend fits two out-of-process torch-NUTS models plus an INLA
# reference. Gated `skip_on_cran()` + `skip_on_ci()`, and skipped when the
# backend's runnable version (or its dev source) is unavailable. Locally, point
# at a gretaR dev tree with `options(flexyBayes.gretaR_home = <src>)` or
# `Sys.setenv(GRETAR_HOME = <src>)`; no path is hard-coded here (r_style I13).
# =============================================================================

# ---- SBC-calibrated conformance threshold (clause 4) ------------------------
# The agreement statistic is T = max_j W1_j / sigma_j over the common parameters,
# where W1_j is the 1D Wasserstein-1 distance between the two engines' marginals
# and sigma_j is the reference engine's posterior SD for parameter j. The
# triangulation-calibration study put the production Gaussian trio's 95%
# agreement-null threshold at 0.327; the battery asserts
# against a margined ceiling so an honest fit does not flake on MCMC noise while
# a genuinely divergent (wrong-posterior) backend still trips it.
.CONFORMANCE_T_NULL_95 <- 0.327
.CONFORMANCE_T_CEILING <- 0.50

# ---- fixture: a strong-signal Gaussian GLM with a known posterior -----------
# Truth: (Intercept) = 1, x = 2, sigma = 0.5. Two predictors so the re-entrancy
# clause can fit a *second*, structurally different model in the same process.
.conformance_data <- function() {
  set.seed(424242L)
  n <- 150L
  x <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  data.frame(
    y = 1.0 + 2.0 * x + stats::rnorm(n, sd = 0.5),
    x = x,
    x2 = x2
  )
}

.conformance_truth <- c("(Intercept)" = 1.0, "x" = 2.0, "sigma" = 0.5)

# ---- reference engine (clause 4) --------------------------------------------
# INLA is the reference: deterministic (Laplace), fast, and exact for a Gaussian
# GLM -- the cleanest fixed point to triangulate every candidate backend against.
.conformance_reference_available <- function() {
  nzchar(system.file(package = "INLA"))
}

.conformance_fit_reference <- function(formula, data) {
  suppressMessages(flexybayes(
    formula,
    data = data,
    backend = "inla",
    verbose = FALSE,
    mcmc_verbose = FALSE
  ))
}

# ---- the parameterised backend registry -------------------------------------
# A descriptor carries:
#   name        -- the canonical `backend=` string + display label.
#   skip_reason -- function() returning NULL when runnable, else a reason string.
#   prepare     -- function() with any side-effecting setup (e.g. point the
#                  worker at a dev source); returns a zero-arg cleanup thunk.
#   fit         -- function(formula, data) returning a flexyBayes fit.
#   capability  -- function(fb_terms) the backend's capability predicate
#                  (TRUE | <reason>), or NULL to skip the clause-5 check.
#
# A backend joins the backend conformance surface by adding ONE entry.

.conformance_gretaR <- function() {
  floor <- flexyBayes:::.GRETAR_VERSION_FLOOR
  home <- getOption("flexyBayes.gretaR_home", Sys.getenv("GRETAR_HOME", ""))
  installed_ok <- nzchar(system.file(package = "gretaR")) &&
    utils::packageVersion("gretaR") >= floor

  list(
    name = "gretaR",
    skip_reason = function() {
      if (!nzchar(home) && !installed_ok) {
        return(paste0(
          "gretaR >= ",
          floor,
          " not available -- set ",
          "options(flexyBayes.gretaR_home = <src>) or install a newer gretaR"
        ))
      }
      NULL
    },
    prepare = function() {
      if (nzchar(home)) {
        old <- options(flexyBayes.gretaR_home = home)
        return(function() options(old))
      }
      function() invisible(NULL)
    },
    fit = function(formula, data) {
      suppressMessages(flexybayes(
        formula,
        data = data,
        backend = "gretaR",
        n_samples = 1500,
        warmup = 750,
        chains = 4,
        verbose = FALSE,
        mcmc_verbose = FALSE
      ))
    },
    capability = flexyBayes:::.capability_gretaR
  )
}

# The registry. A backend joins the backend conformance surface by adding
# one entry here -- gretaR is active. (koine moved to flexyBayesOrchestra in the
# lean-core split, 2026-06-06; it is no longer a core dormant backend.)
.conformance_backends <- function() {
  list(
    gretaR = .conformance_gretaR()
  )
}

# ---- shared helpers ----------------------------------------------------------

.conformance_on_cran_or_ci <- function() {
  !identical(Sys.getenv("NOT_CRAN"), "true") ||
    isTRUE(as.logical(Sys.getenv("CI", "false")))
}

# Clause-4 statistic from a triangulate_result: max W1 normalised by the
# reference engine's (fit_b's) posterior SD, over finitely-scaled parameters.
.conformance_agreement_T <- function(tri) {
  m <- tri$metrics
  keep <- is.finite(m$wasserstein_1) & is.finite(m$sd_b) & m$sd_b > 0
  if (!any(keep)) {
    return(NA_real_)
  }
  max(m$wasserstein_1[keep] / m$sd_b[keep])
}

# Fit everything a backend needs ONCE, reused across the clause tests: the
# primary model (y ~ x), a second model (y ~ x2) for the re-entrancy clause, and
# the INLA reference (y ~ x) for the agreement clause.
.conformance_fit_all <- function(desc) {
  d <- .conformance_data()
  list(
    data = d,
    primary = desc$fit(y ~ x, d),
    second = desc$fit(y ~ x2, d),
    reference = .conformance_fit_reference(y ~ x, d)
  )
}

# ---- the battery (one parameterised pass per backend descriptor) ------------

for (.backend_key in names(.conformance_backends())) {
  local({
    desc <- .conformance_backends()[[.backend_key]]
    reason <- desc$skip_reason()
    if (is.null(reason) && !.conformance_reference_available()) {
      reason <- "INLA reference engine unavailable"
    }

    # Fit once (eagerly) when runnable; the clause tests read these. The
    # skip_on_* guards inside each test still protect CRAN / CI even though the
    # eager fit is itself guarded here.
    fits <- NULL
    if (is.null(reason) && !.conformance_on_cran_or_ci()) {
      cleanup <- desc$prepare()
      fits <- tryCatch(.conformance_fit_all(desc), error = function(e) e)
      cleanup()
    }

    fit_class <- paste0("flexybayes_", desc$name)

    test_that(
      sprintf("[%s] C1.1 construction: re-entrant IR builder", desc$name),
      {
        skip_on_cran()
        skip_on_ci()
        if (!is.null(reason)) {
          skip(reason)
        }
        expect_false(
          inherits(fits, "error"),
          info = if (inherits(fits, "error")) conditionMessage(fits) else ""
        )
        # A model built from the IR, plus a structurally different SECOND model in
        # the same process -- proves the builder is re-entrant, not a session-
        # global single-model graph (contract clause 1).
        expect_s3_class(fits$primary, fit_class)
        expect_s3_class(fits$second, fit_class)
      }
    )

    test_that(sprintf("[%s] C1.2 output: posterior::draws_array", desc$name), {
      skip_on_cran()
      skip_on_ci()
      if (!is.null(reason)) {
        skip(reason)
      }
      skip_if(inherits(fits, "error"), "backend fit failed (see C1.1)")
      expect_s3_class(fits$primary$draws, "draws_array")
      dl <- fb_as_draws_simple(fits$primary)
      expect_type(dl, "list")
      expect_true(length(dl) >= 1L)
      expect_false(is.null(names(dl)))
      expect_true(all(vapply(dl, is.numeric, logical(1L))))
    })

    test_that(
      sprintf("[%s] C1.3 canonical names + posterior recovery", desc$name),
      {
        skip_on_cran()
        skip_on_ci()
        if (!is.null(reason)) {
          skip(reason)
        }
        skip_if(inherits(fits, "error"), "backend fit failed (see C1.1)")
        nm <- names(fb_as_draws_simple(fits$primary))
        expect_true(all(c("(Intercept)", "x", "sigma") %in% nm))
        co <- stats::coef(fits$primary)
        # Strong-signal GLM, n = 150: posterior means recover truth comfortably.
        expect_equal(
          co[["(Intercept)"]],
          .conformance_truth[["(Intercept)"]],
          tolerance = 0.30
        )
        expect_equal(co[["x"]], .conformance_truth[["x"]], tolerance = 0.30)
        expect_equal(
          co[["sigma"]],
          .conformance_truth[["sigma"]],
          tolerance = 0.25
        )
      }
    )

    test_that(sprintf("[%s] C1.4 agreement within SBC threshold", desc$name), {
      skip_on_cran()
      skip_on_ci()
      if (!is.null(reason)) {
        skip(reason)
      }
      skip_if(inherits(fits, "error"), "backend fit failed (see C1.1)")
      tri <- triangulate(fits$primary, fits$reference)
      # The fixed effects + sigma must align across the candidate and the
      # reference (INLA); abstention on a parameter is fine, silence is not.
      expect_gte(tri$n_common, 2L)
      stat <- .conformance_agreement_T(tri)
      expect_true(is.finite(stat))
      # Within the SBC-calibrated null (0.327) plus margin -- a wrong posterior
      # would land far outside this band (the calibration study's whole point).
      expect_lt(stat, .CONFORMANCE_T_CEILING)
    })

    test_that(
      sprintf(
        "[%s] C1.5 capability declaration + catch-up boundary",
        desc$name
      ),
      {
        skip_on_cran()
        skip_on_ci()
        if (!is.null(reason)) {
          skip(reason)
        }
        if (is.null(desc$capability)) {
          skip("descriptor declares no capability predicate")
        }
        d <- .conformance_data()
        # Admits a plain GLM ...
        fb_glm <- fb_from_asreml(fixed = y ~ x, data = d, family = "gaussian")
        expect_true(isTRUE(desc$capability(fb_glm)))
        # ... and refuses its declared catch-up boundary (structured covariance)
        # with a structured reason, never a silent wrong answer.
        d$g <- factor(rep(seq_len(5L), length.out = nrow(d)))
        fb_vm <- suppressWarnings(tryCatch(
          fb_from_asreml(
            fixed = y ~ x,
            random = ~ vm(g),
            data = d,
            family = "gaussian"
          ),
          error = function(e) NULL
        ))
        skip_if(
          is.null(fb_vm),
          "structured-cov IR not constructible in this build"
        )
        cap <- desc$capability(fb_vm)
        expect_type(cap, "character")
        expect_match(cap, "structured_cov")
      }
    )
  })
}
