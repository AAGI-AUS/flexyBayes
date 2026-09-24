# tools/capability_grid_matrix.R -- the capability-matrix CSV builder.
#
# Both tools/execution_grid.R (the live grid, which drives its M-section
# test cells off the same list) and tools/generate_capability_matrix.R (the
# installed-evidence copy, which needs only the CSV) write
# inst/validation/execution_grid/capability_matrix.csv. One function here
# builds the `cell, model_class, spelling, backend, verdict` rows for both,
# so the two writers cannot drift back into the state
# tests/testthat/test-installed-evidence-parity.R exists to catch: a shipped
# CSV four model classes and one release behind the registry it claims to
# describe (last written by 4c81984 at 0.9.3, found 2026-09-16).
#
# `.fb_capability_grid_m_calls()` maps each model class the M section can
# execute to its minimal live call, one list entry per class (a second
# entry for a second spelling of the same class, e.g. `vm()` / `ped()`,
# shares the class's row and does not add one). Row order in the built CSV
# follows this list's order, not the registry's own row order -- two rows
# the grid added after the original table (`Per-trial separable AR1 field`,
# `Grouped temporal AR1 field`) sit later in the list than in the table.
# `.fb_capability_grid_m_executed_elsewhere()` names the registry rows no M
# cell can reach (the joint shrinkage-prior rows, each routed through
# `fb_select()` rather than a bare `flexybayes()` call); every registry row
# must appear in one list or the other, checked by
# `.fb_capability_grid_check_drift()` before either script builds a cell or
# a row from them.

.fb_capability_grid_slug <- function(x) {
  s <- tolower(gsub("[^a-z0-9]+", "-", tolower(x)))
  gsub("(^-|-$)", "", s)
}

# --- the M section's model-class-to-call list --------------------------- #

.fb_capability_grid_m_calls <- function() {
  list(
    list(
      model_class = "Gaussian LMM, simple random intercept",
      inla = 'flexybayes(y ~ env, random = ~ gen, data = d_cap)',
      brms = 'flexybayes(y ~ env, random = ~ gen, data = d_cap)'
    ),
    list(
      model_class = paste(
        "GLMM (binomial, Poisson, negative binomial, gamma, beta),",
        "simple random effect"
      ),
      # One representative family here; S3b runs every admitted spelling
      # on both engines.
      inla = 'flexybayes(bin ~ env, random = ~ gen, data = d_cap, family = "binomial")',
      brms = 'flexybayes(bin ~ env, random = ~ gen, data = d_cap, family = "binomial")'
    ),
    list(
      model_class = "Hurdle gamma (zero mass plus a positive gamma part)",
      inla = 'flexybayes(hg ~ env, random = ~ gen, data = d_cap, family = "hurdle_gamma")',
      brms = 'flexybayes(hg ~ env, random = ~ gen, data = d_cap, family = "hurdle_gamma")'
    ),
    list(
      model_class = "Uncorrelated random slope",
      inla = 'flexybayes(y ~ x + (x || gen), data = d_cap)',
      brms = 'flexybayes(y ~ x + (x || gen), data = d_cap)'
    ),
    list(
      model_class = "Factor-by-numeric fixed interaction",
      inla = 'flexybayes(y ~ env * x, random = ~ gen, data = d_cap)',
      brms = 'flexybayes(y ~ env * x, random = ~ gen, data = d_cap)'
    ),
    list(
      model_class = "Correlated random slope",
      inla = 'flexybayes(y ~ x + (x | gen), data = d_cap)',
      brms = 'flexybayes(y ~ x + (x | gen), data = d_cap)'
    ),
    list(
      model_class = "Nested / interaction random effects, multi-stratum",
      inla = 'flexybayes(y ~ env, random = ~ gen:env, data = d_cap)',
      brms = 'flexybayes(y ~ env, random = ~ gen:env, data = d_cap)'
    ),
    list(
      model_class = "Heterogeneous variance by factor level",
      inla = 'flexybayes(y ~ env, random = ~ diag(env):gen, data = d_cap)',
      brms = 'flexybayes(y ~ env, random = ~ diag(env):gen, data = d_cap)'
    ),
    list(
      model_class = "Unstructured genotype-by-environment covariance",
      inla = 'flexybayes(y ~ env, random = ~ us(env):gen, data = d_cap)',
      brms = 'flexybayes(y ~ env, random = ~ us(env):gen, data = d_cap)'
    ),
    list(
      model_class = "Heterogeneous variances with one shared correlation",
      inla = 'flexybayes(y ~ env, random = ~ corh(env):gen, data = d_cap)',
      brms = 'flexybayes(y ~ env, random = ~ corh(env):gen, data = d_cap)'
    ),
    list(
      model_class = "Heterogeneous residual by factor level",
      inla = paste0('flexybayes(y ~ env, random = ~ gen, ',
                    'residual = ~ dsum(~ units | env), data = d_cap)'),
      brms = paste0('flexybayes(y ~ env, random = ~ gen, ',
                    'residual = ~ dsum(~ units | env), data = d_cap)')
    ),
    list(
      model_class = paste(
        "Baseline heteroscedastic MET: iid genotype, iid",
        "genotype-by-environment and a sectioned residual"
      ),
      inla = paste0('flexybayes(y ~ env, random = ~ gen + gen:env, ',
                    'residual = ~ dsum(~ units | env), data = d_cap)'),
      brms = paste0('flexybayes(y ~ env, random = ~ gen + gen:env, ',
                    'residual = ~ dsum(~ units | env), data = d_cap)')
    ),
    list(
      model_class = "Factor-analytic genotype-by-environment covariance",
      inla = 'flexybayes(y ~ env, random = ~ fa(env, 1):gen, data = d_cap)',
      brms = 'flexybayes(y ~ env, random = ~ fa(env, 1):gen, data = d_cap)'
    ),
    list(
      model_class = "Multi-trait covariance",
      inla = paste0('flexybayes(y ~ env, random = ~ us(trait):vm(gen), ',
                    'data = d_cap, known_matrices = list(K = k_cap))'),
      brms = paste0('flexybayes(y ~ env, random = ~ us(trait):vm(gen), ',
                    'data = d_cap, known_matrices = list(K = k_cap))')
    ),
    list(
      # The row's own note splits the carriers: INLA takes the sparse
      # precision, brms the dense covariance.
      model_class = "Known-covariance genomic / pedigree random effect",
      inla = paste0('flexybayes(y ~ env, ',
                    'random = ~ vm(gen, cov = fb_cov(Q, type = "precision")), ',
                    'data = d_cap, known_matrices = list(Q = q_cap))'),
      brms = paste0('flexybayes(y ~ env, random = ~ vm(gen, K), ',
                    'data = d_cap, known_matrices = list(K = k_cap))')
    ),
    list(
      # `ped()` and `vm()` are two spellings of one known-covariance term
      # -- every branch in the package tests `%in% c("vm", "ped")` -- so the
      # capability is one row and the spelling is a surface on it. The grid
      # executed `vm()` only, and `ped()` was advertised without ever being
      # fitted (audit item 9, 2026-08-27). This cell fits the same row
      # through the other spelling, on the carrier each engine documents:
      # INLA the sparse precision, brms the dense covariance.
      model_class = "Known-covariance genomic / pedigree random effect",
      spelling_variant = "ped",
      inla = paste0('flexybayes(y ~ env, ',
                    'random = ~ ped(gen, cov = fb_cov(Q, type = "precision")), ',
                    'data = d_cap, known_matrices = list(Q = q_cap))'),
      brms = paste0('flexybayes(y ~ env, random = ~ ped(gen, K), ',
                    'data = d_cap, known_matrices = list(K = k_cap))')
    ),
    list(
      model_class = "Separable AR1 spatial field",
      inla = 'flexybayes(y ~ 1, random = ~ ar1(row):ar1(col), data = d_capgrid)',
      brms = 'flexybayes(y ~ 1, random = ~ ar1(row):ar1(col), data = d_capgrid)'
    ),
    list(
      model_class = "Univariate P-spline",
      inla = 'flexybayes(y ~ 1, random = ~ spl(x), data = d_cap)',
      brms = 'flexybayes(y ~ 1, random = ~ spl(x), data = d_cap)'
    ),
    list(
      model_class = "Observation weights (Gaussian, identity link)",
      inla = 'flexybayes(y ~ env, random = ~ gen, data = d_cap, weights = w_cap)',
      brms = 'flexybayes(y ~ env, random = ~ gen, data = d_cap, weights = w_cap)'
    ),
    list(
      model_class = "Exact sufficient-statistic aggregation (binomial, Poisson)",
      inla = paste0('flexybayes(cnt ~ env, random = ~ gen, data = d_agg, ',
                    'aggregate = TRUE, family = "poisson")'),
      brms = paste0('flexybayes(cnt ~ env, random = ~ gen, data = d_agg, ',
                    'aggregate = TRUE, family = "poisson")')
    ),
    list(
      model_class = "Exact sufficient-statistic aggregation (Gaussian)",
      inla = paste0('flexybayes(y ~ env, random = ~ gen, data = d_agg, ',
                    'aggregate = TRUE)'),
      brms = paste0('flexybayes(y ~ env, random = ~ gen, data = d_agg, ',
                    'aggregate = TRUE)')
    ),
    list(
      model_class = "Per-trial separable AR1 field",
      inla = paste0('flexybayes(y ~ 1, random = ~ id(trial):ar1(row):ar1(col), ',
                    'data = d_ar1_at_trial)'),
      brms = paste0('flexybayes(y ~ 1, random = ~ id(trial):ar1(row):ar1(col), ',
                    'data = d_ar1_at_trial)')
    ),
    list(
      model_class = "Grouped temporal AR1 field",
      inla = 'flexybayes(y ~ 1, random = ~ ar1(t):g, data = d_ar1_tg)',
      brms = 'flexybayes(y ~ 1, random = ~ ar1(t):g, data = d_ar1_tg)'
    )
  )
}

# Matrix rows no M cell can execute, each with where it is executed instead.
# A joint shrinkage prior reaches no flexybayes() call: fb_select() is its
# only public route, and its reference fit and folds are several brms fits,
# past the per-cell timeout. Listed by name so the gap prints on every run
# (added 2026-09-15; the rows arrived with 0.11.0 without cells).
.fb_capability_grid_m_executed_elsewhere <- function() {
  c(
    "Ridge shrinkage prior on candidate fixed effects" = paste(
      "fitted live by fb_select() on brms in",
      "review/release_0.12.0/smoke/smoke_installed_package.R;",
      "Stan code asserted in tests/testthat/test-penalty-stancode.R"
    ),
    "Laplace shrinkage prior on candidate fixed effects" =
      "Stan code asserted in tests/testthat/test-penalty-stancode.R; not sampled",
    "Elastic-net shrinkage prior on candidate fixed effects" =
      "Stan code asserted in tests/testthat/test-penalty-stancode.R; not sampled",
    "Regularised horseshoe shrinkage prior on candidate fixed effects" =
      "Stan code asserted in tests/testthat/test-penalty-stancode.R; not sampled",
    "R2D2 shrinkage prior on candidate fixed effects" =
      "Stan code asserted in tests/testthat/test-penalty-stancode.R; not sampled"
  )
}

# --- drift check and CSV rows -------------------------------------------- #

.fb_capability_grid_check_drift <- function(cap, m_calls, elsewhere) {
  named <- vapply(m_calls, function(z) z$model_class, character(1L))
  missing <- setdiff(cap$model_class, c(named, names(elsewhere)))
  extra <- setdiff(named, cap$model_class)
  if (length(missing) > 0L || length(extra) > 0L) {
    stop(
      "The capability matrix and the grid's M section have drifted apart.\n",
      if (length(missing) > 0L) {
        paste0("  matrix rows with no grid cell: ",
               paste(missing, collapse = " | "), "\n")
      },
      if (length(extra) > 0L) {
        paste0("  grid cells naming no matrix row: ",
               paste(extra, collapse = " | "), "\n")
      },
      "Every capability cell is executed, or the gate is not F28.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.fb_capability_grid_cells <- function(cap, m_calls) {
  cells <- data.frame(
    cell = character(0), model_class = character(0),
    spelling = character(0), backend = character(0),
    verdict = character(0), stringsAsFactors = FALSE
  )
  for (spec in m_calls) {
    row <- cap[cap$model_class == spec$model_class, , drop = FALSE]
    for (be in c("inla", "brms")) {
      key <- paste0("cap-", .fb_capability_grid_slug(spec$model_class), "-",
                    be)
      if (!key %in% cells$cell) cells <- rbind(cells, data.frame(
        cell = key,
        model_class = spec$model_class,
        spelling = row$spelling,
        backend = be,
        verdict = row[[be]],
        stringsAsFactors = FALSE
      ))
    }
  }
  cells
}
