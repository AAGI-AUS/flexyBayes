# tools/execution_grid.R -- the F28 execution grid (recipe 62).
#
# A capability claim's oracle is an EXECUTION, not a registry. The
# generated capability matrix proves that every public surface says what
# the R table says; only a live call proves that the table says what the
# code does. This script runs the claim-derived cross product against an
# INSTALLED build and banks the per-cell ledger that gate F28 reads.
#
# Sections, each with the claim surface its expectations come from:
#
#   M    every cell of the generated capability matrix
#        (`.fb_capability_matrix()`, R/capability_matrix.R), model class
#        x backend, one live minimal call per cell.
#   S1   prior route x family x backend, plus the arrival round-trip
#        through each engine's own prior introspection.
#   S2   the malformed-prior corpus at construction (`fb_prior()`).
#   S2F  one fit per corpus entry, separating a construction that is
#        accepted from a fit that must refuse.
#   S3a  the entry allowlist against the installed engines' own rosters.
#   S3b  every admitted family spelling x backend.
#   S3c  the four response families field finding C2 named.
#   S4   structure x family x backend breadth beyond the matrix cells.
#   S5   grammar-surface parity: bar forms and their ASReml siblings.
#   T    the prior-translation table (target x distribution x engine),
#        each cell asserted to translate AND arrive, or to refuse typed.
#   S6   diagnostic structuring on weakly identified fits.
#
# Harness constraints (recipe 62 "Harness pattern"), all load-bearing:
#
#   * runs against an INSTALLED package, never `pkgload::load_all()` --
#     load_all exports internals and resolves symbols the shipped
#     package does not have, so the grid would grade a tree no user
#     will ever hold;
#   * per-cell `parallel::mclapply` over `callr::r()` children with a
#     timeout, because a grid cell can take the engine subprocess down
#     and a crash is a recorded outcome class, not a dead harness;
#   * sampler cells at reachability budgets (1 chain, 300 draws) --
#     this grid asks whether a cell runs, never whether its posterior
#     is right, which is the test suite's job;
#   * one literal seed, one fixed generation order, so every child sees
#     byte-identical fixtures.
#
# Run it (from the package root, after building and installing a
# tarball off the working tree):
#
#   R CMD build --no-build-vignettes .
#   R CMD INSTALL -l <scratch>/lib flexyBayes_<version>.tar.gz
#   FB_GRID_LIB=<scratch>/lib Rscript tools/execution_grid.R
#
# Environment:
#
#   FB_GRID_LIB     required; the library holding the tarball install.
#   FB_GRID_OUT     output directory; defaults to
#                   inst/validation/execution_grid beside this script.
#   FB_GRID_RAW     directory for the wide per-cell record (code,
#                   messages, prior tables). Defaults to FB_GRID_OUT;
#                   point it outside the package to keep the shipped
#                   ledger lean.
#   FB_GRID_DRYRUN  non-empty: expand the grid, parse every snippet,
#                   print the cell list and stop before any engine runs.
#   FB_GRID_ONLY    regular expression on the cell id; runs that subset
#                   and writes under a `_smoke` suffix so a rehearsal
#                   never overwrites a banked run.
#   FB_GRID_CORES   core count, capped at 8.
#
# `tools/stress_corner_to_corner.R` is NOT this gate; see its header.

set.seed(20260819L)

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- configuration and the installed-build guard ---------------------- #

.script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                    value = TRUE)
.script_dir <- if (length(.script_arg) == 1L) {
  dirname(normalizePath(sub("^--file=", "", .script_arg)))
} else {
  normalizePath(getwd())
}
.pkg_root <- normalizePath(file.path(.script_dir, ".."))
if (!file.exists(file.path(.pkg_root, "DESCRIPTION"))) {
  stop("tools/execution_grid.R must sit in <package root>/tools/.",
       call. = FALSE)
}

.lib <- Sys.getenv("FB_GRID_LIB")
if (!nzchar(.lib) || !dir.exists(.lib)) {
  stop(
    "FB_GRID_LIB must point at the library holding the tarball install. ",
    "The grid never runs against a source tree (recipe 62).",
    call. = FALSE
  )
}
.lib <- normalizePath(.lib)
.libPaths(c(.lib, .libPaths()))

.out_dir <- Sys.getenv("FB_GRID_OUT")
if (!nzchar(.out_dir)) {
  .out_dir <- file.path(.pkg_root, "inst", "validation", "execution_grid")
}
dir.create(.out_dir, recursive = TRUE, showWarnings = FALSE)
.out_dir <- normalizePath(.out_dir)

.raw_dir <- Sys.getenv("FB_GRID_RAW")
if (!nzchar(.raw_dir)) {
  .raw_dir <- .out_dir
}
dir.create(.raw_dir, recursive = TRUE, showWarnings = FALSE)
.raw_dir <- normalizePath(.raw_dir)

.cell_timeout <- 300
.n_cores <- min(8L, as.integer(Sys.getenv("FB_GRID_CORES", "8")))

.guard_installed_build <- function() {
  path <- normalizePath(find.package("flexyBayes"))
  if (!startsWith(path, .lib)) {
    stop(
      "Installed-build guard failed: flexyBayes resolved to ", path,
      ", not the library named by FB_GRID_LIB (", .lib, ").",
      call. = FALSE
    )
  }
  if (file.exists(file.path(path, "DESCRIPTION.in"))) {
    stop(
      "Installed-build guard failed: ", path, " looks like a load_all() ",
      "shim rather than an installed build.",
      call. = FALSE
    )
  }
  ver <- as.character(utils::packageVersion("flexyBayes"))
  stamp <- trimws(sub(";.*$", "",
                      utils::packageDescription("flexyBayes")$Packaged))
  cat("guard: flexyBayes", ver, "|", stamp, "|", path, "\n")
  invisible(list(version = ver, stamp = stamp, path = path))
}

.build <- .guard_installed_build()

.head_sha <- tryCatch(
  trimws(system2("git", c("-C", shQuote(.pkg_root), "rev-parse", "--short",
                          "HEAD"),
                 stdout = TRUE, stderr = FALSE)),
  error = function(e) NA_character_
)
if (length(.head_sha) != 1L || !nzchar(.head_sha)) {
  .head_sha <- NA_character_
}
.run_date <- format(Sys.Date(), "%Y-%m-%d")

# --- shared child preamble -------------------------------------------- #
#
# Every child evaluates this before its own cell code. Seeds are literal
# and the generation order is fixed, so the fixtures are byte-identical
# in every child and across runs.

.preamble_body <- '
suppressMessages(library(flexyBayes))
options(flexyBayes.silence_default_prior_note = TRUE)

# Family fixtures: n = 60, one 10-level grouping factor, one covariate.
set.seed(20260818L)
n <- 60L
g <- factor(rep(seq_len(10L), each = 6L))
x <- rnorm(n)
u <- rnorm(10L, 0, 0.7)
eta <- 0.5 + 0.4 * x + u[as.integer(g)]
mu_logit <- plogis(eta)
d_gaussian <- data.frame(y = eta + rnorm(n), x = x, g = g)
d_poisson <- data.frame(y = rpois(n, exp(eta)), x = x, g = g)
d_negative_binomial <- data.frame(
  y = rnbinom(n, mu = exp(eta), size = 2), x = x, g = g
)
d_binomial <- data.frame(y = rbinom(n, 1L, mu_logit), x = x, g = g)
d_binary <- d_binomial
d_negbinom <- d_negative_binomial
d_gamma <- data.frame(
  y = rgamma(n, shape = 2, rate = 2 / exp(eta)), x = x, g = g
)
d_beta <- data.frame(
  y = rbeta(n, mu_logit * 5, (1 - mu_logit) * 5), x = x, g = g
)

# A point mass at zero plus a gamma positive part, for the C2 families.
set.seed(20260821L)
.hz <- rbinom(n, 1L, 0.7)
d_hurdle_gamma <- data.frame(
  y = ifelse(.hz == 1L, rgamma(n, shape = 2, rate = 2 / exp(eta)), 0),
  x = x, g = g
)

# Structure grid: complete 12 x 5 layout, 10 genotypes x 3 environments
# with two replicates per cell, a 4-level block and a 2-level treatment.
# Complete, so ar1(row):ar1(col) has no holes.
set.seed(20260819L)
s_row <- factor(rep(seq_len(12L), each = 5L))
s_col <- factor(rep(seq_len(5L), times = 12L))
s_g <- factor(paste0("g", rep(seq_len(10L), times = 6L)))
s_f <- factor(rep(c("e1", "e2", "e3"), each = 20L))
s_h <- factor(rep(c("b1", "b2", "b3", "b4"), times = 15L))
s_trt <- factor(rep(c("a", "b"), times = 30L))
s_x <- rnorm(60L)
s_eta <- 0.6 + 0.5 * s_x + rnorm(10L, 0, 0.6)[as.integer(s_g)] +
  rnorm(3L, 0, 0.4)[as.integer(s_f)]
d_struct <- data.frame(
  y = s_eta + rnorm(60L, 0, 0.8),
  ycount = rpois(60L, exp(s_eta)),
  x = s_x, g = s_g, f = s_f, h = s_h, trt = s_trt,
  row = s_row, col = s_col
)
d_struct$z <- rnorm(60L)

# Capability-matrix fixtures. The shapes are the ones the matrix rows
# behaviour anchors use (tests/testthat/test-capability-matrix.R), so a
# grid cell and its anchor exercise the same model on the same design.
set.seed(11L)
d_cap <- expand.grid(
  rep = factor(seq_len(2L)),
  gen = factor(seq_len(6L)),
  env = factor(seq_len(3L))
)
d_cap$block <- factor(rep(seq_len(4L), length.out = nrow(d_cap)))
d_cap$trait <- factor(rep(seq_len(2L), length.out = nrow(d_cap)))
d_cap$x <- rnorm(nrow(d_cap))
d_cap$y <- rnorm(nrow(d_cap), 10, 2)
d_cap$bin <- rbinom(nrow(d_cap), 1L, 0.5)
set.seed(21L)
d_cap$hg <- ifelse(
  rbinom(nrow(d_cap), 1L, 0.7) == 1L,
  rgamma(nrow(d_cap), shape = 2, rate = 0.5),
  0
)
w_cap <- as.numeric(seq_len(nrow(d_cap)))

set.seed(12L)
d_capgrid <- expand.grid(
  row = factor(seq_len(6L)),
  col = factor(seq_len(5L))
)
d_capgrid$y <- rnorm(nrow(d_capgrid))

k_cap <- diag(6)
dimnames(k_cap) <- list(as.character(seq_len(6L)),
                        as.character(seq_len(6L)))
q_cap <- Matrix::Matrix(k_cap, sparse = TRUE)
dimnames(q_cap) <- dimnames(k_cap)

# Aggregation fixture: enough replication per cell for the exact
# sufficient-statistic route to be eligible.
set.seed(13L)
d_agg <- expand.grid(
  rep = factor(seq_len(20L)),
  gen = factor(seq_len(6L)),
  env = factor(seq_len(3L))
)
d_agg$y <- rnorm(nrow(d_agg), 10, 2)

# S6 fixtures. d_weak: many clusters of two, so the cluster variance is
# barely identified. d_funnel: three clusters and a near-zero between-
# cluster signal, the centred-parameterisation funnel.
set.seed(20260820L)
w_g <- factor(rep(seq_len(28L), each = 2L))
d_weak <- data.frame(y = rnorm(56L, 0, 1), x = rnorm(56L), g = w_g)
f_g <- factor(rep(c("c1", "c2", "c3"), each = 8L))
d_funnel <- data.frame(
  y = c(rnorm(8L, 0, 0.02), rnorm(8L, 0, 0.02), rnorm(8L, 0, 0.02)),
  x = rnorm(24L), g = f_g
)

# Space-level fixture (C4/FS-26): a factor level containing a space --
# the National Barley Agronomy tables NAPPLIED2 column ("nil N", "low
# N", "high N") is the real-world trigger. Two crossed random-intercept
# groups (TRIAL, N_TRIAL) at 4 and 3 levels with 8 reps per cell and a
# real fixed + random signal baked into the response (not pure noise):
# a 3-level x 2-level version of this same shape hit an INLA
# hyperparameter grid-search numerical instability of its own (a GSL
# interpolation abort, unrelated to the level-name defect this fixture
# exists to exercise) on several seeds. This size and signal is stable
# across five spot-checked seeds and three repeated runs of the grid
# seed used below.
set.seed(14L)
.sp_trial_u <- stats::rnorm(4L, 0, 0.4)
.sp_ntrial_u <- stats::rnorm(3L, 0, 0.3)
d_space_level <- expand.grid(
  NAPPLIED2 = factor(c("nil N", "low N", "high N")),
  TRIAL = factor(seq_len(4L)),
  N_TRIAL = factor(seq_len(3L)),
  rep = seq_len(8L)
)
.sp_napplied_eff <- c("high N" = 1.2, "low N" = 0.6, "nil N" = 0)
d_space_level$GRAIN_YIELD_THA <- 3 +
  .sp_napplied_eff[as.character(d_space_level$NAPPLIED2)] +
  .sp_trial_u[as.integer(d_space_level$TRIAL)] +
  .sp_ntrial_u[as.integer(d_space_level$N_TRIAL)] +
  stats::rnorm(nrow(d_space_level), 0, 0.4)

# Per-trial separable AR1 field fixture (C5/FS-27): three trials, each a
# complete 6 x 8 lattice, row / col labels reused identically across
# trials (the INLA replicate = mechanism expects one field per replicate
# of the SAME node count, not per-trial-distinct labels).
.ar1_mat_fn <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))
.at_trial_field <- function(seed) {
  set.seed(seed)
  grid <- expand.grid(col = seq_len(8L), row = seq_len(6L))
  Lc <- t(chol(1 * kronecker(.ar1_mat_fn(6L, 0.5), .ar1_mat_fn(8L, 0.3))))
  u <- as.numeric(Lc %*% stats::rnorm(48L))
  data.frame(
    y = 5 + u + stats::rnorm(48L, 0, 0.5),
    row = factor(grid$row), col = factor(grid$col)
  )
}
d_ar1_at_trial <- do.call(rbind, lapply(seq_len(3L), function(i) {
  dd <- .at_trial_field(1000L + i)
  dd$trial <- factor(paste0("T", i))
  dd
}))
'

.preamble <- paste0(
  '.libPaths(c("', .lib, '", .libPaths()))\n',
  .preamble_body
)

# --- cell constructor -------------------------------------------------- #

.cells <- list()

.add_cell <- function(id, section, code, expected, expect_src,
                      backend = NA_character_, family = NA_character_,
                      variant = NA_character_, needle = NA_character_,
                      matrix_cell = NA_character_) {
  .cells[[length(.cells) + 1L]] <<- list(
    cell_id = id,
    section = section,
    backend = backend,
    family = family,
    variant = variant,
    matrix_cell = matrix_cell,
    code = code,
    expected = expected,
    expect_src = expect_src,
    needle = needle
  )
  invisible(NULL)
}

# Reachability budget for every sampler cell: one chain, 300 draws.
.brms_args <- "chains = 1L, n_samples = 300L, warmup = 150L, seed = 1L"

.engine_args <- function(backend) {
  if (identical(backend, "brms")) paste0(", ", .brms_args) else ""
}

# --- M: the capability matrix, cell by cell ---------------------------- #
#
# The matrix is read from the INSTALLED package, so the expectations are
# the shipped table's own verdicts rather than a transcription of them.
# `.M_CALLS` maps each model class to the minimal live call that
# exercises it; a matrix row with no entry here, or an entry naming a
# row the matrix does not carry, stops the harness. The matrix is the
# claim surface and this file is the execution -- neither may drift
# past the other silently.

.matrix <- flexyBayes:::.fb_capability_matrix()

.M_CALLS <- list(
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
      "Combined interaction random effects and heterogeneous residual",
      "(full MET)"
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
    model_class = "Exact sufficient-statistic aggregation",
    inla = paste0('flexybayes(y ~ env, random = ~ gen, data = d_agg, ',
                  'aggregate = TRUE)'),
    brms = paste0('flexybayes(y ~ env, random = ~ gen, data = d_agg, ',
                  'aggregate = TRUE)')
  )
)

.m_named <- vapply(.M_CALLS, function(z) z$model_class, character(1L))
.m_missing <- setdiff(.matrix$model_class, .m_named)
.m_extra <- setdiff(.m_named, .matrix$model_class)
if (length(.m_missing) > 0L || length(.m_extra) > 0L) {
  stop(
    "The capability matrix and the grid's M section have drifted apart.\n",
    if (length(.m_missing) > 0L) {
      paste0("  matrix rows with no grid cell: ",
             paste(.m_missing, collapse = " | "), "\n")
    },
    if (length(.m_extra) > 0L) {
      paste0("  grid cells naming no matrix row: ",
             paste(.m_extra, collapse = " | "), "\n")
    },
    "Every capability cell is executed, or the gate is not F28.",
    call. = FALSE
  )
}

# The matrix verdict vocabulary, mapped to an executable expectation.
#
#   fits    -> the call completes and returns a fit
#   emits   -> the same, at emit level (no verdict in the table uses it
#              today; kept so a future `emits` row does not fall through)
#   refuses -> the call raises, with a flexybayes_* condition class
#   n/a     -> the class does not apply to that engine's interface; the
#              executable reading is that the request does not silently
#              produce the structure, so a typed refusal is expected
.m_expectation <- function(verdict) {
  switch(
    verdict,
    fits = "fit",
    emits = "fit",
    refuses = "refuse_typed",
    "n/a" = "refuse_typed",
    stop("unknown capability verdict: ", verdict, call. = FALSE)
  )
}

.m_slug <- function(x) {
  s <- tolower(gsub("[^a-z0-9]+", "-", tolower(x)))
  gsub("(^-|-$)", "", s)
}

.matrix_cells <- data.frame(
  cell = character(0), model_class = character(0),
  spelling = character(0), backend = character(0),
  verdict = character(0), stringsAsFactors = FALSE
)

for (spec in .M_CALLS) {
  row <- .matrix[.matrix$model_class == spec$model_class, , drop = FALSE]
  for (be in c("inla", "brms")) {
    verdict <- row[[be]]
    key <- paste0("cap-", .m_slug(spec$model_class), "-", be)
    # A row whose call already names a family keeps it; every other row
    # runs gaussian. The backend and the reachability budget are
    # appended to the call's own closing parenthesis.
    tail_args <- paste0(
      if (grepl('family = "', spec[[be]], fixed = TRUE)) "" else
        ', family = "gaussian"',
      ', backend = "', be, '"', .engine_args(be), ")"
    )
    code <- sub(")$", tail_args, spec[[be]])
    .add_cell(
      id = paste0("M-", .m_slug(spec$model_class), "-", be),
      section = "M",
      code = code,
      expected = .m_expectation(verdict),
      expect_src = paste0(
        "R/capability_matrix.R .fb_capability_matrix(): row '",
        spec$model_class, "', column ", be, " = ", verdict,
        "; anchor ", row$anchor
      ),
      backend = be,
      family = if (grepl('family = "', spec[[be]], fixed = TRUE)) {
        sub('.*family = "([^"]+)".*', "\\1", spec[[be]])
      } else {
        "gaussian"
      },
      variant = .m_slug(spec$model_class),
      matrix_cell = key
    )
    .matrix_cells <- rbind(.matrix_cells, data.frame(
      cell = key,
      model_class = spec$model_class,
      spelling = row$spelling,
      backend = be,
      verdict = verdict,
      stringsAsFactors = FALSE
    ))
  }
}

# --- M (hand-added): a factor level containing a space (C4/FS-26) ------ #
#
# Not derived from the capability matrix (no `matrix_cell`, so the M/
# matrix drift-check above does not expect a corresponding row) --
# this is a specific regression scenario, not a model-class cross
# product cell. `aggregate = FALSE` isolates the per-row emit_inla()
# fix this cell exercises: the same space-level design is also
# aggregation-eligible, and the aggregated INLA route (a different
# code path, R/emit_gaussian_aggregated.R) carries the identical
# untyped defect and has not been fixed in this slice -- see the WP-C
# report's Handoffs. `backend = "auto"` under the fix now correctly
# reaches INLA on this shape where it fell through silently before
# (FS-26's "Class" note: "a user with spaces in a factor level gets
# brms silently where INLA was available and appropriate").
.add_cell(
  id = "M-gaussian-lmm-space-level-inla",
  section = "M",
  code = paste0(
    'flexybayes(GRAIN_YIELD_THA ~ NAPPLIED2, random = ~ TRIAL + N_TRIAL, ',
    'data = d_space_level, family = "gaussian", backend = "inla", ',
    'aggregate = FALSE', .engine_args("inla"), ')'
  ),
  expected = "fit",
  expect_src = paste0(
    "R/methods_inla.R roxygen, coef.flexybayes_inla: 'A factor with a ",
    "non-syntactic level ... is legalised internally (make.names()) ",
    "before the INLA emit ... this method restores the user's own, ",
    "unlegalised label in the names it returns, as do ranef(), ",
    "summary() and predict(classify = ) on the same fit.' (C4/FS-26)"
  ),
  backend = "inla",
  family = "gaussian",
  variant = "space-level-factor"
)
.add_cell(
  id = "M-gaussian-lmm-space-level-brms",
  section = "M",
  code = paste0(
    'flexybayes(GRAIN_YIELD_THA ~ NAPPLIED2, random = ~ TRIAL + N_TRIAL, ',
    'data = d_space_level, family = "gaussian", backend = "brms"',
    .engine_args("brms"), ')'
  ),
  expected = "fit",
  expect_src = paste0(
    "FS-26 control: the identical model fits on brms regardless of the ",
    "space in the level (brms/Stan formula construction never builds a ",
    "bare symbol from a variable name plus a factor level, so the ",
    "defect was INLA-path-specific throughout)."
  ),
  backend = "brms",
  family = "gaussian",
  variant = "space-level-factor"
)

# --- M (hand-added): a per-trial separable AR1 field (C5/FS-27) -------- #
#
# Not derived from the capability matrix (no `matrix_cell`) -- a
# specific regression scenario. at(trial):ar1(row):ar1(col) lowers to
# INLA only (one separable field per level of `trial`, shared
# hyperparameters via replicate =); brms has no lowering for it at all
# (the pre-existing single-field ar1_spatial refusal, unchanged).
.add_cell(
  id = "M-nested-ar1-field-inla",
  section = "M",
  code = paste0(
    'flexybayes(y ~ 1, random = ~ at(trial):ar1(row):ar1(col), ',
    'data = d_ar1_at_trial, family = "gaussian", backend = "inla"',
    .engine_args("inla"), ')'
  ),
  expected = "fit",
  expect_src = paste0(
    "R/refusal_taxonomy.R, at_field_per_level_hyper_not_representable ",
    "description: 'The supported spelling, at(trial):ar1(row):ar1(col) ",
    "with no level argument, fits one field per level of `trial` on ",
    "INLA via the replicate = mechanism, sharing rho_row, rho_col and ",
    "the field variance across every level ...' (FS-27/C5)"
  ),
  backend = "inla",
  family = "gaussian",
  variant = "nested-ar1-field"
)
.add_cell(
  id = "M-nested-ar1-field-brms",
  section = "M",
  code = paste0(
    'flexybayes(y ~ 1, random = ~ at(trial):ar1(row):ar1(col), ',
    'data = d_ar1_at_trial, family = "gaussian", backend = "brms"',
    .engine_args("brms"), ')'
  ),
  expected = "refuse_typed",
  expect_src = paste0(
    "R/emit_brms.R ar1 / ar1_spatial route: 'An autoregressive latent ",
    "field is emitted by INLA only.' Unchanged by C5 -- the per-trial ",
    "spelling reuses the ar1_spatial term type this refusal already ",
    "keys on (flexybayes_refusal_stan_cannot_represent_ar1_field)."
  ),
  backend = "brms",
  family = "gaussian",
  variant = "nested-ar1-field"
)

# --- S1: prior route x family x backend -------------------------------- #
#
# Claim surface: the capability matrix's simple-random-intercept rows
# (both engines fit every family in this grid) crossed with the four
# documented prior routes plus the auto default. The brms arm reads the
# request back out of `brms::prior_summary()` and the INLA arm out of
# `prior_summary()`: acceptance without arrival was three separate
# silent-wrong defects in the 0.9.1 field engagement.

.s1_routes <- list(
  list(key = "default", args = "", needle = "uniform\\(0,",
       src = "R/flexybayes.R default_prior_active bounded-uniform default"),
  list(key = "vc_sd", args = ", prior_vc_sd = 3",
       needle = "lognormal\\(0, 3\\)",
       src = "man/flexybayes.Rd legacy lognormal(0, prior_vc_sd)"),
  list(key = "fixed_sd", args = ", prior_fixed_sd = 7",
       needle = "normal\\(0, 7\\)",
       src = "man/flexybayes.Rd SD for fixed-effect normal priors"),
  list(key = "fb_full",
       args = paste0(
         ", prior = fb_prior(sigma ~ half_normal(scale = 2), ",
         'sd(group = "g") ~ half_normal(scale = 1.5), ',
         'b("x") ~ normal(mean = 0, sd = 4))'
       ),
       needle = "normal\\(0, 1.5\\)",
       src = "man/fb_prior.Rd full PC-canonical specification"),
  list(key = "fb_partial",
       args = ', prior = fb_prior(sd(group = "g") ~ half_normal(scale = 1.5))',
       needle = "normal\\(0, 1.5\\)",
       src = "man/fb_prior.Rd partial specification, rest defaulted")
)

.s1_families <- c("gaussian", "poisson", "negative_binomial", "binomial",
                  "gamma", "beta")

for (rt in .s1_routes) {
  for (fam in .s1_families) {
    for (be in c("inla", "brms")) {
      .add_cell(
        id = sprintf("S1-%s-%s-%s", rt$key, fam, be),
        section = "S1",
        code = sprintf(
          paste0('flexybayes(y ~ x, random = ~ g, data = d_%s, ',
                 'family = "%s", backend = "%s"%s%s)'),
          fam, fam, be, .engine_args(be), rt$args
        ),
        expected = "fit",
        expect_src = paste0(
          "capability matrix: the simple random-effect GLMM row fits on ",
          "both engines; route source: ", rt$src
        ),
        backend = be, family = fam, variant = rt$key, needle = rt$needle
      )
    }
  }
}

# --- S2 / S2F: the malformed-prior corpus ------------------------------ #
#
# Claim surface: man/fb_prior.Rd -- calls outside the supported set raise
# a structured error naming the supported list. `expected` is "error"
# where the documented contract covers the malformation and "accept" for
# the deliberate controls. Append-only: every field report's
# reproduction joins this corpus permanently (recipe 62 G3).

.s2_cases <- list(
  list("wrong_arg_scale", 'sd(group = "g") ~ normal(scale = 2)', "error",
       "field finding D2; man/fb_prior.Rd normal(mean =, sd =)"),
  list("wrong_arg_mu_sigma", 'sd(group = "g") ~ normal(mu = 0, sigma = 2)',
       "error", "field finding D2; man/fb_prior.Rd normal(mean =, sd =)"),
  list("positional_normal", 'sd(group = "g") ~ normal(0, 2)', "accept",
       "R/fb_prior.R .name_prior_args names positional args (control)"),
  list("wrong_arg_hn_sd", 'sd(group = "g") ~ half_normal(sd = 1.5)', "error",
       "R/fb_prior.R half_normal parameter is `scale`"),
  list("extra_arg_hn", 'sd(group = "g") ~ half_normal(scale = 2, extra = 9)',
       "error", "R/fb_prior.R half_normal takes `scale` only"),
  list("unknown_dist_normla", 'sd(group = "g") ~ normla(0, 2)', "error",
       "R/fb_prior.R unsupported prior distribution refusal"),
  list("unknown_dist_cauchy2", 'sd(group = "g") ~ cauchy2(0, 2)', "error",
       "R/fb_prior.R unsupported prior distribution refusal"),
  list("negative_scale", 'sd(group = "g") ~ half_normal(scale = -1)', "error",
       "sd scale must be positive (the DSL is on the SD scale)"),
  list("zero_scale", 'sd(group = "g") ~ half_normal(scale = 0)', "error",
       "sd scale must be positive (a point mass otherwise)"),
  list("zerolen_scale", 'sd(group = "g") ~ half_normal(scale = numeric(0))',
       "error", "scalar hyperparameter required"),
  list("len2_scale", 'sd(group = "g") ~ half_normal(scale = c(1, 2))',
       "error", "scalar hyperparameter required"),
  list("na_scale", 'sd(group = "g") ~ half_normal(scale = NA)', "error",
       "scalar finite hyperparameter required"),
  list("inf_scale", 'sd(group = "g") ~ half_normal(scale = Inf)', "error",
       "scalar finite hyperparameter required"),
  list("char_scale", 'sd(group = "g") ~ half_normal(scale = "1")', "error",
       "numeric hyperparameter required"),
  list("bad_b_target", 'b("nonexistent_term") ~ normal(mean = 0, sd = 2)',
       "accept", "target mismatch surfaces at fit entry, not construction"),
  list("bad_sd_group",
       'sd(group = "nonexistent_group") ~ half_normal(scale = 1)', "accept",
       "target mismatch surfaces at fit entry, not construction"),
  list("bare_unknown_target", 'phi ~ half_normal(scale = 1)', "error",
       "R/fb_prior.R supported targets: sigma, sd(), b(), cor(), smooth()"),
  list("pc_negative_upper", 'sigma ~ pc(upper = -1, prob = 0.05)', "error",
       "PC upper is an SD-scale bound and must be positive"),
  list("pc_prob_gt1", 'sigma ~ pc(upper = 1, prob = 1.5)', "error",
       "PC prob is a probability"),
  list("pc_missing_prob", 'sigma ~ pc(upper = 1)', "error",
       "R/fb_prior.R pc parameters are (upper, prob)"),
  list("pc_missing_upper", 'sigma ~ pc(prob = 0.05)', "error",
       "R/fb_prior.R pc parameters are (upper, prob)"),
  list("t_negative_df", 'sigma ~ student_t(df = -3, location = 0, scale = 2)',
       "error", "degrees of freedom must be positive"),
  list("t_single_positional", 'sigma ~ student_t(2.5)', "error",
       "R/fb_prior.R student_t parameters are (df, location, scale)"),
  list("exp_negative_rate", 'sigma ~ exponential(rate = -1)', "error",
       "exponential rate must be positive"),
  list("gamma_missing_rate", 'sigma ~ gamma(shape = 1)', "error",
       "R/fb_prior.R gamma parameters are (shape, rate)"),
  list("lkj_negative_eta", 'cor(group = "g") ~ lkj(eta = -1)', "error",
       "LKJ eta must be positive"),
  list("uniform_reversed", 'sigma ~ uniform(lower = 2, upper = 1)', "error",
       "R/fb_prior.R uniform upper must exceed lower (control: is caught)"),
  list("hc_extra_location", 'sigma ~ half_cauchy(location = 0, scale = 2)',
       "error", "R/fb_prior.R half_cauchy takes `scale` only"),
  list("dup_named_arg", 'sigma ~ normal(mean = 0, sd = 2, sd = 3)', "error",
       "duplicate formal in a prior call"),
  list("sd_positional_group", 'sd("g") ~ half_normal(scale = 1)', "accept",
       "R/fb_prior.R .extract_string_arg accepts positional (control)"),
  list("normal_reversed_named", 'sd(group = "g") ~ normal(sd = 2, mean = 0)',
       "accept", "named arguments in any order (control)"),
  list("b_unquoted_target", 'b(x) ~ normal(mean = 0, sd = 2)', "accept",
       "R/fb_prior.R b() deparses its first argument (control)")
)

for (cs in .s2_cases) {
  .add_cell(
    id = paste0("S2-", cs[[1]]), section = "S2",
    code = sprintf("fb_prior(%s)", cs[[2]]),
    expected = cs[[3]], expect_src = cs[[4]], variant = cs[[1]]
  )
}

.s2_shape <- list(
  list("string_not_formula", 'fb_prior("sigma ~ half_normal(scale = 1)")',
       "error", "man/fb_prior.Rd: arguments are two-sided formulas"),
  list("one_sided_formula", "fb_prior(~ half_normal(scale = 1))", "error",
       "R/fb_prior.R two-sided formula required"),
  list("no_args", "fb_prior()", "error",
       "R/fb_prior.R at least one specification required"),
  list("target_not_call", "fb_prior(sigma ~ half_normal)", "error",
       "R/fb_prior.R the distribution must be a call"),
  list("list_arg", "fb_prior(list(sigma ~ half_normal(scale = 1)))", "error",
       "R/fb_prior.R each argument must be a formula")
)

for (cs in .s2_shape) {
  .add_cell(
    id = paste0("S2-", cs[[1]]), section = "S2", code = cs[[2]],
    expected = cs[[3]], expect_src = cs[[4]], variant = cs[[1]]
  )
}

# S2F: one fit per corpus entry. A construction the DSL legitimately
# accepts and a fit that must refuse are not in conflict -- the two
# controls whose target names a term the model does not have refuse at
# fit entry (FS-20), and the three shape controls fit.
.s2f_fit_expected <- c("positional_normal", "normal_reversed_named",
                       "sd_positional_group", "b_unquoted_target")

for (cs in .s2_cases) {
  expected <- if (cs[[1]] %in% .s2f_fit_expected) "fit" else "error"
  .add_cell(
    id = paste0("S2F-", cs[[1]]), section = "S2F",
    code = sprintf(
      paste0('flexybayes(y ~ x, random = ~ g, data = d_gaussian, ',
             'family = "gaussian", backend = "brms", %s, ',
             "prior = fb_prior(%s))"),
      .brms_args, cs[[2]]
    ),
    expected = expected,
    expect_src = paste0(
      "field finding D2 / FS-20: an accepted malformed prior must not ",
      "reach a completed fit under defaults, and a prior naming a term ",
      "the model does not have refuses at fit entry. Construction ",
      "expectation: ", cs[[4]]
    ),
    backend = "brms", family = "gaussian", variant = cs[[1]],
    needle = "___route_specific___"
  )
}

# --- S3a: allowlist against the installed engines' rosters ------------- #
#
# G5. The rosters are read from the INSTALLED engines, never recalled:
# brms families from its own `.family_*` constructor set, INLA
# likelihoods from `names(INLA::inla.models()$likelihood)`.

.add_cell(
  id = "S3a-engine-rosters", section = "S3a",
  code = paste0(
    "{\n",
    '  bf <- ls(asNamespace("brms"), all.names = TRUE)\n',
    '  brms_fams <- sub("^[.]family_", "",\n',
    '                   grep("^[.]family_", bf, value = TRUE))\n',
    '  brms_fams <- setdiff(brms_fams, c("info", "custom"))\n',
    "  inla_fams <- names(INLA::inla.models()$likelihood)\n",
    # The entry allowlist is the `defaults` list inside .resolve_family();
    # it is transcribed here with its source named, and the S3a probe
    # cells below execute every family it leaves out.
    "  fb_allow <- c(\"gaussian\", \"binomial\", \"binary\", \"poisson\",\n",
    "    \"negative_binomial\", \"negbinom\", \"gamma\", \"beta\",\n",
    "    \"hurdle_gamma\")\n",
    "  writeLines(c(\n",
    '    "# flexyBayes entry allowlist (R/utils.R .resolve_family)",\n',
    '    paste(fb_allow, collapse = ", "), "",\n',
    '    paste0("# brms family roster, installed brms ",\n',
    '           as.character(utils::packageVersion("brms")),\n',
    '           " (.family_* constructors, n = ", length(brms_fams), ")"),\n',
    '    paste(sort(brms_fams), collapse = ", "), "",\n',
    '    paste0("# INLA likelihood roster, installed INLA ",\n',
    '           as.character(utils::packageVersion("INLA")),\n',
    '           " (names(inla.models()$likelihood), n = ",\n',
    '           length(inla_fams), ")"),\n',
    '    paste(sort(inla_fams), collapse = ", ")\n',
    "  ), FILE_ENGINE_FAMILIES)\n",
    '  paste0("brms=", length(brms_fams), " inla=", length(inla_fams))\n',
    "}"
  ),
  expected = "artefact",
  expect_src = "recipe 62 G5: rosters read live from the installed engines",
  variant = "rosters"
)

# Families an installed engine carries natively but the entry allowlist
# refuses. Each is expected to refuse; `roster_diff.csv` then records
# whose boundary each refusal is.
#
# `Gamma` is deliberately NOT in this list. `.resolve_family()` lowercases
# before matching, so the base-R spelling is ADMITTED while the refusal
# message's allowlist is lower-case only -- the accepted set is wider than
# the documented one (recorded in the P0 and P3a reports). It gets its own
# fitting cells below, on a response the family can take.
.s3a_probe_families <- c(
  "tweedie", "compound_poisson", "zero_inflated_gamma",
  "hurdle_poisson", "hurdle_negbinomial", "zero_inflated_poisson",
  "zero_inflated_negbinomial", "zero_inflated_binomial", "beta_binomial",
  "student", "lognormal", "weibull", "exponential", "geometric",
  "inverse.gaussian", "skew_normal", "von_mises", "bernoulli",
  "nbinomial", "betabinomial", "gpoisson", "xbinomial"
)

for (fam in .s3a_probe_families) {
  .add_cell(
    id = paste0("S3a-refuse-", fam), section = "S3a",
    code = sprintf(
      paste0('flexybayes(y ~ x, random = ~ g, data = d_gaussian, ',
             'family = "%s", backend = "brms", %s)'),
      fam, .brms_args
    ),
    expected = "refuse_typed",
    expect_src = paste0(
      "R/utils.R .resolve_family unsupported_family refusal; DESCRIPTION ",
      "'typed refusals naming the nearest implemented alternative'"
    ),
    backend = "brms", family = fam, variant = "allowlist_probe"
  )
}

# --- S3b / S3c: family spellings and the C2 families ------------------- #

.s3b_families <- c("gaussian", "binomial", "binary", "poisson",
                   "negative_binomial", "negbinom", "gamma", "beta")

for (fam in .s3b_families) {
  for (be in c("inla", "brms")) {
    .add_cell(
      id = sprintf("S3b-%s-%s", fam, be), section = "S3b",
      code = sprintf(
        paste0('flexybayes(y ~ x, random = ~ g, data = d_%s, ',
               'family = "%s", backend = "%s"%s)'),
        fam, fam, be, .engine_args(be)
      ),
      expected = "fit",
      expect_src = paste0(
        "R/utils.R .resolve_family admits this spelling; the capability ",
        "matrix's GLMM row fits on both engines"
      ),
      backend = be, family = fam, variant = "allowed_family"
    )
  }
}

# The base-R `Gamma` spelling, which the entry allowlist admits through
# `tolower()` while the refusal message lists only `gamma`. The cells
# assert what the code does rather than what the message implies, so the
# gap stays visible until the documented and accepted sets are one set.
for (be in c("inla", "brms")) {
  .add_cell(
    id = sprintf("S3b-Gamma-%s", be), section = "S3b",
    code = sprintf(
      paste0('flexybayes(y ~ x, random = ~ g, data = d_gamma, ',
             'family = "Gamma", backend = "%s"%s)'),
      be, .engine_args(be)
    ),
    expected = "fit",
    expect_src = paste0(
      "R/utils.R .resolve_family lowercases before matching, so the ",
      "base-R spelling is admitted; the refusal message's allowlist is ",
      "lower-case only, which is the documented-versus-accepted gap the ",
      "P0 and P3a reports carry forward"
    ),
    backend = be, family = "Gamma", variant = "allowed_family"
  )
}

.s3c_families <- list(
  list("hurdle_gamma", "fit", "refuse_typed"),
  list("tweedie", "refuse_typed", "refuse_typed"),
  list("zero_inflated_gamma", "refuse_typed", "refuse_typed"),
  list("compound_poisson", "refuse_typed", "refuse_typed")
)

for (fm in .s3c_families) {
  for (be in c("inla", "brms")) {
    .add_cell(
      id = sprintf("S3c-%s-%s", fm[[1]], be), section = "S3c",
      code = sprintf(
        paste0('flexybayes(y ~ x, random = ~ g, data = d_hurdle_gamma, ',
               'family = "%s", backend = "%s"%s)'),
        fm[[1]], be, .engine_args(be)
      ),
      expected = if (identical(be, "brms")) fm[[2]] else fm[[3]],
      expect_src = paste0(
        "field finding C2 read against the installed rosters: ",
        "hurdle_gamma is brms-native; the other three are not, and INLA ",
        "carries none of the four (inst/KNOWN_ISSUES.md)"
      ),
      backend = be, family = fm[[1]], variant = "c2_family"
    )
  }
}

# --- S4: structure x family x backend ---------------------------------- #
#
# The matrix cells run at gaussian in section M; this section carries the
# same structures on a count response, plus the structures whose claim
# lives in a refusal contract rather than in a matrix row.

.s4_structures <- list(
  list(key = "simple", call = "random = ~ g", inla = "fit", brms = "fit",
       src = "capability matrix: simple random intercept fits on both"),
  list(key = "two_terms", call = "random = ~ g + h", inla = "fit",
       brms = "fit",
       src = "capability matrix: two simple random intercepts, both fit"),
  list(key = "us", call = "random = ~ us(f):g", inla = "refuse_typed",
       brms = "fit", src = "capability matrix: ~ us(f):g refuses | fits"),
  list(key = "diag", call = "random = ~ diag(f):g", inla = "refuse_typed",
       brms = "fit", src = "capability matrix: ~ diag(f):g refuses | fits"),
  list(key = "at", call = "random = ~ at(f):g", inla = "refuse_typed",
       brms = "fit",
       src = "capability matrix: at() is a diag() spelling, refuses | fits"),
  list(key = "dsum_residual",
       call = "random = ~ g, residual = ~ dsum(~ units | f)",
       inla = "refuse_typed", brms = "fit",
       # The matrix row's own note: "Refused for families with no
       # residual scale." Poisson has none, so the count arm refuses on
       # both engines and the refusal is the claim, not a miss.
       poisson = list(inla = "refuse_typed", brms = "refuse_typed"),
       src = paste0("capability matrix: heterogeneous residual refuses | ",
                    "fits, and the row's note refuses it for a family ",
                    "with no residual scale")),
  list(key = "ar1", call = "random = ~ ar1(row):ar1(col)", inla = "fit",
       brms = "refuse_typed",
       src = "capability matrix: separable AR1 field fits | refuses"),
  list(key = "spl", call = "random = ~ spl(x)", inla = "fit",
       brms = "refuse_typed",
       # `x` carries the smooth, so it must not also be a fixed term:
       # INLA refuses a key used twice (FS-17), which is why the matrix
       # row's own anchor writes the fixed part as an intercept.
       fixed = "1",
       src = "capability matrix: univariate P-spline fits | refuses"),
  list(key = "smooth_s", call = "fixed_override = y ~ s(x)",
       inla = "refuse_typed", brms = "refuse_typed",
       src = paste0("FS-18: s(x) is the mgcv spelling, so the refusal is ",
                    "typed on both engines and names spl(x)")),
  list(key = "met_combined",
       call = "random = ~ g + g:f, residual = ~ dsum(~ units | f)",
       inla = "refuse_typed", brms = "fit",
       # Carries the dsum residual, so the count arm refuses for the same
       # reason the dsum_residual row above does.
       poisson = list(inla = "refuse_typed", brms = "refuse_typed"),
       src = paste0("capability matrix: full MET refuses | fits, with the ",
                    "heterogeneous-residual row's no-residual-scale note")),
  list(key = "corh", call = "random = ~ corh(f):g", inla = "refuse_typed",
       brms = "refuse_typed",
       src = "capability matrix: ~ corh(f):g refuses | refuses"),
  list(key = "fa", call = "random = ~ fa(f, 1):g", inla = "refuse_typed",
       brms = "refuse_typed",
       src = "capability matrix: ~ fa(env, k):gen refuses | refuses"),
  list(key = "nested", call = "random = ~ g:h", inla = "refuse_typed",
       brms = "fit",
       src = "capability matrix: nested / interaction refuses | fits"),
  list(key = "weights", call = "random = ~ g, weights = runif(60, 0.5, 2)",
       inla = "fit", brms = "fit",
       # C6: weights are lowered for family = "gaussian" with the
       # identity link only. The count arm (poisson, ycount) is exactly
       # the family this row's own note excludes, so it refuses by name
       # (weights_requires_gaussian) on both engines rather than fitting
       # something else.
       poisson = list(inla = "refuse_typed", brms = "refuse_typed"),
       src = paste0("capability matrix: Observation weights (Gaussian, ",
                    "identity link) fits | fits; refused by name off ",
                    "that family (weights_requires_gaussian)"))
)

for (st in .s4_structures) {
  for (be in c("inla", "brms")) {
    for (fam in c("gaussian", "poisson")) {
      resp <- if (identical(fam, "gaussian")) "y" else "ycount"
      if (grepl("^fixed_override", st$call)) {
        fixed_txt <- sub("^fixed_override = y", resp, st$call)
        rest <- ""
      } else {
        fixed_txt <- paste0(resp, " ~ ", st$fixed %||% "x")
        rest <- paste0(", ", st$call)
      }
      # A per-family override where the structure's verdict depends on
      # the response family rather than on the engine alone.
      expected <- if (!is.null(st[[fam]])) st[[fam]][[be]] else st[[be]]
      .add_cell(
        id = sprintf("S4-%s-%s-%s", st$key, fam, be), section = "S4",
        code = sprintf(
          paste0('flexybayes(%s, data = d_struct%s, family = "%s", ',
                 'backend = "%s"%s)'),
          fixed_txt, rest, fam, be, .engine_args(be)
        ),
        expected = expected, expect_src = st$src,
        backend = be, family = fam, variant = st$key
      )
    }
  }
}

# --- S5: grammar-surface parity ---------------------------------------- #
#
# Bar-grammar forms and the ASReml surfaces of the same models. Where a
# bar form refuses and its ASReml sibling fits, the refusal has to name
# the sibling; a refusal that names nothing is a remedy-free refusal
# (miss class M5) and an untyped error is M2.

.s5_forms <- list(
  list("bar_intercept", "y ~ x + (1 | g)", NA, "fit", "fit",
       "capability matrix: (1 | g) fits | fits"),
  list("bar_slope_cor", "y ~ x + (x | g)", NA, "refuse_typed",
       "refuse_typed", "capability matrix: (x | g) refuses | refuses"),
  list("bar_slope_uncor", "y ~ x + (x || g)", NA, "refuse_typed", "fit",
       "capability matrix: (x || g) refuses | fits"),
  list("bar_factor_slope", "y ~ trt + (trt | g)", NA, "refuse_typed",
       "refuse_typed",
       "FS-6 / U1: typed refusal naming us(trt):g as the ASReml surface"),
  list("bar_nested", "y ~ x + (1 | g:h)", NA, "refuse_typed", "fit",
       "capability matrix: nested / interaction refuses | fits"),
  list("bar_zero_factor_slope", "y ~ trt + (0 + trt | g)", NA,
       "refuse_typed", "refuse_typed",
       "FS-6 / U1: typed refusal naming us(trt):g"),
  list("bar_multi_slope_uncor", "y ~ x + z + (x + z || g)", NA,
       "refuse_typed", "refuse_typed",
       "FS-6: multi-variable uncorrelated slopes refuse typed, claiming no
        alternative"),
  list("asreml_simple", "y ~ x", "~ g", "fit", "fit",
       "capability matrix: random = ~ g fits | fits"),
  list("asreml_us_trt", "y ~ trt", "~ us(trt):g", "refuse_typed", "fit",
       "capability matrix: ~ us(f):g refuses | fits"),
  list("asreml_diag_trt", "y ~ trt", "~ diag(trt):g", "refuse_typed", "fit",
       "capability matrix: ~ diag(f):g refuses | fits"),
  list("asreml_nested", "y ~ x", "~ g:h", "refuse_typed", "fit",
       "capability matrix: nested / interaction refuses | fits"),
  list("fixed_factor_numeric", "y ~ f * x", NA, "refuse_typed", "fit",
       "capability matrix: factor-by-numeric fixed interaction refuses | fits"),
  list("fixed_factor_numeric_colon", "y ~ f:x", NA, "refuse_typed", "fit",
       "capability matrix: factor-by-numeric interaction refuses | fits"),
  list("random_factor_numeric", "y ~ x", "~ f:x", "refuse_typed",
       "refuse_typed",
       "field finding U2: a numeric in a random interaction is not a
        grouping factor, and the refusal names the remedy")
)

for (fm in .s5_forms) {
  for (be in c("inla", "brms")) {
    rnd <- if (is.na(fm[[3]])) "" else paste0(", random = ", fm[[3]])
    .add_cell(
      id = sprintf("S5-%s-%s", fm[[1]], be), section = "S5",
      code = sprintf(
        paste0('flexybayes(%s%s, data = d_struct, family = "gaussian", ',
               'backend = "%s"%s)'),
        fm[[2]], rnd, be, .engine_args(be)
      ),
      expected = if (identical(be, "inla")) fm[[4]] else fm[[5]],
      expect_src = gsub("[[:space:]]+", " ", fm[[6]]),
      backend = be, family = "gaussian", variant = fm[[1]]
    )
  }
}

# --- T: the prior-translation table, executed -------------------------- #
#
# One cell per (target, distribution, engine). The expectation is read
# from the package's own translation table and the child then asks the
# ENGINE whether the prior arrived. A completed fit whose prior never
# reached the engine argument is a silent drop -- the FS-21 class, and
# the reason a table alone is not evidence.

.t_specs <- list(
  normal      = "normal(mean = 0, sd = 2)",
  half_normal = "half_normal(scale = 1.5)",
  half_cauchy = "half_cauchy(scale = 1)",
  cauchy      = "cauchy(location = 0, scale = 1)",
  student_t   = "student_t(df = 3, scale = 1)",
  exponential = "exponential(rate = 2)",
  gamma       = "gamma(shape = 1, rate = 2)",
  lkj         = "lkj(eta = 2)",
  pc          = "pc(upper = 1, prob = 0.05)",
  uniform     = "uniform(lower = 0, upper = 5)"
)

.t_targets <- list(
  sigma  = "sigma",
  sd     = 'sd(group = "g")',
  b      = 'b("x")',
  cor    = 'cor(group = "g")',
  smooth = 'smooth("x", basis = "rw2")'
)

for (tg in names(.t_targets)) {
  for (dist in names(.t_specs)) {
    for (be in c("inla", "brms")) {
      # A smooth() prior needs a model that carries a smooth. The spline
      # itself refuses on brms, so the brms arm runs the ordinary model:
      # the prior TARGET is what is under test, and it refuses either way.
      code <- if (identical(tg, "smooth")) {
        sprintf(
          paste0('flexybayes(y ~ 1, random = ~ spl(x), data = d_struct, ',
                 'family = "gaussian", backend = "%s"%s, ',
                 "prior = fb_prior(%s ~ %s))"),
          be, .engine_args(be), .t_targets[[tg]], .t_specs[[dist]]
        )
      } else {
        sprintf(
          paste0('flexybayes(y ~ x, random = ~ g, data = d_gaussian, ',
                 'family = "gaussian", backend = "%s"%s, ',
                 "prior = fb_prior(%s ~ %s))"),
          be, .engine_args(be), .t_targets[[tg]], .t_specs[[dist]]
        )
      }
      .add_cell(
        id = sprintf("T-%s-%s-%s", tg, dist, be), section = "T",
        code = code, expected = "___from_table___",
        expect_src = paste0(
          "R/prior_translation.R: the (engine, target, distribution) ",
          "verdict, plus arrival in the engine's own argument list"
        ),
        backend = be, family = "gaussian",
        variant = paste0(tg, "|", dist),
        needle = paste0("___t___", tg)
      )
    }
  }
}

# Resolve the T expectations from the installed package's own table.
# Reading the table here rather than re-listing it means a table edit
# with no matching behaviour change surfaces as a divergent cell.
.t_verdict <- function(engine, target, dist) {
  flexyBayes:::.fb_prior_translation_lookup(engine, target, dist)
}

for (i in seq_along(.cells)) {
  if (!identical(.cells[[i]]$section, "T")) next
  parts <- strsplit(.cells[[i]]$variant, "|", fixed = TRUE)[[1L]]
  v <- .t_verdict(.cells[[i]]$backend, parts[[1L]], parts[[2L]])
  .cells[[i]]$expected <- if (identical(v, "translate")) {
    "translated"
  } else {
    "refuse_typed"
  }
  .cells[[i]]$expect_src <- paste0(.cells[[i]]$expect_src,
                                   "; table verdict = ", v)
}

# --- S6: diagnostic structuring ---------------------------------------- #

.add_cell(
  id = "S6-inla-weak", section = "S6",
  code = paste0('flexybayes(y ~ x, random = ~ g, data = d_weak, ',
                'family = "gaussian", backend = "inla")'),
  expected = "fit",
  expect_src = "field finding U3: a structured convergence record",
  backend = "inla", family = "gaussian", variant = "weak_identification"
)

.add_cell(
  id = "S6-inla-weak-poisson", section = "S6",
  code = paste0('flexybayes(y ~ x, random = ~ g, ',
                "data = transform(d_weak, y = rpois(nrow(d_weak), 0.05)), ",
                'family = "poisson", backend = "inla")'),
  expected = "fit",
  expect_src = "field finding U3: a zero-heavy weakly identified count part",
  backend = "inla", family = "poisson", variant = "weak_identification"
)

.add_cell(
  id = "S6-brms-funnel", section = "S6",
  code = paste0('flexybayes(y ~ x, random = ~ g, data = d_funnel, ',
                'family = "gaussian", backend = "brms", chains = 1L, ',
                "n_samples = 300L, warmup = 150L, seed = 1L, ",
                "control = list(adapt_delta = 0.5, max_treedepth = 5))"),
  expected = "fit",
  expect_src = "field finding U3: divergent transitions on the fit object",
  backend = "brms", family = "gaussian", variant = "divergence"
)

.add_cell(
  id = "S6-brms-funnel-tight", section = "S6",
  code = paste0('flexybayes(y ~ x, random = ~ g, data = d_funnel, ',
                'family = "gaussian", backend = "brms", chains = 1L, ',
                "n_samples = 300L, warmup = 50L, seed = 7L, ",
                "control = list(adapt_delta = 0.4, max_treedepth = 4))"),
  expected = "fit",
  expect_src = "field finding U3: divergent transitions on the fit object",
  backend = "brms", family = "gaussian", variant = "divergence"
)

# --- child worker ------------------------------------------------------ #
#
# One cell per process. Warnings and messages are captured and muffled,
# so a cell that only warns still reports as a completed fit.

.child_worker <- function(preamble, code, needle, file_engine_families) {
  FILE_ENGINE_FAMILIES <- file_engine_families
  env <- environment()
  eval(parse(text = preamble), envir = env)

  warn_msgs <- character()
  warn_cls <- character()
  note_msgs <- character()

  out <- withCallingHandlers(
    tryCatch(
      list(ok = TRUE, value = eval(parse(text = code), envir = env)),
      error = function(e) {
        list(ok = FALSE, cond_class = paste(class(e), collapse = "|"),
             message = conditionMessage(e))
      }
    ),
    warning = function(w) {
      warn_msgs <<- c(warn_msgs, conditionMessage(w))
      warn_cls <<- c(warn_cls, paste(class(w), collapse = "|"))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      note_msgs <<- c(note_msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  res <- list(
    warn_class = paste(unique(warn_cls), collapse = " ;; "),
    warn_message = paste(unique(warn_msgs), collapse = " ;; "),
    note_message = paste(unique(note_msgs), collapse = " ;; "),
    prior_table = NA_character_,
    prior_reached = NA,
    converge_json = NA_character_,
    diag_structured = NA
  )

  if (!isTRUE(out$ok)) {
    res$status <- "error"
    res$cond_class <- out$cond_class
    res$message <- out$message
    return(res)
  }

  value <- out$value
  res$cond_class <- paste(class(value), collapse = "|")
  res$message <- ""

  if (inherits(value, "fb_prior")) {
    res$status <- "construct_accept"
    return(res)
  }
  if (!inherits(value, c("flexybayes", "flexybayes_aggregated"))) {
    res$status <- "value"
    res$message <- paste(
      utils::capture.output(utils::str(value, max.level = 1L)),
      collapse = " | "
    )
    return(res)
  }

  res$status <- "fit"

  # Translation-grid arrival: ask the ENGINE, in the argument the
  # translation table says the prior arrives through.
  if (!is.na(needle) && startsWith(needle, "___t___")) {
    tgt <- sub("^___t___", "", needle)
    if (!is.null(value$brms)) {
      ps <- try(as.data.frame(brms::prior_summary(value$brms)),
                silent = TRUE)
      res$prior_reached <- if (inherits(ps, "try-error")) {
        NA
      } else {
        any(ps$source == "user")
      }
      if (!inherits(ps, "try-error")) {
        res$prior_table <- paste(
          paste0(ifelse(nzchar(ps$prior), ps$prior, "<flat>"), " @ ",
                 ps$class,
                 ifelse(nzchar(ps$group), paste0("/", ps$group), ""),
                 " [", ps$source, "]"),
          collapse = " ;; "
        )
      }
    } else {
      form <- paste(deparse(value$inla$.args$formula), collapse = "")
      cf <- value$inla$.args$control.fixed
      cfam <- value$inla$.args$control.family
      hyper_set <- FALSE
      if (!is.null(cfam) && length(cfam)) {
        h <- cfam[[1L]]$hyper
        # INLA stamps an inla.read.only attribute on each field; a slot
        # still carrying the engine default is not an arrival.
        hyper_set <- !is.null(h) && any(vapply(
          h,
          function(e) {
            pv <- as.character(e$prior %||% "")[[1L]]
            !pv %in% c("loggamma", "none", "")
          },
          logical(1L)
        ))
      }
      res$prior_reached <- switch(
        tgt,
        "sigma" = hyper_set,
        "sd" = grepl("hyper", form, fixed = TRUE),
        "smooth" = grepl("hyper", form, fixed = TRUE),
        "b" = !is.null(cf$mean) && length(cf$mean) > 0L,
        "cor" = FALSE,
        NA
      )
      res$prior_table <- substr(paste0(
        "formula=", form, " ;; control.fixed=",
        paste(names(cf), collapse = ","), " ;; family_hyper=", hyper_set
      ), 1L, 1200L)
    }
    return(res)
  }

  # S1 round-trip: the requested prior, read back out of the engine.
  if (!is.null(value$brms)) {
    ps <- try(as.data.frame(brms::prior_summary(value$brms)), silent = TRUE)
    if (!inherits(ps, "try-error")) {
      res$prior_table <- paste(
        paste0(ifelse(nzchar(ps$prior), ps$prior, "<flat>"), " @ ", ps$class,
               ifelse(nzchar(ps$coef), paste0(":", ps$coef), ""),
               ifelse(nzchar(ps$group), paste0("/", ps$group), ""),
               " [", ps$source, "]"),
        collapse = " ;; "
      )
      if (identical(needle, "___route_specific___")) {
        res$prior_reached <- any(
          ps$source == "user" & !grepl("^uniform\\(0,", ps$prior)
        )
      } else if (!is.na(needle)) {
        res$prior_reached <- any(grepl(needle, ps$prior))
      }
    }
  }

  if (is.null(value$brms)) {
    pr <- try(
      utils::capture.output(print(flexyBayes::prior_summary(value))),
      silent = TRUE
    )
    if (!inherits(pr, "try-error")) {
      res$prior_table <- paste(pr, collapse = " ;; ")
      if (!is.na(needle) && !identical(needle, "___route_specific___")) {
        res$prior_reached <- any(grepl(needle, pr))
      }
    }
  }

  sm <- try(summary(value), silent = TRUE)
  if (!inherits(sm, "try-error") && !is.null(sm$converge)) {
    res$converge_json <- paste(
      utils::capture.output(utils::str(sm$converge, max.level = 2L)),
      collapse = " | "
    )
    txt <- tolower(res$converge_json)
    res$diag_structured <- any(grepl(
      "diverg|collapse|warn|fail|not converged|unreliable|ess", txt
    ))
  }

  res
}

# --- dispatch ---------------------------------------------------------- #

.file_engine_families <- file.path(.raw_dir, "engine_families.txt")

.run_one <- function(cell) {
  t0 <- Sys.time()
  res <- tryCatch(
    callr::r(
      .child_worker,
      args = list(preamble = .preamble, code = cell$code,
                  needle = cell$needle,
                  file_engine_families = .file_engine_families),
      timeout = .cell_timeout, show = FALSE, spinner = FALSE
    ),
    error = function(e) {
      cls <- paste(class(e), collapse = "|")
      status <- if (grepl("timeout", cls, ignore.case = TRUE)) {
        "timeout"
      } else {
        "crash"
      }
      list(status = status, cond_class = cls, message = conditionMessage(e),
           warn_class = "", warn_message = "", note_message = "",
           prior_table = NA_character_, prior_reached = NA,
           converge_json = NA_character_, diag_structured = NA)
    }
  )
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[%s] %-46s %-16s %6.1fs\n",
              format(Sys.time(), "%H:%M:%S"), cell$cell_id, res$status,
              wall))
  utils::flush.console()

  data.frame(
    cell_id = cell$cell_id,
    section = cell$section,
    backend = cell$backend,
    family = cell$family,
    variant = cell$variant,
    matrix_cell = cell$matrix_cell,
    expected = cell$expected,
    expect_src = cell$expect_src,
    code = cell$code,
    observed = res$status,
    cond_class = res$cond_class %||% "",
    message = substr(gsub("[\r\n]+", " ", res$message %||% ""), 1L, 300L),
    warn_class = substr(res$warn_class %||% "", 1L, 200L),
    warn_message = substr(gsub("[\r\n]+", " ", res$warn_message %||% ""),
                          1L, 300L),
    note_message = substr(gsub("[\r\n]+", " ", res$note_message %||% ""),
                          1L, 200L),
    prior_reached = res$prior_reached %||% NA,
    prior_table = substr(res$prior_table %||% NA_character_, 1L, 1200L),
    converge_json = substr(
      gsub("[\r\n]+", " ", res$converge_json %||% NA_character_), 1L, 600L
    ),
    diag_structured = res$diag_structured %||% NA,
    wall_sec = round(wall, 2),
    stringsAsFactors = FALSE
  )
}

cat("cells to run:", length(.cells), "\n")
cat("cores:", .n_cores, "| per-cell timeout:", .cell_timeout, "s\n")
cat("out  :", .out_dir, "\nraw  :", .raw_dir, "\n\n")

if (nzchar(Sys.getenv("FB_GRID_DRYRUN"))) {
  for (cl in .cells) {
    ok <- tryCatch({
      parse(text = cl$code)
      "parses"
    }, error = function(e) paste("PARSE-FAIL:", conditionMessage(e)))
    cat(sprintf("%-46s %-8s %-14s %s\n", cl$cell_id, cl$section,
                cl$expected, ok))
  }
  print(table(vapply(.cells, function(z) z$section, character(1L))))
  quit(save = "no")
}

.only <- Sys.getenv("FB_GRID_ONLY")
.suffix <- ""
if (nzchar(.only)) {
  keep <- vapply(.cells, function(z) grepl(.only, z$cell_id), logical(1L))
  .cells <- .cells[keep]
  .suffix <- "_smoke"
  cat("SMOKE MODE: ", length(.cells), " cell(s) matching '", .only, "'\n\n",
      sep = "")
}

.t_start <- Sys.time()
.rows <- parallel::mclapply(.cells, .run_one, mc.cores = .n_cores)
.wall_total <- as.numeric(difftime(Sys.time(), .t_start, units = "mins"))

# --- results assembly -------------------------------------------------- #

.bad <- vapply(.rows, function(r) !is.data.frame(r), logical(1L))
if (any(.bad)) {
  cat("\nWARNING: ", sum(.bad), " cell(s) returned a non-data.frame:\n")
  for (i in which(.bad)) {
    cat("  ", .cells[[i]]$cell_id, ": ",
        paste(utils::capture.output(print(.rows[[i]]))[1L], collapse = ""),
        "\n", sep = "")
  }
  .rows <- .rows[!.bad]
}

results <- do.call(rbind, .rows)

# An error splits into a typed refusal and a raw engine / R error by its
# condition class: every structured refusal carries a `flexybayes_*`
# class (R/refusal_taxonomy.R), so an error without one is untyped
# against the DESCRIPTION contract.
results$class <- ifelse(
  results$observed != "error",
  results$observed,
  ifelse(grepl("flexybayes_", results$cond_class), "refuse_typed",
         "error_untyped")
)

# The translation grid splits a completed fit on arrival. A fit whose
# prior never reached the engine is a silent drop, recorded as such in
# the verdict while the ledger keeps the closed outcome vocabulary.
.t_fit <- results$section == "T" & results$class == "fit"
results$arrived <- results$prior_reached
results$resolved <- results$class
results$resolved[.t_fit] <- ifelse(
  results$prior_reached[.t_fit] %in% TRUE, "translated", "silent_drop"
)

.matches_expectation <- function(expected, resolved) {
  switch(
    expected,
    fit = identical(resolved, "fit"),
    refuse_typed = identical(resolved, "refuse_typed"),
    accept = resolved %in% c("construct_accept", "fit"),
    error = resolved %in% c("refuse_typed", "error_untyped"),
    translated = identical(resolved, "translated"),
    artefact = identical(resolved, "value"),
    NA
  )
}

results$verdict <- ifelse(
  mapply(.matches_expectation, results$expected, results$resolved),
  "as_claimed", "DIVERGENT"
)

# --- artefacts --------------------------------------------------------- #

# (1) The wide record: everything, for the reader chasing one cell.
utils::write.csv(
  results,
  file.path(.raw_dir, paste0("grid_results", .suffix, ".csv")),
  row.names = FALSE
)

# (2) The ledger gate F28 reads. Lean on purpose -- it ships in the
#     package, and the columns the gate names are `cell, expected,
#     observed, class, sha, date`. Engine versions travel with it,
#     since an engine upgrade is the commonest reason a green grid
#     turns red.
.engine_version <- function(pkg) {
  v <- try(as.character(utils::packageVersion(pkg)), silent = TRUE)
  if (inherits(v, "try-error")) NA_character_ else v
}

ledger <- data.frame(
  cell = results$cell_id,
  expected = results$expected,
  observed = results$resolved,
  class = results$class,
  sha = .head_sha,
  date = .run_date,
  section = results$section,
  backend = results$backend,
  family = results$family,
  variant = results$variant,
  matrix_cell = results$matrix_cell,
  arrived = results$arrived,
  verdict = results$verdict,
  wall_sec = results$wall_sec,
  pkg_version = .build$version,
  brms_version = .engine_version("brms"),
  inla_version = .engine_version("INLA"),
  stringsAsFactors = FALSE
)

utils::write.csv(
  ledger, file.path(.out_dir, paste0("ledger", .suffix, ".csv")),
  row.names = FALSE
)

# (3) The machine-readable capability matrix, beside the ledger, so the
#     gate can check ledger coverage of every claimed cell without
#     parsing R.
utils::write.csv(
  .matrix_cells,
  file.path(.out_dir, paste0("capability_matrix", .suffix, ".csv")),
  row.names = FALSE
)

# (4) The misses register (M1-M6). Every open miss from this run, one
#     row each; a clean run leaves the header and no rows.
.miss_class <- function(row) {
  if (identical(row$resolved, "silent_drop")) {
    return("M1 silent-wrong")
  }
  if (identical(row$class, "error_untyped")) {
    return("M2 raw-error")
  }
  if (identical(row$class, "crash") || identical(row$class, "timeout")) {
    return("M2 raw-error")
  }
  if (identical(row$section, "S3a") &&
        identical(row$variant, "allowlist_probe")) {
    return("M4 allowlist-narrower")
  }
  "M3 claims-breach"
}

.open <- results[results$verdict == "DIVERGENT" |
                   results$class == "error_untyped" |
                   results$resolved == "silent_drop", , drop = FALSE]

misses <- if (nrow(.open) == 0L) {
  data.frame(miss_id = character(0), class = character(0),
             cell = character(0), expected = character(0),
             observed = character(0), detail = character(0),
             sha = character(0), date = character(0),
             stringsAsFactors = FALSE)
} else {
  data.frame(
    miss_id = sprintf("EG-%s-%03d", .run_date, seq_len(nrow(.open))),
    class = vapply(seq_len(nrow(.open)),
                   function(i) .miss_class(.open[i, ]), character(1L)),
    cell = .open$cell_id,
    expected = .open$expected,
    observed = .open$resolved,
    detail = substr(paste0(.open$message, " | ", .open$expect_src),
                    1L, 400L),
    sha = .head_sha, date = .run_date,
    stringsAsFactors = FALSE
  )
}

utils::write.csv(
  misses, file.path(.out_dir, paste0("misses", .suffix, ".csv")),
  row.names = FALSE
)

# (5) The allowlist-versus-engine roster diff (G5). One row per family
#     an installed engine carries and the entry allowlist refuses, with
#     the layer that owns the boundary.
.roster_probe <- results[results$section == "S3a" &
                           results$variant == "allowlist_probe", ,
                         drop = FALSE]
.engine_txt <- if (file.exists(.file_engine_families)) {
  readLines(.file_engine_families, warn = FALSE)
} else {
  character(0)
}
.roster_of <- function(marker) {
  idx <- grep(marker, .engine_txt)
  if (length(idx) == 0L) {
    return(character(0))
  }
  trimws(strsplit(.engine_txt[idx[1L] + 1L], ",", fixed = TRUE)[[1L]])
}
.brms_roster <- .roster_of("^# brms family roster")
.inla_roster <- .roster_of("^# INLA likelihood roster")

roster_diff <- data.frame(
  family = .roster_probe$family,
  refused_by_flexybayes = .roster_probe$class == "refuse_typed",
  brms_native = .roster_probe$family %in% .brms_roster,
  inla_native = .roster_probe$family %in% .inla_roster,
  stringsAsFactors = FALSE
)
roster_diff$boundary_owner <- ifelse(
  roster_diff$brms_native | roster_diff$inla_native,
  "flexyBayes (an engine carries it; the package has no emit)",
  "engine (no active engine carries it)"
)
roster_diff$boundary_doc <- ifelse(
  roster_diff$brms_native | roster_diff$inla_native,
  "inst/KNOWN_ISSUES.md, Response families: where the boundary is",
  "inst/KNOWN_ISSUES.md, Response families: where the boundary is"
)
roster_diff <- roster_diff[
  roster_diff$brms_native | roster_diff$inla_native, , drop = FALSE
]

utils::write.csv(
  roster_diff, file.path(.out_dir, paste0("roster_diff", .suffix, ".csv")),
  row.names = FALSE
)

# (6) The session record, captured in a child that has ATTACHED the
#     packages under test, so it names the engines the cells ran
#     against rather than the dispatcher's bare session.
writeLines(
  callr::r(function(lib) {
    .libPaths(c(lib, .libPaths()))
    suppressMessages({
      library(flexyBayes)
      library(brms)
      library(INLA)
    })
    out <- c(
      paste0("flexyBayes Packaged: ",
             sub(";.*$", "",
                 utils::packageDescription("flexyBayes")$Packaged)),
      paste0("flexyBayes version : ",
             as.character(utils::packageVersion("flexyBayes"))),
      "",
      utils::capture.output(utils::sessionInfo())
    )
    # The artefact ships, so it records WHICH build was under test and
    # never WHO ran it: the home directory collapses to `~` and the login
    # name to a placeholder, and the library path is not written at all
    # (a scratch install path is outside the home directory, so the `~`
    # substitution alone would not have caught the name inside it).
    out <- sub(path.expand("~"), "~", out, fixed = TRUE)
    user <- Sys.info()[["user"]]
    if (is.character(user) && nzchar(user)) {
      out <- gsub(user, "<user>", out, fixed = TRUE)
    }
    out
  }, args = list(lib = .lib), timeout = 180, show = FALSE, spinner = FALSE),
  file.path(.out_dir, paste0("sessioninfo", .suffix, ".txt"))
)

# --- console summary --------------------------------------------------- #

cat("\n===== SUMMARY =====\n")
cat("total cells :", nrow(results), "\n")
cat("wall (min)  :", round(.wall_total, 2), "\n")
cat("sha / date  :", .head_sha, "/", .run_date, "\n\n")
cat("-- outcome class histogram --\n")
print(table(results$resolved))
cat("\n-- outcome by section --\n")
print(table(results$section, results$resolved))
cat("\n-- verdict by section --\n")
print(table(results$section, results$verdict))

cat("\n-- error_untyped cells --\n")
.eu <- results[results$class == "error_untyped", , drop = FALSE]
if (nrow(.eu) > 0L) {
  for (i in seq_len(nrow(.eu))) {
    cat(sprintf("  %-46s %s\n", .eu$cell_id[i],
                substr(.eu$message[i], 1L, 110L)))
  }
} else {
  cat("  none\n")
}

cat("\n-- silent drops (fit completed, prior never arrived) --\n")
.sd <- results[results$resolved == "silent_drop", , drop = FALSE]
if (nrow(.sd) > 0L) {
  for (i in seq_len(nrow(.sd))) cat("  ", .sd$cell_id[i], "\n")
} else {
  cat("  none\n")
}

cat("\n-- DIVERGENT cells --\n")
.div <- results[results$verdict == "DIVERGENT", , drop = FALSE]
if (nrow(.div) > 0L) {
  for (i in seq_len(nrow(.div))) {
    cat(sprintf("  %-46s expected %-13s observed %-14s :: %s\n",
                .div$cell_id[i], .div$expected[i], .div$resolved[i],
                substr(.div$message[i], 1L, 110L)))
  }
} else {
  cat("  none\n")
}

cat("\n-- capability-matrix cells --\n")
.m <- results[results$section == "M", , drop = FALSE]
cat("  ", nrow(.m), "of", nrow(.matrix_cells),
    "matrix cells executed;", sum(.m$verdict == "as_claimed"),
    "as claimed\n")

cat("\n-- S1 prior arrival by route --\n")
.pr <- results[results$section == "S1" & results$class == "fit", ,
               drop = FALSE]
print(table(.pr$variant, .pr$prior_reached, useNA = "ifany"))

cat("\nartefacts written to", .out_dir, "\n")
cat("wide record written to", .raw_dir, "\n")

if (nrow(.div) > 0L || nrow(.eu) > 0L || nrow(.sd) > 0L) {
  cat("\nGRID: NOT CLEAN -- see the misses register.\n")
} else {
  cat("\nGRID: clean (0 divergent, 0 error_untyped, 0 silent drops).\n")
}
