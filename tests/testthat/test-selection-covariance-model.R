# =============================================================================
# fb_select(covariance = "base_model"): the fold covariance of the benchmark
# and optimiser workflows, read off a refit of the base model on each fold's
# training rows, and whitened block by block where it is block-diagonal.
#
# Each block holds one claim to an oracle the implementation did not write.
#
#   Injection identity. The REML plug-in's own variance components, fed into
#   the base-model route in place of the refit, give the same path, support
#   and scores bit for bit, for nested and for crossed intercepts:
#   everything downstream of the components is the same computation
#   whichever source supplied them.
#
#   Agreement. On one fixture the two sources choose the same support, and
#   their out-of-fold scores differ by less than the standard error of
#   either route's score across scoring units.
#
#   Structures. A random slope, `diag(f):g`, a single-trial
#   `ar1(row):ar1(col)` field and `vm()` / `ped()` precision carriers each run
#   end to end through `fb_select()`. The fold covariance equals V written
#   out element by element from the model's own definition, the recorded
#   components equal a separate refit of the fold, and the conditional
#   covariance of the assessment rows is the Schur complement under
#   `known_cluster` and V_SS under `new_cluster`.
#
#   Blocking. Above the dense bound a block-diagonal V is factored block by
#   block. Each block of V equals the dense V's block exactly, for nested
#   intercepts and for a random slope, and the path agrees with the dense
#   route's to rounding.
#
#   The bound. A covariance that would be formed densely above it refuses by
#   name before any fold runs.
#
# Tolerances, fixed before any result.
#   - Identity and V's blocks: identical(). The same values pass through the
#     same arithmetic in the same order.
#   - V against its element-wise definition, and the Schur complement: 1e-10
#     relative. The two sides sum the same terms in a different order; the
#     rounding of an inner product of at most a hundred terms is below 1e-13.
#   - Dense against block-by-block factorisation: 1.5e-8 (sqrt of machine
#     epsilon) relative. Two Cholesky algorithms on one matrix differ by
#     rounding of order kappa(V) n eps, below 1e-11 on these fixtures.
#   - A recorded component against a separate refit of the same fold: 1e-4
#     relative for a variance, 1e-4 absolute for a correlation. Identical
#     INLA fits on one thread differed by up to 1.1e-5 relative in
#     E[sigma^2] and 1.5e-7 in the posterior mean of a correlation (two
#     probe runs of three fits on each of three fixtures); 1e-4 is ten times
#     the larger.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

fbcm_quiet <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE,
    .local_envir = envir
  )
  if (requireNamespace("INLA", quietly = TRUE)) {
    old <- INLA::inla.getOption("num.threads")
    INLA::inla.setOption(num.threads = "1:1")
    withr::defer(INLA::inla.setOption(num.threads = old), envir = envir)
  }
}

# Twelve sites of eight rows: a protected covariate, five candidates of which
# two carry signal, a site intercept and a site slope on the covariate.
fbcm_site_data <- function(n_site = 12L, per = 8L, seed = 31L) {
  .fb_sel_with_seed(seed, {
    d <- data.frame(
      site = factor(rep(sprintf("s%02d", seq_len(n_site)), each = per)),
      year = factor(rep(c("y1", "y2", "y3"), length.out = n_site * per))
    )
    n <- nrow(d)
    d$x <- stats::rnorm(n)
    for (j in seq_len(5L)) {
      d[[paste0("z", j)]] <- stats::rnorm(n)
    }
    u <- stats::rnorm(n_site, sd = 1.2)
    b <- stats::rnorm(n_site, sd = 0.5)
    d$y <- 2 + 0.6 * d$x + 0.9 * d$z1 - 0.7 * d$z2 +
      u[as.integer(d$site)] + b[as.integer(d$site)] * d$x +
      stats::rnorm(n, sd = 0.6)
    d
  })
}

# A base that carries the parsed model and the data and was never fitted.
fbcm_stub_base <- function(d, fixed = y ~ x, random = ~site) {
  fb <- fb_from_asreml(fixed, random = random, data = d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  structure(
    list(
      fb = fb,
      data = d,
      extras = list(
        call_info = list(fixed = fixed, random = random, family = "gaussian")
      )
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
}

fbcm_penalty <- function() {
  fb_penalty(term(everything()) ~ laplace(scale = 0.15))
}

# The REML plug-in's components for a fold, computed as the REML route
# computes them, standing in for the base-model refit. `calls` counts how
# often the route asked for them.
fbcm_reml_injection <- function(d, protected, random) {
  calls <- new.env()
  calls$labels <- character(0)
  fn <- function(base, split, ...) {
    calls$labels <- c(calls$labels, split$label)
    design <- .fb_sel_design(protected, d, split$train_rows)
    components <- .fb_sel_reml_components(
      as.numeric(d$y),
      design$X,
      .fb_sel_Z_from_vars(d, random),
      split$train_rows,
      label = split$label
    )
    list(
      residual = components$residual,
      terms = components$terms,
      rho = list()
    )
  }
  list(fn = fn, calls = calls)
}

# Fixed components, for comparisons that must not depend on any estimator.
fbcm_fixed_injection <- function(residual, terms) {
  function(base, split, ...) {
    list(residual = residual, terms = terms, rho = list())
  }
}

fbcm_select <- function(base, ..., candidates = ~ z1 + z2 + z3 + z4 + z5,
                        resamples = fb_resamples(unit = ~site, v = 4L)) {
  fb_select(
    base, candidates,
    penalty = fbcm_penalty(),
    workflow = "penalized_optimizer",
    resamples = resamples,
    seed = 1L,
    ...
  )
}

fbcm_rel <- function(a, b) {
  max(abs(a - b), na.rm = TRUE) / max(abs(b), na.rm = TRUE)
}

# E[sigma^2] for every row of a fit's variance-component table.
fbcm_second_moments <- function(fit) {
  vc <- summary(fit)$varcomp
  stats::setNames(vc$estimate^2 + vc$std.error^2, vc$component)
}

# The fold refit, written out: the base's model on the fold's rows with the
# assessment responses missing.
fbcm_hand_refit <- function(d, split, ...) {
  d_fold <- d
  d_fold$y[split$assess_rows] <- NA_real_
  flexybayes(..., data = d_fold, backend = "inla", na_action = "augment",
             verbose = FALSE)
}

fbcm_ar1 <- function(n, rho) {
  rho^abs(outer(seq_len(n), seq_len(n), "-"))
}


# --- (i) injection identity ---------------------------------------------- #

test_that("the plug-in's components fed through base_model give an identical path, support and scores", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("INLA")
  fbcm_quiet()
  d <- fbcm_site_data()
  base <- fbcm_stub_base(d, random = ~ site + site:year)
  reml <- fbcm_select(base)

  injection <- fbcm_reml_injection(d, ~x, c("site", "site:year"))
  local_mocked_bindings(.fb_sel_cov_fold_refit = injection$fn)
  injected <- fbcm_select(base, covariance = "base_model")

  # The route asked for the components once per fold, in fold order.
  expect_identical(injection$calls$labels, c("1", "2", "3", "4"))
  arm_r <- reml$benchmark$glmnet
  arm_i <- injected$benchmark$glmnet
  expect_identical(arm_i$coefficients, arm_r$coefficients)
  expect_identical(arm_i$support, arm_r$support)
  expect_identical(arm_i$selection, arm_r$selection)
  expect_identical(arm_i$clusters_by_model, arm_r$clusters_by_model)
  expect_identical(arm_i$assessment_predictions, arm_r$assessment_predictions)
  expect_identical(injected$selection$chosen_terms, reml$selection$chosen_terms)
  expect_identical(injected$scores, reml$scores)
  expect_identical(injected$clusters, reml$clusters)

  # The records say which source ran; nothing else differs.
  expect_identical(reml$certificate$covariance$source, "reml_plugin")
  expect_identical(injected$certificate$covariance$source, "base_model")
  expect_identical(arm_i$covariance_estimand, "plugin_mean_covariance")
  expect_identical(arm_r$covariance_estimand, "reml_plugin_covariance")

  # The same identity on the unpenalised benchmark alone.
  bench_r <- fb_select(
    base, ~ z1 + z2 + z3 + z4 + z5, workflow = "unpenalized_benchmark",
    resamples = fb_resamples(unit = ~site, v = 4L), seed = 1L
  )
  bench_i <- fb_select(
    base, ~ z1 + z2 + z3 + z4 + z5, workflow = "unpenalized_benchmark",
    resamples = fb_resamples(unit = ~site, v = 4L), seed = 1L,
    covariance = "base_model"
  )
  expect_identical(bench_i$scores$score, bench_r$scores$score)
  expect_identical(bench_i$scores$se, bench_r$scores$se)
  expect_identical(bench_i$clusters, bench_r$clusters)
  expect_identical(bench_i$provider, "gls_base_model")
  expect_identical(bench_r$provider, "gls_reml_plugin")
})


test_that("crossed intercepts: the plug-in's components give an identical path through base_model", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("INLA")
  fbcm_quiet()
  d <- fbcm_site_data()
  base <- fbcm_stub_base(d, random = ~ site + year)
  reml <- fbcm_select(base)
  injection <- fbcm_reml_injection(d, ~x, c("site", "year"))
  local_mocked_bindings(.fb_sel_cov_fold_refit = injection$fn)
  injected <- fbcm_select(base, covariance = "base_model")
  expect_identical(injection$calls$labels, c("1", "2", "3", "4"))
  expect_identical(injected$benchmark$glmnet$coefficients,
                   reml$benchmark$glmnet$coefficients)
  expect_identical(injected$benchmark$glmnet$support,
                   reml$benchmark$glmnet$support)
  expect_identical(injected$scores, reml$scores)
  # No grouping holds both terms: the certificate says the base-model V was
  # factored densely.
  expect_identical(injected$certificate$covariance$builder, "dense_cholesky")
})


# --- (ii) agreement with the plug-in ------------------------------------- #

test_that("the base-model and REML plug-in routes agree on the support and within each route's standard error", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  skip_if_not_long_tier()
  fbcm_quiet()
  d <- fbcm_site_data()
  base <- flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
                     verbose = FALSE)
  reml <- fbcm_select(base)
  model <- fbcm_select(base, covariance = "base_model")
  expect_identical(model$certificate$covariance$source, "base_model")

  expect_setequal(model$selection$chosen_terms, reml$selection$chosen_terms)
  expect_true("z1" %in% model$selection$chosen_terms)
  for (m in c("protected_only", "full")) {
    r <- reml$scores[reml$scores$model == m, ]
    b <- model$scores[model$scores$model == m, ]
    expect_lt(abs(b$score - r$score), min(b$se, r$se))
  }
  # The spread of the base-model route across identical calls: a second run
  # chooses the same support.
  again <- fbcm_select(base, covariance = "base_model")
  expect_setequal(again$selection$chosen_terms, model$selection$chosen_terms)
})


# --- (iii) structures, end to end ---------------------------------------- #

test_that("a random slope runs end to end, with V as the model defines it", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  skip_if_not_long_tier()
  fbcm_quiet()
  d <- fbcm_site_data()
  base <- flexybayes(y ~ x + (x || site), data = d, backend = "inla",
                     verbose = FALSE)
  sel <- fbcm_select(base, covariance = "base_model")
  cov <- sel$certificate$covariance
  expect_identical(cov$source, "base_model")
  expect_identical(cov$engine, "inla")
  expect_identical(cov$estimand, "plugin_mean_covariance")
  expect_false(is.null(sel$benchmark$glmnet))
  expect_length(cov$folds, 4L)
  # INLA is handed every row, the assessment responses missing.
  expect_identical(cov$folds[[1L]]$rows_passed, nrow(d))

  # Fold 1's components against a separate refit of that fold.
  split <- sel$folds[[1L]]
  record <- cov$folds[[1L]]
  hand <- fbcm_second_moments(fbcm_hand_refit(d, split, y ~ x + (x || site)))
  expect_lt(fbcm_rel(record$residual, hand[["sigma"]]), 1e-4)
  expect_lt(fbcm_rel(record$terms[["site"]], hand[["sd_site"]]), 1e-4)
  expect_lt(fbcm_rel(record$terms[["x_site"]], hand[["sd_x_site"]]), 1e-4)

  # V over the fold's rows, element by element, from the recorded values.
  model <- .fb_sel_cov_model(sel$base, seed = 1L)
  spec <- .fb_sel_cov_fold_spec(model, record)
  blocks <- .fb_sel_fold_blocks(spec, split$train_rows, split$assess_rows,
                                "known_cluster")
  same <- outer(d$site, d$site, "==")
  v <- diag(record$residual, nrow(d)) +
    same * (record$terms[["site"]] + record$terms[["x_site"]] * outer(d$x, d$x))
  tr <- split$train_rows
  expect_lt(fbcm_rel(blocks$v_tt, v[tr, tr]), 1e-10)
})

test_that("diag(f):g runs end to end, with one variance per level of f", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  skip_if_not_long_tier()
  fbcm_quiet()
  d <- .fb_sel_with_seed(32L, {
    out <- expand.grid(
      gen = factor(sprintf("g%02d", 1:12)), env = factor(c("E1", "E2", "E3")),
      rep = 1:2
    )
    for (j in 1:4) out[[paste0("z", j)]] <- stats::rnorm(nrow(out))
    effect <- stats::rnorm(36L) * rep(c(0.4, 1, 1.6), each = 12L)
    out$y <- 1 + as.numeric(out$env) + 0.8 * out$z1 +
      effect[(as.integer(out$env) - 1L) * 12L + as.integer(out$gen)] +
      stats::rnorm(nrow(out), sd = 0.4)
    out
  })
  base <- flexybayes(y ~ env, random = ~ diag(env):gen, data = d,
                     backend = "inla", verbose = FALSE)
  sel <- fbcm_select(base, covariance = "base_model",
                     candidates = ~ z1 + z2 + z3 + z4,
                     resamples = fb_resamples(v = 3L),
                     prediction = "known_cluster")
  cov <- sel$certificate$covariance
  expect_identical(cov$source, "base_model")
  levels_env <- c("E1", "E2", "E3")
  slots <- paste0("gen__env", levels_env)
  expect_setequal(names(cov$folds[[1L]]$terms), slots)

  split <- sel$folds[[1L]]
  record <- cov$folds[[1L]]
  hand <- fbcm_second_moments(fbcm_hand_refit(d, split, y ~ env,
                                              random = ~ diag(env):gen))
  for (s in slots) {
    expect_lt(fbcm_rel(record$terms[[s]], hand[[paste0("sd_", s)]]), 1e-4)
  }

  model <- .fb_sel_cov_model(sel$base, seed = 1L)
  spec <- .fb_sel_cov_fold_spec(model, record)
  blocks <- .fb_sel_fold_blocks(spec, split$train_rows, split$assess_rows,
                                "known_cluster")
  variance <- unlist(record$terms[paste0("gen__env", as.character(d$env))])
  same <- outer(d$gen, d$gen, "==") & outer(d$env, d$env, "==")
  v <- diag(record$residual, nrow(d)) + same * variance
  tr <- split$train_rows
  expect_lt(fbcm_rel(blocks$v_tt, v[tr, tr]), 1e-10)
})

test_that("a single-trial ar1(row):ar1(col) field runs end to end, and V_S|T is the Schur complement", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  skip_if_not_long_tier()
  fbcm_quiet()
  n_row <- 8L
  n_col <- 6L
  d <- .fb_sel_with_seed(33L, {
    grid <- expand.grid(col = seq_len(n_col), row = seq_len(n_row))
    out <- data.frame(row = factor(grid$row), col = factor(grid$col))
    for (j in 1:4) out[[paste0("z", j)]] <- stats::rnorm(nrow(out))
    field <- t(chol(kronecker(fbcm_ar1(n_row, 0.7), fbcm_ar1(n_col, 0.4))))
    out$y <- 5 + 0.8 * out$z1 + as.numeric(field %*% stats::rnorm(nrow(out))) +
      stats::rnorm(nrow(out), sd = 0.4)
    out
  })
  base <- flexybayes(y ~ 1, random = ~ ar1(row):ar1(col), data = d,
                     backend = "inla", verbose = FALSE)
  for (prediction in c("known_cluster", "new_cluster")) {
    # The field on this fixture runs to the boundary, and the fit says so.
    # Known issue 9; the expectation holds the package to saying it.
    expect_warning(
      sel <- fbcm_select(base, covariance = "base_model",
                         candidates = ~ z1 + z2 + z3 + z4,
                         resamples = fb_resamples(v = 3L),
                         prediction = prediction),
      "the autoregressive field in this fit did not identify"
    )
    expect_identical(sel$certificate$covariance$source, "base_model")
    expect_identical(sel$certificate$resampling$prediction, prediction)
  }

  split <- sel$folds[[1L]]
  record <- sel$certificate$covariance$folds[[1L]]
  hand_fit <- fbcm_hand_refit(d, split, y ~ 1, random = ~ ar1(row):ar1(col))
  hand <- fbcm_second_moments(hand_fit)
  vc_hand <- summary(hand_fit)$varcomp
  expect_lt(fbcm_rel(record$residual, hand[["sigma"]]), 1e-4)
  expect_lt(fbcm_rel(record$terms[["spatial"]], hand[["sd_spatial"]]), 1e-4)
  for (r in c("row", "col")) {
    expect_lt(
      abs(record$rho$spatial[[r]] -
            vc_hand$estimate[vc_hand$component == paste0("rho_", r)]),
      1e-4
    )
  }

  rho <- record$rho$spatial
  v <- diag(record$residual, nrow(d)) + record$terms[["spatial"]] *
    outer(as.integer(d$row), as.integer(d$row),
          function(a, b) rho[["row"]]^abs(a - b)) *
    outer(as.integer(d$col), as.integer(d$col),
          function(a, b) rho[["col"]]^abs(a - b))
  model <- .fb_sel_cov_model(sel$base, seed = 1L)
  spec <- .fb_sel_cov_fold_spec(model, record)
  tr <- split$train_rows
  as_ <- split$assess_rows
  known <- .fb_sel_fold_blocks(spec, tr, as_, "known_cluster")
  expect_lt(fbcm_rel(known$v_tt, v[tr, tr]), 1e-10)
  schur <- v[as_, as_] - v[as_, tr] %*% solve(v[tr, tr], v[tr, as_])
  expect_lt(fbcm_rel(known$v_s_given_t, schur), 1e-10)
  new <- .fb_sel_fold_blocks(spec, tr, as_, "new_cluster")
  expect_lt(fbcm_rel(new$v_s_given_t, v[as_, as_]), 1e-10)
})

test_that("vm() and ped() precision carriers run end to end, with V = s^2 Z K Z'", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  skip_if_not_long_tier()
  fbcm_quiet()
  levs <- sprintf("g%02d", 1:16)
  fixture <- .fb_sel_with_seed(34L, {
    m <- matrix(stats::rnorm(16L * 24L), 16L)
    k <- tcrossprod(m) / 24 + diag(0.2, 16L)
    dimnames(k) <- list(levs, levs)
    out <- data.frame(gen = factor(rep(levs, each = 4L)))
    for (j in 1:4) out[[paste0("z", j)]] <- stats::rnorm(nrow(out))
    u <- as.numeric(t(chol(k)) %*% stats::rnorm(16L))
    out$y <- 3 + 0.8 * out$z1 + u[as.integer(out$gen)] +
      stats::rnorm(nrow(out), sd = 0.5)
    list(d = out, k = k)
  })
  d <- fixture$d
  q <- solve(fixture$k)
  carriers <- list(
    vm = ~ vm(gen, cov = fb_cov(q, type = "precision")),
    ped = ~ ped(gen, q, use_sparse_precision = TRUE)
  )
  for (nm in names(carriers)) {
    base <- flexybayes(y ~ 1, random = carriers[[nm]], data = d,
                       known_matrices = list(q = q), backend = "inla",
                       verbose = FALSE)
    sel <- fbcm_select(base, covariance = "base_model",
                       candidates = ~ z1 + z2 + z3 + z4,
                       resamples = fb_resamples(v = 3L),
                       prediction = "known_cluster")
    expect_identical(sel$certificate$covariance$source, "base_model")

    split <- sel$folds[[1L]]
    record <- sel$certificate$covariance$folds[[1L]]
    hand <- fbcm_second_moments(fbcm_hand_refit(
      d, split, y ~ 1, random = carriers[[nm]], known_matrices = list(q = q)
    ))
    expect_lt(fbcm_rel(record$terms[["gen"]], hand[["sd_gen_id"]]), 1e-4)

    v <- diag(record$residual, nrow(d)) +
      record$terms[["gen"]] * fixture$k[as.character(d$gen), as.character(d$gen)]
    dimnames(v) <- NULL
    model <- .fb_sel_cov_model(sel$base, seed = 1L)
    spec <- .fb_sel_cov_fold_spec(model, record)
    tr <- split$train_rows
    blocks <- .fb_sel_fold_blocks(spec, tr, split$assess_rows, "known_cluster")
    expect_lt(fbcm_rel(blocks$v_tt, v[tr, tr]), 1e-10)
  }
})


# --- (iv) block by block ----------------------------------------------------- #

test_that("above the dense bound a block-diagonal V is factored block by block, exactly V's blocks", {
  skip_if_not_installed("INLA")
  d <- fbcm_site_data()
  slope <- fbcm_stub_base(d)
  slope$fb <- suppressWarnings(fb_from_brms(y ~ x + (x || site), data = d))
  slope$fb$family <- "gaussian"
  slope$fb$link <- "identity"
  cases <- list(
    nested = list(
      base = fbcm_stub_base(d, random = ~ site + site:year),
      terms = c(site = 1.3, `site:year` = 0.3)
    ),
    slope = list(base = slope, terms = c(site = 1.3, x_site = 0.2))
  )
  folds <- .fb_sel_folds(fb_resamples(v = 3L), d)
  split <- folds[[1L]]
  tr <- split$train_rows
  as_ <- split$assess_rows
  for (case in cases) {
    base <- .fb_sel_bayes_normalise_base(case$base)
    components <- list(residual = 0.4, terms = case$terms, rho = list())
    withr::local_options(flexyBayes.selection_dense_gb = 0.5)
    dense_model <- .fb_sel_cov_model(base, seed = 1L)
    expect_identical(dense_model$whitening$mode, "dense")
    withr::local_options(flexyBayes.selection_dense_gb = 1e-6)
    block_model <- .fb_sel_cov_model(base, seed = 1L)
    expect_identical(block_model$whitening$mode, "by_block")
    expect_identical(block_model$whitening$group, "site")
    expect_identical(block_model$whitening$n_blocks, 12L)

    dense_spec <- .fb_sel_cov_fold_spec(dense_model, components)
    block_spec <- .fb_sel_cov_fold_spec(block_model, components)
    for (prediction in c("known_cluster", "new_cluster")) {
      dense <- .fb_sel_fold_blocks(dense_spec, tr, as_, prediction)
      by_block <- .fb_sel_fold_blocks_by_group(block_spec, tr, as_, prediction)
      # Each block of V is the dense V's block, value for value.
      for (b in by_block$factor_tt$blocks) {
        rows <- tr[b$positions]
        expect_identical(
          .fb_sel_block_covariance(block_spec, rows, rows),
          dense$v_tt[b$positions, b$positions, drop = FALSE]
        )
      }
      expect_lt(abs(by_block$logdet_v_tt - dense$logdet_v_tt) /
                  abs(dense$logdet_v_tt), 1.5e-8)
      y <- matrix(d$y[tr], ncol = 1L)
      expect_lt(
        fbcm_rel(.fb_sel_whiten_train(by_block$factor_tt, y),
                 forwardsolve(t(dense$factor_tt), y)),
        1.5e-8
      )
      key <- as.character(d$site[as_])
      for (k in unique(key)) {
        index <- which(key == k)
        expect_lt(
          fbcm_rel(.fb_sel_conditional_rows(by_block$v_s_given_t, index),
                   dense$v_s_given_t[index, index]),
          1.5e-8
        )
      }
    }
  }
})

test_that("a base without random terms is factored row by row above the bound", {
  skip_if_not_installed("INLA")
  d <- fbcm_site_data()
  base <- .fb_sel_bayes_normalise_base(fbcm_stub_base(d, random = NULL))
  components <- list(residual = 0.4, terms = numeric(0), rho = list())
  withr::local_options(flexyBayes.selection_dense_gb = 1e-6)
  model <- .fb_sel_cov_model(base, seed = 1L)
  expect_identical(model$whitening$mode, "by_block")
  expect_identical(model$whitening$n_blocks, nrow(d))
  split <- .fb_sel_folds(fb_resamples(v = 3L), d)[[1L]]
  blocks <- .fb_sel_fold_blocks_by_group(
    .fb_sel_cov_fold_spec(model, components),
    split$train_rows, split$assess_rows, "known_cluster"
  )
  y <- matrix(d$y[split$train_rows], ncol = 1L)
  expect_lt(
    fbcm_rel(.fb_sel_whiten_train(blocks$factor_tt, y), y / sqrt(0.4)),
    1.5e-8
  )
  expect_identical(
    .fb_sel_conditional_variance(blocks$v_s_given_t),
    rep(0.4, length(split$assess_rows))
  )
})

test_that("the block-by-block route agrees with the dense route through fb_select()", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("INLA")
  fbcm_quiet()
  d <- fbcm_site_data()
  base <- fbcm_stub_base(d)
  local_mocked_bindings(
    .fb_sel_cov_fold_refit = fbcm_fixed_injection(0.4, c(site = 1.3))
  )
  dense <- fbcm_select(base, covariance = "base_model",
                       prediction = "known_cluster",
                       resamples = fb_resamples(v = 4L))
  withr::local_options(flexyBayes.selection_dense_gb = 1e-6)
  by_block <- fbcm_select(base, covariance = "base_model",
                          prediction = "known_cluster",
                          resamples = fb_resamples(v = 4L))
  expect_identical(dense$certificate$covariance$whitening$mode, "dense")
  expect_identical(by_block$certificate$covariance$whitening$mode, "by_block")
  expect_identical(by_block$selection$chosen_terms, dense$selection$chosen_terms)
  # Membership at every point below the null anchor. At s = 1 the first
  # column to enter sits exactly on its entry threshold, so whether it is
  # reported in the support there is decided by rounding on either route
  # (probe_05_block_dense_history.log: the only two differing rows of 400).
  hb <- by_block$benchmark$glmnet$support$history
  hd <- dense$benchmark$glmnet$support$history
  expect_identical(hb[, c("fold", "lambda_index", "column")],
                   hd[, c("fold", "lambda_index", "column")])
  below <- hd$s < 1
  expect_gt(sum(below), 0L)
  expect_identical(hb$in_support[below], hd$in_support[below])
  expect_lt(fbcm_rel(by_block$scores$score, dense$scores$score), 1.5e-8)
  cb <- by_block$benchmark$glmnet$coefficients
  cd <- dense$benchmark$glmnet$coefficients
  expect_identical(cb$term, cd$term)
  expect_identical(is.na(cb$estimate), is.na(cd$estimate))
  expect_lt(fbcm_rel(cb$estimate, cd$estimate), 1.5e-8)
})


# --- (v) the dense bound ------------------------------------------------- #

test_that("a dense covariance above the bound refuses by name before any fold runs", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("INLA")
  fbcm_quiet()
  d <- fbcm_site_data()
  withr::local_options(flexyBayes.selection_dense_gb = 1e-6)
  # The REML plug-in estimates on the dense matrix itself.
  local_mocked_bindings(.fb_sel_reml_components = function(...) {
    stop("the REML plug-in ran")
  })
  expect_error(
    fbcm_select(fbcm_stub_base(d)),
    class = "flexybayes_refusal_selection_dense_covariance_too_large"
  )
  # Crossed intercepts: no grouping holds every term, so V has no blocks.
  local_mocked_bindings(.fb_sel_cov_fold_refit = function(...) {
    stop("the base model was refitted")
  })
  err <- expect_error(
    fbcm_select(fbcm_stub_base(d, random = ~ site + year),
                covariance = "base_model"),
    class = "flexybayes_refusal_selection_dense_covariance_too_large"
  )
  expect_match(conditionMessage(err), "flexyBayes.selection_dense_gb",
               fixed = TRUE)
  expect_match(conditionMessage(err), "not block-diagonal", fixed = TRUE)
})


# --- (vi) what the base-model route does not build ------------------------ #

test_that("a structure the base-model route does not build refuses by name", {
  skip_if_not_installed("INLA")
  d <- fbcm_site_data()
  d$block <- factor(rep(c("b1", "b2"), length.out = nrow(d)))
  base <- fbcm_stub_base(d, random = ~ site:year:block)
  expect_error(
    fbcm_select(base, covariance = "base_model"),
    class = "flexybayes_refusal_selection_covariance_not_buildable"
  )
})


# --- (vii) the argument ---------------------------------------------------- #

test_that("covariance defaults to the REML plug-in, and a legacy estimand string warns and runs it", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("INLA")
  fbcm_quiet()
  withr::local_options(flexyBayes.deprecation_warnings = TRUE)
  expect_identical(eval(formals(fb_select)$covariance),
                   c("reml_plugin", "base_model"))
  d <- fbcm_site_data()
  expect_warning(base <- fbcm_stub_base(d), class = "flexybayes_deprecated")
  default <- fbcm_select(base)
  expect_warning(
    legacy <- fbcm_select(base, covariance = "plugin_mean_covariance"),
    class = "flexybayes_deprecated"
  )
  expect_identical(legacy$scores, default$scores)
  expect_identical(legacy$certificate$covariance$source, "reml_plugin")
})

test_that("covariance = \"base_model\" on a route that plugs in no covariance is an error", {
  skip_if_not_installed("INLA")
  d <- fbcm_site_data()
  expect_error(
    fb_select(fbcm_stub_base(d), ~ z1 + z2, covariance = "base_model",
              seed = 1L),
    "integrates the covariance"
  )
})


# --- (viii) the record -------------------------------------------------- #

test_that("the certificate records the source, engine, point estimate, blocking and equivalence level", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("INLA")
  fbcm_quiet()
  d <- fbcm_site_data()
  base <- fbcm_stub_base(d)
  local_mocked_bindings(
    .fb_sel_cov_fold_refit = fbcm_fixed_injection(0.4, c(site = 1.3))
  )
  sel <- fbcm_select(base, covariance = "base_model")
  cov <- sel$certificate$covariance
  expect_identical(cov$source, "base_model")
  expect_identical(cov$engine, "inla")
  expect_identical(cov$engine_version,
                   as.character(utils::packageVersion("INLA")))
  expect_match(cov$point_estimate, "posterior mean of each variance",
               fixed = TRUE)
  expect_identical(cov$whitening$mode, "dense")
  expect_identical(cov$equivalence_level, "map_exact_fixed_covariance")
  expect_identical(cov$folds[[2L]]$terms, c(site = 1.3))
  expect_identical(sel$certificate$benchmark$gls$covariance_source,
                   "base_model")
  expect_output(print(sel$certificate), "refit of the base model")

  # The plan counts the refits.
  plan <- fbcm_select(base, covariance = "base_model", plan = TRUE)
  expect_identical(plan$selection$covariance, "base_model")
  # One refit per fold, on four folds, beside the base fit itself.
  expect_identical(plan$selection$cost$n_fits, 5L)
})


# --- (ix) a brms base ----------------------------------------------------- #

test_that("a brms base is refitted on brms, the training rows only", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("glmnet")
  skip_if_not_long_tier()
  fbcm_quiet()
  d <- fbcm_site_data()
  base <- suppressWarnings(flexybayes(
    y ~ x, random = ~site, data = d, backend = "brms", chains = 1,
    n_samples = 400, warmup = 400, seed = 1L, verbose = FALSE,
    mcmc_verbose = FALSE
  ))
  sel <- suppressWarnings(fbcm_select(
    base, covariance = "base_model",
    engine = list(chains = 1, n_samples = 400, warmup = 400)
  ))
  cov <- sel$certificate$covariance
  expect_identical(cov$engine, "brms")
  expect_length(cov$folds, 4L)
  expect_identical(
    vapply(cov$folds, function(f) f$rows_passed, integer(1L)),
    vapply(sel$folds, function(s) length(s$train_rows), integer(1L))
  )
})
