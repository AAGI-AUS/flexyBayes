# A base fitted with observation weights: fb_select() uses them or refuses.

weights_fixture <- function() {
  set.seed(11L)
  n_site <- 12L
  n_per <- 20L
  d <- data.frame(site = factor(rep(seq_len(n_site), each = n_per)))
  d$x1 <- stats::rnorm(nrow(d))
  d$x2 <- stats::rnorm(nrow(d))
  d$x3 <- stats::rnorm(nrow(d))
  # Residual variances 0.2 and 3 on alternate sites: weights 15-fold apart.
  s2 <- rep(c(0.2, 3), length.out = n_site)[as.integer(d$site)]
  d$y <- 1 + 0.6 * d$x1 + stats::rnorm(n_site)[as.integer(d$site)] +
    stats::rnorm(nrow(d), sd = sqrt(s2))
  d$w <- stats::median(s2) / s2
  d
}

weights_base <- function(d) {
  suppressMessages(flexybayes(
    y ~ 1,
    random = ~ site,
    data = d,
    weights = d$w,
    backend = "inla",
    verbose = FALSE
  ))
}

weights_select <- function(base, d, ...) {
  fb_select(
    base = base,
    candidates = ~ x1 + x2 + x3,
    resamples = fb_resamples(unit = ~ site, v = 6L),
    seed = 1L,
    ...
  )
}

# Sites of 6 to 30 rows, and weights that vary within a site as well as
# between sites, so no weighting error cancels by balance.
weights_unbalanced <- function() {
  set.seed(23L)
  n_per <- rep(c(6L, 30L, 11L, 19L, 8L, 24L), 2L)
  d <- data.frame(site = factor(rep(seq_along(n_per), n_per)))
  n <- nrow(d)
  d$x1 <- stats::rnorm(n)
  d$x2 <- stats::rnorm(n)
  d$x3 <- stats::rnorm(n)
  d$w <- stats::runif(n, 0.1, 4)
  d$y <- 1 + 0.6 * d$x1 +
    stats::rnorm(length(n_per), sd = 0.8)[as.integer(d$site)] +
    stats::rnorm(n, sd = sqrt(0.5 / d$w))
  d
}

# The fold's GLS written out: V = diag(sigma^2 / w) + sigma_site^2 Z Z', with
# held-out rows at w = 1, and the held-out rows predicted by their fixed
# part (new_cluster).
weights_gls_by_hand <- function(d, split, record, formula) {
  tr <- split$train_rows
  te <- split$assess_rows
  w <- d$w
  w[te] <- 1
  z <- stats::model.matrix(~ site - 1, d)
  v <- diag(record$residual / w) + record$terms[["site"]] * tcrossprod(z)
  x <- stats::model.matrix(formula, d)
  v_inv <- solve(v[tr, tr])
  beta <- solve(
    t(x[tr, ]) %*% v_inv %*% x[tr, ],
    t(x[tr, ]) %*% v_inv %*% d$y[tr]
  )
  residual <- d$y[te] - drop(x[te, ] %*% beta)
  site <- as.character(d$site[te])
  by_site <- split(seq_along(te), factor(site, levels = unique(site)))
  data.frame(
    cluster = names(by_site),
    mse = vapply(by_site, function(i) mean(residual[i]^2), numeric(1L)),
    nlpd = vapply(by_site, function(i) {
      v_s <- v[te[i], te[i], drop = FALSE]
      0.5 * (length(i) * log(2 * pi) + c(determinant(v_s)$modulus) +
        sum(residual[i] * solve(v_s, residual[i]))) / length(i)
    }, numeric(1L))
  )
}

test_that("a weighted base refuses the projection decision by name", {
  skip_if_not_installed("INLA")
  d <- weights_fixture()
  base <- weights_base(d)
  expect_error(
    weights_select(base, d, workflow = "bayes"),
    class = "flexybayes_refusal_selection_weights_not_carried"
  )
})

# A weighted base that carries its parsed model and data and was never fitted.
weights_stub_base <- function(d, class) {
  fb <- fb_from_asreml(y ~ 1, random = ~ site, data = d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  structure(
    list(
      fb = fb,
      data = d,
      extras = list(call_info = list(
        fixed = y ~ 1, random = ~ site, family = "gaussian", weights = d$w
      ))
    ),
    class = c(class, "flexybayes", "list")
  )
}

test_that("each fold refit is handed the weights of the rows it fits", {
  d <- weights_unbalanced()
  split <- list(label = "1", train_rows = which(!d$site %in% c("2", "9")),
                assess_rows = which(d$site %in% c("2", "9")))
  for (provider in c("inla", "brms")) {
    base <- flexyBayes:::.fb_sel_bayes_normalise_base(
      weights_stub_base(d, paste0("flexybayes_", provider))
    )
    seen <- NULL
    local_mocked_bindings(flexybayes = function(...) {
      seen <<- list(...)
      stop("captured")
    })
    expect_error(
      flexyBayes:::.fb_sel_cov_fold_refit(
        base, split, provider = provider, slots = list()
      ),
      "captured"
    )
    rows <- if (provider == "inla") seq_len(nrow(d)) else split$train_rows
    expect_identical(seen$weights, d$w[rows])
    expect_false(".fb_sel_weight" %in% names(seen$data))
  }
})

test_that("a fold refit of a weighted base refuses data without the weights", {
  d <- weights_unbalanced()
  base <- weights_stub_base(d, "flexybayes_inla")
  expect_error(
    flexyBayes:::.fb_sel_bayes_refit(base, d, y ~ 1, "inla"),
    class = "flexybayes_refusal_selection_weights_not_carried"
  )
})

test_that("the residual diagonal divides by each row's weight", {
  z <- list()
  spec <- flexyBayes:::.fb_sel_cov_spec(
    z, 2, stats::setNames(numeric(0), character(0)),
    weights = c(1, 4, 0.5, 2)
  )
  expect_equal(
    flexyBayes:::.fb_sel_residual_diagonal(spec, c(2L, 3L)),
    c(0.5, 4)
  )
  both <- flexyBayes:::.fb_sel_cov_spec(
    z, 2, stats::setNames(numeric(0), character(0)),
    known_variance = c(0.1, 0.2, 0.3, 0.4),
    weights = c(1, 4, 0.5, 2)
  )
  expect_equal(
    flexyBayes:::.fb_sel_residual_diagonal(both, 1:4),
    c(2.1, 0.7, 4.3, 1.4)
  )
})

test_that("the weighted REML plug-in matches lme4 with the same weights", {
  skip_if_not_installed("lme4")
  d <- weights_unbalanced()
  train <- which(!d$site %in% c("2", "9"))
  x <- stats::model.matrix(~ x1, d)
  z <- flexyBayes:::.fb_sel_Z_from_vars(d, "site")
  got <- flexyBayes:::.fb_sel_reml_components(
    d$y, x, z, train, weights = d$w
  )
  oracle <- lme4::lmer(
    y ~ x1 + (1 | site),
    data = d[train, ],
    weights = w,
    REML = TRUE
  )
  vc <- as.data.frame(lme4::VarCorr(oracle))
  expect_lt(abs(got$residual / stats::sigma(oracle)^2 - 1), 1e-4)
  expect_lt(abs(got$terms[[1L]] / vc$vcov[vc$grp == "site"] - 1), 1e-3)

  unweighted <- flexyBayes:::.fb_sel_reml_components(d$y, x, z, train)
  expect_gt(abs(unweighted$residual / got$residual - 1), 0.05)
})

test_that("fb_select() on a weighted base scores each fold as the weighted GLS written out", {
  skip_if_not_installed("INLA")
  d <- weights_unbalanced()
  base <- weights_base(d)
  sel <- fb_select(
    base, ~ x1 + x2 + x3,
    workflow = "unpenalized_benchmark",
    resamples = fb_resamples(unit = ~ site, v = 4L),
    seed = 1L
  )
  cov <- sel$certificate$covariance
  expect_true(cov$weights)
  expect_output(print(sel$certificate), "divided by each training row's weight")

  clusters <- as.data.frame(sel$clusters)
  labels <- vapply(cov$folds, function(r) r$fold, character(1L))
  formulas <- list(protected_only = ~ 1, full = ~ x1 + x2 + x3)
  for (split in sel$folds) {
    record <- cov$folds[[match(split$label, labels)]]
    for (model in names(formulas)) {
      hand <- weights_gls_by_hand(d, split, record, formulas[[model]])
      got <- clusters[clusters$model == model &
                        clusters$fold == split$label, , drop = FALSE]
      got <- got[match(hand$cluster, got$cluster), , drop = FALSE]
      expect_false(anyNA(got$mse))
      expect_lt(max(abs(got$mse / hand$mse - 1)), 1e-6)
      expect_lt(max(abs(got$nlpd / hand$nlpd - 1)), 1e-6)
    }
  }

  # The same data fitted without weights selects on other components.
  unweighted <- suppressMessages(flexybayes(
    y ~ 1, random = ~ site, data = d, backend = "inla", verbose = FALSE
  ))
  plain <- fb_select(
    unweighted, ~ x1 + x2 + x3,
    workflow = "unpenalized_benchmark",
    resamples = fb_resamples(unit = ~ site, v = 4L),
    seed = 1L
  )
  expect_false(plain$certificate$covariance$weights)
  residual_of <- function(s) {
    vapply(s$certificate$covariance$folds, function(r) r$residual, 1)
  }
  expect_gt(max(abs(residual_of(sel) / residual_of(plain) - 1)), 0.05)
})

test_that("every fold and every stability subsample of the optimiser is weighted", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  d <- weights_unbalanced()
  base <- weights_base(d)
  original <- flexyBayes:::.fb_sel_reml_components
  seen <- logical(0)
  local_mocked_bindings(
    .fb_sel_reml_components = function(..., weights = NULL) {
      seen <<- c(seen, !is.null(weights))
      original(..., weights = weights)
    }
  )
  sel <- fb_select(
    base, ~ x1 + x2 + x3,
    penalty = fb_penalty(
      term(everything()) ~ elastic_net(lambda = 1, alpha = 1)
    ),
    workflow = "penalized_optimizer",
    decision = fb_decision(
      method = "full", metric = "rmse", rule = "one_se",
      stability_subsamples = 3L
    ),
    resamples = fb_resamples(unit = ~ site, v = 4L),
    seed = 1L
  )
  expect_true(sel$certificate$covariance$weights)
  expect_length(seen, 4L + 3L)
  expect_true(all(seen))
})

test_that("covariance = \"base_model\" on a weighted base scores each fold as the weighted GLS written out", {
  skip_if_not_installed("INLA")
  # A held-out locality leaves a third of a fold's rows without a response.
  withr::local_options(flexyBayes.silence_high_missingness_warning = TRUE)
  d <- weights_unbalanced()
  base <- weights_base(d)
  run <- function(b) {
    fb_select(
      b, ~ x1 + x2 + x3,
      workflow = "unpenalized_benchmark",
      covariance = "base_model",
      resamples = fb_resamples(unit = ~ site, v = 4L),
      seed = 1L
    )
  }
  sel <- run(base)
  cov <- sel$certificate$covariance
  expect_true(cov$weights)
  clusters <- as.data.frame(sel$clusters)
  labels <- vapply(cov$folds, function(r) r$fold, character(1L))
  formulas <- list(protected_only = ~ 1, full = ~ x1 + x2 + x3)
  for (split in sel$folds) {
    record <- cov$folds[[match(split$label, labels)]]
    for (model in names(formulas)) {
      hand <- weights_gls_by_hand(d, split, record, formulas[[model]])
      got <- clusters[clusters$model == model &
                        clusters$fold == split$label, , drop = FALSE]
      got <- got[match(hand$cluster, got$cluster), , drop = FALSE]
      expect_lt(max(abs(got$mse / hand$mse - 1)), 1e-6)
      expect_lt(max(abs(got$nlpd / hand$nlpd - 1)), 1e-6)
    }
  }

  # Unweighted, the same folds read other components off their refits.
  plain <- run(suppressMessages(flexybayes(
    y ~ 1, random = ~ site, data = d, backend = "inla", verbose = FALSE
  )))
  expect_false(plain$certificate$covariance$weights)
  residual_of <- function(s) {
    vapply(s$certificate$covariance$folds, function(r) r$residual, 1)
  }
  expect_gt(max(abs(residual_of(sel) / residual_of(plain) - 1)), 0.05)
})

test_that("every Bayesian fold refit of a weighted base is fitted with its weights", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_high_missingness_warning = TRUE)
  d <- weights_unbalanced()
  base <- weights_base(d)
  original <- flexybayes
  seen <- list()
  local_mocked_bindings(flexybayes = function(...) {
    args <- list(...)
    seen[[length(seen) + 1L]] <<- args$weights
    original(...)
  })
  sel <- fb_select(
    base, ~ x1 + x2 + x3,
    workflow = "bayes",
    decision = fb_decision(method = "full", metric = "rmse"),
    resamples = fb_resamples(unit = ~ site, v = 4L),
    seed = 1L
  )
  expect_s3_class(sel, "fb_select")
  # The timed cost probe, then the protected baseline and the full model on
  # each of 4 folds.
  expect_length(seen, 1L + 8L)
  expect_true(all(vapply(seen, identical, logical(1L), d$w)))

  seen <- list()
  opt <- fb_select(
    base, ~ x1 + x2 + x3,
    penalty = fb_penalty(
      term(everything()) ~ elastic_net(lambda = 1, alpha = 1)
    ),
    workflow = "penalized_optimizer",
    decision = fb_decision(
      method = "full", metric = "rmse", rule = "one_se",
      stability_subsamples = 0L
    ),
    refit = list(method = "bayes"),
    resamples = fb_resamples(unit = ~ site, v = 4L),
    seed = 1L
  )
  expect_gt(length(seen), 0L)
  expect_true(all(vapply(seen, identical, logical(1L), d$w)))
})
