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

test_that("a weighted base refuses on every route that does not carry weights", {
  skip_if_not_installed("INLA")
  d <- weights_fixture()
  base <- weights_base(d)
  penalty <- fb_penalty(term(everything()) ~ elastic_net(lambda = 1, alpha = 1))
  decision <- fb_decision(
    method = "full", metric = "rmse", rule = "one_se",
    stability_subsamples = 0L
  )
  routes <- list(
    list(workflow = "penalized_optimizer", covariance = "base_model",
         penalty = penalty, decision = decision),
    list(workflow = "penalized_optimizer", covariance = "reml_plugin",
         penalty = penalty, decision = decision,
         refit = list(method = "bayes")),
    list(workflow = "bayes", covariance = "reml_plugin")
  )
  for (r in routes) {
    expect_error(
      do.call(weights_select, c(list(base = base, d = d), r)),
      class = "flexybayes_refusal_selection_weights_not_carried"
    )
  }
})

test_that("the per-fold Bayesian refit refuses weights by name", {
  base <- list(extras = list(call_info = list(weights = c(1, 2, 1, 2))))
  expect_error(
    flexyBayes:::.fb_sel_bayes_refit(
      base, data.frame(y = 1:4), y ~ 1, "inla"
    ),
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
