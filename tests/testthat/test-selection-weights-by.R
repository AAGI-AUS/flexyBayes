# fb_resamples(weights_by = ): observation weights re-estimated in each fold.

# Sites of 6 to 30 rows with residual variances from 0.2 to 3, and base
# weights median(s2) / s2 by site.
weights_by_fixture <- function() {
  set.seed(31L)
  n_per <- c(6L, 30L, 11L, 19L, 8L, 24L, 14L, 9L)
  s2_site <- c(0.2, 3, 0.5, 1.5, 0.3, 2, 0.8, 1)
  d <- data.frame(site = factor(rep(seq_along(n_per), n_per)))
  n <- nrow(d)
  d$x1 <- stats::rnorm(n)
  d$x2 <- stats::rnorm(n)
  d$x3 <- stats::rnorm(n)
  s2 <- s2_site[as.integer(d$site)]
  d$y <- 1 + 0.6 * d$x1 +
    stats::rnorm(length(n_per), sd = 0.8)[as.integer(d$site)] +
    stats::rnorm(n, sd = sqrt(s2))
  d$w <- stats::median(s2_site) / s2
  d
}

weights_by_base <- function(d, weighted = TRUE) {
  args <- list(
    y ~ 1,
    random = ~ site,
    data = d,
    backend = "inla",
    verbose = FALSE
  )
  if (weighted) {
    args$weights <- d$w
  }
  suppressMessages(do.call(flexybayes, args))
}

weights_by_benchmark <- function(base, covariance = "base_model", ...) {
  fb_select(
    base, ~ x1 + x2 + x3,
    workflow = "unpenalized_benchmark",
    covariance = covariance,
    resamples = fb_resamples(unit = ~ site, v = 4L, weights_by = ~ site),
    seed = 1L,
    ...
  )
}

# One fold's weights per row from its certificate record: the site weight on
# a site with training rows, 1 on held-out rows and on every other row.
weights_by_rows <- function(d, split, record) {
  w <- rep(1, nrow(d))
  site <- as.character(d$site)
  hit <- site %in% names(record$weights_by$weights)
  w[hit] <- record$weights_by$weights[site[hit]]
  w[split$assess_rows] <- 1
  w
}

# The fold's GLS written out with the fold's weights: V = diag(sigma^2 / w) +
# sigma_site^2 Z Z', held-out rows predicted by their fixed part.
weights_by_gls <- function(d, split, record, w, formula) {
  tr <- split$train_rows
  te <- split$assess_rows
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
    mse = vapply(by_site, function(i) mean(residual[i]^2), numeric(1L))
  )
}

test_that("fb_resamples() records weights_by and prints it", {
  r <- fb_resamples(unit = ~ site, v = 4L, weights_by = ~ site)
  expect_identical(r$weights_by, "site")
  expect_null(fb_resamples(unit = ~ site, v = 4L)$weights_by)
  expect_output(print(r), "weights : re-estimated in each fold by site")
})

test_that("each fold's weights are median(s2) / s2 from an unweighted refit on its training rows", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_high_missingness_warning = TRUE)
  d <- weights_by_fixture()
  base <- weights_by_base(d)
  original <- flexybayes
  calls <- list()
  local_mocked_bindings(flexybayes = function(...) {
    fit <- original(...)
    calls[[length(calls) + 1L]] <<- list(args = list(...), fit = fit)
    fit
  })
  sel <- weights_by_benchmark(base)

  unweighted <- which(vapply(calls, function(k) {
    w <- k$args$weights
    is.null(w) || all(w == w[[1L]])
  }, logical(1L)))
  expect_length(unweighted, 4L)
  cov <- sel$certificate$covariance
  labels <- vapply(cov$folds, function(r) r$fold, character(1L))
  # The rule by hand on two folds, from the fitted values of the unweighted
  # refit the package made.
  for (k in unweighted[1:2]) {
    refit <- calls[[k]]
    held <- is.na(refit$args$data$y)
    expect_true(any(held))
    train <- which(!held)
    r2 <- (d$y[train] - stats::fitted(refit$fit)[train])^2
    s2 <- tapply(r2, as.character(d$site[train]), mean)
    w_site <- stats::median(s2) / s2
    expected <- rep(1, nrow(d))
    expected[train] <- w_site[as.character(d$site[train])]
    # The next call is the weighted refit of the same fold.
    expect_lt(max(abs(calls[[k + 1L]]$args$weights - expected)), 1e-12)
    split <- sel$folds[[which(vapply(
      sel$folds,
      function(s) setequal(s$assess_rows, which(held)),
      logical(1L)
    ))]]
    record <- cov$folds[[match(split$label, labels)]]
    expect_lt(max(abs(record$weights_by$s2[names(s2)] - s2)), 1e-12)
    expect_lt(max(abs(split$weights - expected)), 1e-12)
  }
})

test_that("held-out rows of a group with training rows get weight 1", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_high_missingness_warning = TRUE)
  d <- weights_by_fixture()
  base <- flexyBayes:::.fb_sel_bayes_normalise_base(weights_by_base(d))
  assess <- which(d$site %in% c("2", "4"))[c(TRUE, FALSE, FALSE)]
  split <- list(
    label = "1",
    train_rows = setdiff(seq_len(nrow(d)), assess),
    assess_rows = assess
  )
  got <- flexyBayes:::.fb_sel_fold_weights(base, split, "site", "inla")
  expect_setequal(names(got$record$weights), levels(d$site))
  expect_true(all(got$weights[assess] == 1))
  train_site <- as.character(d$site[split$train_rows])
  expect_identical(
    got$weights[split$train_rows],
    unname(got$record$weights[train_site])
  )
})

test_that("each fold is scored as the GLS written out with its own weights", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_high_missingness_warning = TRUE)
  d <- weights_by_fixture()
  base <- weights_by_base(d)
  for (covariance in c("base_model", "reml_plugin")) {
    sel <- weights_by_benchmark(base, covariance = covariance)
    cov <- sel$certificate$covariance
    expect_identical(cov$weights_by, "site")
    clusters <- as.data.frame(sel$clusters)
    labels <- vapply(cov$folds, function(r) r$fold, character(1L))
    formulas <- list(protected_only = ~ 1, full = ~ x1 + x2 + x3)
    for (split in sel$folds) {
      record <- cov$folds[[match(split$label, labels)]]
      w <- weights_by_rows(d, split, record)
      # The fold's weights are not the base's.
      expect_gt(max(abs(w[split$train_rows] / d$w[split$train_rows] - 1)), 0.01)
      for (model in names(formulas)) {
        hand <- weights_by_gls(d, split, record, w, formulas[[model]])
        got <- clusters[clusters$model == model &
                          clusters$fold == split$label, , drop = FALSE]
        got <- got[match(hand$cluster, got$cluster), , drop = FALSE]
        expect_false(anyNA(got$mse))
        expect_lt(max(abs(got$mse / hand$mse - 1)), 1e-6)
      }
    }
  }
  expect_output(print(sel$certificate), "re-estimated in each fold by site")
})

test_that("each fold's s2 agrees with lme4's residuals on the same training rows", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  withr::local_options(flexyBayes.silence_high_missingness_warning = TRUE)
  d <- weights_by_fixture()
  sel <- weights_by_benchmark(weights_by_base(d))
  labels <- vapply(
    sel$certificate$covariance$folds,
    function(r) r$fold,
    character(1L)
  )
  for (split in sel$folds) {
    record <- sel$certificate$covariance$folds[[match(split$label, labels)]]
    train <- d[split$train_rows, , drop = FALSE]
    train$site <- droplevels(train$site)
    lmer_fit <- lme4::lmer(y ~ 1 + (1 | site), data = train, REML = TRUE)
    s2 <- tapply(stats::residuals(lmer_fit)^2, as.character(train$site), mean)
    got <- record$weights_by$s2[names(s2)]
    expect_lt(max(abs(log(got / s2))), 0.1)
  }
})

test_that("the Bayesian refit of each fold is fitted with that fold's weights", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  withr::local_options(flexyBayes.silence_high_missingness_warning = TRUE)
  d <- weights_by_fixture()
  base <- weights_by_base(d)
  original <- flexybayes
  seen <- list()
  local_mocked_bindings(flexybayes = function(...) {
    seen[[length(seen) + 1L]] <<- list(...)$weights
    original(...)
  })
  sel <- fb_select(
    base, ~ x1 + x2 + x3,
    penalty = fb_penalty(
      term(everything()) ~ elastic_net(lambda = 1, alpha = 1)
    ),
    workflow = "penalized_optimizer",
    covariance = "base_model",
    decision = fb_decision(
      method = "full", metric = "rmse", rule = "one_se",
      stability_subsamples = 0L
    ),
    refit = list(method = "bayes"),
    resamples = fb_resamples(unit = ~ site, v = 4L, weights_by = ~ site),
    seed = 1L
  )
  fold_w <- lapply(sel$folds, function(s) s$weights)
  expect_length(fold_w, 4L)
  weighted <- Filter(function(w) !is.null(w) && any(w != w[[1L]]), seen)
  # Per fold: the covariance refit and the Bayesian refit of the support.
  expect_length(weighted, 8L)
  matched <- vapply(weighted, function(w) {
    sum(vapply(fold_w, identical, logical(1L), w))
  }, integer(1L))
  expect_true(all(matched == 1L))
  expect_false(any(vapply(weighted, identical, logical(1L), d$w)))
})

test_that("weights_by refuses where it is not carried", {
  skip_if_not_installed("INLA")
  d <- weights_by_fixture()
  expect_error(
    weights_by_benchmark(weights_by_base(d, weighted = FALSE)),
    class = "flexybayes_refusal_selection_weights_by_unweighted_base"
  )
  base <- weights_by_base(d)
  expect_error(
    fb_select(
      base, ~ x1 + x2 + x3,
      workflow = "bayes",
      decision = fb_decision(method = "full", metric = "rmse"),
      resamples = fb_resamples(unit = ~ site, v = 4L, weights_by = ~ site),
      seed = 1L
    ),
    class = "flexybayes_refusal_selection_weights_not_carried"
  )
  skip_if_not_installed("glmnet")
  expect_error(
    fb_select(
      base, ~ x1 + x2 + x3,
      penalty = fb_penalty(
        term(everything()) ~ elastic_net(lambda = 1, alpha = 1)
      ),
      workflow = "penalized_optimizer",
      decision = fb_decision(
        method = "full", metric = "rmse", rule = "one_se",
        stability_subsamples = 3L
      ),
      resamples = fb_resamples(unit = ~ site, v = 4L, weights_by = ~ site),
      seed = 1L
    ),
    class = "flexybayes_refusal_selection_weights_not_carried"
  )
})
