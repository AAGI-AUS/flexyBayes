# =============================================================================
# The penalty path grid set on fb_decision(n_lambda, lambda_min_ratio).
#
# The grid is built in one place, .fb_sel_glmnet_grid(). The tests mock it
# with a wrapper that records the settings it receives and then builds the
# real grid, so a value that fb_decision() stores but no route passes on
# fails here on the route that dropped it.
# =============================================================================

fbpg_data <- function(n_cluster = 8L, per = 6L, seed = 51L) {
  .fb_sel_with_seed(seed, {
    n <- n_cluster * per
    site <- factor(rep(seq_len(n_cluster), each = per))
    u <- stats::rnorm(n_cluster, sd = 1.2)
    d <- data.frame(
      site = site,
      x = stats::rnorm(n),
      f = factor(rep(c("a", "b", "c"), length.out = n)),
      g = factor(rep(c("p", "q", "r"), each = 2L, length.out = n))
    )
    for (j in seq_len(4L)) {
      d[[paste0("z", j)]] <- stats::rnorm(n)
    }
    d$y <- 3 + 0.8 * d$x + 1.3 * d$z1 - 1.0 * d$z2 +
      1.1 * as.integer(d$f == "b") +
      u[as.integer(site)] + stats::rnorm(n, sd = 0.5)
    d
  })
}


# A base object with the parsed model and data the optimiser route reads;
# the optimiser never samples it.
fbpg_stub_base <- function(d) {
  fb <- fb_from_asreml(y ~ x, random = ~site, data = d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  structure(
    list(
      fb = fb,
      data = d,
      extras = list(
        call_info = list(fixed = y ~ x, random = ~site, family = "gaussian")
      )
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
}


fbpg_penalties <- function() {
  list(
    elastic_net = fb_penalty(
      term(everything()) ~ elastic_net(lambda = 0.1, alpha = 0.5)
    ),
    laplace = fb_penalty(term(everything()) ~ laplace(scale = 0.7)),
    ridge = fb_penalty(term(everything()) ~ ridge(scale = 0.7)),
    group_lasso = fb_penalty(
      term(everything()) ~ group_lasso(scale = 0.7)
    ),
    sparse_group_lasso = fb_penalty(
      term(everything()) ~ sparse_group_lasso(scale = 0.7, alpha = 0.35)
    ),
    exclusive_lasso = fb_penalty(
      term(everything()) ~ exclusive_lasso(scale = 0.7)
    ),
    group_exclusive = fb_penalty(
      term(everything()) ~ group_exclusive(
        group_scale = 0.7, exclusive_scale = 0.9
      )
    )
  )
}


fbpg_select <- function(d, penalty, decision, ...) {
  fb_select(
    fbpg_stub_base(d),
    candidates = ~ f + g + z1 + z2 + z3 + z4,
    penalty = penalty,
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~site, v = 2L),
    decision = decision,
    seed = 91L,
    ...
  )
}


fbpg_real_grid <- .fb_sel_glmnet_grid


# Replaces the grid builder with one that records what it was handed.
fbpg_record_grid <- function(env = parent.frame()) {
  seen <- new.env()
  seen$calls <- list()
  testthat::local_mocked_bindings(
    .fb_sel_glmnet_grid = function(parts_by_fold, penalty, decision,
                                   groups_by_fold = NULL) {
      seen$calls[[length(seen$calls) + 1L]] <- list(
        n_lambda = decision$n_lambda,
        lambda_min_ratio = decision$lambda_min_ratio
      )
      fbpg_real_grid(parts_by_fold, penalty, decision, groups_by_fold)
    },
    .env = env
  )
  seen
}


test_that("fb_decision() stores the path grid, with the 0.14.0 values as defaults", {
  d <- fb_decision()
  expect_identical(d$n_lambda, 20L)
  expect_identical(d$lambda_min_ratio, 1e-3)
  d <- fb_decision(n_lambda = 7, lambda_min_ratio = 0.05)
  expect_identical(d$n_lambda, 7L)
  expect_identical(d$lambda_min_ratio, 0.05)
  expect_output(print(d), "path    : 7 points")
})


test_that("fb_decision() refuses an invalid grid by the argument's name", {
  for (bad in list(1L, 0L, -3L, 2.5, NA_real_, Inf, "5", c(5L, 6L), NULL)) {
    expect_error(
      fb_decision(n_lambda = bad),
      regexp = "`n_lambda`",
      class = "flexybayes_refusal_selection_decision_invalid"
    )
  }
  for (bad in list(0, 1, -0.1, 2, NA_real_, Inf, "0.01", c(0.1, 0.01), NULL)) {
    expect_error(
      fb_decision(lambda_min_ratio = bad),
      regexp = "`lambda_min_ratio`",
      class = "flexybayes_refusal_selection_decision_invalid"
    )
  }
})


test_that("a non-default grid reaches the grid builder on every penalised family", {
  skip_if_not_installed("glmnet")
  d <- fbpg_data()
  for (family in names(fbpg_penalties())) {
    seen <- fbpg_record_grid()
    fit <- fbpg_select(
      d, fbpg_penalties()[[family]],
      fb_decision(
        method = "full", metric = "rmse",
        n_lambda = 5L, lambda_min_ratio = 0.05
      )
    )
    expect_gt(length(seen$calls), 0L)
    expect_true(
      all(vapply(seen$calls, function(x) identical(x$n_lambda, 5L),
                 logical(1L))),
      info = family
    )
    expect_true(
      all(vapply(seen$calls, function(x) identical(x$lambda_min_ratio, 0.05),
                 logical(1L))),
      info = family
    )
    grid <- fit$certificate$benchmark$glmnet$path_grid
    expect_identical(grid$n_lambda, 5L, info = family)
    expect_identical(grid$lambda_min_ratio, 0.05, info = family)
    scores <- fit$benchmark$glmnet$selection$scores
    if (identical(grid$tuning_scale, "s = lambda / lambda_max") ||
          identical(grid$tuning_scale, "s = multiplier / multiplier_max")) {
      # A family with a finite lambda_max walks s from 1 down to the ratio.
      expect_equal(grid$s, exp(seq(0, log(0.05), length.out = 5L)),
                   info = family)
      expect_equal(sort(unique(scores$s), decreasing = TRUE), grid$s,
                   info = family)
    } else {
      # Ridge and the exclusive penalty have no finite lambda_max: n_lambda
      # decades either side of the declared point, plus that point.
      expect_true(all(is.na(grid$s)), info = family)
      expect_lte(nrow(scores), 6L)
      expect_gte(nrow(scores), 5L)
    }
  }
})


test_that("the stability subsamples and the ablation walk the declared grid", {
  skip_if_not_installed("glmnet")
  d <- fbpg_data()
  seen <- fbpg_record_grid()
  fit <- fbpg_select(
    d, fbpg_penalties()$laplace,
    fb_decision(
      method = "full", metric = "rmse",
      n_lambda = 4L, lambda_min_ratio = 0.02,
      stability_subsamples = 2L, stability_unit = ~site
    )
  )
  # One cross-validated path and one path per subsample.
  expect_identical(length(seen$calls), 3L)
  expect_true(all(vapply(seen$calls, function(x) {
    identical(x$n_lambda, 4L) && identical(x$lambda_min_ratio, 0.02)
  }, logical(1L))))

  seen <- fbpg_record_grid()
  fit <- fbpg_select(
    d, fbpg_penalties()$group_lasso,
    fb_decision(
      method = "full", metric = "rmse",
      n_lambda = 4L, lambda_min_ratio = 0.02
    ),
    ablation = TRUE
  )
  expect_gt(length(seen$calls), 1L)
  expect_true(all(vapply(seen$calls, function(x) {
    identical(x$n_lambda, 4L) && identical(x$lambda_min_ratio, 0.02)
  }, logical(1L))))
})


test_that("the plan counts and prints the declared grid", {
  d <- fbpg_data()
  plan <- fbpg_select(
    d, fbpg_penalties()$group_lasso,
    fb_decision(
      n_lambda = 6L, lambda_min_ratio = 0.01,
      stability_subsamples = 3L, stability_unit = ~site
    ),
    plan = TRUE
  )
  expect_identical(plan$selection$path_grid$n_lambda, 6L)
  expect_identical(plan$selection$path_grid$lambda_min_ratio, 0.01)
  expect_identical(plan$selection$stability$fits_per_path, 7L)
  expect_output(print(plan), "Path grid:    6 points")
  expect_output(print(plan), "0.01")

  plan <- fbpg_select(
    d, fbpg_penalties()$group_lasso,
    fb_decision(n_lambda = 6L, lambda_min_ratio = 0.01),
    ablation = TRUE,
    plan = TRUE
  )
  expect_identical(plan$selection$ablation$fits_per_path, 7L)
})


test_that("the Bayesian plan states the grid of its benchmark arm", {
  skip_if_not_installed("brms")
  d <- fbpg_data()
  args <- list(
    fbpg_stub_base(d),
    candidates = ~ z1 + z2,
    penalty = fbpg_penalties()$laplace,
    workflow = "bayes",
    resamples = fb_resamples(unit = ~site, v = 2L),
    decision = fb_decision(method = "full", n_lambda = 9L,
                           lambda_min_ratio = 0.02),
    seed = 1L,
    plan = TRUE
  )
  plan <- do.call(fb_select, c(args, benchmark = TRUE))
  expect_identical(plan$selection$path_grid$n_lambda, 9L)
  expect_output(print(plan), "Path grid:    9 points")
  plan <- do.call(fb_select, c(args, benchmark = FALSE))
  expect_null(plan$selection$path_grid)
})


test_that("the fit prints its grid, including where the ratio is not read", {
  skip_if_not_installed("glmnet")
  d <- fbpg_data()
  fit <- fbpg_select(
    d, fbpg_penalties()$laplace,
    fb_decision(
      method = "full", metric = "rmse",
      n_lambda = 5L, lambda_min_ratio = 0.05
    )
  )
  expect_output(print(fit), "path grid: 5 points")
  expect_output(print(fit), "0.05")
  fit <- fbpg_select(
    d, fbpg_penalties()$exclusive_lasso,
    fb_decision(
      method = "full", metric = "rmse",
      n_lambda = 5L, lambda_min_ratio = 0.05
    )
  )
  expect_output(print(fit), "lambda_min_ratio is not read")
})


test_that("the default grid is the one 0.14.0 walked", {
  skip_if_not_installed("glmnet")
  d <- fbpg_data()
  seen <- fbpg_record_grid()
  fit <- fbpg_select(
    d, fbpg_penalties()$elastic_net,
    fb_decision(method = "full", metric = "rmse")
  )
  expect_true(all(vapply(seen$calls, function(x) {
    identical(x$n_lambda, 20L) && identical(x$lambda_min_ratio, 1e-3)
  }, logical(1L))))
  expect_identical(
    fit$certificate$benchmark$glmnet$path_grid$s,
    exp(seq(0, log(1e-3), length.out = 20L))
  )
})


test_that("elastic_net(lambda = ) does not move the optimiser's chosen point", {
  skip_if_not_installed("glmnet")
  d <- fbpg_data()
  fits <- lapply(c(0.05, 5), function(lambda) {
    fbpg_select(
      d,
      fb_penalty(
        term(everything()) ~ elastic_net(lambda = lambda, alpha = 0.5)
      ),
      fb_decision(method = "full", metric = "rmse", n_lambda = 6L)
    )
  })
  sel <- lapply(fits, function(f) f$benchmark$glmnet$selection)
  expect_identical(sel[[1L]]$selected_s, sel[[2L]]$selected_s)
  expect_identical(sel[[1L]]$selected_terms, sel[[2L]]$selected_terms)
  expect_identical(sel[[1L]]$scores, sel[[2L]]$scores)
})
