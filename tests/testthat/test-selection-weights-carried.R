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
    list(workflow = "penalized_optimizer", covariance = "reml_plugin",
         penalty = penalty, decision = decision),
    list(workflow = "penalized_optimizer", covariance = "base_model",
         penalty = penalty, decision = decision),
    list(workflow = "unpenalized_benchmark", covariance = "reml_plugin"),
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
  base <- list(extras = list(call_info = list(weights = rep(1, 4L))))
  expect_error(
    flexyBayes:::.fb_sel_bayes_refit(
      base, data.frame(y = 1:4), y ~ 1, "inla"
    ),
    class = "flexybayes_refusal_selection_weights_not_carried"
  )
})
