# =============================================================================
# fb_select(): the certificate says how the hyperparameters were treated and
# how the data were reused by the selection.
#
# Why the fields exist. Ridge, Laplace and elastic-net reference fits are
# Bayesian conditional on penalty hyperparameters the caller fixed; they do
# not learn the global penalty strength jointly. The regularised horseshoe and
# R2D2 sample their global and local scales with the coefficients. An
# optimiser's lambda is tuned on the folds. The four axes did not record that
# distinction, and "fully Bayesian" read the same for all of them.
#
# The expected treatments come from that definition applied to each family's
# declared hyperparameter mode in the family table, not from the mapping under
# test. No sampler is involved.
# =============================================================================

fbhp_data <- function() {
  set.seed(8L)
  d <- data.frame(site = factor(rep(sprintf("S%02d", seq_len(6L)), each = 4L)))
  for (v in c("x1", "x2", "x3", "x4")) {
    d[[v]] <- stats::rnorm(nrow(d))
  }
  d$y <- stats::rnorm(nrow(d))
  d
}

fbhp_blocks <- function(penalty) {
  .fb_pen_resolve(penalty, .fb_sel_term_map(~ x1 + x2 + x3 + x4, fbhp_data()))
}

fbhp_certificate <- function(blocks,
                             workflow = "bayes",
                             prior = "ridge",
                             decision = "full",
                             projection = NULL) {
  base <- structure(
    list(
      fb = list(
        family = "gaussian", link = "identity", source = "test",
        fixed_terms = list(), random_terms = list(), response = "y"
      ),
      data = fbhp_data()
    ),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
  .fb_select_certificate(
    base = base,
    axes = list(
      workflow = workflow, provider = "brms", prior = prior,
      decision = decision
    ),
    cv = list(
      folds = list(), aggregated = NULL, prediction = "new_cluster",
      seed = 1L, metric = "elpd", covariance_estimand = "integrated"
    ),
    resamples = fb_resamples(unit = ~site, v = 2L),
    candidates = ~ x1 + x2 + x3 + x4,
    blocks = blocks,
    projection = projection
  )
}

fbhp_vocabulary <- c(
  "fixed_user", "calibrated_rule", "tuned_external", "hierarchical"
)


test_that("fixed-hyperparameter families are certified conditional on the caller's values", {
  for (penalty in list(
    fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    fb_penalty(term(everything()) ~ laplace(scale = 0.4)),
    fb_penalty(term(everything()) ~ elastic_net(lambda = 0.1, alpha = 0.5))
  )) {
    blocks <- fbhp_blocks(penalty)
    # The family table's own declaration, which the certificate must follow.
    expect_identical(blocks[[1L]]$hyperparameter_mode, "fixed")
    cert <- fbhp_certificate(blocks, prior = blocks[[1L]]$family)
    expect_identical(cert$hyperparameter_treatment, "fixed_user")
    expect_identical(
      cert$prior_blocks[[1L]]$hyperparameter_treatment, "fixed_user"
    )
  }
})


test_that("global-local families are certified hierarchical", {
  for (penalty in list(
    fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)),
    fb_penalty(term(everything()) ~ r2d2())
  )) {
    blocks <- fbhp_blocks(penalty)
    expect_identical(blocks[[1L]]$hyperparameter_mode, "learned")
    cert <- fbhp_certificate(blocks, prior = blocks[[1L]]$family)
    expect_identical(cert$hyperparameter_treatment, "hierarchical")
    expect_true(all(cert$hyperparameter_treatment %in% fbhp_vocabulary))
  }
})


test_that("an optimiser's lambda is tuned, and a run with no prior has no treatment", {
  blocks <- fbhp_blocks(fb_penalty(term(everything()) ~ ridge(scale = 0.5)))
  tuned <- fbhp_certificate(
    blocks, workflow = "penalized_optimizer", decision = "optimizer_support"
  )
  expect_identical(tuned$hyperparameter_treatment, "tuned_external")

  none <- fbhp_certificate(NULL, prior = "none")
  expect_identical(none$hyperparameter_treatment, character(0))
})


test_that("the data reuse names what the decision was chosen on", {
  blocks <- fbhp_blocks(fb_penalty(term(everything()) ~ laplace(scale = 0.4)))

  full <- fbhp_certificate(blocks, prior = "laplace", decision = "full")
  expect_identical(full$selection_data_reuse$level, "none")

  threshold <- fbhp_certificate(
    blocks, prior = "laplace", decision = "posterior_threshold"
  )
  expect_identical(
    threshold$selection_data_reuse$level, "decision_reads_reference_posterior"
  )

  searched <- fbhp_certificate(
    blocks, prior = "laplace", decision = "projection",
    projection = list(settings = list(validate_search = TRUE))
  )
  expect_identical(
    searched$selection_data_reuse$level, "decision_on_scored_folds"
  )
  expect_match(
    searched$selection_data_reuse$reason, "not a nested estimate", fixed = TRUE
  )
  optimistic <- fbhp_certificate(
    blocks, prior = "laplace", decision = "projection",
    projection = list(settings = list(validate_search = FALSE))
  )
  expect_match(
    optimistic$selection_data_reuse$reason, "optimistic", fixed = TRUE
  )

  tuned <- fbhp_certificate(
    blocks, workflow = "penalized_optimizer", decision = "optimizer_support"
  )
  expect_identical(tuned$selection_data_reuse$level, "tuning_on_scored_folds")

  protected <- fbhp_certificate(NULL, prior = "none")
  expect_identical(protected$selection_data_reuse$level, "none")

  out <- paste(utils::capture.output(print(searched)), collapse = "\n")
  expect_match(out, "hyperparameters: fixed_user", fixed = TRUE)
  expect_match(out, "decision_on_scored_folds", fixed = TRUE)
})
