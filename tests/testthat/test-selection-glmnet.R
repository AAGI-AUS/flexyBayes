# =============================================================================
# The penalised-optimiser provider: measured conventions and the path
# contract.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 2, 8.1, 8.2,
# 8.3 and 10; SPEC_penalised_selection_v3_2026-09-11.md sections 7 and 8.
#
# Nothing here reads a convention out of `glmnet`'s documentation. The
# provider's lambda maps are measured against a Cholesky closed form and a
# subgradient residual this package writes, and the measurement is asserted
# in both directions: the rule that reproduces the closed form is pinned, and
# the rules that do not are asserted to fail. A pin that only ever sees the
# right answer cannot tell that it is the right answer.
#
# Injected defects that turn these gates red, all standing in the file:
#   D4  the naive L2 map, `lambda_canonical / n` with no response scale,
#       which leaves a KKT residual four orders above the measured map
#   D6  a lasso path on the same folds as the ridge path, so the ridge
#       assertion "no support, empty selected set" is not vacuous
#   D7  `glmnet` reported absent, through the one binding that asks
# =============================================================================

# --- fixtures ---------------------------------------------------------- #

fb_glmnet_data <- function(n_cluster = 10L, per = 6L, seed = 51L) {
  .fb_sel_with_seed(seed, {
    n <- n_cluster * per
    site <- factor(rep(seq_len(n_cluster), each = per))
    u <- stats::rnorm(n_cluster, sd = 1.4)
    d <- data.frame(
      site = site,
      x = stats::rnorm(n)
    )
    for (j in seq_len(5L)) {
      d[[paste0("z", j)]] <- stats::rnorm(n)
    }
    d$y <- 3 +
      0.8 * d$x +
      1.1 * d$z1 -
      0.9 * d$z2 +
      u[as.integer(site)] +
      stats::rnorm(n, sd = 0.5)
    d
  })
}


fb_glmnet_substrates <- function(data = fb_glmnet_data(), v = 4L, seed = 77L) {
  folds <- .fb_sel_folds(fb_resamples(unit = ~site, v = v), data, seed = seed)
  substrates <- lapply(folds, function(split) {
    .fb_sel_fold_substrate(
      data,
      split,
      response = "y",
      protected = ~x,
      candidates = ~ z1 + z2 + z3 + z4 + z5,
      random = "site"
    )
  })
  names(substrates) <- vapply(folds, function(f) f$label, character(1L))
  list(data = data, folds = folds, substrates = substrates)
}


fb_glmnet_block <- function(family, hyper, columns = 1:5) {
  list(list(family = family, columns = columns, hyperparameters = hyper))
}


fb_glmnet_matrix <- function(n = 55L, p = 6L, seed = 8L) {
  .fb_sel_with_seed(seed, {
    A <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    colnames(A) <- paste0("g", seq_len(p))
    z <- 1.5 +
      drop(A %*% c(1.8, -1, 0, 0, 0.7, 0)) +
      stats::rnorm(n, sd = 0.6)
    list(A = A, z = z, n = n, p = p)
  })
}


# --- the measured convention -------------------------------------------- #

test_that("the convention measurement pins both lambda maps", {
  skip_if_not_installed("glmnet")
  convention <- .fb_sel_glmnet_convention(refresh = TRUE)

  # The L1 map is `lambda_canonical / n` with no response scale on it. The
  # measurement is two-sided: the "none" rule wins, and the rules that put a
  # response scale on the L1 part are asserted to lose.
  expect_identical(convention$l1, "lambda_canonical / n")
  expect_identical(convention$l1_scale_rule, "none")
  expect_lt(convention$residuals$l1_kkt, 1e-6)
  for (rule in c("rms", "sd_population", "sd_sample")) {
    expect_gt(convention$residuals$l1_kkt_by_rule[[rule]], 1e-2)
  }

  # The L2 map carries the response scale once, and the scale is the
  # UNCENTRED root mean square of the response, not its standard deviation.
  # The standard-deviation rules stay four orders of magnitude worse on a
  # response with a non-zero mean, which is what the internal design gives it.
  expect_identical(convention$response_scale_rule, "rms")
  expect_identical(convention$l2, "lambda_canonical * rms(response) / n")
  expect_lt(convention$residuals$l2_max_abs, 1e-8)
  expect_gt(convention$residuals$l2_max_abs_by_rule[["sd_population"]], 1e-6)
  expect_gt(convention$residuals$l2_max_abs_by_rule[["sd_sample"]], 1e-6)
  # D4. The naive map -- no response scale at all -- is off by a margin no
  # tolerance in this module would absorb.
  expect_gt(convention$residuals$l2_max_abs_by_rule[["none"]], 1e-3)

  expect_identical(
    convention$glmnet_version,
    as.character(utils::packageVersion("glmnet"))
  )
  expect_false(convention$settings$standardize)
  expect_false(convention$settings$intercept)
})


test_that("the convention is cached per session and per glmnet version", {
  skip_if_not_installed("glmnet")
  first <- .fb_sel_glmnet_convention(refresh = TRUE)
  second <- .fb_sel_glmnet_convention()
  expect_identical(first, second)
})


test_that("glmnet at the canonical map is the proximal mode, for lasso", {
  skip_if_not_installed("glmnet")
  fixture <- fb_glmnet_matrix()
  convention <- .fb_sel_convention_at(
    .fb_sel_glmnet_convention(),
    fixture$z
  )

  for (scale in c(0.05, 0.2, 0.8)) {
    penalty <- .fb_sel_canonical_penalty(
      "laplace",
      list(scale = scale),
      n = fixture$n
    )
    provider <- .fb_sel_glmnet_at(
      fixture$A,
      fixture$z,
      penalty$lambda_composite,
      penalty$alpha_composite,
      convention
    )
    internal <- .fb_sel_prox_mode(
      fixture$A,
      fixture$z,
      "laplace",
      list(scale = scale),
      tol = 1e-10
    )
    expect_identical(provider$route, "glmnet")
    expect_lt(provider$kkt, 1e-6)
    expect_equal(
      unname(provider$beta),
      unname(internal$beta),
      tolerance = 1e-6
    )
    # The L1 map is exactly `lambda_canonical / n`.
    expect_equal(
      provider$lambda_glmnet,
      penalty$lambda_composite / fixture$n
    )
  }
})


test_that("glmnet at the corrected map is the proximal mode, for elastic net", {
  skip_if_not_installed("glmnet")
  fixture <- fb_glmnet_matrix()
  convention <- .fb_sel_convention_at(
    .fb_sel_glmnet_convention(),
    fixture$z
  )

  for (alpha in c(0.25, 0.5, 0.85)) {
    penalty <- .fb_sel_canonical_penalty(
      "elastic_net",
      list(lambda = 0.15, alpha = alpha),
      n = fixture$n
    )
    provider <- .fb_sel_glmnet_at(
      fixture$A,
      fixture$z,
      penalty$lambda_composite,
      penalty$alpha_composite,
      convention
    )
    internal <- .fb_sel_prox_mode(
      fixture$A,
      fixture$z,
      "elastic_net",
      list(lambda = 0.15, alpha = alpha),
      tol = 1e-10
    )
    expect_lt(provider$kkt, 1e-6)
    expect_equal(
      unname(provider$beta),
      unname(internal$beta),
      tolerance = 1e-6
    )

    # D4 again, on the mixed penalty: the naive map keeps `alpha` and divides
    # by `n`, and the elastic net inherits the ridge quirk on its L2 part.
    naive <- .fb_sel_glmnet_fit(
      fixture$A,
      fixture$z,
      alpha = penalty$alpha_composite,
      lambda = penalty$lambda_composite / fixture$n,
      thresh = 1e-16,
      maxit = 1e6
    )
    naive_kkt <- .fb_sel_kkt_residual(
      fixture$A,
      fixture$z,
      naive,
      "elastic_net",
      list(lambda = 0.15, alpha = alpha)
    )
    expect_gt(naive_kkt, 1e-2)
    expect_gt(naive_kkt / provider$kkt, 1e3)
  }
})


# --- P6: elastic_net(alpha = 1), the lasso endpoint ---------------------- #
#
# Final specification B6. The oracle is glmnet's own KKT residual at the
# returned point, never a Bayesian summary. `alpha = 1` puts lambda_l2 at
# exactly zero, the same (l1, l2) coordinate laplace() already reaches at
# alpha_composite = 1 -- so the assertion below is the same shape as "glmnet
# at the canonical map is the proximal mode, for lasso" above, entered from
# elastic_net() instead of laplace().
test_that("glmnet at the canonical map is the proximal mode, for the elastic-net lasso endpoint (P6)", {
  skip_if_not_installed("glmnet")
  fixture <- fb_glmnet_matrix()
  convention <- .fb_sel_convention_at(
    .fb_sel_glmnet_convention(),
    fixture$z
  )

  for (lambda in c(0.05, 0.15, 0.4)) {
    penalty <- .fb_sel_canonical_penalty(
      "elastic_net",
      list(lambda = lambda, alpha = 1),
      n = fixture$n
    )
    expect_equal(penalty$lambda_l2, 0)
    expect_equal(penalty$alpha_composite, 1)

    provider <- .fb_sel_glmnet_at(
      fixture$A,
      fixture$z,
      penalty$lambda_composite,
      penalty$alpha_composite,
      convention
    )
    internal <- .fb_sel_prox_mode(
      fixture$A,
      fixture$z,
      "elastic_net",
      list(lambda = lambda, alpha = 1),
      tol = 1e-10
    )
    expect_identical(provider$route, "glmnet")
    expect_identical(provider$alpha_glmnet, 1)
    expect_lt(provider$kkt, 1e-6)
    expect_equal(
      unname(provider$beta),
      unname(internal$beta),
      tolerance = 1e-6
    )
  }
})


test_that("the lasso endpoint fits and certifies through the full glmnet path, against glmnet's own optimality conditions (P6)", {
  skip_if_not_installed("glmnet")
  built <- fb_glmnet_substrates()
  path <- .fb_sel_glmnet_path(
    built$substrates,
    fb_glmnet_block("elastic_net", list(lambda = 0.2, alpha = 1)),
    folds = built$folds,
    decision = list(n_lambda = 6L),
    seed = 11L
  )

  expect_identical(path$workflow, "penalized_optimizer")
  expect_identical(path$provider, "glmnet")
  expect_false(path$fully_bayesian)
  # The equivalence level is reached exactly as it is for any other alpha of
  # this family: alpha is a hyperparameter of elastic_net, not a distinct
  # family, and the level is decided by family and provider, never by the
  # value of alpha.
  expect_identical(
    path$equivalence_level$level,
    "map_exact_fixed_covariance"
  )
  expect_identical(path$equivalence_level$prior_family, "elastic_net")

  # The certificate: glmnet's own KKT residual at the returned point, on
  # every fold, never a comparison against a Bayesian fit.
  fixed_kkt <- path$kkt_by_lambda$kkt[path$kkt_by_lambda$arm == "fixed"]
  expect_length(fixed_kkt, length(built$folds))
  expect_lt(max(fixed_kkt), 1e-6)
  expect_true(all(is.finite(path$canonical_objective$canonical_objective)))
  expect_true(all(is.finite(path$provider_objective$provider_objective)))

  # The fixed arm's own coordinate is the declared alpha, unconditionally on
  # which route solved it: alpha = 1 on every fold.
  fixed_rows <- path$path_coordinates[path$path_coordinates$arm == "fixed", ]
  expect_true(nrow(fixed_rows) > 0)
  expect_true(all(fixed_rows$alpha_canonical == 1))
})


test_that("elastic_net(alpha = 1) names its Bayesian counterpart in the shipped documentation (P6)", {
  # Acceptance (final specification B6): the certified glmnet endpoint and
  # the Bayesian lasso() / laplace() route are never conflated. The help
  # text is the surface a reader of fb_penalty() actually sees, so this
  # reads the installed .Rd rather than the roxygen source.
  rd_path <- testthat::test_path("..", "..", "man", "fb_penalty.Rd")
  skip_if(!file.exists(rd_path), "man/fb_penalty.Rd is not in the tree")
  squish <- function(x) gsub("[[:space:]]+", " ", paste(x, collapse = " "))
  text <- squish(readLines(rd_path, warn = FALSE))

  expect_match(text, "lasso corner of the elastic", fixed = TRUE)
  expect_match(text, "laplace()", fixed = TRUE)
  expect_match(
    text,
    "posterior mean of a Laplace prior is not a lasso estimate",
    fixed = TRUE
  )
})


test_that("glmnet at the corrected map is the ridge closed form", {
  skip_if_not_installed("glmnet")
  fixture <- fb_glmnet_matrix()
  convention <- .fb_sel_convention_at(
    .fb_sel_glmnet_convention(),
    fixture$z
  )
  for (scale in c(0.4, 1.2)) {
    penalty <- .fb_sel_canonical_penalty(
      "ridge",
      list(scale = scale),
      n = fixture$n
    )
    provider <- .fb_sel_glmnet_at(
      fixture$A,
      fixture$z,
      penalty$lambda_composite,
      0,
      convention
    )
    closed <- .fb_sel_bridge_ridge(fixture$A, fixture$z, scale)
    expect_equal(unname(provider$beta), unname(closed), tolerance = 1e-6)
    expect_equal(
      provider$lambda_glmnet,
      penalty$lambda_composite * convention$response_scale / fixture$n
    )
  }
})


# --- the provider and family gates --------------------------------------- #

test_that("an absent provider refuses by name (defect D7)", {
  # D7. The provider probe is one binding, so the refusal can be exercised
  # without uninstalling anything and the registered code has a triggering
  # test rather than a description.
  testthat::local_mocked_bindings(
    .fb_sel_provider_installed = function(package) FALSE
  )
  expect_error(
    .fb_sel_require_glmnet(),
    class = "flexybayes_refusal_selection_provider_absent"
  )
  condition <- tryCatch(
    .fb_sel_require_glmnet(),
    flexybayes_refusal_selection_provider_absent = function(e) e
  )
  expect_identical(condition$provider, "glmnet")
  expect_true(grepl("install.packages", condition$message, fixed = TRUE))
})


test_that("a family with no optimiser provider refuses by name", {
  # The class-wide global-local families are Bayesian-only, and `adaptive()`
  # has no provider at all. The four grouped families left this list at
  # 0.13.0, when their proximal provider and its gates shipped
  # (tests/testthat/test-selection-grouped-penalties.R).
  for (family in c("adaptive", "regularized_horseshoe", "r2d2")) {
    expect_error(
      .fb_sel_glmnet_assert_family(family),
      class = "flexybayes_refusal_selection_penalty_family_withdrawn"
    )
  }
  for (family in c(
    "ridge", "laplace", "elastic_net",
    "group_lasso", "sparse_group_lasso", "exclusive_lasso", "group_exclusive"
  )) {
    expect_identical(.fb_sel_glmnet_assert_family(family), family)
  }
})


test_that("blocks no single glmnet problem expresses refuse by name", {
  skip_if_not_installed("glmnet")
  built <- fb_glmnet_substrates()
  mixed <- list(
    list(family = "ridge", columns = 1:2, hyperparameters = list(scale = 1)),
    list(
      family = "laplace",
      columns = 3:5,
      hyperparameters = list(scale = 1)
    )
  )
  expect_error(
    .fb_sel_glmnet_path(built$substrates, mixed, folds = built$folds),
    class = "flexybayes_refusal_selection_penalty_blocks_not_expressible"
  )

  same_family_two_scales <- list(
    list(
      family = "laplace",
      columns = 1:2,
      hyperparameters = list(scale = 1)
    ),
    list(
      family = "laplace",
      columns = 3:5,
      hyperparameters = list(scale = 0.2)
    )
  )
  expect_error(
    .fb_sel_glmnet_path(
      built$substrates,
      same_family_two_scales,
      folds = built$folds
    ),
    class = "flexybayes_refusal_selection_penalty_blocks_not_expressible"
  )

  # Two blocks that resolve to the same penalty are one problem, and are not
  # refused.
  agreeing <- list(
    list(
      family = "laplace",
      columns = 1:2,
      hyperparameters = list(scale = 0.5)
    ),
    list(
      family = "laplace",
      columns = 3:5,
      hyperparameters = list(scale = 0.5)
    )
  )
  expect_s3_class(
    .fb_sel_glmnet_path(
      built$substrates,
      agreeing,
      folds = built$folds,
      decision = list(benchmark = FALSE)
    ),
    "fb_sel_canonical_path"
  )
})


# --- the path on the built substrate -------------------------------------- #

test_that("a ridge path reports no support, and the assertion is not vacuous", {
  skip_if_not_installed("glmnet")
  built <- fb_glmnet_substrates()

  ridge <- .fb_sel_glmnet_path(
    built$substrates,
    fb_glmnet_block("ridge", list(scale = 0.6)),
    folds = built$folds,
    decision = list(n_lambda = 6L),
    seed = 5L
  )
  expect_identical(ridge$support_by_resolution$semantics, "none")
  expect_false(ridge$support_by_resolution$selects)
  expect_length(ridge$support_by_resolution$selected, 0L)
  expect_true(is.na(ridge$support_by_resolution$selection_resolution))
  expect_true(all(
    ridge$path_coordinates$path_parameter == "log_lambda"
  ))
  expect_true(all(ridge$path_coordinates$null_anchor == "infinite"))
  # Every coefficient a ridge fit returns is non-zero: the "no support"
  # statement is about the geometry, not about an accident of the grid.
  ridge_candidates <- ridge$coefficients_original_scale$estimate[
    ridge$coefficients_original_scale$role == "candidate"
  ]
  expect_true(all(ridge_candidates != 0))

  # D6. The same folds under a lasso do return a support, so the ridge
  # assertion above is a property of the family rather than of a path that
  # happened to select nothing.
  lasso <- .fb_sel_glmnet_path(
    built$substrates,
    fb_glmnet_block("laplace", list(scale = 0.15)),
    folds = built$folds,
    decision = list(n_lambda = 8L),
    seed = 5L
  )
  expect_identical(lasso$support_by_resolution$semantics, "coefficient")
  expect_true(lasso$support_by_resolution$selects)
  expect_gt(length(lasso$support_by_resolution$selected), 0L)
  expect_true(any(
    lasso$coefficients_original_scale$estimate[
      lasso$coefficients_original_scale$role == "candidate"
    ] ==
      0
  ))
})


test_that("the provider certificate carries every field the contract names", {
  skip_if_not_installed("glmnet")
  built <- fb_glmnet_substrates()
  path <- .fb_sel_glmnet_path(
    built$substrates,
    fb_glmnet_block("elastic_net", list(lambda = 0.2, alpha = 0.6)),
    folds = built$folds,
    decision = list(n_lambda = 6L),
    seed = 11L
  )

  required <- c(
    "coefficients_original_scale",
    "fitted_values",
    "support_by_resolution",
    "canonical_objective",
    "provider_objective",
    "kkt_by_lambda",
    "path_coordinates",
    "provider_version",
    "provider_lambda_convention",
    "workflow",
    "provider",
    "fully_bayesian",
    "equivalence_level"
  )
  expect_true(all(required %in% names(path)))
  for (field in required) {
    expect_false(is.null(path[[field]]))
  }

  expect_identical(path$workflow, "penalized_optimizer")
  expect_identical(path$provider, "glmnet")
  expect_false(path$fully_bayesian)
  expect_identical(
    path$provider_version,
    as.character(utils::packageVersion("glmnet"))
  )
  expect_identical(
    path$provider_lambda_convention$response_scale_rule,
    "rms"
  )
  expect_identical(
    path$equivalence_level$level,
    "map_exact_fixed_covariance"
  )
  expect_identical(path$contract, "SPEC_bayesian_selection_v5_2026-09-12")
  expect_identical(path$covariance_estimand, "reml_plugin_covariance")
  expect_identical(path$prediction, "new_cluster")

  # Every point on every fold carries its own residual, and the fixed arm --
  # the canonical point the declared hyperparameters name -- is inside the
  # tolerance the exact-limit tier asks for.
  fixed_kkt <- path$kkt_by_lambda$kkt[path$kkt_by_lambda$arm == "fixed"]
  expect_length(fixed_kkt, length(built$folds))
  expect_lt(max(fixed_kkt), 1e-6)
  expect_true(all(is.finite(path$canonical_objective$canonical_objective)))
  expect_true(all(is.finite(path$provider_objective$provider_objective)))
})


test_that("the benchmark arm decides at the fold unit and reports its rule", {
  skip_if_not_installed("glmnet")
  built <- fb_glmnet_substrates()
  path <- .fb_sel_glmnet_path(
    built$substrates,
    fb_glmnet_block("laplace", list(scale = 0.2)),
    folds = built$folds,
    decision = list(rule = "one_se", metric = "rmse", n_lambda = 8L),
    seed = 3L
  )
  selection <- path$selection
  expect_true(selection$ran)
  expect_identical(selection$rule, "one_se")
  expect_identical(selection$unit, "cluster")
  expect_equal(unique(selection$scores$n_units), length(built$folds))
  # The one-standard-error rule never picks a model more complex than the
  # minimum, which is the property that makes it a rule rather than a label.
  minimum_row <- which(selection$scores$model == selection$minimum_model)
  chosen_row <- which(selection$scores$model == selection$selected_model)
  expect_lte(
    selection$scores$complexity[[chosen_row]],
    selection$scores$complexity[[minimum_row]]
  )

  # Turning the benchmark off leaves the fixed arm and says so.
  fixed_only <- .fb_sel_glmnet_path(
    built$substrates,
    fb_glmnet_block("laplace", list(scale = 0.2)),
    folds = built$folds,
    decision = list(benchmark = FALSE)
  )
  expect_false(fixed_only$selection$ran)
  expect_true(all(fixed_only$kkt_by_lambda$arm == "fixed"))
})


test_that("the path never lets a protected column reach the provider", {
  skip_if_not_installed("glmnet")
  built <- fb_glmnet_substrates()
  # Two candidate columns are named by no block, so they are unpenalised and
  # join the protected side of the reduction. The provider payload's own
  # assertion is what holds this, and the coefficient table has to show the
  # roles it produced.
  path <- .fb_sel_glmnet_path(
    built$substrates,
    fb_glmnet_block("laplace", list(scale = 0.3), columns = c("z1", "z2")),
    folds = built$folds,
    decision = list(benchmark = FALSE)
  )
  candidates <- unique(path$coefficients_original_scale$term[
    path$coefficients_original_scale$role == "candidate"
  ])
  protected <- unique(path$coefficients_original_scale$term[
    path$coefficients_original_scale$role == "protected"
  ])
  expect_setequal(candidates, c("z1", "z2"))
  expect_true(all(c("x", "z3", "z4", "z5") %in% protected))
  expect_false("x" %in% candidates)
})


test_that("the path refuses folds that do not match its substrates", {
  skip_if_not_installed("glmnet")
  built <- fb_glmnet_substrates()
  expect_error(
    .fb_sel_glmnet_path(
      built$substrates,
      fb_glmnet_block("ridge", list(scale = 1)),
      folds = built$folds[1:2]
    ),
    "must describe the same folds",
    fixed = TRUE
  )
  expect_error(
    .fb_sel_glmnet_path(
      list(),
      fb_glmnet_block("ridge", list(scale = 1))
    ),
    "must be an <fb_sel_substrate>",
    fixed = TRUE
  )
})


test_that("the decision settings refuse a name they do not carry", {
  expect_error(
    .fb_sel_glmnet_decision(list(rule = "one_se", stability = TRUE)),
    "unknown decision setting",
    fixed = TRUE
  )
  expect_error(
    .fb_sel_glmnet_decision(list(n_lambda = 1L)),
    "at least two grid points",
    fixed = TRUE
  )
  defaults <- .fb_sel_glmnet_decision()
  expect_identical(defaults$rule, "one_se")
  expect_identical(defaults$metric, "rmse")
  expect_true(defaults$benchmark)
})
