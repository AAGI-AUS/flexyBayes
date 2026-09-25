# =============================================================================
# fb_select() on a base whose data carry a missing response.
#
# The unpenalised benchmark, the penalised optimiser and the comparison
# `benchmark = TRUE` adds leave a row with no response out of every fold's
# training and scored rows. For a Gaussian response that is the exact
# marginal over the row. The row stays in the data, so a field is still
# indexed over the whole padded lattice.
# =============================================================================

fbmr_data <- function(missing_rows = c(4L, 19L)) {
  set.seed(20260919L)
  J <- 12L
  m <- 6L
  d <- data.frame(site = factor(rep(sprintf("S%02d", seq_len(J)), each = m)))
  d$clay <- rep(stats::rnorm(J), each = m)
  d$rain <- rep(stats::rnorm(J), each = m)
  d$plot_z <- stats::rnorm(nrow(d))
  d$yield <- 3 + 0.8 * d$clay + 0.2 * d$plot_z +
    rep(stats::rnorm(J, sd = 0.6), each = m) + stats::rnorm(nrow(d), sd = 0.7)
  d$yield[missing_rows] <- NA
  d
}

fbmr_base <- function(d) {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  flexybayes(yield ~ plot_z, random = ~ site, data = d, backend = "inla",
             verbose = FALSE)
}

fbmr_folds <- function() fb_resamples(unit = ~ site, v = 3L)

# Six trials on a 4 x 6 lattice. Eleven cells have no response, two of them
# inside the lattice, so a covariance indexed over the observed cells alone
# would place every later cell of those trials at the wrong position.
fbmr_field_data <- function() {
  set.seed(20260925L)
  d <- expand.grid(ROW = factor(1:6), RANGE = factor(1:4), TRIAL = factor(1:6))
  d$x1 <- rep(stats::rnorm(6L), each = 24L)
  d$x2 <- rep(stats::rnorm(6L), each = 24L)
  d$HALF <- factor(paste0(d$TRIAL, ":", as.integer(d$RANGE) <= 2L))
  d$y <- 1 + 0.5 * d$x1 + stats::rnorm(6L)[as.integer(d$TRIAL)] +
    stats::rnorm(nrow(d))
  d$y[c(3L, 10L, 30L, 47L, 48L, 60L, 75L, 99L, 100L, 120L, 131L)] <- NA
  d
}

fbmr_field_base <- function(d) {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  suppressWarnings(suppressMessages(flexybayes(
    y ~ 1,
    random = ~ TRIAL + id(TRIAL):ar1(RANGE):ar1(ROW),
    data = d, backend = "inla", verbose = FALSE
  )))
}

# The components every fold refit reports, fixed so the covariance is known.
fbmr_components <- function(slots) {
  terms <- numeric(0)
  rho <- list()
  for (s in slots) {
    if (identical(s$structure$type, "ar1xar1")) {
      terms[[s$name]] <- 0.6
      rho[[s$name]] <- c(row = 0.7, col = 0.4)
    } else {
      terms[[s$name]] <- 0.8
    }
  }
  list(residual = 0.5, terms = terms, rho = rho)
}

# V written from the lattice coordinates, outside the package.
fbmr_closed_form_v <- function(d, rows_a, rows_b) {
  same <- outer(d$TRIAL[rows_a], d$TRIAL[rows_b], "==")
  r <- as.integer(d$RANGE)
  cc <- as.integer(d$ROW)
  field <- 0.7^abs(outer(r[rows_a], r[rows_b], "-")) *
    0.4^abs(outer(cc[rows_a], cc[rows_b], "-"))
  0.5 * outer(rows_a, rows_b, "==") + same * (0.8 + 0.6 * field)
}

fbmr_oracle_fold <- function(d, train, assess, columns) {
  x <- cbind(1, as.matrix(d[, columns, drop = FALSE]))
  v_tt <- fbmr_closed_form_v(d, train, train)
  v_st <- fbmr_closed_form_v(d, assess, train)
  v_ss <- fbmr_closed_form_v(d, assess, assess)
  a <- solve(v_tt)
  beta <- solve(t(x[train, , drop = FALSE]) %*% a %*% x[train, , drop = FALSE],
                t(x[train, , drop = FALSE]) %*% a %*% d$y[train])
  resid_t <- d$y[train] - drop(x[train, , drop = FALSE] %*% beta)
  mu <- drop(x[assess, , drop = FALSE] %*% beta) + drop(v_st %*% a %*% resid_t)
  list(residual = d$y[assess] - mu, v = v_ss - v_st %*% a %*% t(v_st))
}

fbmr_run_cv <- function(base, unit, prediction) {
  cm <- flexyBayes:::.fb_sel_cov_model(base, seed = 1L)
  components <- fbmr_components(cm$slots)
  testthat::local_mocked_bindings(
    .fb_sel_cov_fold_refit = function(base, split, provider, engine = NULL,
                                      seed = NULL, slots = list()) {
      components
    },
    .package = "flexyBayes"
  )
  flexyBayes:::.fb_sel_cv(
    data = base$data,
    response = base$fb$response,
    protected = ~ 1,
    candidates = ~ x1 + x2,
    random = flexyBayes:::.fb_sel_grouping_columns(base$fb)$columns,
    resamples = fb_resamples(unit = unit, v = 3L),
    prediction = prediction,
    metric = "rmse",
    cov_model = cm
  )
}

fbmr_expect_oracle <- function(cv, d, unit) {
  observed <- which(!is.na(d$y))
  rows <- as.data.frame(cv$rows)
  clusters <- as.data.frame(cv$clusters)
  for (split in cv$folds) {
    expect_true(all(split$train_rows %in% observed))
    expect_true(all(split$assess_rows %in% observed))
    for (model in c("protected_only", "full")) {
      columns <- if (identical(model, "full")) c("x1", "x2") else character(0)
      want <- fbmr_oracle_fold(d, split$train_rows, split$assess_rows, columns)
      got <- rows[rows$fold == split$label & rows$model == model, ]
      expect_equal(got$residual, want$residual, tolerance = 1e-10)
      got_clusters <- clusters[
        clusters$fold == split$label & clusters$model == model,
      ]
      key <- as.character(d[[unit]])[split$assess_rows]
      for (k in seq_len(nrow(got_clusters))) {
        i <- which(key == got_clusters$cluster[[k]])
        b <- want$v[i, i, drop = FALSE]
        r <- want$residual[i]
        nlpd <- 0.5 * (length(i) * log(2 * pi) +
          as.numeric(determinant(b)$modulus) +
          drop(t(r) %*% solve(b, r))) / length(i)
        expect_equal(got_clusters$nlpd[[k]], nlpd, tolerance = 1e-10)
      }
    }
  }
}


test_that("held-out trials on a padded field match the closed-form GLS", {
  d <- fbmr_field_data()
  base <- fbmr_field_base(d)
  cv <- fbmr_run_cv(base, ~ TRIAL, "new_cluster")
  expect_identical(nrow(cv$rows), 2L * sum(!is.na(d$y)))
  expect_setequal(
    unlist(lapply(cv$folds, `[[`, "unobserved_rows")),
    which(is.na(d$y))
  )
  fbmr_expect_oracle(cv, d, "TRIAL")
})


test_that("within-trial folds on a padded field match the closed-form GLS", {
  # Half of each trial is held out, so every prediction borrows from the
  # trial's training cells through the field, at their lattice positions.
  d <- fbmr_field_data()
  base <- fbmr_field_base(d)
  cv <- fbmr_run_cv(base, ~ HALF, "known_cluster")
  expect_identical(nrow(cv$rows), 2L * sum(!is.na(d$y)))
  fbmr_expect_oracle(cv, d, "HALF")
})


test_that("the REML plug-in scores a padded base as the base without those rows", {
  padded <- fbmr_base(fbmr_data())
  kept <- fbmr_data()
  kept <- kept[!is.na(kept$yield), ]
  base_kept <- fbmr_base(kept)
  pen <- fb_penalty(term(everything()) ~ ridge(scale = 1))
  for (workflow in c("unpenalized_benchmark", "penalized_optimizer")) {
    if (identical(workflow, "penalized_optimizer")) {
      skip_if_not_installed("glmnet")
    }
    run <- function(base) {
      fb_select(base, candidates = ~ clay + rain,
                penalty = if (identical(workflow, "penalized_optimizer")) pen,
                workflow = workflow, resamples = fbmr_folds(), seed = 1L)
    }
    a <- run(padded)
    b <- run(base_kept)
    expect_equal(a$scores$score, b$scores$score, tolerance = 1e-10)
    expect_equal(a$scores$se, b$scores$se, tolerance = 1e-10)
  }
})


test_that("the study route runs on a padded field and scores observed rows only", {
  skip_if_not_installed("glmnet")
  d <- fbmr_field_data()
  base <- fbmr_field_base(d)
  sel <- suppressWarnings(fb_select(
    base = base,
    candidates = ~ x1 + x2,
    penalty = fb_penalty(term(everything()) ~ elastic_net(lambda = 1, alpha = 1)),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~ TRIAL, v = 3L),
    prediction = "new_cluster",
    decision = fb_decision(
      method = "full", metric = "rmse", rule = "one_se",
      stability_subsamples = 0L
    ),
    covariance = "base_model",
    refit = list(method = "bayes"),
    seed = 1L
  ))
  expect_s3_class(sel, "fb_select")
  scored <- unlist(lapply(sel$folds, `[[`, "assess_rows"))
  expect_setequal(scored, which(!is.na(d$y)))
  expect_true(all(is.finite(sel$scores$score)))
})


test_that("a fold whose held-out units have no response refuses by name", {
  d <- fbmr_data(missing_rows = 1:6)
  folds <- flexyBayes:::.fb_sel_folds(
    fb_resamples(unit = ~ site, v = 12L), d
  )
  expect_error(
    flexyBayes:::.fb_sel_observed_folds(folds, d$yield),
    class = "flexybayes_refusal_selection_fold_not_usable"
  )
})


test_that("a base with every response observed is left as it was", {
  d <- fbmr_data(missing_rows = integer(0))
  folds <- flexyBayes:::.fb_sel_folds(fbmr_folds(), d)
  expect_identical(flexyBayes:::.fb_sel_observed_folds(folds, d$yield), folds)
})
