# =============================================================================
# flexybayes(known_variance = K) with K a known covariance matrix: stage two
# of a two-stage analysis whose first stage returned the covariance of its
# trial means, Var(m) = sigma^2 I + sigma_l^2 Z Z' + K. INLA fits K as a
# generic0 effect with precision K^-1; brms refuses. fb_select(covariance =
# "base_model") whitens by the full V.
#
# Oracles the implementation did not write.
#   - GLS written out: V = s2 I + l2 Z Z' + K formed densely, the candidate
#     columns centred and scaled on each fold's training rows with mean() and
#     sd(), beta = solve(X' V^-1 X, X' V^-1 m) by solve().
#   - The whitened and protected-removed design written out: L = t(chol(V_TT)),
#     L^-1 by solve(), Q0 = qr.Q(qr(L^-1 X0)), y - Q0 Q0' y.
#   - The vector route: a diagonal K is the vector model under a second
#     lowering, and the two fits agree.
#   - REML written out for the INLA fit's variances.
#
# Tolerances, fixed before any result.
#   - 1e-8 relative for the written-out GLS and design (the acceptance bound
#     of the two-stage plan). Whitened QR against solve() on the normal
#     equations differ by rounding of order kappa n eps.
#   - 1e-4 relative between the diagonal-matrix and vector INLA fits (ten
#     times the largest spread between identical INLA fits, U-COV-MODEL
#     probe 04).
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

fbkc_quiet <- function(envir = parent.frame()) {
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

# Forty trial means in twelve location clusters, six trial-level candidates
# of which x1 carries signal, and a known covariance with off-diagonal
# entries of the size a first stage with shared varieties gives.
fbkc_trials <- function(seed = 7L, n = 40L, n_cluster = 12L) {
  .fb_sel_with_seed(seed, {
    cluster <- factor(sprintf("c%02d", rep(seq_len(n_cluster),
                                             length.out = n)))
    x <- matrix(stats::rnorm(n * 6L), n, 6L,
                dimnames = list(NULL, paste0("x", 1:6)))
    a <- matrix(stats::rnorm(n * n, 0, 0.05), n, n)
    k <- diag(stats::runif(n, 0.05, 0.3)) + tcrossprod(a)
    l <- stats::rnorm(n_cluster, 0, 0.6)[as.integer(cluster)]
    m <- 2 + 0.8 * x[, 1L] + l + stats::rnorm(n, 0, 0.4) +
      drop(t(chol(k)) %*% stats::rnorm(n))
    list(d = data.frame(m = m, cluster = cluster, x), k = k)
  })
}

fbkc_rel <- function(a, b) {
  max(abs(a - b)) / max(abs(b), .Machine$double.eps)
}

# The fold covariance of every fold, fixed, so V is known exactly.
fbkc_components <- list(residual = 0.3, terms = c(cluster = 0.5))

fbkc_fixed_fold_refit <- function(base, split, provider, engine = NULL,
                                  seed = NULL, slots = list()) {
  list(
    residual = fbkc_components$residual,
    terms = fbkc_components$terms,
    rho = list(),
    rows_passed = nrow(base$data)
  )
}

fbkc_v <- function(d, k) {
  fbkc_components$residual * diag(nrow(d)) +
    fbkc_components$terms[["cluster"]] * outer(d$cluster, d$cluster, "==") +
    k
}


# --- (i) the argument ------------------------------------------------------ #

test_that("a matrix that is not a covariance over the rows refuses by name", {
  fbkc_quiet()
  tr <- fbkc_trials()
  d <- tr$d
  k <- tr$k
  fit_with <- function(value) {
    flexybayes(m ~ 1, random = ~cluster, data = d, known_variance = value,
               backend = "inla", return_code = TRUE)
  }
  asymmetric <- k
  asymmetric[1L, 2L] <- asymmetric[1L, 2L] + 0.01
  linked_zero <- k
  linked_zero[3L, 3L] <- 0
  bad <- list(
    k[, -1L],
    k[-1L, -1L],
    replace(k, 5L, NA_real_),
    asymmetric,
    k - diag(2, nrow(k)),
    linked_zero,
    matrix(as.character(k), nrow(k))
  )
  for (value in bad) {
    expect_error(fit_with(value),
                 class = "flexybayes_refusal_known_variance_invalid")
  }
})

test_that("a matrix is carried as its diagonal and a row index, and follows its rows", {
  tr <- fbkc_trials()
  fb <- fb_from_asreml(m ~ 1, random = ~cluster, data = tr$d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  attached <- .fb_known_variance_attach(fb, tr$d, tr$k, aggregate = "auto")
  expect_identical(attached$data$flexybayes_known_variance, diag(tr$k))
  expect_identical(attached$data$flexybayes_known_row, seq_len(nrow(tr$d)))
  expect_identical(
    .fb_known_covariance_values(attached$fb, attached$data),
    tr$k
  )
  # Dropping and reordering rows takes the matching rows and columns; a row
  # the design grid adds has no index and no known part.
  rows <- c(9L, 2L, 30L)
  kept <- attached$data[rows, ]
  expect_identical(
    .fb_known_covariance_values(attached$fb, kept),
    tr$k[rows, rows]
  )
  grown <- rbind(kept, kept[1L, ])
  grown <- .fb_known_variance_fill(attached$fb, grown, 1L)
  got <- .fb_known_covariance_values(attached$fb, grown)
  expect_identical(got[1:3, 1:3], tr$k[rows, rows])
  expect_true(all(got[4L, ] == 0) && all(got[, 4L] == 0))
  expect_identical(
    .fb_known_variance_reissue(attached$fb, kept),
    tr$k[rows, rows]
  )
  # A vector re-issues as the vector.
  vec <- .fb_known_variance_attach(fb, tr$d, diag(tr$k), aggregate = "auto")
  expect_identical(.fb_known_variance_reissue(vec$fb, vec$data), diag(tr$k))
})

test_that("INLA fits a generic0 effect over K^-1, and brms refuses a matrix", {
  skip_if_not_installed("INLA")
  fbkc_quiet()
  tr <- fbkc_trials()
  code <- flexybayes(m ~ 1, random = ~cluster, data = tr$d,
                     known_variance = tr$k, backend = "inla",
                     return_code = TRUE)
  form <- gsub("\\s+", " ", paste(deparse(code$formula), collapse = ""))
  expect_match(
    form,
    paste0("f(flexybayes_known_index, model = \"generic0\", Cmatrix = ",
           "flexybayes_known_precision, hyper = list(prec = ",
           "list(initial = 0, fixed = TRUE)))"),
    fixed = TRUE
  )
  skip_if_not_installed("brms")
  expect_error(
    flexybayes(m ~ 1, random = ~cluster, data = tr$d,
               known_variance = tr$k, backend = "brms", return_code = TRUE),
    class = "flexybayes_refusal_known_variance_not_supported"
  )
})


# --- (ii) the fit ---------------------------------------------------------- #

test_that("a diagonal matrix fits the vector model, and update() re-issues the matrix", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbkc_quiet()
  tr <- fbkc_trials()
  fit <- function(value) {
    flexybayes(m ~ 1, random = ~cluster, data = tr$d,
               known_variance = value, backend = "inla", verbose = FALSE)
  }
  by_vector <- fit(diag(tr$k))
  by_matrix <- fit(diag(diag(tr$k)))
  vc_v <- summary(by_vector)$varcomp
  vc_m <- summary(by_matrix)$varcomp
  expect_lt(fbkc_rel(vc_m$estimate, vc_v$estimate), 1e-4)
  expect_lt(fbkc_rel(vc_m$std.error, vc_v$std.error), 1e-4)
  expect_lt(fbkc_rel(coef(by_matrix), coef(by_vector)), 1e-4)
  # fitted() leaves the known part out on both lowerings.
  expect_lt(fbkc_rel(fitted(by_matrix), fitted(by_vector)), 1e-4)

  full <- fit(tr$k)
  again <- update(full, verbose = FALSE)
  expect_identical(again$extras$call_info$known_variance, tr$k)
  expect_lt(fbkc_rel(coef(again), coef(full)), 1e-4)
})

test_that("the INLA fit with a full K recovers the REML variances of sigma^2 I + l2 Z Z' + K", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbkc_quiet()
  tr <- fbkc_trials()
  d <- tr$d
  zz <- outer(d$cluster, d$cluster, "==") * 1
  x <- matrix(1, nrow(d), 1L)
  nll <- function(th) {
    v <- exp(th[1L]) * diag(nrow(d)) + exp(th[2L]) * zz + tr$k
    r <- chol(v)
    xw <- backsolve(r, x, transpose = TRUE)
    yw <- backsolve(r, d$m, transpose = TRUE)
    a <- crossprod(xw)
    b <- solve(a, crossprod(xw, yw))
    0.5 * (2 * sum(log(diag(r))) + as.numeric(determinant(a)$modulus) +
             sum((yw - xw %*% b)^2))
  }
  o <- stats::optim(c(log(0.2), log(0.3)), nll, method = "BFGS",
                    control = list(reltol = 1e-14, maxit = 1000L))
  expect_identical(o$convergence, 0L)
  reml <- exp(o$par)

  fit <- flexybayes(m ~ 1, random = ~cluster, data = d,
                    known_variance = tr$k, backend = "inla", verbose = FALSE)
  vc <- summary(fit)$varcomp
  low <- vc$conf.low^2
  high <- vc$conf.high^2
  expect_true(all(reml >= low[c(1L, 2L)] & reml <= high[c(1L, 2L)]))
  # Fitted with K's off-diagonal part ignored, the residual absorbs it: the
  # two fits differ, so the generic0 effect is not a diagonal one.
  diagonal <- flexybayes(m ~ 1, random = ~cluster, data = d,
                         known_variance = diag(tr$k), backend = "inla",
                         verbose = FALSE)
  expect_gt(fbkc_rel(coef(fit), coef(diagonal)), 1e-3)
})


# --- (iii) fb_select() on a known covariance ------------------------------- #

test_that("each fold refit is handed the known covariance of every row", {
  tr <- fbkc_trials()
  fb <- fb_from_asreml(m ~ 1, random = ~cluster, data = tr$d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  attached <- .fb_known_variance_attach(fb, tr$d, tr$k, aggregate = "auto")
  base <- structure(
    list(
      fb = attached$fb,
      data = attached$data,
      extras = list(call_info = list(
        fixed = m ~ 1, random = ~cluster, family = "gaussian"
      ))
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
  split <- list(label = "1", train_rows = 1:30, assess_rows = 31:40)
  seen <- NULL
  local_mocked_bindings(flexybayes = function(...) {
    seen <<- list(...)
    stop("captured")
  })
  expect_error(
    .fb_sel_cov_fold_refit(base, split, provider = "inla", slots = list()),
    "captured"
  )
  expect_identical(seen$known_variance, tr$k)
  expect_false(any(.fb_known_variance_columns() %in% names(seen$data)))
  expect_true(all(is.na(seen$data$m[split$assess_rows])))
})

test_that("fb_select() whitens by s2 I + l2 Z Z' + K and its GLS refit and glmnet design equal the written-out ones to 1e-8", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  fbkc_quiet()
  tr <- fbkc_trials()
  d <- tr$d
  base <- flexybayes(m ~ 1, random = ~cluster, data = d,
                     known_variance = tr$k, backend = "inla", verbose = FALSE)
  local_mocked_bindings(.fb_sel_cov_fold_refit = fbkc_fixed_fold_refit)
  candidates <- paste0("x", 1:6)
  sel <- fb_select(
    base, ~ x1 + x2 + x3 + x4 + x5 + x6,
    penalty = fb_penalty(term(everything()) ~
                           elastic_net(lambda = 1, alpha = 0.6)),
    workflow = "penalized_optimizer", covariance = "base_model",
    decision = fb_decision(method = "full", metric = "rmse", n_lambda = 10L),
    resamples = fb_resamples(unit = ~cluster, v = 12L), seed = 1L,
    keep_fold_data = TRUE
  )
  expect_identical(sel$certificate$covariance$whitening$mode, "dense")
  expect_true(sel$certificate$covariance$known_variance)

  v <- fbkc_v(d, tr$k)
  folds <- fb_select_fold_data(sel)
  supports <- list(character(0), c("x1", "x3"), candidates)
  refits <- lapply(supports, function(s) fb_select_refit(sel, s))
  gap_design <- 0
  gap_beta <- 0
  gap_score <- 0
  for (fold in folds) {
    tr_rows <- fold$train_rows
    as_rows <- fold$assess_rows
    # The candidates standardised on the training rows.
    raw <- as.matrix(d[, candidates])
    centre <- colMeans(raw[tr_rows, ])
    spread <- apply(raw[tr_rows, ], 2L, stats::sd)
    x_std <- sweep(sweep(raw, 2L, centre, "-"), 2L, spread, "/")

    # Plan 4.3: whitening and the protected block removed.
    l_inv <- solve(t(chol(v[tr_rows, tr_rows])))
    y_w <- drop(l_inv %*% d$m[tr_rows])
    x0_w <- l_inv %*% matrix(1, length(tr_rows), 1L)
    x_w <- l_inv %*% x_std[tr_rows, ]
    q0 <- qr.Q(qr(x0_w))
    y_perp <- y_w - drop(q0 %*% crossprod(q0, y_w))
    x_perp <- x_w - q0 %*% crossprod(q0, x_w)
    gap_design <- max(
      gap_design,
      fbkc_rel(fold$whitened$y, y_w),
      fbkc_rel(fold$whitened$x_candidate, x_w),
      fbkc_rel(fold$reduced$y, y_perp),
      fbkc_rel(fold$reduced$x, x_perp)
    )

    # Plan 4.4: the unpenalised refit of each support, and its scores.
    v_tt_inv <- solve(v[tr_rows, tr_rows])
    for (i in seq_along(supports)) {
      x_s <- cbind(1, x_std[, supports[[i]], drop = FALSE])
      beta <- solve(
        crossprod(x_s[tr_rows, ], v_tt_inv %*% x_s[tr_rows, ]),
        crossprod(x_s[tr_rows, ], v_tt_inv %*% d$m[tr_rows])
      )
      got <- refits[[i]]$coefficients
      got <- got$estimate[got$fold == fold$fold]
      gap_beta <- max(gap_beta, fbkc_rel(got, drop(beta)))
      r <- d$m[as_rows] - drop(x_s[as_rows, , drop = FALSE] %*% beta)
      scores <- refits[[i]]$scores
      scores <- scores[scores$fold == fold$fold, ]
      for (key in scores$cluster) {
        at <- which(d$cluster[as_rows] == key)
        v_c <- v[as_rows[at], as_rows[at], drop = FALSE]
        nlpd <- 0.5 * (length(at) * log(2 * pi) +
                         as.numeric(determinant(v_c)$modulus) +
                         sum(r[at] * solve(v_c, r[at]))) / length(at)
        row <- scores[scores$cluster == key, ]
        gap_score <- max(
          gap_score,
          fbkc_rel(row$mse, mean(r[at]^2)),
          fbkc_rel(row$nlpd, nlpd)
        )
      }
    }
  }
  expect_lt(gap_design, 1e-8)
  expect_lt(gap_beta, 1e-8)
  expect_lt(gap_score, 1e-8)

  # The path solved again at the run's alpha is the run's path.
  path <- fb_select_path(sel)
  pkg <- as.data.frame(sel$benchmark$glmnet$support$by_lambda)
  pkg <- pkg[pkg$arm == "benchmark", ]
  pkg$support <- vapply(
    strsplit(pkg$terms, ", ", fixed = TRUE),
    function(t) paste(sort(t[nzchar(t)]), collapse = ";"),
    character(1L)
  )
  both <- merge(path, pkg[, c("fold", "lambda_index", "s", "support")],
                by = c("fold", "lambda_index"))
  expect_identical(nrow(both), nrow(path))
  expect_identical(both$support.x, both$support.y)
  expect_identical(both$s.x, both$s.y)
  expect_identical(attr(path, "alpha"), 0.6)
})

test_that("fb_select() with the fold refits run whitens by the refit's variances plus K", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_long_tier()
  fbkc_quiet()
  tr <- fbkc_trials()
  d <- tr$d
  base <- flexybayes(m ~ 1, random = ~cluster, data = d,
                     known_variance = tr$k, backend = "inla", verbose = FALSE)
  sel <- fb_select(base, ~ x1 + x2 + x3 + x4 + x5 + x6,
                   workflow = "unpenalized_benchmark",
                   covariance = "base_model",
                   resamples = fb_resamples(unit = ~cluster, v = 4L),
                   seed = 1L, keep_fold_data = TRUE)
  cov <- sel$certificate$covariance
  labels <- vapply(cov$folds, function(r) r$fold, character(1L))
  for (fold in fb_select_fold_data(sel)) {
    record <- cov$folds[[match(fold$fold, labels)]]
    v <- record$residual * diag(nrow(d)) +
      record$terms[["cluster"]] * outer(d$cluster, d$cluster, "==") + tr$k
    expect_lt(
      fbkc_rel(fold$v_assess, v[fold$assess_rows, fold$assess_rows]),
      1e-12
    )
    l_inv <- solve(t(chol(v[fold$train_rows, fold$train_rows])))
    expect_lt(
      fbkc_rel(fold$whitened$y, drop(l_inv %*% d$m[fold$train_rows])),
      1e-8
    )
  }
})


# --- (iv) the Bayesian refit of a fold's support ---------------------------- #

test_that("the INLA refit of a support adds each held-out row's known variance to its predictive variance", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbkc_quiet()
  tr <- fbkc_trials()
  base <- flexybayes(m ~ 1, random = ~cluster, data = tr$d,
                     known_variance = tr$k, backend = "inla", verbose = FALSE)
  held <- tr$d$cluster %in% c("c10", "c11", "c12")
  split <- list(
    label = "1", repeat_id = 1L,
    train_rows = which(!held), assess_rows = which(held),
    not_estimable_rows = integer(0), unobserved_rows = integer(0)
  )
  inputs <- .fb_sel_bayes_fold_inputs(base$data, split, response = "m",
                                      protected = ~1, random = "cluster")
  fit <- .fb_sel_bayes_refit(base, inputs$engine_data, m ~ 1, "inla",
                             na_action = "augment",
                             carry_known_variance = TRUE)
  # The refit carries the known covariance of every engine row.
  expect_identical(
    .fb_known_covariance_values(fit$fb %||% fit$extras$fb_terms, fit$data),
    tr$k
  )
  block <- list(design = matrix(numeric(0), length(split$train_rows), 0L),
                key = NULL, engine_names = list(inla = character(0)),
                design_columns = character(0))
  score <- function(inputs) {
    .fb_sel_value_fold_inla(fit, inputs, split, ~1, character(0),
                            "new_cluster", 400L, 11L, block, NA_character_)
  }
  # Both scorings read one set of posterior draws.
  sample_draws <- .fb_inla_posterior_sample
  kept <- NULL
  local_mocked_bindings(.fb_inla_posterior_sample = function(...) {
    if (is.null(kept)) {
      kept <<- sample_draws(...)
    }
    kept
  })
  with_known <- score(inputs)
  expect_identical(with_known$method, "analytic_gaussian")
  without <- inputs
  without$engine_data$flexybayes_known_variance <- NULL
  without_known <- score(without)
  # Same draws, so the predictive variances differ by K_ii exactly.
  expect_lt(
    fbkc_rel(with_known$scored$sd^2 - without_known$scored$sd^2,
             diag(tr$k)[split$assess_rows]),
    1e-10
  )
  expect_identical(with_known$scored$mean, without_known$scored$mean)
})

test_that("fb_select(refit = bayes) runs on an INLA base with a known covariance", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  skip_if_not_long_tier()
  fbkc_quiet()
  tr <- fbkc_trials()
  base <- flexybayes(m ~ 1, random = ~cluster, data = tr$d,
                     known_variance = tr$k, backend = "inla", verbose = FALSE)
  sel <- fb_select(
    base, ~ x1 + x2 + x3 + x4 + x5 + x6,
    penalty = fb_penalty(term(everything()) ~
                           elastic_net(lambda = 1, alpha = 1)),
    workflow = "penalized_optimizer", covariance = "base_model",
    decision = fb_decision(method = "full", metric = "rmse", n_lambda = 10L),
    refit = list(method = "bayes", n_draws = 400L),
    resamples = fb_resamples(unit = ~cluster, v = 4L), seed = 1L
  )
  rows <- sel$value$rows
  expect_identical(nrow(rows), nrow(tr$d))
  expect_true(all(rows$method == "analytic_gaussian"))
  # The known part is a floor under each row's predictive variance.
  expect_true(all(rows$predictive_sd^2 > diag(tr$k)[rows$row]))
})
