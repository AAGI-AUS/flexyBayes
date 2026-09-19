# =============================================================================
# flexybayes(known_variance = v): a known variance per row added to the
# residual, Var(y_i) = sigma^2 + v_i, the second stage of a two-stage
# analysis. brms fits y | se(sqrt(v), sigma = TRUE); INLA a row effect of
# fixed unit precision scaled by sqrt(v_i). fb_select(covariance =
# "base_model") carries v into the fold covariance's diagonal.
#
# Oracles the implementation did not write.
#   - REML. A restricted likelihood written out below and maximised with
#     optim(); on this fixture it agrees with metafor::rma.mv to 1e-5
#     (review/units/U-STAGE2-VAR_2026-09-19, probe 02).
#   - GLS. beta = (X' V^-1 X)^-1 X' V^-1 y with V = diag(sigma^2 + v) +
#     sigma_t^2 Z Z', formed densely and solved with solve(); the per-cluster
#     squared error and negative log predictive density are computed from it.
#   - The other engine. brms and INLA reach the model by different lowerings.
#
# Tolerances, fixed before any result.
#   - The simulated values and the REML estimate lie inside each engine's
#     95% credible interval of each variance: the interval is the engine's
#     own measured posterior spread.
#   - INLA against brms: |E_inla - E_brms| <= z * sqrt(mcse^2 + (1e-4 E)^2),
#     mcse the Monte Carlo standard error of the brms mean measured from its
#     draws. Across five brms seeds its spread matched the within-fit MCSE
#     (probe 03). 1e-4 relative is ten times the largest spread measured
#     between identical INLA fits (1.1e-5, U-COV-MODEL probe 04; 2.8e-7 on
#     this fixture, probe 03). z = 4 for each set of ten scalar comparisons,
#     5 for the per-row fitted values (220 rows).
#   - A vector of zeros: identical engine inputs; INLA estimates within 1e-4
#     relative (above); brms draws identical under one seed.
#   - GLS through fb_select(): 1.5e-8 relative (sqrt of machine epsilon).
#     Whitened QR against solve() on the normal equations differ by rounding
#     of order kappa n eps.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

fbkv_quiet <- function(envir = parent.frame()) {
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

# Variety means from 30 trials of 8 varieties, 20 cells absent: variety
# fixed, trial random (variance t2), and each mean with its own known
# variance, zero on 24 cells.
fbkv_stage2 <- function(seed = 20260919L, n_trial = 30L, n_var = 8L,
                        s2 = 0.25, t2 = 1, n_drop = 20L, n_zero = 24L) {
  .fb_sel_with_seed(seed, {
    d <- expand.grid(
      variety = paste0("v", seq_len(n_var)),
      trial = sprintf("t%02d", seq_len(n_trial)),
      stringsAsFactors = FALSE
    )
    alpha <- seq(-1, 1, length.out = n_var)
    s_t <- stats::rnorm(n_trial, 0, sqrt(t2))
    v <- stats::runif(nrow(d), 0.05, 1.5)
    v[sample.int(nrow(d), n_zero)] <- 0
    d$ybar <- alpha[match(d$variety, paste0("v", seq_len(n_var)))] +
      s_t[match(d$trial, sprintf("t%02d", seq_len(n_trial)))] +
      stats::rnorm(nrow(d), 0, sqrt(s2 + v))
    d$v <- v
    keep <- sort(sample.int(nrow(d), nrow(d) - n_drop))
    d <- d[keep, ]
    rownames(d) <- NULL
    d$variety <- factor(d$variety)
    d$trial <- factor(d$trial)
    d
  })
}

# Variety means with four trial-level candidates, one of which carries
# signal, for fb_select().
fbkv_select_data <- function(seed = 2L, n_trial = 24L, n_var = 6L) {
  .fb_sel_with_seed(seed, {
    d <- expand.grid(
      variety = paste0("v", seq_len(n_var)),
      trial = sprintf("t%02d", seq_len(n_trial)),
      stringsAsFactors = FALSE
    )
    z <- matrix(stats::rnorm(n_trial * 4L), n_trial)
    colnames(z) <- paste0("z", 1:4)
    d <- cbind(d, z[match(d$trial, sprintf("t%02d", seq_len(n_trial))), ])
    d$variety <- factor(d$variety)
    d$trial <- factor(d$trial)
    v <- stats::runif(nrow(d), 0.05, 1.5)
    v[sample.int(nrow(d), 12L)] <- 0
    d$ybar <- 0.3 * as.integer(d$variety) + 0.8 * d$z1 +
      stats::rnorm(n_trial, 0, 0.7)[as.integer(d$trial)] +
      stats::rnorm(nrow(d), 0, sqrt(0.25 + v))
    d$v <- v
    d
  })
}

# The restricted likelihood of y = X b + Z s + e, s ~ N(0, t2 I),
# e ~ N(0, diag(s2 / w + v)), maximised over (log s2, log t2).
fbkv_reml <- function(y, X, Z, v, w = 1) {
  nll <- function(th) {
    V <- diag(exp(th[1]) / w + v) + exp(th[2]) * tcrossprod(Z)
    L <- chol(V)
    Xw <- backsolve(L, X, transpose = TRUE)
    yw <- backsolve(L, y, transpose = TRUE)
    A <- crossprod(Xw)
    b <- solve(A, crossprod(Xw, yw))
    r <- yw - Xw %*% b
    0.5 * (2 * sum(log(diag(L))) +
             as.numeric(determinant(A)$modulus) + sum(r^2))
  }
  o <- stats::optim(c(log(0.3), log(0.8)), nll, method = "BFGS",
                    control = list(reltol = 1e-14, maxit = 1000L))
  stopifnot(o$convergence == 0L)
  c(sigma = exp(o$par[1]), sd_trial = exp(o$par[2]))
}

# Each variance's posterior mean and 95% interval from a fit's table.
fbkv_variances <- function(fit) {
  vc <- summary(fit)$varcomp
  data.frame(
    component = vc$component,
    mean = vc$estimate^2 + vc$std.error^2,
    low = vc$conf.low^2,
    high = vc$conf.high^2,
    stringsAsFactors = FALSE
  )
}

# The per-cluster scores of one fold, written out: GLS on the training rows
# with V = diag(s2 + v) + t2 (same trial), predictions X_S beta at new
# clusters, and each assessment trial's mean squared error and negative log
# predictive density under V_SS.
fbkv_gls_clusters <- function(d, v, split, record, formula) {
  s2 <- record$residual
  t2 <- record$terms[["trial"]]
  X <- stats::model.matrix(formula, d)
  vcov_rows <- function(rows) {
    diag(s2 + v[rows], length(rows)) +
      t2 * outer(d$trial[rows], d$trial[rows], "==")
  }
  tr <- split$train_rows
  as_ <- split$assess_rows
  vi_x <- solve(vcov_rows(tr), X[tr, , drop = FALSE])
  beta <- solve(crossprod(X[tr, , drop = FALSE], vi_x),
                crossprod(vi_x, d$ybar[tr]))
  r <- d$ybar[as_] - drop(X[as_, , drop = FALSE] %*% beta)
  clusters <- as.character(d$trial[as_])
  out <- lapply(unique(clusters), function(key) {
    at <- which(clusters == key)
    rows <- as_[at]
    V <- vcov_rows(rows)
    n <- length(rows)
    data.frame(
      cluster = key,
      mse = mean(r[at]^2),
      nlpd = 0.5 * (n * log(2 * pi) +
                      as.numeric(determinant(V)$modulus) +
                      sum(r[at] * solve(V, r[at]))) / n,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

fbkv_penalty <- function() {
  fb_penalty(term(everything()) ~ laplace(scale = 0.15))
}

fbkv_rel <- function(a, b) {
  max(abs(a - b) / pmax(abs(b), .Machine$double.eps))
}

# Compare fb_select()'s per-cluster scores with the written-out GLS, fold by
# fold, for both models the benchmark scores.
fbkv_expect_gls_scores <- function(sel, d) {
  cov <- sel$certificate$covariance
  clusters <- as.data.frame(sel$clusters)
  formulas <- list(
    protected_only = ~variety,
    full = ~ variety + z1 + z2 + z3 + z4
  )
  labels <- vapply(cov$folds, function(r) r$fold, character(1L))
  for (split in sel$folds) {
    record <- cov$folds[[match(split$label, labels)]]
    for (model in names(formulas)) {
      hand <- fbkv_gls_clusters(d, d$v, split, record, formulas[[model]])
      got <- clusters[clusters$model == model &
                        clusters$fold == split$label, , drop = FALSE]
      got <- got[match(hand$cluster, got$cluster), , drop = FALSE]
      expect_false(anyNA(got$mse))
      expect_lt(fbkv_rel(got$mse, hand$mse), 1.5e-8)
      expect_lt(fbkv_rel(got$nlpd, hand$nlpd), 1.5e-8)
    }
  }
}

# A base that carries the parsed model, a known variance and the data, and
# was never fitted.
fbkv_stub_base <- function(d, class = "flexybayes_inla") {
  fb <- fb_from_asreml(ybar ~ variety, random = ~trial, data = d)
  fb$family <- "gaussian"
  fb$link <- "identity"
  attached <- .fb_known_variance_attach(fb, d, d$v, aggregate = "auto")
  structure(
    list(
      fb = attached$fb,
      data = attached$data,
      extras = list(call_info = list(
        fixed = ybar ~ variety, random = ~trial, family = "gaussian"
      ))
    ),
    class = c(class, "flexybayes", "list")
  )
}


# --- (i) the argument ------------------------------------------------------ #

test_that("known_variance that is not one finite, non-negative variance per row refuses by name", {
  d <- fbkv_stage2()
  fit_with <- function(value) {
    flexybayes(ybar ~ variety, random = ~trial, data = d,
               known_variance = value, backend = "inla", return_code = TRUE)
  }
  bad <- list(
    as.character(d$v),
    d$v[-1L],
    replace(d$v, 3L, NA_real_),
    replace(d$v, 3L, Inf),
    replace(d$v, 3L, -0.1),
    matrix(d$v, ncol = 1L)
  )
  for (value in bad) {
    expect_error(fit_with(value),
                 class = "flexybayes_refusal_known_variance_invalid")
  }
})

test_that("a model that cannot carry a known variance refuses by name", {
  d <- fbkv_stage2()
  d$count <- .fb_sel_with_seed(3L, stats::rpois(nrow(d), 3))
  d$hit <- as.integer(d$ybar > 0)
  d$pos <- exp(d$ybar)
  refused <- "flexybayes_refusal_known_variance_not_supported"
  expect_error(
    flexybayes(count ~ variety, random = ~trial, data = d, family = "poisson",
               known_variance = d$v, backend = "inla", return_code = TRUE),
    class = refused
  )
  expect_error(
    flexybayes(hit ~ variety, random = ~trial, data = d, family = "binomial",
               known_variance = d$v, backend = "brms", return_code = TRUE),
    class = refused
  )
  expect_error(
    flexybayes(pos ~ variety, random = ~trial, data = d, link = "log",
               known_variance = d$v, backend = "brms", return_code = TRUE),
    class = refused
  )
  expect_error(
    flexybayes(ybar ~ variety, random = ~trial, data = d, aggregate = TRUE,
               known_variance = d$v, backend = "inla", return_code = TRUE),
    class = refused
  )
  d$flexybayes_known_variance <- 1
  expect_error(
    flexybayes(ybar ~ variety, random = ~trial, data = d,
               known_variance = d$v, backend = "inla", return_code = TRUE),
    class = "flexybayes_refusal_nested_index_name_collision"
  )
})


# --- (ii) the lowering on each engine ------------------------------------ #

test_that("INLA fits a unit-precision row effect scaled by sqrt(v), and zeros change nothing", {
  skip_if_not_installed("INLA")
  d <- fbkv_stage2()
  code <- function(value) {
    flexybayes(ybar ~ variety, random = ~trial, data = d,
               known_variance = value, backend = "inla", return_code = TRUE)
  }
  with_v <- code(d$v)
  form <- paste(deparse(with_v$formula), collapse = "")
  expect_match(
    gsub("\\s+", " ", form),
    paste0("f(flexybayes_known_index, flexybayes_known_scale, ",
           "model = \"iid\", hyper = list(prec = list(initial = 0, ",
           "fixed = TRUE)))"),
    fixed = TRUE
  )
  zero <- code(rep(0, nrow(d)))
  none <- code(NULL)
  expect_identical(deparse(zero$formula), deparse(none$formula))
  expect_identical(zero[names(zero) != "formula"],
                   none[names(none) != "formula"])
})

test_that("brms fits se(sqrt(v), sigma = TRUE), and zeros change nothing", {
  skip_if_not_installed("brms")
  d <- fbkv_stage2()
  code <- function(value) {
    flexybayes(ybar ~ variety, random = ~trial, data = d,
               known_variance = value, backend = "brms", return_code = TRUE)
  }
  expect_match(code(d$v), "normal_lpdf(Y | mu, sqrt(square(sigma) + se2))",
               fixed = TRUE)
  expect_identical(code(rep(0, nrow(d))), code(NULL))
})


# --- (iii) recovery, agreement and the REML oracle ------------------------ #

test_that("INLA and brms recover the variances and agree with each other and with REML", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  skip_if_not_installed("posterior")
  fbkv_quiet()
  d <- fbkv_stage2()
  n <- nrow(d)
  fit_inla <- flexybayes(ybar ~ variety, random = ~trial, data = d,
                         known_variance = d$v, backend = "inla",
                         verbose = FALSE)
  fit_brms <- suppressWarnings(flexybayes(
    ybar ~ variety, random = ~trial, data = d, known_variance = d$v,
    backend = "brms", chains = 2, n_samples = 1000, warmup = 500, seed = 1L,
    verbose = FALSE, mcmc_verbose = FALSE
  ))

  reml <- fbkv_reml(d$ybar, stats::model.matrix(~variety, d),
                    stats::model.matrix(~ 0 + trial, d), d$v)
  truth <- c(sigma = 0.25, sd_trial = 1)
  vi <- fbkv_variances(fit_inla)
  vb <- fbkv_variances(fit_brms)
  for (k in c("sigma", "sd_trial")) {
    for (tab in list(vi, vb)) {
      row <- tab[tab$component == k, ]
      expect_true(reml[[k]] > row$low && reml[[k]] < row$high, info = k)
      expect_true(truth[[k]] > row$low && truth[[k]] < row$high, info = k)
    }
  }

  # The engines, within the measured Monte Carlo spread of the brms means.
  draws <- posterior::as_draws_matrix(fit_brms$brms)
  within <- function(e_inla, e_brms, mcse, z) {
    all(abs(e_inla - e_brms) <= z * sqrt(mcse^2 + (1e-4 * e_inla)^2))
  }
  sigma2 <- as.numeric(draws[, "sigma"])^2
  trial2 <- as.numeric(draws[, "sd_trial__Intercept"])^2
  expect_true(within(vi$mean[vi$component == "sigma"], mean(sigma2),
                     posterior::mcse_mean(sigma2), 4))
  expect_true(within(vi$mean[vi$component == "sd_trial"], mean(trial2),
                     posterior::mcse_mean(trial2), 4))
  b_inla <- coef(fit_inla)
  b_names <- paste0("b_", sub("(Intercept)", "Intercept", names(b_inla),
                              fixed = TRUE))
  b_mcse <- apply(draws[, b_names, drop = FALSE], 2L, function(x) {
    posterior::mcse_mean(as.numeric(x))
  })
  expect_true(within(unname(b_inla), unname(coef(fit_brms)[names(b_inla)]),
                     unname(b_mcse), 4))

  # Fitted values leave the known part out on both engines. INLA's linear
  # predictor carries the row effect, up to 2.4 away on this fixture.
  epred <- brms::posterior_epred(fit_brms$brms)
  f_mcse <- apply(epred, 2L, function(x) posterior::mcse_mean(as.numeric(x)))
  expect_true(within(fitted(fit_inla), colMeans(epred), f_mcse, 5))
  expect_gt(
    max(abs(fit_inla$inla$summary.fitted.values$mean[seq_len(n)] -
              fitted(fit_inla))),
    1
  )
  expect_identical(names(ranef(fit_inla)), "trial")

  # brms predicts the mean of new rows without their known variance.
  nd <- d[1:5, c("variety", "trial")]
  expect_equal(predict(fit_brms, newdata = nd), predict(fit_brms)[1:5],
               tolerance = 1e-12)
})

test_that("a vector of zeros reproduces the fit without it on both engines", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  fbkv_quiet()
  d <- fbkv_stage2()
  zeros <- rep(0, nrow(d))
  fit <- function(backend, value) {
    suppressWarnings(flexybayes(
      ybar ~ variety, random = ~trial, data = d, known_variance = value,
      backend = backend, chains = 1, n_samples = 300, warmup = 300,
      seed = 7L, verbose = FALSE, mcmc_verbose = FALSE
    ))
  }
  i0 <- fit("inla", zeros)
  i1 <- fit("inla", NULL)
  expect_lt(fbkv_rel(fbkv_variances(i0)$mean, fbkv_variances(i1)$mean), 1e-4)
  expect_lt(fbkv_rel(coef(i0), coef(i1)), 1e-4)
  expect_identical(i0$extras$call_info$known_variance, zeros)
  b0 <- fit("brms", zeros)
  b1 <- fit("brms", NULL)
  expect_identical(as.matrix(b0$brms), as.matrix(b1$brms))
})


test_that("weights and a known variance compose as sigma^2 / w + v on both engines", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  skip_if_not_installed("posterior")
  fbkv_quiet()
  d <- fbkv_stage2()
  w <- .fb_sel_with_seed(11L, stats::runif(nrow(d), 0.5, 2))
  fit_inla <- flexybayes(ybar ~ variety, random = ~trial, data = d,
                         weights = w, known_variance = d$v, backend = "inla",
                         verbose = FALSE)
  fit_brms <- suppressWarnings(flexybayes(
    ybar ~ variety, random = ~trial, data = d, weights = w,
    known_variance = d$v, backend = "brms", chains = 2, n_samples = 1000,
    warmup = 500, seed = 1L, verbose = FALSE, mcmc_verbose = FALSE
  ))
  reml <- fbkv_reml(d$ybar, stats::model.matrix(~variety, d),
                    stats::model.matrix(~ 0 + trial, d), d$v, w)
  for (tab in list(fbkv_variances(fit_inla), fbkv_variances(fit_brms))) {
    for (k in c("sigma", "sd_trial")) {
      row <- tab[tab$component == k, ]
      expect_true(reml[[k]] > row$low && reml[[k]] < row$high, info = k)
    }
  }
  draws <- posterior::as_draws_matrix(fit_brms$brms)
  b_inla <- coef(fit_inla)
  b_names <- paste0("b_", sub("(Intercept)", "Intercept", names(b_inla),
                              fixed = TRUE))
  b_mcse <- apply(draws[, b_names, drop = FALSE], 2L, function(x) {
    posterior::mcse_mean(as.numeric(x))
  })
  b_brms <- coef(fit_brms)[names(b_inla)]
  expect_true(all(abs(b_inla - b_brms) <=
                    4 * sqrt(b_mcse^2 + (1e-4 * b_inla)^2)))
})


test_that("backend = \"auto\", review_code then proceed(), and a separate process carry it", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  fbkv_quiet()
  d <- fbkv_stage2()
  inla <- flexybayes(ybar ~ variety, random = ~trial, data = d,
                     known_variance = d$v, backend = "inla", verbose = FALSE)
  auto <- flexybayes(ybar ~ variety, random = ~trial, data = d,
                     known_variance = d$v, verbose = FALSE)
  expect_s3_class(auto, "flexybayes_inla")
  expect_lt(fbkv_rel(fbkv_variances(auto)$mean, fbkv_variances(inla)$mean),
            1e-4)

  args <- list(ybar ~ variety, random = ~trial, data = d,
               known_variance = d$v, backend = "brms", chains = 1,
               n_samples = 300, warmup = 300, seed = 3L)
  review <- do.call(flexybayes, c(args, review_code = TRUE))
  proceeded <- suppressWarnings(proceed(review))
  direct <- suppressWarnings(do.call(
    flexybayes, c(args, verbose = FALSE, mcmc_verbose = FALSE)
  ))
  expect_identical(as.matrix(proceeded$brms), as.matrix(direct$brms))

  skip_if_not_installed("callr")
  isolated <- withr::with_options(
    list(flexyBayes.isolation = "process"),
    flexybayes(ybar ~ variety, random = ~trial, data = d,
               known_variance = d$v, backend = "inla", verbose = FALSE)
  )
  expect_identical(isolated$extras$isolation$mode, "process")
  expect_lt(
    fbkv_rel(fbkv_variances(isolated)$mean, fbkv_variances(inla)$mean),
    1e-4
  )
})


# --- (iv) rows through the missing-data layer ----------------------------- #

test_that("a known variance follows its rows through omit and augment", {
  skip_if_not_installed("INLA")
  fbkv_quiet()
  d <- fbkv_stage2()
  miss <- c(4L, 50L, 111L)
  d$ybar[miss] <- NA_real_
  omit <- flexybayes(ybar ~ variety, random = ~trial, data = d,
                     known_variance = d$v, na_action = "omit",
                     backend = "inla", verbose = FALSE)
  expect_identical(omit$data$flexybayes_known_variance, d$v[-miss])
  augment <- flexybayes(ybar ~ variety, random = ~trial, data = d,
                        known_variance = d$v, backend = "inla",
                        verbose = FALSE)
  expect_identical(augment$data$flexybayes_known_variance, d$v)
  skip_if_not_installed("brms")
  code <- flexybayes(ybar ~ variety, random = ~trial, data = d,
                     known_variance = d$v, backend = "brms",
                     return_code = TRUE)
  expect_match(code, "normal_lpdf(Yl | mu, sqrt(square(sigma) + se2))",
               fixed = TRUE)

})

test_that("a plot the design grid adds carries no known variance", {
  # A plot absent from a field's grid is added with no response. The data
  # the engine would receive are read where flexybayes() hands them over.
  g <- .fb_sel_with_seed(5L, {
    grid <- expand.grid(col = 1:5, row = 1:6)
    out <- data.frame(row = factor(grid$row), col = factor(grid$col))
    out$y <- stats::rnorm(nrow(out))
    out$v <- stats::runif(nrow(out), 0.1, 0.5)
    out[-7L, ]
  })
  local_mocked_bindings(.dispatch_backend = function(fb, data, ...) data)
  handed <- flexybayes(y ~ 1, random = ~ ar1(row):ar1(col), data = g,
                       known_variance = g$v, backend = "inla",
                       verbose = FALSE)
  added <- is.na(handed$y)
  expect_identical(sum(added), 1L)
  expect_identical(handed$flexybayes_known_variance[added], 0)
  key <- paste(handed$row, handed$col)[!added]
  expect_identical(
    handed$flexybayes_known_variance[!added],
    g$v[match(key, paste(g$row, g$col))]
  )
})

test_that("update() re-issues the known variance with the rows it refits", {
  skip_if_not_installed("INLA")
  fbkv_quiet()
  d <- fbkv_stage2()
  fit <- flexybayes(ybar ~ variety, random = ~trial, data = d,
                    known_variance = d$v, backend = "inla", verbose = FALSE)
  again <- update(fit, verbose = FALSE)
  expect_identical(again$extras$call_info$known_variance, d$v)
  expect_lt(fbkv_rel(fbkv_variances(again)$mean, fbkv_variances(fit)$mean),
            1e-4)
})


# --- (v) fb_select(covariance = "base_model") ----------------------------- #

test_that("fb_select() whitens with sigma^2 + v, dense and block by block, and scores equal GLS written out", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbkv_quiet()
  d <- fbkv_select_data()
  base <- flexybayes(ybar ~ variety, random = ~trial, data = d,
                     known_variance = d$v, backend = "inla", verbose = FALSE)
  run <- function() {
    fb_select(base, ~ z1 + z2 + z3 + z4, workflow = "unpenalized_benchmark",
              covariance = "base_model",
              resamples = fb_resamples(unit = ~trial, v = 3L), seed = 1L)
  }
  dense <- run()
  expect_identical(dense$certificate$covariance$whitening$mode, "dense")
  expect_true(dense$certificate$covariance$known_variance)
  expect_output(print(dense$certificate), "known variance added")
  fbkv_expect_gls_scores(dense, d)

  withr::local_options(flexyBayes.selection_dense_gb = 1e-6)
  by_block <- run()
  expect_identical(by_block$certificate$covariance$whitening$mode, "by_block")
  fbkv_expect_gls_scores(by_block, d)

  # The penalised optimiser runs on the same folds and covariance.
  skip_if_not_installed("glmnet")
  opt <- fb_select(base, ~ z1 + z2 + z3 + z4,
                   penalty = fbkv_penalty(),
                   workflow = "penalized_optimizer", covariance = "base_model",
                   resamples = fb_resamples(unit = ~trial, v = 3L), seed = 1L)
  expect_true(opt$certificate$covariance$known_variance)
  expect_identical(opt$certificate$covariance$whitening$mode, "by_block")
})

test_that("each fold refit is handed the known variance of the rows it fits", {
  d <- fbkv_select_data()
  split <- list(label = "1", train_rows = 1:60, assess_rows = 61:90)
  for (provider in c("inla", "brms")) {
    base <- fbkv_stub_base(d, paste0("flexybayes_", provider))
    seen <- NULL
    local_mocked_bindings(flexybayes = function(...) {
      seen <<- list(...)
      stop("captured")
    })
    expect_error(
      .fb_sel_cov_fold_refit(base, split, provider = provider, slots = list()),
      "captured"
    )
    rows <- if (provider == "inla") seq_len(nrow(d)) else split$train_rows
    expect_identical(seen$known_variance, d$v[rows])
    expect_false(any(.fb_known_variance_columns() %in% names(seen$data)))
    if (provider == "inla") {
      expect_true(all(is.na(seen$data$ybar[split$assess_rows])))
    }
  }
})

test_that("a brms base carries its known variance through fb_select() end to end", {
  skip_on_cran()
  skip_if_not_installed("brms")
  fbkv_quiet()
  d <- fbkv_select_data()
  base <- suppressWarnings(flexybayes(
    ybar ~ variety, random = ~trial, data = d, known_variance = d$v,
    backend = "brms", chains = 1, n_samples = 400, warmup = 400, seed = 1L,
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  sel <- suppressWarnings(fb_select(
    base, ~ z1 + z2 + z3 + z4, workflow = "unpenalized_benchmark",
    covariance = "base_model",
    resamples = fb_resamples(unit = ~trial, v = 3L), seed = 1L,
    engine = list(chains = 1, n_samples = 400, warmup = 400)
  ))
  expect_identical(sel$certificate$covariance$engine, "brms")
  expect_true(sel$certificate$covariance$known_variance)
  fbkv_expect_gls_scores(sel, d)
})

test_that("routes that would leave a known variance out refuse by name", {
  skip_if_not_installed("INLA")
  fbkv_quiet()
  d <- fbkv_select_data()
  base <- flexybayes(ybar ~ variety, random = ~trial, data = d,
                     known_variance = d$v, backend = "inla", verbose = FALSE)
  refused <- "flexybayes_refusal_known_variance_not_supported"
  resamples <- fb_resamples(unit = ~trial, v = 3L)
  # Each refuses before any fold is refitted or any fold covariance built.
  local({
    local_mocked_bindings(
      .fb_sel_bayes_refit = function(...) stop("a fold was refitted"),
      .fb_sel_fold_substrate = function(...) stop("a fold was compiled")
    )
    expect_error(fb_select(base, ~ z1 + z2, resamples = resamples, seed = 1L),
                 class = refused)
    expect_error(fb_select(base, ~ z1 + z2, workflow = "unpenalized_benchmark",
                           resamples = resamples, seed = 1L),
                 class = refused)
    expect_error(fb_select(base, ~ z1 + z2, workflow = "unpenalized_benchmark",
                           resamples = resamples, seed = 1L, plan = TRUE),
                 class = refused)
    expect_error(
      fb_select(base, ~ z1 + z2,
                penalty = fbkv_penalty(),
                workflow = "penalized_optimizer", covariance = "base_model",
                refit = list(method = "bayes"), resamples = resamples,
                seed = 1L),
      class = refused
    )
  })

  # The loci beneath, each refusing where it would otherwise build a model
  # without the known variance.
  stub <- fbkv_stub_base(d)
  split <- list(label = "1", train_rows = 1:60, assess_rows = 61:90)
  expect_error(
    .fb_sel_bayes_refit(stub, stub$data, ybar ~ variety, "inla"),
    class = refused
  )
  expect_error(
    .fb_sel_fold_substrate(stub$data, split, response = "ybar",
                           protected = ~variety, random = "trial"),
    class = refused
  )
  expect_error(.fb_sel_variance_components(stub), class = refused)
  expect_error(
    .fb_sel_reference_fit(stub, ~ z1 + z2, prior = NULL, seed = 1L),
    class = refused
  )
})

test_that("the grouped ablation of a known-variance base scores the protected model with sigma^2 + v", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("glmnet")
  fbkv_quiet()
  d <- fbkv_select_data()
  base <- flexybayes(ybar ~ variety, random = ~trial, data = d,
                     known_variance = d$v, backend = "inla", verbose = FALSE)
  sel <- fb_select(base, ~ z1 + z2 + z3 + z4, penalty = fbkv_penalty(),
                   workflow = "penalized_optimizer", covariance = "base_model",
                   resamples = fb_resamples(unit = ~trial, v = 3L), seed = 1L,
                   decision = fb_decision(method = "full", metric = "rmse"),
                   ablation = TRUE)
  expect_true(sel$certificate$covariance$known_variance)
  # Dropping the only declared block leaves the protected model, so each
  # trial's ablated loss is the root mean squared error of the GLS fit with
  # V = diag(sigma^2 + v) + sigma_t^2 Z Z', written out above.
  cov <- sel$certificate$covariance
  labels <- vapply(cov$folds, function(r) r$fold, character(1L))
  hand <- NULL
  for (split in sel$folds) {
    record <- cov$folds[[match(split$label, labels)]]
    hand <- rbind(hand, fbkv_gls_clusters(d, d$v, split, record, ~variety))
  }
  hand_rmse <- stats::setNames(sqrt(hand$mse), hand$cluster)
  units <- sel$ablation$units
  expect_setequal(units$unit, names(hand_rmse))
  expect_lt(fbkv_rel(units$loss_ablated, hand_rmse[units$unit]), 1.5e-8)
})
