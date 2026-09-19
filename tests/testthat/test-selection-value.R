# =============================================================================
# fb_select(refit = ): the Bayesian refit of each fold's frozen support, its
# held-out log predictive density per assessment row and per held-out unit,
# and the share of between-trial variance a declared block explains.
#
# Each block holds the valuation to an oracle written in this file without
# package code.
#
#   The closed-form predictive. With both variance components pinned at one
#   by lognormal(0, 0.001) priors (`prior_vc_sd = 0.001`), which have no hard
#   edge for INLA's optimiser to meet, and the engine's near-flat default
#   prior on the fixed effects, the posterior predictive of an
#   assessment row is Gaussian and has a closed form: the mixed-model
#   equations over the fold's training rows give the joint posterior of the
#   fixed effects and the site effects; a new cluster adds the site variance
#   to the fixed part's posterior variance, and a known cluster conditions on
#   the site's posterior. The stored density is a Monte Carlo average over
#   posterior draws. The seeded sampler repeats itself on this build (two
#   calls with one seed agreed to 7e-5), so three runs with one seed would
#   measure nothing of that Monte Carlo error. The three runs are identical
#   except for the seed that drives the posterior sampler: the folds are
#   explicit, so the allocation and every frozen support stay the same. The
#   tolerance is five times the largest per-row standard deviation across
#   them, and is never widened.
#
#   Paired differencing. Two valuations on one value contract hold out the
#   same rows in the same folds, so their unit tables pair unit by unit and
#   the paired difference and its standard error are recomputed here from
#   the two stored tables.
#
#   The block share. Three groups of forty trials are built so that in
#   every union of two groups the signal covariate and the trial effect have
#   mean zero, equal sample variance and zero cross-product: the share of
#   between-trial variance the signal explains is 0.5 in every training fold
#   by construction. A noise covariate orthogonal to both is near zero.
#
#   Known issue 19. A fold whose training rows cannot tell a numeric
#   candidate from one column of a factor is refused by the formula refits;
#   the refit of a frozen support writes each support column as its own
#   design column, so that fold is refitted on the columns it kept.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

fbv_quiet <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE,
    .local_envir = envir
  )
}

# Twelve sites of eight rows: one protected covariate, four candidates of
# which two carry signal, a site effect and a residual both of SD one.
fbv_data <- function(seed = 101L) {
  .fb_sel_with_seed(seed, {
    n_site <- 12L
    per <- 8L
    d <- data.frame(
      site = factor(rep(sprintf("S%02d", seq_len(n_site)), each = per))
    )
    n <- nrow(d)
    d$x <- stats::rnorm(n)
    for (j in 1:4) {
      d[[paste0("z", j)]] <- stats::rnorm(n)
    }
    d$y <- 2 + 0.5 * d$x + 0.9 * d$z1 - 0.7 * d$z2 +
      stats::rnorm(n_site)[as.integer(d$site)] + stats::rnorm(n)
    # Three folds of rows, fixed here, so a site has rows on both sides of
    # every fold and the seed moves only the posterior sampler.
    d$fold <- sample(rep(1:3, length.out = n))
    d
  })
}

fbv_penalty <- function() {
  fb_penalty(term(everything()) ~ laplace(scale = 0.5))
}

fbv_select <- function(base,
                       prediction,
                       candidates = ~ z1 + z2 + z3 + z4,
                       refit = list(method = "bayes"),
                       seed = 11L) {
  fb_select(
    base,
    candidates,
    penalty = fbv_penalty(),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(fold_id = base$data$fold),
    prediction = prediction,
    seed = seed,
    refit = refit
  )
}

# INLA fits and valuations are the slow part of the file, so each is built
# once.
fbv_cache <- new.env(parent = emptyenv())

fbv_cached <- function(key, expr) {
  if (!exists(key, envir = fbv_cache, inherits = FALSE)) {
    assign(key, expr, envir = fbv_cache)
  }
  get(key, envir = fbv_cache, inherits = FALSE)
}

# Both standard deviations pinned at one: lognormal(0, 0.001) on each, which
# puts 99.7% of the prior within 0.3% of one.
fbv_pinned_base <- function() {
  fbv_cached("pinned_base", flexybayes(
    y ~ x,
    random = ~site,
    data = fbv_data(),
    backend = "inla",
    prior_vc_sd = 0.001,
    verbose = FALSE
  ))
}

# Three runs identical except for the seed of the posterior sampler.
fbv_runs <- function(prediction) {
  fbv_cached(paste0("runs_", prediction), lapply(11:13, function(seed) {
    fbv_select(fbv_pinned_base(), prediction, seed = seed)
  }))
}

# The closed-form posterior predictive of one fold's assessment rows under a
# Gaussian random-intercept model with known variances and a flat prior on
# the fixed effects. `columns` are the numeric columns of the fixed part
# besides the intercept.
fbv_closed_form <- function(d, split, columns, sigma_e, sigma_u, target) {
  X <- cbind(1, as.matrix(d[, columns, drop = FALSE]))
  train <- split$train_rows
  site <- as.character(d$site)
  levels_train <- sort(unique(site[train]))
  Xt <- X[train, , drop = FALSE]
  Zt <- 1 * outer(site[train], levels_train, "==")
  p <- ncol(X)
  q <- length(levels_train)
  Q <- rbind(
    cbind(crossprod(Xt), crossprod(Xt, Zt)),
    cbind(crossprod(Zt, Xt), crossprod(Zt))
  ) / sigma_e^2
  Q[p + seq_len(q), p + seq_len(q)] <-
    Q[p + seq_len(q), p + seq_len(q)] + diag(1 / sigma_u^2, q)
  S <- solve(Q)
  m <- S %*% c(crossprod(Xt, d$y[train]), crossprod(Zt, d$y[train])) /
    sigma_e^2
  out <- lapply(split$assess_rows, function(i) {
    w <- c(X[i, ], rep(0, q))
    unseen <- sigma_u^2
    if (identical(target, "known_cluster") && site[[i]] %in% levels_train) {
      w[p + match(site[[i]], levels_train)] <- 1
      unseen <- 0
    }
    mean <- sum(w * m)
    variance <- drop(crossprod(w, S %*% w)) + unseen + sigma_e^2
    data.frame(
      row = i,
      mean = mean,
      variance = variance,
      lpd = stats::dnorm(d$y[[i]], mean, sqrt(variance), log = TRUE)
    )
  })
  do.call(rbind, out)
}

fbv_expected <- function(sel, d, target) {
  supports <- sel$value$supports
  out <- do.call(rbind, lapply(sel$folds, function(split) {
    kept <- supports$column[supports$fold == split$label]
    fbv_closed_form(d, split, c("x", kept), 1, 1, target)
  }))
  out <- out[order(out$row), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Three groups of forty trials of five plots. Within each group the signal
# `s1`, the trial effect `u` and the noise `n1` are the columns of an
# orthonormal basis orthogonal to the constant, scaled to sample variance
# one; so in any union of groups they keep mean zero, equal sample variance
# and zero cross-products.
fbv_share_data <- function(seed = 202L) {
  .fb_sel_with_seed(seed, {
    n_group <- 3L
    per_group <- 40L
    plots <- 5L
    n_trial <- n_group * per_group
    group <- rep(seq_len(n_group), each = per_group)
    s1 <- numeric(n_trial)
    u <- numeric(n_trial)
    n1 <- numeric(n_trial)
    for (g in seq_len(n_group)) {
      idx <- which(group == g)
      raw <- matrix(stats::rnorm(per_group * 3L), per_group)
      basis <- qr.Q(qr(cbind(1, raw)))
      basis <- basis[, 2:4] * sqrt(per_group - 1)
      s1[idx] <- basis[, 1L]
      u[idx] <- basis[, 2L]
      n1[idx] <- basis[, 3L]
    }
    d <- data.frame(
      trial = factor(rep(sprintf("T%02d", seq_len(n_trial)), each = plots)),
      fold = rep(group, each = plots)
    )
    tr <- as.integer(d$trial)
    n <- nrow(d)
    d$s1 <- s1[tr]
    d$n1 <- n1[tr]
    d$u <- u[tr]
    d$x <- stats::rnorm(n)
    d$w1 <- stats::rnorm(n)
    d$y <- 1 + 0.3 * d$x + d$s1 + d$u + stats::rnorm(n, sd = 0.5)
    d
  })
}

# The block-share valuations: the signal block `s1` is a candidate, the noise
# block `n1` is protected and so in every fold's model.
fbv_share_run <- function(which = c("signal", "noise")) {
  which <- match.arg(which)
  d <- fbv_share_data()
  base <- fbv_cached("share_base", flexybayes(
    y ~ x + n1, random = ~trial, data = d, backend = "inla", verbose = FALSE
  ))
  fbv_cached(paste0("share_", which), fb_select(
    base, ~ s1 + w1,
    penalty = fbv_penalty(),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~trial, fold_id = d$fold),
    decision = fb_decision(rule = "min"),
    seed = 7L,
    refit = list(
      method = "bayes",
      block = if (identical(which, "signal")) ~s1 else ~n1
    )
  ))
}

# Six locations of two trials of six plots, folds of two locations. `k`
# equals the indicator of soil C on every row outside location L1, so fold 1,
# which holds L1 out, cannot tell `k` from the soil C column.
fbv_drop_data <- function(seed = 303L) {
  .fb_sel_with_seed(seed, {
    locations <- sprintf("L%d", 1:6)
    soil_of <- c(L1 = "A", L2 = "A", L3 = "B", L4 = "B", L5 = "C", L6 = "C")
    fold_of <- c(L1 = 1L, L3 = 1L, L2 = 2L, L5 = 2L, L4 = 3L, L6 = 3L)
    trials_per <- 2L
    plots <- 6L
    n_trial <- length(locations) * trials_per
    d <- data.frame(
      trial = factor(rep(sprintf("T%02d", seq_len(n_trial)), each = plots)),
      location = factor(rep(rep(locations, each = trials_per), each = plots))
    )
    n <- nrow(d)
    loc <- as.character(d$location)
    d$soil <- factor(unname(soil_of[loc]), levels = c("A", "B", "C"))
    d$fold <- unname(fold_of[loc])
    d$k <- as.numeric(d$soil == "C") + as.numeric(loc == "L1")
    d$blk <- factor(rep(c("b1", "b2"), length.out = n))
    d$x <- stats::rnorm(n)
    tr <- as.integer(d$trial)
    d$y <- 1 + 0.4 * d$x + 1.5 * (d$soil == "B") + 0.2 * (d$blk == "b2") +
      stats::rnorm(n_trial, sd = 0.4)[tr] + stats::rnorm(n, sd = 0.5)
    d
  })
}

# A base carrying the parsed model, the data and the call record, with no
# engine behind it: the refusals below fire before any fit.
fbv_stub_base <- function(d) {
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


# --- (i) the request is validated before anything is refitted ------------ #

test_that("a refit on a route with no frozen support refuses by name", {
  fbv_quiet()
  d <- fbv_data()
  stub <- fbv_stub_base(d)
  # The Bayesian route with no penalty holds every candidate: no support.
  err <- tryCatch(
    fb_select(
      stub, ~ z1 + z2,
      resamples = fb_resamples(v = 3L), seed = 1L,
      refit = list(method = "bayes")
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_selection_refit_support_undefined")
  # A ridge path shrinks every coefficient and zeroes none, so it has no
  # support to freeze.
  err <- tryCatch(
    fb_select(
      stub, ~ z1 + z2,
      penalty = fb_penalty(term(everything()) ~ ridge(scale = 1)),
      workflow = "penalized_optimizer",
      resamples = fb_resamples(v = 3L), seed = 1L,
      refit = list(method = "bayes")
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_selection_refit_support_undefined")
  expect_match(conditionMessage(err), "ridge", fixed = TRUE)
})

test_that("a malformed refit request stops with the part it is missing", {
  fbv_quiet()
  d <- fbv_data()
  stub <- fbv_stub_base(d)
  call_with <- function(refit) {
    fb_select(
      stub, ~ z1 + z2,
      penalty = fbv_penalty(), workflow = "penalized_optimizer",
      resamples = fb_resamples(v = 3L), seed = 1L, refit = refit
    )
  }
  bayes <- function(...) list(method = "bayes", ...)
  expect_error(call_with("bayes"), "named list", fixed = TRUE)
  expect_error(call_with(list(method = "frequentist")), "method", fixed = TRUE)
  expect_error(call_with(bayes(draws = 10L)), "draws", fixed = TRUE)
  expect_error(call_with(bayes(n_draws = 0L)), "n_draws", fixed = TRUE)
  expect_error(call_with(bayes(group = "nope")), "site", fixed = TRUE)
  expect_error(call_with(bayes(block = ~nope)), "nope", fixed = TRUE)
  expect_error(call_with(bayes(block = "z1")), "one-sided", fixed = TRUE)
})


test_that("a fold refit carries the base's scalar priors", {
  d <- fbv_data()
  capture <- function(supplied) {
    base <- fbv_stub_base(d)
    base$fb$prior_scalars <- list(
      fixed_sd = 10,
      vc_sd = 0.001,
      supplied = c(fixed_sd = supplied, vc_sd = supplied)
    )
    seen <- NULL
    local_mocked_bindings(
      flexybayes = function(...) {
        seen <<- list(...)
        structure(
          list(extras = list()),
          class = c("flexybayes_inla", "flexybayes", "list")
        )
      }
    )
    .fb_sel_bayes_refit(
      base, d, .fb_sel_bayes_fixed("y", "x", "z1"), "inla",
      na_action = "augment"
    )
    seen
  }
  stated <- capture(TRUE)
  expect_identical(stated$prior_vc_sd, 0.001)
  expect_identical(stated$prior_fixed_sd, 10)
  # A base fitted without them keeps the package default.
  unstated <- capture(FALSE)
  expect_false("prior_vc_sd" %in% names(unstated))
  expect_false("prior_fixed_sd" %in% names(unstated))
})


test_that("scoring a block of rows at a time gives the scores of all rows at once", {
  draws <- 50L
  rows <- 23L
  mean <- .fb_sel_with_seed(5L, matrix(stats::rnorm(draws * rows), draws, rows))
  variance <- .fb_sel_with_seed(6L, stats::rexp(draws) + 0.5)
  y <- .fb_sel_with_seed(7L, stats::rnorm(rows))
  whole <- .fb_sel_value_score(y, mean, variance)
  blocks <- .fb_sel_value_score_blocks(
    y,
    function(columns) {
      list(mean = mean[, columns, drop = FALSE], variance = variance)
    },
    size = 5L
  )
  expect_identical(blocks, whole)
  # The density of each row, written out here: the log of the mean over
  # draws of the Gaussian density.
  by_hand <- vapply(seq_len(rows), function(j) {
    log(mean(stats::dnorm(y[[j]], mean[, j], sqrt(variance))))
  }, numeric(1L))
  expect_equal(whole$lpd, by_hand, tolerance = 1e-12)
})


# --- (ii) the stored density against the closed form --------------------- #

test_that("the stored per-row log predictive density equals the closed-form Gaussian predictive", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbv_quiet()
  d <- fbv_data()
  for (target in c("new_cluster", "known_cluster")) {
    runs <- fbv_runs(target)
    tabs <- lapply(runs, augment, type = "value_rows")
    # One allocation and one frozen support in every run.
    for (r in 2:3) {
      expect_identical(tabs[[r]]$row, tabs[[1L]]$row)
      expect_identical(runs[[r]]$value$supports, runs[[1L]]$value$supports)
    }
    lpd <- vapply(tabs, function(t) t$lpd, numeric(nrow(tabs[[1L]])))
    spread <- apply(lpd, 1L, stats::sd)
    # A zero spread would mean the seed does not reach the sampler, and the
    # tolerance would measure nothing.
    expect_gt(max(spread), 0)
    tolerance <- 5 * max(spread)
    expected <- fbv_expected(runs[[1L]], d, target)
    expect_identical(tabs[[1L]]$row, expected$row)
    expect_lte(max(abs(tabs[[1L]]$lpd - expected$lpd)), tolerance)
    expect_lte(
      max(abs(tabs[[1L]]$predictive_mean - expected$mean)),
      5 * max(apply(
        vapply(tabs, function(t) t$predictive_mean, numeric(nrow(d))),
        1L, stats::sd
      ))
    )
  }
  new_rows <- augment(fbv_runs("new_cluster")[[1L]], type = "value_rows")
  known_rows <- augment(fbv_runs("known_cluster")[[1L]], type = "value_rows")
  # The unseen site's effect is integrated in closed form per draw; a known
  # site is read off the draws.
  expect_true(all(new_rows$method == "analytic_gaussian"))
  expect_true(all(known_rows$method == "posterior_predictive_draws"))
})

test_that("the new-cluster predictive is wider than the known-cluster one by the closed-form amount", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbv_quiet()
  d <- fbv_data()
  new_runs <- fbv_runs("new_cluster")
  known_runs <- fbv_runs("known_cluster")
  widen <- vapply(1:3, function(r) {
    a <- augment(new_runs[[r]], type = "value_rows")
    b <- augment(known_runs[[r]], type = "value_rows")
    a$predictive_sd^2 - b$predictive_sd^2
  }, numeric(nrow(d)))
  spread <- apply(widen, 1L, stats::sd)
  expect_gt(max(spread), 0)
  tolerance <- 5 * max(spread)
  expected <- fbv_expected(new_runs[[1L]], d, "new_cluster")$variance -
    fbv_expected(known_runs[[1L]], d, "known_cluster")$variance
  expect_true(all(expected > 0))
  expect_true(all(widen[, 1L] > 0))
  expect_lte(max(abs(widen[, 1L] - expected)), tolerance)
})


# --- (iii) paired differencing ------------------------------------------- #

test_that("two valuations on one contract difference unit by unit", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbv_quiet()
  a <- fbv_runs("new_cluster")[[1L]]
  b <- fbv_cached(
    "noise_only",
    fbv_select(fbv_pinned_base(), "new_cluster", candidates = ~ z3 + z4)
  )
  ga <- glance(a)
  gb <- glance(b)
  expect_identical(ga$value_contract_hash, gb$value_contract_hash)
  ua <- augment(a, type = "value_units")
  ub <- augment(b, type = "value_units")
  # The pairing is exact: the same units, held out in the same folds, with
  # the same rows.
  expect_identical(ua[, c("repeat_id", "unit", "fold", "n_rows")],
                   ub[, c("repeat_id", "unit", "fold", "n_rows")])
  paired <- augment(a, type = "value_units", against = b)
  expect_identical(paired$elpd, ua$elpd)
  expect_identical(paired$elpd_against, ub$elpd)
  expect_identical(paired$difference, ua$elpd - ub$elpd)
  # The difference and its standard error, from the paired unit differences.
  delta <- ua$elpd - ub$elpd
  g <- glance(a, against = b)
  expect_equal(g$value_difference, sum(delta), tolerance = 1e-12)
  expect_equal(
    g$value_difference_se, sqrt(length(delta)) * stats::sd(delta),
    tolerance = 1e-12
  )
  expect_equal(g$value_difference, ga$value_elpd - gb$value_elpd,
               tolerance = 1e-10)
  # The support that carries signal predicts the held-out rows better.
  expect_gt(g$value_difference, 2 * g$value_difference_se)
  # Two valuations on different contracts do not pair.
  known <- fbv_runs("known_cluster")[[1L]]
  err <- tryCatch(glance(a, against = known), error = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_selection_value_contract_mismatch")
  expect_match(conditionMessage(err), "prediction", fixed = TRUE)
})


# --- (iv) the block share ------------------------------------------------- #

test_that("the block share covers a known share by construction, and a noise block gives near zero", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbv_quiet()
  signal <- fbv_share_run("signal")
  noise <- fbv_share_run("noise")

  # The signal is in every fold's support, so every fold's share has a
  # numerator; the construction makes the true share 0.5 in each.
  supports <- signal$value$supports
  expect_setequal(supports$fold[supports$column == "s1"], c("1", "2", "3"))
  share <- signal$value$block_share
  expect_identical(share$group, "trial")
  expect_identical(nrow(share$by_fold), 3L)
  expect_true(all(share$by_fold$lower <= 0.5 & share$by_fold$upper >= 0.5))
  expect_true(share$pooled$lower <= 0.5 && share$pooled$upper >= 0.5)
  g <- glance(signal)
  expect_identical(g$block_share_median, share$pooled$median)
  expect_identical(g$block_share_lower, share$pooled$lower)
  expect_identical(g$block_share_upper, share$pooled$upper)

  # The noise covariate is protected, so it is in every fold's model. On T
  # training trials with m plots each, its coefficient is centred at zero
  # with variance (sigma_t^2 + sigma_e^2 / m) / ((T - 1) V(n1)), so each
  # draw's share is near (1 + sigma_e^2 / (m sigma_t^2)) chi-squared(1) /
  # (T - 1). "Near zero" is a pooled median below five times that
  # distribution's median: with T = 80, m = 5, sigma_e = 0.5, sigma_t = 1,
  # 5 * 1.05 * 0.455 / 79 = 0.030.
  near_zero <- 5 * (1 + 0.5^2 / 5) * stats::qchisq(0.5, 1) / 79
  expect_lt(noise$value$block_share$pooled$median, near_zero)
  expect_lt(
    noise$value$block_share$pooled$median,
    share$pooled$lower
  )
})


# --- (v) known issue 19 --------------------------------------------------- #

test_that("a fold that leaves part of a factor unidentified is refitted on its support, never dropped", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbv_quiet()
  d <- fbv_drop_data()
  resamples <- fb_resamples(unit = ~location, fold_id = d$fold)
  base <- fbv_cached("drop_base", flexybayes(
    y ~ x + blk, random = ~trial, data = d, backend = "inla", verbose = FALSE
  ))
  # The formula refits refuse the fold, naming it (known issue 19).
  err <- tryCatch(
    fb_select(base, ~ k + soil, resamples = resamples, seed = 1L),
    error = function(e) e
  )
  expect_s3_class(
    err, "flexybayes_refusal_selection_fold_candidate_drop_not_representable"
  )
  expect_match(conditionMessage(err), "fold 1", fixed = TRUE)

  sel <- fbv_cached("drop_run", fb_select(
    base, ~ k + soil,
    penalty = fbv_penalty(),
    workflow = "penalized_optimizer",
    resamples = resamples,
    seed = 1L,
    refit = list(method = "bayes")
  ))
  rows <- augment(sel, type = "value_rows")
  # Every assessment row of every fold is scored: no fold is dropped.
  expect_identical(sort(rows$row), seq_len(nrow(d)))
  expect_setequal(unique(rows$fold), c("1", "2", "3"))
  expect_true(all(is.finite(rows$lpd)))
  # Fold 1 dropped soilC and kept soilB, so its support holds one column of
  # the factor: the shape the formula refits cannot write.
  space <- sel$certificate$design$candidate_space$dropped
  expect_true(any(space$fold == "1" & space$column == "soilC"))
  fold_one <- sel$value$supports$column[sel$value$supports$fold == "1"]
  expect_true("soilB" %in% fold_one)
  expect_false("soilC" %in% fold_one)
  # The protected factor is read by name, so the unseen trial is integrated
  # in closed form on every fold.
  expect_true(all(rows$method == "analytic_gaussian"))
})


# --- (vi) the accessors ---------------------------------------------------- #

test_that("glance, tidy and augment carry the valuation", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbv_quiet()
  d <- fbv_share_data()
  sel <- fbv_share_run("signal")

  rows <- augment(sel, type = "value_rows")
  expect_true(all(c(
    "repeat_id", "fold", "row", "row_key", "unit", "observed",
    "predictive_mean", "predictive_sd", "lpd", "method", "prediction",
    "n_draws", "contract_hash"
  ) %in% names(rows)))
  expect_identical(rows$observed, d$y[rows$row])
  units <- augment(sel, type = "value_units")
  expect_true(all(c(
    "repeat_id", "fold", "unit", "n_rows", "elpd", "contract_hash"
  ) %in% names(units)))
  # Each unit's density is the sum of its rows'.
  by_unit <- tapply(rows$lpd, rows$unit, sum)
  expect_equal(units$elpd, as.numeric(by_unit[units$unit]), tolerance = 1e-12)
  expect_identical(units$n_rows, as.integer(table(rows$unit)[units$unit]))
  expect_identical(nrow(units), nlevels(d$trial))

  g <- glance(sel)
  expect_equal(g$value_elpd, sum(units$elpd), tolerance = 1e-12)
  expect_equal(
    g$value_elpd_se, sqrt(nrow(units)) * stats::sd(units$elpd),
    tolerance = 1e-12
  )
  expect_identical(g$value_n_units, nlevels(d$trial))
  expect_identical(g$value_unit, "trial")
  expect_identical(g$value_contract_hash, unique(rows$contract_hash))
  expect_true(nzchar(g$value_formula))

  # The refit does not reselect: each fold's frozen support is the
  # optimiser's support on that fold at the chosen point. Read on the run
  # whose folds hold different supports, so a support pooled over folds, or
  # chosen again, would differ from it.
  reselect <- fbv_runs("new_cluster")[[1L]]
  arm <- reselect$benchmark$glmnet
  path <- as.data.frame(arm$support$by_lambda)
  chosen <- path[
    path$arm == "benchmark" &
      path$lambda_index == arm$selection$selected_lambda_index, ,
    drop = FALSE
  ]
  held <- lapply(chosen$terms, function(t) {
    parts <- strsplit(t, ", ", fixed = TRUE)[[1L]]
    parts[nzchar(parts)]
  })
  expect_gt(length(unique(lapply(held, sort))), 1L)
  supports <- reselect$value$supports
  for (i in seq_len(nrow(chosen))) {
    expect_setequal(
      supports$column[supports$fold == as.character(chosen$fold[[i]])],
      held[[i]]
    )
  }

  td <- tidy(sel)
  expect_true(all(c("value_folds_in_support", "value_in_block") %in% names(td)))
  s1 <- td[td$term == "s1", , drop = FALSE]
  expect_identical(
    s1$value_folds_in_support,
    sum(sel$value$supports$column == "s1")
  )
  expect_true(s1$value_in_block)
  expect_false(td$value_in_block[td$term == "w1"])

  expect_output(print(sel), "valuation", fixed = TRUE)

  # Without a refit the accessors keep their schema and say nothing ran.
  plain <- fb_select(
    fbv_stub_base(fbv_data()), ~ z1 + z2,
    penalty = fbv_penalty(), workflow = "penalized_optimizer",
    resamples = fb_resamples(v = 3L), seed = 1L
  )
  expect_null(plain$value)
  expect_true(is.na(glance(plain)$value_elpd))
  expect_identical(nrow(augment(plain, type = "value_rows")), 0L)
  expect_true(all(names(rows) %in% names(augment(plain, type = "value_rows"))))
  expect_true(all(is.na(tidy(plain)$value_folds_in_support)))
})


# --- (vii) a brms base ------------------------------------------------------ #

test_that("a brms base is refitted on brms and scored the same way", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  fbv_quiet()
  d <- fbv_data()
  d <- d[as.integer(d$site) <= 6L, , drop = FALSE]
  d$site <- droplevels(d$site)
  base <- flexybayes(
    y ~ x, random = ~site, data = d, backend = "brms",
    chains = 2L, n_samples = 400L, warmup = 400L, seed = 1L,
    verbose = FALSE, mcmc_verbose = FALSE
  )
  sel <- suppressWarnings(fb_select(
    base, ~ z1 + z2,
    penalty = fbv_penalty(),
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~site, v = 2L),
    seed = 1L,
    refit = list(method = "bayes")
  ))
  expect_identical(sel$value$engine, "brms")
  rows <- augment(sel, type = "value_rows")
  expect_identical(sort(rows$row), seq_len(nrow(d)))
  expect_true(all(is.finite(rows$lpd)))
  expect_true(all(rows$method == "analytic_gaussian"))
  expect_length(unique(rows$n_draws), 1L)
  expect_gte(rows$n_draws[[1L]], 400L)
  units <- augment(sel, type = "value_units")
  expect_equal(glance(sel)$value_elpd, sum(units$elpd), tolerance = 1e-12)
  expect_identical(nrow(sel$value$block_share$by_fold), 2L)
})
