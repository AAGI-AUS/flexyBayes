# =============================================================================
# fb_select() fold compilation: every quantity a fold learns reads its
# training rows, and `prediction = "new_cluster"` relabels every grouping
# column the emitters index a random effect by, or refuses the term by name.
#
# Each block holds one defect to an oracle the implementation did not write.
#
#   Default priors on INLA. A fold refit carries its assessment rows with the
#   response set to NA, and the covariate-scaled defaults (the slope bound
#   divided by sd(x), the fixed-effect prior divided by sd(x_j)) were read
#   over all of those rows. The oracle is invariance: an extreme covariate
#   placed only in the assessment rows leaves every default bit-identical,
#   and each default equals its formula evaluated on the training rows.
#
#   Relabelling. The grouping of a random slope, of an `ide()` term and of a
#   multi-way interaction was left alone, and the outer factor of
#   `diag(f):g` was relabelled, which brms refuses as a new factor level.
#   The oracle is the emitters' own output: brms's parse of the emitted
#   formula, and the index of every f() term in the data an INLA fold refit
#   was handed.
#
#   The fold refit of a brms-grammar base. The random terms of
#   `y ~ x + (x || g)` are written in the model formula, and the refit
#   rebuilt that formula from the fixed labels alone, so it fitted no random
#   effect.
#
#   The benchmark route. A transformed response was read as a column name,
#   and a random slope or a three-way interaction reached the covariance
#   builder as nothing, or as its factors fitted as separate intercepts.
#
# Comparisons are `identical()`: every quantity here is computed by the same
# arithmetic on the same values in the same order, so a difference is a
# change of input rows, not rounding.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

# Nine sites of four rows, a slope covariate, two candidates, an offset, three
# crossed factors for an interaction term and an environment factor that every
# site sees. The response is positive so that log(y) is defined.
fbfi_data <- function(seed = 7L) {
  .fb_sel_with_seed(seed, {
    n_site <- 9L
    per <- 4L
    d <- data.frame(
      site = factor(rep(sprintf("S%02d", seq_len(n_site)), each = per))
    )
    n <- nrow(d)
    d$x <- stats::rnorm(n)
    d$z1 <- stats::rnorm(n)
    d$z2 <- stats::rnorm(n)
    d$o <- stats::runif(n, 0, 2)
    d$a <- factor(rep(c("a1", "a2"), length.out = n))
    d$b <- factor(rep(c("b1", "b2", "b3"), each = 2L, length.out = n))
    d$cc <- factor(rep(c("c1", "c2"), each = 3L, length.out = n))
    d$env <- factor(rep(c("E1", "E2"), times = n / 2L))
    d$y <- exp(
      1 + 0.3 * d$x +
        stats::rnorm(n_site, sd = 0.4)[as.integer(d$site)] +
        stats::rnorm(n, sd = 0.3)
    )
    d
  })
}

fbfi_folds <- function(d) {
  .fb_sel_folds(fb_resamples(unit = ~site, v = 3L), d, seed = 1L)
}

# A base carrying the parsed model, the data and the call record the fold
# refit reads, with no engine behind it.
fbfi_stub <- function(d, fixed, random = NULL, class = "flexybayes_inla") {
  fb <- if (is.null(random)) {
    fb_from_brms(fixed, data = d, family = "gaussian")
  } else {
    fb_from_asreml(fixed, random = random, data = d)
  }
  fb$family <- "gaussian"
  fb$link <- fb$link %||% "identity"
  structure(
    list(
      fb = fb,
      data = d,
      extras = list(
        call_info = list(fixed = fixed, random = random, family = "gaussian")
      )
    ),
    class = c(class, "flexybayes", "list")
  )
}

# The columns whose assessment-row values the fold builder changed, other
# than the response it sets to NA.
fbfi_relabelled <- function(d, inputs, response = "y") {
  rows <- inputs$assess_rows
  changed <- vapply(
    names(d),
    function(nm) {
      !identical(
        as.character(d[[nm]][rows]),
        as.character(inputs$engine_data[[nm]][rows])
      )
    },
    logical(1L)
  )
  setdiff(names(d)[changed], response)
}

# The data columns brms groups the emitted formula's random terms by, read
# from brms's own parse of that formula.
fbfi_brms_groups <- function(fb) {
  terms <- brms::brmsterms(stats::as.formula(.fb_to_brms_formula(fb)))
  groups <- unlist(lapply(terms$dpars$mu$re$gcall, function(g) g$groups))
  unique(unlist(strsplit(groups, ":", fixed = TRUE)))
}

# The index column of every f() term in an INLA formula.
fbfi_f_indices <- function(formula) {
  vars <- as.list(attr(stats::terms(formula), "variables"))[-1L]
  out <- character(0)
  for (v in vars) {
    if (is.call(v) && identical(v[[1L]], as.name("f"))) {
      out <- c(out, as.character(v[[2L]]))
    }
  }
  out
}

# INLA fits are the slow part of the file, so each is built once.
fbfi_cache <- new.env(parent = emptyenv())

fbfi_cached <- function(key, expr) {
  if (!exists(key, envir = fbfi_cache, inherits = FALSE)) {
    assign(key, expr, envir = fbfi_cache)
  }
  get(key, envir = fbfi_cache, inherits = FALSE)
}

fbfi_quiet <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE,
    .local_envir = envir
  )
}

fbfi_inla_refit <- function(base, engine_data, candidate = "z1") {
  spec <- .fb_sel_bayes_base_spec(base)
  .fb_sel_bayes_refit(
    base,
    engine_data,
    .fb_sel_bayes_fixed(
      spec$response, spec$protected_labels, candidate, spec$offsets
    ),
    "inla",
    na_action = "augment"
  )
}


# --- (i) default priors on INLA read the training rows ------------------- #

test_that("the INLA default functions read a fold's training rows from its engine data", {
  d <- fbfi_data()
  split <- fbfi_folds(d)[[1L]]
  fb <- fb_from_brms(y ~ x + z1 + (x || site), data = d, family = "gaussian")
  fb$link <- "identity"
  engine <- .fb_sel_bayes_fold_inputs(
    d, split, "y", ~x, ~z1, random = "site"
  )$engine_data
  extreme <- engine
  extreme$x[split$assess_rows[[1L]]] <- 1e3
  extreme$z1[split$assess_rows[[2L]]] <- -1e3
  train <- d[split$train_rows, , drop = FALSE]

  for (f in list(.fb_inla_slope_covariate_sd, .fb_inla_default_fixed_sd)) {
    # No spread under identical input, so identity is the right comparison.
    expect_identical(f(fb, engine), f(fb, engine))
    expect_identical(f(fb, extreme), f(fb, engine))
    expect_identical(f(fb, engine), f(fb, train))
  }
  # Each default against its formula on the training rows, written out here.
  expect_identical(
    .fb_inla_slope_covariate_sd(fb, engine),
    c(x_site = stats::sd(train$x))
  )
  expect_identical(
    .fb_inla_default_fixed_sd(fb, engine),
    c(
      x = 1e5 * stats::sd(train$y) / stats::sd(train$x),
      z1 = 1e5 * stats::sd(train$y) / stats::sd(train$z1)
    )
  )
})


test_that("an extreme covariate only in the assessment rows leaves an INLA fold refit's defaults unchanged", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbfi_quiet()
  d <- fbfi_data()
  split <- fbfi_folds(d)[[1L]]
  base <- .fb_sel_bayes_normalise_base(fbfi_cached(
    "slope_base",
    fb(y ~ x + (x || site), data = d, backend = "inla", verbose = FALSE)
  ))
  spec <- .fb_sel_bayes_base_spec(base)
  inputs <- .fb_sel_bayes_fold_inputs(
    base$data, split, "y", spec$protected, ~z1, random = spec$random_vars
  )
  extreme <- inputs$engine_data
  extreme$x[split$assess_rows[[1L]]] <- 25
  extreme$z1[split$assess_rows[[2L]]] <- -25

  clean_fit <- fbfi_cached(
    "slope_refit", fbfi_inla_refit(base, inputs$engine_data)
  )
  extreme_fit <- fbfi_inla_refit(base, extreme)
  clean <- clean_fit$extras$prior_completion
  moved <- extreme_fit$extras$prior_completion

  train <- d[split$train_rows, , drop = FALSE]
  # The slope reached the refit, and its bound was divided by the training
  # rows' sd(x).
  expect_identical(clean$covariate_sd, c(x_site = stats::sd(train$x)))
  expect_identical(moved$covariate_sd, clean$covariate_sd)
  expect_identical(moved$fixed_default_sd, clean$fixed_default_sd)
  expect_identical(
    clean$fixed_default_sd[["x"]],
    1e5 * stats::sd(train$y) / stats::sd(train$x)
  )
})


# --- (ii) the working response in fold compilation ----------------------- #

test_that("an offset model's fold default bound is five SDs of the training rows' modelled response", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbfi_quiet()
  d <- fbfi_data()
  split <- fbfi_folds(d)[[1L]]
  base <- .fb_sel_bayes_normalise_base(fbfi_cached(
    "offset_base",
    fb(y ~ x + offset(o) + (1 | site), data = d, backend = "inla",
       verbose = FALSE)
  ))
  spec <- .fb_sel_bayes_base_spec(base)
  uppers <- function(frame) {
    inputs <- .fb_sel_bayes_fold_inputs(
      frame, split, "y", spec$protected, ~z1, random = spec$random_vars
    )
    fit <- fbfi_inla_refit(base, inputs$engine_data)
    vapply(
      fit$fb$priors$specs,
      function(s) as.numeric(s$spec$args$upper),
      numeric(1L)
    )
  }
  extreme <- d
  extreme$y[split$assess_rows[[1L]]] <- 1e4
  train <- d[split$train_rows, , drop = FALSE]

  clean <- uppers(d)
  expect_true(all(clean == 5 * stats::sd(train$y - train$o)))
  expect_identical(uppers(extreme), clean)
})


test_that("the benchmark substrate reads a transformed response through the fit's resolver", {
  d <- fbfi_data()
  split <- fbfi_folds(d)[[1L]]
  extreme <- d
  extreme$y[split$assess_rows[[1L]]] <- 1e4
  for (response in c("log(y)", "I(10 * y)")) {
    modelled <- as.numeric(eval(str2lang(response), d))
    substrate <- .fb_sel_fold_substrate(
      d, split, response, ~x, ~z1, random = "site"
    )
    expect_identical(substrate$response, modelled)
    moved <- .fb_sel_fold_substrate(
      extreme, split, response, ~x, ~z1, random = "site"
    )
    expect_identical(moved$learned, substrate$learned)
  }
})


test_that("fb_select() scores a log(y) base on the benchmark route and refuses it on the Bayesian one", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbfi_quiet()
  d <- fbfi_data()
  base <- fbfi_cached(
    "log_base",
    flexybayes(log(y) ~ x, random = ~site, data = d, backend = "inla",
               verbose = FALSE)
  )
  sel <- fb_select(
    base, ~ z1 + z2,
    workflow = "unpenalized_benchmark",
    resamples = fb_resamples(unit = ~site, v = 3L),
    seed = 1L
  )
  expect_s3_class(sel, "fb_select")
  expect_true(all(is.finite(sel$scores$score)))
  # The Bayesian route sets the assessment rows' response to NA, and an
  # expression has no column to set; that refusal is unchanged.
  expect_error(
    fb_select(
      base, ~z1,
      resamples = fb_resamples(unit = ~site, v = 3L),
      seed = 1L
    ),
    "is an expression",
    fixed = TRUE
  )
})


# --- (iii) relabelling reaches every grouping column --------------------- #

test_that("new_cluster relabels exactly the grouping columns brms reads off the emitted formula", {
  skip_if_not_installed("brms")
  d <- fbfi_data()
  # Folds on rows split every site across training and assessment rows, so an
  # unrelabelled grouping level is shared between them. On folds that hold
  # whole sites out, relabelling changes no index, and the disjointness
  # check below could not see a column that was skipped.
  split <- .fb_sel_folds(fb_resamples(v = 3L), d, seed = 1L)[[1L]]
  expect_true(any(d$site[split$assess_rows] %in% d$site[split$train_rows]))
  cases <- list(
    intercept = fbfi_stub(d, y ~ x + (1 | site)),
    slope = fbfi_stub(d, y ~ x + (x || site)),
    three_way_brms = fbfi_stub(d, y ~ x + (1 | a:b:cc)),
    three_way = fbfi_stub(d, y ~ x, ~ a:b:cc),
    nested = fbfi_stub(d, y ~ x, ~ site:a),
    ide = fbfi_stub(d, y ~ x, ~ ide(site)),
    heterogeneous = fbfi_stub(d, y ~ x, ~ diag(env):site)
  )
  for (name in names(cases)) {
    base <- cases[[name]]
    spec <- .fb_sel_bayes_base_spec(base)
    inputs <- .fb_sel_bayes_fold_inputs(
      d, split, "y", spec$protected, ~z1,
      random = spec$random_vars, prediction = "new_cluster"
    )
    relabelled <- fbfi_relabelled(d, inputs)
    expect_identical(
      sort(relabelled),
      sort(fbfi_brms_groups(base$fb)),
      label = paste0(name, ": relabelled columns")
    )
    for (column in relabelled) {
      shared <- intersect(
        as.character(inputs$engine_data[[column]][split$assess_rows]),
        as.character(inputs$engine_data[[column]][split$train_rows])
      )
      expect_length(shared, 0L)
    }
    known <- .fb_sel_bayes_fold_inputs(
      d, split, "y", spec$protected, ~z1,
      random = spec$random_vars, prediction = "known_cluster"
    )
    expect_length(fbfi_relabelled(d, known), 0L)
  }
})


test_that("a known-covariance term refuses by name under new_cluster, before anything is fitted", {
  skip_if_not_installed("brms")
  d <- fbfi_data()
  K <- diag(nlevels(d$site))
  dimnames(K) <- list(levels(d$site), levels(d$site))
  base <- fbfi_stub(d, y ~ x, ~ vm(site, K), class = "flexybayes_brms")
  reached <- 0L
  refit <- function(...) {
    reached <<- reached + 1L
    stop("the refit was reached", call. = FALSE)
  }
  run <- function(prediction) {
    tryCatch(
      .fb_sel_bayes_cv(
        base, ~z1,
        resamples = fb_resamples(unit = ~site, v = 3L),
        prediction = prediction,
        metric = "elpd",
        seed = 1L,
        refit_fn = refit,
        probe_seconds = 1
      ),
      error = function(e) e
    )
  }

  err <- run("new_cluster")
  expect_s3_class(
    err, "flexybayes_refusal_selection_new_cluster_not_representable"
  )
  expect_match(conditionMessage(err), "vm(site)", fixed = TRUE)
  expect_identical(reached, 0L)

  # known_cluster keeps the levels, so the same base reaches the refit.
  err_known <- run("known_cluster")
  expect_identical(conditionMessage(err_known), "the refit was reached")
  expect_identical(reached, 1L)
})


test_that("a smooth of a covariate carries no grouping column and passes new_cluster", {
  d <- fbfi_data()
  base <- fbfi_stub(d, y ~ 1, ~ spl(x))
  expect_identical(.fb_sel_grouping_columns(base$fb)$columns, character(0))
  expect_null(.fb_sel_require_relabelable(base$fb, "new_cluster"))
})


test_that("an INLA fold refit indexes no assessment row's effect by a training row's level", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbfi_quiet()
  d <- fbfi_data()
  unit_split <- fbfi_folds(d)[[1L]]
  # Folds on rows split the sites, so a grouping left unrelabelled shares its
  # index values between training and assessment rows; see the brms block.
  row_split <- .fb_sel_folds(fb_resamples(v = 3L), d, seed = 1L)[[1L]]
  expect_true(any(d$site[row_split$assess_rows] %in% d$site[row_split$train_rows]))
  # Adversarial: a level of `cc`, and so every interaction cell carrying it,
  # exists only in one site, which the unit fold holds out.
  lone <- d$site == d$site[unit_split$assess_rows[[1L]]]
  d_lone <- d
  d_lone$cc <- factor(as.character(d_lone$cc), levels = c("c1", "c2", "c3"))
  d_lone$cc[lone] <- "c3"

  slope <- .fb_sel_bayes_normalise_base(fbfi_cached(
    "slope_base",
    fb(y ~ x + (x || site), data = d, backend = "inla", verbose = FALSE)
  ))
  three_way <- .fb_sel_bayes_normalise_base(fbfi_cached(
    "three_way_base",
    flexybayes(y ~ x, random = ~ a:b:cc, data = d, backend = "inla",
               verbose = FALSE)
  ))
  three_way_lone <- .fb_sel_bayes_normalise_base(fbfi_cached(
    "three_way_lone_base",
    flexybayes(y ~ x, random = ~ a:b:cc, data = d_lone, backend = "inla",
               verbose = FALSE)
  ))
  refit_on <- function(base, split) {
    spec <- .fb_sel_bayes_base_spec(base)
    inputs <- .fb_sel_bayes_fold_inputs(
      base$data, split, "y", spec$protected, ~z1, random = spec$random_vars
    )
    list(fit = fbfi_inla_refit(base, inputs$engine_data), split = split)
  }
  cases <- list(
    slope = refit_on(slope, row_split),
    three_way = refit_on(three_way, row_split),
    three_way_lone_level = refit_on(three_way_lone, unit_split)
  )

  for (name in names(cases)) {
    fit <- cases[[name]]$fit
    split <- cases[[name]]$split
    indices <- fbfi_f_indices(fit$extras$formula)
    expect_gt(length(indices), 0L, label = paste0(name, ": f() terms"))
    data_inla <- fit$inla$.args$data
    for (index in indices) {
      shared <- intersect(
        stats::na.omit(data_inla[[index]][split$assess_rows]),
        stats::na.omit(data_inla[[index]][split$train_rows])
      )
      expect_length(shared, 0L)
    }
  }
})


# --- (iv) the fold refit carries a brms-grammar base's random terms ------ #

test_that("a fold refit of a brms-grammar base re-issues the base's random terms", {
  d <- fbfi_data()
  capture <- function(base) {
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

  slope <- capture(fbfi_stub(d, y ~ x + (x || site)))
  expect_identical(
    attr(stats::terms(slope$fixed), "term.labels"),
    c("x", "z1", "x || site")
  )
  expect_null(slope$random)

  asreml <- capture(fbfi_stub(d, y ~ x, ~site))
  expect_identical(attr(stats::terms(asreml$fixed), "term.labels"), c("x", "z1"))
  expect_identical(deparse(asreml$random), "~site")
})


# --- (v) the benchmark route's random terms and learned quantities ------- #

test_that("the benchmark covariance builder receives the base's random terms, or the route refuses by name", {
  d <- fbfi_data()
  seen <- new.env(parent = emptyenv())
  local_mocked_bindings(
    .fb_sel_cv = function(data, response, protected, candidates = NULL,
                          random = character(0), ...) {
      seen$random <- random
      stop("the covariance builder was reached", call. = FALSE)
    }
  )
  run <- function(base) {
    tryCatch(
      fb_select(
        base, ~z1,
        workflow = "unpenalized_benchmark",
        resamples = fb_resamples(unit = ~site, v = 3L),
        seed = 1L
      ),
      error = function(e) e
    )
  }

  seen$random <- NULL
  err <- run(fbfi_stub(d, y ~ x, ~ site:a))
  expect_identical(conditionMessage(err), "the covariance builder was reached")
  # One interaction term, not two crossed intercepts.
  expect_identical(seen$random, "site:a")

  for (base in list(
    fbfi_stub(d, y ~ x + (x || site)),
    fbfi_stub(d, y ~ x, ~ a:b:cc)
  )) {
    seen$random <- NULL
    err <- run(base)
    expect_s3_class(
      err, "flexybayes_refusal_selection_covariance_not_buildable"
    )
    expect_null(seen$random)
  }
})


test_that("an extreme value only in the assessment rows moves nothing the benchmark and optimiser folds learn", {
  d <- fbfi_data()
  split <- fbfi_folds(d)[[1L]]
  extreme <- d
  extreme$y[split$assess_rows[[1L]]] <- 1e4
  extreme$z1[split$assess_rows[[2L]]] <- -1e4
  build <- function(frame) {
    .fb_sel_fold_substrate(frame, split, "y", ~x, ~ z1 + z2, random = "site")
  }
  clean <- build(d)
  moved <- build(extreme)
  # Centring and scaling, the variance components, the whitening factor and
  # the pilot fit.
  expect_identical(moved$learned, clean$learned)
  expect_identical(moved$blocks$factor_tt, clean$blocks$factor_tt)
  # The optimiser's whitened, protected-residualised payload.
  expect_identical(
    .fb_sel_glmnet_reduction(moved, c(1L, 2L))$payload,
    .fb_sel_glmnet_reduction(clean, c(1L, 2L))$payload
  )
})


test_that("the lasso grid starts at the largest of the folds' own training-row anchors", {
  skip_if_not_installed("glmnet")
  d <- fbfi_data()
  cv <- .fb_sel_cv(
    d, "y", ~x, ~ z1 + z2, random = "site",
    resamples = fb_resamples(unit = ~site, v = 3L), seed = 1L
  )
  blocks <- list(list(
    family = "elastic_net",
    columns = c(1L, 2L),
    hyperparameters = list(lambda = 1e-3, alpha = 0.5)
  ))
  path <- .fb_sel_glmnet_path(
    cv$substrates, blocks, folds = cv$folds,
    decision = list(n_lambda = 4L), seed = 1L
  )
  alpha <- path$penalty$alpha_composite
  anchors <- vapply(
    cv$substrates,
    function(s) {
      payload <- .fb_sel_glmnet_reduction(s, c(1L, 2L))$payload
      max(abs(drop(crossprod(payload$x, payload$y)))) / alpha
    },
    numeric(1L)
  )
  grid <- path$path_coordinates$lambda_canonical[
    path$path_coordinates$arm == "benchmark"
  ]
  # Non-vacuous: the data anchor, not the declared lambda, tops the grid.
  expect_gt(max(anchors), path$penalty$lambda_composite)
  # The grid is evenly spaced in log lambda, so its first point is the
  # anchor through exp(log()), the arithmetic the grid itself uses.
  expect_identical(max(grid), exp(log(max(anchors))))
})


# --- (vi) adversarial folds complete or refuse by name ------------------- #

test_that("a fixed-factor level only in the assessment rows refuses by name on both fold builders", {
  d <- fbfi_data()
  split <- fbfi_folds(d)[[1L]]
  lone <- d$site == d$site[split$assess_rows[[1L]]]
  d$f1 <- factor(ifelse(lone, "lone", "common"))
  expect_error(
    .fb_sel_bayes_fold_inputs(d, split, "y", ~x, ~f1, random = "site"),
    class = "flexybayes_refusal_selection_fold_level_absent"
  )
  expect_error(
    .fb_sel_fold_substrate(d, split, "y", ~x, ~f1, random = "site"),
    class = "flexybayes_refusal_selection_fold_level_absent"
  )
})


test_that("an interaction cell only in the assessment rows completes on the benchmark substrate", {
  d <- fbfi_data()
  split <- fbfi_folds(d)[[1L]]
  lone <- d$site == d$site[split$assess_rows[[1L]]]
  # The cell (a2, b3) exists only in the lone site; both levels occur in
  # training rows.
  d$b[!lone & d$a == "a2" & d$b == "b3"] <- "b1"
  d$a[lone] <- "a2"
  d$b[lone] <- "b3"
  expect_true(any(d$a[split$train_rows] == "a2"))
  expect_true(any(d$b[split$train_rows] == "b3"))
  expect_false(any(d$a[split$train_rows] == "a2" & d$b[split$train_rows] == "b3"))

  substrate <- .fb_sel_fold_substrate(d, split, "y", ~x, ~ a:b, random = "site")
  expect_true(all(is.finite(substrate$learned$pilot_beta)))
  extreme <- d
  extreme$y[split$assess_rows[[1L]]] <- 1e4
  moved <- .fb_sel_fold_substrate(extreme, split, "y", ~x, ~ a:b, random = "site")
  expect_identical(moved$learned, substrate$learned)
})
