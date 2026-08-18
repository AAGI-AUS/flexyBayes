# =============================================================================
# The design-preserving missing-response layer.
#
# A missing response removes an observation, not a design node. Where the
# model carries a covariance indexed by the design -- a separable AR1(row) x
# AR1(col) residual over a field trial -- deleting the row of a lost plot
# changes the index set the covariance is built over, so what gets fitted is
# no longer the model that was written down.
#
# na_action = "augment" keeps the design intact: rows whose response is
# missing are retained and carried as latent quantities, and design cells
# absent from the data frame are reinstated when their model variables are
# determinable. Under ignorability the parameter posterior is the same as
# under deletion -- augmentation preserves the representation, not
# information -- but the representation is what makes the fit possible at
# all here, since deletion breaks the grid and the emit refuses.
#
# References: Yates (1933) missing-plot technique; Houtman and Speed (1984);
# Verbyla and Cullis (1992); Rubin (1976) on ignorability; Tanner and Wong
# (1987) on data augmentation.
# =============================================================================

.ar1_m <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))

# One observation per node on a complete n_row x n_col grid with a
# separable AR1 field plus nugget.
.grid_sim <- function(n_row = 8L, n_col = 6L, rho_r = 0.6, rho_c = 0.3,
                      sd_s = 1.2, sd_e = 0.4, seed = 2026L) {
  set.seed(seed)
  g <- expand.grid(col = seq_len(n_col), row = seq_len(n_row))
  L <- t(chol(sd_s^2 * kronecker(.ar1_m(n_row, rho_r), .ar1_m(n_col, rho_c))))
  g$y <- 20 + as.numeric(L %*% stats::rnorm(n_row * n_col)) +
    stats::rnorm(n_row * n_col, 0, sd_e)
  g$row <- factor(g$row)
  g$col <- factor(g$col)
  g
}

.hyper <- function(fit) {
  s <- fit$inla$summary.hyperpar
  stats::setNames(s[, "mean"], rownames(s))
}

# ---------------------------------------------------------------- #
# 1. The layer itself, without fitting                              #
# ---------------------------------------------------------------- #

test_that("augment reinstates absent design cells with a missing response", {
  g <- .grid_sim(4L, 3L)
  ir <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(),
    random_terms = list(list(
      type = "ar1_spatial", row_var = "row", col_var = "col"
    )),
    residual_terms = list(list(type = "units")),
    source = "asreml"
  )
  short <- g[-c(2L, 7L), ]
  res <- flexyBayes:::.fb_apply_na_action(ir, short, "augment")

  expect_equal(nrow(res$data), nrow(g))
  expect_equal(res$meta$n_cells_completed, 2L)
  expect_equal(sum(is.na(res$data$y)), 2L)
  # Every grid node present exactly once -- the condition the emit checks.
  expect_equal(
    nrow(unique(res$data[, c("row", "col")])),
    nlevels(g$row) * nlevels(g$col)
  )
  expect_false(anyDuplicated(paste(res$data$row, res$data$col)) > 0L)
})

test_that("omit drops missing-response rows", {
  g <- .grid_sim(4L, 3L)
  g$y[c(1L, 5L)] <- NA
  ir <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(), random_terms = list(),
    residual_terms = list(list(type = "units")), source = "asreml"
  )
  res <- flexyBayes:::.fb_apply_na_action(ir, g, "omit")
  expect_equal(nrow(res$data), nrow(g) - 2L)
  expect_false(anyNA(res$data$y))
  expect_equal(res$meta$n_missing_response, 2L)
})

test_that("fail refuses a missing response", {
  g <- .grid_sim(4L, 3L)
  g$y[1L] <- NA
  ir <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(), random_terms = list(),
    residual_terms = list(list(type = "units")), source = "asreml"
  )
  expect_error(
    flexyBayes:::.fb_apply_na_action(ir, g, "fail"),
    class = "flexybayes_missing_response_refused"
  )
})

test_that("a missing covariate is refused under every na_action", {
  # The device carries a missing RESPONSE as a latent quantity. A filled-in
  # predictor is a fabricated observation, and no ignorability argument
  # makes it safe -- so this refuses even under augment.
  g <- .grid_sim(4L, 3L)
  g$x <- stats::rnorm(nrow(g))
  g$x[3L] <- NA
  ir <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(list(type = "continuous", var = "x")),
    random_terms = list(),
    residual_terms = list(list(type = "units")), source = "asreml"
  )
  for (act in c("augment", "omit", "fail")) {
    expect_error(
      flexyBayes:::.fb_apply_na_action(ir, g, act),
      class = "flexybayes_missing_covariate_not_supported",
      label = act
    )
  }
})

test_that("an absent cell with an unrecorded design factor refuses", {
  # Reinstating the cell would mean inventing a genotype for a plot whose
  # assignment nobody wrote down. Refuse and say what to supply instead.
  g <- .grid_sim(4L, 3L)
  g$geno <- factor(rep(letters[1:4], length.out = nrow(g)))
  ir <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(),
    random_terms = list(
      list(type = "simple", var = "geno"),
      list(type = "ar1_spatial", row_var = "row", col_var = "col")
    ),
    residual_terms = list(list(type = "units")),
    source = "asreml"
  )
  expect_error(
    flexyBayes:::.fb_apply_na_action(ir, g[-2L, ], "augment"),
    class = "flexybayes_augment_cell_not_determinable"
  )
})

test_that("complete data is untouched by any na_action", {
  g <- .grid_sim(4L, 3L)
  ir <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(),
    random_terms = list(list(
      type = "ar1_spatial", row_var = "row", col_var = "col"
    )),
    residual_terms = list(list(type = "units")),
    source = "asreml"
  )
  for (act in c("augment", "omit", "fail")) {
    res <- flexyBayes:::.fb_apply_na_action(ir, g, act)
    expect_equal(nrow(res$data), nrow(g), label = act)
  }
})

# ---------------------------------------------------------------- #
# 2. End to end: the design survives, and deletion does not          #
# ---------------------------------------------------------------- #

test_that("augment fits a holed grid that omit refuses", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .grid_sim()
  lost <- c(5L, 19L, 33L)

  na_row <- g
  na_row$y[lost] <- NA
  absent <- g[-lost, ]

  for (dat in list(na_row, absent)) {
    expect_error(
      suppressMessages(flexybayes(
        y ~ 1, random = ~ ar1(row):ar1(col), data = dat,
        backend = "inla", na_action = "omit", verbose = FALSE
      )),
      class = "flexybayes_refusal_ar1_spatial_requires_complete_grid"
    )
    expect_no_error(suppressMessages(flexybayes(
      y ~ 1, random = ~ ar1(row):ar1(col), data = dat,
      backend = "inla", na_action = "augment", verbose = FALSE
    )))
  }
})

test_that("the two routes to a completed grid give the same analysis", {
  # Marking a present row NA and reinstating an absent row are two ways of
  # describing the same missing plot, so they must produce the same fit.
  # This is the layer's internal consistency check; agreement with an
  # engine-independent oracle is a separate matter.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .grid_sim()
  lost <- c(5L, 19L, 33L)
  na_row <- g
  na_row$y[lost] <- NA

  f_na <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(row):ar1(col), data = na_row,
    backend = "inla", na_action = "augment", verbose = FALSE
  ))
  f_absent <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(row):ar1(col), data = g[-lost, ],
    backend = "inla", na_action = "augment", verbose = FALSE
  ))

  h1 <- .hyper(f_na)
  h2 <- .hyper(f_absent)
  expect_identical(names(h1), names(h2))

  # The two routes are compared at the precision each parameter actually
  # supports, which is not the same for all four.
  #
  # The Gaussian observation precision is weakly identified in this model:
  # the simulated nugget is 10% of the total, so the likelihood is nearly
  # flat in that direction and the optimiser stops at different points on the
  # ridge. Measured over three identical repeats of the SAME data, this
  # parameter alone spans 2.6425 / 2.6627 / 2.7464 -- a 3.9% run-to-run
  # spread, with no missing-data difference involved at all. The other three
  # hyperparameters repeat to between 0.04% and 0.4%.
  #
  # A single 2% tolerance across all four therefore fails at random on the
  # first one (it did, in a full-suite run) while being far looser than the
  # rest deserve. Splitting the comparison keeps the tight check where it is
  # meaningful and states plainly why one parameter cannot carry it.
  obs_prec <- grepl("Precision for the Gaussian observations", names(h1),
                    fixed = TRUE)
  if (any(obs_prec)) {
    expect_equal(unname(h1[obs_prec]), unname(h2[obs_prec]), tolerance = 0.08)
  }
  expect_equal(unname(h1[!obs_prec]), unname(h2[!obs_prec]), tolerance = 0.02)
})

test_that("the fit records what the layer did", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .grid_sim()
  g$y[c(5L, 19L)] <- NA
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(row):ar1(col), data = g,
    backend = "inla", na_action = "augment", verbose = FALSE
  ))
  m <- fit$extras$na_action
  expect_identical(m$na_action, "augment")
  expect_equal(m$n_missing_response, 2L)
  expect_setequal(m$design_index_vars, c("row", "col"))
})

# =============================================================================
# High missing fraction: the estimand degrades, and the fit says so (S8)
#
# The augmentation identity is algebraic and holds at any missing fraction.
# What weakens past roughly a third unobserved is the restricted likelihood
# it targets -- variance components reach the boundary, the surface
# flattens, and two correct implementations legitimately stop at different
# points. The warning states that once per session so a user comparing
# against ASReml or lme4 does not read a disagreement as a defect.
# =============================================================================

.missing_frac_data <- function(n = 100L, frac_missing = 0.35, seed = 606L) {
  set.seed(seed)
  d <- data.frame(
    y = stats::rnorm(n),
    x = stats::rnorm(n),
    g = factor(rep(seq_len(10L), length.out = n))
  )
  n_na <- as.integer(round(frac_missing * n))
  d$y[seq_len(n_na)] <- NA_real_
  d
}

.apply_augment <- function(d) {
  fb <- flexyBayes:::.build_ir_polymorphic(
    fixed = y ~ x,
    random = ~g,
    residual = NULL,
    data = d,
    family = "gaussian",
    link = NULL,
    weights = NULL,
    known_matrices = list(),
    prior = NULL,
    prior_fixed_sd = 100,
    prior_vc_sd = 1,
    syntax = "asreml"
  )
  flexyBayes:::.fb_apply_na_action(fb, d, "augment")
}

test_that("augmenting past 30 percent missing warns about the estimand", {
  local_clean_emit_state()
  d <- .missing_frac_data(frac_missing = 0.35)

  expect_warning(
    out <- .apply_augment(d),
    "carry no observed response"
  )
  expect_equal(out$meta$missing_fraction, 0.35, tolerance = 1e-8)
})

test_that("the high-missingness warning names the boundary, not the device", {
  local_clean_emit_state()
  d <- .missing_frac_data(frac_missing = 0.35)

  w <- tryCatch(.apply_augment(d), warning = function(w) w)
  expect_s3_class(w, "warning")
  expect_match(conditionMessage(w), "boundary")
  expect_match(conditionMessage(w), "property of the estimand")
  expect_match(conditionMessage(w), "augmentation identity itself is")
})

test_that("augmenting at 20 percent missing is silent", {
  local_clean_emit_state()
  d <- .missing_frac_data(frac_missing = 0.20)

  expect_warning(out <- .apply_augment(d), NA)
  expect_equal(out$meta$missing_fraction, 0.20, tolerance = 1e-8)
})

test_that("the high-missingness warning fires once and is silenceable", {
  local_clean_emit_state()
  d <- .missing_frac_data(frac_missing = 0.35)

  expect_warning(.apply_augment(d), "carry no observed response")
  # Once per session, in the manner of the routing notes: a simulation
  # study fitting the same design a thousand times says this once.
  expect_warning(.apply_augment(d), NA)

  local_clean_emit_state()
  withr::local_options(flexyBayes.silence_high_missingness_warning = TRUE)
  expect_warning(.apply_augment(d), NA)
})


# =============================================================================
# The record the layer leaves behind (0.9.1)
#
# One builder writes the record on all three response policies, so the
# fields cannot differ between branches. Two numbers carry one name each
# across the record, summary(fit) and nobs(type = ): `n_design` is how
# many rows the engine was handed and `n_observed` how many of those
# carried an observed response. `missing_fraction` is derived from the
# pair on every path, which is what makes it comparable between them.
# =============================================================================

test_that("the record carries the same fields on every response policy", {
  g <- .grid_sim(4L, 3L)
  ir <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(), random_terms = list(),
    residual_terms = list(list(type = "units")), source = "asreml"
  )
  fields <- c(
    "na_action", "na_covariate", "n_missing_response", "n_cells_completed",
    "design_index_vars", "n_design", "n_observed",
    "n_covariate_rows_dropped", "missing_fraction"
  )
  for (act in c("augment", "omit", "fail")) {
    meta <- flexyBayes:::.fb_apply_na_action(ir, g, act)$meta
    expect_identical(names(meta), fields, label = act)
    expect_identical(meta$na_action, act, label = act)
  }
})

test_that("n_design and n_observed count the augmented design", {
  g <- .grid_sim(4L, 3L)
  g$y[c(2L, 5L)] <- NA
  ir <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(),
    random_terms = list(list(
      type = "ar1_spatial", row_var = "row", col_var = "col"
    )),
    residual_terms = list(list(type = "units")), source = "asreml"
  )

  # Two lost yields on a complete grid: the design keeps all 12 cells,
  # ten of which were observed.
  meta <- flexyBayes:::.fb_apply_na_action(ir, g, "augment")$meta
  expect_equal(meta$n_design, 12L)
  expect_equal(meta$n_observed, 10L)
  expect_equal(meta$missing_fraction, 2 / 12, tolerance = 1e-12)

  # A cell absent from the data frame entirely is reinstated, so it
  # counts in the design and not in the observations.
  meta2 <- flexyBayes:::.fb_apply_na_action(ir, g[-1L, ], "augment")$meta
  expect_equal(meta2$n_cells_completed, 1L)
  expect_equal(meta2$n_design, 12L)
  expect_equal(meta2$n_observed, 9L)
})

test_that("missing_fraction is present and zero on the omit path", {
  # It used to be absent there, so a consumer expecting the five
  # documented fields on every path got NULL on one of them.
  g <- .grid_sim(4L, 3L)
  g$y[c(2L, 5L)] <- NA
  ir <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(), random_terms = list(),
    residual_terms = list(list(type = "units")), source = "asreml"
  )
  meta <- flexyBayes:::.fb_apply_na_action(ir, g, "omit")$meta
  expect_false(is.null(meta$missing_fraction))
  expect_equal(meta$missing_fraction, 0)
  expect_equal(meta$n_missing_response, 2L)
  expect_equal(meta$n_design, 10L)
  expect_equal(meta$n_observed, 10L)
})

test_that("fail leaves the design grid alone", {
  # `fail` refuses a missing response. Completing the grid would create
  # the missing responses it exists to refuse, so it does not run.
  g <- .grid_sim(4L, 3L)
  ir <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(),
    random_terms = list(list(
      type = "ar1_spatial", row_var = "row", col_var = "col"
    )),
    residual_terms = list(list(type = "units")), source = "asreml"
  )
  res <- flexyBayes:::.fb_apply_na_action(ir, g[-1L, ], "fail")
  expect_equal(nrow(res$data), nrow(g) - 1L)
  expect_identical(res$meta$na_action, "fail")
  expect_equal(res$meta$n_cells_completed, 0L)
})

test_that("the fit records the requested policy in its call record", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .grid_sim()
  g$y[c(5L, 19L)] <- NA
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(row):ar1(col), data = g,
    backend = "inla", na_action = "augment", verbose = FALSE
  ))
  expect_identical(fit$extras$call_info$na_action, "augment")
  expect_equal(fit$extras$na_action$n_design, nrow(g))
  expect_equal(fit$extras$na_action$n_observed, nrow(g) - 2L)
})

test_that("an asreml na.method() value fits the same model as the string", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .grid_sim()
  g$y[c(5L, 19L)] <- NA
  as_list <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(row):ar1(col), data = g, backend = "inla",
    na_action = .asreml_na_method_recorded("explicit"), verbose = FALSE
  ))
  expect_identical(as_list$extras$na_action$na_action, "augment")
  expect_identical(as_list$extras$call_info$na_action, "augment")
  expect_equal(stats::nobs(as_list), nrow(g))
})


# =============================================================================
# The recorded call (0.9.1)
#
# update() rebuilds the original call from `$extras$call_info`, and
# refuses when the record is short, because a re-fit that substituted a
# default for something the user set would return a different model under
# the same name. The INLA emit recorded six of the fifteen arguments brms
# records, and that gap -- nothing about the engine -- was the whole of
# the refusal.
# =============================================================================

test_that("an INLA fit records the same call fields a brms fit does", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .grid_sim()
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(row):ar1(col), data = g,
    backend = "inla", verbose = FALSE
  ))

  expected <- c(
    "fixed", "random", "residual", "data_name", "family", "link",
    "known_matrices", "weights", "n_samples", "warmup", "chains",
    "seed", "control", "prior_fixed_sd", "prior_vc_sd", "na_action",
    "backend", "aggregate", "verbose"
  )
  expect_identical(names(fit$extras$call_info), expected)

  # The sampler settings are recorded as what INLA used, which is
  # nothing: a nested Laplace approximation runs no chains and discards
  # no warmup, so recording the numbers the call happened to carry would
  # put a sampler on the record of a fit that never ran one.
  expect_null(fit$extras$call_info$n_samples)
  expect_null(fit$extras$call_info$warmup)
  expect_null(fit$extras$call_info$chains)
  expect_identical(fit$extras$call_info$prior_fixed_sd, 100)
})

test_that("update() re-fits an INLA model rather than refusing it", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .grid_sim()
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(row):ar1(col), data = g,
    backend = "inla", verbose = FALSE
  ))
  refit <- suppressMessages(stats::update(fit, random = ~ ar1(row):ar1(col)))
  expect_s3_class(refit, "flexybayes_inla")
  expect_identical(refit$extras$call_info$na_action, "augment")
})


# =============================================================================
# What the refusals say (0.9.1)
#
# Each of the three messages a user meets on a trial with lost plots now
# names what ASReml does in the same situation, and what to type instead.
# A refusal that only says no leaves the reader to guess whether the
# package is stricter than ASReml or simply less capable.
# =============================================================================

test_that("the broken-grid refusal cites ASReml and names the remedy", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .grid_sim()
  g$y[c(5L, 19L)] <- NA

  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ 1, random = ~ ar1(row):ar1(col), data = g,
      backend = "inla", na_action = "omit", verbose = FALSE
    )),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_ar1_spatial_refusal")
  msg <- conditionMessage(err)

  # It counts what is actually there rather than describing the class of
  # problem in the abstract.
  expect_match(msg, "rows for the", fixed = TRUE)
  expect_match(msg, "nodes", fixed = TRUE)
  # It says ASReml refuses the same data for the same reason, so a
  # reader does not read this as flexyBayes being the weaker tool.
  expect_match(msg, "ASReml", fixed = TRUE)
  expect_match(msg, "residual = ~ ar1:ar1", fixed = TRUE)
  # And it names a remedy the reader can type.
  expect_match(msg, "na_action = \"augment\"", fixed = TRUE)
  expect_match(msg, "field-book rows", fixed = TRUE)
})

test_that("the undeterminable-cell refusal names LANCER and the field book", {
  g <- .grid_sim(4L, 3L)
  g$geno <- factor(rep(letters[1:4], length.out = nrow(g)))
  ir <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(),
    random_terms = list(
      list(type = "simple", var = "geno"),
      list(type = "ar1_spatial", row_var = "row", col_var = "col")
    ),
    residual_terms = list(list(type = "units")),
    source = "asreml"
  )
  err <- tryCatch(
    flexyBayes:::.fb_apply_na_action(ir, g[-2L, ], "augment"),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_augment_cell_not_determinable")
  msg <- conditionMessage(err)
  expect_match(msg, "geno", fixed = TRUE)
  expect_match(msg, "LANCER", fixed = TRUE)
  expect_match(msg, "nin89", fixed = TRUE)
  expect_match(msg, "field-book rows", fixed = TRUE)
  expect_match(msg, "na_action = \"omit\"", fixed = TRUE)
})


# ---------------------------------------------------------------- #
# 7. The unobserved cells reach an accessor (WP-E2)                 #
# ---------------------------------------------------------------- #

test_that("summary()$missing names the lost plots on an INLA field fit", {
  # ASReml's `mv` factor. The engine already computes a posterior for each
  # unobserved cell -- INLA treats an NA response as a latent prediction
  # target -- so the test is that the accessor reports what the fit holds,
  # at the right rows, on the right scale.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .grid_sim()
  lost <- c(5L, 19L)
  g$y[lost] <- NA
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(row):ar1(col), data = g,
    backend = "inla", na_action = "augment", verbose = FALSE
  ))

  captured <- utils::capture.output(s <- summary(fit))
  mv <- s$missing
  expect_s3_class(mv, "data.frame")
  expect_identical(
    names(mv)[1:5],
    c("row", "estimate", "std.error", "conf.low", "conf.high")
  )
  expect_identical(nrow(mv), 2L)
  expect_identical(
    nrow(mv),
    as.integer(stats::nobs(fit) - stats::nobs(fit, type = "observed"))
  )
  expect_identical(mv$row, lost)
  expect_true(all(is.finite(mv$estimate)))
  expect_true(all(mv$std.error > 0))
  expect_true(all(mv$conf.low <= mv$conf.high))

  # The design index the covariance is built over is carried alongside, so
  # a reader can see which plot each row is.
  expect_true(all(c("row.index", "col") %in% names(mv)))
  expect_identical(as.character(mv$row.index), as.character(g$row[lost]))
  expect_identical(as.character(mv$col), as.character(g$col[lost]))

  # One construction behind both accessors.
  expect_identical(stats::coef(fit, what = "missing"), mv)

  # And the values are INLA's own fitted-value marginals, read from the
  # full table rather than the head() copy on $extras.
  expect_equal(
    mv$estimate,
    as.numeric(fit$inla$summary.fitted.values[["mean"]][lost]),
    tolerance = 1e-12
  )
})

test_that("a complete grid reports a typed zero-row missing table", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  g <- .grid_sim()
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(row):ar1(col), data = g,
    backend = "inla", na_action = "augment", verbose = FALSE
  ))
  mv <- stats::coef(fit, what = "missing")
  expect_s3_class(mv, "data.frame")
  expect_identical(nrow(mv), 0L)
  expect_identical(
    names(mv),
    c("row", "estimate", "std.error", "conf.low", "conf.high")
  )
})

test_that("a non-identity link withholds the table and says why", {
  # INLA computes a fitted value for an NA response on the identity link
  # unless control.predictor names another, so on a Poisson fit those rows
  # would be on the linear-predictor scale while the rows around them are
  # on the response scale. The table is withheld rather than mixed, and
  # the reason is not left to be discovered.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  local_clean_emit_state()

  set.seed(707L)
  d <- data.frame(
    g = factor(rep(letters[1:6], each = 8L)),
    x = stats::rnorm(48L)
  )
  d$y <- stats::rpois(48L, exp(0.6 + 0.3 * d$x))
  d$y[c(2L, 17L)] <- NA
  fit <- suppressMessages(suppressWarnings(flexybayes(
    y ~ x, random = ~g, data = d, family = "poisson",
    backend = "inla", na_action = "augment", verbose = FALSE
  )))

  expect_warning(
    mv <- stats::coef(fit, what = "missing"),
    "linear-predictor scale"
  )
  expect_identical(nrow(mv), 0L)
  # The count is still reachable, which is the point of saying so.
  expect_equal(
    stats::nobs(fit) - stats::nobs(fit, type = "observed"),
    2
  )
  # Once per session, not once per call.
  expect_no_warning(stats::coef(fit, what = "missing"))
})
