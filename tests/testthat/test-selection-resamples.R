# =============================================================================
# Selection resampling: oracles 4, 7 and 9 of the correctness ledger.
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md sections 3.3, 4.1, 4.2,
# 5.4, 9 and 12.
#
#   Oracle 4  no fold leakage. Folds fall on original units before any
#             rotation, and every quantity a fold declares as learned is a
#             function of its training rows alone. The gate is a falsifier:
#             the fold is rebuilt with its assessment rows corrupted, and
#             anything that moves was read from them.
#   Oracle 7  the prediction estimand. The conditional mean and the
#             conditional covariance match dense multivariate-normal algebra;
#             a genuinely new cluster has V_ST = 0; and new_cluster and
#             known_cluster differ in the direction that makes new_cluster
#             the honest default.
#   Oracle 9  the scoring unit. Aggregation happens at the declared
#             independent cluster, a row-level standard error understates it
#             on correlated data, and one-standard-error selection reproduces
#             across reruns under the same seed.
#
# Every gate here ships with the defect that makes it red. Two of the gates
# audit A25 and A27 examined passed for every input because nobody had ever
# made them fail on purpose, and a gate nobody has seen fail is not evidence.
# The leakage gate carries two of its falsifiers as standing tests -- a
# builder that learns its scaling on every row, and one that learns its
# covariance on every row -- because a leakage test that has never been red is
# worth nothing at all.
#
# Almost nothing here fits a model. The prediction algebra is exercised
# against a hand-built dense covariance, and the variance components come from
# the module's own deterministic REML plugin, so the file runs in seconds and
# its numbers do not move between runs.
# =============================================================================

# --- fixtures ---------------------------------------------------------- #

# A strongly clustered design: twelve sites, six rows each, a between-site
# standard deviation of 2 against a within-site 0.3, and four site-constant
# covariates with no true effect. The ratio is what makes the difference
# between a unit holdout and a row holdout visible; the null covariates are
# what a wrongly-held-out cluster lets a model exploit.
fb_sel_rs_data <- function(
  n_cluster = 12L,
  per = 6L,
  sd_u = 2,
  sd_e = 0.3,
  n_static = 4L,
  seed = 11L
) {
  .fb_sel_with_seed(seed, {
    n <- n_cluster * per
    site <- factor(rep(seq_len(n_cluster), each = per))
    x <- stats::rnorm(n)
    u <- stats::rnorm(n_cluster, sd = sd_u)
    y <- 1 + 0.5 * x + u[as.integer(site)] + stats::rnorm(n, sd = sd_e)
    d <- data.frame(y = y, x = x, site = site)
    static <- matrix(stats::rnorm(n_cluster * n_static), nrow = n_cluster)
    for (j in seq_len(n_static)) {
      d[[paste0("z", j)]] <- static[as.integer(site), j]
    }
    d
  })
}


# The builder the module ships, wrapped so the leakage gate can call it.
fb_sel_rs_builder <- function(data, split) {
  .fb_sel_fold_substrate(
    data,
    split,
    response = "y",
    protected = ~x,
    candidates = ~ z1 + z2 + z3 + z4,
    random = "site"
  )
}


# --- fb_resamples(): a unit is not a specification --------------------- #

test_that("a grouping unit without a scheme refuses", {
  # Contract section 3.3. `unit = ~ LOCATION` says what may not be split
  # across a fold boundary; it does not say how the splitting is done, and
  # letting the default `v` supply the answer would make a grouped design
  # silently ten-fold.
  expect_error(
    fb_resamples(unit = ~site),
    "no resampling scheme",
    fixed = TRUE
  )
})


test_that("a grouping unit with a scheme is a specification", {
  expect_s3_class(fb_resamples(unit = ~site, v = 4L), "fb_resamples")
  expect_s3_class(fb_resamples(unit = ~site, repeats = 3L), "fb_resamples")
  expect_s3_class(
    fb_resamples(unit = ~site, fold_id = rep(1:2, 4)),
    "fb_resamples"
  )
  # No unit named: the unit-free defaults are a specification on their own.
  expect_s3_class(fb_resamples(), "fb_resamples")
})


test_that("an explicit partition and a scheme are not both accepted", {
  expect_error(
    fb_resamples(unit = ~site, v = 3L, fold_id = rep(1:2, 4)),
    "two answers to the same question",
    fixed = TRUE
  )
  expect_error(
    fb_resamples(fold_id = rep(1:2, 4), repeats = 3L),
    "report the repetition as independent evidence",
    fixed = TRUE
  )
})


test_that("the scheme arguments are checked before anything uses them", {
  expect_error(fb_resamples(v = 1L), "at least 2", fixed = TRUE)
  expect_error(fb_resamples(v = 2.5), "whole number", fixed = TRUE)
  expect_error(fb_resamples(repeats = 0L), "at least 1", fixed = TRUE)
  expect_error(fb_resamples(unit = 42), "one-sided formula", fixed = TRUE)
})


# --- ORACLE 4a: folds fall on original units --------------------------- #

test_that("a unit is never split across a fold boundary", {
  d <- fb_sel_rs_data()
  folds <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)

  expect_length(folds, 4L)
  for (split in folds) {
    train_sites <- unique(as.character(d$site)[split$train_rows])
    assess_sites <- unique(as.character(d$site)[split$assess_rows])
    # The whole content of "folds on original units": no site contributes
    # rows to both sides of the same split.
    expect_length(intersect(train_sites, assess_sites), 0L)
  }
})


test_that("the assessment sets partition the rows exactly once", {
  d <- fb_sel_rs_data()
  folds <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)
  assessed <- unlist(lapply(folds, `[[`, "assess_rows"))
  expect_setequal(assessed, seq_len(nrow(d)))
  expect_identical(length(assessed), nrow(d))
})


test_that("indices index the ORIGINAL rows, in the original order", {
  d <- fb_sel_rs_data()
  folds <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)
  for (split in folds) {
    expect_identical(split$train_rows, sort(split$train_rows))
    expect_true(all(split$train_rows >= 1L & split$train_rows <= nrow(d)))
    expect_length(
      intersect(split$train_rows, split$assess_rows),
      0L
    )
  }
})


test_that("an explicit fold_id that splits a unit refuses", {
  # This is the leakage the `unit` argument exists to prevent, arriving by
  # the back door. It has to be caught at the same place.
  d <- fb_sel_rs_data()
  fid <- rep(1:2, length.out = nrow(d)) # alternates inside every site
  expect_error(
    .fb_sel_folds(fb_resamples(unit = ~site, fold_id = fid), d),
    "splits unit",
    fixed = TRUE
  )
})


test_that("a malformed explicit fold_id refuses by shape", {
  d <- fb_sel_rs_data()
  expect_error(
    .fb_sel_folds(fb_resamples(fold_id = rep(1:2, 3)), d),
    "entries and the data has",
    fixed = TRUE
  )
  fid <- rep(1:2, length.out = nrow(d))
  fid[3L] <- NA
  expect_error(
    .fb_sel_folds(fb_resamples(fold_id = fid), d),
    "missing value",
    fixed = TRUE
  )
  expect_error(
    .fb_sel_folds(fb_resamples(fold_id = rep(1L, nrow(d))), d),
    "there is no holdout",
    fixed = TRUE
  )
})


test_that("a fold variable absent from the data refuses, and names it", {
  d <- fb_sel_rs_data()
  expect_error(
    .fb_sel_folds(fb_resamples(unit = ~paddock, v = 3L), d),
    "absent from the data: paddock",
    fixed = TRUE
  )
})


test_that("a stratifying variable that varies inside a unit refuses", {
  d <- fb_sel_rs_data()
  d$block <- rep(c("a", "b"), length.out = nrow(d))
  expect_error(
    .fb_sel_folds(fb_resamples(unit = ~site, v = 3L, strata = ~block), d),
    "is not constant within unit",
    fixed = TRUE
  )
})


test_that("strata are balanced across folds", {
  d <- fb_sel_rs_data()
  # Six wet sites and six dry ones, constant within site.
  d$season <- ifelse(as.integer(d$site) <= 6L, "wet", "dry")
  folds <- .fb_sel_folds(
    fb_resamples(unit = ~site, v = 3L, strata = ~season),
    d
  )
  for (split in folds) {
    seasons <- unique(d$season[split$assess_rows])
    expect_setequal(seasons, c("wet", "dry"))
  }
})


test_that("a design with one unit refuses rather than score nothing", {
  d <- fb_sel_rs_data(n_cluster = 1L, per = 8L)
  expect_error(
    .fb_sel_folds(fb_resamples(unit = ~site, v = 2L), d),
    "training and",
    fixed = TRUE
  )
})


# --- ORACLE 4b: the seed rule ------------------------------------------ #

test_that("repeats without a seed refuses, and no seed is invented", {
  d <- fb_sel_rs_data()
  expect_error(
    .fb_sel_folds(fb_resamples(unit = ~site, v = 4L, repeats = 3L), d),
    "No seed is invented here",
    fixed = TRUE
  )
})


test_that("the same seed twice gives bit-identical fold assignments", {
  d <- fb_sel_rs_data()
  spec <- fb_resamples(unit = ~site, v = 4L, repeats = 3L)
  a <- .fb_sel_folds(spec, d, seed = 2026L)
  b <- .fb_sel_folds(spec, d, seed = 2026L)
  expect_identical(
    lapply(a, `[`, c("fold_id", "repeat_id", "train_rows", "assess_rows")),
    lapply(b, `[`, c("fold_id", "repeat_id", "train_rows", "assess_rows"))
  )
  # And a different seed is a different assignment, or the seed is decorative.
  c_folds <- .fb_sel_folds(spec, d, seed = 7L)
  expect_false(identical(
    lapply(a, `[[`, "assess_rows"),
    lapply(c_folds, `[[`, "assess_rows")
  ))
})


test_that("an unrepeated unseeded deal is systematic and touches no stream", {
  d <- fb_sel_rs_data()
  spec <- fb_resamples(unit = ~site, v = 4L)
  set.seed(1L)
  before <- get(".Random.seed", envir = globalenv())
  a <- .fb_sel_folds(spec, d)
  after <- get(".Random.seed", envir = globalenv())
  b <- .fb_sel_folds(spec, d)

  expect_identical(attr(a, "assignment"), "systematic")
  expect_identical(before, after)
  expect_identical(
    lapply(a, `[[`, "assess_rows"),
    lapply(b, `[[`, "assess_rows")
  )
})


test_that("a seeded deal puts the caller's random stream back", {
  d <- fb_sel_rs_data()
  set.seed(3L)
  before <- get(".Random.seed", envir = globalenv())
  .fb_sel_folds(fb_resamples(unit = ~site, v = 4L, repeats = 2L), d, seed = 9L)
  expect_identical(before, get(".Random.seed", envir = globalenv()))
})


# --- frozen levels and contrasts --------------------------------------- #

test_that("contrast columns are frozen across folds", {
  d <- fb_sel_rs_data()
  d$variety <- factor(rep(c("a", "b", "c", "d"), length.out = nrow(d)))
  folds <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)

  designs <- lapply(folds, function(split) {
    .fb_sel_fold_substrate(
      d,
      split,
      response = "y",
      protected = ~ x + variety,
      random = "site"
    )$x_protected
  })
  widths <- vapply(designs, ncol, integer(1L))
  names_all <- lapply(designs, colnames)

  expect_true(all(widths == widths[[1L]]))
  for (nm in names_all) {
    expect_identical(nm, names_all[[1L]])
  }
  expect_identical(widths[[1L]], 5L)
})


test_that("an assessment level absent from training is a named outcome", {
  # Contract section 3.3: never a silent drop. The refusal names the term and
  # the level, because "a level was missing" is not actionable on its own.
  d <- fb_sel_rs_data()
  d$variety <- factor(ifelse(as.integer(d$site) == 1L, "rare", "common"))
  fid <- ifelse(as.integer(d$site) == 1L, "1", "2")
  folds <- .fb_sel_folds(fb_resamples(unit = ~site, fold_id = fid), d)
  held_out_first <- folds[[which(vapply(
    folds,
    function(s) identical(s$fold_id, "1"),
    logical(1L)
  ))]]

  expect_error(
    .fb_sel_fold_substrate(
      d,
      held_out_first,
      response = "y",
      protected = ~ x + variety,
      random = "site"
    ),
    "absent from its training rows: rare",
    fixed = TRUE
  )
})


# --- ORACLE 4c: the leakage gate, and the defects that make it red ------ #

test_that("the module's own fold substrate passes the leakage gate", {
  d <- fb_sel_rs_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]
  report <- .fb_sel_leakage_gate(fb_sel_rs_builder, d, split)

  expect_true(report$clean)
  expect_length(report$leaked, 0L)
  # Nothing inert: every declared quantity moved when the TRAINING rows moved,
  # so the gate is looking at something rather than at a set of constants.
  expect_length(report$inert, 0L)
  expect_setequal(
    report$checked,
    c(
      "components_residual",
      "components_terms",
      "preprocessing_centre",
      "preprocessing_scale",
      "pilot_beta",
      "logdet_v_tt"
    )
  )
  expect_silent(.fb_sel_certify_no_leakage(fb_sel_rs_builder, d, split))
})


test_that("FALSIFIER: scaling learned on all rows turns the gate red", {
  # The single most important falsifier in this file. The defect is one
  # argument: the preprocessing is learned from every row instead of from the
  # training rows. Nothing else changes, no test of a coefficient or a score
  # would notice, and the gate must.
  d <- fb_sel_rs_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]

  leaky <- function(data, split) {
    substrate <- fb_sel_rs_builder(data, split)
    all_rows <- seq_len(nrow(data))
    candidates <- .fb_sel_design(~ z1 + z2 + z3 + z4, data)$X
    candidates <- candidates[,
      colnames(candidates) != "(Intercept)",
      drop = FALSE
    ]
    leaked <- .fb_sel_preprocessing(candidates, all_rows)
    substrate$learned$preprocessing_centre <- leaked$centre
    substrate$learned$preprocessing_scale <- leaked$scale
    substrate
  }

  report <- .fb_sel_leakage_gate(leaky, d, split)
  expect_false(report$clean)
  expect_setequal(
    report$leaked,
    c("preprocessing_centre", "preprocessing_scale")
  )
  expect_error(
    .fb_sel_certify_no_leakage(leaky, d, split),
    "were not learned on its training rows",
    fixed = TRUE
  )
})


test_that("FALSIFIER: a covariance learned on all rows turns the gate red", {
  d <- fb_sel_rs_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]

  leaky <- function(data, split) {
    substrate <- fb_sel_rs_builder(data, split)
    all_rows <- seq_len(nrow(data))
    leaked <- .fb_sel_reml_components(
      data$y,
      .fb_sel_design(~x, data)$X,
      .fb_sel_Z_from_vars(data, "site"),
      all_rows
    )
    substrate$learned$components_residual <- leaked$residual
    substrate$learned$components_terms <- leaked$terms
    substrate
  }

  report <- .fb_sel_leakage_gate(leaky, d, split)
  expect_false(report$clean)
  expect_setequal(
    report$leaked,
    c("components_residual", "components_terms")
  )
  expect_error(
    .fb_sel_certify_no_leakage(leaky, d, split),
    "components_residual",
    fixed = TRUE
  )
})


test_that("FALSIFIER: an inert declared quantity is reported, not passed", {
  # A gate that only asks "did it move with the assessment rows" stays green
  # on a hard-coded constant, which is how a gate comes to pass vacuously.
  d <- fb_sel_rs_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]

  inert <- function(data, split) {
    substrate <- fb_sel_rs_builder(data, split)
    substrate$learned$pilot_beta <- c(0, 0, 0, 0, 0)
    substrate
  }

  report <- .fb_sel_leakage_gate(inert, d, split)
  expect_false(report$clean)
  expect_length(report$leaked, 0L)
  expect_identical(report$inert, "pilot_beta")
  expect_error(
    .fb_sel_certify_no_leakage(inert, d, split),
    "the gate cannot see them",
    fixed = TRUE
  )
})


test_that("the gate refuses a builder that claims nothing", {
  d <- fb_sel_rs_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]
  expect_error(
    .fb_sel_leakage_gate(
      function(data, split) list(learned = list()),
      d,
      split
    ),
    "would pass on an empty claim",
    fixed = TRUE
  )
})


test_that("the REML plugin is bit-identical across identical calls", {
  # The gate compares with identical() and carries no tolerance at all, which
  # is only legitimate because the estimator is deterministic. Measured here
  # rather than assumed: five calls, spread exactly zero.
  d <- fb_sel_rs_data()
  Z <- .fb_sel_Z_from_vars(d, "site")
  X <- .fb_sel_design(~x, d)$X
  rows <- seq_len(nrow(d))
  runs <- replicate(
    5L,
    {
      out <- .fb_sel_reml_components(d$y, X, Z, rows)
      c(out$residual, out$terms)
    },
    simplify = "array"
  )
  expect_identical(max(apply(runs, 1L, function(v) diff(range(v)))), 0)
})


# --- ORACLE 4d: the rotation order (contract section 4.2) --------------- #

test_that("splitting rotated rows gives an optimistic score, and by how much", {
  # A rotated row is a contrast among the rows of its covariance block, not an
  # observation, so the index set that names it does not name a unit. Fitting
  # on the remaining rotated rows and mapping the residual back to the
  # original scale, e = L r*, gives the prediction error the rotated-space fit
  # actually implies -- and it is small, because most of what it left to
  # predict had already been reconstructed from cluster mates that stayed in
  # training.
  #
  # Measured on this fixture: correct 2.5911, rotated 1.2168, ratio 0.4696.
  # Across ten independent draws of the same design the ratio ran 0.30 to
  # 0.64, so the assertion below sits inside the observed range by a margin
  # and it is the direction, not the value, that is being fixed. The
  # assertion is not satisfied by any two numbers: two correctly-ordered
  # routes differing only in their fold draw give a ratio of 1.06 on this
  # same fixture, and it fails on them.
  d <- fb_sel_rs_data(seed = 21L)
  design <- model.matrix(~ x + z1 + z2 + z3 + z4, d)
  Z <- .fb_sel_Z_from_vars(d, "site")
  n <- nrow(d)

  correct <- .fb_sel_cv(
    d,
    response = "y",
    protected = ~x,
    candidates = ~ z1 + z2 + z3 + z4,
    random = "site",
    resamples = fb_resamples(unit = ~site, v = 4L),
    prediction = "new_cluster",
    cluster = "site"
  )
  rmse_correct <- correct$aggregated$score[
    correct$aggregated$model == "full"
  ]

  # The wrong order: covariance on every row, rotate, then split.
  components <- .fb_sel_reml_components(
    d$y,
    .fb_sel_design(~x, d)$X,
    Z,
    seq_len(n)
  )
  V <- diag(components$residual, n) +
    components$terms[["site"]] * tcrossprod(Z$site)
  L <- t(chol(V))
  y_rotated <- forwardsolve(L, d$y)
  x_rotated <- forwardsolve(L, design)
  fold <- .fb_sel_with_seed(5L, sample(rep(seq_len(4L), length.out = n)))
  rotated_residual <- numeric(n)
  for (f in seq_len(4L)) {
    assess <- which(fold == f)
    train <- which(fold != f)
    beta <- qr.coef(qr(x_rotated[train, , drop = FALSE]), y_rotated[train])
    beta[is.na(beta)] <- 0
    rotated_residual[train] <- y_rotated[train] -
      drop(x_rotated[train, , drop = FALSE] %*% beta)
    rotated_residual[assess] <- y_rotated[assess] -
      drop(x_rotated[assess, , drop = FALSE] %*% beta)
  }
  implied <- drop(L %*% rotated_residual)
  rmse_rotated <- mean(vapply(
    split(implied, d$site),
    function(r) sqrt(mean(r^2)),
    numeric(1L)
  ))

  expect_lt(rmse_rotated, 0.75 * rmse_correct)
  expect_gt(rmse_correct, 1)
})


# --- ORACLE 7: the prediction estimand ---------------------------------- #

# A dense case with two crossed groupings, small enough to condition by brute
# force. No model is fitted: the coefficients are given, so what is under test
# is the conditioning algebra and nothing else.
fb_sel_rs_dense <- function(seed = 7L) {
  .fb_sel_with_seed(seed, {
    n <- 40L
    d <- data.frame(
      g1 = factor(rep(1:8, each = 5)),
      g2 = factor(rep(1:5, times = 8)),
      x = stats::rnorm(n)
    )
    Z <- .fb_sel_Z_from_vars(d, c("g1", "g2"))
    residual <- 0.37
    terms <- c(g1 = 1.3, g2 = 0.6)
    V <- diag(residual, n)
    for (nm in names(Z)) {
      V <- V + terms[[nm]] * tcrossprod(Z[[nm]])
    }
    X <- model.matrix(~x, d)
    beta <- c(1.2, -0.8)
    y <- drop(X %*% beta) + drop(t(chol(V)) %*% stats::rnorm(n))
    list(
      data = d,
      Z = Z,
      V = V,
      X = X,
      beta = beta,
      y = y,
      cov_spec = .fb_sel_cov_spec(Z, residual, terms),
      train = sort(sample(n, 28L)),
      assess = NULL
    )
  })
}


test_that("oracle 7: the conditional prediction equals dense MVN algebra", {
  case <- fb_sel_rs_dense()
  train <- case$train
  assess <- setdiff(seq_len(nrow(case$data)), train)

  blocks <- .fb_sel_fold_blocks(case$cov_spec, train, assess, "known_cluster")
  mu <- .fb_sel_fold_predict(
    blocks,
    case$beta,
    case$X[train, , drop = FALSE],
    case$y[train],
    case$X[assess, , drop = FALSE]
  )

  v_tt <- case$V[train, train]
  v_st <- case$V[assess, train]
  v_ss <- case$V[assess, assess]
  mu_reference <- drop(case$X[assess, ] %*% case$beta) +
    drop(
      v_st %*%
        solve(
          v_tt,
          case$y[train] -
            drop(case$X[train, ] %*% case$beta)
        )
    )
  conditional_reference <- v_ss - v_st %*% solve(v_tt, t(v_st))

  # Measured on this case: 2.36e-14 relative on the mean and 1.72e-15 on the
  # covariance, against the contract's bar of relative 1e-8 for deterministic
  # linear algebra. The two routes share no code: the reference conditions
  # with solve(), the module with Cholesky solves.
  expect_lt(max(abs(mu - mu_reference) / abs(mu_reference)), 1e-8)
  expect_lt(
    max(abs(blocks$v_s_given_t - conditional_reference)) /
      max(abs(conditional_reference)),
    1e-8
  )
  expect_identical(max(abs(blocks$v_tt - v_tt)), 0)
})


test_that("oracle 7: a genuinely new cluster has V_ST = 0", {
  case <- fb_sel_rs_dense()
  train <- case$train
  assess <- setdiff(seq_len(nrow(case$data)), train)

  blocks <- .fb_sel_fold_blocks(case$cov_spec, train, assess, "new_cluster")
  expect_true(all(blocks$v_st == 0))
  expect_true(all(blocks$gain == 0))
  # And with V_ST zero the conditional covariance is the marginal one, so the
  # fold is scored on the fixed part alone.
  expect_identical(max(abs(blocks$v_s_given_t - case$V[assess, assess])), 0)

  mu <- .fb_sel_fold_predict(
    blocks,
    case$beta,
    case$X[train, , drop = FALSE],
    case$y[train],
    case$X[assess, , drop = FALSE]
  )
  expect_identical(mu, drop(case$X[assess, ] %*% case$beta))
})


test_that("oracle 7: folds on units make the cross-covariance exactly zero", {
  # The structural cross-covariance, not the estimand's zeroing of it. When
  # the fold boundary follows the grouping factor, the held-out rows share no
  # random-effect level with the training rows and V_ST is zero because of
  # the design, which is what "genuinely new cluster" means.
  d <- fb_sel_rs_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]
  Z <- .fb_sel_Z_from_vars(d, "site")
  spec <- .fb_sel_cov_spec(Z, 0.1, c(site = 4))
  blocks <- .fb_sel_fold_blocks(
    spec,
    split$train_rows,
    split$assess_rows,
    "known_cluster"
  )
  expect_true(all(blocks$v_st_structural == 0))
})


test_that("oracle 7: the fold blocks agree with the dense reference", {
  # Cross-check against the module's existing oracle rather than against a
  # second copy of the same arithmetic. `.fb_sel_V_dense()` assembles V over
  # every row, elementwise and deliberately naively; the fold blocks are the
  # same operator restricted to the training rows, the assessment rows and the
  # cross of the two, and if they disagree one of the two is wrong.
  #
  # The one fit in this file. john.alpha is the nested alpha design the
  # covariance oracles use, so the two files exercise the same base.
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  d <- agridat::john.alpha
  d$gen <- factor(d$gen)
  d$rep <- factor(d$rep)
  d$block <- factor(d$block)
  fit <- suppressMessages(flexybayes(
    fixed = yield ~ gen,
    random = ~ rep + rep:block,
    data = d,
    backend = "inla",
    verbose = FALSE
  ))

  reference <- .fb_sel_V_dense(fit, "plugin_mean_covariance")
  components <- .fb_sel_variance_components(fit, "plugin_mean_covariance")
  spec <- .fb_sel_cov_spec(
    .fb_sel_random_designs(fit),
    components$residual,
    components$terms,
    components$estimand
  )

  split <- .fb_sel_folds(fb_resamples(unit = ~rep, v = 3L), fit$data)[[1L]]
  blocks <- .fb_sel_fold_blocks(
    spec,
    split$train_rows,
    split$assess_rows,
    "known_cluster"
  )

  expect_equal(
    blocks$v_tt,
    reference[split$train_rows, split$train_rows],
    tolerance = 1e-12
  )
  expect_equal(
    blocks$v_ss,
    reference[split$assess_rows, split$assess_rows],
    tolerance = 1e-12
  )
  expect_equal(
    blocks$v_st_structural,
    reference[split$assess_rows, split$train_rows],
    tolerance = 1e-12
  )
  # Folds on the replicate: the held-out replicate shares no level of either
  # random term with the training rows, so the cross-block is exactly zero
  # and the base is a genuinely new cluster.
  expect_identical(max(abs(blocks$v_st_structural)), 0)
})


test_that("oracle 7: V_TT is factored and V^-1 is never formed", {
  case <- fb_sel_rs_dense()
  train <- case$train
  assess <- setdiff(seq_len(nrow(case$data)), train)
  blocks <- .fb_sel_fold_blocks(case$cov_spec, train, assess, "known_cluster")

  # Structural, not aspirational. The object holds an upper-triangular factor
  # of V_TT and an S-by-T gain, and nothing that is or contains an inverse.
  expect_true(all(blocks$factor_tt[lower.tri(blocks$factor_tt)] == 0))
  expect_equal(crossprod(blocks$factor_tt), blocks$v_tt, tolerance = 1e-12)
  expect_identical(dim(blocks$gain), c(length(assess), length(train)))
  expect_false(any(vapply(
    blocks,
    function(part) {
      is.matrix(part) &&
        identical(dim(part), dim(blocks$v_tt)) &&
        isTRUE(all.equal(part, solve(blocks$v_tt), tolerance = 1e-6))
    },
    logical(1L)
  )))
})


test_that("oracle 7: new_cluster scores worse than known_cluster, materially", {
  # The direction is the whole argument for the default. Conditioning on a
  # held-out location's own random effect hands the model the thing it was
  # asked to do without, and the score it earns is not the score of a
  # prediction. Row-level folds here, so every assessment row's cluster is
  # present in training and the two estimands genuinely differ.
  #
  # Measured on this fixture, for the full model: known 0.3100, new 1.0692, a
  # ratio of 0.290 on RMSE; known 0.3946 against new 1.3792 on the negative
  # log predictive density. Across eight independent draws the RMSE ratio
  # stayed between 0.176 and 0.278.
  d <- fb_sel_rs_data()
  spec <- fb_resamples(v = 4L)
  known <- .fb_sel_cv(
    d,
    "y",
    ~x,
    ~ z1 + z2 + z3 + z4,
    random = "site",
    resamples = spec,
    prediction = "known_cluster",
    cluster = "site"
  )
  new <- .fb_sel_cv(
    d,
    "y",
    ~x,
    ~ z1 + z2 + z3 + z4,
    random = "site",
    resamples = spec,
    prediction = "new_cluster",
    cluster = "site"
  )
  index <- which(known$aggregated$model == "full")
  expect_lt(
    known$aggregated$score[index],
    0.5 * new$aggregated$score[index]
  )

  known_nlpd <- .fb_sel_cv(
    d,
    "y",
    ~x,
    ~ z1 + z2 + z3 + z4,
    random = "site",
    resamples = spec,
    prediction = "known_cluster",
    cluster = "site",
    metric = "nlpd"
  )
  new_nlpd <- .fb_sel_cv(
    d,
    "y",
    ~x,
    ~ z1 + z2 + z3 + z4,
    random = "site",
    resamples = spec,
    prediction = "new_cluster",
    cluster = "site",
    metric = "nlpd"
  )
  expect_lt(
    known_nlpd$aggregated$score[index],
    new_nlpd$aggregated$score[index] - 0.5
  )
})


test_that("a training covariance with no Cholesky factor refuses by name", {
  # A variance component reported negative by an engine is the realistic way
  # this happens. No jitter, no repair: the fold refuses.
  d <- fb_sel_rs_data()
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]
  spec <- .fb_sel_cov_spec(.fb_sel_Z_from_vars(d, "site"), 0.1, c(site = -5))
  expect_error(
    .fb_sel_fold_blocks(spec, split$train_rows, split$assess_rows),
    "No jitter was applied",
    fixed = TRUE
  )
})


# --- ORACLE 9: the scoring unit ----------------------------------------- #

test_that("oracle 9: the reported standard error is the cluster one", {
  # Same point estimate, different standard error. The per-row squared errors
  # inside a site share that site's random effect, so treating them as
  # independent evidence divides the standard error by a factor the data never
  # earned.
  #
  # Measured on this fixture: point estimate 6.509953394 both ways, agreeing
  # to 8.9e-16; cluster standard error 3.5308 against a row-level 1.4072, a
  # ratio of 2.509. Across eight independent draws the ratio ran 2.26 to 2.53.
  d <- fb_sel_rs_data()
  cv <- .fb_sel_cv(
    d,
    "y",
    ~x,
    ~ z1 + z2 + z3 + z4,
    random = "site",
    resamples = fb_resamples(unit = ~site, v = 4L),
    prediction = "new_cluster",
    metric = "mse"
  )
  clusters <- cv$by_model$full
  rows <- cv$rows
  keep <- rows$model == "full"

  per_cluster <- vapply(
    split(clusters$mse, clusters$cluster),
    mean,
    numeric(1L)
  )
  se_cluster <- stats::sd(per_cluster) / sqrt(length(per_cluster))
  se_row <- stats::sd(rows$squared_error[keep]) / sqrt(sum(keep))

  # The two are estimating the standard error of the SAME number.
  expect_equal(
    mean(per_cluster),
    mean(rows$squared_error[keep]),
    tolerance = 1e-12
  )
  expect_gt(se_cluster / se_row, 2)

  reported <- cv$aggregated
  index <- which(reported$model == "full")
  expect_identical(reported$se_unit[[index]], "cluster")
  expect_equal(reported$se[[index]], se_cluster, tolerance = 1e-12)
  expect_identical(reported$n_units[[index]], 12L)
})


test_that("oracle 9: a cluster measured again is not new evidence", {
  # With two repeats every site is scored twice. The independent unit is still
  # the site, so the repeats average within the site before the standard error
  # is taken and the unit count does not double.
  d <- fb_sel_rs_data()
  cv <- .fb_sel_cv(
    d,
    "y",
    ~x,
    ~ z1 + z2 + z3 + z4,
    random = "site",
    resamples = fb_resamples(unit = ~site, v = 4L, repeats = 2L),
    seed = 101L,
    prediction = "new_cluster"
  )
  expect_identical(nrow(cv$by_model$full), 24L)
  expect_true(all(cv$aggregated$n_units == 12L))
})


test_that("oracle 9: the cluster density is the joint one, not the rows'", {
  # The negative log predictive density of a correlated block is not the sum
  # of its rows' marginal densities. A common site effect explains every row
  # of the site at once, so the joint density is the higher one, and it is
  # the density of what was actually observed.
  #
  # Measured on this fixture: joint 11.8489 against a row-marginal 29.3951, a
  # ratio of 0.4031. Across six independent draws the ratio ran 0.19 to 0.40.
  d <- fb_sel_rs_data(seed = 401L)
  split <- .fb_sel_folds(fb_resamples(unit = ~site, v = 4L), d)[[1L]]
  substrate <- fb_sel_rs_builder(d, split)
  design <- cbind(substrate$x_protected, substrate$x_candidate)
  fit <- .fb_sel_gls_fit(
    substrate$blocks$factor_tt,
    design[split$train_rows, , drop = FALSE],
    substrate$response[split$train_rows]
  )
  mu <- .fb_sel_fold_predict(
    substrate$blocks,
    fit$beta,
    design[split$train_rows, , drop = FALSE],
    substrate$response[split$train_rows],
    design[split$assess_rows, , drop = FALSE]
  )
  scored <- .fb_sel_score_fold(
    substrate$assessment$y,
    mu,
    substrate$assessment$v_s_given_t,
    as.character(d$site)[split$assess_rows]
  )
  joint <- sum(scored$clusters$nlpd * scored$clusters$n)
  marginal <- sum(scored$rows$nlpd_row)
  expect_lt(joint, 0.6 * marginal)

  # Where the block IS diagonal the two must agree exactly, or one of them is
  # not a Gaussian density.
  independent <- .fb_sel_score_fold(
    substrate$assessment$y,
    mu,
    diag(diag(substrate$assessment$v_s_given_t)),
    as.character(d$site)[split$assess_rows]
  )
  expect_equal(
    sum(independent$clusters$nlpd * independent$clusters$n),
    sum(independent$rows$nlpd_row),
    tolerance = 1e-10
  )
})


test_that("a scoring unit with fewer than two clusters refuses", {
  d <- fb_sel_rs_data()
  split <- .fb_sel_folds(fb_resamples(v = 2L), d)[[1L]]
  substrate <- fb_sel_rs_builder(d, split)
  design <- cbind(substrate$x_protected, substrate$x_candidate)
  fit <- .fb_sel_gls_fit(
    substrate$blocks$factor_tt,
    design[split$train_rows, , drop = FALSE],
    substrate$response[split$train_rows]
  )
  mu <- .fb_sel_fold_predict(
    substrate$blocks,
    fit$beta,
    design[split$train_rows, , drop = FALSE],
    substrate$response[split$train_rows],
    design[split$assess_rows, , drop = FALSE]
  )
  scored <- .fb_sel_score_fold(
    substrate$assessment$y,
    mu,
    substrate$assessment$v_s_given_t,
    rep("one_and_only", length(split$assess_rows))
  )
  expect_error(
    .fb_sel_aggregate_scores(scored$clusters),
    "A standard error needs at least two",
    fixed = TRUE
  )
})


# --- ORACLE 9b: the decision rule --------------------------------------- #

test_that("one-standard-error selection uses the error at the minimum", {
  # A hand-built path, so the rule is under test and not the fit that produced
  # it. The minimum sits at the most complex model; a much simpler one is
  # inside one standard error of it, and that is the one the rule takes.
  path <- data.frame(
    model = c("m0", "m1", "m2", "m3"),
    score = c(1.40, 1.02, 1.05, 1.00),
    se = c(0.20, 0.10, 0.10, 0.10),
    complexity = c(0L, 1L, 2L, 3L)
  )
  expect_identical(.fb_sel_select_model(path, "min")$selected, "m3")
  expect_identical(.fb_sel_select_model(path, "one_se")$selected, "m1")
  expect_equal(.fb_sel_select_model(path, "one_se")$threshold, 1.10)

  # With a standard error that is too small -- which is exactly what a
  # row-level unit gives on correlated data -- the rule collapses to the
  # minimum and stops being a rule.
  tight <- path
  tight$se <- tight$se / 10
  expect_identical(.fb_sel_select_model(tight, "one_se")$selected, "m3")
})


test_that("one-standard-error selection reproduces across reruns", {
  d <- fb_sel_rs_data()
  sets <- list(m0 = integer(0), m1 = 1L, m2 = 1:2, m4 = 1:4)
  run <- function() {
    .fb_sel_cv(
      d,
      "y",
      ~x,
      ~ z1 + z2 + z3 + z4,
      random = "site",
      resamples = fb_resamples(unit = ~site, v = 4L, repeats = 2L),
      seed = 505L,
      prediction = "new_cluster",
      candidate_sets = sets,
      cluster = "site"
    )
  }
  first <- run()
  second <- run()

  expect_identical(first$aggregated$model, second$aggregated$model)
  expect_identical(first$aggregated$score, second$aggregated$score)
  expect_identical(first$aggregated$se, second$aggregated$se)
  expect_identical(
    .fb_sel_select_model(first$aggregated, "one_se"),
    .fb_sel_select_model(second$aggregated, "one_se")
  )
  # A different seed is a different set of folds, or the seed is decorative.
  different <- .fb_sel_cv(
    d,
    "y",
    ~x,
    ~ z1 + z2 + z3 + z4,
    random = "site",
    resamples = fb_resamples(unit = ~site, v = 4L, repeats = 2L),
    seed = 606L,
    prediction = "new_cluster",
    candidate_sets = sets,
    cluster = "site"
  )
  expect_false(identical(first$aggregated$score, different$aggregated$score))
})


test_that("the protected-only baseline is scored alongside the full model", {
  # Contract section 11: the increment falls out of the substrate rather than
  # needing a second call, because the protected model is always one of the
  # scored candidate sets.
  d <- fb_sel_rs_data()
  cv <- .fb_sel_cv(
    d,
    "y",
    ~x,
    ~ z1 + z2 + z3 + z4,
    random = "site",
    resamples = fb_resamples(unit = ~site, v = 4L),
    prediction = "new_cluster"
  )
  expect_setequal(cv$aggregated$model, c("protected_only", "full"))
  expect_identical(
    cv$aggregated$complexity[cv$aggregated$model == "protected_only"],
    0L
  )
  expect_identical(cv$aggregated$complexity[cv$aggregated$model == "full"], 4L)
  # Four site-constant covariates with no true effect, scored on held-out
  # sites: the increment is negative, and the substrate says so.
  expect_gt(
    cv$aggregated$score[cv$aggregated$model == "full"],
    cv$aggregated$score[cv$aggregated$model == "protected_only"]
  )
})


# --- refusals ------------------------------------------------------------ #

test_that("every refusal this module can throw is registered", {
  # The package raises a contract violation when an unregistered code is
  # thrown, which means an unregistered code surfaces as the wrong error and a
  # bare expect_error() passes on it. This test is the direct check, and every
  # code below has a test above that fires it.
  codes <- c(
    "selection_seed_required",
    "selection_fold_level_absent",
    "selection_resamples_underspecified",
    "selection_fold_id_invalid",
    "selection_fold_variable_not_in_data",
    "selection_strata_not_constant_within_unit",
    "selection_fold_not_usable",
    "selection_fold_covariance_not_factorisable",
    "selection_leakage_detected",
    "selection_scoring_unit_undeclared"
  )
  expect_true(all(codes %in% fb_refusals()$reason_code))
})


test_that("the refusals carry their reason code as a condition class", {
  d <- fb_sel_rs_data()
  condition <- tryCatch(
    .fb_sel_folds(fb_resamples(unit = ~site, v = 4L, repeats = 2L), d),
    flexybayes_refusal_selection_seed_required = function(e) e
  )
  expect_s3_class(condition, "flexybayes_refusal_selection_seed_required")
})
