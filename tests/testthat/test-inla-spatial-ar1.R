# =============================================================================
# Separable AR1(row) x AR1(col) spatial latent field on INLA.
#
# A designed field trial's spatial structure is written on the RANDOM side --
# random = ~ ar1(row):ar1(col) -- and emitted on INLA as the grouped-AR1
# latent field
#   f(<row>_id, model = "ar1", group = <col>_id,
#     control.group = list(model = "ar1"))
# plus the Gaussian observation nugget. Faithful under one observation per
# grid node; an incomplete or replicated grid refuses.
#
# 0.9.0 respelling. The same field used to be reachable through the RESIDUAL
# spelling, which is ASReml's name for a three-parameter nugget-free process
# -- a different model from the four-parameter field-plus-nugget the emit
# builds. The residual spelling is now refused by name and the refusal points
# at the random-side spelling; the pinning test is at the foot of this file.
#
# The FAITHFULNESS clause (OO-spatial) is an OPEN oracle in the I4 sense: INLA
# is validated against an INDEPENDENT hand-built GLS/REML implementation of the
# SAME AR1(row)xAR1(col) model (Kronecker-eigenbasis profile likelihood) that
# does NOT call INLA -- so agreement means the emit is a faithful
# implementation, not self-consistency. Mirrors rebuild/wp16_inla_spatial_oracle.R.
# =============================================================================

.ar1_mat <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))

# One observation per node on a complete n_row x n_col grid, with a spatial
# field ~ N(0, sd_s^2 (R_row(rho_r) (x) R_col(rho_c))) + nugget ~ N(0, sd_e^2).
# Row is the OUTER index, col the INNER index, matching kronecker(R_row, R_col).
.spatial_sim <- function(n_row, n_col, rho_r, rho_c, sd_s, sd_e, mu = 5) {
  grid <- expand.grid(col = seq_len(n_col), row = seq_len(n_row)) # col fast
  L <- t(chol(sd_s^2 * kronecker(.ar1_mat(n_row, rho_r), .ar1_mat(n_col, rho_c))))
  u <- as.numeric(L %*% stats::rnorm(n_row * n_col))
  data.frame(
    y = mu + u + stats::rnorm(n_row * n_col, 0, sd_e),
    row = factor(grid$row), col = factor(grid$col)
  )
}

# INDEPENDENT open oracle: REML for y ~ 1, V = ss2 (Rr (x) Rc) + se2 I, solved
# in the Kronecker eigenbasis. Does NOT call INLA. Returns the ML/REML point
# estimates of (rho_row, rho_col, sd_s, sd_e).
.spatial_reml_oracle <- function(y, n_row, n_col) {
  X <- matrix(1, length(y), 1L)
  vy <- stats::var(y)
  nll <- function(par) {
    rr <- tanh(par[1]); rc <- tanh(par[2]); ss2 <- exp(par[3]); se2 <- exp(par[4])
    er <- eigen(.ar1_mat(n_row, rr), symmetric = TRUE)
    ec <- eigen(.ar1_mat(n_col, rc), symmetric = TRUE)
    U <- kronecker(er$vectors, ec$vectors)
    dv <- as.numeric(ss2 * kronecker(er$values, ec$values) + se2)
    if (any(dv <= 0)) return(1e10)
    yt <- crossprod(U, y); Xt <- crossprod(U, X)
    A <- crossprod(Xt, Xt / dv)
    beta <- solve(A, crossprod(Xt, yt / dv))
    rt <- yt - Xt %*% beta
    0.5 * (sum(log(dv)) + sum((rt^2) / dv) +
             determinant(A, logarithm = TRUE)$modulus[1])
  }
  best <- NULL
  for (s in list(c(atanh(0.3), atanh(0.3), log(vy / 2), log(vy / 2)),
                 c(atanh(0.6), atanh(0.2), log(vy / 2), log(vy / 4)))) {
    o <- tryCatch(stats::optim(s, nll, method = "Nelder-Mead",
                               control = list(maxit = 800L, reltol = 1e-10)),
                  error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  p <- best$par
  list(rho_row = tanh(p[1]), rho_col = tanh(p[2]),
       sd_s = sqrt(exp(p[3])), sd_e = sqrt(exp(p[4])))
}

.inla_spatial_hyper <- function(fit) {
  hp <- fit$inla$summary.hyperpar
  gv <- function(p) {
    i <- grep(p, rownames(hp)); if (length(i)) hp$mean[i[1]] else NA_real_
  }
  list(
    rho_row = gv("^Rho for row_id"),
    rho_col = gv("GroupRho for row_id"),
    sd_s = 1 / sqrt(gv("^Precision for row_id")),
    sd_e = 1 / sqrt(gv("Gaussian observations"))
  )
}

# ---- emit: the grouped-AR1 f() string ---------------------------------------

test_that("emit: separable ar1(row):ar1(col) -> INLA grouped-AR1 f()", {
  term <- list(type = "ar1_spatial", row_var = "row", col_var = "col",
               col_ar1 = TRUE)
  expect_identical(
    flexyBayes:::.inla_ar1_field_term(term),
    paste0("f(row_id, model = \"ar1\", group = col_id, ",
           "control.group = list(model = \"ar1\"))")
  )
})

test_that("emit: ar1(row):id(col) -> AR1 x iid group model", {
  term <- list(type = "ar1_spatial", row_var = "row", col_var = "col",
               col_ar1 = FALSE)
  expect_match(
    flexyBayes:::.inla_ar1_field_term(term),
    "control.group = list(model = \"iid\")", fixed = TRUE
  )
})

test_that("emit: 1D ar1(t) -> f(t_id, model = 'ar1')", {
  expect_identical(
    flexyBayes:::.inla_ar1_field_term(list(type = "ar1", var = "t")),
    "f(t_id, model = \"ar1\")"
  )
})

# ---- gate: ar1 / ar1_spatial accepted on INLA as a RANDOM term --------------

test_that("gate: a random-side ar1_spatial is accepted on INLA", {
  d <- .spatial_sim(6L, 5L, 0.7, 0.4, 1, 0.5)
  fb <- fb_from_asreml(y ~ 1, random = ~ ar1(row):ar1(col), data = d)
  expect_true(flexyBayes:::.lgm_check_random_term_inla_support(fb)$pass)
  expect_true(flexyBayes:::.lgm_check_residual_term_inla_support(fb)$pass)
})

test_that("gate: the INLA residual allowlist is the homogeneous form only", {
  # The allowlist is the gate's copy of what the emit builds. With the field
  # on the random side, a structured residual has no INLA emit at all, and
  # the allowlist must say so or the gate will accept what emit refuses.
  expect_identical(flexyBayes:::.inla_residual_term_type_allowlist(), "units")
})

test_that("gate: the hyperparameter budget counts the field's own parameters", {
  # An under-count buys a fit that should have been refused: the budget is
  # what decides whether INLA's numerical integration stays tractable. The
  # 1-D field carries a precision and a correlation; the separable form adds
  # the group correlation; the observation nugget is the likelihood's own
  # precision, counted on the family side.
  fb1 <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(),
    random_terms = list(list(type = "ar1", var = "t")),
    residual_terms = list(list(type = "units")), source = "asreml"
  )
  expect_identical(flexyBayes:::.lgm_count_hyperparams(fb1), 3L)
  fb2 <- flexyBayes:::new_fb_terms(
    response = "y", family = "gaussian", link = "identity",
    fixed_terms = list(),
    random_terms = list(list(
      type = "ar1_spatial", row_var = "row", col_var = "col", col_ar1 = TRUE
    )),
    residual_terms = list(list(type = "units")), source = "asreml"
  )
  expect_identical(flexyBayes:::.lgm_count_hyperparams(fb2), 4L)
})

# ---- fit + OPEN-ORACLE faithfulness -----------------------------------------

test_that("INLA grouped-AR1 recovers the separable structure vs an independent REML oracle", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_convergence_warning = TRUE
  )
  set.seed(20260725L)
  # Distinct row/col correlations make recovery a defect-catcher: an iid emit,
  # a single-AR1, or a swapped-Kronecker orientation cannot match both.
  n_row <- 12L; n_col <- 10L
  rho_r <- 0.8; rho_c <- 0.4; sd_s <- 1.0; sd_e <- 0.5
  d <- .spatial_sim(n_row, n_col, rho_r, rho_c, sd_s, sd_e)

  fit <- flexybayes(y ~ 1, random = ~ ar1(row):ar1(col), data = d,
                    backend = "inla")
  inla_h <- .inla_spatial_hyper(fit)
  orc <- .spatial_reml_oracle(d$y, n_row, n_col)

  # Orientation: two DISTINCT correlations in the right order (row > col).
  expect_gt(inla_h$rho_row, inla_h$rho_col)
  # INLA agrees with the INDEPENDENT REML oracle (the faithfulness claim). The
  # tolerance accommodates INLA's default-prior shrinkage + single-realisation
  # noise (the spike measured |INLA - oracle| ~ 0.09 at this grid size).
  expect_lt(abs(inla_h$rho_row - orc$rho_row), 0.15)
  expect_lt(abs(inla_h$rho_col - orc$rho_col), 0.15)
  expect_lt(abs(inla_h$sd_s - orc$sd_s), 0.25)
  # Both correlations recovered well away from zero (not collapsed to iid).
  expect_gt(inla_h$rho_row, 0.4)
  expect_gt(inla_h$rho_col, 0.1)
})

test_that("summary and print report all four field parameters", {
  # The four parameters ARE the difference between what is fitted and the
  # ASReml residual it used to be spelled as, so a reader must be able to
  # see all of them without translating INLA's precision-scale names. The
  # numbers are checked against the simulation, not merely counted.
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(11L)
  d <- .spatial_sim(11L, 22L, 0.7, 0.5, 1, 0.5)
  fit <- flexybayes(y ~ 1, random = ~ ar1(row):ar1(col), data = d,
                    backend = "inla")

  tab <- flexyBayes:::.inla_spatial_hyper_table(fit)
  expect_equal(nrow(tab), 4L)
  expect_identical(
    tab$parameter,
    c(
      "correlation along row (rho_row)",
      "correlation along col (rho_col)",
      "field SD (sigma_field)", "nugget SD (sigma_e)"
    )
  )
  # Posterior medians near the values the grid was simulated from, and every
  # interval ordered.
  expect_lt(abs(tab$median[[1L]] - 0.7), 0.2)
  expect_lt(abs(tab$median[[2L]] - 0.5), 0.2)
  expect_lt(abs(tab$median[[3L]] - 1.0), 0.3)
  expect_lt(abs(tab$median[[4L]] - 0.5), 0.3)
  expect_true(all(tab$lower <= tab$median & tab$median <= tab$upper))

  for (txt in list(
    utils::capture.output(print(fit)),
    utils::capture.output(summary(fit))
  )) {
    joined <- paste(txt, collapse = "\n")
    expect_match(joined, "correlation along row (rho_row)", fixed = TRUE)
    expect_match(joined, "correlation along col (rho_col)", fixed = TRUE)
    expect_match(joined, "field SD (sigma_field)", fixed = TRUE)
    expect_match(joined, "nugget SD (sigma_e)", fixed = TRUE)
  }
})

test_that("a 1-D field reports its correlation, field SD and nugget SD", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(4L)
  d1 <- data.frame(t = factor(seq_len(60L)))
  d1$y <- as.numeric(stats::arima.sim(list(ar = 0.7), 60L)) +
    stats::rnorm(60L, 0, 0.3)
  fit <- flexybayes(y ~ 1, random = ~ ar1(t), data = d1, backend = "inla")
  tab <- flexyBayes:::.inla_spatial_hyper_table(fit)
  expect_equal(nrow(tab), 3L)
  expect_match(tab$parameter[[1L]], "correlation along t")
  expect_match(tab$parameter[[3L]], "nugget SD")
})

# ---- canonical names --------------------------------------------------------

test_that("canonical names: rho_row / rho_col / sd_spatial resolve", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(3L)
  d <- .spatial_sim(7L, 6L, 0.7, 0.4, 1, 0.5)
  fit <- flexybayes(y ~ 1, random = ~ ar1(row):ar1(col), data = d,
                    backend = "inla")
  co <- flexyBayes:::.fb_canonical_draws(fit)
  expect_true(all(c("rho_row", "rho_col", "sd_spatial") %in% names(co)))
})

# ---- routing + refusal ------------------------------------------------------

test_that("auto routes a designed spatial trial to INLA", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(5L)
  d <- .spatial_sim(6L, 5L, 0.7, 0.4, 1, 0.5)
  fit <- flexybayes(y ~ 1, random = ~ ar1(row):ar1(col), data = d,
                    backend = "auto")
  expect_identical(fit$extras$backend_decision$backend, "inla")
})

test_that("the residual spelling refuses and names the random-side field", {
  # D-A. ASReml's residual ar1(row):ar1(col) is a three-parameter nugget-free
  # process; the INLA emit builds a latent field plus the observation nugget,
  # four parameters. The package used to fit the second under the first's
  # name. It now refuses, on every engine, and says where the field lives.
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(5L)
  d <- .spatial_sim(6L, 5L, 0.7, 0.4, 1, 0.5)
  refusal <- tryCatch(
    suppressMessages(flexybayes(
      y ~ 1, residual = ~ ar1(row):ar1(col), data = d, backend = "inla"
    )),
    error = function(e) e
  )
  expect_s3_class(
    refusal, "flexybayes_refusal_ar1_residual_not_representable"
  )
  msg <- conditionMessage(refusal)
  expect_match(msg, "random = ~ ar1(row):ar1(col)", fixed = TRUE)
  expect_match(msg, "four parameters", fixed = TRUE)
  expect_match(msg, "ASReml", fixed = TRUE)
  # The 1-D residual spelling refuses on the same grounds.
  expect_error(
    suppressMessages(flexybayes(
      y ~ 1, residual = ~ ar1(row), data = d, backend = "inla"
    )),
    class = "flexybayes_refusal_ar1_residual_not_representable"
  )
})

test_that("a replicated / incomplete grid refuses (fidelity condition)", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(7L)
  d <- .spatial_sim(6L, 5L, 0.7, 0.4, 1, 0.5)
  d2 <- rbind(d, d[1:4, ]) # 4 replicated nodes
  expect_error(
    flexybayes(y ~ 1, random = ~ ar1(row):ar1(col), data = d2,
               backend = "inla"),
    class = "flexybayes_refusal_ar1_spatial_requires_complete_grid"
  )
})

test_that("auto refuses (never silently drops the residual) on an incomplete grid", {
  # Regression test, revised 2026-08-14. An incomplete/replicated grid is a
  # DATA-dependent runtime failure inside emit_inla(), not a structural
  # lgm_gate() refusal, so it reaches a separate auto-fallback branch in
  # .dispatch_backend().
  #
  # This test previously asserted that the fallback to brms SUCCEEDED, and
  # checked only the returned backend. That certified a model-fidelity bug:
  # brms has no residual-structure lowering, so the emitted Stan program was
  # an intercept-only iid Gaussian -- the requested AR1xAR1 residual vanished
  # with no error. Checking the wrapper class could never have caught it.
  #
  # The correct outcome is a refusal: INLA will not fit an incomplete grid,
  # brms cannot represent the residual, and greta is quarantined, so no
  # active backend can faithfully fit the model.
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(7L)
  d <- .spatial_sim(6L, 5L, 0.7, 0.4, 1, 0.5)
  d2 <- rbind(d, d[1:4, ]) # 4 replicated nodes
  expect_error(
    suppressMessages(flexybayes(
      y ~ 1,
      random = ~ ar1(row):ar1(col),
      data = d2,
      backend = "auto"
    )),
    class = "flexybayes_refusal_auto_no_active_route"
  )
})
