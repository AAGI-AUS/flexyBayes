# =============================================================================
# The grouped temporal AR1 field: random = ~ ar1(t):g.
#
# One AR1 series along `t` within each level of `g`, independent across
# levels, one correlation and one field SD shared by every level, plus the
# Gaussian observation nugget -- three hyperparameters. It is the
# units-by-occasions repeated-measures covariance of Verbyla and Cullis
# (1992) written on the random side, and it is emitted on INLA as
#   f(<t>_id, model = "ar1", group = <g>_id,
#     control.group = list(model = "iid")).
#
# The field has been reachable since the grouped-AR1 emit landed (ar1(t):g
# folds to the same node as ar1(row):id(col)), but it was described on
# every print surface as `ar1(t):ar1(g)` with four parameters, and the
# commuted spelling `g:ar1(t)` -- ASReml's usual order for a set of plots
# followed over seasons -- fell to an untyped interaction. This file pins
# the spelling, the emitted term, the description and the fit.
#
# FAITHFULNESS is checked against an INDEPENDENT oracle that does not call
# INLA: REML for y ~ 1 with V = ss2 (I_g (x) R_t(rho)) + se2 I, profiled in
# the eigenbasis of R_t (the same device test-inla-spatial-ar1.R uses for
# the separable grid, with the group correlation fixed at zero).
# =============================================================================

skip_if_no_inla <- function() skip_if_not_installed("INLA")

.tg_ar1_mat <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))

# n_g units, each a complete series of n_t occasions: AR1(rho) along t with
# field SD sd_s, independent across units, plus a nugget with SD sd_e.
# `t` varies fastest, so unit i occupies rows (i - 1) * n_t + 1 .. i * n_t.
.tg_sim <- function(n_t, n_g, rho, sd_s, sd_e, mu = 5, seed = 31L) {
  set.seed(seed)
  grid <- expand.grid(t = seq_len(n_t), g = seq_len(n_g))
  Lt <- t(chol(sd_s^2 * .tg_ar1_mat(n_t, rho)))
  u <- as.numeric(vapply(
    seq_len(n_g),
    function(i) as.numeric(Lt %*% stats::rnorm(n_t)),
    numeric(n_t)
  ))
  data.frame(
    y = mu + u + stats::rnorm(n_t * n_g, 0, sd_e),
    t = factor(grid$t),
    g = factor(grid$g)
  )
}

# INDEPENDENT oracle: REML point estimates of (rho, sd_s, sd_e) for
# y ~ 1, V = ss2 (I_g (x) R_t(rho)) + se2 I. The Kronecker with an identity
# on the unit side means V's eigenvalues are ss2 * lambda_k(R_t) + se2, each
# repeated n_g times, and the eigenvectors are I_g (x) U_t.
.tg_reml_oracle <- function(y, n_t, n_g) {
  X <- matrix(1, length(y), 1L)
  vy <- stats::var(y)
  nll <- function(par) {
    rr <- tanh(par[1]); ss2 <- exp(par[2]); se2 <- exp(par[3])
    et <- eigen(.tg_ar1_mat(n_t, rr), symmetric = TRUE)
    U <- kronecker(diag(n_g), et$vectors)
    dv <- as.numeric(ss2 * rep(et$values, times = n_g) + se2)
    if (any(dv <= 0)) return(1e10)
    yt <- crossprod(U, y); Xt <- crossprod(U, X)
    A <- crossprod(Xt, Xt / dv)
    beta <- solve(A, crossprod(Xt, yt / dv))
    rt <- yt - Xt %*% beta
    0.5 * (sum(log(dv)) + sum((rt^2) / dv) +
             determinant(A, logarithm = TRUE)$modulus[1])
  }
  best <- NULL
  for (s in list(c(atanh(0.3), log(vy / 2), log(vy / 2)),
                 c(atanh(0.7), log(vy / 2), log(vy / 4)))) {
    o <- tryCatch(stats::optim(s, nll, method = "Nelder-Mead",
                               control = list(maxit = 800L, reltol = 1e-10)),
                  error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  p <- best$par
  list(rho = tanh(p[1]), sd_s = sqrt(exp(p[2])), sd_e = sqrt(exp(p[3])))
}

.tg_inla_hyper <- function(fit) {
  hp <- fit$inla$summary.hyperpar
  gv <- function(p) {
    i <- grep(p, rownames(hp)); if (length(i)) hp$mean[i[1]] else NA_real_
  }
  list(
    rho = gv("^Rho for t_id"),
    sd_s = 1 / sqrt(gv("^Precision for t_id")),
    sd_e = 1 / sqrt(gv("Gaussian observations"))
  )
}

# ---- parse: three spellings, one node -----------------------------------------

test_that("ar1(t):g, g:ar1(t) and ar1(t):id(g) parse to the same grouped field", {
  d <- .tg_sim(6L, 4L, 0.5, 1, 0.5)
  node <- function(rhs) {
    fb <- fb_from_asreml(y ~ 1, random = rhs, data = d)
    expect_length(fb$random_terms, 1L)
    fb$random_terms[[1L]]
  }
  a <- node(~ ar1(t):g)
  b <- node(~ g:ar1(t))
  c3 <- node(~ ar1(t):id(g))
  for (term in list(a, b, c3)) {
    expect_identical(term$type, "ar1_spatial")
    expect_identical(term$row_var, "t")
    expect_identical(term$col_var, "g")
    expect_false(isTRUE(term$col_ar1))
    expect_null(term$at_var)
  }
  expect_identical(a[c("type", "row_var", "col_var", "col_ar1")],
                   b[c("type", "row_var", "col_var", "col_ar1")])
})

test_that("the commuted spelling does not disturb the separable grid or units", {
  d <- .tg_sim(6L, 4L, 0.5, 1, 0.5)
  d$row <- d$t; d$col <- d$g
  sep <- fb_from_asreml(y ~ 1, random = ~ ar1(row):ar1(col), data = d)
  expect_true(isTRUE(sep$random_terms[[1L]]$col_ar1))
  # A plain nested effect is still a nested effect.
  nest <- fb_from_asreml(y ~ 1, random = ~ g:t, data = d)
  expect_identical(nest$random_terms[[1L]]$type, "nested")
})

# ---- emit: the grouped-AR1 f() with an iid group model ----------------------

test_that("emit: the grouped field is AR1 along t grouped iid over g", {
  term <- list(type = "ar1_spatial", row_var = "t", col_var = "g",
               col_ar1 = FALSE)
  expect_identical(
    flexyBayes:::.inla_ar1_field_term(term),
    paste0("f(t_id, model = \"ar1\", group = g_id, ",
           "control.group = list(model = \"iid\"))")
  )
  expect_identical(flexyBayes:::.ar1_term_spelling(term), "ar1(t):g")
})

# ---- description: the model line says what the field is ---------------------

test_that("the model line and the ASReml note describe the grouped field as itself", {
  d <- .tg_sim(6L, 4L, 0.5, 1, 0.5)
  ir <- fb_from_asreml(y ~ 1, random = ~ g:ar1(t), data = d)
  line <- flexyBayes:::.fb_model_line(ir)
  expect_match(line, "ar1(t):g field, AR1 along t and independent across g", fixed = TRUE)
  expect_match(line, "3 parameters", fixed = TRUE)
  expect_false(grepl("ar1(g)", line, fixed = TRUE))
  expect_false(grepl("4 parameters", line, fixed = TRUE))
  expect_identical(
    flexyBayes:::.fb_asreml_residual_note(ir),
    "ASReml residual = ~ g:ar1(t) (2 parameters)"
  )
  # The separable grid keeps its four-parameter description.
  d$row <- d$t; d$col <- d$g
  sep <- fb_from_asreml(y ~ 1, random = ~ ar1(row):ar1(col), data = d)
  expect_match(flexyBayes:::.fb_model_line(sep),
               "ar1(row):ar1(col) field + nugget (4 parameters)", fixed = TRUE)
  expect_identical(flexyBayes:::.fb_asreml_residual_note(sep),
                   "ASReml residual = ~ ar1:ar1 (3 parameters)")
  one <- fb_from_asreml(y ~ 1, random = ~ ar1(t), data = d)
  expect_identical(flexyBayes:::.fb_asreml_residual_note(one),
                   "ASReml residual = ~ ar1(t) (2 parameters)")
})

# ---- gate and refusals -----------------------------------------------------------

test_that("the residual spelling in either order refuses and names the random-side field", {
  d <- .tg_sim(6L, 4L, 0.5, 1, 0.5)
  for (rhs in list(~ g:ar1(t), ~ ar1(t):g)) {
    res <- tryCatch(
      suppressMessages(flexybayes(y ~ 1, residual = rhs, data = d,
                                  backend = "inla", return_code = TRUE,
                                  verbose = FALSE)),
      error = function(e) e
    )
    expect_true(inherits(res, "flexybayes_refusal_ar1_residual_not_representable"))
    expect_match(conditionMessage(res), "random = ~ ar1(t):g", fixed = TRUE)
  }
})

test_that("brms refuses the grouped field by name", {
  skip_if_not_installed("brms")
  d <- .tg_sim(6L, 4L, 0.5, 1, 0.5)
  res <- tryCatch(
    suppressMessages(flexybayes(y ~ 1, random = ~ g:ar1(t), data = d,
                                backend = "brms", return_code = TRUE,
                                verbose = FALSE)),
    error = function(e) e
  )
  expect_true(inherits(res, "flexybayes_refusal_stan_cannot_represent_ar1_field"))
})

test_that("a non-Gaussian family carries the grouped field through the INLA gate", {
  skip_if_no_inla()
  d <- .tg_sim(6L, 4L, 0.5, 1, 0.5)
  d$y <- exp(d$y / 5)
  res <- tryCatch(
    suppressMessages(flexybayes(y ~ 1, random = ~ ar1(t):g, data = d,
                                backend = "inla", family = "gamma",
                                link = "log", return_code = TRUE,
                                verbose = FALSE)),
    error = function(e) e
  )
  expect_false(inherits(res, "error"))
})

# ---- fit: INLA recovers the grouped structure against the REML oracle ---------

test_that("INLA recovers the grouped temporal field against an independent REML oracle", {
  skip_if_no_inla()
  skip_on_cran()
  n_t <- 24L; n_g <- 12L
  d <- .tg_sim(n_t, n_g, rho = 0.6, sd_s = 1.0, sd_e = 0.5, seed = 41L)
  oracle <- .tg_reml_oracle(d$y, n_t, n_g)
  fit <- suppressMessages(flexybayes(y ~ 1, random = ~ g:ar1(t), data = d,
                                     backend = "inla", verbose = FALSE))
  got <- .tg_inla_hyper(fit)
  expect_lt(abs(got$rho - oracle$rho), 0.15)
  expect_lt(abs(got$sd_s - oracle$sd_s) / oracle$sd_s, 0.30)
  expect_lt(abs(got$sd_e - oracle$sd_e) / oracle$sd_e, 0.30)
  # The oracle itself is near the truth on a series this long.
  expect_lt(abs(oracle$rho - 0.6), 0.2)
  # No second correlation: the group model is iid, so INLA reports no
  # GroupRho for this field.
  expect_false(any(grepl("GroupRho", rownames(fit$inla$summary.hyperpar))))
})

test_that("an occasion a unit missed is carried as a design row with NA", {
  skip_if_no_inla()
  d <- .tg_sim(10L, 6L, 0.6, 1, 0.5, seed = 43L)
  d$y[c(3L, 27L, 50L)] <- NA
  fit <- suppressMessages(flexybayes(y ~ 1, random = ~ ar1(t):g, data = d,
                                     backend = "inla", verbose = FALSE))
  printed <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_match(printed, "60 design rows, 57 observed responses", fixed = TRUE)
  expect_match(printed, "Grouped AR1 field (AR1 along t, independent across g)",
               fixed = TRUE)
  expect_match(printed, "ar1(t):g field, AR1 along t and independent across g",
               fixed = TRUE)
  expect_match(printed, "ASReml residual = ~ g:ar1(t) (2 parameters)", fixed = TRUE)
  expect_false(grepl("ar1(t):ar1(g)", printed, fixed = TRUE))
  expect_false(grepl("Separable AR1 field", printed, fixed = TRUE))
})

test_that("a replicated node still refuses: the grid must be one row per (t, g)", {
  skip_if_no_inla()
  d <- .tg_sim(6L, 4L, 0.5, 1, 0.5)
  d2 <- rbind(d, d[1:2, ])
  # The lattice check runs on the fit path, not under return_code = TRUE,
  # so this is a live call: it refuses before INLA is reached.
  res <- tryCatch(
    suppressMessages(flexybayes(y ~ 1, random = ~ ar1(t):g, data = d2,
                                backend = "inla", verbose = FALSE)),
    error = function(e) e
  )
  expect_true(inherits(res, "flexybayes_refusal_ar1_spatial_requires_complete_grid"))
  expect_match(conditionMessage(res), "ar1(t):g", fixed = TRUE)
})
