# Tests for code generation (indexing approach, no model.matrix)
# Migrated from 04_test_translator.R

# Helper function
check_code <- function(
  code_str,
  must_contain = character(0),
  must_not_contain = character(0)
) {
  ok <- all(vapply(
    must_contain,
    function(p) grepl(p, code_str, fixed = TRUE),
    logical(1)
  ))
  ok <- ok &
    !any(vapply(
      must_not_contain,
      function(p) grepl(p, code_str, fixed = TRUE),
      logical(1)
    ))
  ok
}

# Create a minimal test dataset
make_test_data <- function(N = 120, n_geno = 10, n_env = 4, n_rep = 3) {
  dat <- expand.grid(
    geno = factor(paste0("G", seq_len(n_geno))),
    env = factor(paste0("E", seq_len(n_env))),
    rep = factor(paste0("R", seq_len(n_rep)))
  )
  dat$yield <- rnorm(nrow(dat), 50, 5)
  dat$x_cov <- rnorm(nrow(dat))
  dat$row <- factor(rep(seq_len(10), length.out = nrow(dat)))
  dat$col <- factor(rep(seq_len(12), length.out = nrow(dat)))
  dat$block <- factor(rep(seq_len(6), length.out = nrow(dat)))
  dat$bin_y <- rbinom(nrow(dat), 1, 0.5)
  dat$count_y <- rpois(nrow(dat), 10)
  dat
}

test_that("CG-1: intercept-only generates mu_atg, no model.matrix", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ 1,
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(is.character(code))
  expect_true(check_code(
    code,
    must_contain = c("mu_atg"),
    must_not_contain = c("model.matrix")
  ))
})

test_that("CG-1b: fixed factor uses indexing", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ 0 + env,
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(
    code,
    must_contain = c("alpha_env", "env_id"),
    must_not_contain = c("X_env %*%", "model.matrix")
  ))
})

test_that("CG-2: random geno uses u_geno[geno_id]", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ 1,
    random = ~geno,
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(
    code,
    must_contain = c("u_geno[geno_id]", "sigma_geno"),
    must_not_contain = c("model.matrix")
  ))
})

test_that("CG-3: vm(geno, Gmat) uses Cholesky", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  n_geno <- 10
  G_mat <- diag(n_geno) + 0.1
  code <- flexybayes(
    fixed = yield ~ 1,
    random = ~ vm(geno, Gmat),
    data = dat,
    known_matrices = list(Gmat = G_mat),
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(
    code,
    must_contain = c("t(chol(Gmat))", "sigma_geno"),
    must_not_contain = c("model.matrix")
  ))
})

test_that("CG-4a: at(env):geno generates DIAG structure", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ 1,
    random = ~ at(env):geno,
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(
    code,
    must_contain = c("sigma_geno_env", "cbind(geno_id, env_id)")
  ))
})

test_that("CG-4b: us(env):id(geno) generates LKJ correlation", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ 1,
    random = ~ us(env):id(geno),
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(
    code,
    must_contain = c("lkj_correlation", "cbind(geno_id, env_id)")
  ))
})

test_that("CG-4c: fa(env,2):id(geno) generates FA structure", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ 1,
    random = ~ fa(env, 2):id(geno),
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(code, must_contain = c("Lambda_", "psi_", "F_")))
})

test_that("CG-4c-id: fa loadings are identified (lower-tri, positive diagonal)", {
  # B8 regression guard. The factor-analytic loadings Lambda must be emitted with
  # the lower-triangular + positive-diagonal identification (Lopes & West 2004),
  # NOT as a free normal() matrix -- a free Lambda is rotation/sign/label-switch
  # unidentified for k > 1, so its posterior summaries are meaningless. The
  # strict-lower mask plus the positive-diagonal sweep ARE the identification.
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ 1,
    random = ~ fa(env, 2):id(geno),
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(
    code,
    must_contain = c(
      ".fa_lmask_", ".fa_dmask_", "Lambda_dvec_", "sweep(.fa_dmask_"
    )
  ))
  # the assembled Lambda is the masked construction, not a free normal()
  expect_false(grepl("Lambda_geno_env_fa2 <- normal\\(0, 1, dim", code))
})

test_that("CG-5: ar1(row):id(col) generates spatial terms", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ 1,
    random = ~ ar1(row):id(col),
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(
    code,
    must_contain = c("row_id", "col_id", "sigma_sp_")
  ))
})

test_that("CG-6: spl(x_cov) generates B-spline basis", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ 1,
    random = ~ spl(x_cov),
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(code, must_contain = c("splines::bs", "sigma_spl_")))
})

test_that("CG-8: at(env):units generates heterogeneous residual", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ 1,
    rcov = ~ at(env):units,
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(code, must_contain = c("sigma_e_atg", "env_id")))
})

test_that("CG-10: nested block:rep uses pre-computed ID", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ 1,
    random = ~ block:rep,
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(
    code,
    must_contain = c("nested_id_"),
    must_not_contain = c("factor(paste")
  ))
})

test_that("return_code=TRUE returns string, not fit", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  result <- flexybayes(
    fixed = yield ~ 1,
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(is.character(result))
  expect_true(nchar(result) > 50)
})

test_that("binomial family generates bernoulli distribution", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = bin_y ~ env,
    data = dat,
    family = "binomial",
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(code, must_contain = c("bernoulli", "ilogit")))
})

test_that("poisson family generates poisson distribution", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = count_y ~ env,
    data = dat,
    family = "poisson",
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(code, must_contain = c("poisson", "exp")))
})

test_that("probit link uses iprobit", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = bin_y ~ env,
    data = dat,
    family = "binomial",
    link = "probit",
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(code, must_contain = c("iprobit")))
})

test_that("unknown family raises error", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  expect_error(
    flexybayes(
      fixed = yield ~ 1,
      data = dat,
      family = "tweedie",
      return_code = TRUE,
      verbose = FALSE
    ),
    "Unsupported"
  )
})

test_that("missing known_matrix raises error", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  expect_error(
    flexybayes(
      fixed = yield ~ 1,
      random = ~ vm(geno, Gmat),
      data = dat,
      return_code = TRUE,
      verbose = FALSE
    ),
    "not found in known_matrices"
  )
})

test_that("missing response variable raises error", {
  skip_if_greta_backend_unusable()
  dat <- data.frame(x = rnorm(10))
  expect_error(
    flexybayes(fixed = y ~ x, data = dat, return_code = TRUE, verbose = FALSE),
    "not found"
  )
})

test_that("crossed random effects generate both terms", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ 1,
    random = ~ geno + env,
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(
    code,
    must_contain = c(
      "sigma_geno",
      "sigma_env",
      "u_geno[geno_id]",
      "u_env[env_id]"
    )
  ))
})

test_that("continuous covariate generates beta coefficient", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ x_cov,
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(
    code,
    must_contain = c("beta_x_cov", "as_data(x_cov)")
  ))
})

test_that("factor + covariate model has correct structure", {
  skip_if_greta_backend_unusable()
  dat <- make_test_data()
  code <- flexybayes(
    fixed = yield ~ env + x_cov,
    data = dat,
    return_code = TRUE,
    verbose = FALSE
  )
  expect_true(check_code(
    code,
    must_contain = c("mu_atg", "tau_env", "beta_x_cov")
  ))
})
