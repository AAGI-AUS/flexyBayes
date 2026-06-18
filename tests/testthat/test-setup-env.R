# Tests for environment setup

test_that("setup_env creates id vectors for factors", {
  dat <- data.frame(
    y = rnorm(20),
    geno = factor(rep(paste0("G", 1:5), 4)),
    env = factor(rep(paste0("E", 1:4), each = 5))
  )
  fixed_info <- flexyBayes:::.parse_fixed(y ~ env, dat)
  ev <- new.env(parent = emptyenv())

  flexyBayes:::.setup_env(
    ev,
    fixed_info,
    list(),
    list(list(type = "units")),
    dat,
    list(),
    NULL
  )

  expect_true("env_id" %in% ls(ev))
  expect_true("n_env" %in% ls(ev))
  expect_equal(ev$n_env, 4)
  expect_equal(length(ev$env_id), 20)
  expect_true(all(ev$env_id %in% 1:4))
})

test_that("setup_env creates random effect ids", {
  dat <- data.frame(
    y = rnorm(20),
    geno = factor(rep(paste0("G", 1:5), 4))
  )
  fixed_info <- flexyBayes:::.parse_fixed(y ~ 1, dat)
  random_terms <- flexyBayes:::.parse_formula(~geno, dat)
  ev <- new.env(parent = emptyenv())

  flexyBayes:::.setup_env(
    ev,
    fixed_info,
    random_terms,
    list(list(type = "units")),
    dat,
    list(),
    NULL
  )

  expect_true("geno_id" %in% ls(ev))
  expect_true("n_geno" %in% ls(ev))
  expect_equal(ev$n_geno, 5)
})

test_that("setup_env creates nested ids", {
  dat <- data.frame(
    y = rnorm(12),
    block = factor(rep(1:3, 4)),
    rep = factor(rep(1:2, 6))
  )
  fixed_info <- flexyBayes:::.parse_fixed(y ~ 1, dat)
  random_terms <- flexyBayes:::.parse_formula(~ block:rep, dat)
  ev <- new.env(parent = emptyenv())

  flexyBayes:::.setup_env(
    ev,
    fixed_info,
    random_terms,
    list(list(type = "units")),
    dat,
    list(),
    NULL
  )

  expect_true("nested_id_rep_in_block" %in% ls(ev))
  expect_true("n_rep_in_block" %in% ls(ev))
})

test_that("setup_env stores known matrices", {
  dat <- data.frame(
    y = rnorm(10),
    geno = factor(paste0("G", 1:10))
  )
  G_mat <- diag(10)
  fixed_info <- flexyBayes:::.parse_fixed(y ~ 1, dat)
  random_terms <- flexyBayes:::.parse_formula(~ vm(geno, Gmat), dat)
  ev <- new.env(parent = emptyenv())

  flexyBayes:::.setup_env(
    ev,
    fixed_info,
    random_terms,
    list(list(type = "units")),
    dat,
    list(Gmat = G_mat),
    NULL
  )

  expect_true("Gmat" %in% ls(ev))
  expect_equal(ev$Gmat, G_mat)
})

test_that("setup_env errors on missing matrix", {
  dat <- data.frame(
    y = rnorm(10),
    geno = factor(paste0("G", 1:10))
  )
  fixed_info <- flexyBayes:::.parse_fixed(y ~ 1, dat)
  random_terms <- flexyBayes:::.parse_formula(~ vm(geno, Gmat), dat)
  ev <- new.env(parent = emptyenv())

  expect_error(
    flexyBayes:::.setup_env(
      ev,
      fixed_info,
      random_terms,
      list(list(type = "units")),
      dat,
      list(),
      NULL
    ),
    "not found in known_matrices"
  )
})

test_that("setup_env creates factor interaction ids", {
  dat <- data.frame(
    y = rnorm(12),
    A = factor(rep(1:2, 6)),
    B = factor(rep(1:3, 4))
  )
  fixed_info <- flexyBayes:::.parse_fixed(y ~ A:B, dat)
  ev <- new.env(parent = emptyenv())

  flexyBayes:::.setup_env(
    ev,
    fixed_info,
    list(),
    list(list(type = "units")),
    dat,
    list(),
    NULL
  )

  expect_true("A_x_B_id" %in% ls(ev))
  expect_true("n_A_x_B" %in% ls(ev))
})

test_that("setup_env stores weights", {
  dat <- data.frame(y = rnorm(10))
  fixed_info <- flexyBayes:::.parse_fixed(y ~ 1, dat)
  wts <- runif(10, 0.5, 2)
  ev <- new.env(parent = emptyenv())

  flexyBayes:::.setup_env(
    ev,
    fixed_info,
    list(),
    list(list(type = "units")),
    dat,
    list(),
    wts
  )

  expect_true("wt_atg" %in% ls(ev))
  expect_equal(ev$wt_atg, wts)
})

test_that("setup_env handles at_units rcov", {
  dat <- data.frame(
    y = rnorm(20),
    env = factor(rep(1:4, 5))
  )
  fixed_info <- flexyBayes:::.parse_fixed(y ~ 1, dat)
  rcov_terms <- flexyBayes:::.parse_formula(~ at(env):units, dat)
  ev <- new.env(parent = emptyenv())

  flexyBayes:::.setup_env(ev, fixed_info, list(), rcov_terms, dat, list(), NULL)

  expect_true("env_id" %in% ls(ev))
  expect_true("n_env" %in% ls(ev))
})
