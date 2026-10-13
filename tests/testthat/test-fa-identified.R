# The identified covariance of a fa(f, k):g fit, G = Lambda Lambda' + Psi,
# and the agreement of several fits on it.
#
# G and the genetic correlations are the quantities compared across fits and
# against ASReml, because they do not change when the loadings are rotated
# or change sign. The loadings themselves do, so no check here reads one.

test_that(".fb_fa_g() is the closed form and ignores rotation and sign", {
  mat_l <- cbind(c(1.0, 1.4, 0.6, -0.8), c(0.0, 1.0, 0.7, 0.4))
  vec_sf <- c(0.9, 0.6)
  vec_psi <- c(0.25, 0.3, 0.2, 0.35)
  mat_ls <- mat_l %*% diag(vec_sf)
  expect_equal(.fb_fa_g(mat_l, vec_sf, vec_psi),
               mat_ls %*% t(mat_ls) + diag(vec_psi^2))
  # Rotate the scaled loadings by 0.7 radians and flip the sign of factor 2.
  rot <- matrix(c(cos(0.7), sin(0.7), -sin(0.7), cos(0.7)), 2L)
  mat_rot <- mat_ls %*% rot %*% diag(c(1, -1))
  expect_equal(.fb_fa_g(mat_rot, c(1, 1), vec_psi),
               .fb_fa_g(mat_l, vec_sf, vec_psi))
  # A changed loading changes G, so the comparison above is not vacuous.
  mat_bad <- mat_l
  mat_bad[3L, 1L] <- 0.9
  expect_gt(max(abs(.fb_fa_g(mat_bad, vec_sf, vec_psi) -
                      .fb_fa_g(mat_l, vec_sf, vec_psi))), 0.1)
})

test_that(".fb_fa_identified() reads the anchors and the zero triangle in fa order", {
  # Three levels, k = 2, anchors at levels B then A. Factor 1 is entered at
  # B, A, C, and factor 2 at A and C only.
  fit <- list(fb = list(random_terms = list(list(
    type = "fa_gxe", outer = "env", inner = "gen", k = 2L,
    outer_levels = c("A", "B", "C"), fa_order = c(2L, 1L, 3L)
  ))))
  draws <- data.frame(
    sd_fa1 = c(2, 2), sd_fa2 = c(1, 1),
    fa1_envA = c(0.5, 0.5), fa1_envC = c(-1, -1),
    fa2_envC = c(0.8, 0.8),
    sd_gen__envA = c(0.1, 0.1), sd_gen__envB = c(0.2, 0.2),
    sd_gen__envC = c(0.3, 0.3),
    check.names = FALSE
  )
  local_mocked_bindings(fb_variance_draws = function(object, n, seed) draws)
  got <- .fb_fa_identified(fit)
  mat_l <- rbind(A = c(0.5, 1), B = c(1, 0), C = c(-1, 0.8))
  want <- .fb_fa_g(mat_l, c(2, 1), c(0.1, 0.2, 0.3))
  dimnames(want) <- list(c("A", "B", "C"), c("A", "B", "C"))
  expect_equal(got$G, want)
  expect_equal(got$cor, stats::cov2cor(want))
  expect_identical(got$anchors, c("B", "A"))
})

test_that(".fb_fa_identified() refuses a fit with no fa() term", {
  fit <- list(fb = list(random_terms = list(list(type = "diag_gxe"))))
  expect_error(.fb_fa_identified(fit), "no fa\\(f, k\\):g term")
})

test_that(".fb_fa_agreement() is zero on identical fits and measures a spread", {
  mk <- function(g, r) {
    mat_g <- diag(g)
    mat_g[1L, 2L] <- mat_g[2L, 1L] <- r * sqrt(g[[1L]] * g[[2L]])
    dimnames(mat_g) <- list(c("E1", "E2"), c("E1", "E2"))
    list(G = mat_g, cor = stats::cov2cor(mat_g))
  }
  a <- mk(c(1, 2), 0.5)
  same <- .fb_fa_agreement(list(a, a, a))
  expect_identical(same$g_rel_max, 0)
  expect_identical(same$cor_abs_max, 0)
  expect_identical(same$n_fits, 3L)

  got <- .fb_fa_agreement(list(a, mk(c(1.1, 2), 0.62), mk(c(1.05, 2), 0.55)))
  expect_equal(got$g_rel_max, 0.1)
  expect_equal(unname(got$g_rel), c(0.1, 0))
  expect_equal(got$cor_abs_max, 0.12)
  expect_identical(names(got$cor_abs), "E1:E2")
})

test_that(".fb_fa_agreement() refuses one fit and fits over different levels", {
  mat_g <- diag(2)
  dimnames(mat_g) <- list(c("E1", "E2"), c("E1", "E2"))
  a <- list(G = mat_g, cor = mat_g)
  b <- a
  dimnames(b$G) <- dimnames(b$cor) <- list(c("E2", "E1"), c("E2", "E1"))
  expect_error(.fb_fa_agreement(list(a)), "at least two fits")
  expect_error(.fb_fa_agreement(list(a, b)), "same levels")
})

# ---------------------------------------------------------------- #
# The stored ASReml 4.2 reference fits                              #
# ---------------------------------------------------------------- #

.fa_asreml_record <- function(name) {
  path <- system.file("validation", "asreml-4.2", paste0(name, ".rds"),
                      package = "flexyBayes")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "validation", "asreml-4.2",
                                paste0(name, ".rds"))
  }
  readRDS(path)
}

.fa_asreml_names <- c("20_fa_besag_met", "21_fa_besag_met_30",
                      "22_fa_dasilva_maize")

test_that("each stored ASReml G is Lambda Lambda' + Psi of its own variance table", {
  # A hand edit to G, the correlations or the variance table would break the
  # agreement between them.
  for (nm in .fa_asreml_names) {
    rec <- .fa_asreml_record(nm)
    levs <- levels(rec$data$env)
    for (fk in names(rec$fits)) {
      x <- rec$fits[[fk]]
      vc <- stats::setNames(x$varcomp$component, x$varcomp$term)
      pre <- paste0("fa(env, ", x$k, "):gen!")
      mat_l <- vapply(seq_len(x$k), function(j) {
        unname(vc[paste0(pre, levs, "!fa", j)])
      }, numeric(length(levs)))
      mat_l <- matrix(mat_l, nrow = length(levs))
      vec_psi <- unname(vc[paste0(pre, levs, "!var")])
      want <- mat_l %*% t(mat_l) + diag(vec_psi)
      expect_equal(unname(x$G), want, label = paste(nm, fk))
      expect_equal(x$cor, stats::cov2cor(x$G), label = paste(nm, fk))
      expect_true(is.finite(x$loglik), label = paste(nm, fk))
    }
  }
})

test_that("the stored records say which ASReml fits converged", {
  conv <- vapply(.fa_asreml_names, function(nm) {
    rec <- .fa_asreml_record(nm)
    c(rec$fits$fa1$converged, rec$fits$fa2$converged)
  }, logical(2L))
  # besag.met at 30 genotypes does not converge with two factors.
  expect_identical(
    unname(conv),
    matrix(c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE), 2L)
  )
})

# ---------------------------------------------------------------- #
# The k-fit agreement gate on the INLA fa route                     #
# ---------------------------------------------------------------- #

# An fa(env, 1):gen term over three levels anchored at B, and a stand-in INLA
# fit carrying its joint mode: log precisions and copy scalings.
.fag_term <- function() {
  list(type = "fa_gxe", outer = "env", inner = "gen", k = 1L,
       outer_levels = c("A", "B", "C"), fa_order = c(2L, 1L, 3L))
}
.fag_fake <- function(mlik, beta_c = -1) {
  theta <- c(
    "Log precision for the Gaussian observations" = 0,
    "Log precision for gen_fa1_env_B" = -2 * log(2),
    "Log precision for gen_psi_env_A" = -2 * log(0.1),
    "Log precision for gen_psi_env_B" = -2 * log(0.2),
    "Log precision for gen_psi_env_C" = -2 * log(0.3),
    "Beta_intern for gen_fa1_env_A" = 0.5,
    "Beta_intern for gen_fa1_env_C" = beta_c
  )
  list(mlik = matrix(c(mlik, mlik), 2L, 1L),
       mode = list(mode.status = 0L, theta = theta))
}

# The k - 1 fresh-process fits replaced by `fits`, numbered by start. Fit j
# finishes `delays[[j]]` seconds after it starts.
.fag_mock <- function(fits, delays = 0, env = parent.frame()) {
  delays <- rep_len(delays, length(fits))
  n <- 0L
  testthat::local_mocked_bindings(
    .fb_inla_isolated_start = function(args) {
      n <<- n + 1L
      list(i = n, done_at = Sys.time() + delays[[n]])
    },
    .fb_inla_isolated_running = function(handle) Sys.time() < handle$done_at,
    .fb_inla_isolated_collect = function(handle) {
      list(ok = TRUE, value = fits[[handle$i]],
           seed = get0(".Random.seed", globalenv()))
    },
    .fb_inla_isolated_stop = function(handle) invisible(NULL),
    .fb_isolation_callr_available = function() TRUE,
    .fb_free_memory_bytes = function() 64 * 1024^3,
    .fb_inla_fit_bytes = function(fb, data) 1024^3,
    .env = env
  )
}

test_that(".fb_fa_mode_identified() is G at the joint mode, anchors read in fa order", {
  got <- .fb_fa_mode_identified(.fag_fake(-10), .fag_term())
  mat_l <- matrix(c(0.5, 1, -1), ncol = 1L)
  want <- .fb_fa_g(mat_l, 2, c(0.1, 0.2, 0.3))
  dimnames(want) <- list(c("A", "B", "C"), c("A", "B", "C"))
  expect_equal(got$G, want)
  expect_equal(got$cor, stats::cov2cor(want))
})

test_that("the gate records agreement when the k fits land on one mode", {
  .fag_mock(list(.fag_fake(-10), .fag_fake(-10)))
  rec <- .fb_inla_best_of(.fag_fake(-10), list(), 3L,
                          fa_term = .fag_term())$record
  gate <- rec$fa_gate
  expect_true(gate$agree)
  expect_identical(gate$g_rel_max, 0)
  expect_identical(gate$cor_abs_max, 0)
  expect_identical(gate$g_tol, 0.10)
  expect_identical(gate$cor_tol, 0.05)
  expect_identical(gate$n_fits, 3L)
})

test_that("the gate records disagreement when one fit lands elsewhere", {
  # Fit 3 puts the loading at C at -0.5 instead of -1: G_CC moves from 4.09
  # to 1.09.
  .fag_mock(list(.fag_fake(-10), .fag_fake(-11, beta_c = -0.5)))
  rec <- .fb_inla_best_of(.fag_fake(-10), list(), 3L,
                          fa_term = .fag_term())$record
  gate <- rec$fa_gate
  expect_false(gate$agree)
  expect_gt(gate$g_rel_max, 0.10)
  expect_identical(names(which.max(gate$g_rel)), "C")
})

test_that("the gate is recorded when the k fits ran in parallel", {
  # Launched before the session fit, run beside it, and finishing in reverse
  # order of start. Fit 3 lands on another mode.
  .fag_mock(list(.fag_fake(-10), .fag_fake(-11, beta_c = -0.5)),
            delays = c(0.3, 0.05))
  pool <- .fb_inla_best_of_launch(list(), 3L, NULL, NULL)
  rec <- .fb_inla_best_of(.fag_fake(-10), list(), 3L, pool,
                          fa_term = .fag_term())$record
  expect_true(rec$parallel)
  expect_identical(rec$fa_gate$n_fits, 3L)
  expect_false(rec$fa_gate$agree)
  expect_identical(names(which.max(rec$fa_gate$g_rel)), "C")
  # The same k fits run one at a time give the same gate.
  .fag_mock(list(.fag_fake(-10), .fag_fake(-11, beta_c = -0.5)))
  seq_rec <- .fb_inla_best_of(.fag_fake(-10), list(), 3L,
                              fa_term = .fag_term())$record
  expect_false(seq_rec$parallel)
  expect_identical(seq_rec$fa_gate, rec$fa_gate)
})

test_that("a fit without an fa() term records no gate", {
  .fag_mock(list(.fag_fake(-10), .fag_fake(-10)))
  rec <- .fb_inla_best_of(.fag_fake(-10), list(), 3L)$record
  expect_null(rec$fa_gate)
})

test_that("a live fa fit with k = 3 records the gate on its reproducibility record", {
  skip_if_not_long_tier()
  skip_if_not_installed("INLA")
  skip_if_not_installed("callr")
  withr::local_options(
    flexyBayes.inla_threads = "1:1",
    flexyBayes.inla_best_of = 3L,
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE
  )
  set.seed(51L)
  vec_f <- rnorm(30L, 0, 0.9)
  tab <- expand.grid(gen = factor(1:30), env = factor(1:4), rep = 1:2)
  tab$y <- 5 + outer(vec_f, c(1, 1.6, 0.5, -0.9))[cbind(tab$gen, tab$env)] +
    rnorm(nrow(tab), 0, 0.3)
  fit <- .fa_fit_past_refusal(y ~ env, ~ fa(env, 1):gen, tab)
  gate <- fit$extras$reproducibility$best_of$fa_gate
  expect_identical(gate$n_fits, 3L)
  expect_type(gate$agree, "logical")
  expect_setequal(names(gate$g_rel), as.character(1:4))
})
