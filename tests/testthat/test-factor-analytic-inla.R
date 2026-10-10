# Factor-analytic genotype-by-environment covariance: fa(f, k):g.
#
# ASReml's headline multi-environment-trial model. INLA has no latent model
# called "factor analytic" and does not need one: each common factor is a
# single latent field entered once per level of the outer factor as a SCALED
# COPY, and `copy=`'s scaling parameter IS the loading. Factor j is anchored
# at level j and entered only from level j onward, the zero upper triangle
# ASReml uses to identify the loading matrix against rotation.
#
# That emit works and is still here. What is withdrawn at 0.13.0 is the claim
# attached to its output. The structure is refused on both engines
# (factor_analytic_not_repeatable) because repeated identical calls did not
# agree on the covariance, and brms has no factor-analytic lowering to check
# it against.
#
# The measurement, made against a tolerance fixed before the runs and not
# moved afterwards: the reference is the identified covariance
# G = Lambda Lambda' + diag(psi), which is invariant to the rotation and sign
# the raw loadings are free in. Eight identical calls per design, tolerance a
# ratio of 1.10 on each diagonal entry of G and a range of 0.05 on each
# genetic correlation. The two simulated designs this file used to assert
# recovery on repeated to 1.003 and 1.007, and agridat::dasilva.maize to
# 1.0001; agridat::besag.met at 30 genotypes, where the genetic signal is
# weak, spanned 1.23 and 0.12 on INLA's default threads and 1.20 and 0.09 on
# one thread. Every call reported mode.status = 0 and a passing numerical
# confirm. Evidence: review/units/U-FA-1_2026-09-16/.
#
# The live-fit recovery blocks that stood here (loadings at k = 1 including a
# negative one, the k = 2 anchoring pattern, the specific variance under the
# package prior, and the summary-label guard) were removed with the fitting
# route: each fitted the structure, and none of them asked whether a second
# identical call returned the same answer. The emit-level blocks below are
# kept, because the emit is unchanged and still reachable through
# emit_inla(), as the us() branch is.

.fa_fixture <- function(seed = 51L, n_gen = 80L, n_rep = 6L,
                        loadings = c(1.0, 1.6, 0.5, -0.9),
                        sd_f = 0.9, sd_spec = 0.25, sd_resid = 0.30) {
  set.seed(seed)
  n_env <- length(loadings)
  vec_f <- rnorm(n_gen, 0, sd_f)
  mat_u <- outer(vec_f, loadings) +
    matrix(rnorm(n_gen * n_env, 0, sd_spec), ncol = n_env)
  tab <- expand.grid(
    gen = factor(seq_len(n_gen)),
    env = factor(seq_len(n_env)),
    rep = seq_len(n_rep)
  )
  tab$y <- 5 + mat_u[cbind(as.integer(tab$gen), as.integer(tab$env))] +
    rnorm(nrow(tab), 0, sd_resid)
  tab
}

# ---------------------------------------------------------------- #
# The emit is unchanged                                             #
# ---------------------------------------------------------------- #

test_that("the INLA gate still admits fa(f, k):g", {
  # The gate answers "can INLA represent this", and it can. The refusal is
  # not here: it is dispatch-wide (.refuse_unrepresentable_structures()),
  # because the limit is not INLA's representation but that the fitted
  # numbers do not repeat, and brms has nothing to route to. Keeping the
  # gate accurate about representability is what leaves the emit reachable
  # through a direct emit_inla() call.
  tab <- .fa_fixture(n_gen = 12L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ fa(env, 1):gen,
    data = tab,
    family = "gaussian"
  )
  expect_true("lgm_compatible" %in% lgm_gate(fb)$capabilities)
})

test_that("the emitted formula carries scaled copies and specific variances", {
  tab <- .fa_fixture(n_gen = 12L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ fa(env, 1):gen,
    data = tab,
    family = "gaussian"
  )
  txt <- gsub("\\s+", " ", paste(deparse(.build_inla_formula(fb)),
                                collapse = " "))
  anchor <- paste0("gen_fa1_env_", fb$random_terms[[1L]]$fa_order[[1L]])
  # The anchor is an ordinary iid field.
  expect_match(txt, paste0("f(", anchor, ", model = \"iid\""), fixed = TRUE)
  # Every other level is a scaled copy of it -- that scaling is the loading.
  expect_match(txt, paste0("copy = \"", anchor, "\", fixed = FALSE"),
               fixed = TRUE)
  # And one specific variance per level.
  expect_match(txt, "f(gen_psi_env_4, model = \"iid\"", fixed = TRUE)
})

# ---------------------------------------------------------------- #
# The anchor levels are chosen from the data                        #
# ---------------------------------------------------------------- #

# A k = 2 fixture whose first level loads on neither factor, so factor
# order would anchor both factors where the true loadings are near zero.
.fa_fixture_k2_weak_first <- function(seed = 77L, n_gen = 100L, n_rep = 4L) {
  set.seed(seed)
  mat_l <- cbind(
    c(0.05, 1.6, 0.6, -0.8, 1.1),
    c(0.02, 0.0, 1.2, 0.4, -0.9)
  )
  n_env <- nrow(mat_l)
  mat_u <- cbind(rnorm(n_gen, 0, 0.9), rnorm(n_gen, 0, 0.8)) %*% t(mat_l) +
    matrix(rnorm(n_gen * n_env, 0, 0.25), ncol = n_env)
  tab <- expand.grid(
    gen = factor(seq_len(n_gen)),
    env = factor(seq_len(n_env)),
    rep = seq_len(n_rep)
  )
  tab$y <- 5 + mat_u[cbind(as.integer(tab$gen), as.integer(tab$env))] +
    rnorm(nrow(tab), 0, 0.3)
  tab
}

test_that("a level with a near-zero loading is not chosen as the anchor", {
  # Level 1 has true loading 0.05. Anchored there, every other loading is a
  # ratio to almost zero.
  tab <- .fa_fixture(n_gen = 60L, n_rep = 3L,
                     loadings = c(0.05, 1.6, 0.5, -0.9))
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ fa(env, 1):gen,
    data = tab,
    family = "gaussian"
  )
  term <- fb$random_terms[[1L]]
  # Level 2 carries the largest loading.
  expect_identical(term$fa_order[[1L]], 2L)
  expect_setequal(term$fa_order, seq_len(4L))
  # The level order the rest of the package reads is unchanged.
  expect_identical(term$outer_levels, as.character(1:4))
  txt <- gsub("\\s+", " ", paste(deparse(.build_inla_formula(fb)),
                                collapse = " "))
  expect_match(txt, "f(gen_fa1_env_2, model = \"iid\"", fixed = TRUE)
  expect_match(txt, "f(gen_fa1_env_1, copy = \"gen_fa1_env_2\"", fixed = TRUE)
  expect_false(grepl("f(gen_fa1_env_1, model", txt, fixed = TRUE))
})

test_that("with k = 2 the second anchor is a different level carrying factor 2", {
  tab <- .fa_fixture_k2_weak_first()
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ fa(env, 2):gen,
    data = tab,
    family = "gaussian"
  )
  ord <- fb$random_terms[[1L]]$fa_order
  # Levels 2 and 3 carry the largest loadings on the two factors; level 1,
  # which factor order would anchor first, carries almost none.
  expect_identical(ord[1:2], c(2L, 3L))
  expect_setequal(ord, seq_len(5L))
  txt <- gsub("\\s+", " ", paste(deparse(.build_inla_formula(fb)),
                                collapse = " "))
  expect_match(txt, "f(gen_fa1_env_2, model = \"iid\"", fixed = TRUE)
  expect_match(txt, "f(gen_fa2_env_3, model = \"iid\"", fixed = TRUE)
  # Factor 2 is not entered at the first anchor: the zero upper triangle.
  expect_false(grepl("gen_fa2_env_2", txt, fixed = TRUE))
  expect_match(txt, "f(gen_fa2_env_1, copy = \"gen_fa2_env_3\"", fixed = TRUE)
})

test_that("the index columns follow the chosen anchors at k = 2", {
  # With k = 2 the data must carry a factor-2 column at every level after
  # the second anchor, in anchor order. Columns built in factor order would
  # leave the formula naming a column that does not exist.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.inla_threads = "1:1",
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE
  )
  fit <- .fa_fit_past_refusal(y ~ env, ~ fa(env, 2):gen,
                              .fa_fixture_k2_weak_first(n_gen = 60L))
  expect_s3_class(fit, "flexybayes_inla")
  expect_identical(.fb_fa_identified(fit)$anchors, c("2", "3"))
})

test_that("the chosen anchor recovers G on a fixture with a near-zero first loading", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.inla_threads = "1:1",
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE
  )
  vec_l <- c(0.05, 1.6, 0.5, -0.9)
  tab <- .fa_fixture(n_gen = 80L, n_rep = 4L, loadings = vec_l)
  mat_true <- 0.9^2 * outer(vec_l, vec_l) + diag(0.25^2, 4L)
  fit <- .fa_fit_past_refusal(y ~ env, ~ fa(env, 1):gen, tab)
  got <- .fb_fa_identified(fit)
  expect_identical(got$anchors, "2")
  rel <- abs(diag(got$G) - diag(mat_true)) / diag(mat_true)
  expect_lt(max(rel), 0.30)
  r_true <- stats::cov2cor(mat_true)
  off <- upper.tri(r_true)
  expect_lt(max(abs(got$cor[off] - r_true[off])), 0.10)
})

# ---------------------------------------------------------------- #
# The fitting route refuses, on every backend                       #
# ---------------------------------------------------------------- #

test_that("fa(f, k):g is refused on inla, brms and auto alike", {
  # One reason code for all three. Before 0.13.0 backend = "inla" fitted,
  # backend = "brms" gave that engine's own brms_cannot_represent_term, and
  # backend = "auto" fitted on INLA -- three different answers to one
  # question about whether the package will report this covariance.
  tab <- .fa_fixture(n_gen = 8L, n_rep = 2L)
  for (be in c("inla", "brms", "auto")) {
    err <- tryCatch(
      suppressMessages(flexybayes(
        y ~ env,
        random = ~ fa(env, 1):gen,
        data = tab,
        backend = be,
        verbose = FALSE
      )),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_refusal_factor_analytic_not_repeatable")
    expect_identical(err$reason_code, "factor_analytic_not_repeatable")
  }
})

test_that("the refusal fires before any engine runs", {
  # Code generation never reaches an engine, so a refusal raised on the
  # return_code = TRUE path is raised before the fit. This also means the
  # refusal does not depend on INLA or brms being installed.
  tab <- .fa_fixture(n_gen = 8L, n_rep = 2L)
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ env,
      random = ~ fa(env, 2):gen,
      data = tab,
      backend = "inla",
      return_code = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_factor_analytic_not_repeatable")
})

test_that("the refusal names the rank, the term and what to write instead", {
  tab <- .fa_fixture(n_gen = 8L, n_rep = 2L)
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ env,
      random = ~ fa(env, 2):gen,
      data = tab,
      backend = "inla",
      verbose = FALSE
    )),
    error = function(e) e
  )
  msg <- conditionMessage(err)
  expect_match(msg, "fa(env, 2):gen", fixed = TRUE)
  expect_match(msg, "diag(env):gen", fixed = TRUE)
  expect_match(msg, "us(env):gen", fixed = TRUE)
  expect_match(msg, "marks some variances the data barely identify as `weakly identified`",
               fixed = TRUE)
  # It says why, with the coverage the simulation test required.
  expect_match(msg, "134 of 150", fixed = TRUE)
  expect_match(msg, "against 136 required", fixed = TRUE)
  # And it does not offer the other engine as a way to get the fit, because
  # brms has no lowering for the structure at all.
  expect_false(grepl("backend = \"brms\"", msg, fixed = TRUE))
})

test_that("the refusal is targeted: the same data still fits without fa()", {
  # A blanket refusal that caught neighbouring models would be a regression
  # dressed as a fix. diag(env):gen is the structure the message recommends,
  # on the same rows.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  tab <- .fa_fixture(n_gen = 12L, n_rep = 3L)
  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~ diag(env):gen,
    data = tab,
    backend = "inla",
    verbose = FALSE
  ))
  expect_s3_class(fit, "flexybayes")
  expect_true(any(grepl("^sd_gen__env", summary(fit)$varcomp$component)))
})

test_that("the bare fa(v, k) spelling keeps its own refusal", {
  # fa(v, k) with no grouping factor has no covariance to apply and no emit
  # anywhere; it is a different refusal and must not be swallowed by the new
  # one.
  tab <- .fa_fixture(n_gen = 8L, n_rep = 2L)
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ env,
      random = ~ fa(env, 1),
      data = tab,
      backend = "inla",
      verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_fa_not_representable")
})

test_that("fa(f, k):g still parses, so the catalogue keeps the term", {
  # The refusal is about what the package will report, not about what it can
  # read. fb_terms() stays a formula catalogue.
  tab <- .fa_fixture(n_gen = 8L, n_rep = 2L)
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ fa(env, 2):gen,
    data = tab,
    family = "gaussian"
  )
  expect_identical(fb$random_terms[[1L]]$type, "fa_gxe")
  expect_identical(fb$random_terms[[1L]]$k, 2L)
  expect_identical(fb$random_terms[[1L]]$outer, "env")
  expect_identical(fb$random_terms[[1L]]$inner, "gen")
})
