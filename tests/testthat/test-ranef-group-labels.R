# ranef() names the group the user wrote, not the index INLA was given.
#
# emit_inla() builds an index column per structured random term,
# `rep_block_id` for `rep:block`, `gen_by_env_id` for an unstructured
# covariance, and INLA names each summary.random entry after the column it
# was indexed by. The variance components were mapped back to user labels
# when that machinery was built. The random-effect table was not, so
# ranef(fit) came back with a grouping factor called `rep_block_id`. Same
# defect class as the label leak closed for summary()$varcomp, reached
# through a different accessor.

test_that("a nested term keeps its own name in ranef()", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  data(john.alpha, package = "agridat")
  fit <- suppressMessages(flexybayes(
    yield ~ rep, random = ~ rep:block + gen,
    data = john.alpha, backend = "inla", verbose = FALSE
  ))
  r <- ranef(fit)
  expect_true("rep:block" %in% names(r))
  expect_true("gen" %in% names(r))
  expect_false(any(grepl("_id$", names(r))))
  # And inside the table, not only on the list name.
  expect_identical(unique(r[["rep:block"]]$group), "rep:block")
})

test_that("a plain grouping factor is left alone", {
  # The map must not rename a group whose index IS the variable the user
  # named, which is the common case and the one a blanket rewrite would
  # break.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  data(john.alpha, package = "agridat")
  fit <- suppressMessages(flexybayes(
    yield ~ rep, random = ~gen, data = john.alpha,
    backend = "inla", verbose = FALSE
  ))
  expect_identical(names(ranef(fit)), "gen")
})

test_that("the label map covers every structured term type it claims", {
  # Unit-level, so it runs without INLA: build the IR for each structured
  # term and check the generated index name maps to something a reader
  # would recognise.
  d <- expand.grid(rep = factor(1:2), gen = factor(1:6), env = factor(1:3))
  d$y <- stats::rnorm(nrow(d))
  d$x <- stats::rnorm(nrow(d))

  cases <- list(
    list(random = ~ rep:gen, index = "rep_gen_id", label = "rep:gen"),
    list(random = ~ diag(env):gen, index = "gen_at_env_1",
         label = "gen (env = 1)"),
    list(random = ~ us(env):gen, index = "gen_by_env_id", label = "gen:env")
  )
  for (cs in cases) {
    fb <- fb_from_asreml(y ~ 1, random = cs$random, data = d)
    map <- .fb_inla_group_labels(fb)
    expect_true(cs$index %in% names(map), info = cs$index)
    expect_identical(unname(map[[cs$index]]), cs$label)
  }
})

test_that("the correlated slope's two index blocks take distinct labels", {
  # They used to share one. The iid2d field is laid out intercepts then
  # slopes, so the two index names address two different blocks of
  # deviations, and mapping both to `x | g` returned a list with two
  # elements of the same name, `ranef(fit)[["x | g"]]` reached only the
  # first, and no caller could tell which block that was.
  d <- data.frame(y = stats::rnorm(40), x = stats::rnorm(40),
                  g = factor(rep(letters[1:8], each = 5)))
  fb <- fb_from_brms(y ~ x + (x | g), data = d)
  map <- .fb_inla_group_labels(fb)
  expect_identical(unname(map[["g_cor_x_id"]]), "(Intercept) | g")
  expect_identical(unname(map[["g_cor_x_slope_id"]]), "x | g")
})


# --- the levels inside the table, not only the name on it -------------- #
#
# Finding S01's defect D2. `.fb_random_inla()` restored a level label only
# where INLA's `ID` column came back CHARACTER, which is what a plain `~ g`
# produces. Every structured term is indexed by an integer column
# emit_inla() built, so `ID` is 1..n and the guard never fired: `rep:block`
# on john.alpha reported levels "1" to "18" against ASReml's
# `rep_R1:block_B1 ...`, and `(Days | Subject)` reported "1" to "36" with
# the subject identifiers sitting unused on the fit. Two of the four BLUP
# comparisons in the S01 programme could not be made at all, and aligning by
# row order was refused because `interaction()` orders `rep` fastest while
# ASReml orders it slowest, a positional match would have been silently
# wrong on exactly the row where it looked easiest.

test_that("an interaction group reports the interaction's own levels", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  e <- new.env()
  utils::data("john.alpha", package = "agridat", envir = e)
  d <- e$john.alpha

  fit <- suppressMessages(flexybayes(
    yield ~ rep, random = ~ rep:block + gen,
    data = d, backend = "inla", verbose = FALSE
  ))
  got <- ranef(fit)[["rep:block"]]$level

  # The expected labels are written out here instead of taken from
  # `interaction()`, which is the call the emitter itself makes: an oracle
  # that shares the code under test proves self-consistency. The index is
  # lexicographic with the FIRST variable slowest, so the order is
  # `order(rep, block)` over the observed cells.
  cells <- unique(d[c("rep", "block")])
  cells <- cells[order(
    as.integer(factor(cells$rep, levels = levels(factor(d$rep)))),
    as.integer(factor(cells$block, levels = levels(factor(d$block))))
  ), ]
  expected <- paste0(cells$rep, ":", cells$block)

  expect_identical(got, expected)
  expect_length(got, 18L)
  # And nothing that looks like a bare index survived.
  expect_false(any(grepl("^[0-9]+$", got)))
})


# --- every other generated index (review finding N5) -------------------- #
#
# The level map covered the interaction and the correlated slope and nothing
# else, so the terms below reported levels "1", "2", ... Each test takes the
# expected group and levels from the data (or from lme4 where lme4 fits the
# model), and then checks the levels sit on the right ROWS with an oracle the
# labelling code never touches.
#
# The oracle. For a Gaussian model the posterior mean of the linear predictor
# is the sum of the posterior means of its parts, and INLA reports that
# predictor per observation, from the observation's own row of data. So each
# observation's random effect, read out of ranef() BY ITS OWN LABELS and added
# to the fixed part, has to reproduce it. A level attached to the wrong row of
# INLA's table breaks that sum by the size of the effects themselves, which
# is why the fixtures carry strong effects against a small residual. The
# reproduction is not exact. INLA's predictor summary is computed on its
# own, and differed from the sum by about a thousandth when this was measured,
# so the gap is read relative to the spread of the effects.

.rgl_effect <- function(r, group, level) {
  tab <- r[[group]]
  if (is.null(tab)) {
    return(rep(NA_real_, length(level)))
  }
  tab$estimate[match(as.character(level), tab$level)]
}

.rgl_placement_gap <- function(fit, fixed, contrib) {
  eta <- fit$inla$summary.linear.predictor$mean
  max(abs(eta - fixed - contrib)) / stats::sd(contrib)
}

.rgl_intercept <- function(fit) {
  fit$inla$summary.fixed["(Intercept)", "mean"]
}

test_that("an uncorrelated slope reports lme4's levels, on the right rows", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- lme4::sleepstudy
  fit <- suppressMessages(flexybayes(
    Reaction ~ Days + (Days || Subject),
    data = d, backend = "inla", verbose = FALSE
  ))
  r <- ranef(fit)
  oracle <- rownames(lme4::ranef(
    lme4::lmer(Reaction ~ Days + (Days || Subject), data = d)
  )$Subject)

  expect_true("Days | Subject" %in% names(r))
  expect_identical(r[["Days | Subject"]]$level, oracle)
  expect_identical(r[["Subject"]]$level, oracle)

  fixed <- .rgl_intercept(fit) + fit$inla$summary.fixed["Days", "mean"] * d$Days
  contrib <- .rgl_effect(r, "Subject", d$Subject) +
    d$Days * .rgl_effect(r, "Days | Subject", d$Subject)
  expect_lt(.rgl_placement_gap(fit, fixed, contrib), 0.1)
})

test_that("diag(f):g reports each level of f as a group over the levels of g", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(301L)
  d <- expand.grid(
    gen = factor(sprintf("G%02d", 1:12)),
    env = factor(c("dry", "wet", "irrigated")),
    rep = 1:3
  )
  mat_u <- matrix(stats::rnorm(36L, sd = rep(c(0.8, 1.6, 2.4), each = 12L)),
                  12L, 3L)
  d$y <- 10 + mat_u[cbind(as.integer(d$gen), as.integer(d$env))] +
    stats::rnorm(nrow(d), 0, 0.3)
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ diag(env):gen, data = d, backend = "inla", verbose = FALSE
  ))
  r <- ranef(fit)

  groups <- paste0("gen (env = ", levels(d$env), ")")
  expect_setequal(names(r), groups)
  for (g in groups) {
    expect_identical(r[[g]]$level, levels(d$gen), info = g)
  }
  contrib <- vapply(seq_len(nrow(d)), function(i) {
    .rgl_effect(r, paste0("gen (env = ", d$env[i], ")"), d$gen[i])
  }, numeric(1L))
  expect_lt(.rgl_placement_gap(fit, .rgl_intercept(fit), contrib), 0.1)
})

test_that("us(f):g fitted by emit_inla() reports its cells as g:f", {
  # Refused at the gate, so the fit is made the one way it still can be.
  skip_if_not_installed("INLA")
  skip_on_cran()
  set.seed(302L)
  d <- expand.grid(
    gen = factor(sprintf("G%02d", 1:10)),
    env = factor(c("north", "south")),
    rep = 1:3
  )
  mat_u <- matrix(stats::rnorm(20L, sd = 1.5), 10L, 2L)
  d$y <- 10 + mat_u[cbind(as.integer(d$gen), as.integer(d$env))] +
    stats::rnorm(nrow(d), 0, 0.3)
  fb <- fb_from_asreml(y ~ 1, random = ~ us(env):gen, data = d)
  expect_s3_class(lgm_gate(fb), "lgm_refusal")
  fit <- suppressMessages(emit_inla(fb, data = d, verbose = FALSE))
  r <- ranef(fit)

  expect_named(r, "gen:env")
  cell <- paste0(d$gen, ":", d$env)
  expect_setequal(r[["gen:env"]]$level, unique(cell))
  expect_length(r[["gen:env"]]$level, length(unique(cell)))
  contrib <- .rgl_effect(r, "gen:env", cell)
  expect_lt(.rgl_placement_gap(fit, .rgl_intercept(fit), contrib), 0.1)
})

test_that("fa(f, k):g has no random-effect table to label, because it refuses", {
  # This block used to fit fa(env, 1):gen and check that ranef() labelled its
  # factor and specific groups over gen. The structure left the fitting route
  # at 0.13.0 (factor_analytic_not_repeatable): repeated identical calls did
  # not agree on the covariance where the genetic signal is weak, so there is
  # no fit to take a random-effect table from. The labelling code is
  # unchanged and still covered for every structure that does fit. What is
  # asserted here is that the refusal arrives instead of a fit, so this file
  # records the reason its fa() coverage went away instead of losing it
  # silently.
  set.seed(303L)
  d <- expand.grid(
    gen = factor(sprintf("G%02d", 1:15)),
    env = factor(c("E1", "E2", "E3")),
    rep = 1:2
  )
  vec_f <- stats::rnorm(15L, 0, 1.2)
  mat_u <- outer(vec_f, c(1, 1.5, -0.8)) +
    matrix(stats::rnorm(45L, 0, 0.3), 15L, 3L)
  d$y <- 10 + mat_u[cbind(as.integer(d$gen), as.integer(d$env))] +
    stats::rnorm(nrow(d), 0, 0.3)
  err <- tryCatch(
    suppressMessages(flexybayes(
      y ~ 1, random = ~ fa(env, 1):gen, data = d,
      backend = "inla", verbose = FALSE
    )),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_factor_analytic_not_repeatable")
})

# A relationship matrix among n levels, positive definite, named by level.
.rgl_kinship <- function(levs, seed) {
  set.seed(seed)
  n <- length(levs)
  mat_m <- matrix(stats::rnorm(n * 40L), n, 40L)
  out <- tcrossprod(scale(mat_m, scale = FALSE)) / 40 + diag(n) * 0.2
  dimnames(out) <- list(levs, levs)
  out
}

test_that("vm() reports the genotypes it was given, on the right rows", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  levs <- sprintf("line%02d", 1:10)
  mat_k <- .rgl_kinship(levs, 304L)
  d <- data.frame(geno = factor(rep(levs, 4L)))
  d$y <- 10 + stats::rnorm(10L, 0, 1.5)[as.integer(d$geno)] +
    stats::rnorm(nrow(d), 0, 0.3)
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ vm(geno, cov = fb_cov(Qk, type = "precision")),
    data = d, known_matrices = list(Qk = solve(mat_k)),
    backend = "inla", verbose = FALSE
  ))
  r <- ranef(fit)

  expect_named(r, "geno")
  expect_identical(r[["geno"]]$level, levels(d$geno))
  contrib <- .rgl_effect(r, "geno", d$geno)
  expect_lt(.rgl_placement_gap(fit, .rgl_intercept(fit), contrib), 0.1)
})

test_that("ped() reports the animals it was given, on the right rows", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  levs <- paste0("A", 101:110)
  mat_a <- .rgl_kinship(levs, 305L)
  d <- data.frame(animal = factor(rep(levs, 3L)))
  d$y <- 10 + stats::rnorm(10L, 0, 1.5)[as.integer(d$animal)] +
    stats::rnorm(nrow(d), 0, 0.3)
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ped(animal, Ainv, use_sparse_precision = TRUE),
    data = d, known_matrices = list(Ainv = solve(mat_a)),
    backend = "inla", verbose = FALSE
  ))
  r <- ranef(fit)

  expect_named(r, "animal")
  expect_identical(r[["animal"]]$level, levels(d$animal))
  contrib <- .rgl_effect(r, "animal", d$animal)
  expect_lt(.rgl_placement_gap(fit, .rgl_intercept(fit), contrib), 0.1)
})

test_that("ar1(t) reports the occasions it was fitted over", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(309L)
  n_t <- 30L
  d <- data.frame(t = factor(sprintf("wk%02d", seq_len(n_t))))
  d$y <- 10 + as.numeric(stats::arima.sim(list(ar = 0.8), n_t, sd = 1)) +
    stats::rnorm(n_t, 0, 0.2)
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ ar1(t), data = d, backend = "inla", verbose = FALSE
  ))
  r <- ranef(fit)

  expect_named(r, "ar1(t)")
  expect_identical(r[["ar1(t)"]]$level, levels(d$t))
  contrib <- .rgl_effect(r, "ar1(t)", d$t)
  expect_lt(.rgl_placement_gap(fit, .rgl_intercept(fit), contrib), 0.1)
})

test_that("one AR1 field per trial reports its plots as trial:row:col", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  set.seed(310L)
  n_row <- 5L
  n_col <- 4L
  ar1_cor <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))
  mat_l <- t(chol(1.5^2 * kronecker(ar1_cor(n_col, 0.5), ar1_cor(n_row, 0.7))))
  d <- expand.grid(
    row = factor(paste0("R", seq_len(n_row))),
    col = factor(paste0("C", seq_len(n_col))),
    trial = factor(c("Roseworthy", "Turretfield"))
  )
  field <- c(mat_l %*% stats::rnorm(n_row * n_col),
             mat_l %*% stats::rnorm(n_row * n_col))
  d$y <- 10 + field + stats::rnorm(nrow(d), 0, 0.3)
  # Rows in an order unrelated to the field's own layout, so the table's
  # order cannot be read off the data's.
  d <- d[sample(nrow(d)), ]
  fit <- suppressMessages(flexybayes(
    y ~ 1, random = ~ id(trial):ar1(row):ar1(col), data = d,
    backend = "inla", verbose = FALSE
  ))
  r <- ranef(fit)

  group <- "id(trial):ar1(row):ar1(col)"
  expect_named(r, group)
  plot_id <- paste(d$trial, d$row, d$col, sep = ":")
  expect_setequal(r[[group]]$level, plot_id)
  expect_length(r[[group]]$level, nrow(d))
  contrib <- .rgl_effect(r, group, plot_id)
  expect_lt(.rgl_placement_gap(fit, .rgl_intercept(fit), contrib), 0.1)
})

test_that("relabelling moves no number: nested and correlated-slope tables", {
  # ranef() may put names on INLA's table and pick one block out of a longer
  # field. It must not change a value, so each column is compared with the
  # engine's own table, row for row.
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  skip_if_not_installed("lme4")
  skip_on_cran()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  engine_cols <- c(estimate = "mean", std.error = "sd",
                   conf.low = "0.025quant", conf.high = "0.975quant")
  same_numbers <- function(got, tab, rows) {
    for (k in names(engine_cols)) {
      expect_identical(got[[k]], as.numeric(tab[[engine_cols[[k]]]][rows]),
                       info = k)
    }
  }

  e <- new.env()
  utils::data("john.alpha", package = "agridat", envir = e)
  fit <- suppressMessages(flexybayes(
    yield ~ rep, random = ~ rep:block + gen,
    data = e$john.alpha, backend = "inla", verbose = FALSE
  ))
  tab <- fit$inla$summary.random$rep_block_id
  same_numbers(ranef(fit)[["rep:block"]], tab, seq_len(nrow(tab)))

  # (x | g) is refused at the gate and fitted here by emit_inla() directly.
  set.seed(311L)
  n_g <- 12L
  d <- expand.grid(obs = 1:6, g = factor(sprintf("plot%02d", seq_len(n_g))))
  d$x <- stats::runif(nrow(d), -2, 2)
  mat_u <- matrix(stats::rnorm(2L * n_g), n_g, 2L) %*% diag(c(1.2, 0.6))
  d$y <- 5 + 1.5 * d$x + mat_u[as.integer(d$g), 1L] +
    mat_u[as.integer(d$g), 2L] * d$x + stats::rnorm(nrow(d), 0, 0.3)
  fb <- fb_from_brms(y ~ x + (x | g), data = d, family = "gaussian")
  fit <- suppressMessages(emit_inla(fb, data = d, verbose = FALSE))
  r <- ranef(fit)
  oracle <- rownames(lme4::ranef(lme4::lmer(y ~ x + (x | g), data = d))$g)

  expect_identical(r[["(Intercept) | g"]]$level, oracle)
  expect_identical(r[["x | g"]]$level, oracle)
  same_numbers(r[["(Intercept) | g"]], fit$inla$summary.random$g_cor_x_id,
               seq_len(n_g))
  same_numbers(r[["x | g"]], fit$inla$summary.random$g_cor_x_slope_id,
               n_g + seq_len(n_g))
})
