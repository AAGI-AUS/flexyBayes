# =============================================================================
# plot(type = "variogram") and the engine-aware plot() default.
#
# Two things land here. The first is the display an ASReml user expects
# after fitting a spatial model: the empirical semivariance of the
# residuals against separation along the design index. The second is a live
# 0.9.0 defect the display sits next to -- plot() decided whether a fit had
# sampler draws by testing the since-withdrawn native engine's draws slot
# alone (see NEWS.md, 0.9.3), which is NULL on a brms fit, so
# plot(brms_fit) declined to draw its diagnostics while naming brms as a
# supported backend in the same sentence.
#
# The variogram uses OBSERVED ROWS ONLY. A row carried as a latent design
# cell has no observed response and therefore no residual; pairing an NA
# into a squared difference would drop whole lags rather than one pair, and
# the picture would quietly rest on less of the array than it appeared to.
# The subtitle prints both counts for the same reason.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


.pv_ar1_m <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))

.pv_grid <- function(n_row = 8L, n_col = 6L, rho_r = 0.6, rho_c = 0.3,
                     sd_s = 1.2, sd_e = 0.4, seed = 2026L) {
  set.seed(seed)
  g <- expand.grid(col = seq_len(n_col), row = seq_len(n_row))
  chol_l <- t(chol(
    sd_s^2 * kronecker(.pv_ar1_m(n_row, rho_r), .pv_ar1_m(n_col, rho_c))
  ))
  g$y <- 20 + as.numeric(chol_l %*% stats::rnorm(n_row * n_col)) +
    stats::rnorm(n_row * n_col, 0, sd_e)
  g$row <- factor(g$row)
  g$col <- factor(g$col)
  g
}

.pv_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    .local_envir = parent.frame()
  )
}

.pv_cache <- new.env(parent = emptyenv())

.pv_field_fit <- function() {
  if (is.null(.pv_cache$field)) {
    .pv_cache$field <- suppressMessages(flexybayes(
      y ~ 1,
      random = ~ ar1(row):ar1(col),
      data = .pv_grid(),
      backend = "inla",
      na_action = "augment",
      verbose = FALSE
    ))
  }
  .pv_cache$field
}

.pv_holed_fit <- function() {
  if (is.null(.pv_cache$holed)) {
    g <- .pv_grid()
    g$y[c(5L, 19L)] <- NA
    .pv_cache$holed <- suppressMessages(flexybayes(
      y ~ 1,
      random = ~ ar1(row):ar1(col),
      data = g,
      backend = "inla",
      na_action = "augment",
      verbose = FALSE
    ))
  }
  .pv_cache$holed
}

.pv_iid_fit <- function() {
  if (is.null(.pv_cache$iid)) {
    set.seed(88L)
    d <- data.frame(
      g = factor(rep(letters[1:6], each = 8L)),
      x = stats::rnorm(48L)
    )
    d$y <- 2 + 0.5 * d$x + stats::rnorm(48L)
    .pv_cache$iid <- suppressMessages(flexybayes(
      y ~ x,
      random = ~g,
      data = d,
      backend = "inla",
      aggregate = FALSE,
      verbose = FALSE
    ))
  }
  .pv_cache$iid
}

.pv_null_device <- function(.local_envir = parent.frame()) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  withr::defer(
    {
      grDevices::dev.off()
      unlink(path)
    },
    envir = .local_envir
  )
  invisible(path)
}


# ---------------------------------------------------------------- #
# 1. The display                                                    #
# ---------------------------------------------------------------- #

test_that("a complete field fit draws a variogram over the array", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pv_silence()
  .pv_null_device()
  fit <- .pv_field_fit()

  expect_no_error(tab <- plot(fit, type = "variogram"))
  expect_s3_class(tab, "data.frame")
  expect_identical(names(tab), c("lag_row", "lag_col", "semivariance",
    "n_pairs"))
  # Every lag of an 8 x 6 array except (0, 0), which is a plot paired with
  # itself and contributes no separation.
  expect_identical(nrow(tab), 8L * 6L - 1L)
  expect_identical(max(tab$lag_row), 7)
  expect_identical(max(tab$lag_col), 5)
  expect_identical(sum(tab$lag_row == 0 & tab$lag_col == 0), 0L)
  expect_true(all(tab$semivariance >= 0))
  expect_true(all(tab$n_pairs > 0L))
  # Pair counts add to every unordered pair of observed rows.
  expect_identical(sum(tab$n_pairs), as.integer(48L * 47L / 2L))
  # Along the row direction at zero column lag the semivariance rises with
  # separation, which is what an AR1 field looks like.
  strip <- tab[tab$lag_col == 0, ]
  strip <- strip[order(strip$lag_row), ]
  expect_lt(strip$semivariance[[1L]], strip$semivariance[[nrow(strip)]])
})

test_that("the variogram uses observed rows only (F-p2)", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pv_silence()
  .pv_null_device()
  fit <- .pv_holed_fit()

  expect_identical(stats::nobs(fit), 48L)
  expect_identical(as.integer(stats::nobs(fit, type = "observed")), 46L)

  tab <- plot(fit, type = "variogram")
  expect_false(anyNA(tab$semivariance))
  # 46 observed rows, so 46 * 45 / 2 pairs -- not 48 * 47 / 2.
  expect_identical(sum(tab$n_pairs), as.integer(46L * 45L / 2L))
  expect_identical(attr(tab, "n_observed"), 46L)
  expect_identical(attr(tab, "n_design"), 48L)
})

test_that("the titles say what the picture is and what it rests on", {
  # The title and subtitle are carried on the returned table as well as
  # drawn, from one construction, so asserting on them is asserting on
  # what the panel says.
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pv_silence()
  .pv_null_device()
  fit <- .pv_holed_fit()

  tab <- plot(fit, type = "variogram")
  expect_identical(attr(tab, "title"), "Empirical residual variogram")
  expect_identical(attr(tab, "subtitle"), "46 observed / 48 design rows")
  # No fitted-variogram overlay is claimed on the panel or in the sources.
  panel_text <- paste(
    attr(tab, "title"), attr(tab, "subtitle")
  )
  expect_false(grepl("fitted", panel_text, ignore.case = TRUE))
  expect_false(grepl("ASReml", panel_text, fixed = TRUE))
})


# ---------------------------------------------------------------- #
# 2. The refusal                                                    #
# ---------------------------------------------------------------- #

test_that("a fit with no design index refuses by name", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pv_silence()
  .pv_null_device()
  fit <- .pv_iid_fit()

  err <- tryCatch(
    plot(fit, type = "variogram"),
    condition = function(e) e
  )
  expect_s3_class(err, "flexybayes_variogram_requires_design_index")
  msg <- conditionMessage(err)
  expect_match(msg, "design index", fixed = TRUE)
  expect_match(msg, "ar1(row):ar1(col)", fixed = TRUE)
  expect_match(msg, "type = \"residuals\"", fixed = TRUE)
})

test_that("the refusal code is registered in the taxonomy", {
  entries <- ls(flexyBayes:::.refusal_registry)
  expect_true("variogram_requires_design_index" %in% entries)
})


# ---------------------------------------------------------------- #
# 3. The engine-aware default (F-p3, AD-6)                          #
# ---------------------------------------------------------------- #

test_that("the sampler-draws predicate reads the engine, not one slot", {
  # The 0.9.0 predicate tested a since-withdrawn native engine's draws
  # slot alone (see NEWS.md, 0.9.3). A brms fit carries its draws inside
  # the brmsfit, and an INLA fit carries none -- with that engine gone,
  # the predicate is now unconditionally `inherits(x$brms, "brmsfit")`.
  brms_like <- structure(
    list(brms = structure(list(), class = "brmsfit")),
    class = c("flexybayes_brms", "flexybayes", "list")
  )
  inla_like <- structure(
    list(inla = list()),
    class = c("flexybayes_inla", "flexybayes", "list")
  )

  expect_true(flexyBayes:::.fb_has_sampler_draws(brms_like))
  expect_false(flexyBayes:::.fb_has_sampler_draws(inla_like))

  expect_identical(flexyBayes:::.fb_default_plot_type(brms_like),
    "diagnostics")
  expect_identical(flexyBayes:::.fb_default_plot_type(inla_like), "residuals")
})

test_that("plot(inla_fit) draws residuals rather than declining traces", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pv_silence()
  .pv_null_device()
  fit <- .pv_iid_fit()

  # No message: the default is a display this fit can actually draw.
  expect_no_message(plot(fit))
  # And the explicit request still says why it cannot be honoured.
  expect_message(plot(fit, type = "diagnostics"), "no sampler draws")
})

test_that("plot(brms_fit) draws its diagnostics (AD-6, live 0.9.0 bug)", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  .pv_silence()
  .pv_null_device()

  set.seed(19L)
  d <- data.frame(
    g = factor(rep(letters[1:5], each = 8L)),
    x = stats::rnorm(40L)
  )
  d$y <- 1 + 0.4 * d$x + stats::rnorm(40L)
  fit <- suppressMessages(suppressWarnings(flexybayes(
    y ~ x, random = ~g, data = d, backend = "brms", aggregate = FALSE,
    n_samples = 200L, warmup = 200L, chains = 1L,
    verbose = FALSE, mcmc_verbose = FALSE
  )))

  # Before the fix this emitted a message naming brms as a supported
  # engine for trace/density plots while refusing to draw them on this
  # exact brms fit.
  expect_no_message(plot(fit))
  expect_no_message(plot(fit, type = "diagnostics"))
})

test_that("the explicit types still work and an unknown one refuses", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  .pv_silence()
  .pv_null_device()
  fit <- .pv_iid_fit()

  expect_no_error(plot(fit, type = "residuals"))
  expect_no_error(plot(fit, type = "effects"))
  expect_no_error(plot(fit, type = "variance"))
  expect_error(plot(fit, type = "traceplot"), "should be one of")
})
