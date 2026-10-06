# DEF-004: a variance component the data barely identify is flagged, and
# neither that flag nor the boundary flag depends on a covariate's units.
#
# The weak-identification flag reads how much of the prior's probability the
# posterior's 95% interval still spans. Measured on the default half-normal
# prior, 10 seeds each: 3 levels 0.55 to 0.80, a random copy of a fixed
# factor 0.88, 5 levels 0.19 to 0.36, 30 levels 0.07; the threshold is 0.4.

.wv_silent_prior_note <- function(env = parent.frame()) {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE, .local_envir = env)
  threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")
  withr::defer(INLA::inla.setOption(num.threads = threads), envir = env)
}

# Fits and returns the fit with the warnings it raised.
.wv_fit <- function(expr) {
  msgs <- character(0)
  fit <- withCallingHandlers(
    expr,
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = msgs)
}

.wv_note <- function(fit, component) {
  vc <- summary(fit)$varcomp
  vc$note[vc$component == component]
}

.wv_intercepts <- function(levels, rows, seed) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(levels), each = rows)))
  d$y <- stats::rnorm(levels)[d$g] + stats::rnorm(nrow(d))
  d
}

.wv_slopes <- function(levels, seed, multiplier) {
  set.seed(seed)
  d <- data.frame(
    g = factor(rep(seq_len(levels), each = 10)),
    x = stats::rnorm(levels * 10)
  )
  d$y <- 1 + 0.5 * d$x + stats::rnorm(levels)[d$g] +
    stats::rnorm(levels, sd = 0.5)[d$g] * d$x + stats::rnorm(nrow(d))
  d$x <- d$x * multiplier
  d
}

test_that("a random copy of a fixed factor is flagged as weakly identified", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  .wv_silent_prior_note()
  set.seed(1)
  d <- data.frame(blk = factor(rep(1:6, each = 10)))
  # `plot` splits the rows exactly as `blk` does, under another name.
  d$plot <- factor(paste0("p", d$blk))
  d$y <- stats::rnorm(6)[d$blk] + stats::rnorm(nrow(d))
  res <- .wv_fit(flexybayes(
    y ~ blk, random = ~ plot, data = d, backend = "inla", verbose = FALSE
  ))
  expect_true(any(grepl("weakly identified", res$warnings, fixed = TRUE)))
  expect_true(any(grepl("sd_plot", res$warnings, fixed = TRUE)))
  expect_equal(.wv_note(res$fit, "sd_plot"), "weakly identified")
})

test_that("a grouping factor with three levels is flagged", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  .wv_silent_prior_note()
  for (seed in 1:3) {
    res <- .wv_fit(fb(
      y ~ 1 + (1 | g), data = .wv_intercepts(3, 10, seed),
      backend = "inla", verbose = FALSE
    ))
    expect_true(
      any(grepl("weakly identified", res$warnings, fixed = TRUE)),
      label = paste("seed", seed)
    )
    expect_equal(.wv_note(res$fit, "sd_g"), "weakly identified",
                 label = paste("seed", seed))
  }
})

test_that("a well-identified grouping factor is not flagged", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  .wv_silent_prior_note()
  for (seed in 1:3) {
    res <- .wv_fit(fb(
      y ~ 1 + (1 | g), data = .wv_intercepts(30, 8, seed),
      backend = "inla", verbose = FALSE
    ))
    expect_false(
      any(grepl("weakly identified", res$warnings, fixed = TRUE)),
      label = paste("seed", seed)
    )
    expect_equal(.wv_note(res$fit, "sd_g"), "", label = paste("seed", seed))
  }
})

test_that("rescaling a slope covariate by 1000 changes neither flag", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  .wv_silent_prior_note()
  for (levels in c(3, 20)) {
    notes <- character(0)
    weak <- logical(0)
    boundary <- logical(0)
    for (multiplier in c(1, 1000)) {
      res <- .wv_fit(fb(
        y ~ x + (1 + x || g), data = .wv_slopes(levels, 1, multiplier),
        backend = "inla", verbose = FALSE
      ))
      notes <- c(notes, .wv_note(res$fit, "sd_x_g"))
      weak <- c(weak, any(grepl("sd_x_g", res$warnings, fixed = TRUE) &
        grepl("weakly identified", res$warnings, fixed = TRUE)))
      boundary <- c(boundary, any(grepl("at the boundary", res$warnings,
                                        fixed = TRUE)))
    }
    label <- paste(levels, "levels")
    expect_identical(notes[[1L]], notes[[2L]], label = label)
    expect_identical(weak[[1L]], weak[[2L]], label = label)
    expect_identical(boundary, c(FALSE, FALSE), label = label)
    # Three levels must fire and twenty must not, at both scales.
    expect_identical(weak, rep(levels == 3, 2L), label = label)
  }
})

test_that("a slope SD is compared with the residual SD in response units", {
  vc <- data.frame(
    component = c("sigma", "sd_x_g"),
    estimate = c(1, 0.0005),
    sd = c(0.1, 0.0001),
    q2.5 = c(0.9, 0.0002),
    q97.5 = c(1.1, 0.0008),
    stringsAsFactors = FALSE
  )
  attr(vc, "posterior_median") <- c(sigma = 1, sd_x_g = 0.0005)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rep(c(-1, 1), 50) * 1000)
  d$y <- 0
  fit <- structure(
    list(
      fb = list(
        response = "y",
        random_terms = list(list(
          type = "simple_slope_uncor", var = "g", slope_var = "x"
        ))
      ),
      data = d,
      extras = list(variance_comps = vc)
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
  # 0.0008 per unit of x is 0.8 per standard deviation of x.
  expect_length(.fb_boundary_collapse_reasons(fit), 0L)
  expect_equal(flexyBayes:::.fb_summary_varcomp(fit)$note, c("", ""))

  # The falsifier: a slope that is negligible in response units is flagged.
  fit$extras$variance_comps$q97.5[[2L]] <- 1e-8
  expect_length(.fb_boundary_collapse_reasons(fit), 1L)
  expect_equal(flexyBayes:::.fb_summary_varcomp(fit)$note, c("", "collapsed"))
})

test_that("the weak-identification warning can be silenced by option", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  .wv_silent_prior_note()
  withr::local_options(flexyBayes.silence_weak_identification_warning = TRUE)
  res <- .wv_fit(fb(
    y ~ 1 + (1 | g), data = .wv_intercepts(3, 10, 1),
    backend = "inla", verbose = FALSE
  ))
  expect_false(any(grepl("weakly identified", res$warnings, fixed = TRUE)))
  # The table still says so.
  expect_equal(.wv_note(res$fit, "sd_g"), "weakly identified")
})
