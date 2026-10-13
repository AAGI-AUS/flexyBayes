# DEF-004: a variance component the data barely identify is flagged, and
# neither that flag nor the boundary flag depends on a covariate's units.
#
# The weak-identification flag reads how much of the prior's probability the
# posterior's 95% interval still spans. Measured on the default half-normal
# prior, 10 seeds each: 3 levels 0.55 to 0.80, a random copy of a fixed
# factor 0.88, 5 levels 0.19 to 0.36, 30 levels 0.07. The threshold is 0.4.

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

test_that("DEF-004: a random copy of a fixed factor is flagged as weakly identified", {
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

test_that("DEF-004: a grouping factor with three levels is flagged", {
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
    ratio <- attr(res$fit$extras$variance_comps, "interval_ratio")
    expect_gte(ratio[["sd_g"]], .FB_WEAK_IDENTIFICATION_RATIO)
  }
})

test_that("DEF-004: a well-identified grouping factor is not flagged", {
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

test_that("DEF-004: rescaling a slope covariate by 1000 changes neither flag", {
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

test_that("DEF-004: a slope SD is compared with the residual SD in response units", {
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

test_that("DEF-004: the weak-identification warning can be silenced by option", {
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

# The interval ratio, (q97.5 - q2.5) / median on the SD, and its threshold.

.wv_calibration <- function() {
  path <- system.file("extdata", "weak_identification_calibration.csv",
                      package = "flexyBayes")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "extdata",
                                "weak_identification_calibration.csv")
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

test_that("the interval-ratio threshold is the 99th percentile of its calibration rows", {
  x <- .wv_calibration()
  expect_identical(nrow(x), 80L)
  expect_equal(x$ratio, (x$q97.5 - x$q2.5) / x$median)
  calib <- x$levels >= 5L & x$design != "copy_of_fixed_6"
  expect_identical(sum(calib), 50L)
  q99 <- stats::quantile(x$ratio[calib], 0.99, names = FALSE)
  expect_equal(.FB_WEAK_IDENTIFICATION_RATIO, q99, tolerance = 0.005 / q99)
  # Every three-level row is above it and every four-level row below.
  expect_true(all(x$ratio[x$levels == 3L] >= .FB_WEAK_IDENTIFICATION_RATIO))
  expect_true(all(x$ratio[x$levels == 4L] < .FB_WEAK_IDENTIFICATION_RATIO))
  # The random copy of a fixed factor is below it, and the prior share
  # catches it.
  copy <- x$design == "copy_of_fixed_6"
  expect_true(all(x$ratio[copy] < .FB_WEAK_IDENTIFICATION_RATIO))
  expect_true(all(x$prior_share[copy] >= .FB_WEAK_IDENTIFICATION_SHARE))
})

test_that("either reading flags a component, and only an sd_ row", {
  cmp <- c("sigma", "sd_a", "sd_b", "sd_c", "sd_spatial")
  share <- c(sigma = 0.9, sd_a = 0.5, sd_b = 0.1, sd_c = 0.1, sd_spatial = 0.9)
  ratio <- c(sigma = 9, sd_a = 1, sd_b = 5, sd_c = 1, sd_spatial = 9)
  expect_identical(
    .fb_weakly_identified(cmp, share, ratio),
    c(FALSE, TRUE, TRUE, FALSE, FALSE)
  )
  expect_identical(
    .fb_weakly_identified(cmp, numeric(0), numeric(0)),
    rep(FALSE, 5L)
  )
})

# DEF-201. besag.met, random = ~ diag(env):gen, fitted by
# data-raw/besag_diag_gxe_varcomps.R. ASReml 4.2 puts the C2 and C5 genotype
# variances at 3.28 and 3.77 with standard errors 23.6 and 23.7.
.wv_besag <- function() {
  vc <- readRDS(testthat::test_path("fixtures", "besag_diag_gxe_varcomps.rds"))
  attr(vc, "interval_ratio") <- .fb_interval_ratio(vc)
  vc
}

test_that("DEF-201: besag C2 and C5 read on the SD scale stay below the ratio threshold", {
  vc <- .wv_besag()
  med <- attr(vc, "posterior_median")
  # The posterior medians are SDs. Squared, they are the variances 19.5 and
  # 19.6 that DEFECTS.tsv records.
  expect_equal(unname(med[c("sd_gen__envC2", "sd_gen__envC5")])^2,
               c(19.5, 19.6), tolerance = 0.02)
  ratio <- attr(vc, "interval_ratio")
  expect_equal(unname(ratio[c("sd_gen__envC2", "sd_gen__envC5")]),
               c(1.69, 1.75), tolerance = 0.02)
  # DEF-201 stays open: neither reading flags them.
  notes <- .fb_varcomp_notes(
    vc$component, vc$q97.5, med,
    prior_share = attr(vc, "prior_share"),
    interval_ratio = ratio
  )
  expect_identical(notes, rep("", nrow(vc)))
  fit <- list(extras = list(variance_comps = vc))
  expect_identical(.fb_weak_identification_reasons(fit), character(0))

  # The falsifier: C2's interval widened to 4.3 times its median is flagged.
  i <- match("sd_gen__envC2", vc$component)
  vc$q97.5[[i]] <- vc$q2.5[[i]] + 4.3 * med[["sd_gen__envC2"]]
  attr(vc, "interval_ratio") <- .fb_interval_ratio(vc)
  notes <- .fb_varcomp_notes(
    vc$component, vc$q97.5, med,
    prior_share = attr(vc, "prior_share"),
    interval_ratio = attr(vc, "interval_ratio")
  )
  expect_identical(vc$component[notes == "weakly identified"], "sd_gen__envC2")
  fit <- list(extras = list(variance_comps = vc))
  reasons <- .fb_weak_identification_reasons(fit)
  expect_length(reasons, 1L)
  expect_match(reasons, "sd_gen__envC2 .*4\\.30 times its median")
})
