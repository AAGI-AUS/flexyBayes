# =============================================================================
# The default variance scale reads the response the model is written on.
#
# The default is half_normal(scale = s) with s = 5 x the SD of the modelled
# response, over the adopted divisor. The scale (the bound of the uniform
# default before 0.13.0) was read from data[[response]], which is NULL for a
# transformed response such as log(y) or I(10 * y), so it fell back to 1 on
# both engines whatever the units. It is now read through the resolver the
# rest of the fit uses, after the offset policy.
#
# A NaN or infinite response is not a missing observation. INLA fitted such a
# model without a word, and brms failed with a message that named neither the
# row nor the value. The fit refuses by name before any engine is called.
# =============================================================================

.resp_data <- function(seed = 9L, n = 200L) {
  set.seed(seed)
  d <- data.frame(
    x = stats::rnorm(n),
    g = factor(rep(seq_len(10L), length.out = n)),
    o = stats::runif(n, 0, 2)
  )
  y <- 3 + 2 * d$x + stats::rnorm(n, 0, 1.5)
  # positive, so log(y) is defined
  d$y <- y - min(y) + 1
  d
}

# One entry per response form: the formula and the quantity the residual is
# written on, computed here from the data, not by the package.
.resp_cases <- function(d) {
  d_na <- d
  d_na$y[c(3L, 50L)] <- NA
  list(
    list(label = "y", fixed = y ~ x, data = d, modelled = d$y),
    list(label = "I(10 * y)", fixed = I(10 * y) ~ x, data = d, modelled = 10 * d$y),
    list(label = "log(y)", fixed = log(y) ~ x, data = d, modelled = log(d$y)),
    list(
      label = "log(y) with an offset",
      fixed = log(y) ~ x + offset(o),
      data = d,
      modelled = log(d$y) - d$o
    ),
    list(
      label = "log(y) with missing responses",
      fixed = log(y) ~ x,
      data = d_na,
      modelled = log(d_na$y)
    )
  )
}

.five_sd <- function(v) {
  5 * stats::sd(v[is.finite(v)]) / flexyBayes:::.FB_DEFAULT_SCALE_DIVISOR
}

.brms_sigma_bound <- function(code) {
  code <- paste(as.character(code), collapse = "\n")
  hit <- regmatches(
    code,
    regexpr("normal_lpdf\\(sigma \\| 0, [0-9.eE+-]+\\)", code)
  )
  if (length(hit) != 1L) {
    return(NA_real_)
  }
  as.numeric(sub("^normal_lpdf\\(sigma \\| 0, ([0-9.eE+-]+)\\)$", "\\1", hit))
}

# ---------------------------------------------------------------- #
# 1. The default scale on each engine                                #
# ---------------------------------------------------------------- #

test_that("the default residual scale on brms is five SDs of the modelled response", {
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  cases <- .resp_cases(.resp_data())

  for (cs in cases) {
    # The missing response is deliberate. brms reports the expected omission
    # during both prior reachability and Stan-code construction. Muffle only
    # that exact dependency warning and leave every other warning visible.
    code <- withCallingHandlers(
      suppressMessages(flexybayes(
        fixed = cs$fixed,
        random = ~g,
        data = cs$data,
        backend = "brms",
        return_code = TRUE
      )),
      warning = function(w) {
        if (identical(conditionMessage(w), "Rows containing NAs were excluded from the model.")) {
          invokeRestart("muffleWarning")
        }
      }
    )
    expect_equal(
      .brms_sigma_bound(code),
      .five_sd(cs$modelled),
      tolerance = 1e-6,
      label = paste0("brms sigma scale, ", cs$label)
    )
  }
})

test_that("the default residual and group scales on INLA are five SDs of the modelled response", {
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  cases <- .resp_cases(.resp_data())

  for (cs in cases) {
    code <- suppressMessages(flexybayes(
      fixed = cs$fixed,
      random = ~g,
      data = cs$data,
      backend = "inla",
      return_code = TRUE
    ))
    scale <- code$response_standardisation$scale
    for (key in c("sigma", "g")) {
      expect_identical(code$hyper[[key]]$meta$family, "half_normal")
      expect_equal(
        code$hyper[[key]]$meta$scale * scale,
        .five_sd(cs$modelled),
        tolerance = 1e-6,
        label = paste0("INLA ", key, " scale on the original scale, ", cs$label)
      )
    }
  }
})

test_that("a scale INLA fills beside a stated prior reads the modelled response", {
  skip_if_not_installed("INLA")
  d <- .resp_data()
  code <- suppressMessages(flexybayes(
    fixed = log(y) ~ x,
    random = ~g,
    data = d,
    backend = "inla",
    prior = fb_prior(sigma ~ uniform(lower = 0, upper = 10)),
    return_code = TRUE
  ))
  expect_equal(
    code$hyper[["g"]]$meta$scale * code$response_standardisation$scale,
    .five_sd(log(d$y)),
    tolerance = 1e-6
  )
})

# ---------------------------------------------------------------- #
# 2. A NaN or infinite response refuses by name                      #
# ---------------------------------------------------------------- #

test_that("a NaN or infinite response refuses by name before any engine runs", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  local_mocked_bindings(
    .dispatch_backend = function(...) stop("the engine was reached")
  )
  d <- .resp_data()
  d_inf <- d
  d_inf$y[5L] <- Inf
  d_nan <- d
  d_nan$y[5L] <- NaN
  d_zero <- d
  d_zero$y[5L] <- 0

  cases <- list(
    list(fixed = y ~ x, data = d_inf, shown = "row 5: Inf"),
    list(fixed = y ~ x, data = d_nan, shown = "row 5: NaN"),
    list(fixed = log(y) ~ x, data = d_zero, shown = "row 5: -Inf")
  )

  for (cs in cases) {
    for (engine in c("inla", "brms")) {
      err <- tryCatch(
        suppressMessages(flexybayes(
          fixed = cs$fixed,
          random = ~g,
          data = cs$data,
          backend = engine,
          verbose = FALSE
        )),
        error = function(e) e
      )
      expect_s3_class(err, "flexybayes_refusal_response_not_finite")
      expect_s3_class(err, "flexybayes_response_refusal")
      expect_match(conditionMessage(err), cs$shown, fixed = TRUE)
    }
  }
})

test_that("a missing response is left to na_action and reaches the engine", {
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  local_mocked_bindings(
    .dispatch_backend = function(...) stop("the engine was reached")
  )
  d <- .resp_data()
  d$y[5L] <- NA

  expect_error(
    suppressMessages(flexybayes(
      fixed = log(y) ~ x,
      random = ~g,
      data = d,
      backend = "inla",
      verbose = FALSE
    )),
    "the engine was reached",
    fixed = TRUE
  )
})
