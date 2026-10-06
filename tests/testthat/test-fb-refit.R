# =============================================================================
# fb_refit(): the same model refitted on new data, with added fixed columns,
# extra brms prior rows and Stan code, a formula hook, and a reduced readout.
#
# Each carried setting is read back from the refit's own record (the call
# record, the parsed model and the model manifest), never assumed.
# =============================================================================

fbrf_quiet <- function(env = parent.frame()) {
  withr::local_options(
    list(
      flexyBayes.inla_threads = "1:1",
      flexyBayes.silence_default_prior_note = TRUE,
      flexyBayes.silence_uniform_inla_approx = TRUE,
      flexyBayes.silence_high_missingness_warning = TRUE
    ),
    .local_envir = env
  )
}

fbrf_data <- function(seed = 20261006L, n_site = 12L, per = 10L) {
  withr::with_seed(seed, {
    n <- n_site * per
    d <- data.frame(
      site = factor(rep(sprintf("S%02d", seq_len(n_site)), each = per)),
      x = stats::rnorm(n),
      z = stats::rnorm(n),
      o = stats::runif(n),
      v = stats::runif(n, 0.05, 0.2),
      w = rep(c(1, 0.6), length.out = n)
    )
    u <- stats::rnorm(n_site, sd = 0.8)
    d$y <- 2 + 0.5 * d$x + 0.3 * d$z + d$o + u[as.integer(d$site)] +
      stats::rnorm(n, sd = 0.5)
    d
  })
}

# The posterior medians of the fixed effects and hyperparameters, by name.
fbrf_medians <- function(fit) {
  c(
    stats::setNames(
      fit$inla$summary.fixed[["0.5quant"]],
      rownames(fit$inla$summary.fixed)
    ),
    stats::setNames(
      fit$inla$summary.hyperpar[["0.5quant"]],
      rownames(fit$inla$summary.hyperpar)
    )
  )
}

# The largest relative difference between two sets of medians.
fbrf_rel <- function(a, b) {
  max(abs(a - b) / abs(b))
}

# The value of `code`, and the arguments of the first INLA call it makes.
fbrf_capture <- function(code) {
  seen <- new.env()
  seen$args <- NULL
  real <- .inla_call
  value <- testthat::with_mocked_bindings(
    code,
    .inla_call = function(args) {
      if (is.null(seen$args)) {
        seen$args <- args
      }
      real(args)
    }
  )
  list(value = value, args = seen$args)
}

# Two INLA calls with the same formula, data, likelihood, priors and
# settings.
fbrf_expect_same_call <- function(a, b) {
  expect_identical(sort(names(a)), sort(names(b)))
  expect_identical(deparse(a$formula), deparse(b$formula))
  for (nm in setdiff(names(a), "formula")) {
    expect_identical(a[[nm]], b[[nm]], info = nm)
  }
}

# Identical flexybayes() calls on these fixtures return medians that differ
# by up to 2.4e-6 relative (3.7e-4 with a known variance), measured over six
# repeats on INLA at "1:1"; the tolerances below are about ten times that.
FBRF_TOL <- 5e-5
FBRF_TOL_KNOWN <- 5e-3

# The prior each engine slot was sent, and what INLA reports it received.
fbrf_slots <- function(fit) {
  s <- fit$extras$manifest$slots
  s <- s[s$section != "fixed", c("slot_id", "emitted", "readback"),
         drop = FALSE]
  s[order(s$slot_id), , drop = FALSE]
}


test_that("fb_refit(fit) with nothing changed reproduces the fit", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbrf_quiet()
  d <- fbrf_data()
  original <- fbrf_capture(
    flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
               verbose = FALSE)
  )
  fit <- original$value
  again <- fbrf_capture(fb_refit(fit))
  refit <- again$value

  expect_s3_class(refit, "flexybayes_inla")
  expect_identical(refit$extras$call_info$fixed, fit$extras$call_info$fixed)
  fbrf_expect_same_call(again$args, original$args)
  a <- fbrf_medians(fit)
  b <- fbrf_medians(refit)
  expect_identical(names(b), names(a))
  expect_lt(fbrf_rel(b, a), FBRF_TOL)
  expect_identical(fbrf_slots(refit), fbrf_slots(fit))
})


test_that("add_fixed refits the model flexybayes() fits with the column written in", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbrf_quiet()
  d <- fbrf_data()
  fit <- flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
                    verbose = FALSE)
  direct <- fbrf_capture(
    flexybayes(y ~ x + z, random = ~site, data = d, backend = "inla",
               verbose = FALSE)
  )
  refit <- fbrf_capture(fb_refit(fit, add_fixed = ~z))
  by_formula <- refit$value
  by_label <- fb_refit(fit, add_fixed = "z")

  fbrf_expect_same_call(refit$args, direct$args)

  expect_identical(
    deparse(by_formula$extras$call_info$fixed),
    deparse(y ~ x + z)
  )
  expect_identical(rownames(by_formula$inla$summary.fixed),
                   c("(Intercept)", "x", "z"))
  expect_lt(fbrf_rel(fbrf_medians(by_formula), fbrf_medians(direct$value)),
            FBRF_TOL)
  expect_lt(fbrf_rel(fbrf_medians(by_label), fbrf_medians(direct$value)),
            FBRF_TOL)
})


test_that("the refit carries priors, weights, offsets, known variance and engine settings", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbrf_quiet()
  d <- fbrf_data()
  priors <- fb_prior(
    sigma ~ uniform(lower = 0, upper = 20),
    sd(group = "site") ~ uniform(lower = 0, upper = 20)
  )
  fit <- suppressMessages(flexybayes(
    y ~ x + offset(o), random = ~site, data = d, backend = "inla",
    prior = priors, weights = d$w, known_variance = d$v, na_action = "omit",
    verbose = FALSE
  ))
  refit <- suppressMessages(fb_refit(fit, add_fixed = ~z))

  # Priors: the user's fb_prior() is the one the refit ran under, and every
  # variance slot was sent, and received, the prior the fit's slot was.
  expect_equal(refit$extras$fb_terms$priors, priors)
  expect_identical(prior_summary(refit)$default_origin, "user")
  expect_identical(fbrf_slots(refit), fbrf_slots(fit))
  # Weights, offsets and the known variance, read off the refit's parsed
  # model.
  expect_equal(.fb_ir_weights(refit$fb), .fb_ir_weights(fit$fb))
  expect_identical(refit$fb$offsets, "offset(o)")
  expect_equal(
    .fb_known_variance_term(refit$fb)$values,
    .fb_known_variance_term(fit$fb)$values
  )
  # Engine and engine settings, read off the refit's call record.
  for (field in c("backend", "aggregate", "na_action", "family", "link",
                  "random", "residual", "control", "seed")) {
    expect_identical(
      refit$extras$call_info[[field]],
      fit$extras$call_info[[field]],
      info = field
    )
  }
  expect_identical(
    deparse(refit$extras$call_info$fixed),
    deparse(y ~ x + offset(o) + z)
  )

  # A package default prior is rebuilt over the refitted model, and on the
  # same data every slot gets the prior the fit's slot got.
  fit_default <- flexybayes(y ~ x, random = ~site, data = d,
                            backend = "inla", verbose = FALSE)
  refit_default <- fb_refit(fit_default, add_fixed = ~z)
  expect_identical(fbrf_slots(refit_default), fbrf_slots(fit_default))
  expect_identical(prior_summary(refit_default)$default_origin, "auto")
})


test_that("weights and a known variance follow new rows only when they can", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbrf_quiet()
  d <- fbrf_data()
  fit <- flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
                    weights = d$w, verbose = FALSE)
  rows <- 1:100
  expect_error(
    fb_refit(fit, data = d[rows, ]),
    class = "flexybayes_refusal_refit_not_reconstructable"
  )
  refit <- fb_refit(fit, data = d[rows, ], weights = d$w[rows])
  expect_equal(.fb_ir_weights(refit$fb), d$w[rows])

  kv <- flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
                   known_variance = d$v, verbose = FALSE)
  # The fit's own data carries the known variance in a column; other data
  # does not, and the refit refuses before fitting.
  expect_error(
    fb_refit(kv, data = d[rows, ]),
    class = "flexybayes_refusal_known_variance_not_supported"
  )
  refit_kv <- fbrf_capture(fb_refit(kv, data = kv$data[rows, ]))
  expect_equal(.fb_known_variance_term(refit_kv$value$fb)$values, d$v[rows])
  refit_kv2 <- fbrf_capture(
    fb_refit(kv, data = d[rows, ], known_variance = d$v[rows])
  )
  fbrf_expect_same_call(refit_kv2$args, refit_kv$args)
  expect_lt(
    fbrf_rel(fbrf_medians(refit_kv2$value), fbrf_medians(refit_kv$value)),
    FBRF_TOL_KNOWN
  )
})


test_that("fb_refit() refuses a fit it cannot rebuild and arguments that change the model", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbrf_quiet()
  d <- fbrf_data()
  fit <- flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
                    verbose = FALSE)
  bare <- fit
  bare$extras$call_info$fixed <- NULL
  expect_error(
    fb_refit(bare),
    class = "flexybayes_refusal_refit_not_reconstructable"
  )
  expect_error(fb_refit(fit, random = ~x), "update\\(\\)")
  expect_error(fb_refit(fit, add_fixed = y ~ z), "one-sided")
  expect_error(fb_refit(fit, readout = "everything"))
  expect_error(
    fb_refit(fit, formula_hook = function(formula, prior) NULL),
    class = "flexybayes_refusal_refit_hooks_require_brms"
  )
})


test_that("formula_hook, extra_prior and stanvars reach the brms model", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  fbrf_quiet()
  d <- fbrf_data()
  fit <- flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
                    verbose = FALSE)
  seen <- new.env()
  hook <- function(formula, prior) {
    seen$formula <- formula
    seen$prior <- prior
    prior$prior[prior$class == "sigma"] <- "normal(0, 3)"
    list(formula = brms::bf(formula, decomp = "QR"), prior = prior)
  }
  code <- fb_refit(
    fit,
    add_fixed = ~z,
    extra_prior = brms::set_prior("normal(0, 0.25)", class = "b", coef = "z"),
    stanvars = brms::stanvar(scode = "real fbrefit_marker = 1;",
                             block = "tdata"),
    formula_hook = hook,
    backend = "brms",
    return_code = TRUE,
    verbose = FALSE
  )

  # The hook receives the linear formula the core built from the fit's
  # model and the added column.
  expect_identical(
    deparse(seen$formula),
    deparse(y ~ x + z + (1 | site))
  )
  expect_s3_class(seen$prior, "brmsprior")
  expect_type(code, "character")
  # The hook's formula: the QR decomposition of the design.
  expect_match(code, "XQ", fixed = TRUE)
  # The hook's prior row.
  expect_match(code, "normal_lpdf(sigma | 0, 3)", fixed = TRUE)
  # The extra prior row on the added column, on its QR coefficient.
  expect_match(code, "normal_lpdf\\(bQ\\[[0-9]+\\] \\| 0, 0\\.25\\)")
  # The Stan code.
  expect_match(code, "real fbrefit_marker = 1;", fixed = TRUE)

  # Without a hook the core's own linear formula is fitted.
  plain <- fb_refit(fit, add_fixed = ~z, backend = "brms",
                    extra_prior = brms::set_prior("normal(0, 0.25)",
                                                  class = "b", coef = "z"),
                    return_code = TRUE, verbose = FALSE)
  expect_false(grepl("XQ", plain, fixed = TRUE))
  expect_match(plain, "normal_lpdf\\(b\\[[0-9]+\\] \\| 0, 0\\.25\\)")
})


test_that("readout asks INLA for that readout only and returns the same estimates", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  fbrf_quiet()
  d <- fbrf_data()
  fit <- flexybayes(y ~ x, random = ~site, data = d, backend = "inla",
                    verbose = FALSE)
  fd <- fit$data
  fd$y[fd$site %in% c("S01", "S02", "S03")] <- NA
  seen <- new.env()
  seen$calls <- list()
  real <- .inla_call
  testthat::local_mocked_bindings(
    .inla_call = function(args) {
      seen$calls[[length(seen$calls) + 1L]] <- args
      real(args)
    }
  )

  full <- fb_refit(fit, data = fd, na_action = "augment")
  expect_identical(
    seen$calls[[1L]]$control.compute,
    list(config = TRUE, return.marginals = TRUE, dic = TRUE, waic = TRUE)
  )
  for (readout in c("varcomp", "fitted", "draws")) {
    seen$calls <- list()
    lean <- fb_refit(fit, data = fd, na_action = "augment", readout = readout)
    expect_identical(
      seen$calls[[1L]]$control.compute,
      list(
        config = identical(readout, "draws"),
        return.marginals = FALSE,
        dic = FALSE,
        waic = FALSE
      ),
      info = readout
    )
    expect_null(lean$inla$marginals.fixed)
    expect_lt(fbrf_rel(fbrf_medians(lean), fbrf_medians(full)), FBRF_TOL)
  }
  expect_null(.fb_fold_refit_active())
})
