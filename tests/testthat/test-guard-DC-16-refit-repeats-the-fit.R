# Class guard DC-16: a re-fit does not repeat the fit it re-issues.
#
# update() with no change must return the fit it was called on: the same
# class, route, recorded call, prior, data and estimates. The INLA fits
# below cover the routes and arguments the class's register rows broke:
# the aggregated count route (DEF-062, DEF-067), the auto-default and a
# stated prior (DEF-063, DEF-068), an explicit engine and representation
# (DEF-072), the missing-response policy (DEF-062), a known variance
# (DEF-169, DEF-175), weights, and an `omit` fit. The brms block runs
# review_code then proceed() (DEF-176) and sits in the long tier.

.g16_data <- function() {
  set.seed(16)
  d <- data.frame(
    g = factor(rep(seq_len(6L), each = 10L)),
    f = factor(rep(c("a", "b"), 30L))
  )
  d$y <- 1 + 0.5 * (d$f == "b") + stats::rnorm(6L)[d$g] + stats::rnorm(60L)
  d$n <- stats::rpois(60L, 3 + (d$f == "b"))
  d$v <- stats::runif(60L, 0.1, 0.4)
  d$w <- stats::runif(60L, 0.5, 2)
  d
}

# Each fit names its data, as a user's call does, so the recorded data
# name is a name.
.g16_fits <- function(d) {
  list(
    gaussian_default = list(y ~ f, random = ~g, data = quote(d)),
    gaussian_prior_vc_sd = list(
      y ~ f,
      random = ~g,
      data = quote(d),
      prior_vc_sd = 2
    ),
    gaussian_stated_prior = list(
      y ~ f,
      random = ~g,
      data = quote(d),
      prior = fb_prior(
        sigma ~ uniform(lower = 0, upper = 20),
        sd(group = "g") ~ uniform(lower = 0, upper = 20)
      )
    ),
    gaussian_weights_inla = list(
      y ~ f,
      random = ~g,
      data = quote(d),
      weights = d$w,
      backend = "inla"
    ),
    gaussian_known_variance = list(
      y ~ f,
      random = ~g,
      data = quote(d),
      known_variance = d$v
    ),
    gaussian_omit = list(
      y ~ f,
      random = ~g,
      data = quote(d_na),
      na_action = "omit"
    ),
    poisson_aggregated = list(
      n ~ f,
      random = ~g,
      data = quote(d),
      family = "poisson"
    ),
    poisson_per_row = list(
      n ~ f,
      random = ~g,
      data = quote(d),
      family = "poisson",
      aggregate = FALSE
    )
  )
}

# What a re-fit must repeat, read from one fit.
.g16_record <- function(fit) {
  ci <- fit$extras$call_info
  ci$verbose <- NULL
  s <- suppressMessages(summary(fit))
  list(
    class = class(fit),
    route = fit$extras$manifest$route$route,
    call = ci,
    prior = utils::capture.output(print(prior_summary(fit))),
    data = fit$glm$data %||% fit$data,
    fixed = s$fixed$estimate,
    varcomp = s$varcomp$estimate
  )
}

.g16_compare <- function(label, a, b) {
  out <- character(0)
  for (nm in c("class", "route", "call", "prior", "data")) {
    if (!identical(a[[nm]], b[[nm]])) {
      out <- c(out, paste0(label, ": ", nm, " differs"))
    }
  }
  # Identical INLA calls on these fits differ by up to 2.3e-5 relative.
  for (nm in c("fixed", "varcomp")) {
    same <- length(a[[nm]]) == length(b[[nm]]) &&
      isTRUE(all.equal(a[[nm]], b[[nm]], tolerance = 1e-4))
    if (!same) {
      out <- c(
        out,
        paste0(
          label,
          ": ",
          nm,
          " estimates differ (",
          paste(signif(a[[nm]], 7L), collapse = ", "),
          " against ",
          paste(signif(b[[nm]], 7L), collapse = ", "),
          ")"
        )
      )
    }
  }
  out
}

test_that("DC-16 guard: update() with no change repeats every INLA fit", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  d <- .g16_data()
  d_na <- d
  d_na$y[c(3L, 17L)] <- NA
  specs <- .g16_fits(d)
  diffs <- character(0)
  routes <- character(0)
  for (nm in names(specs)) {
    fit <- suppressWarnings(suppressMessages(
      do.call(flexybayes, c(specs[[nm]], list(verbose = FALSE)))
    ))
    refit <- suppressWarnings(suppressMessages(update(fit, verbose = FALSE)))
    a <- .g16_record(fit)
    routes <- c(routes, a$route)
    diffs <- c(diffs, .g16_compare(nm, a, .g16_record(refit)))
  }
  # The sweep reaches both the per-row and the aggregated count route.
  expect_true(all(c("aggregated_count", "per_row") %in% routes))
  expect_identical(diffs, character(0))
})

test_that("DC-16 guard: proceed() and update() repeat a brms fit", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_long_tier()
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE
  )
  d <- .g16_data()
  args <- list(
    y ~ f,
    random = ~g,
    data = d,
    known_variance = d$v,
    backend = "brms",
    chains = 1L,
    n_samples = 200L,
    warmup = 200L,
    seed = 16L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )
  fit <- suppressWarnings(suppressMessages(do.call(flexybayes, args)))
  review <- suppressWarnings(suppressMessages(
    do.call(flexybayes, c(args, list(review_code = TRUE)))
  ))
  proceeded <- suppressWarnings(suppressMessages(proceed(review)))
  refit <- suppressWarnings(suppressMessages(update(fit, verbose = FALSE)))
  a <- .g16_record(fit)
  diffs <- c(
    .g16_compare("proceed", a, .g16_record(proceeded)),
    .g16_compare("update", a, .g16_record(refit))
  )
  expect_identical(diffs, character(0))
})
