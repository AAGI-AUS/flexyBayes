# Class guard DC-04: every route of the fitting verb carries the whole input.
#
# One model, one data set and one set of arguments are sent down every route
# that fits it (flexybayes() per row, flexybayes() aggregated,
# flexybayes_stream() and the fb_inla() pin) for each of several inputs:
# the default prior, a stated prior, a missing response and observation
# weights. A route either refuses with a classed refusal or returns the same
# fit as the others: the same prior in every engine slot (read back from the
# fit), the same number of observed responses, and the same estimates. A
# route that drops part of the input shows up as a difference. Past members:
# the streamed route passed no prior (e5527fbe), skipped the missing-response
# layer (9f0fd89f), the aggregated Gaussian route fitted its own residual
# prior (0d3b3892), fb_select() ignored weights on one route (9ca962fc).

.g04_data <- function() {
  set.seed(4)
  d <- data.frame(
    site = factor(rep(sprintf("S%02d", 1:8), each = 12)),
    trt = factor(rep(c("a", "b", "c"), 32))
  )
  u <- stats::rnorm(8, sd = 0.5)
  d$cnt <- stats::rpois(96, exp(0.5 + 0.3 * as.integer(d$trt) + u[d$site]))
  d$y <- 2 + 0.4 * as.integer(d$trt) + u[d$site] + stats::rnorm(96, sd = 0.6)
  d$w <- stats::runif(96, 0.5, 2)
  d
}

# Each route as a function of the shared arguments. A route that takes no
# argument of that name receives it all the same: dropping it is the defect.
.g04_routes <- function() {
  list(
    per_row = function(a) do.call(flexybayes, c(a, list(backend = "inla", aggregate = FALSE))),
    aggregated = function(a) do.call(flexybayes, c(a, list(backend = "inla", aggregate = TRUE))),
    stream = function(a) {
      a$source <- a$data
      a$data <- NULL
      do.call(flexybayes_stream, a)
    },
    fb_inla_pin = function(a) do.call(fb_inla, c(a, list(aggregate = FALSE)))
  )
}

.g04_run <- function(route, a) {
  tryCatch(
    suppressWarnings(suppressMessages(route(c(a, list(verbose = FALSE))))),
    error = function(e) e
  )
}

# What a fit received and returned, in comparable form.
.g04_profile <- function(fit) {
  sl <- prior_summary(fit)$slots
  sl <- sl[sl$section != "fixed", c("slot_id", "readback", "source"), drop = FALSE]
  sl <- sl[order(sl$slot_id), , drop = FALSE]
  rownames(sl) <- NULL
  s <- summary(fit)
  list(
    priors = sl,
    n_observed = s$n_observed,
    fixed = stats::setNames(s$fixed$estimate, rownames(s$fixed)),
    varcomp = stats::setNames(s$varcomp$estimate, s$varcomp$component)
  )
}

test_that("DC-04 guard: every route of the fitting verb carries the whole input", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  # The Poisson missing-response input empties the unobserved-cell table.
  # That warning is tested in test-na-action-augment.R.
  local_clean_emit_state()
  set_emit_state("missing_link_scale_warning")
  d <- .g04_data()
  d_na <- d
  d_na$cnt[5] <- NA
  d_na$y[5] <- NA
  stated <- fb_prior(sd(group = "site") ~ half_normal(scale = 2))

  inputs <- list(
    count_default_prior = list(fixed = cnt ~ trt, random = ~site, data = d, family = "poisson"),
    count_stated_prior = list(fixed = cnt ~ trt, random = ~site, data = d, family = "poisson", prior = stated),
    count_missing_response = list(fixed = cnt ~ trt, random = ~site, data = d_na, family = "poisson"),
    gaussian_stated_prior = list(fixed = y ~ trt, random = ~site, data = d, family = "gaussian", prior = stated),
    gaussian_weights = list(fixed = y ~ trt, random = ~site, data = d, family = "gaussian", weights = d$w)
  )
  routes <- .g04_routes()

  for (input_nm in names(inputs)) {
    fits <- lapply(routes, .g04_run, a = inputs[[input_nm]])
    refused <- vapply(fits, inherits, logical(1L), what = "flexybayes_refusal")
    errored <- vapply(fits, inherits, logical(1L), what = "error") & !refused
    expect_identical(names(fits)[errored], character(0),
      label = paste0(input_nm, ": routes that stopped without a refusal (",
        paste(vapply(fits[errored], conditionMessage, ""), collapse = " | "), ")")
    )
    fitted <- names(fits)[!refused & !errored]
    # Every input fits on at least the per-row route, so a comparison is made.
    expect_true("per_row" %in% fitted, label = paste0(input_nm, ": the per-row route fits"))
    if (length(fitted) < 2L) next

    ref <- .g04_profile(fits[[fitted[1L]]])
    for (route_nm in fitted[-1L]) {
      got <- .g04_profile(fits[[route_nm]])
      where <- paste0(input_nm, ", ", route_nm, " against ", fitted[1L])
      expect_identical(got$priors, ref$priors, label = paste0(where, ": the prior in every engine slot"))
      expect_identical(got$n_observed, ref$n_observed, label = paste0(where, ": observed responses"))
      expect_identical(names(got$fixed), names(ref$fixed), label = paste0(where, ": fixed-effect names"))
      expect_equal(unname(got$fixed), unname(ref$fixed), tolerance = 0.01,
        label = paste0(where, ": fixed-effect estimates")
      )
      expect_equal(unname(got$varcomp), unname(ref$varcomp), tolerance = 0.02,
        label = paste0(where, ": variance-component estimates")
      )
    }
  }
})
