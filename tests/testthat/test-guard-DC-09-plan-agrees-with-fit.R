# Class guard DC-09: preflight disagrees with the fitting gate.
#
# A corpus of model shapes goes through fb_plan() and through flexybayes()
# under backend = "auto" and "inla". The plan's verdict (refused, per-row
# INLA, aggregated INLA or brms) must be the route the fit takes. The brms
# emitter is replaced by a stub, so a model routed to brms is recorded
# without sampling. Every INLA route runs.
#
# The register rows behind it: the plan refused models the fit accepts
# (DEF-003, DEF-028, DEF-088), read the INLA gate as the verdict under
# auto (DEF-042), and announced the aggregated route where the gate or the
# aggregator refused (DEF-051, DEF-055).
#
# `.g09_open` lists the case and backend pairs where the two disagree
# today, with their register rows. The list must stay exact.

.g09_data <- function() {
  set.seed(9)
  d <- expand.grid(
    rowf = factor(1:6), colf = factor(1:4), env = factor(c("e1", "e2", "e3"))
  )
  d$gen <- factor(sample(sprintf("v%d", 1:8), nrow(d), replace = TRUE))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 2 + as.integer(d$env) + stats::rnorm(8L)[d$gen] + stats::rnorm(nrow(d))
  d$n <- stats::rpois(nrow(d), 3)
  d$pos <- exp(d$y / 3)
  d
}

.g09_corpus <- list(
  gaussian_iid = list(y ~ env, random = ~gen),
  poisson_iid = list(n ~ env, random = ~gen, family = "poisson"),
  poisson_continuous = list(n ~ x, random = ~gen, family = "poisson"),
  interaction = list(y ~ env, random = ~ gen + env:gen),
  dsum_residual = list(y ~ env, random = ~gen, residual = ~ dsum(~ units | env)),
  diag = list(y ~ env, random = ~ diag(env):gen),
  us = list(y ~ env, random = ~ us(env):gen),
  correlated_slope = list(y ~ x + (x | gen)),
  uncorrelated_slope = list(y ~ x + (x || gen)),
  aggregate_gaussian = list(y ~ env, random = ~gen, aggregate = TRUE),
  aggregate_continuous = list(
    n ~ x, random = ~gen, family = "poisson", aggregate = TRUE
  ),
  ar1_field = list(y ~ env, random = ~ gen + id(env):ar1(rowf):ar1(colf)),
  lognormal = list(pos ~ env, random = ~gen, family = "lognormal"),
  spline = list(y ~ env, random = ~ gen + spl(x))
)

# Pairs where the plan and the fit disagree today.
.g09_open <- stats::setNames(character(0), character(0))

.g09_plan_verdict <- function(p) {
  if (inherits(p, "error")) {
    return("refused")
  }
  if (!isTRUE(p$will_fit)) {
    return("refused")
  }
  if (identical(p$backend_chosen, "brms")) {
    return("brms")
  }
  if (identical(p$path, "aggregated_inla")) {
    return("aggregated")
  }
  "per_row"
}

.g09_fit_verdict <- function(fit) {
  if (inherits(fit, "error")) {
    return("refused")
  }
  if (inherits(fit, "g09_brms")) {
    return("brms")
  }
  if (inherits(fit, "flexybayes_aggregated")) {
    return("aggregated")
  }
  "per_row"
}

test_that("DC-09 guard: fb_plan() names the route the fit takes", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  d <- .g09_data()
  emit <- .backend_emit_fn
  local_mocked_bindings(.backend_emit_fn = function(name) {
    if (identical(name, "brms")) {
      return(function(...) structure(list(), class = "g09_brms"))
    }
    emit(name)
  })
  quiet <- function(expr) {
    tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) e)
  }
  disagree <- stats::setNames(character(0), character(0))
  verdicts <- character(0)
  for (nm in names(.g09_corpus)) {
    for (backend in c("auto", "inla")) {
      a <- c(.g09_corpus[[nm]], list(backend = backend))
      pa <- a
      pa$formula <- pa[[1L]]
      pa[[1L]] <- NULL
      plan <- .g09_plan_verdict(quiet(do.call(fb_plan, c(pa, list(data = d)))))
      fit <- .g09_fit_verdict(quiet(
        do.call(flexybayes, c(a, list(data = d, verbose = FALSE)))
      ))
      verdicts <- c(verdicts, fit)
      if (!identical(plan, fit)) {
        disagree[[paste(nm, backend)]] <- paste0("plan ", plan, ", fit ", fit)
      }
    }
  }
  # The corpus reaches every route.
  expect_setequal(unique(verdicts), c("refused", "per_row", "aggregated", "brms"))
  expect_identical(
    setdiff(names(disagree), names(.g09_open)),
    character(0),
    label = "plan and fit disagree"
  )
  expect_identical(
    setdiff(names(.g09_open), names(disagree)),
    character(0),
    label = "open-list pairs that now agree (remove them from .g09_open)"
  )
})
