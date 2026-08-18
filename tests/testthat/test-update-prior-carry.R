# =============================================================================
# The prior a re-fit runs under.
#
# The auto-default bounded uniform on the SD scale fires only when the call
# supplies neither a `prior` nor a `prior_vc_sd`. update() re-issues the
# recorded `prior_vc_sd`, so the gate saw an argument that was present and
# the default never fired on a re-fit: the model fell through to the legacy
# scalar bridge and the engine's own hyperprior. An identity update() --
# nothing changed -- therefore returned a different model, and said so
# nowhere. Measured before the repair, on the fixture below:
#
#   route                 sd component   original   identity update()
#   per-row INLA          sd_Block        2.2785     1.1443
#   aggregated Gaussian   tau             1.5581     0.9798
#   brms                  sd_Block        2.1877     (engine default)
#
# with `summary(fit)$varcomp$prior` reading `uniform(lower=0, upper=...)`
# before and `engine default` after.
#
# The repair re-issues the resolved prior the fit carries
# (`fit$extras$fb_terms$priors`) rather than inferring the default's
# applicability from argument missingness. These tests hold the three
# active routes to the invariant, and hold the two override spellings to
# outranking the record.
# =============================================================================

suppressPackageStartupMessages(library(testthat))


.upc_data <- function(seed = 11L, n = 72L) {
  set.seed(seed)
  d <- data.frame(
    Trt = factor(rep(c("a", "b", "c"), length.out = n)),
    Block = factor(rep(paste0("B", 1:6), length.out = n)),
    # A continuous covariate makes every design cell distinct, so the
    # Gaussian route stays per-row instead of auto-routing to the
    # aggregated emit. Dropped in the fixtures that want the aggregated
    # representation.
    x = stats::runif(n)
  )
  b_block <- stats::rnorm(6L, sd = 1.0)[as.integer(d$Block)]
  d$y <- 10 + 1.5 * (d$Trt == "b") - 0.7 * (d$Trt == "c") +
    0.5 * d$x + b_block + stats::rnorm(n, sd = 0.8)
  d
}

.upc_silence <- function() {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    flexyBayes.silence_auto_fallback_note = TRUE,
    flexyBayes.silence_auto_inla_missing_note = TRUE,
    .local_envir = parent.frame()
  )
}

# summary() prints as a side effect and returns invisibly, so every read of
# a slot goes through this.
.upc_summary <- function(fit) {
  out <- NULL
  invisible(utils::capture.output(out <- summary(fit)))
  out
}

# The variance components as a named vector, so the two fits are compared
# component by name rather than by row position.
.upc_varcomp <- function(sm) {
  stats::setNames(sm$varcomp$estimate, sm$varcomp$component)
}

.upc_prior_cells <- function(sm) {
  stats::setNames(sm$varcomp$prior, sm$varcomp$component)
}


# ---------------------------------------------------------------- #
# 1. The three active routes                                        #
# ---------------------------------------------------------------- #

test_that("an identity update() on a per-row INLA fit keeps its prior", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .upc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ Trt + x, random = ~Block, data = .upc_data(), backend = "inla",
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  # The premise: this fixture fires the auto-default, which is the prior
  # the re-fit has to reproduce.
  expect_identical(prior_summary(fit)$kind, "fb_prior")
  expect_identical(prior_summary(fit)$default_origin, "auto")

  refit <- suppressMessages(stats::update(
    fit, verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_identical(prior_summary(refit)$kind, "fb_prior")
  expect_identical(prior_summary(refit)$default_origin, "auto")
  expect_identical(
    prior_summary(refit)$default_scale,
    prior_summary(fit)$default_scale
  )

  sm0 <- .upc_summary(fit)
  sm1 <- .upc_summary(refit)
  # Cell for cell: the prior column of the variance-component table is the
  # surface a user reads, and it is what changed.
  expect_identical(.upc_prior_cells(sm1), .upc_prior_cells(sm0))
  expect_true(all(grepl("^uniform\\(", .upc_prior_cells(sm1))))
  # The band is set well inside the defect it guards and outside the
  # engine's own reproducibility. Two runs of the same INLA model on this
  # machine have differed in the third decimal of sd_Block (2.2779 against
  # 2.2837, a quarter of a per cent), while the defect halved it -- 2.2785
  # to 1.1443.
  v0 <- .upc_varcomp(sm0)
  v1 <- .upc_varcomp(sm1)
  expect_identical(names(v1), names(v0))
  expect_equal(unname(v1), unname(v0), tolerance = 0.05)
})

test_that("an identity update() on an aggregated fit keeps its prior", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .upc_silence()
  d <- .upc_data()
  d$x <- NULL
  fit <- suppressMessages(flexybayes(
    y ~ Trt, random = ~Block, data = d, verbose = FALSE, mcmc_verbose = FALSE
  ))
  # The premise, asserted rather than assumed: this is the aggregated
  # representation, which returns the same unified summary object as
  # every other route.
  expect_s3_class(fit, "flexybayes_aggregated")
  expect_identical(prior_summary(fit)$default_origin, "auto")

  refit <- suppressMessages(stats::update(
    fit, verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_s3_class(refit, "flexybayes_aggregated")
  expect_identical(prior_summary(refit)$kind, "fb_prior")
  expect_identical(prior_summary(refit)$default_origin, "auto")
  # The resolved prior object itself, not a summary of it.
  expect_equal(
    refit$extras$fb_terms$priors,
    fit$extras$fb_terms$priors
  )

  sm0 <- .upc_summary(fit)
  sm1 <- .upc_summary(refit)
  # Read off the unified table, component by name, exactly as on the
  # per-row route. Same band, and for the same reason: the defect took
  # tau from 1.558 to 0.980.
  v0 <- .upc_varcomp(sm0)
  v1 <- .upc_varcomp(sm1)
  expect_identical(names(v1), names(v0))
  expect_equal(unname(v1), unname(v0), tolerance = 0.05)
  expect_identical(.upc_prior_cells(sm1), .upc_prior_cells(sm0))
  # The raw aggregated pieces are still where they were, so nothing that
  # read them through `$extras$summary` moved.
  expect_equal(
    refit$extras$summary$sigma_means,
    fit$extras$summary$sigma_means,
    tolerance = 0.05
  )
  expect_equal(
    refit$extras$summary$tau_means,
    fit$extras$summary$tau_means,
    tolerance = 0.05
  )
})

test_that("an identity update() on a brms fit keeps its prior", {
  skip_if_not_installed("brms")
  skip_on_cran()
  .upc_silence()
  d <- .upc_data()
  d$x <- NULL
  # Stan's own sampler warnings on a 500-draw fixture are not what this
  # file is measuring, and they would otherwise be the suite's only
  # warnings.
  fit <- suppressMessages(suppressWarnings(flexybayes(
    y ~ Trt, random = ~Block, data = d, backend = "brms",
    n_samples = 500L, warmup = 500L, chains = 2L, seed = 4321L,
    verbose = FALSE, mcmc_verbose = FALSE
  )))
  expect_identical(prior_summary(fit)$default_origin, "auto")

  # `backend` is not part of the recorded call, so a bare update() of a
  # brms fit re-routes to the auto backend. The engine is named here so
  # this test measures the prior carry and not the routing.
  refit <- suppressMessages(suppressWarnings(stats::update(
    fit, backend = "brms", verbose = FALSE, mcmc_verbose = FALSE
  )))
  expect_s3_class(refit, "flexybayes_brms")
  expect_identical(prior_summary(refit)$kind, "fb_prior")
  expect_identical(prior_summary(refit)$default_origin, "auto")

  sm0 <- .upc_summary(fit)
  sm1 <- .upc_summary(refit)
  expect_identical(.upc_prior_cells(sm1), .upc_prior_cells(sm0))
  # The recorded seed is re-issued too, so two runs of the same model
  # under the same prior agree to sampler noise at worst.
  v0 <- .upc_varcomp(sm0)
  v1 <- .upc_varcomp(sm1)
  expect_identical(names(v1), names(v0))
  expect_equal(unname(v1), unname(v0), tolerance = 0.05)
})


# ---------------------------------------------------------------- #
# 2. A user-supplied prior survives the re-fit                      #
# ---------------------------------------------------------------- #

test_that("a user fb_prior() is the prior the re-fit runs under", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .upc_silence()
  # The bound sits clear of the posterior mass. A uniform whose upper
  # bound cuts into it (upper = 3 against an sd_Block near 2.27) makes
  # the inla binary abort, which is an upstream matter and not what this
  # file is measuring.
  priors <- fb_prior(
    sigma ~ uniform(lower = 0, upper = 6),
    sd(group = "Block") ~ uniform(lower = 0, upper = 6)
  )
  fit <- suppressMessages(flexybayes(
    y ~ Trt + x, random = ~Block, data = .upc_data(), backend = "inla",
    prior = priors, verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_identical(prior_summary(fit)$default_origin, "user")

  refit <- suppressMessages(stats::update(
    fit, verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_identical(prior_summary(refit)$kind, "fb_prior")
  expect_identical(prior_summary(refit)$default_origin, "user")
  expect_equal(refit$extras$fb_terms$priors, priors)
  expect_identical(
    .upc_prior_cells(.upc_summary(refit)),
    .upc_prior_cells(.upc_summary(fit))
  )
})


# ---------------------------------------------------------------- #
# 3. What the caller writes outranks the record                     #
# ---------------------------------------------------------------- #

test_that("an explicit prior_vc_sd in update() overrides the recorded prior", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .upc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ Trt + x, random = ~Block, data = .upc_data(), backend = "inla",
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_identical(prior_summary(fit)$default_origin, "auto")

  refit <- suppressMessages(stats::update(
    fit, prior_vc_sd = 3, verbose = FALSE, mcmc_verbose = FALSE
  ))
  # The scalar spelling asks for the legacy bridge. Re-issuing the
  # recorded prior alongside it would let the record win over the number
  # the caller just wrote.
  ps <- prior_summary(refit)
  expect_identical(ps$kind, "legacy_scalar")
  expect_identical(ps$vc_sd, 3)
  expect_false(inherits(refit$extras$fb_terms$priors, "fb_prior"))
})

test_that("an explicit prior in update() overrides the recorded prior", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .upc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ Trt + x, random = ~Block, data = .upc_data(), backend = "inla",
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  new_prior <- fb_prior(
    sigma ~ uniform(lower = 0, upper = 12),
    sd(group = "Block") ~ uniform(lower = 0, upper = 12)
  )
  refit <- suppressMessages(stats::update(
    fit, prior = new_prior, verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_identical(prior_summary(refit)$default_origin, "user")
  expect_equal(refit$extras$fb_terms$priors, new_prior)
  expect_true(all(grepl(
    "upper=12", .upc_prior_cells(.upc_summary(refit)),
    fixed = TRUE
  )))
})


# ---------------------------------------------------------------- #
# 4. The two records that are not an fb_prior                       #
# ---------------------------------------------------------------- #

test_that("a legacy-scalar fit re-fits on its recorded scalars", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .upc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ Trt + x, random = ~Block, data = .upc_data(), backend = "inla",
    prior_vc_sd = 2, verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_identical(prior_summary(fit)$kind, "legacy_scalar")

  refit <- suppressMessages(stats::update(
    fit, verbose = FALSE, mcmc_verbose = FALSE
  ))
  ps <- prior_summary(refit)
  expect_identical(ps$kind, "legacy_scalar")
  expect_identical(ps$vc_sd, 2)
})

test_that("a fit carrying no prior record re-fits as it did before", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .upc_silence()
  fit <- suppressMessages(flexybayes(
    y ~ Trt + x, random = ~Block, data = .upc_data(), backend = "inla",
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  # A pre-IR object: the record the repair reads is simply absent, and
  # the re-fit falls back to the recorded scalars rather than refusing.
  fit$extras$fb_terms$priors <- NULL
  expect_identical(prior_summary(fit)$kind, "no_prior_recorded")

  refit <- suppressMessages(stats::update(
    fit, verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_s3_class(refit, "flexybayes")
  expect_identical(prior_summary(refit)$kind, "legacy_scalar")
})


# ---------------------------------------------------------------- #
# 5. A policy re-fires, a bespoke prior carries                     #
# ---------------------------------------------------------------- #
#
# An update() that ADDS a random term exposes the difference between the
# two things the prior record can hold. Carrying the auto-default's
# object re-applies yesterday's output to today's model, and the added
# term falls to the engine while its siblings keep the bounded uniform --
# a mixed state in one table, visible on the vignette's own first page.
# The auto-default is a policy, so it is re-fired over the updated model
# instead. A user's fb_prior() names its terms and is carried verbatim,
# and a term it never named keeps following the engine, which is what it
# did on the first fit too.

# The fixture the section needs: a second grouping factor to add.
.upc_data_two_groups <- function(seed = 11L, n = 72L) {
  d <- .upc_data(seed = seed, n = n)
  d$Site <- factor(rep(paste0("S", 1:4), length.out = nrow(d)))
  d
}

test_that("adding a term under the auto-default priors every component", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .upc_silence()
  d <- .upc_data_two_groups()

  fit <- suppressMessages(flexybayes(
    y ~ Trt + x, random = ~Block, data = d, backend = "inla",
    aggregate = FALSE, verbose = FALSE, mcmc_verbose = FALSE
  ))
  ps <- prior_summary(fit)
  expect_identical(ps$kind, "fb_prior")
  expect_identical(ps$default_origin, "auto")

  refit <- suppressMessages(stats::update(fit, random = ~ Block + Site))
  cells <- .upc_prior_cells(.upc_summary(refit))

  # The added term is present, and priored like the rest.
  expect_true(all(c("sigma", "sd_Block", "sd_Site") %in% names(cells)))
  expect_true(all(grepl("^uniform\\(", cells)))
  expect_false(any(cells %in% "engine default"))
  # One bound, shared: the policy read it off the same data for every
  # component, so a mixed table is what this test exists to rule out.
  expect_identical(length(unique(unname(cells))), 1L)
  # And the re-fit reports the policy as the policy, not as a user prior.
  expect_identical(prior_summary(refit)$default_origin, "auto")
})

test_that("adding a term under a user prior leaves it to the engine", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .upc_silence()
  d <- .upc_data_two_groups()

  user_prior <- fb_prior(
    sigma ~ uniform(0, 10),
    sd(group = "Block") ~ uniform(0, 10)
  )
  fit <- suppressMessages(flexybayes(
    y ~ Trt + x, random = ~Block, data = d, backend = "inla",
    aggregate = FALSE, prior = user_prior,
    verbose = FALSE, mcmc_verbose = FALSE
  ))
  expect_identical(prior_summary(fit)$default_origin, "user")

  refit <- suppressMessages(stats::update(fit, random = ~ Block + Site))
  cells <- .upc_prior_cells(.upc_summary(refit))

  # The named terms keep the user's words, verbatim.
  expect_match(cells[["sigma"]], "^uniform\\(")
  expect_match(cells[["sd_Block"]], "^uniform\\(")
  expect_identical(unname(cells[["sigma"]]), unname(cells[["sd_Block"]]))
  # The term the user never named follows the engine, exactly as it would
  # have on a first fit that added it.
  expect_identical(unname(cells[["sd_Site"]]), "engine default")
  # The prior is still theirs, not re-decided into a policy.
  expect_identical(prior_summary(refit)$default_origin, "user")
})
