# data-raw/inla_spline_constrained.R
#
# An INLA fit of y ~ trt, random = ~ site + spl(x) whose stored precision at
# one or more hyperparameter configurations is singular to working
# precision, with the posterior moments INLA's own sampler,
# INLA::inla.posterior.sample(), gives for the same fit. They are the oracle
# for tests/testthat/test-inla-config-constrained.R.
#
# Run from the package root:
#
#   Rscript --vanilla data-raw/inla_spline_constrained.R
#
# The rw2 smooth and the intercept share the constant direction, which only
# the smooth's sum-to-zero constraint fixes. Whether the Cholesky factor of a
# configuration's precision fails depends on rounding inside the INLA run,
# so the script refits on new simulated data until one fails, and stores
# that fit.
pkgload::load_all(".", quiet = TRUE)
options(
  flexyBayes.inla_threads = "1:1",
  flexyBayes.silence_default_prior_note = TRUE
)

not_positive_definite <- function(q) {
  q <- Matrix::forceSymmetric(q, uplo = "U")
  out <- tryCatch(
    suppressWarnings(Matrix::Cholesky(q, LDL = FALSE)),
    error = function(e) e
  )
  inherits(out, "error") && grepl("not positive", conditionMessage(out))
}

simulate <- function(seed) {
  set.seed(seed)
  d <- data.frame(
    site = factor(rep(1:6, each = 12)),
    trt = factor(rep(c("a", "b", "c"), 24)),
    x = stats::rnorm(72)
  )
  d$y <- 1 + c(0, 0.5, -0.3)[d$trt] + 0.4 * d$x +
    stats::rnorm(6)[d$site] + stats::rnorm(72)
  d
}

seed <- 9L
repeat {
  d <- simulate(seed)
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ trt,
    random = ~ site + spl(x),
    data = d,
    backend = "inla",
    verbose = FALSE
  )))
  failing <- which(vapply(
    fit$inla$misc$configs$config,
    function(k) not_positive_definite(k$Q),
    logical(1L)
  ))
  if (length(failing) > 0L) {
    break
  }
  seed <- seed + 1L
}

# INLA's sampler, on the standardised scale the INLA fit is built on.
n_draws <- 100000L
draws <- INLA::inla.posterior.sample(
  n_draws,
  fit$inla,
  seed = 20261008L,
  num.threads = "1:1"
)
latent_rows <- rownames(draws[[1L]]$latent)
fixed_rows <- paste0(rownames(fit$inla$summary.fixed), ":1")
spline_rows <- grep("^x:", latent_rows, value = TRUE)
fixed <- t(vapply(
  draws,
  function(z) z$latent[fixed_rows, 1L],
  numeric(length(fixed_rows))
))
spline <- t(vapply(
  draws,
  function(z) z$latent[spline_rows, 1L],
  numeric(length(spline_rows))
))

# The parts of the fit that vcov(), predict() and the draws read; the
# INLA call arguments, log and marginal densities are dropped for size.
fit$inla$.args <- NULL
fit$inla$logfile <- NULL
fit$inla$marginals.random <- NULL
fit$extras$formula <- NULL

inla_spline_constrained <- list(
  stamp = list(
    date = as.character(Sys.Date()),
    inla = as.character(utils::packageVersion("INLA")),
    r = R.version.string,
    script = "data-raw/inla_spline_constrained.R"
  ),
  seed = seed,
  failing_configs = failing,
  data = d,
  fit = fit,
  oracle = list(
    n_draws = n_draws,
    fixed_mean = colMeans(fixed),
    fixed_cov = stats::cov(fixed),
    spline_mean = colMeans(spline),
    spline_sd = apply(spline, 2L, stats::sd),
    spline_sum_max = max(abs(rowSums(spline)))
  )
)

saveRDS(
  inla_spline_constrained,
  "tests/testthat/fixtures/inla_spline_constrained.rds",
  compress = "xz"
)
cat(
  "written tests/testthat/fixtures/inla_spline_constrained.rds: seed", seed,
  "failing configurations", paste(failing, collapse = ", "), "\n"
)
