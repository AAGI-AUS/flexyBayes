# Helpers for the Bayesian selection tests: one simulated dataset, one
# base fit, and a stand-in for the joint-prior compiler.
#
# The compiler (`.fb_pen_compile_brms()`) is built in a sibling unit. What
# this file provides is the SHAPE of its output -- a list of prior blocks
# plus the brms prior data frame and stanvars they compile to -- so the
# reference fit and the projection can be exercised against the interface
# rather than against an implementation. The prior rows are written with
# `brms::set_prior()` directly, which is what the compiler emits.

# CmdStan is the Stan toolchain of record on this machine, and the fits in
# this file are small enough that compile time dominates; the rstan backend
# would recompile the same programs against an older Stan.
#
# The option is set only when cmdstanr is actually reachable. A helper is
# sourced for EVERY test file, so an unconditional setting leaks into files
# that have nothing to do with this module: under `R CMD check`, where
# cmdstanr is not on the library path, it turned seven unrelated brms fits
# into "Please install the 'cmdstanr' package".
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  options(brms.backend = "cmdstanr")
}


# --- data -------------------------------------------------------------- #

# A clustered Gaussian dataset with two real candidate signals and two
# noise columns, plus one protected within-cluster covariate. The signals
# are within-cluster so that a cluster-level fold does not remove them.
fb_sel_bayes_data <- function(seed = 20260912L, J = 12L, m = 6L) {
  set.seed(seed)
  sites <- sprintf("S%02d", seq_len(J))
  d <- data.frame(site = factor(rep(sites, each = m)))
  n <- nrow(d)
  for (v in c("p1", "x1", "x2", "z1", "z2")) {
    d[[v]] <- stats::rnorm(n)
  }
  u <- stats::setNames(stats::rnorm(J, sd = 0.7), sites)
  d$y <- 1 + 0.5 * d$p1 + 1.4 * d$x1 - 1.1 * d$x2 +
    as.numeric(u[as.character(d$site)]) + stats::rnorm(n, sd = 0.5)
  d
}

# The protected model: one covariate and one grouping factor. Small
# sampler budget on purpose -- these tests are about the module's contract,
# not about the posterior's precision.
fb_sel_bayes_base <- function(d,
                              chains = 2L,
                              n_samples = 400L,
                              warmup = 400L,
                              seed = 20260912L) {
  fb_brms(
    y ~ p1 + (1 | site),
    data = d,
    chains = chains,
    n_samples = n_samples,
    warmup = warmup,
    seed = seed,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )
}

# One fit per file rather than one per test: the same posterior answers
# every question asked of it here, and a Stan compile per test would put
# the file well past its time budget.
fb_sel_bayes_cache <- new.env(parent = emptyenv())

fb_sel_bayes_memo <- function(key, value) {
  if (is.null(fb_sel_bayes_cache[[key]])) {
    assign(key, value, envir = fb_sel_bayes_cache)
  }
  get(key, envir = fb_sel_bayes_cache)
}


# --- the compiler double ------------------------------------------------- #

# `.fb_sel_test_blocks()` --- what `.fb_pen_compile_brms()` returns.
#
# Two families, and they take the two routes section 4.4 distinguishes:
#
#   regularized_horseshoe   class_wide  one row on the whole `b` class,
#                                       because brms refuses a
#                                       coefficient-specific row beside a
#                                       class-wide special prior
#   ridge                   coef_rows   one normal row per candidate
#                                       coefficient, which targets the
#                                       block exactly
#
# `par_ratio` is brms's own calibration of the global scale: p0 / (p - p0),
# the prior guess at the number of non-zero coefficients over the rest.
.fb_sel_test_blocks <- function(coef_names,
                                family = c("regularized_horseshoe", "ridge"),
                                p0 = 2,
                                scale = 0.5) {
  family <- match.arg(family)
  p <- length(coef_names)
  if (p < 2L) {
    stop("the test double needs at least two candidate coefficients.")
  }

  if (identical(family, "regularized_horseshoe")) {
    par_ratio <- p0 / (p - p0)
    prior <- brms::set_prior(
      paste0("horseshoe(par_ratio = ", format(par_ratio), ")"),
      class = "b"
    )
    route <- "class_wide"
    hyper <- list(p0 = p0, slab_scale = 2, slab_df = 4, df = 1, df_global = 1)
  } else {
    prior <- do.call(
      rbind,
      lapply(coef_names, function(nm) {
        brms::set_prior(
          paste0("normal(0, ", format(scale), ")"),
          class = "b",
          coef = nm
        )
      })
    )
    prior <- structure(
      prior,
      class = class(brms::set_prior("normal(0, 1)", class = "b"))
    )
    route <- "coef_rows"
    hyper <- list(scale = scale)
  }

  blocks <- list(list(
    id = 1L,
    family = family,
    target = "term(everything())",
    columns = seq_along(coef_names),
    coef_names = coef_names,
    protected = FALSE,
    hyperparameters = hyper,
    hyperparameter_mode = if (identical(family, "ridge")) "fixed" else "learned",
    scaling = "none",
    route = route
  ))

  list(
    blocks = blocks,
    prior = prior,
    stanvars = NULL,
    param_map = stats::setNames(
      as.list(paste0("b_", coef_names)),
      coef_names
    ),
    route = route
  )
}

# The column map the reference fit carries as an attribute. The real one
# comes from `.fb_pen_resolve()`; this is its shape.
.fb_sel_test_column_map <- function(coef_names, protected = "p1") {
  data.frame(
    column = c(protected, coef_names),
    role = c(
      rep("protected", length(protected)),
      rep("candidate", length(coef_names))
    ),
    block = c(rep(NA_integer_, length(protected)), rep(1L, length(coef_names))),
    centre = 0,
    scale = 1,
    stringsAsFactors = FALSE
  )
}
