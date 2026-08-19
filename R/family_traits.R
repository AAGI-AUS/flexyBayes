# family_traits -- per-engine facts about response families, one table each.
#
# A family trait is a fact about an engine, not about flexyBayes: whether
# brms parameterises a residual `sigma` for a family, and which
# hyperparameter keyword INLA's likelihood declares. Each fact is written
# down once here and read wherever it is needed.
#
# The rule exists because it was broken. Two hand-maintained lists in this
# package disagreed about which families carry a residual scale --
# `.brms_family_has_sigma()` in R/priors_to_brms.R listed gamma and beta,
# the heteroscedastic-residual gate in R/emit_brms.R did not -- and the
# wrong one governed the prior emit. Every prior route on gamma and beta,
# the plain default included, therefore sent brms a `sigma` prior for a
# parameter the model does not have, and brms refused the fit: two of the
# six advertised families were unreachable on the Stan backend.
#
# Each table names the engine call that grounds it, and a test in
# tests/testthat/test-family-traits.R re-reads that call from the
# installed engine and diffs it against the table. The engine is the
# oracle; the table is a cache of the engine's answer, not a second
# opinion about it.

# ---------------------------------------------------------------- #
# brms: which families carry a residual sigma                       #
# ---------------------------------------------------------------- #

# The flexyBayes family spellings whose brms counterpart parameterises a
# residual standard deviation `sigma`.
#
# Ground truth is `brms::brmsfamily(<brms family>)$dpars`, brms's own
# declaration of a family's distributional parameters, read through the
# spelling map in `.fb_family_to_brms()` (R/emit_brms.R). At brms 2.23.0:
# gaussian and lognormal carry `mu, sigma`; Gamma and negbinomial carry
# `mu, shape`; Beta carries `mu, phi`; bernoulli and poisson carry `mu`
# alone. Dispersion is not sigma -- a prior written for `sigma` cannot be
# retargeted onto `shape` or `phi`, because they are not on the
# standard-deviation scale and are not the same parameter.
#
# A family outside this vector still takes priors on its fixed effects
# and its variance components; only the residual-scale row is dropped,
# exactly as `.build_inla_control_family()` drops it on the other engine.
.fb_brms_families_with_sigma <- function() {
  c("gaussian", "lognormal")
}

# Predicate form, for the emit paths. A NULL family is the gaussian
# default, which does carry a sigma.
.fb_family_has_brms_sigma <- function(fam) {
  if (is.null(fam)) {
    return(TRUE)
  }
  fam <- if (inherits(fam, "family") || inherits(fam, "brmsfamily")) {
    fam$family
  } else {
    fam
  }
  tolower(as.character(fam)) %in% .fb_brms_families_with_sigma()
}


# ---------------------------------------------------------------- #
# INLA: which likelihood declares which residual hyperparameter     #
# ---------------------------------------------------------------- #

# The keyword INLA's `control.family$hyper` expects for the residual-
# scale hyperparameter of each likelihood flexyBayes emits, keyed by the
# INLA family name (the resolved spelling, not the flexyBayes one).
#
# Ground truth is `INLA::inla.models()$likelihood[[<family>]]$hyper`,
# whose entries carry the keyword as `short.name`. At INLA 25.10.19:
# gaussian declares `prec` (and `precoffset`), lognormal, logistic, t and
# gamma declare `prec`, and beta declares **`phi`** -- which is why a
# residual prior on a beta likelihood used to fail with INLA's own
# "Unknown keyword in `hyper` ` prec `" rather than fitting. The engine's
# error message quoted the right alternatives verbatim; the table did
# not.
#
# A likelihood absent from this vector carries no residual-scale
# hyperparameter this package sets: poisson and binomial declare none at
# all, and nbinomial's `size` and betabinomial's `rho` are overdispersion
# parameters rather than scales, so an SD-scale prior does not describe
# them. The sigma specification is dropped for those, as it is on the
# brms side for the families with no residual sigma.
.fb_inla_residual_hyper <- function() {
  c(
    gaussian = "prec",
    lognormal = "prec",
    logistic = "prec",
    t = "prec",
    gamma = "prec",
    beta = "phi"
  )
}

# The keyword for one INLA family name, or NULL when the likelihood
# carries no residual-scale hyperparameter.
.fb_inla_hyper_keyword <- function(family) {
  keywords <- .fb_inla_residual_hyper()
  key <- tolower(as.character(family %||% "")[[1L]])
  if (!key %in% names(keywords)) {
    return(NULL)
  }
  unname(keywords[[key]])
}
