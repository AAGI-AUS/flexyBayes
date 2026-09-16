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
# brms: which families take a missing response through mi()        #
# ---------------------------------------------------------------- #

# The flexyBayes family spellings for which `na_action = "augment"` on brms
# samples a missing response through brms's `mi()` addition term.
#
# Ground truth for what brms accepts is brms's own family table: at brms
# 2.23.0 `mi()` is accepted for gaussian, lognormal, student and gamma, and
# refused for poisson and bernoulli. A family is admitted here only once a
# live fit shows the missing rows carried and the augmented posterior
# matching the complete-case one (tests/testthat/test-brms-missing-response.R).
# Student-t is accepted by brms but not admitted, because emit_brms() has no
# translation for it. That this vector equals the residual-sigma roster above
# is a coincidence of the two validations, not a shared rule.
.fb_brms_families_augmentable <- function() {
  c("gaussian", "lognormal")
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

# ---------------------------------------------------------------- #
# Response support: which families require which response domain   #
# ---------------------------------------------------------------- #

# The response domain each flexyBayes family's likelihood is defined
# over, keyed by the flexyBayes family spelling. A family absent from
# this table is unrestricted (gaussian, binomial, binary, poisson,
# negative_binomial, negbinom): the response gate has no domain check for
# it here, and the existing binomial-specific guard
# (.refuse_non_binary_binomial_response(), R/dispatch.R) covers the
# 0/1 constraint that family carries.
#
# Ground truth is each engine's own likelihood, read off a live failure
# rather than assumed from a family's textbook support: at brms 2.23.0,
# brm() stops with "Family '<x>' requires response greater than 0" for
# gamma and lognormal, "... greater than or equal to 0" for hurdle_gamma,
# and both "... greater than 0" and "... smaller than 1" for beta; at
# INLA 25.10.19 the inla-program dies on a zero or negative gamma
# observation ("(m)Gamma data[i] (y) = <v> ... is void") and on a beta
# observation at the boundary ("Beta data: N observations are either 0
# or 1"). Before this table, admitting a family onto the roster
# (.fb_family_links(), R/utils.R) was necessary and sufficient to reach
# an engine; nothing checked the response against what the family can
# represent, so lognormal (added b979c5d) and every other bounded-support
# family here reached brms::brm() or the INLA subprocess and died
# unrefused (U-FAM-SUPPORT, 2026-09-16).
#
#   "positive"      -- strictly greater than 0 (gamma, lognormal)
#   "nonnegative"   -- 0 or greater, zero admitted (hurdle_gamma: the
#                      hurdle mass point is exactly zero, the continuous
#                      part is gamma)
#   "unit_interval" -- strictly between 0 and 1, both bounds excluded
#                      (beta)
.fb_family_response_domain <- function() {
  list(
    gamma = "positive",
    hurdle_gamma = "nonnegative",
    lognormal = "positive",
    beta = "unit_interval"
  )
}

# Row indices of `y` that violate `domain`. NA entries are never flagged:
# a comparison against NA is NA, and which() reads only TRUE, so a
# response left NA for na_action to carry as a latent quantity is passed
# over rather than read as an out-of-domain value.
.fb_response_domain_violation_rows <- function(domain, y) {
  switch(
    domain,
    positive = which(y <= 0),
    nonnegative = which(y < 0),
    unit_interval = which(y <= 0 | y >= 1),
    integer(0)
  )
}

# The domain, worded for the refusal message.
.fb_response_domain_describe <- function(domain) {
  switch(
    domain,
    positive = "strictly greater than 0",
    nonnegative = "0 or greater",
    unit_interval = "strictly between 0 and 1",
    domain
  )
}

# The nearest supported alternative for a family whose response domain
# check just failed. gamma's own hurdle variant is the closer fit when
# every offending value is an exact zero and none is negative -- the
# hurdle mass point admits the zero, and the gamma part is unaffected;
# every other case falls back to gaussian, the one family in the roster
# with no domain restriction at all.
.fb_response_domain_alternative <- function(fam, offending) {
  if (
    identical(fam, "gamma") &&
      length(offending) > 0L &&
      all(offending == 0)
  ) {
    return(paste0(
      "family = \"hurdle_gamma\" admits an exact zero as a hurdle mass ",
      "point plus a positive-gamma continuous part for the rest; use it ",
      "if the zeros are genuine observations."
    ))
  }
  paste0(
    "family = \"gaussian\" places no restriction on the response and is ",
    "the nearest supported alternative."
  )
}
