# A prior for tests that fit agridat::john.alpha as `rep + rep:block` on INLA.
#
# Under the default uniform bound INLA's reported interval for the three-level
# `rep` runs past the bound (0.17 to 3.65 against 3.05), and the support audit
# refuses the fit. Tests that fit this design for another purpose state a
# half-normal prior on each standard deviation, with scale the default bound
# itself, five standard deviations of the response. It has no bound to cross,
# it moves with the response as the default does, so a rescaling check still
# holds, and it is weak enough that `rep` stays as weakly identified as under
# the default: a scale of one standard deviation identified it (posterior CV
# 0.60), which the theta-oracle and covariance snapshots record as a change.
.fb_test_john_alpha_prior <- function(d, response = "yield") {
  s <- 5 * stats::sd(d[[response]])
  fb_prior(
    sigma ~ half_normal(scale = s),
    sd(group = "rep") ~ half_normal(scale = s),
    sd(group = "rep:block") ~ half_normal(scale = s)
  )
}
