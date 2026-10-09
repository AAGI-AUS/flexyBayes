# testthat setup -- runs once before the suite.
#
# Silence the fit-time convergence warning across the suite. Almost every
# fit here uses a deliberately tiny MCMC budget for speed, not for
# convergence, so the warning would be pervasive noise and would trip the
# expect_silent() / expect_no_warning() probes that target other
# behaviour. The warning itself is exercised on purpose in
# test-convergence-warning.R, which re-enables it locally.
options(flexyBayes.silence_convergence_warning = TRUE)

# Each INLA fit in the suite runs once. A fit run three times starts two
# further R processes, which do not see a test's mocked bindings; the
# best-of-k tests set the option themselves (test-inla-best-of.R).
options(flexyBayes.inla_best_of = 1L)
