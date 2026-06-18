# testthat setup -- runs once before the suite.
#
# Silence the fit-time convergence warning across the suite. Almost every
# fit here uses a deliberately tiny MCMC budget for speed, not for
# convergence, so the warning would be pervasive noise and would trip the
# expect_silent() / expect_no_warning() probes that target other
# behaviour. The warning itself is exercised on purpose in
# test-convergence-warning.R, which re-enables it locally.
options(flexyBayes.silence_convergence_warning = TRUE)
