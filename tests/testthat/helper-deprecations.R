# The suite deliberately exercises compatibility wrappers throughout. Their
# warning behaviour is asserted explicitly in test-deprecation-conditions.R;
# elsewhere, mute deprecations so legacy-equivalence tests remain readable.
# All other warning classes still fail normally.
# The dedicated warning-contract files opt back in locally.
options(flexyBayes.deprecation_warnings = FALSE)
