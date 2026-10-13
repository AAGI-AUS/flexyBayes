# Targeted muffle for the low effective-sample-size warnings of Stan and
# posterior, for structural assertions on small brms fits only.
#
# Some brms fits run at a tiny budget (n_samples = warmup = 200, chains = 1)
# because the property under test is structural: which backend a fit routes
# to, or whether a level name round-trips. rstan's throw_sampler_warnings()
# then raises the bulk and tail ESS warnings, and posterior's .ess(), called
# by summary(), adds "The ESS has been capped to avoid unstable estimates."
#
# Only those warnings are muffled. Divergent transitions, max-treedepth, low
# BFMI and R-hat warnings still propagate, which a blanket
# suppressWarnings() would hide. Where the property under test is
# quantitative, raise the budget instead: muffling would hide a small-sample
# bias in the estimate.

.ess_warning_pattern <- paste0(
  "mc-stan\\.org/misc/warnings\\.html#(bulk|tail)-ess",
  "|ESS has been capped to avoid unstable estimates"
)

.muffle_ess_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl(.ess_warning_pattern, conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
