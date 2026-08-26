# helper-ess-warnings.R -- targeted muffle for Stan/posterior's
# low-effective-sample-size warnings, for structural brms fit
# assertions only.
#
# WP-C review F1: three brms fits across test-plan-tells-truth.R and
# test-space-level-factor.R run at a deliberately tiny sampler budget
# (n_samples = warmup = 200, chains = 1) because the property under
# test is structural (which backend a fit routes to; whether a random-
# effect level round-trips unmangled) and does not depend on effective
# sample size. rstan's own `throw_sampler_warnings()`
# (`asNamespace("rstan")$throw_sampler_warnings`) raises these two
# warnings, verbatim, whenever a chain's bulk or tail ESS falls below
# 100 * n_chains:
#
#   "Bulk Effective Samples Size (ESS) is too low, indicating
#   posterior means and medians may be unreliable. Running the chains
#   for more iterations may help. See
#   https://mc-stan.org/misc/warnings.html#bulk-ess"
#
#   "Tail Effective Samples Size (ESS) is too low, indicating
#   posterior variances and tail quantiles may be unreliable. Running
#   the chains for more iterations may help. See
#   https://mc-stan.org/misc/warnings.html#tail-ess"
#
# A second, related warning fires from `posterior:::.ess()` itself --
# the shared internal estimator BOTH `posterior::ess_bulk()` and
# `ess_tail()` call, which flexyBayes's own post-fit convergence
# summary (the "Min ESS (bulk):" / "Min ESS (tail):" lines printed by
# `summary()` / `print()`) invokes on every brms fit, independent of
# rstan's own sampler-warning pass:
#
#   "The ESS has been capped to avoid unstable estimates."
#
# This fires when the autocorrelation-time estimate itself is too
# unstable at very few draws (`tau_hat < 1/log10(ess)`, verified
# against the installed `posterior` package's own `.ess()` source) --
# the same root cause (a deliberately tiny structural-test sampler
# budget), a different call site (flexyBayes's own summary
# construction, not rstan's post-sampling check), so it does not carry
# rstan's mc-stan.org/misc/warnings.html URL. Folded into the same
# pattern below rather than a second helper.
#
# `.muffle_ess_warnings()` mirrors `base::suppressWarnings()`'s own
# withCallingHandlers() + invokeRestart("muffleWarning") pattern
# exactly, but filters on the ESS-warning pattern only -- divergent
# transitions, max-treedepth, low BFMI and R-hat warnings (rstan's
# other four sampler warnings, each with its own
# mc-stan.org/misc/warnings.html anchor) still propagate normally, so
# a genuinely broken fit is not silenced by using this helper. Never
# use a blanket suppressWarnings() around a fit call: that would also
# hide those four.
#
# Where the property under test IS quantitative (a numeric agreement
# or a metamorphic ratio), the correct fix is to raise the sampler
# budget until the warning stops firing on its own, not to muffle it
# -- muffling would hide a genuine small-sample bias in the estimate
# the test is asserting on.

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
