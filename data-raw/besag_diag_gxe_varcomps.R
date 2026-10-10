# data-raw/besag_diag_gxe_varcomps.R
#
# The variance-component table of an INLA fit of agridat::besag.met with
# random = ~ diag(env):gen (shape 03_diag_gxe of the ASReml parity suite),
# under the package defaults: best of 3, one thread, default priors.
# tests/testthat/test-weak-variance-flag.R reads it to check the
# weak-identification flag on counties C2 and C5, which ASReml 4.2 estimates
# at 3.28 (std.error 23.6) and 3.77 (23.7) (inst/extdata/asreml_parity.rds).
#
# Run from the package root:
#
#   Rscript --no-init-file data-raw/besag_diag_gxe_varcomps.R
pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-asreml-parity-data.R")
options(
  flexyBayes.inla_threads = "1:1",
  flexyBayes.silence_default_prior_note = TRUE
)
message("best of ", getOption("flexyBayes.inla_best_of", 3L))
d <- .ap_data("03_diag_gxe")
fit <- withCallingHandlers(
  flexybayes(yield ~ env, random = ~ diag(env):gen, data = d,
             backend = "inla", verbose = FALSE),
  warning = function(w) {
    message("warning: ", conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
vc <- fit$extras$variance_comps
med <- attr(vc, "posterior_median")[vc$component]
print(data.frame(
  component = vc$component,
  sd_q2.5 = vc$q2.5,
  sd_median = unname(med),
  sd_q97.5 = vc$q97.5,
  var_median = unname(med)^2,
  ratio = (vc$q97.5 - vc$q2.5) / unname(med),
  prior_share = unname(attr(vc, "prior_share")[vc$component])
), digits = 4)
print(summary(fit)$varcomp)
print(fit$extras$reproducibility$best_of)
saveRDS(vc, "tests/testthat/fixtures/besag_diag_gxe_varcomps.rds", version = 2)
