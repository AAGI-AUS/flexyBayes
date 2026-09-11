# scenarios.R -- the paired scenario definitions of SPEC.md.
#
# Each scenario supplies the data, the two calls, and the invariants to
# compare. The calls sit side by side deliberately: a reader checking that the
# pairing is fair should not have to look in two files to do it.

SEED <- 20260911L

scenario_registry <- function() {
  list(

    # ---- 01 fixed effects + one iid random intercept -------------------- #
    "01_iid_intercept" = list(
      label = "Fixed genotype effects with an iid replicate effect",
      registry_row = "Simple iid random intercept",
      factor_cols = c("gen", "rep"),
      data = function() {
        e <- new.env()
        utils::data("john.alpha", package = "agridat", envir = e)
        d <- e$john.alpha
        d$gen <- factor(d$gen); d$rep <- factor(d$rep); d$block <- factor(d$block)
        d
      },
      asreml = function(d) {
        asreml::asreml(fixed = yield ~ gen, random = ~ rep,
                       data = d, trace = FALSE)
      },
      flexybayes = function(d) {
        flexyBayes::flexybayes(fixed = yield ~ gen, random = ~ rep,
                               data = d, backend = "inla", verbose = FALSE)
      },
      ratios = list(list(asreml = "rep", flexybayes = "sd_rep")),
      blup = list(list(asreml = "rep", flexybayes = "rep", label = "rep"))
    ),

    # ---- 02 nested / multi-stratum -------------------------------------- #
    "02_nested_stratum" = list(
      label = "Alpha design: replicate and block-within-replicate strata",
      registry_row = "Nested and interacted iid terms",
      factor_cols = c("gen", "rep", "block"),
      data = function() {
        e <- new.env()
        utils::data("john.alpha", package = "agridat", envir = e)
        d <- e$john.alpha
        d$gen <- factor(d$gen); d$rep <- factor(d$rep); d$block <- factor(d$block)
        d
      },
      asreml = function(d) {
        asreml::asreml(fixed = yield ~ gen, random = ~ rep + rep:block,
                       data = d, trace = FALSE)
      },
      flexybayes = function(d) {
        flexyBayes::flexybayes(fixed = yield ~ gen, random = ~ rep + rep:block,
                               data = d, backend = "inla", verbose = FALSE)
      },
      ratios = list(list(asreml = "rep", flexybayes = "sd_rep")),
      blup = list(list(asreml = "rep", flexybayes = "rep", label = "rep"))
    ),

    # ---- 03 diagonal heterogeneous G by E ------------------------------- #
    "03_diag_gxe" = list(
      label = "One genotype variance per environment, no covariance",
      registry_row = "Diagonal heterogeneous variance over a factor",
      factor_cols = c("gen", "env"),
      data = function() {
        e <- new.env()
        utils::data("besag.met", package = "agridat", envir = e)
        d <- e$besag.met
        d$env <- factor(d$county); d$gen <- factor(d$gen); d$rep <- factor(d$rep)
        d <- d[!is.na(d$yield), ]
        d <- d[order(d$env, d$gen), ]
        rownames(d) <- NULL
        d
      },
      asreml = function(d) {
        asreml::asreml(fixed = yield ~ env, random = ~ diag(env):gen,
                       data = d, trace = FALSE)
      },
      flexybayes = function(d) {
        flexyBayes::flexybayes(fixed = yield ~ env, random = ~ diag(env):gen,
                               data = d, backend = "inla", verbose = FALSE)
      },
      ratios = list(),
      blup = list()
    ),

    # ---- 04 sectioned residual ------------------------------------------ #
    "04_dsum_residual" = list(
      label = "One residual variance per environment",
      registry_row = "Sectioned residual, one variance per section",
      factor_cols = c("gen", "env"),
      data = function() {
        e <- new.env()
        utils::data("yan.winterwheat", package = "agridat", envir = e)
        d <- e$yan.winterwheat
        d$gen <- factor(d$gen); d$env <- factor(d$env)
        d <- d[order(d$env, d$gen), ]
        rownames(d) <- NULL
        d
      },
      asreml = function(d) {
        asreml::asreml(fixed = yield ~ env, random = ~ gen,
                       residual = ~ dsum(~ units | env),
                       data = d, trace = FALSE)
      },
      flexybayes = function(d) {
        flexyBayes::flexybayes(fixed = yield ~ env, random = ~ gen,
                               residual = ~ dsum(~ units | env),
                               data = d, backend = "inla", verbose = FALSE)
      },
      ratios = list(list(asreml = "gen", flexybayes = "sd_gen")),
      blup = list(list(asreml = "gen", flexybayes = "gen", label = "gen"))
    ),

    # ---- 05 combined MET ------------------------------------------------- #
    "05_combined_met" = list(
      label = "Heterogeneous G by E together with a sectioned residual",
      registry_row = "Sectioned residual, one variance per section",
      factor_cols = c("gen", "env"),
      data = function() {
        e <- new.env()
        utils::data("besag.met", package = "agridat", envir = e)
        d <- e$besag.met
        d$env <- factor(d$county); d$gen <- factor(d$gen); d$rep <- factor(d$rep)
        d <- d[!is.na(d$yield), ]
        d <- d[order(d$env, d$gen), ]
        rownames(d) <- NULL
        d
      },
      asreml = function(d) {
        asreml::asreml(fixed = yield ~ env, random = ~ diag(env):gen,
                       residual = ~ dsum(~ units | env),
                       data = d, trace = FALSE)
      },
      flexybayes = function(d) {
        flexyBayes::flexybayes(fixed = yield ~ env, random = ~ diag(env):gen,
                               residual = ~ dsum(~ units | env),
                               data = d, backend = "inla", verbose = FALSE)
      },
      ratios = list(),
      blup = list()
    )
  )
}
