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

    # ---- 06 / 07 unstructured covariance, dimensions 2 and 3 ------------- #
    #
    # us() is where the two systems' parameterisations differ most: ASReml
    # estimates the covariance directly, INLA carries an iidKd field under a
    # joint Wishart. The comparison is therefore on the reconstructed
    # covariance and the correlations derived from it, never on raw
    # parameters, which is what SPEC.md's covariance tolerance is for.
    "06_us_dim2" = list(
      label = "Unstructured genotype covariance across two environments",
      registry_row = "Unstructured covariance over a factor, on INLA",
      factor_cols = c("gen", "env"),
      data = function() {
        e <- new.env()
        utils::data("besag.met", package = "agridat", envir = e)
        d <- e$besag.met
        d$env <- factor(d$county); d$gen <- factor(d$gen); d$rep <- factor(d$rep)
        d <- d[!is.na(d$yield), ]
        d <- d[d$env %in% levels(d$env)[seq_len(2)], ]
        d$env <- droplevels(d$env); d$gen <- droplevels(d$gen)
        d <- d[order(d$env, d$gen), ]
        rownames(d) <- NULL
        d
      },
      asreml = function(d) {
        asreml::asreml(fixed = yield ~ env, random = ~ us(env):gen,
                       data = d, trace = FALSE)
      },
      flexybayes = function(d) {
        flexyBayes::flexybayes(fixed = yield ~ env, random = ~ us(env):gen,
                               data = d, backend = "inla", verbose = FALSE)
      },
      ratios = list(),
      blup = list(),
      # The quantity us() claims to represent, and the one the scenario
      # existed to test. Without this the run banked "1 pass" -- the
      # intercept -- while the genetic correlation was 0.14 against 0.93,
      # which is the instrument answering a question nobody asked. ASReml
      # reports the covariance directly, so the correlation is reconstructed
      # from its three parameters rather than read off a row.
      correlations = list(
        list(
          label = "gen, C1 vs C2",
          asreml = c(v1 = "env:gen!env_C1:C1", v2 = "env:gen!env_C2:C2",
                     cov = "env:gen!env_C2:C1"),
          flexybayes = "cor_gen__envC1__envC2"
        )
      )
    ),

    "07_us_dim3" = list(
      label = "Unstructured genotype covariance across three environments",
      registry_row = "Unstructured covariance over a factor, on INLA",
      factor_cols = c("gen", "env"),
      data = function() {
        e <- new.env()
        utils::data("besag.met", package = "agridat", envir = e)
        d <- e$besag.met
        d$env <- factor(d$county); d$gen <- factor(d$gen); d$rep <- factor(d$rep)
        d <- d[!is.na(d$yield), ]
        d <- d[d$env %in% levels(d$env)[seq_len(3)], ]
        d$env <- droplevels(d$env); d$gen <- droplevels(d$gen)
        d <- d[order(d$env, d$gen), ]
        rownames(d) <- NULL
        d
      },
      asreml = function(d) {
        asreml::asreml(fixed = yield ~ env, random = ~ us(env):gen,
                       data = d, trace = FALSE)
      },
      flexybayes = function(d) {
        flexyBayes::flexybayes(fixed = yield ~ env, random = ~ us(env):gen,
                               data = d, backend = "inla", verbose = FALSE)
      },
      ratios = list(),
      blup = list(),
      correlations = list(
        list(
          label = "gen, C1 vs C2",
          asreml = c(v1 = "env:gen!env_C1:C1", v2 = "env:gen!env_C2:C2",
                     cov = "env:gen!env_C2:C1"),
          flexybayes = "cor_gen__envC1__envC2"
        ),
        list(
          label = "gen, C1 vs C3",
          asreml = c(v1 = "env:gen!env_C1:C1", v2 = "env:gen!env_C3:C3",
                     cov = "env:gen!env_C3:C1"),
          flexybayes = "cor_gen__envC1__envC3"
        ),
        list(
          label = "gen, C2 vs C3",
          asreml = c(v1 = "env:gen!env_C2:C2", v2 = "env:gen!env_C3:C3",
                     cov = "env:gen!env_C3:C2"),
          flexybayes = "cor_gen__envC2__envC3"
        )
      )
    ),

    # ---- 15 Gaussian relative weights ------------------------------------ #
    #
    # Both systems claim variance proportional to sigma^2 / w. The check is
    # that they agree on the fixed effects, which is where a mis-scaled
    # weight shows up first.
    #
    # Numbered 15 to match SPEC.md's matrix, not in run order. The number is
    # the registry row's identity and the fixtures are named by it, so a
    # scenario that drifts from its number in the specification is a claim
    # surface out of step with its own record.
    "15_weights" = list(
      label = "Gaussian relative weights on a replicated trial",
      registry_row = "Gaussian relative weights",
      factor_cols = c("gen", "rep"),
      data = function() {
        e <- new.env()
        utils::data("john.alpha", package = "agridat", envir = e)
        d <- e$john.alpha
        d$gen <- factor(d$gen); d$rep <- factor(d$rep); d$block <- factor(d$block)
        set.seed(20260911L)
        # Weights varying by a factor of four, so a dropped or inverted
        # weight cannot pass unnoticed.
        d$wt <- rep(c(1, 2, 4, 2), length.out = nrow(d))
        d
      },
      asreml = function(d) {
        asreml::asreml(fixed = yield ~ gen, random = ~ rep, weights = wt,
                       data = d, trace = FALSE)
      },
      flexybayes = function(d) {
        flexyBayes::flexybayes(fixed = yield ~ gen, random = ~ rep,
                               weights = d$wt, data = d,
                               backend = "inla", verbose = FALSE)
      },
      ratios = list(list(asreml = "rep", flexybayes = "sd_rep")),
      blup = list(list(asreml = "rep", flexybayes = "rep", label = "rep"))
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
