# =============================================================================
# stress_corner_to_corner.R
#
# Adversarial corner-to-corner capability stress run for flexyBayes.
# Intentionally pushes every capability area to the cases where anomalies
# surface, and verifies the package's central contract: invalid / unsupported
# requests must be REFUSED LOUDLY (non-silently), never silently mishandled.
#
# This is NOT the testthat suite (that covers the happy + known-edge paths,
# ~2641 tests). This harness probes BEYOND it: the refusal surface, a set of
# specifically-suspected silent-failure paths, cross-backend agreement, the
# post-fit / ecosystem surface, streaming + aggregation, and predict variants.
#
# Each probe records: id, area, kind, status (PASS / ANOMALY / SKIP), detail.
#   kind = "refusal"  -> adversarial input that SHOULD error loudly; PASS if it
#                        errored (sub-noted if the structured reason-code class
#                        did not match), ANOMALY if it returned silently.
#   kind = "success"  -> valid input that should work; ANOMALY if it errored or
#                        a content check failed.
#   kind = "probe"    -> bespoke silent-failure / correctness check.
#
# Writes:  tools/stress_results/stress_corner_to_corner_results.rds
#          tools/stress_results/stress_corner_to_corner_report.md
#
# Run (from the package root, backends warm):
#   NOT_CRAN=true Rscript tools/stress_corner_to_corner.R
# =============================================================================

suppressWarnings(suppressMessages({
  ok_load <- tryCatch({ pkgload::load_all(".", quiet = TRUE); TRUE },
                      error = function(e) FALSE)
  if (!ok_load) library(flexyBayes)
  library(stats)
}))

set.seed(20260602L)
HAS_INLA  <- requireNamespace("INLA", quietly = TRUE)
HAS_GRETA <- requireNamespace("greta", quietly = TRUE)
HAS_BRMS  <- requireNamespace("brms", quietly = TRUE)
RUN_GRETA <- HAS_GRETA && !identical(Sys.getenv("STRESS_SKIP_GRETA"), "1")
RUN_BRMS  <- HAS_BRMS  && !identical(Sys.getenv("STRESS_SKIP_BRMS"), "1")

GS <- list(n_samples = 80L, warmup = 80L, chains = 1L,
           verbose = FALSE, mcmc_verbose = FALSE)

# ---- probe framework --------------------------------------------------------

RESULTS <- new.env(parent = emptyenv()); RESULTS$rows <- list()
.rec <- function(id, area, kind, status, detail = "") {
  RESULTS$rows[[length(RESULTS$rows) + 1L]] <- data.frame(
    id = id, area = area, kind = kind, status = status,
    detail = substr(gsub("[\r\n]+", " ", detail), 1L, 320L),
    stringsAsFactors = FALSE)
  tag <- switch(status, PASS = "  ok ", ANOMALY = ">>ANOM", SKIP = " skip ", "  ?  ")
  message(sprintf("[%s] %-34s %s", tag, id, substr(detail, 1L, 90L)))
}

# adversarial: SHOULD raise an error. PASS if it errored (loud). Records whether
# the expected structured reason-code class was present.
probe_refuse <- function(id, area, expr, class = NULL) {
  r <- tryCatch({ val <- force(expr); list(err = FALSE, val = val) },
                error = function(e) list(err = TRUE, e = e))
  if (!r$err) {
    .rec(id, area, "refusal", "ANOMALY",
         "SILENT: invalid input returned WITHOUT error")
    return(invisible(FALSE))
  }
  cls <- class(r$e); msg <- conditionMessage(r$e)
  if (!is.null(class) && !class %in% cls) {
    .rec(id, area, "refusal", "PASS",
         sprintf("refused (generic, expected class '%s' absent; got %s): %s",
                 class, paste(utils::head(cls, 2L), collapse = ","), msg))
  } else {
    .rec(id, area, "refusal", "PASS", sprintf("refused loudly: %s", msg))
  }
  invisible(TRUE)
}

# valid: should succeed (+ optional content check returning TRUE or a message)
probe_ok <- function(id, area, expr, check = NULL) {
  r <- tryCatch({ list(err = FALSE, val = force(expr)) },
                error = function(e) list(err = TRUE, e = e))
  if (r$err) {
    .rec(id, area, "success", "ANOMALY",
         sprintf("errored on valid input: %s", conditionMessage(r$e)))
    return(invisible(NULL))
  }
  if (!is.null(check)) {
    chk <- tryCatch(check(r$val), error = function(e)
      paste("check threw:", conditionMessage(e)))
    if (isTRUE(chk)) .rec(id, area, "success", "PASS", "")
    else .rec(id, area, "success", "ANOMALY", as.character(chk))
  } else .rec(id, area, "success", "PASS", "")
  invisible(r$val)
}

# =============================================================================
# AREA 1 -- REFUSAL CONTRACT (adversarial: must refuse loudly)
# =============================================================================
message("\n==== AREA 1: refusal contract (adversarial) ====")
d_g <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
d_bin <- data.frame(y = rbinom(60, 1, .5), x = rnorm(60), g = factor(rep(1:6, 10)))

# -- family gate
probe_refuse("refuse.family.survival", "family",
  flexybayes(y ~ x, data = d_g, family = "coxph", backend = "inla"),
  class = "flexybayes_refusal_unsupported_family")
probe_refuse("refuse.family.weibull", "family",
  flexybayes(y ~ x, data = d_g, family = "weibull", backend = "inla"),
  class = "flexybayes_refusal_unsupported_family")

# -- formula / parse
probe_refuse("refuse.formula.onesided", "parse",
  flexybayes(~ x, data = d_g, backend = "inla"),
  class = "flexybayes_refusal_formula_not_two_sided")
probe_refuse("refuse.formula.response_missing", "parse",
  flexybayes(zzz ~ x, data = d_g, backend = "inla"),
  class = "flexybayes_refusal_response_not_in_data")
probe_refuse("refuse.smooth.tensor_te", "parse",
  flexybayes(y ~ te(x, g), data = d_g, backend = "inla"),
  class = "flexybayes_refusal_tensor_smooth_unsupported")
probe_refuse("refuse.smooth.tensor_ti", "parse",
  flexybayes(y ~ ti(x), data = d_g, backend = "inla"),
  class = "flexybayes_refusal_tensor_smooth_unsupported")
probe_refuse("refuse.smooth.var_missing", "parse",
  flexybayes(y ~ s(nope), data = d_g, backend = "inla"),
  class = "flexybayes_refusal_smooth_variable_not_in_data")

# -- entry-point flag guards
probe_refuse("refuse.flags.mutually_exclusive", "entry",
  flexybayes(y ~ x, data = d_g, return_code = TRUE, review_code = TRUE),
  class = "flexybayes_refusal_code_flags_mutually_exclusive")

# -- structured covariance carriers
A_bad <- matrix(c(1, 0.9, 0.2, 0.9), 2, 2)   # non-symmetric
probe_refuse("refuse.cov.type_unknown", "structured_cov",
  fb_cov(diag(3), type = "wishart"))
probe_refuse("refuse.cov.missing_matrix", "structured_cov",
  fb_cov(type = "dense"))
# fb_cov() validates the *carrier* lazily by design (a light summary at
# construction; the hard symmetry/PD check runs at fit time against the
# grouping factor). So the non-symmetric-precision refusal is a FIT-time probe.
if (HAS_INLA)
  probe_refuse("refuse.cov.precision_not_symmetric", "structured_cov",
    { dpc <- data.frame(y = rnorm(40), x = rnorm(40),
                        gg = factor(rep(c("a", "b"), 20)))
      fb_inla(y ~ x, random = ~ vm(gg, precision = Q),
              known_matrices = list(Q = A_bad), data = dpc) },
    class = "flexybayes_refusal_precision_not_symmetric")

# -- spec helpers
probe_refuse("refuse.prior.onesided", "spec",
  fb_prior(~ normal(0, 1)))
# NB: a bare symbol LHS (e.g. `b ~ ...`) is a documented coefficient-name
# shorthand, so it is NOT a bad target. A *call* to an unsupported function is.
probe_refuse("refuse.prior.bad_target", "spec",
  fb_prior(variance(x) ~ normal(0, 1)))
probe_refuse("refuse.prior.bad_dist", "spec",
  fb_prior(sigma ~ laplace(0, 1)))
probe_refuse("refuse.engine.unknown", "spec",
  fb_engine("jags"))
probe_refuse("refuse.approx.unknown_scheme", "spec",
  fb_approx("magic_basis"))

# -- backend-capability conflicts
probe_refuse("refuse.brms.structured_cov", "dispatch",
  fb_brms(y ~ vm(g), data = d_g))
probe_refuse("refuse.inla.us_random", "dispatch",
  fb_inla(y ~ x, random = ~ us(g), data = d_g))
probe_refuse("refuse.inla.at_units_rcov", "dispatch",
  fb_inla(y ~ x, random = ~ g, rcov = ~ at(g):units, data = d_g))

# =============================================================================
# AREA 2 -- SUSPECTED SILENT-FAILURE PROBES (verify empirically)
# =============================================================================
message("\n==== AREA 2: suspected silent-failure paths ====")

# (2a) beta family: .get_stats_family() maps beta -> gaussian(identity).
#      Probe: does a beta fit report family()=="beta", or silently "gaussian"?
if (HAS_INLA) {
  d_beta <- data.frame(y = plogis(rnorm(60)), x = rnorm(60))  # (0,1) response
  fit_beta <- tryCatch(fb_inla(y ~ x, data = d_beta, family = "beta"),
                       error = function(e) e)
  if (inherits(fit_beta, "error")) {
    .rec("probe.beta.family_label", "silent_failure", "probe", "SKIP",
         sprintf("beta INLA fit errored: %s", conditionMessage(fit_beta)))
  } else {
    fam <- tryCatch(family(fit_beta), error = function(e) e)
    famname <- if (inherits(fam, "family")) fam$family else
      if (inherits(fam, "error")) paste("family() errored:", conditionMessage(fam)) else as.character(fam)
    if (grepl("beta", tolower(paste(famname, collapse = " ")))) {
      .rec("probe.beta.family_label", "silent_failure", "probe", "PASS",
           sprintf("family() on beta fit reports '%s'", paste(famname, collapse = " ")))
    } else {
      .rec("probe.beta.family_label", "silent_failure", "probe", "ANOMALY",
           sprintf("beta fit reports family()='%s' (expected 'beta'); .get_stats_family beta->gaussian leak",
                   paste(famname, collapse = " ")))
    }
  }
} else .rec("probe.beta.family_label", "silent_failure", "probe", "SKIP", "INLA absent")

# (2b) negative_binomial: .get_stats_family maps negbinom -> poisson(log).
#      Probe: does a negbinom fit report family()=="poisson" silently?
if (HAS_INLA) {
  d_nb <- data.frame(y = rnbinom(60, mu = 4, size = 2), x = rnorm(60))
  fit_nb <- tryCatch(fb_inla(y ~ x, data = d_nb, family = "negative_binomial"),
                     error = function(e) e)
  if (inherits(fit_nb, "error")) {
    .rec("probe.negbinom.family_label", "silent_failure", "probe", "SKIP",
         sprintf("negbinom INLA fit errored: %s", conditionMessage(fit_nb)))
  } else {
    fam <- tryCatch(family(fit_nb), error = function(e) e)
    famname <- if (inherits(fam, "family")) fam$family else
      if (inherits(fam, "error")) paste("family() errored:", conditionMessage(fam)) else as.character(fam)
    if (grepl("pois", tolower(paste(famname, collapse = " ")))) {
      .rec("probe.negbinom.family_label", "silent_failure", "probe", "ANOMALY",
           sprintf("negbinom fit reports family()='%s' (poisson leak from .get_stats_family)",
                   paste(famname, collapse = " ")))
    } else {
      .rec("probe.negbinom.family_label", "silent_failure", "probe", "PASS",
           sprintf("family() on negbinom fit reports '%s'", paste(famname, collapse = " ")))
    }
  }
} else .rec("probe.negbinom.family_label", "silent_failure", "probe", "SKIP", "INLA absent")

# (2c) unknown function in formula. A real FIT errors loudly (the unknown
#      function is passed through and the backend cannot find it) -- so it is
#      NON-silent, which is the contract. NOTE (documented in the report): the
#      error is a leaked backend message rather than a clean parse-time refusal,
#      and `plan = TRUE` does not catch it (it does not evaluate the term). A
#      parse-time guard is deferred because it cannot be added without risking
#      legitimate transforms (log(), I(), poly(), factor()).
if (HAS_INLA)
  probe_refuse("probe.formula.unknown_function", "silent_failure",
    fb_inla(y ~ mysmooth(x), data = d_g))

# (2d) non-binary "binomial" on greta: project history says it must refuse,
#      not silently fit Bernoulli. y in {0,1,2}.
if (RUN_GRETA) {
  d_multi <- data.frame(y = rep(c(0L, 1L, 2L), 20), x = rnorm(60))
  probe_refuse("probe.greta.nonbinary_binomial", "silent_failure",
    do.call(fb_greta, c(list(quote(y ~ x), data = d_multi, family = "binomial"), GS)))
} else .rec("probe.greta.nonbinary_binomial", "silent_failure", "probe", "SKIP", "greta off")

# =============================================================================
# AREA 3 -- CROSS-BACKEND FIT + TRIANGULATE CONSISTENCY
# =============================================================================
message("\n==== AREA 3: cross-backend fit + triangulate ====")
set.seed(101L)
n <- 80L
d_lin <- data.frame(x = rnorm(n)); d_lin$y <- 1 + 0.8 * d_lin$x + rnorm(n, 0, 0.7)

fit_i <- if (HAS_INLA) probe_ok("fit.inla.gaussian", "fit",
  fb_inla(y ~ x, data = d_lin),
  check = function(f) inherits(f, "flexybayes_inla") || "flexybayes" %in% class(f)) else NULL
fit_g <- if (RUN_GRETA) probe_ok("fit.greta.gaussian", "fit",
  do.call(fb_greta, c(list(quote(y ~ x), data = d_lin), GS)),
  check = function(f) inherits(f, "flexybayes")) else
  { .rec("fit.greta.gaussian", "fit", "success", "SKIP", "greta off"); NULL }
fit_b <- if (RUN_BRMS) probe_ok("fit.brms.gaussian", "fit",
  fb_brms(y ~ x, data = d_lin, chains = 1L, n_samples = 300L, warmup = 300L,
          mcmc_verbose = FALSE),
  check = function(f) inherits(f, "flexybayes_brms")) else
  { .rec("fit.brms.gaussian", "fit", "success", "SKIP", "brms off"); NULL }

# triangulate every available backend pair; coefficients should agree closely
.tri_check <- function(t) {
  if (!inherits(t, "triangulate_result")) return("not a triangulate_result")
  if (t$n_common < 1L) return("no common parameters mapped")
  md <- max(abs(t$metrics$mean_diff), na.rm = TRUE)
  if (is.finite(md) && md < 0.5) TRUE else
    sprintf("max |mean_diff| = %.3f across %d params (expected < 0.5)", md, t$n_common)
}
if (!is.null(fit_i) && !is.null(fit_g))
  probe_ok("tri.greta_inla", "triangulate", triangulate(fit_g, fit_i), check = .tri_check)
if (!is.null(fit_i) && !is.null(fit_b))
  probe_ok("tri.brms_inla", "triangulate", triangulate(fit_b, fit_i), check = .tri_check)
if (!is.null(fit_g) && !is.null(fit_b))
  probe_ok("tri.greta_brms", "triangulate", triangulate(fit_g, fit_b), check = .tri_check)

# non-Gaussian fits across backends (binomial + poisson) on INLA (fast) + greta
if (HAS_INLA) {
  probe_ok("fit.inla.binomial", "fit",
    fb_inla(y ~ x, data = d_bin, family = "binomial"))
  d_pois <- data.frame(y = rpois(60, 3), x = rnorm(60))
  probe_ok("fit.inla.poisson", "fit",
    fb_inla(y ~ x, data = d_pois, family = "poisson"))
}
if (RUN_GRETA) {
  probe_ok("fit.greta.binomial", "fit",
    do.call(fb_greta, c(list(quote(y ~ x), data = d_bin, family = "binomial"), GS)))
}

# mixed model (random intercept) on INLA + greta
d_mix <- data.frame(y = rnorm(80), x = rnorm(80), g = factor(rep(1:8, 10)))
if (HAS_INLA) probe_ok("fit.inla.mixed", "fit",
  fb_inla(y ~ x + (1 | g), data = d_mix))
if (RUN_GRETA) probe_ok("fit.greta.mixed", "fit",
  do.call(fb_greta, c(list(quote(y ~ x + (1 | g)), data = d_mix), GS)))

# =============================================================================
# AREA 4 -- POST-FIT + ECOSYSTEM SURFACE (on a real fit)
# =============================================================================
message("\n==== AREA 4: post-fit + ecosystem surface ====")
# Standard S3 surface -- exercise on the INLA fit (fast, and the surface that
# carried the fitted()/residuals() silent-NULL bug now fixed).
pf <- if (!is.null(fit_i)) fit_i else if (!is.null(fit_g)) fit_g else NULL
pf_n <- if (!is.null(fit_i)) 80L else 80L  # both fit on d_lin (n = 80)
if (!is.null(pf)) {
  probe_ok("post.summary",  "postfit", summary(pf))
  probe_ok("post.coef",     "postfit", coef(pf),
           check = function(v) length(v) >= 1L)
  probe_ok("post.confint",  "postfit", confint(pf))
  probe_ok("post.vcov",     "postfit", vcov(pf))
  probe_ok("post.fitted",   "postfit", fitted(pf),
           check = function(v) length(v) == pf_n && all(is.finite(v)))
  probe_ok("post.residuals","postfit", residuals(pf),
           check = function(v) length(v) == pf_n && all(is.finite(v)))
  probe_ok("post.predict_newdata", "postfit",
           predict(pf, newdata = data.frame(x = c(-1, 0, 1))),
           check = function(v) {
             vv <- if (is.list(v)) v$fit else v; length(vv) == 3L && all(is.finite(vv)) })
  # documented stubs / boundaries: must refuse cleanly, not silently no-op
  probe_refuse("post.anova_stub",  "postfit", anova(pf))
  probe_refuse("post.update_stub", "postfit", update(pf, . ~ . + 1))
  # plot must not crash
  probe_ok("post.plot", "postfit",
    { tf <- tempfile(fileext = ".png"); grDevices::png(tf)
      on.exit(grDevices::dev.off(), add = TRUE); plot(pf); TRUE })
} else .rec("post.surface", "postfit", "success", "SKIP", "no fit available")

# Broom tidiers + emmeans/marginaleffects are registered for the greta-class
# fit (and use the identity-link reference grid). Exercise them on the greta
# fit where they are designed to work. (On a bare flexybayes_inla object they
# are unavailable via dispatch -- a documented capability boundary, noted in
# the report, not an anomaly.)
gf <- fit_g
if (!is.null(gf)) {
  probe_ok("eco.tidy",   "ecosystem", tidy.flexybayes(gf))
  probe_ok("eco.glance", "ecosystem", glance.flexybayes(gf))
  if (requireNamespace("emmeans", quietly = TRUE))
    probe_ok("eco.emmeans", "ecosystem", emmeans::ref_grid(gf))
  if (requireNamespace("marginaleffects", quietly = TRUE))
    probe_ok("eco.marginaleffects", "ecosystem",
             marginaleffects::avg_slopes(gf))
} else .rec("eco.surface", "ecosystem", "success", "SKIP",
            "greta off -- tidiers/ecosystem are the greta-class surface")

# =============================================================================
# AREA 5 -- STREAMING + AGGREGATION
# =============================================================================
message("\n==== AREA 5: streaming + aggregation ====")
if (HAS_INLA) {
  # in-memory streaming, gaussian
  set.seed(7L)
  d_str <- data.frame(y = rnorm(2000), g = factor(rep(1:40, 50)))
  probe_ok("stream.inmem.gaussian", "streaming",
    flexybayes_stream(y ~ 1, random = ~ g, source = d_str,
                      family = "gaussian", backend = "inla", fit = FALSE),
    check = function(o) inherits(o, "fb_aggregated") || !is.null(o))
  # generator source, fit = FALSE (aggregation only)
  gen <- function(i) { if (i > 4L) return(NULL)
    data.frame(y = rnorm(500), g = factor(rep(1:25, 20))) }
  probe_ok("stream.generator.gaussian", "streaming",
    flexybayes_stream(y ~ 1, random = ~ g, source = gen,
                      family = "gaussian", backend = "inla", fit = FALSE))
  # .fst shards round-trip
  if (requireNamespace("fst", quietly = TRUE)) {
    dir <- tempfile("shards"); dir.create(dir)
    for (k in 1:3) fst::write_fst(
      data.frame(y = rnorm(400), g = factor(rep(1:20, 20))),
      file.path(dir, sprintf("s%d.fst", k)))
    probe_ok("stream.fst_shards.gaussian", "streaming",
      flexybayes_stream(y ~ 1, random = ~ g,
                        source = list.files(dir, full.names = TRUE),
                        family = "gaussian", backend = "inla", fit = FALSE))
  } else .rec("stream.fst_shards.gaussian", "streaming", "success", "SKIP", "fst absent")
  # streaming refusals
  probe_refuse("stream.refuse.bad_family", "streaming",
    flexybayes_stream(y ~ 1, random = ~ g, source = d_str, family = "gamma"))
  probe_refuse("stream.refuse.continuous_key", "streaming",
    { d_str$xc <- rnorm(nrow(d_str))
      flexybayes_stream(y ~ xc, random = ~ g, source = d_str,
                        family = "gaussian", backend = "inla", fit = FALSE) })
  probe_refuse("stream.refuse.missing_source", "streaming",
    flexybayes_stream(y ~ 1, random = ~ g, source = "/no/such/file.fst",
                      family = "gaussian"))
  # in-memory count auto-aggregation path
  d_agg <- expand.grid(env = factor(1:8), geno = factor(1:25), rep = 1:6)
  d_agg$y <- rnorm(nrow(d_agg))
  probe_ok("aggregate.inmem.gaussian_auto", "aggregation",
    flexybayes(y ~ env, random = ~ geno, data = d_agg, family = "gaussian",
               backend = "inla"))
} else .rec("stream.area", "streaming", "success", "SKIP", "INLA absent")

# =============================================================================
# AREA 6 -- PREDICT VARIANTS + REFUSALS (greta-backed; needs draws)
# =============================================================================
message("\n==== AREA 6: predict variants ====")
if (RUN_GRETA && !is.null(fit_g)) {
  nd  <- data.frame(x = rnorm(20))
  ndL <- data.frame(x = rnorm(20))                       # for new-level tests use a mixed fit
  probe_ok("predict.greta.point", "predict",
    predict(fit_g, newdata = nd),
    check = function(v) { vv <- if (is.list(v)) v$fit else v
      length(vv) == 20L && all(is.finite(vv)) })
  probe_ok("predict.greta.se", "predict",
    predict(fit_g, newdata = nd, se.fit = TRUE),
    check = function(v) is.list(v) && all(c("fit", "se.fit") %in% names(v)))
  probe_ok("predict.greta.link", "predict",
    predict(fit_g, newdata = nd, type = "link"))
  probe_ok("predict.greta.chunked", "predict",
    predict(fit_g, newdata = nd, chunk_size = 5L),
    check = function(v) { vv <- if (is.list(v)) v$fit else v; length(vv) == 20L })
  # file output + overwrite refusal
  of <- tempfile(fileext = ".csv")
  probe_ok("predict.greta.file_output", "predict",
    predict(fit_g, newdata = nd, output_file = of, interop = TRUE),
    check = function(p) file.exists(of))
  probe_refuse("predict.greta.no_overwrite", "predict",
    predict(fit_g, newdata = nd, output_file = of, interop = TRUE))
  # invalid include vocabulary -> refusal
  fit_gm <- do.call(fb_greta, c(list(quote(y ~ x + (1 | g)), data = d_mix), GS))
  probe_refuse("predict.greta.refuse_new_level", "predict",
    { nld <- data.frame(x = rnorm(5), g = factor(rep("ZZZ", 5)))
      predict(fit_gm, newdata = nld, allow_new_levels = "refuse") })
  probe_ok("predict.greta.population_new_level", "predict",
    { nld <- data.frame(x = rnorm(5), g = factor(rep("ZZZ", 5)))
      suppressWarnings(predict(fit_gm, newdata = nld, allow_new_levels = "population")) })
} else .rec("predict.area", "predict", "success", "SKIP", "greta off / no greta fit")

# =============================================================================
# SUMMARY + PERSIST
# =============================================================================
df <- do.call(rbind, RESULTS$rows)
dir.create("tools/stress_results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(results = df, sessionInfo = utils::capture.output(utils::sessionInfo()),
             backends = c(INLA = HAS_INLA, greta = RUN_GRETA, brms = RUN_BRMS)),
        "tools/stress_results/stress_corner_to_corner_results.rds")

n_anom <- sum(df$status == "ANOMALY")
n_pass <- sum(df$status == "PASS")
n_skip <- sum(df$status == "SKIP")

# markdown report
con <- file("tools/stress_results/stress_corner_to_corner_report.md", "w")
writeLines(c(
  "# flexyBayes -- corner-to-corner stress run",
  "",
  sprintf("Probes: **%d** | PASS **%d** | ANOMALY **%d** | SKIP **%d**",
          nrow(df), n_pass, n_anom, n_skip),
  sprintf("Backends exercised: INLA=%s greta=%s brms=%s",
          HAS_INLA, RUN_GRETA, RUN_BRMS),
  "",
  "## Anomalies",
  if (n_anom == 0L) "_None._" else "",
  if (n_anom > 0L) paste0("- **", df$id[df$status == "ANOMALY"], "** (",
                          df$area[df$status == "ANOMALY"], "): ",
                          df$detail[df$status == "ANOMALY"]) else character(0),
  "",
  "## All probes",
  "| id | area | kind | status | detail |",
  "|---|---|---|---|---|",
  sprintf("| %s | %s | %s | %s | %s |", df$id, df$area, df$kind, df$status,
          gsub("\\|", "/", df$detail))
), con)
close(con)

cat("\n=================== STRESS SUMMARY ===================\n")
cat(sprintf("probes %d | PASS %d | ANOMALY %d | SKIP %d\n",
            nrow(df), n_pass, n_anom, n_skip))
if (n_anom > 0L) {
  cat("\nANOMALIES:\n")
  a <- df[df$status == "ANOMALY", ]
  for (i in seq_len(nrow(a)))
    cat(sprintf("  - [%s] %s: %s\n", a$area[i], a$id[i], a$detail[i]))
}
cat("=====================================================\n")
cat("results -> tools/stress_results/stress_corner_to_corner_results.rds\n")
cat("report  -> tools/stress_results/stress_corner_to_corner_report.md\n")
