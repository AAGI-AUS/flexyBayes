# run_one.R -- run one paired scenario and write its manifest row.
#
#   FB_LIB=<lib> FB_SCENARIO=<id> FB_OUT=<dir> R --file=run_one.R
#
# One scenario per process, deliberately: these fits are memory-heavy and a
# single long-lived session accumulating all of them is what killed two earlier
# whole-suite runs on this machine.

.libPaths(path.expand(c(Sys.getenv("FB_LIB"),
                        "~/Library/R/arm64/4.5/library", .libPaths())))
suppressMessages(library(asreml))
suppressMessages(library(flexyBayes))
# Scenario 19 reads its substrate from lme4. Loaded unconditionally so a
# missing dependency fails at the top of the run rather than inside one
# scenario, where it would look like a modelling failure.
suppressMessages(library(lme4))

here <- Sys.getenv("FB_ORACLE_DIR")
source(file.path(here, "compare.R"))
source(file.path(here, "scenarios.R"))

id <- Sys.getenv("FB_SCENARIO")
out_dir <- Sys.getenv("FB_OUT")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sc <- scenario_registry()[[id]]
if (is.null(sc)) stop("unknown scenario: ", id, call. = FALSE)

d <- sc$data()
t0 <- Sys.time()

asr <- tryCatch(sc$asreml(d), error = function(e) e)
fb <- tryCatch(sc$flexybayes(d), error = function(e) e)

record <- list(
  scenario = id,
  label = sc$label,
  registry_row = sc$registry_row,
  run_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  data_rows = nrow(d),
  data_sha256 = data_digest(d),
  factor_levels = lapply(
    stats::setNames(sc$factor_cols, sc$factor_cols),
    function(cn) levels(factor(d[[cn]]))
  ),
  versions = list(
    r = R.version.string,
    asreml = as.character(utils::packageVersion("asreml")),
    flexybayes = as.character(utils::packageVersion("flexyBayes")),
    inla = if (requireNamespace("INLA", quietly = TRUE))
      as.character(utils::packageVersion("INLA")) else NA_character_,
    agridat = as.character(utils::packageVersion("agridat"))
  ),
  source_sha = Sys.getenv("FB_SOURCE_SHA"),
  tarball_sha256 = Sys.getenv("FB_TARBALL_SHA"),
  # Review A21: a future reader has to be able to establish which code
  # produced this number without trusting the sentence that shipped with it.
  # The commit and the tarball hash pin the package; these pin the
  # comparison itself, the numerical library underneath it, and the command.
  provenance = list(
    scenarios_sha256 = Sys.getenv("FB_SCENARIOS_SHA"),
    compare_sha256 = Sys.getenv("FB_COMPARE_SHA"),
    runner_sha256 = Sys.getenv("FB_RUNNER_SHA"),
    dirty = identical(Sys.getenv("FB_DIRTY"), "TRUE"),
    run_command = Sys.getenv("FB_RUN_CMD"),
    blas = tryCatch(
      extSoftVersion()[["BLAS"]],
      error = function(e) NA_character_
    ),
    lapack = tryCatch(La_library(), error = function(e) NA_character_),
    rng_kind = paste(RNGkind(), collapse = " / "),
    platform = R.version$platform
  )
)

if (inherits(asr, "error") || inherits(fb, "error")) {
  record$status <- "not_compared"
  record$reason <- paste0(
    if (inherits(asr, "error")) {
      paste0("asreml did not fit: ", conditionMessage(asr))
    } else {
      ""
    },
    if (inherits(fb, "error")) {
      paste0("flexybayes did not fit: ", conditionMessage(fb))
    } else {
      ""
    }
  )
  saveRDS(record, file.path(out_dir, paste0(id, ".rds")))
  cat(id, "\tNOT_COMPARED\t", record$reason, "\n", sep = "")
  quit(save = "no", status = 0)
}

s <- summary(fb)
record$asreml_converged <- isTRUE(asr$converge)
# Tri-state, and never defaulted. This line used to read
# isTRUE(s$converge$ok %||% TRUE) against a field that did not exist, and
# NULL %||% TRUE is TRUE, so it was a CONSTANT: no input produced FALSE and
# all eight banked fixtures recorded a convergence that had not been
# measured (audit A25). `ok` now exists and is TRUE / FALSE / NA.
record$flexybayes_converged <- s$converge$ok
record$asreml_loglik <- asr$loglik

# ---- Layer 1 ---- #
record$structural <- structural_checks(asr, s, d, sc$factor_cols)

# ---- Layer 2 ---- #
v_asr <- asr_varcomp(asr)
se_asr <- asr_varcomp_se(asr)
v_fb <- fb_varcomp(s)
record$raw <- list(asreml_varcomp = v_asr, asreml_se = se_asr,
                   flexybayes_varcomp_as_variance = v_fb)

checks <- list()
for (r in sc$ratios) {
  if (r$asreml %in% names(v_asr) && r$flexybayes %in% names(v_fb)) {
    checks[[length(checks) + 1L]] <- cmp_variance_ratio(r, v_asr, v_fb, se_asr)
  }
}

fixed_map <- lapply(
  intersect(rownames(summary(asr, coef = TRUE)$coef.fixed), s$fixed$term),
  function(tm) list(asreml = tm, flexybayes = tm)
)
checks <- c(checks, cmp_fixed(asr, s, fixed_map))

for (b in sc$blup) {
  ra <- asreml::summary.asreml(asr, coef = TRUE)$coef.random
  keep <- grep(paste0("^", b$asreml, "_"), rownames(ra))
  blup_asr <- stats::setNames(
    ra[keep, "solution"],
    sub(paste0("^", b$asreml, "_"), "", rownames(ra)[keep])
  )
  # ranef() returns a named list of data frames, one per group, not a
  # single frame. Read the contract rather than the shape it happened to
  # print in.
  rf_all <- ranef(fb)
  rf <- if (is.data.frame(rf_all)) {
    rf_all[rf_all$group == b$flexybayes, , drop = FALSE]
  } else {
    rf_all[[b$flexybayes]]
  }
  blup_fb <- if (is.null(rf) || !nrow(rf)) {
    stats::setNames(numeric(0), character(0))
  } else {
    stats::setNames(rf$estimate, as.character(rf$level))
  }
  checks[[length(checks) + 1L]] <- cmp_blup(blup_asr, blup_fb, b$label)
}

# Correlations. ASReml parameterises us() as variances and a covariance, so
# the correlation is reconstructed from the three named parameters; flexyBayes
# reports it directly. Comparing the reconstruction against the reported value
# is the point -- it is the one quantity invariant to the two systems'
# different parameterisations, and the one a breeder reads.
for (r in sc$correlations %||% list()) {
  nm <- r$asreml
  if (all(nm %in% names(v_asr)) && r$flexybayes %in% names(v_fb)) {
    r_asr <- unname(v_asr[[nm[["cov"]]]] /
                      sqrt(v_asr[[nm[["v1"]]]] * v_asr[[nm[["v2"]]]]))
    checks[[length(checks) + 1L]] <- cmp_correlation(
      r_asr, unname(v_fb[[r$flexybayes]]), r$label
    )
  } else {
    checks[[length(checks) + 1L]] <- list(
      quantity = paste0("correlation: ", r$label), pass = NA,
      note = "not asserted -- a named component was absent on one side"
    )
  }
}

# ---- Layer 3 ---- #
#
# SPEC.md Amendment 5. Runs only for a scenario that declares the pair, so
# every earlier scenario's record is unchanged in shape and in content.
if (is.function(sc$asreml_protected)) {
  asr_p <- tryCatch(sc$asreml_protected(d), error = function(e) e)
  fb_p <- tryCatch(sc$flexybayes_protected(d), error = function(e) e)
  v1 <- list()
  if (inherits(asr_p, "error") || inherits(fb_p, "error")) {
    v1$status <- "not_compared"
    v1$reason <- paste0(
      if (inherits(asr_p, "error")) {
        paste0("asreml protected fit failed: ", conditionMessage(asr_p))
      } else {
        ""
      },
      if (inherits(fb_p, "error")) {
        paste0("flexybayes protected fit failed: ", conditionMessage(fb_p))
      } else {
        ""
      }
    )
  } else {
    s_p <- summary(fb_p)
    moments_full <- .fb_inla_variance_moments(fb)
    moments_prot <- .fb_inla_variance_moments(fb_p)
    v1$status <- "compared"
    v1$inla_moments_full <- moments_full
    v1$inla_moments_protected <- moments_prot
    v1$protected_converged <- s_p$converge$ok
    v1$protected_asreml_converged <- isTRUE(asr_p$converge)
    v1$protected_asreml_loglik <- asr_p$loglik

    v1_checks <- c(
      cmp_fixed_v1(asr_p, s_p, "protected_only"),
      cmp_fixed_v1(asr, s, "protected_plus_candidates")
    )
    if (!is.null(sc$varcomp_map)) {
      v1_checks <- c(
        v1_checks,
        cmp_varcomp_v1(asr_p, moments_prot, sc$varcomp_map, "protected_only"),
        cmp_varcomp_v1(asr, moments_full, sc$varcomp_map,
                       "protected_plus_candidates")
      )
    }
    for (b in sc$blup_v1 %||% list()) {
      v1_checks[[length(v1_checks) + 1L]] <- cmp_blup_v1(
        .fb_blup_asreml(asr, b$asreml),
        .fb_blup_flexybayes(fb, b$flexybayes),
        paste0(b$label, " (protected_plus_candidates)")
      )
    }
    v1$checks <- v1_checks
  }
  record$layer3 <- v1
}

record$checks <- checks
record$elapsed_sec <-
  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
record$status <- "compared"

saveRDS(record, file.path(out_dir, paste0(id, ".rds")))

np <- sum(vapply(checks, function(c) isTRUE(c$pass), logical(1)))
nf <- sum(vapply(checks, function(c) identical(c$pass, FALSE), logical(1)))
nn <- sum(vapply(checks, function(c) is.na(c$pass), logical(1)))
ns <- sum(vapply(record$structural, function(c) isTRUE(c$pass), logical(1)))
nsf <- sum(vapply(
  record$structural,
  function(c) identical(c$pass, FALSE),
  logical(1)
))
l3 <- record$layer3$checks %||% list()
l3p <- sum(vapply(l3, function(c) isTRUE(c$pass), logical(1)))
l3f <- sum(vapply(l3, function(c) identical(c$pass, FALSE), logical(1)))
l3n <- sum(vapply(l3, function(c) is.na(c$pass), logical(1)))
cat(sprintf(
  paste0(
    "%s\tstructural %d pass / %d fail\tnumeric %d pass / %d fail / %d ",
    "not-asserted\tlayer3 %d / %d / %d\t%.0fs\n"
  ),
  id, ns, nsf, np, nf, nn, l3p, l3f, l3n, record$elapsed_sec
))
