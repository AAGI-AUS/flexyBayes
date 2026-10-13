# data-raw/asreml_fa_reference.R
#
# ASReml 4.2 reference fits of the factor-analytic genotype-by-environment
# model `fa(env, k):gen` for k = 1 and k = 2, on three public trials, stored
# so the tests can compare flexyBayes with ASReml without an ASReml licence.
#
#   20_fa_besag_met      agridat::besag.met, all 64 genotypes, 6 counties
#   21_fa_besag_met_30   the same, first 30 genotypes by level order
#   22_fa_dasilva_maize  agridat::dasilva.maize, 55 genotypes, 9 environments,
#                        one residual variance per environment
#
# Each record holds, per k, ASReml's variance components, the genetic
# covariance G = Lambda Lambda' + Psi, the genetic correlations, the REML
# log-likelihood and whether ASReml converged. G and the correlations do not
# depend on how the loadings are rotated or signed; the loadings do.
#
# Run from the package root on a clean tree, with a licensed ASReml:
#
#   FLEXYBAYES_ASREML_LIB=/path/to/library \
#     Rscript --vanilla data-raw/asreml_fa_reference.R
#
# The records name the commit they were run from, and that commit is
# declared in inst/validation/asreml-4.2/ANCHORS.md
# (tests/testthat/test-oracle-fixture-provenance.R). `tarball_sha256` is the
# SHA-256 of `git archive --format=tar <commit>`, the source the run read.
# This one file holds the data, the fits and the extraction, so it is hashed
# into all three of the scenario, comparison and runner fields.
#
# FLEXYBAYES_ALLOW_DIRTY=1 runs on a dirty tree and marks the records dirty,
# which the provenance test rejects.

asreml_lib <- Sys.getenv("FLEXYBAYES_ASREML_LIB", unset = "")
if (nzchar(asreml_lib)) {
  if (!dir.exists(asreml_lib)) {
    stop(
      "FLEXYBAYES_ASREML_LIB is set to a directory that does not exist: ",
      asreml_lib,
      call. = FALSE
    )
  }
  .libPaths(c(asreml_lib, .libPaths()))
}
suppressMessages(library(asreml))

script_path <- "data-raw/asreml_fa_reference.R"
out_dir <- "inst/validation/asreml-4.2"
stopifnot(file.exists(script_path), dir.exists(out_dir))

# --- 1. provenance ------------------------------------------------------- #

git <- function(...) {
  system2("git", c(...), stdout = TRUE)
}
source_sha <- git("rev-parse", "HEAD")
dirty <- length(git("status", "--porcelain")) > 0L
if (dirty && !identical(Sys.getenv("FLEXYBAYES_ALLOW_DIRTY"), "1")) {
  stop(
    "The working tree is not clean, so the records could not name the code ",
    "they were produced by. Commit first.",
    call. = FALSE
  )
}
archive_sha256 <- system(
  paste0("git archive --format=tar ", source_sha, " | shasum -a 256"),
  intern = TRUE
)
archive_sha256 <- sub(" .*$", "", archive_sha256)
script_sha256 <- digest::digest(file = script_path, algo = "sha256")

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
provenance <- list(
  scenarios_sha256 = script_sha256,
  compare_sha256 = script_sha256,
  runner_sha256 = script_sha256,
  dirty = dirty,
  run_command = paste0("Rscript --vanilla ", script_path),
  blas = extSoftVersion()[["BLAS"]],
  lapack = La_library(),
  rng_kind = paste(RNGkind(), collapse = " / "),
  platform = R.version$platform
)
versions <- list(
  r = R.version.string,
  asreml = as.character(utils::packageVersion("asreml")),
  agridat = as.character(utils::packageVersion("agridat"))
)

# --- 2. the data --------------------------------------------------------- #

agridat_data <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "agridat", envir = e)
  e[[name]]
}

besag_met <- function(n_gen = NULL) {
  d <- agridat_data("besag.met")
  d$env <- factor(d$county)
  d$gen <- factor(d$gen)
  d <- d[!is.na(d$yield), ]
  if (!is.null(n_gen)) {
    d <- d[d$gen %in% levels(d$gen)[seq_len(n_gen)], ]
    d$env <- droplevels(d$env)
    d$gen <- droplevels(d$gen)
  }
  d <- d[order(d$env, d$gen), c("env", "gen", "rep", "yield")]
  rownames(d) <- NULL
  d
}

dasilva_maize <- function() {
  d <- agridat_data("dasilva.maize")
  d <- d[order(d$env, d$gen), c("env", "gen", "rep", "yield")]
  rownames(d) <- NULL
  d
}

designs <- list(
  "20_fa_besag_met" = list(
    label = "agridat::besag.met, all genotypes",
    data = besag_met(),
    sectioned = FALSE
  ),
  "21_fa_besag_met_30" = list(
    label = "agridat::besag.met, first 30 genotypes by level order",
    data = besag_met(30L),
    sectioned = FALSE
  ),
  "22_fa_dasilva_maize" = list(
    label = "agridat::dasilva.maize, one residual variance per environment",
    data = dasilva_maize(),
    sectioned = TRUE
  )
)

# --- 3. one ASReml fit and its identified covariance --------------------- #

fit_fa <- function(d, k, sectioned) {
  f_rand <- stats::as.formula(paste0("~ fa(env, ", k, "):gen"))
  warn <- character(0)
  t0 <- proc.time()
  m <- withCallingHandlers(
    if (sectioned) {
      asreml(
        fixed = yield ~ env,
        random = f_rand,
        residual = ~ dsum(~ units | env),
        data = d,
        trace = FALSE,
        maxit = 200
      )
    } else {
      asreml(
        fixed = yield ~ env,
        random = f_rand,
        data = d,
        trace = FALSE,
        maxit = 200
      )
    },
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # One further pass from the converged values, as is usual for an fa model,
  # where the last iteration often still moves a component by more than 1%.
  warn_first <- warn
  warn <- character(0)
  m <- withCallingHandlers(
    update(m),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  elapsed <- unname((proc.time() - t0)[["elapsed"]])

  vc <- summary(m)$varcomp
  levs <- levels(d$env)
  prefix <- paste0("fa(env, ", k, "):gen!")
  mat_l <- matrix(
    0, length(levs), k,
    dimnames = list(levs, paste0("fa", seq_len(k)))
  )
  vec_psi <- stats::setNames(numeric(length(levs)), levs)
  for (lev in levs) {
    vec_psi[[lev]] <- vc[paste0(prefix, lev, "!var"), "component"]
    for (j in seq_len(k)) {
      mat_l[lev, j] <- vc[paste0(prefix, lev, "!fa", j), "component"]
    }
  }
  mat_g <- mat_l %*% t(mat_l) + diag(vec_psi, nrow = length(levs))
  dimnames(mat_g) <- list(levs, levs)

  list(
    k = k,
    converged = isTRUE(m$converge),
    loglik = m$loglik,
    warnings_first_pass = warn_first,
    warnings = warn,
    elapsed_sec = elapsed,
    varcomp = data.frame(
      term = rownames(vc),
      component = vc$component,
      se = vc$std.error,
      bound = as.character(vc$bound),
      stringsAsFactors = FALSE
    ),
    loadings = mat_l,
    psi = vec_psi,
    G = mat_g,
    cor = stats::cov2cor(mat_g)
  )
}

# --- 4. the records ------------------------------------------------------ #

for (nm in names(designs)) {
  des <- designs[[nm]]
  d <- des$data
  record <- list(
    scenario = nm,
    label = des$label,
    run_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    data = d,
    data_rows = nrow(d),
    data_sha256 = digest::digest(d, algo = "sha256"),
    model = list(
      fixed = "yield ~ env",
      random = "~ fa(env, k):gen",
      residual = if (des$sectioned) "~ dsum(~ units | env)" else "~ units"
    ),
    versions = versions,
    source_sha = source_sha,
    tarball_sha256 = archive_sha256,
    provenance = provenance,
    fits = list(
      fa1 = fit_fa(d, 1L, des$sectioned),
      fa2 = fit_fa(d, 2L, des$sectioned)
    )
  )
  saveRDS(record, file.path(out_dir, paste0(nm, ".rds")))
  print(paste0(
    nm, ": fa1 converged = ", record$fits$fa1$converged,
    ", loglik = ", round(record$fits$fa1$loglik, 3),
    "; fa2 converged = ", record$fits$fa2$converged,
    ", loglik = ", round(record$fits$fa2$loglik, 3)
  ))
}
print(paste0("source ", source_sha, ", dirty = ", dirty))
