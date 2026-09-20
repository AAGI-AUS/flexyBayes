# tools/check_effectsurf_surfaces.R -- a response surface that effectsurf
# draws from a flexyBayes fit must be that fit's posterior.
#
# effectsurf draws the surface of population-level predictions over two
# covariates of a fitted model through marginaleffects, which reaches a
# flexyBayes fit through the package's get_coef(), get_vcov() and
# get_predict() methods (R/marginaleffects_support.R). effectsurf is not on
# CRAN, so this is a release-gate script rather than a package test, and
# flexyBayes does not depend on it.
#
# The script simulates a known surface, y = 3 + 0.8 x1 - 1.2 x2 + 0.15 x1 x2
# with a site effect, fits it on each engine, asks effectsurf for the
# surface over (x1, x2), and compares every grid point with the posterior of
# the linear predictor computed here: with a design matrix written here, and
# from the engine's own output rather than from the package's summaries.
#
#   INLA  runs on a standardised response. The relation between the data
#         and the coordinates INLA ran on is read off the engine's own data
#         (`$inla$.args$data`), and the posterior mean of the fixed effects
#         on those coordinates is the weight-averaged conditional mean over
#         INLA's hyperparameter configurations; the posterior spread comes
#         from INLA's own draws on the same coordinates.
#   brms  the draws of the fixed effects held by the brmsfit.
#
# Checks, per engine:
#   grid   the surface is the requested 25 x 20 lattice over x1 in [0, 10]
#          and x2 in [0, 5];
#   mean   its estimate equals the posterior mean of the linear predictor
#          at every grid point;
#   band   its interval half-width equals 1.96 posterior standard
#          deviations of the linear predictor at every grid point;
#   truth  the simulated surface lies within 4 posterior standard
#          deviations of the estimate at every grid point.
# The departure of the normal band from the posterior's own 2.5 and 97.5
# per cent quantiles is reported and not tested: effectsurf's interval is a
# normal approximation by construction.
#
# Tolerances, from the spread measured on this fixture before the checks
# were written:
#   INLA mean   1e-3 posterior sd. effectsurf's estimate is the package's
#               stored posterior mean, which INLA integrates numerically; it
#               sat 2.7e-5 sd from the configuration mean.
#   INLA band   relative, 4 x sqrt(r^2 + 1 / (2 (n - 1))), where r is the
#               largest relative spread of effectsurf's half-width over five
#               identical calls in this run (its covariance is a Monte Carlo
#               estimate from 2,000 draws) and n = 4,000 draws are taken
#               here. Measured: r = 0.019, half-widths within 0.024.
#   brms mean   1e-8 posterior sd: the package's coefficients are the mean
#               of the same draws, so any difference is floating point.
#   brms band   1e-6, relative. The package's covariance is the covariance
#               of the same draws, but marginaleffects forms the delta
#               method's Jacobian by finite differences, whose rounding error
#               over the step reaches the half-width: measured 2.4e-8.
#
# Run from the package root:
#
#   Rscript tools/check_effectsurf_surfaces.R                     # INLA and brms
#   Rscript tools/check_effectsurf_surfaces.R --engines=inla      # one engine
#   Rscript tools/check_effectsurf_surfaces.R --installed         # the installed flexyBayes
#   Rscript tools/check_effectsurf_surfaces.R --out=<dir>         # where the grid tables go
#   Rscript tools/check_effectsurf_surfaces.R --isolation=process # fits in a separate process
#
# Without --installed the package is loaded from the working tree with
# pkgload. Exit status 0 when every check passes, 1 when a check fails, and
# 2 when the run cannot be made (a missing package or engine), which is
# never a pass.

if (!file.exists("DESCRIPTION")) {
  stop("Run this script from the package root.", call. = FALSE)
}

# ---- arguments and requirements ------------------------------------- #

args_vec <- commandArgs(trailingOnly = TRUE)
.arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args_vec, value = TRUE)
  if (length(hit) == 0L) default else sub(paste0("^--", name, "="), "", hit[[1L]])
}
engines <- strsplit(.arg_value("engines", "inla,brms"), ",", fixed = TRUE)[[1L]]
out_dir <- .arg_value("out", file.path(tempdir(), "effectsurf-gate"))
isolation <- .arg_value("isolation", "none")
installed <- "--installed" %in% args_vec

.cannot_run <- function(...) {
  cat("[effectsurf-gate] CANNOT RUN: ", ..., "\n", sep = "")
  quit(status = 2L, save = "no")
}

needed <- c("effectsurf", "marginaleffects", if ("inla" %in% engines) "INLA",
            if ("brms" %in% engines) c("brms", "posterior"),
            if (!installed) "pkgload")
absent <- needed[!vapply(needed, requireNamespace, logical(1L), quietly = TRUE)]
if (length(absent) > 0L) {
  .cannot_run("not installed: ", paste(absent, collapse = ", "), ".")
}
if (!all(engines %in% c("inla", "brms"))) {
  .cannot_run("--engines takes inla, brms or both.")
}

if (installed) {
  suppressPackageStartupMessages(library(flexyBayes))
} else {
  suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
}
options(
  flexyBayes.silence_default_prior_note = TRUE,
  flexyBayes.silence_uniform_inla_approx = TRUE
)
if (!identical(isolation, "none")) {
  options(flexyBayes.isolation = isolation)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("[effectsurf-gate] flexyBayes ", as.character(utils::packageVersion("flexyBayes")),
    " from ", find.package("flexyBayes"), "; effectsurf ",
    as.character(utils::packageVersion("effectsurf")), "; marginaleffects ",
    as.character(utils::packageVersion("marginaleffects")), "; isolation ",
    getOption("flexyBayes.isolation", "none"), "\n", sep = "")

# ---- the known surface ----------------------------------------------- #

.truth <- function(x1, x2) 3 + 0.8 * x1 - 1.2 * x2 + 0.15 * x1 * x2

set.seed(20260919L)
n_site <- 12L
per_site <- 10L
dat <- data.frame(
  site = factor(rep(sprintf("s%02d", seq_len(n_site)), each = per_site)),
  x1 = stats::runif(n_site * per_site, 0, 10),
  x2 = stats::runif(n_site * per_site, 0, 5)
)
dat$y <- .truth(dat$x1, dat$x2) +
  stats::rnorm(n_site, sd = 0.7)[as.integer(dat$site)] +
  stats::rnorm(nrow(dat), sd = 0.5)

x1_grid <- seq(0, 10, length.out = 25L)
x2_grid <- seq(0, 5, length.out = 20L)
terms_here <- c("(Intercept)", "x1", "x2", "x1:x2")

# The design of the fixed effects at the grid, written here.
.design <- function(x1, x2) {
  cbind(1, x1, x2, x1 * x2)
}

.surface <- function(fit) {
  es <- effectsurf::surf_prediction(
    fit, x = "x1", y = "x2", x_range = c(0, 10), y_range = c(0, 5),
    x_length = 25L, y_length = 20L, level = 0.95
  )
  out <- as.data.frame(effectsurf::surf_data(es))
  out[order(out$x2, out$x1), , drop = FALSE]
}

# ---- the posterior of the linear predictor, per engine --------------- #

# How the data map onto the coordinates the engine ran on: y = a + b y*, and
# each covariate divided by a constant. Read off the engine's own data.
.engine_map <- function(engine_data, response_name) {
  y_star <- as.numeric(engine_data[[response_name]])
  ab <- stats::coef(stats::lm(dat$y ~ y_star))
  gap <- max(abs(dat$y - ab[[1L]] - ab[[2L]] * y_star))
  divisor <- vapply(c("x1", "x2"), function(v) {
    ratio <- dat[[v]] / as.numeric(engine_data[[v]])
    if (diff(range(ratio)) > 1e-12 * max(abs(ratio))) {
      .cannot_run("covariate ", v, " is not a constant multiple of the data.")
    }
    ratio[[1L]]
  }, numeric(1L))
  if (gap > 1e-8 * stats::sd(dat$y)) {
    .cannot_run("the engine's response is not an affine map of the data's.")
  }
  list(a = ab[[1L]], b = ab[[2L]], divisor = divisor)
}

# The engine-coordinate design at the grid.
.engine_design <- function(x1, x2, map) {
  .design(x1 / map$divisor[["x1"]], x2 / map$divisor[["x2"]])
}

.posterior_inla <- function(fit, grid) {
  engine_data <- fit$inla$.args$data
  response_name <- all.vars(fit$inla$.args$formula[[2L]])
  map <- .engine_map(engine_data, response_name)
  X_star <- .engine_design(grid$x1, grid$x2, map)

  # The posterior mean: configurations weighted by their posterior mass.
  configs <- fit$inla$misc$configs
  offset <- sum(configs$contents$length) - length(configs$config[[1L]]$mean)
  at <- configs$contents$start[match(terms_here, configs$contents$tag)] - offset
  if (anyNA(at)) {
    .cannot_run("INLA's configurations do not hold every fixed effect.")
  }
  log_post <- vapply(configs$config, function(k) k$log.posterior, numeric(1L))
  weight <- exp(log_post - max(log_post))
  weight <- weight / sum(weight)
  beta_star <- Reduce(`+`, Map(
    function(k, w) w * k$improved.mean[at], configs$config, weight
  ))
  mean_eta <- map$a + map$b * drop(X_star %*% beta_star)

  # The spread: INLA's own draws on the same coordinates.
  n_draws <- 4000L
  draws <- INLA::inla.posterior.sample(
    n_draws, fit$inla, num.threads = "1:1", seed = 1L
  )
  rows <- match(paste0(terms_here, ":1"), rownames(draws[[1L]]$latent))
  beta_draws <- vapply(draws, function(s) s$latent[rows, 1L], numeric(4L))
  eta_draws <- map$a + map$b * (X_star %*% beta_draws)
  list(
    mean = mean_eta,
    sd = apply(eta_draws, 1L, stats::sd),
    q = t(apply(eta_draws, 1L, stats::quantile, probs = c(0.025, 0.975))),
    n_draws = n_draws,
    map = map
  )
}

.posterior_brms <- function(fit, grid) {
  engine_data <- fit$brms$data
  map <- .engine_map(engine_data, "y")
  X_star <- .engine_design(grid$x1, grid$x2, map)
  draws <- posterior::as_draws_matrix(fit$brms)
  beta_draws <- t(unclass(draws)[, paste0("b_", c("Intercept", "x1", "x2", "x1:x2")),
                                 drop = FALSE])
  eta_draws <- map$a + map$b * (X_star %*% beta_draws)
  list(
    mean = rowMeans(eta_draws),
    sd = apply(eta_draws, 1L, stats::sd),
    q = t(apply(eta_draws, 1L, stats::quantile, probs = c(0.025, 0.975))),
    n_draws = ncol(beta_draws),
    map = map
  )
}

# ---- the checks --------------------------------------------------------- #

results <- list()
.record <- function(engine, check, pass, detail) {
  results[[length(results) + 1L]] <<- data.frame(
    engine = engine, check = check, pass = pass, detail = detail,
    stringsAsFactors = FALSE
  )
  cat(sprintf("[effectsurf-gate] %-5s %-6s %s  %s\n", engine, check,
              if (pass) "PASS" else "FAIL", detail))
}

# A warning from a fit or from effectsurf is printed, not swallowed.
.noting_warnings <- function(expr, what) {
  withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      cat("[effectsurf-gate] warning from ", what, ": ", conditionMessage(w),
          "\n", sep = "")
      invokeRestart("muffleWarning")
    }
  )
}

for (engine in engines) {
  fit <- .noting_warnings(
    if (identical(engine, "inla")) {
      flexybayes(y ~ x1 * x2, random = ~site, data = dat, backend = "inla",
                 verbose = FALSE)
    } else {
      flexybayes(y ~ x1 * x2, random = ~site, data = dat, backend = "brms",
                 chains = 2L, n_samples = 1000L, warmup = 500L, seed = 1L,
                 verbose = FALSE)
    },
    paste0("the ", engine, " fit")
  )
  surface <- .noting_warnings(.surface(fit), "effectsurf")

  # grid
  lattice <- expand.grid(x1 = x1_grid, x2 = x2_grid)
  grid_ok <- nrow(surface) == nrow(lattice) &&
    isTRUE(all.equal(surface$x1, lattice$x1, tolerance = 1e-12)) &&
    isTRUE(all.equal(surface$x2, lattice$x2, tolerance = 1e-12))
  .record(engine, "grid", grid_ok, sprintf(
    "%d points; x1 %s to %s, x2 %s to %s", nrow(surface),
    format(min(surface$x1)), format(max(surface$x1)),
    format(min(surface$x2)), format(max(surface$x2))
  ))
  if (!grid_ok) {
    next
  }

  post <- if (identical(engine, "inla")) {
    .posterior_inla(fit, surface)
  } else {
    .posterior_brms(fit, surface)
  }

  # mean
  mean_gap <- max(abs(surface$estimate - post$mean) / post$sd)
  mean_tol <- if (identical(engine, "inla")) 1e-3 else 1e-8
  .record(engine, "mean", mean_gap <= mean_tol, sprintf(
    "max |estimate - posterior mean| = %.3g posterior sd (tolerance %.0e); engine map y = %.6g + %.6g y*",
    mean_gap, mean_tol, post$map$a, post$map$b
  ))

  # band
  half <- (surface$conf.high - surface$conf.low) / 2
  half_here <- stats::qnorm(0.975) * post$sd
  band_gap <- max(abs(half / half_here - 1))
  band_tol <- if (identical(engine, "inla")) {
    repeat_half <- vapply(seq_len(5L), function(i) {
      again <- .noting_warnings(.surface(fit), "effectsurf")
      (again$conf.high - again$conf.low) / 2
    }, numeric(nrow(surface)))
    spread <- max(apply(repeat_half, 1L, stats::sd) / rowMeans(repeat_half))
    4 * sqrt(spread^2 + 1 / (2 * (post$n_draws - 1)))
  } else {
    1e-6
  }
  quantile_gap <- max(abs(cbind(surface$conf.low, surface$conf.high) - post$q) / post$sd)
  .record(engine, "band", band_gap <= band_tol, sprintf(
    "max |half-width / (1.96 posterior sd) - 1| = %.3g (tolerance %.3g); bounds sit up to %.3g sd from the posterior quantiles (%d draws)",
    band_gap, band_tol, quantile_gap, post$n_draws
  ))

  # truth
  truth <- .truth(surface$x1, surface$x2)
  truth_gap <- max(abs(surface$estimate - truth) / post$sd)
  .record(engine, "truth", truth_gap <= 4, sprintf(
    "max |estimate - simulated surface| = %.3g posterior sd; %.1f%% of the grid inside the band",
    truth_gap, 100 * mean(truth >= surface$conf.low & truth <= surface$conf.high)
  ))

  table_out <- data.frame(
    x1 = surface$x1, x2 = surface$x2, estimate = surface$estimate,
    conf_low = surface$conf.low, conf_high = surface$conf.high,
    posterior_mean = post$mean, posterior_sd = post$sd,
    posterior_q025 = post$q[, 1L], posterior_q975 = post$q[, 2L],
    simulated = truth
  )
  utils::write.csv(
    table_out, file.path(out_dir, paste0("effectsurf_surface_", engine, ".csv")),
    row.names = FALSE
  )
}

summary_tab <- do.call(rbind, results)
utils::write.csv(summary_tab, file.path(out_dir, "effectsurf_gate_summary.csv"),
                 row.names = FALSE)
n_fail <- sum(!summary_tab$pass)
cat("[effectsurf-gate] ", nrow(summary_tab), " checks, ", n_fail, " failed; tables in ",
    normalizePath(out_dir), "\n", sep = "")
quit(status = if (n_fail == 0L) 0L else 1L, save = "no")
