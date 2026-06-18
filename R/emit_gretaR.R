# =============================================================================
# emit_gretaR.R -- the gretaR backend (torch NUTS), driven OUT OF PROCESS.
#
# gretaR (torch) and greta (TensorFlow) export the same symbols and cannot
# co-load, so flexyBayes never loads gretaR in its own session: emit_gretaR()
# lowers the flexyBayes IR (fb_terms) into a worker spec, runs the worker
# (inst/gretaR_worker.R, gretaR-only process, model_from_arrays + torch NUTS),
# and wraps the returned posterior::draws_array into a `flexybayes_gretaR` fit.
#
# Activation contract: inst/backend_contract.md. Scope: GLM (gaussian/binomial/
# poisson) + a single random intercept -- the classes gretaR is validated for
# (greta-benchmark/). Structured covariance (vm/ped/fa/us/ar1) is refused by
# .capability_gretaR(): gretaR's catch-up boundary. A performance advisory fires
# on the hierarchical path (gretaR is markedly slower there; see benchmark).
# =============================================================================

# Minimum gretaR version carrying model_from_arrays() (the re-entrant,
# deparse-free builder) + the GB1 NUTS fix. Below this floor the backend is
# dormant (the slot probe). Dev override: options(flexyBayes.gretaR_home=<src>).
.GRETAR_VERSION_FLOOR <- "0.3.0.9000"

# ---- capability predicate ----------------------------------------------------
# fb is an fb_terms IR. Allow GLM + random intercept; refuse structured cov.
.capability_gretaR <- function(fb) {
  rt <- fb$random_terms %||% list()
  has_structured <- any(vapply(
    rt,
    function(t) (t$type %||% "") %in% .STRUCTURED_COV_TYPES,
    logical(1L)
  ))
  if (has_structured) {
    return("gretaR_cannot_represent_structured_cov")
  }
  ok_random <- vapply(
    rt,
    function(t) (t$type %||% "") %in% c("simple", "ide", "id"),
    logical(1L)
  )
  if (length(rt) && !all(ok_random)) {
    return("gretaR_random_term_type_unsupported")
  }
  TRUE
}

# ---- the emit function (signature mirrors emit_greta) ------------------------
emit_gretaR <- function(
  fb,
  data,
  known_matrices = list(),
  weights = NULL,
  n_samples = 1000,
  warmup = 500,
  chains = 4,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  verbose = TRUE,
  mcmc_verbose = TRUE,
  return_code = FALSE,
  the_call = NULL,
  fixed = NULL,
  random = NULL,
  rcov = NULL,
  family = NULL,
  link = NULL,
  data_name = NULL,
  seed = NULL,
  ...
) {
  .gretaR_resolve_or_stop() # version floor / dev-home probe

  fam <- fb$family
  if (!fam %in% c("gaussian", "binomial", "binary", "poisson")) {
    stop(.fb_refusal_condition(
      reason_code = "gretaR_family_unsupported",
      message = paste0(
        "backend = \"gretaR\" supports gaussian / binomial / ",
        "poisson; got \"",
        fam,
        "\"."
      ),
      family_class = "flexybayes_gretaR_refusal"
    ))
  }

  # --- lower the IR into a numeric worker spec --------------------------------
  y <- as.numeric(data[[fb$response]])
  ff <- stats::as.formula(fixed %||% paste(fb$response, "~ 1"))
  if (length(ff) == 3L) {
    ff <- .brms_fixed_only_formula(ff)
  }
  X <- stats::model.matrix(ff, data = data)
  spec <- list(
    family = if (fam == "binary") "binomial" else fam,
    X = X,
    y = y,
    random = NULL,
    canonical_names = colnames(X),
    priors = list(
      fixed_sd = prior_fixed_sd,
      sigma_scale = max(prior_vc_sd * 5, 5),
      tau_scale = max(prior_vc_sd * 2, 2)
    ),
    mcmc = list(
      n_samples = as.integer(n_samples),
      warmup = as.integer(warmup),
      chains = as.integer(chains),
      seed = as.integer(seed %||% 1L),
      n_threads = 1L
    )
  )

  rt <- fb$random_terms %||% list()
  if (length(rt) >= 1L) {
    gv <- rt[[1L]]$var
    if (is.null(gv) || !gv %in% names(data)) {
      stop(.fb_refusal_condition(
        reason_code = "gretaR_random_group_not_in_data",
        message = paste0(
          "gretaR random intercept needs grouping factor '",
          gv,
          "' in the data."
        ),
        family_class = "flexybayes_gretaR_refusal"
      ))
    }
    g <- as.factor(data[[gv]])
    J <- nlevels(g)
    Z <- matrix(0, length(g), J)
    Z[cbind(seq_along(g), as.integer(g))] <- 1
    spec$random <- list(Z = Z, J = J, group = gv)
    if (!isTRUE(getOption("flexyBayes.silence_gretaR_perf_note", FALSE))) {
      message(
        "flexyBayes: the gretaR backend is markedly slower than greta / ",
        "INLA on hierarchical models (see the engine benchmark); for a ",
        "random-effects fit greta or INLA is usually the better choice. ",
        "Silence via options(flexyBayes.silence_gretaR_perf_note = TRUE)."
      )
    }
  }

  if (isTRUE(return_code)) {
    return(structure(
      list(spec = spec, worker = .gretaR_worker_path()),
      class = "flexybayes_gretaR_code"
    ))
  }

  # --- run the out-of-process worker ------------------------------------------
  res <- .gretaR_run_worker(spec)
  if (!identical(res$status, "ok")) {
    stop(
      "gretaR worker failed: ",
      res$message %||% "unknown error",
      call. = FALSE
    )
  }
  draws <- res$draws
  summ <- tryCatch(
    as.data.frame(posterior::summarise_draws(
      draws,
      mean = mean,
      sd = stats::sd,
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk
    )),
    error = function(e) NULL
  )

  structure(
    list(
      draws = draws,
      fb = fb,
      data = data,
      extras = list(
        summary = summ,
        model_info = list(
          n_obs = length(y),
          family = fam,
          link = fb$link,
          n_random = length(rt),
          gretaR_version = res$gretaR_version
        ),
        call_info = list(
          fixed = fixed,
          random = random,
          rcov = rcov,
          data_name = data_name,
          family = family,
          link = link
        ),
        parse_info = list(family = list(family = fb$family, link = fb$link)),
        param_names = posterior::variables(draws),
        run_time = res$t_total
      )
    ),
    class = c("flexybayes_gretaR", "list")
  )
}

# ---- worker plumbing ---------------------------------------------------------
.gretaR_worker_path <- function() {
  p <- system.file("gretaR_worker.R", package = "flexyBayes")
  if (!nzchar(p)) {
    stop("gretaR_worker.R not found in the installed package.", call. = FALSE)
  }
  p
}

.gretaR_run_worker <- function(spec) {
  sp <- tempfile(fileext = ".rds")
  op <- tempfile(fileext = ".rds")
  on.exit(unlink(c(sp, op)), add = TRUE)
  saveRDS(spec, sp)
  env <- c("OMP_NUM_THREADS=1", "NOT_CRAN=true")
  home <- getOption("flexyBayes.gretaR_home", "")
  if (nzchar(home)) {
    env <- c(env, paste0("GRETAR_HOME=", home))
  }
  # The exit code is intentionally ignored: success is defined by the worker
  # having written its output file, checked below.
  suppressWarnings(system2(
    "Rscript",
    c(shQuote(.gretaR_worker_path()), shQuote(sp), shQuote(op)),
    env = env,
    stdout = FALSE,
    stderr = FALSE
  ))
  if (!file.exists(op)) {
    return(list(
      status = "error",
      message = paste0(
        "worker produced no output ",
        "(Rscript / gretaR / torch unavailable?)"
      )
    ))
  }
  readRDS(op)
}

# Version-floor / dev-home probe. Errors (dormant) if neither an installed
# gretaR at >= the floor nor a dev source home is available.
.gretaR_resolve_or_stop <- function() {
  home <- getOption("flexyBayes.gretaR_home", "")
  if (nzchar(home)) {
    return(invisible(TRUE))
  }
  if (nzchar(system.file(package = "gretaR"))) {
    v <- tryCatch(utils::packageVersion("gretaR"), error = function(e) NULL)
    if (!is.null(v) && v >= .GRETAR_VERSION_FLOOR) {
      return(invisible(TRUE))
    }
    stop(.fb_refusal_condition(
      reason_code = "gretaR_below_version_floor",
      message = paste0(
        "the gretaR backend needs gretaR >= ",
        .GRETAR_VERSION_FLOOR,
        " (model_from_arrays + the NUTS ",
        "fix); installed: ",
        as.character(v),
        ". Install a newer gretaR or set ",
        "options(flexyBayes.gretaR_home = <source dir>)."
      ),
      family_class = "flexybayes_gretaR_refusal"
    ))
  }
  stop(.fb_refusal_condition(
    reason_code = "gretaR_not_installed",
    message = paste0(
      "the gretaR backend needs the gretaR package (>= ",
      .GRETAR_VERSION_FLOOR,
      ") or ",
      "options(flexyBayes.gretaR_home = <source dir>)."
    ),
    family_class = "flexybayes_gretaR_refusal"
  ))
}

# ---- canonical-name mapper (draws are already canonical from the worker) -----
.mapper_gretaR_real <- function(fit, fb_terms) {
  # model_from_arrays(names=) put canonical tokens straight onto the draws;
  # no relabel and no precision->SD transform are needed.
  list(map = character(0), transform = list(), source = "registry")
}

# ---- triangulate() draws extractor ------------------------------------------
#' @export
fb_as_draws_simple.flexybayes_gretaR <- function(fit, ...) {
  if (is.null(fit$draws) || !inherits(fit$draws, "draws_array")) {
    stop(
      "flexybayes_gretaR fit carries no posterior::draws_array.",
      call. = FALSE
    )
  }
  m <- posterior::as_draws_matrix(fit$draws)
  cols <- colnames(m)
  stats::setNames(
    lapply(seq_len(ncol(m)), function(j) as.numeric(m[, j])),
    cols
  )
}

# ---- minimal post-fit surface ------------------------------------------------
#' @export
print.flexybayes_gretaR <- function(x, ...) {
  mi <- x$extras$model_info
  cat("-- flexyBayes: gretaR fit (torch NUTS, out-of-process) -----\n")
  cat(
    "  family: ",
    mi$family,
    " | obs: ",
    mi$n_obs,
    " | gretaR ",
    as.character(mi$gretaR_version),
    "\n",
    sep = ""
  )
  cat(
    "  params: ",
    paste(x$extras$param_names, collapse = ", "),
    "\n",
    sep = ""
  )
  cat("------------------------------------------------------------\n")
  invisible(x)
}

#' @export
summary.flexybayes_gretaR <- function(object, ...) {
  cat("flexyBayes gretaR fit -- posterior summary\n")
  print(object$extras$summary)
  invisible(object$extras$summary)
}

#' @export
coef.flexybayes_gretaR <- function(object, ...) {
  s <- object$extras$summary
  stats::setNames(s$mean, s$param %||% s$variable)
}
