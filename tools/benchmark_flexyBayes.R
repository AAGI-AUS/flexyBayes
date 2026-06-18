#!/usr/bin/env Rscript

# flexyBayes assessment / benchmark harness
# =========================================
#
# Purpose
# -------
# This script benchmarks the current local flexyBayes source tree with
# a pragmatic "trustworthy instrument" focus:
#
#   1. package/runtime environment;
#   2. parse / IR / gate / prior-translation speed and correctness;
#   3. backend smoke fits through INLA, greta, and brms where available;
#   4. triangulation behaviour;
#   5. selected code-generation and review-code surfaces;
#   6. representation benchmarks (per the design spec):
#      generated_code_size_bytes / pre_fit_object_size_bytes /
#      code_generation_time_sec per indexed term class;
#   7. an honest report that records failures, skips, and timeouts.
#
# The filename intentionally carries no version stamp; the harness
# runs against whichever version is currently installed (read from
# DESCRIPTION at run time and recorded in environment.csv).
#
# Style constraints
# -----------------
# The implementation deliberately uses base R + data.table style:
# no tidyverse, no future framework, no benchmark-specific package
# dependency beyond data.table. Multicore execution uses base
# parallel::mclapply() and child Rscript processes.
#
# Safety / time budgeting
# -----------------------
# Heavy Bayesian backends can stall during Python/TensorFlow or Stan
# setup. Each backend task runs in a separate R process via system2()
# with a per-task timeout. That way a stalled greta / brms / INLA
# run is recorded as a timeout rather than blocking the whole study.
#
# Typical runtime on the local Mac with all dependencies installed:
#   - smoke profile:    5-15 minutes
#   - standard profile: 15-45 minutes
#   - full profile:     30-120 minutes, depending mostly on brms/Stan
#
# The main entry point below defaults to the standard profile. For
# the requested run, use:
#
#   Rscript tools/benchmark_flexyBayes.R \
#     --profile full --max-cores 8 --budget-min 180 --include-brms-fit true
#
# Outputs
# -------
# A timestamped directory under benchmark_results/ containing:
#   - environment.csv
#   - microbench_results.csv
#   - child_task_results.csv
#   - summary_by_group.csv
#   - raw_results.rds
#   - benchmark_report.md

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
      "Package 'data.table' is required for this benchmark harness.",
      call. = FALSE
    )
  }
})

dt <- data.table::data.table
`%chin%` <- data.table::`%chin%`
`%||%` <- function(x, y) if (is.null(x)) y else x

# Keep per-process numerical libraries from oversubscribing all cores
# while the harness itself uses multicore parallelism.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  TF_CPP_MIN_LOG_LEVEL = "2"
)

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) {
      i <- i + 1L
      next
    }
    nm <- sub("^--", "", key)
    val <- TRUE
    if (i < length(x) && !startsWith(x[[i + 1L]], "--")) {
      val <- x[[i + 1L]]
      i <- i + 1L
    }
    out[[nm]] <- val
    i <- i + 1L
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  if (is.logical(x)) {
    return(isTRUE(x))
  }
  tolower(as.character(x)) %in% c("true", "t", "yes", "y", "1")
}

now_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

load_local_flexybayes <- function(pkg_dir = "flexyBayes") {
  # Load the local source tree without installing. devtools is already
  # a development dependency in this workspace; pkgload is the fallback.
  if (requireNamespace("devtools", quietly = TRUE)) {
    suppressMessages(devtools::load_all(pkg_dir, quiet = TRUE))
  } else if (requireNamespace("pkgload", quietly = TRUE)) {
    suppressMessages(pkgload::load_all(pkg_dir, quiet = TRUE))
  } else {
    stop(
      "Need devtools or pkgload to load local flexyBayes source.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

pkg_version_from_description <- function(pkg_dir = "flexyBayes") {
  dcf <- read.dcf(file.path(pkg_dir, "DESCRIPTION"))
  list(package = dcf[1, "Package"], version = dcf[1, "Version"])
}

installed_table <- function() {
  ip <- installed.packages()[, c("Package", "Version"), drop = FALSE]
  pkgs <- c(
    "devtools",
    "pkgload",
    "data.table",
    "testthat",
    "greta",
    "INLA",
    "brms",
    "posterior",
    "lme4",
    "mgcv",
    "coda",
    "Matrix",
    "sn",
    "vdiffr"
  )
  dt(
    package = pkgs,
    installed = pkgs %in% ip[, "Package"],
    version = ifelse(
      pkgs %in% ip[, "Package"],
      ip[match(pkgs, ip[, "Package"]), "Version"],
      NA_character_
    )
  )
}

make_gaussian_ri <- function(n = 80L, n_group = 8L, seed = 1L) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_group), length.out = n))
  x <- stats::rnorm(n)
  u <- stats::rnorm(n_group, 0, 0.7)[as.integer(g)]
  y <- 2 + 0.6 * x + u + stats::rnorm(n, 0, 0.6)
  data.frame(y = y, x = x, g = g)
}

make_poisson_ri <- function(n = 100L, n_group = 10L, seed = 2L) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_group), length.out = n))
  x <- stats::rnorm(n)
  eta <- 0.2 + 0.35 * x + stats::rnorm(n_group, 0, 0.35)[as.integer(g)]
  y <- stats::rpois(n, lambda = exp(eta))
  data.frame(y = y, x = x, g = g)
}

make_structured_data <- function(seed = 3L) {
  set.seed(seed)
  n_geno <- 12L
  n_env <- 3L
  dat <- expand.grid(
    geno = factor(seq_len(n_geno)),
    env = factor(seq_len(n_env)),
    rep = seq_len(2L)
  )
  dat$x <- stats::rnorm(nrow(dat))
  dat$y <- 10 +
    as.numeric(dat$env) +
    0.5 * dat$x +
    stats::rnorm(nrow(dat), 0, 1)
  A <- crossprod(matrix(stats::rnorm(n_geno * n_geno), n_geno, n_geno))
  diag(A) <- diag(A) + 0.1
  list(dat = dat, Gmat = A)
}

safe_eval <- function(expr) {
  start <- proc.time()[["elapsed"]]
  out <- tryCatch(
    {
      value <- force(expr)
      list(status = "ok", value = value, error = NA_character_)
    },
    error = function(e) {
      list(status = "error", value = NULL, error = conditionMessage(e))
    },
    warning = function(w) {
      invokeRestart("muffleWarning")
    }
  )
  out$elapsed_sec <- proc.time()[["elapsed"]] - start
  out
}

bench_expr <- function(id, group, label, reps, expr) {
  expr_sub <- substitute(expr)
  expr_env <- parent.frame()
  gc()
  start <- proc.time()[["elapsed"]]
  err <- NA_character_
  status <- "ok"
  for (i in seq_len(reps)) {
    ok <- tryCatch(
      {
        eval(expr_sub, envir = expr_env)
        TRUE
      },
      error = function(e) {
        err <<- conditionMessage(e)
        FALSE
      }
    )
    if (!ok) {
      status <- "error"
      break
    }
  }
  elapsed <- proc.time()[["elapsed"]] - start
  dt(
    task_id = id,
    group = group,
    label = label,
    reps = reps,
    status = status,
    elapsed_sec = elapsed,
    sec_per_rep = elapsed / max(1L, reps),
    error = err
  )
}

run_microbenchmarks <- function(
  pkg_dir = "flexyBayes",
  reps = 100L,
  cores = 4L
) {
  load_local_flexybayes(pkg_dir)

  fb_from_brms <- getFromNamespace("fb_from_brms", "flexyBayes")
  fb_from_asreml <- getFromNamespace("fb_from_asreml", "flexyBayes")
  lgm_gate <- getFromNamespace("lgm_gate", "flexyBayes")
  priors_to_inla <- getFromNamespace("priors_to_inla", "flexyBayes")
  priors_to_legacy <- getFromNamespace("priors_to_legacy", "flexyBayes")
  priors_to_brms_specs <- getFromNamespace(
    ".priors_to_brms_specs",
    "flexyBayes"
  )

  d_g <- make_gaussian_ri(80L, 8L, 11L)
  d_p <- make_poisson_ri(100L, 10L, 12L)
  sdat <- make_structured_data(13L)

  fb_ri <- fb_from_brms(y ~ x + (1 | g), data = d_g)
  p <- fb_prior(
    sigma ~ uniform(lower = 0, upper = 5),
    sd(group = "g") ~ pc(upper = 2, prob = 0.05),
    b("x") ~ normal(mean = 0, sd = 10)
  )

  jobs <- list(
    function() {
      bench_expr(
        "parse_brms_fixed",
        "parse_ir",
        "fb_from_brms fixed-only gaussian",
        reps,
        fb_from_brms(y ~ x, data = d_g)
      )
    },
    function() {
      bench_expr(
        "parse_brms_ri",
        "parse_ir",
        "fb_from_brms gaussian random intercept",
        reps,
        fb_from_brms(y ~ x + (1 | g), data = d_g)
      )
    },
    function() {
      bench_expr(
        "parse_brms_crossed",
        "parse_ir",
        "fb_from_brms crossed random intercepts",
        reps,
        {
          dd <- d_g
          dd$h <- factor(rep(1:4, length.out = nrow(dd)))
          fb_from_brms(y ~ x + (1 | g) + (1 | h), data = dd)
        }
      )
    },
    function() {
      bench_expr(
        "parse_asreml_ri",
        "parse_ir",
        "fb_from_asreml random intercept",
        reps,
        fb_from_asreml(y ~ x, random = ~g, data = d_g)
      )
    },
    function() {
      bench_expr(
        "parse_asreml_vm_code_ir",
        "parse_ir",
        "fb_from_asreml vm(geno, Gmat)",
        max(10L, reps %/% 5L),
        fb_from_asreml(
          y ~ env + x,
          random = ~ vm(geno, Gmat),
          data = sdat$dat,
          known_matrices = list(Gmat = sdat$Gmat)
        )
      )
    },
    function() {
      bench_expr(
        "lgm_gate_accept",
        "gate",
        "lgm_gate accepted gaussian RI",
        reps,
        lgm_gate(fb_ri)
      )
    },
    function() {
      bench_expr(
        "lgm_gate_refuse_us",
        "gate",
        "lgm_gate structured covariance refusal",
        max(10L, reps %/% 5L),
        {
          fb_us <- fb_from_asreml(
            y ~ env + x,
            random = ~ us(env):id(geno),
            data = sdat$dat
          )
          lgm_gate(fb_us)
        }
      )
    },
    function() {
      bench_expr(
        "prior_parse",
        "priors",
        "fb_prior DSL parse",
        reps,
        fb_prior(
          sigma ~ uniform(0, 5),
          sd(group = "g") ~ pc(2, 0.05),
          b("x") ~ normal(0, 10)
        )
      )
    },
    function() {
      bench_expr(
        "prior_to_inla",
        "priors",
        "priors_to_inla translation",
        reps,
        priors_to_inla(p)
      )
    },
    function() {
      bench_expr(
        "prior_to_legacy",
        "priors",
        "priors_to_legacy translation",
        reps,
        priors_to_legacy(p, fixed_sd_default = 100)
      )
    },
    function() {
      bench_expr(
        "prior_to_brms_specs",
        "priors",
        "priors_to_brms_specs translation",
        reps,
        priors_to_brms_specs(p, fb_ri)
      )
    }
  )

  # Fork only pure-R / parse-level work here. Heavy backend fits are
  # handled as isolated child Rscript processes below.
  ans <- parallel::mclapply(jobs, function(f) f(), mc.cores = cores)
  data.table::rbindlist(ans, fill = TRUE)
}

# --------------------------------------------------------------------- #
# Representation benchmarks (the design spec)              #
# --------------------------------------------------------------------- #
#
# Three measurements per (term_class, n) cell:
#
#   generated_code_size_bytes   nchar(return_code(fit)) at representative n
#   pre_fit_object_size_bytes   object.size() of the largest R object
#                               the codegen path creates BEFORE
#                               greta::mcmc() is called -- the design
#                               matrix or, on the indexed-emit path,
#                               its indexed equivalent
#   code_generation_time_sec    wall time for flexybayes(..., return_code = TRUE)
#
# The point of these rows is the n-scaling check: under the indexed
# emit path, code size stays nearly constant in n because the design
# matrix is bound into the model environment by reference rather than
# inlined as a literal. The "dense baseline" cells are not measured
# here because no dense-baseline path is exposed at the R surface in
# the indexed-codegen package; the n-scaling check below is the
# operational version of the analytic memory comparison.
#
# Term classes covered (post-Stage-0):
#   fixed_factor_simple    y ~ g          (g a 5-level factor)
#   random_intercept       y ~ x + (1|g)
#   smooth                 y ~ s(x)
#
# Two additional classes scheduled for v0.2.6 (Stage 1A + 1B) are not
# exercised here yet -- they are added once those merges land.

make_repr_data <- function(n, n_group = 8L, seed = 71L) {
  set.seed(seed)
  g <- factor(rep(seq_len(n_group), length.out = n))
  x <- stats::rnorm(n)
  y <- 1 +
    0.4 * x +
    stats::rnorm(n_group, 0, 0.3)[as.integer(g)] +
    stats::rnorm(n, 0, 0.5)
  data.frame(y = y, x = x, g = g)
}

# Estimate the largest R object the codegen path creates BEFORE
# greta::mcmc() runs. For the supported indexed-emit term classes this
# is either the integer level-index vector (random intercept, simple
# fixed factor) or the basis matrix (smooth s(x)). The estimate is
# constructed deterministically from data + the IR shape, matching what
# the emit path actually allocates -- not what model.matrix() would
# allocate on a hypothetical dense fall-back.
.repr_pre_fit_object_size <- function(term_class, data, fb = NULL) {
  switch(
    term_class,
    fixed_factor_simple = {
      # Indexed emit binds the factor level vector (integer per row).
      utils::object.size(as.integer(data$g))
    },
    random_intercept = {
      # Indexed emit binds the level-index vector for the random group;
      # the per-row matrix is never materialised.
      utils::object.size(as.integer(data$g))
    },
    smooth = {
      # Smooth emit binds the basis matrix B_s_x; the largest pre-fit
      # object is therefore that matrix, not a model.matrix() expansion.
      if (
        !is.null(fb) &&
          length(fb$random_terms) >= 1L &&
          !is.null(fb$random_terms[[1L]]$smooth_obj)
      ) {
        smobj <- fb$random_terms[[1L]]$smooth_obj
        if (!is.null(smobj$X)) {
          return(utils::object.size(smobj$X))
        }
      }
      # Defensive fallback: synthesise the basis using mgcv directly.
      if (requireNamespace("mgcv", quietly = TRUE)) {
        sm <- mgcv::smoothCon(
          mgcv::s(x, k = 6L),
          data = data,
          knots = NULL,
          absorb.cons = TRUE
        )[[1L]]
        utils::object.size(sm$X)
      } else {
        NA_real_
      }
    },
    NA_real_
  )
}

.repr_one_cell <- function(term_class, n, pkg_dir) {
  d <- make_repr_data(n)
  formula_pieces <- switch(
    term_class,
    fixed_factor_simple = list(fixed = y ~ g, random = NULL),
    random_intercept = list(fixed = y ~ x, random = ~g),
    smooth = list(fixed = y ~ 1, random = ~ s(x, k = 6L)),
    stop("unknown term_class: ", term_class, call. = FALSE)
  )

  # Code-generation wall time: return_code = TRUE never enters mcmc().
  t0 <- proc.time()[["elapsed"]]
  code <- tryCatch(
    suppressMessages(flexybayes(
      fixed = formula_pieces$fixed,
      random = formula_pieces$random,
      data = d,
      return_code = TRUE,
      verbose = FALSE,
      mcmc_verbose = FALSE
    )),
    error = function(e) structure(NA_character_, error = conditionMessage(e))
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  if (is.na(code) || is.null(code) || !is.character(code)) {
    return(dt(
      task_id = paste("repr", term_class, n, sep = "_"),
      group = "representation",
      term_class = term_class,
      n = n,
      status = "error",
      generated_code_size_bytes = NA_real_,
      pre_fit_object_size_bytes = NA_real_,
      code_generation_time_sec = elapsed,
      error = attr(code, "error") %||% "codegen returned non-character"
    ))
  }

  # Re-derive the IR off the same call so the smooth Smooth object is
  # available for the basis-size estimate without paying the codegen
  # cost again. Pulled via getFromNamespace() so this function does
  # not depend on the symbol being assigned into its enclosing env.
  fb_from_asreml_ns <- getFromNamespace("fb_from_asreml", "flexyBayes")
  fb_for_size <- tryCatch(
    suppressMessages(fb_from_asreml_ns(
      fixed = formula_pieces$fixed,
      random = formula_pieces$random,
      data = d
    )),
    error = function(e) NULL
  )

  dt(
    task_id = paste("repr", term_class, n, sep = "_"),
    group = "representation",
    term_class = term_class,
    n = n,
    status = "ok",
    generated_code_size_bytes = as.numeric(nchar(code, type = "bytes")),
    pre_fit_object_size_bytes = as.numeric(
      .repr_pre_fit_object_size(term_class, d, fb_for_size)
    ),
    code_generation_time_sec = elapsed,
    error = NA_character_
  )
}

run_representation_benchmarks <- function(
  pkg_dir = "flexyBayes",
  n_grid = c(50L, 500L, 5000L),
  term_classes = c("fixed_factor_simple", "random_intercept", "smooth"),
  cores = 2L
) {
  # Loaded greta is required only because the greta-emit codegen path
  # itself touches greta's namespace at code-build time. We skip the
  # representation block cleanly if greta is not installed rather than
  # erroring -- the rest of the harness should still run.
  if (!requireNamespace("greta", quietly = TRUE)) {
    return(dt(
      task_id = "representation_skipped_no_greta",
      group = "representation",
      term_class = NA_character_,
      n = NA_integer_,
      status = "skipped",
      generated_code_size_bytes = NA_real_,
      pre_fit_object_size_bytes = NA_real_,
      code_generation_time_sec = NA_real_,
      error = "greta not installed; representation block skipped"
    ))
  }

  load_local_flexybayes(pkg_dir)

  cells <- expand.grid(
    term_class = term_classes,
    n = n_grid,
    stringsAsFactors = FALSE
  )

  worker <- function(i) {
    .repr_one_cell(cells$term_class[i], cells$n[i], pkg_dir)
  }

  ans <- parallel::mclapply(
    seq_len(nrow(cells)),
    worker,
    mc.cores = max(1L, cores)
  )
  data.table::rbindlist(ans, fill = TRUE)
}

child_result <- function(
  task_id,
  group,
  status,
  elapsed_sec,
  details = list(),
  error = NA_character_
) {
  structure(
    list(
      task_id = task_id,
      group = group,
      status = status,
      elapsed_sec = elapsed_sec,
      error = error,
      details = details
    ),
    class = "fb_bench_child_result"
  )
}

run_child_task <- function(task_id, pkg_dir = "flexyBayes") {
  t0 <- proc.time()[["elapsed"]]
  load_local_flexybayes(pkg_dir)

  finish <- function(
    group,
    status = "ok",
    details = list(),
    error = NA_character_
  ) {
    child_result(
      task_id,
      group,
      status,
      proc.time()[["elapsed"]] - t0,
      details,
      error
    )
  }

  tryCatch(
    {
      if (identical(task_id, "review_code_greta")) {
        d <- make_gaussian_ri(60L, 6L, 101L)
        rev <- suppressMessages(flexybayes(
          fixed = y ~ x,
          random = ~g,
          data = d,
          review_code = TRUE,
          verbose = FALSE,
          mcmc_verbose = FALSE
        ))
        code <- capture.output(cat_code(rev))
        return(finish(
          "code_review",
          details = list(
            class = paste(class(rev), collapse = "/"),
            code_lines = length(code),
            code_chars = sum(nchar(code))
          )
        ))
      }

      if (identical(task_id, "structured_cov_return_code")) {
        sdat <- make_structured_data(102L)
        code <- suppressMessages(flexybayes(
          fixed = y ~ env + x,
          random = ~ vm(geno, Gmat) + at(env):geno,
          data = sdat$dat,
          known_matrices = list(Gmat = sdat$Gmat),
          return_code = TRUE,
          verbose = FALSE,
          mcmc_verbose = FALSE
        ))
        return(finish(
          "code_generation",
          details = list(
            code_lines = length(strsplit(code, "\n", fixed = TRUE)[[1L]]),
            code_chars = nchar(code),
            has_vm = grepl("L_G_geno", code, fixed = TRUE),
            has_at = grepl("sigma_geno_env", code, fixed = TRUE)
          )
        ))
      }

      if (identical(task_id, "inla_gaussian_ri")) {
        d <- make_gaussian_ri(100L, 10L, 201L)
        fit <- suppressMessages(fb_brms(
          y ~ x + (1 | g),
          data = d,
          backend = "inla",
          verbose = FALSE
        ))
        cn <- canonical_names(fit)
        bd <- backend_decision(fit)
        return(finish(
          "backend_fit",
          details = list(
            backend = bd$backend,
            path = bd$path,
            n_fixed = nrow(fit$inla$summary.fixed),
            n_hyper = nrow(fit$inla$summary.hyperpar),
            canonical_n = length(cn$map),
            numerical_confirm = isTRUE(fit$num_check$pass)
          )
        ))
      }

      if (identical(task_id, "inla_poisson_ri")) {
        d <- make_poisson_ri(120L, 10L, 202L)
        fit <- suppressMessages(fb_brms(
          y ~ x + (1 | g),
          data = d,
          family = "poisson",
          backend = "inla",
          verbose = FALSE
        ))
        bd <- backend_decision(fit)
        return(finish(
          "backend_fit",
          details = list(
            backend = bd$backend,
            path = bd$path,
            n_fixed = nrow(fit$inla$summary.fixed),
            n_hyper = nrow(fit$inla$summary.hyperpar),
            numerical_confirm = isTRUE(fit$num_check$pass)
          )
        ))
      }

      if (identical(task_id, "greta_gaussian_fixed")) {
        d <- make_gaussian_ri(50L, 5L, 301L)
        fit <- suppressMessages(fb_brms(
          y ~ x,
          data = d,
          backend = "greta",
          n_samples = 60L,
          warmup = 60L,
          chains = 1L,
          verbose = FALSE,
          mcmc_verbose = FALSE
        ))
        bd <- backend_decision(fit)
        return(finish(
          "backend_fit",
          details = list(
            backend = bd$backend,
            path = bd$path,
            n_params = fit$extras$model_info$n_params,
            max_rhat = suppressWarnings(max(
              fit$extras$convergence$gelman$psrf[, "Point est."],
              na.rm = TRUE
            ))
          )
        ))
      }

      if (identical(task_id, "greta_gaussian_ri")) {
        d <- make_gaussian_ri(60L, 6L, 302L)
        fit <- suppressMessages(fb_brms(
          y ~ x + (1 | g),
          data = d,
          backend = "greta",
          n_samples = 60L,
          warmup = 60L,
          chains = 1L,
          verbose = FALSE,
          mcmc_verbose = FALSE
        ))
        bd <- backend_decision(fit)
        ps <- prior_summary(fit)
        return(finish(
          "backend_fit",
          details = list(
            backend = bd$backend,
            path = bd$path,
            n_params = fit$extras$model_info$n_params,
            prior_kind = ps$kind,
            prior_origin = ps$default_origin
          )
        ))
      }

      if (identical(task_id, "brms_stancode")) {
        d <- make_gaussian_ri(60L, 6L, 401L)
        rev <- suppressMessages(fb_brms(
          y ~ x + (1 | g),
          data = d,
          backend = "brms",
          review_code = TRUE,
          verbose = FALSE,
          mcmc_verbose = FALSE
        ))
        return(finish(
          "stan_passthrough",
          details = list(
            class = paste(class(rev), collapse = "/"),
            backend = rev$backend,
            code_chars = nchar(rev$code),
            has_parameters_block = grepl("parameters", rev$code, fixed = TRUE)
          )
        ))
      }

      if (identical(task_id, "brms_gaussian_fit")) {
        d <- make_gaussian_ri(50L, 5L, 402L)
        fit <- suppressMessages(suppressWarnings(fb_brms(
          y ~ x + (1 | g),
          data = d,
          backend = "brms",
          n_samples = 80L,
          warmup = 80L,
          chains = 1L,
          verbose = FALSE,
          mcmc_verbose = FALSE
        )))
        bd <- backend_decision(fit)
        return(finish(
          "backend_fit",
          details = list(
            backend = bd$backend,
            path = bd$path,
            n_params = fit$extras$model_info$n_params,
            coef_names = paste(names(coef(fit)), collapse = ",")
          )
        ))
      }

      if (identical(task_id, "triangulate_greta_inla")) {
        d <- make_gaussian_ri(60L, 6L, 501L)
        fit_g <- suppressMessages(fb_brms(
          y ~ x + (1 | g),
          data = d,
          backend = "greta",
          n_samples = 60L,
          warmup = 60L,
          chains = 1L,
          verbose = FALSE,
          mcmc_verbose = FALSE
        ))
        fit_i <- suppressMessages(fb_brms(
          y ~ x + (1 | g),
          data = d,
          backend = "inla",
          verbose = FALSE
        ))
        tri <- triangulate(fit_g, fit_i, n_samples = 200L)
        return(finish(
          "triangulation",
          details = list(
            source_a = tri$source_a,
            source_b = tri$source_b,
            n_common = tri$n_common,
            common = paste(tri$common, collapse = ","),
            max_wasserstein = if (nrow(tri$metrics)) {
              max(tri$metrics$wasserstein_1, na.rm = TRUE)
            } else {
              NA_real_
            }
          )
        ))
      }

      if (identical(task_id, "smooth_greta_fit")) {
        d <- make_gaussian_ri(60L, 6L, 601L)
        fit <- suppressMessages(flexybayes(
          fixed = y ~ 1,
          random = ~ s(x),
          data = d,
          backend = "greta",
          n_samples = 60L,
          warmup = 60L,
          chains = 1L,
          verbose = FALSE,
          mcmc_verbose = FALSE
        ))
        pred <- predict(fit)
        return(finish(
          "smooths",
          details = list(
            n_pred = length(pred),
            pred_has_na = anyNA(pred),
            n_params = fit$extras$model_info$n_params
          )
        ))
      }

      if (identical(task_id, "testthat_tally")) {
        # This runs the project's own tally script. It can be expensive
        # because NOT_CRAN=true lets selected greta/INLA integration
        # tests run. The parent process enforces the timeout.
        Sys.setenv(NOT_CRAN = "true")
        out <- utils::capture.output(source("tools/tally.R", local = TRUE))
        line <- out[grepl("PASS:", out, fixed = TRUE)]
        return(finish(
          "test_suite",
          details = list(
            tally_line = paste(line, collapse = " | "),
            output_lines = length(out)
          )
        ))
      }

      finish(
        "unknown",
        status = "error",
        error = paste("Unknown child task:", task_id)
      )
    },
    error = function(e) {
      finish("runtime", status = "error", error = conditionMessage(e))
    }
  )
}

child_to_dt <- function(x) {
  det <- x$details %||% list()
  detail_json <- paste(
    sprintf(
      "%s=%s",
      names(det),
      vapply(
        det,
        function(v) {
          if (length(v) == 0L) {
            return("")
          }
          paste(as.character(v), collapse = ";")
        },
        character(1)
      )
    ),
    collapse = " | "
  )
  dt(
    task_id = x$task_id,
    group = x$group,
    status = x$status,
    elapsed_sec = x$elapsed_sec,
    error = x$error,
    details = detail_json
  )
}

run_child_as_process <- function(
  script,
  task_id,
  out_dir,
  timeout_sec,
  pkg_dir = "flexyBayes"
) {
  child_rds <- file.path(out_dir, paste0("child_", task_id, ".rds"))
  child_log <- file.path(out_dir, paste0("child_", task_id, ".log"))
  args <- c(
    "--vanilla",
    normalizePath(script),
    "--child",
    task_id,
    "--out",
    child_rds,
    "--pkg-dir",
    pkg_dir
  )

  t0 <- proc.time()[["elapsed"]]
  warn_msg <- NA_character_
  status <- tryCatch(
    withCallingHandlers(
      system2(
        file.path(R.home("bin"), "Rscript"),
        args = args,
        stdout = child_log,
        stderr = child_log,
        timeout = timeout_sec
      ),
      warning = function(w) {
        warn_msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      structure(125L, error_message = conditionMessage(e))
    }
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  if (file.exists(child_rds)) {
    res <- readRDS(child_rds)
    return(child_to_dt(res))
  }

  log_tail <- if (file.exists(child_log)) {
    paste(
      utils::tail(readLines(child_log, warn = FALSE), 20L),
      collapse = " / "
    )
  } else {
    ""
  }

  dt(
    task_id = task_id,
    group = "child_process",
    status = if (identical(as.integer(status), 124L)) "timeout" else "error",
    elapsed_sec = elapsed,
    error = paste0(
      "child did not produce RDS; system status = ",
      paste(status, collapse = ","),
      if (!is.na(warn_msg)) paste0("; warning: ", warn_msg) else "",
      if (nzchar(log_tail)) paste0("; log tail: ", log_tail) else ""
    ),
    details = paste0("timeout_sec=", timeout_sec)
  )
}

task_plan <- function(profile = "standard", include_brms_fit = FALSE) {
  # timeout_sec is intentionally generous enough to let honest backend
  # setup happen, while short enough to keep the whole run under the
  # requested 3 hour budget.
  base <- dt(
    task_id = c(
      "review_code_greta",
      "structured_cov_return_code",
      "inla_gaussian_ri",
      "inla_poisson_ri",
      "greta_gaussian_fixed",
      "greta_gaussian_ri",
      "brms_stancode",
      "triangulate_greta_inla",
      "smooth_greta_fit"
    ),
    tier = c(
      "light",
      "light",
      "fit",
      "fit",
      "fit",
      "fit",
      "light",
      "fit",
      "fit"
    ),
    timeout_sec = c(300, 300, 600, 600, 1200, 1200, 900, 1800, 1200)
  )

  if (identical(profile, "smoke")) {
    base <- base[
      task_id %chin%
        c(
          "review_code_greta",
          "structured_cov_return_code",
          "inla_gaussian_ri",
          "brms_stancode"
        )
    ]
  }

  if (identical(profile, "full")) {
    base <- rbind(
      base,
      dt(task_id = "testthat_tally", tier = "fit", timeout_sec = 3600),
      fill = TRUE
    )
  }

  if (isTRUE(include_brms_fit)) {
    base <- rbind(
      base,
      dt(task_id = "brms_gaussian_fit", tier = "fit", timeout_sec = 2400),
      fill = TRUE
    )
  }

  base[]
}

run_child_batch <- function(
  plan,
  script,
  out_dir,
  cores,
  pkg_dir,
  budget_deadline
) {
  if (!nrow(plan)) {
    return(dt())
  }

  worker <- function(row_i) {
    row <- plan[row_i]
    remaining <- as.numeric(difftime(
      budget_deadline,
      Sys.time(),
      units = "secs"
    ))
    if (!is.finite(remaining) || remaining <= 60) {
      return(dt(
        task_id = row$task_id,
        group = "budget",
        status = "skipped_budget_exhausted",
        elapsed_sec = 0,
        error = "Global benchmark budget exhausted before task start.",
        details = ""
      ))
    }
    timeout <- min(row$timeout_sec, max(60, floor(remaining - 30)))
    run_child_as_process(
      script,
      row$task_id,
      out_dir,
      timeout,
      pkg_dir = pkg_dir
    )
  }

  ans <- parallel::mclapply(seq_len(nrow(plan)), worker, mc.cores = cores)
  data.table::rbindlist(ans, fill = TRUE)
}

write_report <- function(
  out_dir,
  env_dt,
  micro_dt,
  child_dt,
  summary_dt,
  config,
  started,
  finished,
  repr_dt = NULL
) {
  pkg <- pkg_version_from_description(config$pkg_dir)
  ok_child <- if (nrow(child_dt)) sum(child_dt$status == "ok") else 0L
  fail_child <- if (nrow(child_dt)) sum(child_dt$status != "ok") else 0L
  ok_micro <- if (nrow(micro_dt)) sum(micro_dt$status == "ok") else 0L
  fail_micro <- if (nrow(micro_dt)) sum(micro_dt$status != "ok") else 0L

  lines <- c(
    "# flexyBayes v0.2.x Local Assessment / Benchmark Report",
    "",
    paste0("Started: ", format(started, "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Finished: ", format(finished, "%Y-%m-%d %H:%M:%S %Z")),
    paste0(
      "Elapsed minutes: ",
      round(as.numeric(difftime(finished, started, units = "mins")), 2)
    ),
    "",
    "## Configuration",
    "",
    paste0("- Local package: `", pkg$package, "`"),
    paste0("- DESCRIPTION version: `", pkg$version, "`"),
    paste0("- Requested profile: `", config$profile, "`"),
    paste0("- Max cores requested/used: `", config$max_cores, "`"),
    paste0("- Heavy backend parallelism: `", config$fit_cores, "`"),
    paste0("- Budget minutes: `", config$budget_min, "`"),
    paste0("- Include brms sampling fit: `", config$include_brms_fit, "`"),
    "",
    "Note: this report benchmarks the current local source tree; the",
    "DESCRIPTION version is reported above without rewriting it.",
    "",
    "## Headline Result",
    "",
    paste0("- Microbench tasks OK/error: ", ok_micro, " / ", fail_micro),
    paste0("- Child/backend tasks OK/non-OK: ", ok_child, " / ", fail_child),
    "",
    "## Summary By Group",
    "",
    if (nrow(summary_dt)) {
      c(
        "| group | tasks | ok | non_ok | elapsed_sec_median | elapsed_sec_max |",
        "|---|---:|---:|---:|---:|---:|",
        apply(summary_dt, 1L, function(r) {
          paste0(
            "| ",
            r[["group"]],
            " | ",
            r[["tasks"]],
            " | ",
            r[["ok"]],
            " | ",
            r[["non_ok"]],
            " | ",
            r[["elapsed_sec_median"]],
            " | ",
            r[["elapsed_sec_max"]],
            " |"
          )
        })
      )
    } else {
      "_No summary rows._"
    },
    "",
    "## Non-OK Child Tasks",
    "",
    if (nrow(child_dt[status != "ok"])) {
      c(
        "| task_id | status | elapsed_sec | error |",
        "|---|---|---:|---|",
        apply(child_dt[status != "ok"], 1L, function(r) {
          err <- gsub("\\|", "/", r[["error"]] %||% "")
          paste0(
            "| ",
            r[["task_id"]],
            " | ",
            r[["status"]],
            " | ",
            round(as.numeric(r[["elapsed_sec"]]), 2),
            " | ",
            err,
            " |"
          )
        })
      )
    } else {
      "_None._"
    },
    "",
    "## Slowest Child Tasks",
    "",
    if (nrow(child_dt)) {
      top <- child_dt[order(-elapsed_sec)][seq_len(min(.N, 10L))]
      c(
        "| task_id | group | status | elapsed_sec | details |",
        "|---|---|---|---:|---|",
        apply(top, 1L, function(r) {
          det <- gsub("\\|", "/", r[["details"]] %||% "")
          paste0(
            "| ",
            r[["task_id"]],
            " | ",
            r[["group"]],
            " | ",
            r[["status"]],
            " | ",
            round(as.numeric(r[["elapsed_sec"]]), 2),
            " | ",
            det,
            " |"
          )
        })
      )
    } else {
      "_No child tasks were run._"
    },
    "",
    "## Microbenchmark Notes",
    "",
    "Microbenchmarks are repeated parse/translation operations, so the",
    "`sec_per_rep` column in `microbench_results.csv` is the useful number.",
    "Backend fit tasks are deliberately separated because they include",
    "runtime setup, compilation, or sampler initialization overhead.",
    "",
    "## Representation Benchmarks (the design spec)",
    "",
    if (!is.null(repr_dt) && nrow(repr_dt)) {
      c(
        "Three measurements per (term_class, n) cell:",
        "`generated_code_size_bytes`, `pre_fit_object_size_bytes`,",
        "`code_generation_time_sec`. Indexed-emit paths are expected to",
        "keep generated code size nearly constant in `n`; the only",
        "n-scaling object the codegen path materialises is the per-row",
        "level-index vector (random intercept, simple fixed factor) or",
        "the basis matrix (smooth `s(x)`).",
        "",
        "| task_id | term_class | n | status | code_bytes | pre_fit_bytes | codegen_sec |",
        "|---|---|---:|---|---:|---:|---:|",
        apply(repr_dt, 1L, function(r) {
          paste0(
            "| ",
            r[["task_id"]],
            " | ",
            r[["term_class"]],
            " | ",
            r[["n"]],
            " | ",
            r[["status"]],
            " | ",
            r[["generated_code_size_bytes"]],
            " | ",
            r[["pre_fit_object_size_bytes"]],
            " | ",
            round(as.numeric(r[["code_generation_time_sec"]]), 4),
            " |"
          )
        })
      )
    } else {
      "_Representation block produced no rows (greta may not be installed)._"
    },
    "",
    "## Output Files",
    "",
    "- `environment.csv`",
    "- `microbench_results.csv`",
    "- `representation_results.csv`",
    "- `child_task_results.csv`",
    "- `summary_by_group.csv`",
    "- `raw_results.rds`",
    "- child logs: `child_*.log`"
  )

  writeLines(lines, file.path(out_dir, "benchmark_report.md"))
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  # Child mode: run exactly one task and save a compact RDS result.
  if (!is.null(args$child)) {
    out <- args$out %||%
      stop("--out is required in --child mode", call. = FALSE)
    pkg_dir <- args[["pkg-dir"]] %||% "flexyBayes"
    res <- run_child_task(args$child, pkg_dir = pkg_dir)
    saveRDS(res, out)
    return(invisible(NULL))
  }

  started <- Sys.time()
  pkg_dir <- args[["pkg-dir"]] %||% "flexyBayes"
  profile <- args$profile %||% "standard"
  budget_min <- as.numeric(args[["budget-min"]] %||% 180)
  detected <- parallel::detectCores(logical = TRUE)
  max_cores <- min(
    as.integer(args[["max-cores"]] %||% 8L),
    detected %||% 1L,
    8L
  )
  max_cores <- max(1L, max_cores)
  fit_cores <- min(max_cores, as.integer(args[["fit-cores"]] %||% 4L))
  include_brms_fit <- as_bool(args[["include-brms-fit"]], default = FALSE)
  reps <- as.integer(
    args$reps %||% if (identical(profile, "full")) 150L else 80L
  )

  out_root <- args[["out-root"]] %||% "benchmark_results"
  pkg_meta_early <- pkg_version_from_description(
    args[["pkg-dir"]] %||% "flexyBayes"
  )
  out_dir <- file.path(
    out_root,
    paste0("flexyBayes_", pkg_meta_early$version, "_", now_stamp())
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  config <- list(
    pkg_dir = pkg_dir,
    profile = profile,
    budget_min = budget_min,
    max_cores = max_cores,
    fit_cores = fit_cores,
    include_brms_fit = include_brms_fit,
    reps = reps
  )

  cat("flexyBayes benchmark harness\n")
  cat("  output: ", out_dir, "\n", sep = "")
  cat("  profile: ", profile, "\n", sep = "")
  cat(
    "  cores: ",
    max_cores,
    " total; ",
    fit_cores,
    " heavy-backend workers\n",
    sep = ""
  )
  cat("  budget: ", budget_min, " minutes\n", sep = "")

  env_dt <- installed_table()
  pkg <- pkg_version_from_description(pkg_dir)
  env_dt <- rbind(
    dt(
      package = "LOCAL_SOURCE",
      installed = TRUE,
      version = paste(pkg$package, pkg$version)
    ),
    dt(
      package = "R",
      installed = TRUE,
      version = paste(R.version$major, R.version$minor, sep = ".")
    ),
    dt(
      package = "logical_cores_detected",
      installed = TRUE,
      version = as.character(detected)
    ),
    env_dt,
    fill = TRUE
  )
  data.table::fwrite(env_dt, file.path(out_dir, "environment.csv"))

  budget_deadline <- started + budget_min * 60

  cat("Running parse/gate/prior microbenchmarks...\n")
  micro_dt <- run_microbenchmarks(
    pkg_dir,
    reps = reps,
    cores = min(max_cores, 8L)
  )
  data.table::fwrite(micro_dt, file.path(out_dir, "microbench_results.csv"))

  cat("Running representation benchmarks (the design spec)...\n")
  repr_n_grid <- if (identical(profile, "smoke")) {
    c(50L, 500L)
  } else if (identical(profile, "full")) {
    c(50L, 500L, 5000L, 10000L)
  } else {
    c(50L, 500L, 5000L)
  }
  repr_dt <- run_representation_benchmarks(
    pkg_dir = pkg_dir,
    n_grid = repr_n_grid,
    cores = min(max_cores, 2L)
  )
  data.table::fwrite(repr_dt, file.path(out_dir, "representation_results.csv"))

  plan <- task_plan(profile, include_brms_fit = include_brms_fit)

  # Light child tasks can use the full core cap. Heavy fit tasks use
  # a lower cap to avoid making greta/Stan/INLA compete too hard.
  light_plan <- plan[tier == "light"]
  fit_plan <- plan[tier != "light"]

  ca <- commandArgs(FALSE)
  file_arg <- ca[grep("^--file=", ca)][1]
  script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE)
  if (!file.exists(script_path)) {
    # Fallback for interactive/source() use.
    script_path <- normalizePath("tools/benchmark_flexyBayes.R")
  }

  cat("Running light child tasks...\n")
  if (nrow(light_plan)) {
    light_dt <- run_child_batch(
      light_plan,
      script = script_path,
      out_dir = out_dir,
      cores = max_cores,
      pkg_dir = pkg_dir,
      budget_deadline = budget_deadline
    )
  } else {
    light_dt <- dt()
  }

  cat("Running backend / integration child tasks...\n")
  if (nrow(fit_plan)) {
    fit_dt <- run_child_batch(
      fit_plan,
      script = script_path,
      out_dir = out_dir,
      cores = fit_cores,
      pkg_dir = pkg_dir,
      budget_deadline = budget_deadline
    )
  } else {
    fit_dt <- dt()
  }

  child_dt <- data.table::rbindlist(list(light_dt, fit_dt), fill = TRUE)
  data.table::fwrite(child_dt, file.path(out_dir, "child_task_results.csv"))

  all_elapsed <- data.table::rbindlist(
    list(
      micro_dt[, .(group, status, elapsed_sec)],
      child_dt[, .(group, status, elapsed_sec)]
    ),
    fill = TRUE
  )

  summary_dt <- if (nrow(all_elapsed)) {
    all_elapsed[,
      .(
        tasks = .N,
        ok = sum(status == "ok"),
        non_ok = sum(status != "ok"),
        elapsed_sec_median = round(stats::median(elapsed_sec, na.rm = TRUE), 3),
        elapsed_sec_max = round(max(elapsed_sec, na.rm = TRUE), 3)
      ),
      by = group
    ][order(group)]
  } else {
    dt()
  }
  data.table::fwrite(summary_dt, file.path(out_dir, "summary_by_group.csv"))

  finished <- Sys.time()
  raw <- list(
    config = config,
    environment = env_dt,
    microbench = micro_dt,
    representation = repr_dt,
    child_tasks = child_dt,
    summary_by_group = summary_dt,
    started = started,
    finished = finished
  )
  saveRDS(raw, file.path(out_dir, "raw_results.rds"))

  write_report(
    out_dir,
    env_dt,
    micro_dt,
    child_dt,
    summary_dt,
    config,
    started,
    finished,
    repr_dt = repr_dt
  )

  cat("Benchmark complete.\n")
  cat("  report: ", file.path(out_dir, "benchmark_report.md"), "\n", sep = "")
  cat(
    "  elapsed minutes: ",
    round(as.numeric(difftime(finished, started, units = "mins")), 2),
    "\n",
    sep = ""
  )
  invisible(raw)
}

main()
