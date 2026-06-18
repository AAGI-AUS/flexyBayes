#!/usr/bin/env Rscript
# =============================================================================
# gretaR_worker.R -- the out-of-process gretaR fitting worker.
#
# flexyBayes drives gretaR through this worker because greta (TensorFlow) and
# gretaR (torch) export the same symbols (normal/model/mcmc/distribution/
# as_data) and cannot co-load in one R process. The worker loads ONLY gretaR,
# builds the model via gretaR::model_from_arrays() (re-entrant, deparse-free,
# with out-of-band canonical parameter `names=` that flow straight to the
# draws), samples with the torch NUTS backend, and writes a
# posterior::draws_array back. flexyBayes never loads gretaR in its own session.
#
# Spec (rds at args[1]) -- the lowered flexyBayes IR:
#   list(family, X (named-cols design), y, random = list(Z, J) | NULL,
#        priors = list(fixed_sd, sigma_scale, tau_scale),
#        mcmc = list(n_samples, warmup, chains, seed, n_threads),
#        canonical_names (chr, the X column tokens))
# Output (rds at args[2]):
#   list(status, draws (draws_array, canonical names) | NULL, message, t_total)
#
# gretaR source: GRETAR_HOME env var (load_all the dev tree) or installed pkg.
# Usage: Rscript gretaR_worker.R <spec_rds> <out_rds>
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x
a <- commandArgs(trailingOnly = TRUE)
spec_path <- a[1]
out_path <- a[2]

res <- tryCatch(
  {
    suppressMessages(suppressWarnings({
      home <- Sys.getenv("GRETAR_HOME", "")
      if (nzchar(home)) {
        pkgload::load_all(home, quiet = TRUE)
      } else {
        library(gretaR)
      }
      library(posterior)
    }))
    spec <- readRDS(spec_path)
    nt <- spec$mcmc$n_threads %||% 1L
    if (exists("torch_set_num_threads", where = asNamespace("torch"))) {
      torch::torch_set_num_threads(as.integer(nt))
    }

    reset_gretaR_env()
    t0 <- Sys.time()

    X <- spec$X
    y <- as_data(spec$y)
    p <- ncol(X)
    cn <- spec$canonical_names # X column tokens, e.g. (Intercept), x
    # fixed-effect coefficient vector
    b <- normal(0, spec$priors$fixed_sd, dim = p)
    eta <- as_data(X) %*% b

    targets <- list(b)
    names_list <- list(cn)

    if (!is.null(spec$random)) {
      # non-centred random intercept
      J <- spec$random$J
      tau <- half_normal(spec$priors$tau_scale)
      u_raw <- normal(0, 1, dim = J)
      u <- u_raw * tau
      eta <- eta + as_data(spec$random$Z) %*% u
      targets <- c(targets, list(tau))
      names_list <- c(names_list, list("tau"))
    }

    fam <- spec$family
    if (fam == "gaussian") {
      sigma <- half_normal(spec$priors$sigma_scale)
      distribution(y) <- normal(eta, sigma)
      targets <- c(targets, list(sigma))
      names_list <- c(names_list, list("sigma"))
    } else if (fam == "binomial") {
      distribution(y) <- bernoulli(1 / (1 + exp(-eta)))
    } else if (fam == "poisson") {
      distribution(y) <- poisson_dist(exp(eta))
    } else {
      stop("gretaR_worker: unsupported family '", fam, "'.")
    }

    m <- model_from_arrays(
      targets = targets,
      likelihood = y,
      names = names_list
    )
    fit <- mcmc(
      m,
      n_samples = spec$mcmc$n_samples,
      warmup = spec$mcmc$warmup,
      chains = spec$mcmc$chains,
      sampler = "nuts",
      backend = "torch",
      seed = spec$mcmc$seed,
      verbose = FALSE
    )
    t_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    list(
      status = "ok",
      draws = fit$draws,
      t_total = t_total,
      gretaR_version = as.character(utils::packageVersion("gretaR"))
    )
  },
  error = function(e) list(status = "error", message = conditionMessage(e))
)

saveRDS(res, out_path)
