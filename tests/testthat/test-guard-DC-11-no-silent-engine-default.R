# Class guard DC-11: an engine default taken silently.
#
# 1. Every INLA hyperparameter a fit estimates is listed in its prior
#    record, prior_summary(fit)$slots, one slot each. A standard deviation of a term the
#    default-prior walker covers carries the package's prior, never INLA's
#    own. Only the smoothing precision of spl() is left to the engine, and
#    the record says so. Swept over the term kinds of the class's register
#    rows: nested, diag, slope and AR1 field terms (DEF-047, DEF-129,
#    DEF-189), the residual, a sectioned residual and a multi-way term.
# 2. Every fit records the thread setting it ran with (DEF-019).
# 3. In the package source, every hyperparameter draw passes
#    improve.marginals = TRUE (DEF-002), and every in-package joint draw
#    sets num.threads (DEF-019).

.g11_data <- function() {
  set.seed(11)
  d <- expand.grid(
    rowf = factor(1:6), colf = factor(1:4), env = factor(c("e1", "e2", "e3"))
  )
  d$gen <- factor(sample(sprintf("v%d", 1:8), nrow(d), replace = TRUE))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 2 + as.integer(d$env) + stats::rnorm(8L)[d$gen] + stats::rnorm(nrow(d))
  d$n <- stats::rpois(nrow(d), 3)
  d
}

.g11_corpus <- list(
  iid = list(y ~ env, random = ~gen),
  poisson_aggregated = list(n ~ env, random = ~gen, family = "poisson"),
  nested = list(y ~ env, random = ~ gen + env:gen),
  multi_way = list(y ~ env, random = ~ gen + env:rowf:colf),
  dsum_residual = list(y ~ env, random = ~gen, residual = ~ dsum(~ units | env)),
  diag = list(y ~ env, random = ~ diag(env):gen),
  slope = list(y ~ x + (x || gen)),
  ar1_field = list(y ~ env, random = ~ gen + id(env):ar1(rowf):ar1(colf)),
  spline = list(y ~ env, random = ~ gen + spl(x))
)

# Structures whose standard deviation the package leaves to INLA, by design.
.g11_engine_sd <- "spline"

# The likelihood and latent-field hyperparameters INLA estimated, keyed as
# prior_summary() keys its slots, with whether each is a precision.
.g11_hypers <- function(fit) {
  ah <- fit$inla$all.hyper
  out <- list()
  add <- function(prefix, hyper) {
    for (k in names(hyper)) {
      h <- hyper[[k]]
      if (isTRUE(h$fixed)) next
      out[[length(out) + 1L]] <<- data.frame(
        engine_slot = paste(prefix, k),
        precision = grepl("precision", h$name %||% "", fixed = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  for (j in seq_along(ah$family)) {
    add(paste0("family INLA.Data", j), ah$family[[j]]$hyper)
  }
  for (r in ah$random) {
    add(paste("random", r$hyperid), r$hyper)
  }
  do.call(rbind, out)
}

test_that("DC-11 guard: every INLA hyperparameter is in the prior record, on the package prior", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE,
    flexyBayes.silence_spatial_collapse_warning = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  d <- .g11_data()
  off <- character(0)
  for (nm in names(.g11_corpus)) {
    fit <- suppressWarnings(suppressMessages(do.call(
      flexybayes, c(.g11_corpus[[nm]], list(data = d, backend = "inla", verbose = FALSE))
    )))
    if (!identical(fit$inla$.args$num.threads, "1:1")) {
      off <- c(off, paste0(nm, ": threads not recorded as set"))
    }
    slots <- prior_summary(fit)$slots
    # One slot per estimated hyperparameter, group correlations included.
    n_slots <- sum(slots$section %in% c("family", "random", "group"))
    if (n_slots != nrow(fit$inla$summary.hyperpar)) {
      off <- c(off, sprintf(
        "%s: %d hyperparameters estimated, %d in the prior record",
        nm, nrow(fit$inla$summary.hyperpar), n_slots
      ))
    }
    h <- .g11_hypers(fit)
    for (i in seq_len(nrow(h))) {
      s <- slots[slots$engine_slot == h$engine_slot[i], , drop = FALSE]
      if (!nrow(s)) {
        off <- c(off, paste0(nm, ": ", h$engine_slot[i], " not in the prior record"))
        next
      }
      engine_sd <- h$precision[i] && identical(s$source[1L], "engine_default")
      if (engine_sd && !s$structure[1L] %in% .g11_engine_sd) {
        off <- c(off, paste0(
          nm, ": ", h$engine_slot[i], " (", s$structure[1L],
          ") runs on INLA's default prior ", s$readback[1L]
        ))
      }
    }
  }
  expect_identical(off, character(0))
})

test_that("DC-11 guard: the package's INLA draws set the sampler arguments", {
  root <- testthat::test_path("..", "..")
  skip_if(!dir.exists(file.path(root, "R")), "package source tree not available")
  off <- character(0)
  calls_seen <- c(hyperpar = 0L, posterior = 0L)
  for (f in list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)) {
    walk <- function(x) {
      if (!is.call(x)) {
        if (is.pairlist(x) || is.list(x)) for (a in as.list(x)) if (!missing(a)) walk(a)
        return(invisible())
      }
      fn <- paste(deparse(x[[1L]]), collapse = "")
      args <- as.list(x)[-1L]
      passes_dots <- any(vapply(args, identical, logical(1L), quote(...)))
      if (fn %in% c("INLA::inla.hyperpar.sample", "inla.hyperpar.sample")) {
        calls_seen[["hyperpar"]] <<- calls_seen[["hyperpar"]] + 1L
        if (!isTRUE(x[["improve.marginals"]])) {
          off <<- c(off, paste0(basename(f), ": ", paste(deparse(x), collapse = "")))
        }
      }
      if (fn == ".fb_inla_posterior_sample" && !passes_dots) {
        calls_seen[["posterior"]] <<- calls_seen[["posterior"]] + 1L
        if (is.null(x[["num.threads"]])) {
          off <<- c(off, paste0(basename(f), ": ", paste(deparse(x), collapse = "")))
        }
      }
      for (a in as.list(x)) if (!missing(a)) walk(a)
    }
    for (e in parse(f, keep.source = FALSE)) walk(e)
  }
  expect_true(all(calls_seen > 0L))
  expect_identical(off, character(0))
})
