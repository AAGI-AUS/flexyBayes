# Class guard DC-13: no package code reads a field or an argument by
# partial matching.
#
# The main routes run with warnPartialMatchDollar, warnPartialMatchArgs and
# warnPartialMatchAttr on: INLA fits whose terms lack a field that a longer
# field extends (a:b:cc carries `vars` and no `var`), every S3 method and
# fb_* accessor of those fits, the brms prior and code builders, fb_plan()
# and lgm_gate(). A partial-match warning raised in a flexyBayes function
# fails the guard; one raised inside another package is not counted
# (DEF-198: term$var returned a three-way term's `vars`).

.g13_data <- function() {
  set.seed(3)
  d <- expand.grid(a = factor(1:3), b = factor(1:2), cc = factor(1:2), rep = 1:4)
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(12)[interaction(d$a, d$b, d$cc)] +
    stats::rnorm(nrow(d))
  d$n <- stats::rpois(nrow(d), 3)
  d
}

.g13_models <- list(
  three_way_vc = list(fixed = y ~ x, random = ~ a + a:b + a:b:cc, prior_vc_sd = 2),
  three_way = list(fixed = y ~ x + a * b, random = ~ a:b:cc),
  slope = list(fixed = y ~ x + b + (x || a)),
  diag = list(fixed = y ~ x, random = ~ diag(b):a),
  poisson = list(fixed = n ~ x + a:x, random = ~ b + a:b:cc, family = "poisson")
)

# Side effects (plot) or a second fit (anova, fb_refit).
.g13_skip <- c("plot", "anova", "fb_refit")

# The flexyBayes function that raised a warning, or "" when another
# package's code raised it. A C-level warning reaches its handler through
# .signalSimpleWarning(); the frames below it are searched for the first
# one outside base R, because `x$a` passed as an argument is evaluated
# inside the base function that receives it, such as paste().
.g13_raised_in_package <- function(calls) {
  ns <- asNamespace("flexyBayes")
  sig <- which(vapply(calls, function(cl) {
    identical(cl[[1L]], quote(.signalSimpleWarning))
  }, logical(1L)))
  k <- if (length(sig)) max(sig) - 1L else length(calls) - 2L
  for (i in rev(seq_len(max(k, 0L)))) {
    env <- environment(sys.function(i))
    top <- if (is.null(env)) baseenv() else topenv(env)
    if (identical(top, baseenv()) || identical(top, .BaseNamespaceEnv)) {
      next
    }
    if (identical(top, ns)) {
      return(paste(deparse(calls[[i]][[1L]]), collapse = ""))
    }
    return("")
  }
  ""
}

test_that("DC-13 guard: no partial matching on the main routes", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  d <- .g13_data()
  ns <- asNamespace("flexyBayes")
  methods_inla <- unique(c(
    as.character(utils::.S3methods(class = "flexybayes_inla", envir = ns)),
    as.character(utils::.S3methods(class = "flexybayes", envir = ns))
  ))
  generics <- unique(sub("\\.flexybayes(_inla)?$", "", methods_inla))
  accessors <- grep("^fb_", getNamespaceExports("flexyBayes"), value = TRUE)
  takes_fit <- vapply(accessors, function(nm) {
    first <- names(formals(getExportedValue("flexyBayes", nm)))[1L]
    isTRUE(first %in% c("object", "fit", "x"))
  }, logical(1L))
  calls <- setdiff(c(generics, accessors[takes_fit]), .g13_skip)

  hits <- character(0)
  routes_run <- 0L
  record <- function(route, expr) {
    withCallingHandlers(
      tryCatch(
        suppressMessages(force(expr)),
        error = function(e) NULL
      ),
      warning = function(w) {
        if (grepl("^partial (argument )?match", conditionMessage(w))) {
          fn <- .g13_raised_in_package(sys.calls())
          if (nzchar(fn)) {
            hits <<- c(hits, paste0(
              route, ": ", fn, "(): ", conditionMessage(w), " in ",
              paste(deparse(conditionCall(w)), collapse = "")
            ))
          }
        }
        invokeRestart("muffleWarning")
      }
    )
  }

  withr::local_options(
    warnPartialMatchDollar = TRUE,
    warnPartialMatchArgs = TRUE,
    warnPartialMatchAttr = TRUE
  )
  for (nm in names(.g13_models)) {
    m <- .g13_models[[nm]]
    args <- c(m, list(data = d, backend = "inla", verbose = FALSE))
    fit <- record(paste0(nm, " fit"), do.call(flexybayes, args))
    if (inherits(fit, "flexybayes_inla")) {
      routes_run <- routes_run + 1L
      for (g in calls) {
        fun <- if (exists(g, envir = ns)) get(g, envir = ns) else match.fun(g)
        record(paste0(nm, " ", g), utils::capture.output(do.call(fun, list(fit))))
      }
    }
    fb <- record(paste0(nm, " ir"), .fb_fit_ir(fit))
    record(paste0(nm, " lgm_gate"), lgm_gate(fb))
    record(paste0(nm, " brms priors"), .priors_to_brms(
      NULL, fb, prior_fixed_sd = 100, prior_vc_sd = 2
    ))
    record(paste0(nm, " brms formula"), .fb_to_brms_formula(fb))
    for (be in c("inla", "brms")) {
      record(paste0(nm, " fb_plan ", be), utils::capture.output(print(fb_plan(
        m$fixed, data = d, random = m$random, backend = be,
        family = m$family %||% "gaussian"
      ))))
    }
  }

  expect_identical(routes_run, length(.g13_models))
  expect_identical(hits, character(0))
})
