# Class guard DC-05: two calls with the same seed give the same result.
#
# Every S3 method the package registers for its INLA fit class, and every
# exported fb_* accessor that takes a fit, is called twice on one small fit
# with set.seed(1) before each call. A result that differs between the two
# calls draws from a stream the R seed does not control (DEF-005: vcov()
# took inla.posterior.sample() draws, which set.seed() does not fix).

.g05_fit <- function() {
  set.seed(3)
  d <- data.frame(
    site = factor(rep(1:6, each = 12)),
    trt = factor(rep(c("a", "b", "c"), 24)),
    x = stats::rnorm(72)
  )
  d$y <- 1 +
    c(0, 0.5, -0.3)[d$trt] +
    0.4 * d$x +
    stats::rnorm(6)[d$site] +
    stats::rnorm(72)
  suppressWarnings(suppressMessages(flexybayes(
    y ~ trt + x,
    random = ~site,
    data = d,
    backend = "inla",
    verbose = FALSE
  )))
}

# Side effects (plot, print), a refit (update, fb_refit) or a second fit (anova).
.g05_skip <- c("plot", "print", "update", "fb_refit", "anova")

# Accessors that need more than the fit.
.g05_args <- list(
  vpredict = list(r ~ sd_site^2 / (sd_site^2 + sigma^2), n_draws = 200L),
  fb_variance_draws = list(200L)
)

# Routes known to draw. The guard fails if any of them did not run.
.g05_must_run <- c("vcov", "predict", "vpredict", "fb_variance_draws")

test_that("DC-05 guard: every accessor of an INLA fit repeats under one seed", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  fit <- .g05_fit()

  ns <- asNamespace("flexyBayes")
  methods_inla <- unique(c(
    as.character(utils::.S3methods(class = "flexybayes_inla", envir = ns)),
    as.character(utils::.S3methods(class = "flexybayes", envir = ns))
  ))
  generics <- unique(sub("\\.flexybayes(_inla)?$", "", methods_inla))
  accessors <- grep("^fb_", getNamespaceExports("flexyBayes"), value = TRUE)
  takes_fit <- vapply(
    accessors,
    function(nm) {
      first <- names(formals(getExportedValue("flexyBayes", nm)))[1L]
      isTRUE(first %in% c("object", "fit", "x"))
    },
    logical(1L)
  )
  calls <- setdiff(c(generics, accessors[takes_fit]), .g05_skip)

  ran <- character(0)
  for (nm in calls) {
    fun <- if (exists(nm, envir = ns)) get(nm, envir = ns) else match.fun(nm)
    args <- c(list(fit), .g05_args[[nm]])
    once <- function() {
      set.seed(1)
      tryCatch(
        suppressWarnings(suppressMessages(do.call(fun, args))),
        error = function(e) e
      )
    }
    a <- once()
    if (inherits(a, "error")) {
      next
    }
    b <- once()
    expect_identical(a, b, label = paste0(nm, "() on the second call"))
    ran <- c(ran, nm)
  }
  expect_true(
    all(.g05_must_run %in% ran),
    label = paste0(
      "drawing routes reached (missing: ",
      paste(setdiff(.g05_must_run, ran), collapse = ", "),
      ")"
    )
  )
})
