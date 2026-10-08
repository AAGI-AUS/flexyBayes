# Class guard DC-12: one fit class lacks a method or field the others have.
#
# Swept over every generic the package registers a fit-class method for:
#
# 1. each fit class (per-row INLA, aggregated INLA, brms) reaches a
#    flexyBayes method, never another package's default (DEF-011, DEF-020,
#    DEF-039);
# 2. an argument one class's method takes is not swallowed by another
#    class's `...` without being read (DEF-069: tidy() on INLA fits had no
#    `effects` argument and returned fixed effects);
# 3. where an installed suggested package exports the same generic, that
#    generic dispatches to the package's methods too, so attaching the
#    other package does not mask them (DEF-192);
# 4. per-row and aggregated INLA fits answer summary(), fitted() and
#    residuals() in the same shape, and every fit carries its parsed model
#    and data (DEF-061, DEF-123, DEF-144).

.g12_fit_classes <- c(
  "flexybayes", "flexybayes_inla", "flexybayes_aggregated", "flexybayes_brms"
)

# The class a brms fit carries (R/emit_brms.R); a long-tier block checks
# it against a real fit.
.g12_brms_class <- c("flexybayes_brms", "flexybayes", "list")

.g12_methods <- function() {
  ns <- asNamespace("flexyBayes")
  s3 <- getNamespaceInfo(ns, "S3methods")
  s3 <- s3[s3[, 2L] %in% .g12_fit_classes, , drop = FALSE]
  data.frame(
    generic = sub(".*::", "", s3[, 1L]),
    home = ifelse(grepl("::", s3[, 1L]), sub("::.*", "", s3[, 1L]), ""),
    class = s3[, 2L],
    fun = ifelse(is.na(s3[, 3L]), paste0(s3[, 1L], ".", s3[, 2L]), s3[, 3L]),
    stringsAsFactors = FALSE
  )
}

# The method a class vector reaches for one generic among the package's
# registrations, or NA.
.g12_resolve <- function(m, generic, cls) {
  for (k in cls) {
    i <- which(m$generic == generic & m$class == k & m$home == "")
    if (length(i)) {
      return(m$fun[i[1L]])
    }
  }
  NA_character_
}

# Packages whose methods for one generic share a named argument: the
# arguments a method reads, by name or through `...`.
.g12_reads_dots <- function(fn) {
  "..." %in% all.names(body(fn))
}

# Method and argument pairs that swallow a sibling's argument today.
.g12_open <- c(
  "predict.flexybayes_inla re_formula" = "DEF-221",
  "predict.flexybayes_inla summary" = "DEF-221",
  "vcov.flexybayes n_samples" = "DEF-221",
  "fb_as_draws_simple.flexybayes_brms n_samples" = "DEF-221"
)

.g12_data <- function() {
  set.seed(12)
  d <- data.frame(
    g = factor(rep(seq_len(6L), each = 10L)),
    f = factor(rep(c("a", "b"), 30L))
  )
  d$y <- 1 + stats::rnorm(6L)[d$g] + stats::rnorm(60L)
  d$n <- stats::rpois(60L, 3)
  d
}

test_that("DC-12 guard: every fit class reaches a package method", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  d <- .g12_data()
  per_row <- suppressMessages(flexybayes(y ~ f, random = ~g, data = d, verbose = FALSE))
  aggregated <- suppressMessages(flexybayes(
    n ~ f, random = ~g, data = d, family = "poisson", verbose = FALSE
  ))
  classes <- list(
    per_row = class(per_row),
    aggregated = class(aggregated),
    brms = .g12_brms_class
  )
  expect_true(all(vapply(classes, function(k) "flexybayes" %in% k, logical(1L))))
  expect_true("flexybayes_inla" %in% classes$aggregated)

  m <- .g12_methods()
  generics <- unique(m$generic[m$home == ""])
  unreached <- character(0)
  for (g in generics) {
    for (nm in names(classes)) {
      if (is.na(.g12_resolve(m, g, classes[[nm]]))) {
        unreached <- c(unreached, paste0(g, "() on the ", nm, " fit"))
      }
    }
  }
  expect_gt(length(generics), 30L)
  expect_identical(unreached, character(0))
})

test_that("DC-12 guard: no fit-class method swallows an argument a sibling takes", {
  ns <- asNamespace("flexyBayes")
  m <- .g12_methods()
  classes <- list(
    per_row = c("flexybayes_inla", "flexybayes"),
    aggregated = c("flexybayes_aggregated", "flexybayes_inla", "flexybayes"),
    brms = .g12_brms_class
  )
  swallowed <- character(0)
  for (g in unique(m$generic[m$home == ""])) {
    funs <- unique(stats::na.omit(vapply(
      classes, function(k) .g12_resolve(m, g, k), character(1L)
    )))
    args <- lapply(funs, function(f) setdiff(names(formals(get(f, envir = ns))), "..."))
    names(args) <- funs
    every <- unique(unlist(args))
    for (f in funs) {
      lacking <- setdiff(every, args[[f]])
      if (!length(lacking) || .g12_reads_dots(get(f, envir = ns))) next
      swallowed <- c(swallowed, paste(f, lacking))
    }
  }
  new_swallowed <- setdiff(swallowed, names(.g12_open))
  now_read <- setdiff(names(.g12_open), swallowed)
  expect_identical(new_swallowed, character(0), label = "arguments swallowed")
  expect_identical(
    now_read,
    character(0),
    label = "open-list pairs that no longer swallow (remove them from .g12_open)"
  )
})

test_that("DC-12 guard: a generic another package exports still reaches the package methods", {
  m <- .g12_methods()
  generics <- unique(m$generic)
  root <- testthat::test_path("..", "..")
  desc <- if (file.exists(file.path(root, "DESCRIPTION"))) {
    file.path(root, "DESCRIPTION")
  } else {
    system.file("DESCRIPTION", package = "flexyBayes")
  }
  suggested <- trimws(sub(
    "\\(.*", "", strsplit(read.dcf(desc, fields = "Suggests")[1L, 1L], ",")[[1L]]
  ))
  classes <- list(
    per_row = c("flexybayes_inla", "flexybayes"),
    aggregated = c("flexybayes_aggregated", "flexybayes_inla", "flexybayes"),
    brms = .g12_brms_class
  )
  masked <- character(0)
  checked <- 0L
  for (p in suggested) {
    if (!requireNamespace(p, quietly = TRUE)) next
    for (g in intersect(getNamespaceExports(p), generics)) {
      fn <- getExportedValue(p, g)
      if (!is.function(fn) || !"UseMethod" %in% all.names(body(fn))) next
      home <- environment(fn)
      if (identical(home, asNamespace("flexyBayes"))) next
      table <- get(".__S3MethodsTable__.", envir = home)
      checked <- checked + 1L
      for (nm in names(classes)) {
        hit <- any(vapply(classes[[nm]], function(k) {
          exists(paste0(g, ".", k), envir = table, inherits = FALSE)
        }, logical(1L)))
        if (!hit) {
          masked <- c(masked, paste0(p, "::", g, "() on the ", nm, " fit"))
        }
      }
    }
  }
  expect_gt(checked, 5L)
  expect_identical(masked, character(0))
})

test_that("DC-12 guard: per-row and aggregated INLA fits answer in one shape", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE,
    flexyBayes.inla_threads = "1:1"
  )
  d <- .g12_data()
  fits <- list(
    per_row = suppressMessages(flexybayes(
      n ~ f, random = ~g, data = d, family = "poisson", aggregate = FALSE,
      verbose = FALSE
    )),
    aggregated = suppressMessages(flexybayes(
      n ~ f, random = ~g, data = d, family = "poisson", verbose = FALSE
    ))
  )
  expect_true(inherits(fits$aggregated, "flexybayes_aggregated"))
  shapes <- lapply(fits, function(fit) {
    s <- suppressMessages(summary(fit))
    list(
      summary = names(s),
      varcomp = names(s$varcomp),
      fitted = length(stats::fitted(fit)),
      residuals = length(stats::residuals(fit)),
      ir = !is.null(.fb_fit_ir(fit)),
      data = nrow(fit$glm$data %||% fit$data)
    )
  })
  expect_identical(shapes$aggregated, shapes$per_row)
  expect_identical(shapes$per_row$fitted, nrow(d))
  expect_true(shapes$per_row$ir)
})

test_that("DC-12 guard: a brms fit carries the class and slots the guard assumes", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_long_tier()
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  fit <- suppressWarnings(suppressMessages(flexybayes(
    y ~ f, random = ~g, data = .g12_data(), backend = "brms", chains = 1L,
    n_samples = 200L, warmup = 200L, seed = 12L, verbose = FALSE,
    mcmc_verbose = FALSE
  )))
  expect_identical(class(fit), .g12_brms_class)
  expect_false(is.null(.fb_fit_ir(fit)))
  expect_identical(nrow(fit$data), 60L)
})
