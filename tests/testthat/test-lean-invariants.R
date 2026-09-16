# Tier-1 invariants of the package source.
#
# Four properties a security scan found true on 2026-09-15 and that nothing
# kept true: no network call, no evaluation of parsed text, no global option
# or random-number state left changed by an exported call, and no Import
# outside a listed set. The first two walk the syntax tree of every closure in
# the namespace, which is the same object on a source tree and on an installed
# package, so they also hold under R CMD check. Each test was shown to fail on
# an injected violation before it was trusted.

# `pkg::fun(...)` names c(pkg, fun); `fun(...)` names c("", fun); any other
# call head (an anonymous function, `x$f()`) names c("", "").
.lean_call_name <- function(cl) {
  head <- cl[[1L]]
  if (is.symbol(head)) {
    return(c(pkg = "", fun = as.character(head)))
  }
  if (is.call(head) && length(head) == 3L && is.symbol(head[[1L]]) &&
        as.character(head[[1L]]) %in% c("::", ":::")) {
    return(c(pkg = as.character(head[[2L]]), fun = as.character(head[[3L]])))
  }
  c(pkg = "", fun = "")
}

.lean_walk <- function(e, visit) {
  if (is.call(e)) {
    visit(e)
    for (a in as.list(e)) {
      if (!missing(a)) .lean_walk(a, visit)
    }
  }
  invisible(NULL)
}

# Every finding `visit_one()` reports, over the body and formals of every
# closure in the namespace, as "object: finding".
.lean_scan <- function(visit_one) {
  ns <- asNamespace("flexyBayes")
  hits <- character(0)
  for (nm in ls(ns, all.names = TRUE)) {
    obj <- get(nm, envir = ns)
    if (!is.function(obj) || is.primitive(obj)) next
    found <- character(0)
    visit <- function(cl) found <<- c(found, visit_one(cl))
    .lean_walk(body(obj), visit)
    for (a in as.list(formals(obj))) {
      if (!missing(a)) .lean_walk(a, visit)
    }
    if (length(found)) hits <- c(hits, paste0(nm, ": ", found))
  }
  unique(hits)
}

.lean_network_funs <- c(
  "url", "download.file", "curlGetHeaders", "socketConnection",
  "make.socket", "serverSocket", "socketAccept", "browseURL"
)
.lean_network_pkgs <- c("curl", "httr", "httr2", "RCurl", "crul", "websocket")
.lean_url_readers <- c(
  "file", "gzcon", "readLines", "readRDS", "load", "source", "scan",
  "read.table", "read.csv", "read.dcf", "fread"
)

.lean_network_visitor <- function(cl) {
  nm <- .lean_call_name(cl)
  if (nm[["fun"]] %in% .lean_network_funs || nm[["pkg"]] %in% .lean_network_pkgs) {
    return(paste0("calls ", nm[["pkg"]], if (nzchar(nm[["pkg"]])) "::", nm[["fun"]], "()"))
  }
  if (identical(nm[["fun"]], "do.call") && length(cl) > 1L &&
        is.character(cl[[2L]]) && cl[[2L]] %in% .lean_network_funs) {
    return(paste0("do.call()s ", cl[[2L]], "()"))
  }
  if (nm[["fun"]] %in% .lean_url_readers) {
    for (a in as.list(cl)[-1L]) {
      if (!missing(a) && is.character(a) && length(a) == 1L &&
            grepl("^(https?|ftp)://", a)) {
        return(paste0("reads a URL with ", nm[["fun"]], "()"))
      }
    }
  }
  character(0)
}

.lean_parse_visitor <- function(cl) {
  nm <- .lean_call_name(cl)
  if (!(nm[["fun"]] %in% c("eval", "evalq")) || length(cl) < 2L) {
    return(character(0))
  }
  inner <- character(0)
  for (a in as.list(cl)[-1L]) {
    if (!missing(a)) {
      .lean_walk(a, function(c2) inner <<- c(inner, .lean_call_name(c2)[["fun"]]))
    }
  }
  if (any(inner %in% c("parse", "str2lang", "str2expression"))) {
    return(paste0(nm[["fun"]], "() of parsed text"))
  }
  character(0)
}

test_that("no function in the namespace makes a network call", {
  expect_identical(.lean_scan(.lean_network_visitor), character(0))
})

test_that("no function in the namespace evaluates parsed text", {
  expect_identical(.lean_scan(.lean_parse_visitor), character(0))
})

test_that("exported calls that need no engine leave options and the random-number state unchanged", {
  withr::local_preserve_seed()
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  trial <- expand.grid(row = factor(1:4), col = factor(1:3))
  trial$yield <- seq_len(nrow(trial))
  holed <- trial[-c(2L, 7L), ]
  G <- crossprod(matrix(c(2, 1, 0, 1, 3, 1, 0, 1, 4), 3L, 3L))
  # fb_engine("brms") is left out: it loads the brms namespace, whose own
  # load-time effects are not this package's.
  calls <- list(
    fb_prior = quote(fb_prior(sigma ~ pc(upper = 2, prob = 0.05))),
    fb_cov = quote(fb_cov(G, type = "dense")),
    fb_approx = quote(fb_approx("low_rank_smooth", rank = 5L)),
    fb_engine = quote(fb_engine("inla")),
    fb_penalty = quote(fb_penalty(term(everything()) ~ ridge(scale = 1))),
    fb_resamples = quote(fb_resamples(unit = ~ g, v = 3L)),
    fb_decision = quote(fb_decision(rule = "tolerance", tolerance = 0.1)),
    fb_refusals = quote(fb_refusals()),
    fb_backend_status = quote(fb_backend_status()),
    fb_from_brms = quote(fb_from_brms(y ~ x + (1 | g), data = d)),
    fb_from_asreml = quote(fb_from_asreml(y ~ x, random = ~ g, data = d)),
    fb_complete_grid = quote(fb_complete_grid(holed, ~ row * col, response = "yield")),
    fb_plan = quote(fb_plan(y ~ x + (1 | g), data = d)),
    flexybayes_plan = quote(flexybayes(y ~ x + (1 | g), data = d, plan = TRUE, verbose = FALSE))
  )
  for (nm in names(calls)) {
    set.seed(20260916L)
    opt_before <- options()
    seed_before <- get(".Random.seed", envir = globalenv())
    suppressMessages(eval(calls[[nm]]))
    expect_identical(options(), opt_before, label = paste0("options() after ", nm, "()"))
    expect_identical(
      get(".Random.seed", envir = globalenv()), seed_before,
      label = paste0(".Random.seed after ", nm, "()")
    )
  }
})

test_that("Imports and namespace imports stay within the listed set", {
  allowed <- c(
    "cli", "data.table", "generics",
    "stats", "methods", "graphics", "mgcv", "Matrix", "parallel", "utils"
  )
  dcf <- read.dcf(system.file("DESCRIPTION", package = "flexyBayes"), fields = "Imports")
  imports <- trimws(sub("\\(.*\\)", "", strsplit(dcf[1L, "Imports"], ",")[[1L]]))
  imports <- imports[nzchar(imports)]
  expect_identical(setdiff(imports, allowed), character(0))
  ns_imports <- setdiff(names(getNamespaceImports(asNamespace("flexyBayes"))), c("", "base"))
  expect_identical(setdiff(ns_imports, allowed), character(0))
})

test_that("every pkg:: call in the namespace names a declared or base R package", {
  # A package dropped from Imports can survive as a `pkg::` call, which the
  # DESCRIPTION and NAMESPACE reads above do not see.
  dcf <- read.dcf(
    system.file("DESCRIPTION", package = "flexyBayes"),
    fields = c("Imports", "Suggests")
  )
  declared <- trimws(sub("\\(.*\\)", "", unlist(strsplit(dcf[1L, ], ","))))
  declared <- declared[!is.na(declared) & nzchar(declared)]
  base_r <- c(
    "base", "compiler", "datasets", "graphics", "grDevices", "grid",
    "methods", "parallel", "splines", "stats", "stats4", "tcltk", "tools",
    "utils"
  )
  allowed <- c(declared, base_r, "flexyBayes")
  undeclared <- .lean_scan(function(cl) {
    nm <- .lean_call_name(cl)
    if (nzchar(nm[["pkg"]]) && !(nm[["pkg"]] %in% allowed)) {
      return(paste0("calls ", nm[["pkg"]], "::", nm[["fun"]], "()"))
    }
    character(0)
  })
  expect_identical(undeclared, character(0))
})
