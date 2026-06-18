# skip_ledger.R -- catalogue of skip() sites across the test suite.
#
# Parses tests/testthat/test-*.R via parse() + utils::getParseData() and
# classifies each skip site against the audit-named reason taxonomy.
# Writes inst/skip-ledger.md (one Markdown table per reason class).
#
# Recognised call shapes:
#   testthat::skip(<string>)                         skip(<string>)
#   testthat::skip_if(<cond>, <message>)             skip_if(...)
#   testthat::skip_if_not(<cond>, <message>)         skip_if_not(...)
#   testthat::skip_if_not_installed(<pkg>)           skip_if_not_installed(...)
#   testthat::skip_on_cran()                         skip_on_cran()
#   testthat::skip_on_ci()                           skip_on_ci()
#   project-local helpers .skip_if_no_<pkg>()        skip_if_no_<pkg>()
#
# Reason taxonomy (five audit classes from
# feedback_v0.3.7_to_v0.4.0_2026-05-25.md plus two in-portfolio extensions
# ratified at v0.3.8 Phase D for documented flakes):
#   unavailable_optional_backend, cran_ci_time_budget, stress_only,
#   external_toolchain_unavailable, deliberately_deferred_feature,
#   test_infrastructure_flake, stochastic_rng_dependent.
#
# Run as the closing step of tools/tally.R, as the second step of
# tools/profiles/release.R, or standalone via
# `Rscript flexyBayes/tools/skip_ledger.R` from the workspace root.

# Reason taxonomy -------------------------------------------------------

.skip_ledger_reason_classes <- c(
  unavailable_optional_backend = "Unavailable optional backend",
  cran_ci_time_budget = "CRAN / CI time budget",
  stress_only = "Stress-only (env-gated)",
  external_toolchain_unavailable = "External toolchain unavailable",
  deliberately_deferred_feature = "Deliberately deferred feature",
  test_infrastructure_flake = "Test-infrastructure flake",
  stochastic_rng_dependent = "Stochastic / RNG-dependent"
)

.skip_ledger_optional_backends <- c(
  "greta",
  "INLA",
  "inla",
  "brms",
  "rstan",
  "cmdstanr",
  "nimble",
  "JAGS",
  "rjags",
  "lme4",
  "mgcv",
  "Matrix",
  "asreml",
  "fst"
)

# Classifier ------------------------------------------------------------

.skip_classify_site <- function(call_name, reason_text) {
  # Explicit flake-prefix tokens win over every other rule -------------
  if (grepl("[flake-test-infrastructure]", reason_text, fixed = TRUE)) {
    return("test_infrastructure_flake")
  }
  if (grepl("[flake-stochastic-rng]", reason_text, fixed = TRUE)) {
    return("stochastic_rng_dependent")
  }

  # Structural call-name routing ---------------------------------------
  if (call_name %in% c("skip_on_cran", "skip_on_ci")) {
    return("cran_ci_time_budget")
  }
  if (call_name == "skip_if_not_installed") {
    return(.classify_install_skip(call_name, reason_text))
  }
  if (grepl("^skip_if_no_", call_name)) {
    return(.classify_install_skip(call_name, reason_text))
  }

  # Reason-text keyword routing ----------------------------------------
  rt <- tolower(reason_text)
  if (
    grepl(
      "flexybayes_run_stress|flexybayes_run_huge_n|stress-only|stress only",
      rt
    )
  ) {
    return("stress_only")
  }
  if (grepl("cran|on[ _]ci|ci runner|time budget", rt)) {
    return("cran_ci_time_budget")
  }
  if (grepl("toolchain|jags|stan binary|external|pandoc|latex", rt)) {
    return("external_toolchain_unavailable")
  }
  for (pkg in .skip_ledger_optional_backends) {
    if (
      grepl(
        paste0("\\b", pkg, "\\b"),
        reason_text,
        perl = TRUE,
        ignore.case = TRUE
      )
    ) {
      return("unavailable_optional_backend")
    }
  }
  if (grepl("defer|future|todo|not yet|v0\\.[0-9]", rt)) {
    return("deliberately_deferred_feature")
  }

  # Fallback: deliberately deferred (the audit's catch-all class) ------
  "deliberately_deferred_feature"
}

.classify_install_skip <- function(call_name, reason_text) {
  helper_pkg <- sub("^skip_if_no_", "", call_name)
  helper_pkg <- sub("_quiet$", "", helper_pkg)
  if (helper_pkg != call_name && nzchar(helper_pkg)) {
    return("unavailable_optional_backend")
  }
  for (pkg in .skip_ledger_optional_backends) {
    if (
      grepl(
        paste0("\\b", pkg, "\\b"),
        reason_text,
        perl = TRUE,
        ignore.case = TRUE
      )
    ) {
      return("unavailable_optional_backend")
    }
  }
  "external_toolchain_unavailable"
}

# Parser ----------------------------------------------------------------

.parse_one_test_file <- function(path) {
  exprs <- parse(path, keep.source = TRUE)
  pd <- utils::getParseData(exprs)
  if (is.null(pd) || !nrow(pd)) {
    return(data.table::data.table())
  }
  symbols <- pd[pd$token == "SYMBOL_FUNCTION_CALL", ]
  if (!nrow(symbols)) {
    return(data.table::data.table())
  }
  skip_idx <- grep("^skip($|_)", symbols$text)
  if (!length(skip_idx)) {
    return(data.table::data.table())
  }

  out <- data.table::data.table(
    file = basename(path),
    line = symbols$line1[skip_idx],
    call = symbols$text[skip_idx],
    reason = vapply(
      skip_idx,
      function(i) {
        .extract_str_arg(pd, symbols$id[i])
      },
      character(1L)
    )
  )
  out[,
    reason_class := vapply(
      seq_len(.N),
      function(i) {
        .skip_classify_site(call[i], reason[i])
      },
      character(1L)
    )
  ]
  out[]
}

.extract_str_arg <- function(parse_data, call_id) {
  # The SYMBOL_FUNCTION_CALL token's parent is the call expression for
  # bare calls (`skip("x")`) but for namespace-qualified calls
  # (`testthat::skip("x")`) it is the inner namespace expr; the call
  # expression is one level further out.  The call expr is the smallest
  # enclosing expr whose direct children include `(`.
  call_expr <- .find_call_expr(parse_data, call_id)
  if (is.na(call_expr)) {
    return("")
  }
  collected <- integer(0L)
  frontier <- call_expr
  repeat {
    children <- parse_data$id[parse_data$parent %in% frontier]
    if (!length(children)) {
      break
    }
    collected <- c(collected, children)
    frontier <- children
  }
  strs <- parse_data[
    parse_data$id %in% collected & parse_data$token == "STR_CONST",
  ]
  if (!nrow(strs)) {
    return("")
  }
  raw <- strs$text[1L]
  raw <- sub("^['\"]", "", raw)
  raw <- sub("['\"]$", "", raw)
  raw
}

.find_call_expr <- function(parse_data, call_id) {
  cur <- parse_data$parent[parse_data$id == call_id][1L]
  while (!is.na(cur) && cur != 0L) {
    direct_kids <- parse_data[parse_data$parent == cur, ]
    if (any(direct_kids$token == "'('")) {
      return(cur)
    }
    cur <- parse_data$parent[parse_data$id == cur][1L]
  }
  NA_integer_
}

# Markdown writer -------------------------------------------------------

.write_ledger_markdown <- function(dt, output_path) {
  header <- c(
    "# Skip ledger",
    "",
    paste0(
      "Auto-generated by `tools/skip_ledger.R`; refreshed by ",
      "`tools/tally.R` and `tools/profiles/release.R`."
    ),
    "",
    paste0(
      "Reason taxonomy: five audit classes ",
      "(feedback_v0.3.7_to_v0.4.0_2026-05-25.md ",
      "Sec. Test And Benchmark Advice) plus two in-portfolio ",
      "extensions ratified at v0.3.8 Phase D for documented flakes."
    ),
    "",
    paste0(
      "Snapshot ",
      format(Sys.Date(), "%Y-%m-%d"),
      ": ",
      nrow(dt),
      " skip sites across ",
      length(unique(dt$file)),
      " test files."
    ),
    ""
  )
  body <- character(0L)
  for (class in names(.skip_ledger_reason_classes)) {
    sub_dt <- dt[reason_class == class]
    if (!nrow(sub_dt)) {
      next
    }
    body <- c(
      body,
      paste0(
        "## ",
        .skip_ledger_reason_classes[[class]],
        " (",
        nrow(sub_dt),
        ")"
      ),
      "",
      "| File | Line | Call | Reason |",
      "|---|---|---|---|"
    )
    for (i in seq_len(nrow(sub_dt))) {
      reason <- sub_dt$reason[i]
      if (!nzchar(reason)) {
        reason <- "(no string arg)"
      }
      if (nchar(reason) > 80L) {
        reason <- paste0(substr(reason, 1L, 77L), "...")
      }
      reason <- gsub("\\|", "\\\\|", reason)
      body <- c(
        body,
        sprintf(
          "| `%s` | %d | `%s()` | %s |",
          sub_dt$file[i],
          sub_dt$line[i],
          sub_dt$call[i],
          reason
        )
      )
    }
    body <- c(body, "")
  }
  writeLines(c(header, body), output_path)
  invisible(NULL)
}

.write_empty_ledger <- function(output_path) {
  writeLines(c("# Skip ledger", "", "No skip sites detected."), output_path)
  invisible(NULL)
}

# Validation ------------------------------------------------------------

.check_skip_ledger_paths <- function(test_dir, output_path) {
  if (!dir.exists(test_dir)) {
    stop(call. = FALSE, "`test_dir` does not exist: ", test_dir)
  }
  out_dir <- dirname(output_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  invisible(NULL)
}

# Public entry point ----------------------------------------------------

build_skip_ledger <- function(
  test_dir = "tests/testthat",
  output_path = "inst/skip-ledger.md"
) {
  .check_skip_ledger_paths(test_dir, output_path)
  files <- list.files(test_dir, pattern = "^test-.*\\.R$", full.names = TRUE)
  if (!length(files)) {
    .write_empty_ledger(output_path)
    return(invisible(data.table::data.table()))
  }
  parts <- lapply(files, .parse_one_test_file)
  dt <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  if (!nrow(dt)) {
    .write_empty_ledger(output_path)
    return(invisible(dt))
  }
  data.table::setorder(dt, reason_class, file, line)
  .write_ledger_markdown(dt, output_path)
  invisible(dt[])
}

# Taxonomy-lock assertion ----------------------------------------------

# .assert_taxonomy_lock() --- the v0.4.0 section-11 gate. Classifies
# every skip site in `test_dir` and asserts each resolves to exactly one
# of the seven reason classes (the five audit classes plus the two
# in-portfolio extensions). Returns the classified table invisibly on
# success; stops with the offending sites listed on failure. The
# classifier's catch-all (deliberately_deferred_feature) makes an
# out-of-vocabulary class structurally impossible, so this gate is a
# standing guarantee: an edit that broadened the vocabulary without
# ratifying it -- e.g. a classifier branch returning a new class string
# not added to .skip_ledger_reason_classes -- would trip here.
.assert_taxonomy_lock <- function(test_dir) {
  tmp <- tempfile(fileext = ".md")
  on.exit(unlink(tmp), add = TRUE)
  dt <- build_skip_ledger(test_dir = test_dir, output_path = tmp)
  if (!nrow(dt)) {
    return(invisible(dt))
  }
  valid <- names(.skip_ledger_reason_classes)
  bad <- dt[!dt$reason_class %in% valid, ]
  if (nrow(bad)) {
    msg <- paste(
      sprintf(
        "  %s:%d (%s) -> '%s'",
        bad$file,
        bad$line,
        bad$call,
        bad$reason_class
      ),
      collapse = "\n"
    )
    stop(
      call. = FALSE,
      "skip_ledger: ",
      nrow(bad),
      " skip site(s) classify outside the ",
      "seven-class taxonomy:\n",
      msg
    )
  }
  invisible(dt)
}

# Standalone CLI --------------------------------------------------------

if (sys.nframe() == 0L) {
  pkg_root <- if (dir.exists("flexyBayes/tests/testthat")) {
    "flexyBayes"
  } else if (dir.exists("tests/testthat")) {
    "."
  } else {
    stop(
      call. = FALSE,
      "skip_ledger.R: could not locate tests/testthat; ",
      "run from package root or workspace root"
    )
  }
  args <- commandArgs(trailingOnly = TRUE)
  test_dir <- file.path(pkg_root, "tests/testthat")

  if ("--assert-taxonomy-lock" %in% args) {
    dt <- .assert_taxonomy_lock(test_dir)
    counts <- table(factor(
      dt$reason_class,
      levels = names(.skip_ledger_reason_classes)
    ))
    cat(
      "skip_ledger taxonomy-lock: PASS --",
      nrow(dt),
      "sites, all within the seven-class taxonomy.\n"
    )
    for (cls in names(counts)) {
      cat(sprintf("  %-32s %d\n", cls, counts[[cls]]))
    }
  } else {
    out <- build_skip_ledger(
      test_dir = test_dir,
      output_path = file.path(pkg_root, "inst/skip-ledger.md")
    )
    cat(sprintf(
      "skip_ledger: %d sites across %d files -> %s\n",
      nrow(out),
      length(unique(out$file)),
      file.path(pkg_root, "inst/skip-ledger.md")
    ))
  }
}
