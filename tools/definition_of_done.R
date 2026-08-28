# tools/definition_of_done.R - the frozen release contract for 0.10.0.
#
# "Are we ready?" has been a judgement call in every previous release
# cycle, and the answer changed each time because the criteria changed
# each time. This file makes it a command.
#
# Every criterion below is a row with an owner, a cost tier and, where
# the machine can decide, a check that returns TRUE or FALSE. Running
# this script prints the burn-down and exits non-zero while any blocking
# criterion is unmet. Stage 6 of the release plan is not an audit
# somebody performs; it is this script returning green.
#
# THE FREEZE RULE. This file is signed before package code moves. After
# that point new DEFECTS are admitted freely - they are failures of the
# criteria already here. New CRITERIA are not: adding one means editing
# this file, which is a visible act requiring the same sign-off. That is
# the whole anti-cycling device, and it only works if it is honoured.
#
# Cost tiers:
#   quick  - seconds, pure file inspection, runs by default
#   full   - minutes to hours, needs `--full` (gates, lint, check, CI)
#   manual - not the machine's to decide; reported, never auto-passed
#
# Usage:
#   Rscript tools/definition_of_done.R              # quick tier
#   Rscript tools/definition_of_done.R --full       # everything runnable
#   Rscript tools/definition_of_done.R --section B  # one section
#
# Run from the package root. Not shipped (`^tools$` in .Rbuildignore).

TARGET_VERSION <- "0.10.0"
TARGET_TIER <- "V2"
# 31 was the arithmetic of the originally proposed trim (44 - 11 - 2).
# Raised to 34 on 2026-08-28, signed off, after the pre-trim sweep found
# three of those eleven were wrongly grouped: fb_gblup_cv() is taught and
# executed in the met-and-genomics vignette, and genomic_summary() and
# triangulate_genomic() operate on posterior draws. The /rpkg E1 cap of 30
# is therefore missed by four, deliberately and on the record.
MAX_EXPORTS <- 34L

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

.dod_ok <- function(pass, detail = "") {
  list(pass = isTRUE(pass), detail = detail)
}

.dod_read <- function(path) {
  if (!file.exists(path)) {
    return(character(0))
  }
  readLines(path, warn = FALSE)
}

.dod_desc <- function(field) {
  dcf <- read.dcf("DESCRIPTION")
  if (!field %in% colnames(dcf)) {
    return(NA_character_)
  }
  as.character(dcf[1, field])
}

.dod_exports <- function() {
  ns <- .dod_read("NAMESPACE")
  hit <- grep("^export\\(", ns, value = TRUE)
  sub("^export\\((.*)\\)$", "\\1", hit)
}

.dod_r_files <- function() {
  list.files("R", pattern = "\\.R$", full.names = TRUE)
}

.dod_vignette_files <- function() {
  list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
}

# Comment text with the markers stripped and whitespace collapsed. A
# line-based grep on wrapped prose gives false PASSes - the phrase the
# check is looking for straddles two lines and is never seen. Every
# prose check goes through here.
.dod_prose <- function(path) {
  lines <- .dod_read(path)
  lines <- sub("^\\s*#+'?", "", lines)
  gsub("\\s+", " ", paste(lines, collapse = " "))
}

# Non-comment source lines only. A package named in a roxygen block is
# MENTIONED, not used; counting the mention as use is how an unbacked
# Suggests survives an audit.
.dod_code_lines <- function(paths) {
  src <- unlist(lapply(paths, .dod_read))
  src[!grepl("^\\s*#", src)]
}

.dod_last_commit_date <- function() {
  out <- suppressWarnings(tryCatch(
    system2("git", c("log", "-1", "--format=%cd", "--date=short"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)
  ))
  if (!length(out)) NA else as.Date(out[1])
}

# Every `name =` token appearing inside a call written in roxygen prose,
# paired with the function it claims to be calling. This is the check
# that `tools::checkDocFiles()` cannot make: it compares \usage against
# \arguments, so a wrong argument name in a DESCRIPTION paragraph (the
# `fb_plan(fixed = ...)` class) passes R CMD check untouched.
.dod_prose_call_args <- function() {
  out <- list()
  for (f in .dod_r_files()) {
    lines <- .dod_read(f)
    roxy <- grep("^\\s*#'", lines, value = TRUE)
    if (!length(roxy)) next
    calls <- regmatches(
      roxy,
      gregexpr("\\b(flexybayes|fb|fb_[a-z_]+)\\(([^()]*)\\)", roxy)
    )
    for (one in unlist(calls)) {
      fun <- sub("\\(.*$", "", one)
      body_txt <- sub("^[^(]*\\(", "", sub("\\)$", "", one))
      args <- regmatches(
        body_txt,
        gregexpr("[A-Za-z._][A-Za-z0-9._]*(?=\\s*=)", body_txt, perl = TRUE)
      )[[1]]
      if (length(args)) {
        out[[length(out) + 1L]] <- list(file = f, fun = fun, args = args)
      }
    }
  }
  out
}

.dod_formals_of <- function(fun) {
  if (!requireNamespace("flexyBayes", quietly = TRUE)) {
    return(NULL)
  }
  if (!exists(fun, envir = asNamespace("flexyBayes"), inherits = FALSE)) {
    return(NULL)
  }
  obj <- get(fun, envir = asNamespace("flexyBayes"))
  if (!is.function(obj)) {
    return(NULL)
  }
  names(formals(obj))
}

# ---------------------------------------------------------------------
# Section A - declared posture
# ---------------------------------------------------------------------

.dod_section_a <- list(
  list(
    id = "A1", cost = "quick", owner = "claude", blocking = TRUE,
    what = "validationTier is V2, not V3",
    check = function() {
      got <- .dod_desc("Config/rpkg/validationTier")
      .dod_ok(identical(got, TARGET_TIER), paste("found:", got))
    }
  ),
  list(
    id = "A2", cost = "quick", owner = "claude", blocking = TRUE,
    what = "no expired-version waiver text in inst/validation/README.md",
    check = function() {
      txt <- .dod_read("inst/validation/README.md")
      hit <- grep("waiver applies to|applies to \\*\\*0\\.", txt)
      .dod_ok(
        length(hit) == 0L,
        if (length(hit)) paste("line", paste(hit, collapse = ", ")) else "clean"
      )
    }
  ),
  list(
    id = "A3", cost = "quick", owner = "claude", blocking = TRUE,
    what = paste0("DESCRIPTION Version is ", TARGET_VERSION),
    check = function() {
      got <- .dod_desc("Version")
      .dod_ok(identical(got, TARGET_VERSION), paste("found:", got))
    }
  ),
  list(
    id = "A4", cost = "quick", owner = "claude", blocking = TRUE,
    what = "NEWS top section matches DESCRIPTION Version",
    check = function() {
      news <- .dod_read("NEWS.md")
      top <- grep("^# ", news, value = TRUE)[1]
      .dod_ok(
        !is.na(top) && grepl(TARGET_VERSION, top, fixed = TRUE),
        paste("found:", top)
      )
    }
  )
  # A5 RETIRED 2026-08-28, signed off. It required a
  # lifecycle::badge("experimental") on every export. Its purpose was to
  # stop API_STABILITY.md claiming badges that did not exist, and B7
  # already achieves that by stating the lifecycle stage without claiming
  # a badge. The package has no man/figures/, so the criterion would have
  # meant shipping SVG assets to state, for a fifth time, a fact already
  # carried by the DESCRIPTION version, the README, API_STABILITY.md and
  # the startup message. CRAN does not require badges. Retiring it is
  # recorded here rather than deleted silently, because retiring a
  # criterion is the same visible act as adding one.
)

# ---------------------------------------------------------------------
# Section B - claim-surface parity, the anti-drift guards
# ---------------------------------------------------------------------

.dod_section_b <- list(
  list(
    id = "B1", cost = "quick", owner = "claude", blocking = TRUE,
    what = "every argument named in roxygen prose exists in formals()",
    check = function() {
      bad <- character(0)
      n_checked <- 0L
      n_dots <- 0L
      for (one in .dod_prose_call_args()) {
        fmls <- .dod_formals_of(one$fun)
        if (is.null(fmls)) next
        # a function taking `...` accepts arbitrary names, so nothing
        # static can be concluded about them
        if ("..." %in% fmls) {
          n_dots <- n_dots + 1L
          next
        }
        n_checked <- n_checked + 1L
        miss <- setdiff(one$args, fmls)
        if (length(miss)) {
          bad <- c(bad, paste0(
            basename(one$file), ": ", one$fun, "(",
            paste(miss, collapse = ", "), ")"
          ))
        }
      }
      .dod_ok(
        length(bad) == 0L,
        paste0(
          if (length(bad)) paste(unique(bad), collapse = "; ") else "clean",
          " [", n_checked, " calls checked, ", n_dots,
          " skipped for `...`]"
        )
      )
    }
  ),
  list(
    id = "B2", cost = "quick", owner = "claude", blocking = TRUE,
    what = "package-level Rd vignette count matches vignettes/",
    check = function() {
      n_real <- length(.dod_vignette_files())
      txt <- .dod_prose("R/flexyBayes-package.R")
      words <- c("One", "Two", "Three", "Four", "Five", "Six", "Seven",
                 "Eight", "Nine", "Ten", "Eleven", "Twelve", "Thirteen")
      hit <- vapply(words, function(w) {
        grepl(paste0("\\b", w, " vignettes?\\b"), txt, ignore.case = TRUE)
      }, logical(1))
      claimed <- words[hit]
      ok <- length(claimed) == 1L &&
        identical(tolower(claimed), tolower(words[n_real]))
      .dod_ok(ok, paste0("claims '", paste(claimed, collapse = "/"),
                         "', actual ", n_real))
    }
  ),
  list(
    id = "B3", cost = "quick", owner = "claude", blocking = TRUE,
    what = "every package named in DESCRIPTION Description is a real dep",
    check = function() {
      descr <- .dod_desc("Description")
      deps <- paste(.dod_desc("Imports"), .dod_desc("Suggests"))
      quoted <- unique(unlist(regmatches(
        descr, gregexpr("(?<=')[A-Za-z][A-Za-z0-9.]+(?=')", descr, perl = TRUE)
      )))
      named <- unique(c(quoted, "effectsize", "emmeans", "marginaleffects"))
      named <- named[vapply(named, function(p) {
        grepl(p, descr, fixed = TRUE)
      }, logical(1))]
      # CODE lines only: a package named in a roxygen block is mentioned,
      # not used, and counting the mention is how `effectsize` survived
      # three audits while no method existed
      code <- .dod_code_lines(c(
        .dod_r_files(),
        list.files("tests", recursive = TRUE, full.names = TRUE,
                   pattern = "\\.R$")
      ))
      used <- vapply(named, function(p) {
        any(grepl(p, code, fixed = TRUE))
      }, logical(1))
      bad <- named[!used]
      .dod_ok(length(bad) == 0L,
              if (length(bad)) paste("claimed, unused:",
                                     paste(bad, collapse = ", ")) else "clean")
    }
  ),
  list(
    id = "B4", cost = "quick", owner = "claude", blocking = TRUE,
    what = "every Imports package is actually used",
    check = function() {
      imp <- trimws(strsplit(.dod_desc("Imports"), ",")[[1]])
      imp <- sub("\\s*\\(.*\\)$", "", imp)
      # Base packages ship with R and are exempt. `methods` especially:
      # the package reaches S4 through `@` slot access (R/structured_cov.R,
      # `L@uplo`), which no grep for `pkg::` or for an imported symbol
      # can see.
      imp <- setdiff(imp, c("", "stats", "methods", "utils", "graphics"))
      # Two ways this guard used to pass on a dead Import, both fixed
      # here. It read comments as use, and it accepted the mere presence
      # of the package NAME anywhere in the source -- including inside a
      # refusal-message string. `splines` sat in Imports with `bs()`
      # called nowhere and this check reported clean, which is the defect
      # 089a64d fixed by hand for coda, recurring inside the guard
      # written to stop it recurring. A symbol taken by importFrom()
      # now counts only when that symbol is actually called.
      src <- .dod_code_lines(.dod_r_files())
      ns <- .dod_read("NAMESPACE")
      .used <- function(p) {
        if (any(grepl(paste0(p, "::"), src, fixed = TRUE))) {
          return(TRUE)
        }
        tags <- grep(paste0("^importFrom\\(", p, ","), ns, value = TRUE)
        if (!length(tags)) {
          return(FALSE)
        }
        syms <- sub(paste0("^importFrom\\(", p, ",(.*)\\)$"), "\\1", tags)
        syms <- trimws(unlist(strsplit(syms, ",")))
        any(vapply(syms, function(sym) {
          any(grepl(paste0("(^|[^[:alnum:]_.])", sym, "\\s*\\("), src))
        }, logical(1)))
      }
      unused <- imp[!vapply(imp, .used, logical(1))]
      .dod_ok(length(unused) == 0L,
              if (length(unused)) paste("unused:",
                                        paste(unused, collapse = ", "))
              else "clean")
    }
  ),
  list(
    id = "B5", cost = "quick", owner = "claude", blocking = TRUE,
    what = "CITATION.cff / codemeta.json / inst/CITATION agree on version",
    check = function() {
      v <- .dod_desc("Version")
      cff <- any(grepl(paste0("version: ", v), .dod_read("CITATION.cff")))
      cmj <- any(grepl(paste0("\"version\": \"", v),
                       .dod_read("codemeta.json")))
      cit <- any(grepl(v, .dod_read("inst/CITATION"), fixed = TRUE))
      .dod_ok(cff && cmj && cit,
              paste0("cff=", cff, " codemeta=", cmj, " CITATION=", cit))
    }
  ),
  list(
    id = "B6", cost = "quick", owner = "claude", blocking = TRUE,
    what = "CITATION.cff date-released is not before the last commit",
    check = function() {
      cff <- grep("date-released", .dod_read("CITATION.cff"), value = TRUE)
      iso <- regmatches(cff, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", cff))
      d <- if (length(iso)) as.Date(iso[1]) else as.Date(NA)
      last <- .dod_last_commit_date()
      .dod_ok(!is.na(d) && !is.na(last) && d >= last,
              paste0("cff=", d, " last commit=", last))
    }
  ),
  list(
    id = "B7", cost = "quick", owner = "claude", blocking = TRUE,
    what = "API_STABILITY.md makes no claim contradicted by R/",
    check = function() {
      api <- .dod_read("API_STABILITY.md")
      claims_badges <- any(grepl("badge", api, ignore.case = TRUE))
      has_badges <- any(grepl("lifecycle::badge",
                              unlist(lapply(.dod_r_files(), .dod_read))))
      .dod_ok(!claims_badges || has_badges,
              paste0("claims badges=", claims_badges,
                     " badges present=", has_badges))
    }
  ),
  list(
    id = "B8", cost = "quick", owner = "claude", blocking = TRUE,
    what = "no I14 flag words on any shipped surface",
    check = function() {
      pat <- "honesty check|graduation|SA-gate|vibe-check"
      # text only: inst/ carries binaries, and grepl() on those emits an
      # invalid-input warning that would drown a real signal
      # R/ and man/ ship in the tarball too, and were not being read: a
      # flag word in a roxygen comment or a generated Rd page passed a
      # check whose criterion says "any shipped surface".
      files <- c("README.md", "NEWS.md", "DESCRIPTION",
                 list.files("inst", recursive = TRUE, full.names = TRUE,
                            pattern = "\\.(md|R|Rmd|txt|csv|yaml|yml)$"),
                 .dod_r_files(),
                 list.files("man", pattern = "\\.Rd$", full.names = TRUE),
                 .dod_vignette_files())
      files <- files[file.exists(files) & !dir.exists(files)]
      bad <- files[vapply(files, function(f) {
        any(grepl(pat, .dod_read(f), ignore.case = TRUE))
      }, logical(1))]
      .dod_ok(length(bad) == 0L,
              if (length(bad)) paste(bad, collapse = ", ") else "clean")
    }
  ),
  list(
    id = "B9", cost = "quick", owner = "claude", blocking = TRUE,
    what = "README makes no stability claim contradicting experimental",
    check = function() {
      txt <- .dod_read("README.md")
      .dod_ok(
        !any(grepl("stable release", txt, ignore.case = TRUE)),
        if (any(grepl("stable release", txt, ignore.case = TRUE)))
          "README claims a stable release" else "clean"
      )
    }
  ),
  list(
    id = "B10", cost = "quick", owner = "claude", blocking = TRUE,
    what = "workflow files make no false claim about their own gates",
    check = function() {
      yml <- list.files(".github/workflows", full.names = TRUE)
      claims_zero <- any(vapply(yml, function(f) {
        any(grepl("to zero", .dod_read(f), fixed = TRUE))
      }, logical(1)))
      .dod_ok(!claims_zero,
              if (claims_zero) "a workflow comment asserts a lint state"
              else "clean")
    }
  ),
  list(
    id = "B11", cost = "quick", owner = "claude", blocking = TRUE,
    what = "every advertised structure spelling is executed by the grid",
    check = function() {
      # The ledger's `cell` column carries prose slugs, so a spelling is
      # only visible in the grid SOURCE. `ped()` and `vm()` share one
      # matrix cell, which is how `ped()` came to be advertised and
      # never live-fitted.
      grid <- "tools/execution_grid.R"
      if (!file.exists(grid)) {
        return(.dod_ok(FALSE, "no execution_grid.R"))
      }
      # The oracle is the ledger, which records cells that RAN, not the
      # source, which only declares them. Reading the source meant
      # `k_cap <- diag(6)` in fixture setup satisfied the `diag`
      # criterion permanently, and a deleted cell whose call string
      # survived in a comment would still pass. The ledger is committed,
      # so a fresh clone has it; fall back to the source only when it is
      # absent, and say so rather than passing quietly.
      # The oracle is the recorded run, not the source. Reading the
      # source meant `k_cap <- diag(6)` in fixture setup satisfied the
      # `diag` criterion permanently, and a deleted cell whose call
      # string survived in a comment would still pass. grid_results.csv
      # carries a `code` column holding the call each cell actually
      # executed, and it is committed, so a fresh clone has it. Only
      # that column is read: the `message` column quotes spellings back
      # inside refusal text and would satisfy this check by accident.
      results <- "inst/validation/execution_grid/grid_results.csv"
      advertised <- c("vm", "ped", "us", "diag", "dsum", "ar1", "spl", "at")
      if (file.exists(results)) {
        tab <- tryCatch(
          utils::read.csv(results, stringsAsFactors = FALSE),
          error = function(e) NULL
        )
        if (is.null(tab) || !("code" %in% names(tab))) {
          return(.dod_ok(FALSE, "grid_results.csv has no `code` column"))
        }
        src <- tab$code
        where <- "the recorded run"
      } else {
        src <- grep("flexybayes\\(|random =|residual =",
                    .dod_code_lines(grid), value = TRUE)
        where <- "grid source (no run recorded)"
      }
      missing <- advertised[!vapply(advertised, function(a) {
        any(grepl(paste0("\\b", a, "\\("), src))
      }, logical(1))]
      .dod_ok(length(missing) == 0L,
              if (length(missing)) paste("advertised, never executed:",
                                         paste(missing, collapse = ", "))
              else paste("clean, read from the", where))
    }
  ),
  list(
    id = "B12", cost = "quick", owner = "claude", blocking = TRUE,
    what = "every export is exercised by at least one test file",
    check = function() {
      # The count in C1 is a symptom; this is the principle. An export
      # with no test contact is a claim with no oracle, which is how a
      # quarter of the public API came to sit outside the 408-cell grid.
      # Static reference, not a live fit - named for what it checks.
      tests <- list.files("tests", pattern = "\\.R$", recursive = TRUE,
                          full.names = TRUE)
      txt <- unlist(lapply(tests, .dod_read))
      exps <- .dod_exports()
      untouched <- exps[!vapply(exps, function(e) {
        any(grepl(paste0("\\b", e, "\\b"), txt))
      }, logical(1))]
      .dod_ok(
        length(untouched) == 0L,
        if (length(untouched))
          paste0(length(untouched), " untested: ",
                 paste(utils::head(untouched, 10), collapse = ", "))
        else paste0("all ", length(exps), " exports referenced")
      )
    }
  ),
  list(
    id = "B13", cost = "quick", owner = "claude", blocking = TRUE,
    what = "no Rd example calls a function the package does not export",
    check = function() {
      # Withdrawing a function from the API does not touch its examples,
      # and R CMD check runs the examples of internal Rd pages too. The
      # 0.10.0 trim unexported eight functions and left ten example
      # blocks calling them unqualified, which check reported as two
      # ERRORs and no test could see. The rule: a documented-but-
      # unexported function may appear in an example only as
      # `flexyBayes:::name()`.
      rds <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
      if (!length(rds)) {
        return(.dod_ok(FALSE, "no man/ pages"))
      }
      exps <- .dod_exports()
      rd_name <- function(rd) {
        hit <- grep("^\\\\name\\{", .dod_read(rd), value = TRUE)
        if (!length(hit)) return(NA_character_)
        sub("^\\\\name\\{(.*)\\}.*$", "\\1", hit[[1]])
      }
      documented <- vapply(rds, rd_name, character(1), USE.NAMES = FALSE)
      internal <- setdiff(stats::na.omit(documented), exps)
      if (!length(internal)) {
        return(.dod_ok(TRUE, "no documented internals"))
      }
      tmp <- tempfile()
      on.exit(unlink(tmp), add = TRUE)
      offenders <- character(0)
      for (rd in rds) {
        ok <- tryCatch({ tools::Rd2ex(rd, out = tmp); TRUE },
                       error = function(e) FALSE)
        if (!ok || !file.exists(tmp)) next
        ex <- paste(.dod_read(tmp), collapse = "\n")
        unlink(tmp)
        if (!nzchar(trimws(ex))) next
        for (nm in internal) {
          # bare call, not already reached through the namespace
          if (grepl(paste0("(?<![\\w:.])", nm, "\\s*\\("), ex, perl = TRUE)) {
            offenders <- c(offenders, paste0(basename(rd), ":", nm))
          }
        }
      }
      .dod_ok(
        length(offenders) == 0L,
        if (length(offenders))
          paste0(length(offenders), " unqualified: ",
                 paste(utils::head(offenders, 8), collapse = ", "))
        else paste0(length(internal), " documented internals, all qualified")
      )
    }
  ),
  list(
    id = "B14", cost = "quick", owner = "claude", blocking = TRUE,
    what = "every non-internal help page appears in the pkgdown index",
    check = function() {
      # pkgdown's build_reference_index() errors when a topic is neither
      # listed in _pkgdown.yml nor marked internal, and that is a full-
      # tier gate costing twenty minutes to reach. Withdrawing the
      # Dirichlet and GEV constructors left their four S3 method pages
      # indexed with no reachable class, and the failure surfaced only
      # in the release bake. This is the same rule, in the quick tier.
      yml <- "_pkgdown.yml"
      rds <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
      if (!file.exists(yml) || !length(rds)) {
        return(.dod_ok(FALSE, "no _pkgdown.yml or no man/ pages"))
      }
      idx <- trimws(sub("^\\s*-\\s*", "", grep("^\\s*-\\s+[A-Za-z.]",
                                             .dod_read(yml), value = TRUE)))
      missing <- character(0)
      for (rd in rds) {
        txt <- .dod_read(rd)
        if (any(grepl("\\\\keyword\\{internal\\}", txt))) {
          next
        }
        hit <- grep("^\\\\name\\{", txt, value = TRUE)
        if (!length(hit)) {
          next
        }
        nm <- sub("^\\\\name\\{(.*)\\}.*$", "\\1", hit[[1]])
        # reexports.Rd documents borrowed generics and is never indexed
        if (identical(nm, "reexports")) {
          next
        }
        if (!(nm %in% idx)) {
          missing <- c(missing, nm)
        }
      }
      .dod_ok(length(missing) == 0L,
              if (length(missing))
                paste0(length(missing), " unindexed: ",
                       paste(utils::head(missing, 6), collapse = ", "))
              else paste0(length(rds), " pages, index complete"))
    }
  ),
  list(
    id = "B15", cost = "quick", owner = "claude", blocking = TRUE,
    what = "every export is named in API_STABILITY.md",
    check = function() {
      # The stability document is the contract a user reads before
      # depending on a function. Three exports were absent from it --
      # fb_complete_grid(), genomic_summary() and ranef() -- so their
      # stability stage was undeclared while the document presented
      # itself as the inventory.
      doc <- "API_STABILITY.md"
      if (!file.exists(doc)) {
        return(.dod_ok(FALSE, "no API_STABILITY.md"))
      }
      txt <- paste(.dod_read(doc), collapse = "\n")
      exps <- .dod_exports()
      absent <- exps[!vapply(exps, function(e) {
        grepl(paste0("`", e), txt, fixed = TRUE)
      }, logical(1))]
      .dod_ok(length(absent) == 0L,
              if (length(absent))
                paste0(length(absent), " undeclared: ",
                       paste(utils::head(absent, 8), collapse = ", "))
              else paste0("all ", length(exps), " exports declared"))
    }
  )
)

# ---------------------------------------------------------------------
# Section C - consistency invariants (the "functional and consistent"
# half of the brief, expressed as things a user can hit)
# ---------------------------------------------------------------------

.dod_section_c <- list(
  list(
    id = "C1", cost = "quick", owner = "claude", blocking = TRUE,
    what = paste0("export count is at most ", MAX_EXPORTS),
    check = function() {
      n <- length(.dod_exports())
      .dod_ok(n <= MAX_EXPORTS, paste0(n, " exports"))
    }
  ),
  list(
    id = "C2", cost = "quick", owner = "claude", blocking = TRUE,
    what = "no export always abstains (3 known unconditional abstainers)",
    check = function() {
      # Every export whose only reachable outcome is a message and an
      # empty return. fb_structured_cov() joined the list at 0.10.0: no
      # active engine emits an fa() term, so its single return path is
      # invisible(list()). A name is added here only by withdrawing it.
      zombies <- c("fb_met_summary", "fb_log_posterior",
                   "fb_structured_cov")
      present <- intersect(zombies, .dod_exports())
      .dod_ok(length(present) == 0L,
              if (length(present)) paste("still exported:",
                                         paste(present, collapse = ", "))
              else "clean")
    }
  ),
  list(
    id = "C3", cost = "quick", owner = "claude", blocking = TRUE,
    what = "every export has an \\examples{} block",
    check = function() {
      # An export is not always documented on man/<name>.Rd: `fb` is an
      # alias on flexybayes.Rd, and the re-exported generics live on
      # reexports.Rd. Resolve each export to the page that carries its
      # \alias{}, and exempt re-exports, which are documented upstream.
      rds <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
      pages <- lapply(rds, .dod_read)
      names(pages) <- rds
      page_of <- function(e) {
        tag <- paste0("\\alias{", e, "}")
        hit <- names(pages)[vapply(pages, function(x) {
          any(grepl(tag, x, fixed = TRUE))
        }, logical(1))]
        if (length(hit)) hit[[1]] else NA_character_
      }
      exps <- .dod_exports()
      no_ex <- exps[vapply(exps, function(e) {
        pg <- page_of(e)
        if (is.na(pg)) return(TRUE)
        if (basename(pg) == "reexports.Rd") return(FALSE)
        !any(grepl("\\examples", pages[[pg]], fixed = TRUE))
      }, logical(1))]
      .dod_ok(length(no_ex) == 0L,
              if (length(no_ex)) paste0(length(no_ex), " without examples: ",
                                        paste(utils::head(no_ex, 8),
                                              collapse = ", "))
              else "clean")
    }
  ),
  list(
    id = "C4", cost = "full", owner = "claude", blocking = TRUE,
    what = "a boundary-pinned variance component warns, on any Gaussian term",
    command = "testthat::test_file('tests/testthat/test-boundary-collapse.R')",
    check = function() {
      # The oracle here was file.exists(), so the criterion passed while
      # every test in the file could be failing. It runs them now.
      f <- "tests/testthat/test-boundary-collapse.R"
      if (!file.exists(f)) {
        return(.dod_ok(FALSE, "no test asserting the iid case warns"))
      }
      if (!requireNamespace("pkgload", quietly = TRUE) ||
            !requireNamespace("testthat", quietly = TRUE)) {
        return(.dod_ok(FALSE, "pkgload/testthat absent; cannot run C4"))
      }
      suppressMessages(pkgload::load_all(".", quiet = TRUE))
      res <- as.data.frame(testthat::test_file(f, reporter = "silent"))
      bad <- sum(res$failed) + sum(res$error)
      .dod_ok(bad == 0L,
              sprintf("%d failed, %d passed", bad, sum(res$passed)))
    }
  ),
  list(
    id = "C5", cost = "quick", owner = "claude", blocking = TRUE,
    what = "na_action documented behaviour matches the measured behaviour",
    check = function() {
      # Measurement revised 2026-08-28. The original check demanded the
      # "same posterior" sentence be deleted. It should not be: the
      # identity is a correct statement about the posterior. What the
      # measurement showed is that it is not a promise about what an
      # optimiser returns, so the claim must travel with that caveat.
      #
      # Widened 2026-08-28: reading only R/na_action.R is what let the
      # uncaveated identity stand in README.md and in R/emit_brms.R after
      # the correction landed. Every surface that states the identity is
      # checked, because a correction applied to one file is not a
      # correction.
      # Derived, not listed. The hand-written list read na_action.R,
      # emit_brms.R, README and the vignettes, and was still short by
      # R/flexybayes.R -- the entry point, and the most-read statement of
      # the identity of the four. A list is only as good as the person
      # writing it remembering every surface; deriving the set means a
      # new surface cannot escape by not being thought of.
      surfaces <- c(.dod_r_files(), "README.md", .dod_vignette_files())
      surfaces <- surfaces[file.exists(surfaces)]
      # Patterns, not literals. Three hand-copied sentences missed the
      # entry point's phrasing ("posterior FOR THE MODEL PARAMETERS is
      # the same whether"), so widening the surface list changed nothing
      # -- the guard could see the file and still not see the claim. The
      # regex tolerates the words a writer puts between "posterior" and
      # "is the same".
      states <- "posterior[^.]{0,60}is (then )?the same"
      caveats <- "not (a promise )?about what[[:space:]]+an optimiser returns"
      bare <- character(0)
      for (f in surfaces) {
        txt <- .dod_prose(f)
        if (!grepl(states, txt)) {
          next
        }
        if (!grepl(caveats, txt)) {
          bare <- c(bare, f)
        }
      }
      .dod_ok(length(bare) == 0L,
              if (length(bare))
                paste("identity stated without the optimiser caveat:",
                      paste(bare, collapse = ", "))
              else paste(length(surfaces), "surfaces checked"))
    }
  )
)

# ---------------------------------------------------------------------
# Section D - machine evidence
# ---------------------------------------------------------------------

.dod_section_d <- list(
  list(id = "D1", cost = "full", owner = "claude", blocking = TRUE,
       what = "testthat suite: 0 fail, 0 warn",
       command = "review/phase_reports_0100/scripts/bake_and_gate.sh gates"),
  list(id = "D2", cost = "full", owner = "claude", blocking = TRUE,
       what = "R CMD check --as-cran, INLA present: 0/0/<=1",
       command = "R CMD check --as-cran flexyBayes_0.10.0.tar.gz"),
  list(id = "D3", cost = "full", owner = "claude", blocking = TRUE,
       what = "R CMD check --as-cran, INLA absent: 0/0/<=1",
       command = "INLA masked; same check"),
  list(id = "D4", cost = "full", owner = "claude", blocking = TRUE,
       what = "execution grid: 0 divergent, 0 untyped, 0 crash",
       command = "Rscript tools/execution_grid.R"),
  list(id = "D5", cost = "full", owner = "claude", blocking = TRUE,
       what = "lintr::lint_package() returns 0 lints (lint.yaml hard-gates)",
       # with the dev tree loaded, as CI does via local::. -- against a
       # stale installed build object_usage_linter reports phantoms
       command = paste("Rscript -e 'pkgload::load_all(quiet=TRUE);",
                       "print(length(lintr::lint_package()))'")),
  list(id = "D6", cost = "full", owner = "claude", blocking = TRUE,
       what = "pkgdown builds clean",
       command = "Rscript -e 'pkgdown::build_site()'"),
  list(id = "D7", cost = "full", owner = "claude", blocking = TRUE,
       what = "clean-room: tarball built from the tag tree, checked in scratch",
       command = "R CMD build from a fresh clone, check in a scratch dir"),
  list(id = "D8", cost = "full", owner = "claude", blocking = TRUE,
       what = "all four CI workflows green on a local dry run",
       command = "act, or a scratch private repo; none has ever executed"),
  list(id = "D9", cost = "full", owner = "claude", blocking = FALSE,
       what = "coverage measured once with INLA and brms present (recorded)",
       command = "Rscript -e 'covr::package_coverage()'")
)

# ---------------------------------------------------------------------
# Section E - not the machine's to decide. Reported, never auto-passed.
# ---------------------------------------------------------------------

.dod_section_e <- list(
  list(id = "E1", cost = "manual", owner = "Max", blocking = TRUE,
       what = "D8 governance: licence consult, co-author consent, public home"),
  list(id = "E2", cost = "manual", owner = "Max", blocking = TRUE,
       what = "E19 human De-Vibe pass on the release tree"),
  list(id = "E3", cost = "manual", owner = "Max", blocking = TRUE,
       what = "win-builder devel and release archived (outward)"),
  list(id = "E4", cost = "manual", owner = "Max", blocking = TRUE,
       what = "R-hub Linux archived (outward)"),
  list(id = "E5", cost = "manual", owner = "Max", blocking = TRUE,
       what = "CRAN submission and confirmation (outward)")
)

# ---------------------------------------------------------------------
# Explicitly OUT of 0.10.0. Not criteria. Admitting one means editing
# this file, which is the point.
# ---------------------------------------------------------------------

.dod_out_of_scope <- list(
  c("V3 floors: MFD, SBC n_sim >= 1000, null-recovery, metamorphics",
    "0.11"),
  c("Bayesian FA GxE on brms, the incumbent's headline model",
    "0.11+"),
  c("INLA-SPDE Matern; hierarchically shrunk INLA residual; Tweedie",
    "0.11+"),
  c("S3 solver-aware latent-field preflight",
    "0.11+"),
  c("The 11 ML side-track exports (gwas, gblup, gev, dirichlet)",
    "deferred"),
  c("omit + informative missingness: the 3M-evaluation spin",
    "0.11, documented at 0.10.0"),
  c("augment/omit variance disagreement under MCAR: the mechanism",
    "0.11, claim corrected at 0.10.0"),
  c("covr driven to an 85 per cent floor",
    "not planned, measured once"),
  c("Full MET live evidence beyond 120 rows",
    "0.11+")
)

# ---------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------

.dod_all <- function() {
  c(.dod_section_a, .dod_section_b, .dod_section_c,
    .dod_section_d, .dod_section_e)
}

.dod_run_one <- function(item, do_full) {
  if (identical(item$cost, "manual")) {
    return(list(status = "MANUAL", detail = paste("owner:", item$owner)))
  }
  if (identical(item$cost, "full") && !do_full) {
    return(list(status = "SKIP", detail = item$command %||% "full tier"))
  }
  if (is.null(item$check)) {
    # A blocking criterion with no in-process check is UNMET until its
    # evidence is recorded, not passed. It used to report SKIP, and SKIP
    # did not count, so D1-D9 -- the suite, both check runs, the grid,
    # lint, pkgdown, the clean room, CI and coverage -- could every one
    # be red while this script exited 0 and the header claimed it exits
    # non-zero while any blocking criterion is unmet.
    return(list(status = if (isTRUE(item$blocking)) "PENDING" else "SKIP",
                detail = item$command %||% "no check yet"))
  }
  res <- tryCatch(item$check(), error = function(e) {
    .dod_ok(FALSE, paste("check errored:", conditionMessage(e)))
  })
  list(status = if (res$pass) "PASS" else "FAIL", detail = res$detail)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.dod_main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  do_full <- "--full" %in% args
  sect <- NA_character_
  if (any(grepl("^--section", args))) {
    sect <- sub("^--section=?", "",
                grep("^--section", args, value = TRUE)[1])
  }

  if (!file.exists("DESCRIPTION")) {
    stop("run from the package root (no DESCRIPTION here)")
  }

  items <- .dod_all()
  if (!is.na(sect) && nzchar(sect)) {
    items <- items[vapply(items, function(x) {
      startsWith(x$id, toupper(sect))
    }, logical(1))]
  }

  cat("\n")
  cat("flexyBayes definition of done - target ", TARGET_VERSION,
      ", tier ", TARGET_TIER, "\n", sep = "")
  cat("tree version ", .dod_desc("Version"),
      " | tier ", .dod_desc("Config/rpkg/validationTier"),
      " | exports ", length(.dod_exports()),
      " | tier run: ", if (do_full) "full" else "quick", "\n", sep = "")
  cat(strrep("-", 78), "\n", sep = "")

  n_fail <- 0L
  n_pass <- 0L
  n_manual <- 0L
  n_skip <- 0L
  n_pending <- 0L

  for (item in items) {
    res <- .dod_run_one(item, do_full)
    n_fail <- n_fail + (res$status == "FAIL" && isTRUE(item$blocking))
    n_pass <- n_pass + (res$status == "PASS")
    n_manual <- n_manual + (res$status == "MANUAL")
    n_skip <- n_skip + (res$status == "SKIP")
    n_pending <- n_pending + (res$status == "PENDING")
    cat(sprintf("%-4s %-6s %-62s\n", item$id, res$status,
                substr(item$what, 1L, 62L)))
    if (nzchar(res$detail) && res$status != "PASS") {
      cat(sprintf("          %s\n", substr(res$detail, 1L, 100L)))
    }
  }

  cat(strrep("-", 78), "\n", sep = "")
  cat(sprintf(
    "PASS %d | FAIL %d (blocking) | PENDING %d | MANUAL %d | SKIP %d\n",
    n_pass, n_fail, n_pending, n_manual, n_skip))
  if (do_full && n_pending > 0L) {
    cat("PENDING is unmet, not passed: run each item's command and ",
        "record its artefact.\n", sep = "")
  }

  cat("\nExplicitly OUT of ", TARGET_VERSION,
      " - admitting one means editing this file:\n", sep = "")
  for (row in .dod_out_of_scope) {
    cat(sprintf("  - %-62s -> %s\n", row[1], row[2]))
  }
  cat("\n")

  # Quick tier gates on what it can read. The full tier is the release
  # question, so there an unmet blocking criterion -- failed or merely
  # unevidenced -- exits non-zero.
  if (n_fail > 0L || (do_full && n_pending > 0L)) {
    quit(status = 1L)
  }
  invisible(NULL)
}

# Run only when invoked directly, so the file can be sourced into a
# session to inspect an individual check. `--file=` is set by
# `Rscript tools/definition_of_done.R` and not by `Rscript -e source(...)`.
.dod_invoked_directly <- function() {
  any(grepl("^--file=.*definition_of_done\\.R$",
            commandArgs(trailingOnly = FALSE)))
}

if (.dod_invoked_directly()) {
  .dod_main()
}
