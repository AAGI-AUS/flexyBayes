# Class guard DC-08: every call the documentation shows can be made.
#
# Reads the R code a reader is shown -- \code{} and \examples in every help
# page, fenced and inline code in README.md, API_STABILITY.md,
# inst/KNOWN_ISSUES.md and the vignette sources -- and checks each call to a
# flexyBayes function: the function exists, an example calls it as a user
# can (exported, or through :::), and every named argument reaches a formal.
# A named argument that lands in a `...` the function never passes on is
# dropped, so it counts as unmatched. The help pages once showed
# fb_plan(fixed = ...), whose `...` is unused (4e04e3d4), and ten help pages
# called withdrawn functions by bare name in their examples (8a00cc14).

.g08_root <- function() {
  candidate <- testthat::test_path("..", "..")
  if (file.exists(file.path(candidate, "DESCRIPTION")) &&
    dir.exists(file.path(candidate, "man"))) {
    return(normalizePath(candidate))
  }
  NA_character_
}

# The \code{} and \examples blocks of one help page, with their kind.
.g08_rd_snippets <- function(file) {
  rd <- tools::parse_Rd(file, permissive = TRUE)
  text <- character(0)
  kind <- character(0)
  walk <- function(x) {
    tag <- attr(x, "Rd_tag")
    if (identical(tag, "\\code") || identical(tag, "\\examples")) {
      text[[length(text) + 1L]] <<- paste(unlist(x), collapse = "")
      kind[[length(kind) + 1L]] <<- if (identical(tag, "\\examples")) "example" else "code"
      return(invisible())
    }
    if (is.list(x)) {
      for (el in x) walk(el)
    }
  }
  walk(rd)
  data.frame(text = text, kind = kind, stringsAsFactors = FALSE)
}

# Fenced R chunks (run as a user runs them) and inline code spans of one
# Markdown or R Markdown file.
.g08_md_snippets <- function(file) {
  lines <- readLines(file, warn = FALSE)
  fence <- grepl("^\\s*```", lines)
  opens <- grepl("^\\s*```\\s*(\\{r|r\\b|R\\b)", lines)
  chunks <- character(0)
  inside <- FALSE
  buf <- character(0)
  is_r <- FALSE
  for (i in seq_along(lines)) {
    if (fence[i]) {
      if (!inside) {
        inside <- TRUE
        is_r <- opens[i]
        buf <- character(0)
      } else {
        if (is_r) chunks[[length(chunks) + 1L]] <- paste(buf, collapse = "\n")
        inside <- FALSE
      }
      next
    }
    if (inside) buf <- c(buf, lines[i])
  }
  prose <- paste(lines[!fence], collapse = "\n")
  spans <- regmatches(prose, gregexpr("`[^`\n]+`", prose))[[1]]
  spans <- sub("^r ", "", gsub("^`|`$", "", spans))
  data.frame(
    text = c(chunks, spans),
    kind = c(rep("example", length(chunks)), rep("code", length(spans))),
    stringsAsFactors = FALSE
  )
}

# The argument names a function accepts: its formals, plus the formals of
# whatever it passes `...` to; for an S3 generic, the union over its methods.
# NULL means any name may be right (the dots go to a function with dots of
# its own, or the call is rebuilt with match.call() and evaluated again, in
# the body or in the caller's frame handed to a helper).
.g08_accepted <- function(f, ns, fname = "") {
  formal_names <- names(formals(f))
  if (!"..." %in% formal_names) {
    return(formal_names)
  }
  body_calls <- all.names(body(f))
  if ("UseMethod" %in% body_calls) {
    table_s3 <- get(".__S3MethodsTable__.", envir = ns)
    method_names <- grep(paste0("^", gsub(".", "\\.", fname, fixed = TRUE), "\\."),
      ls(table_s3, all.names = TRUE), value = TRUE)
    out <- setdiff(formal_names, "...")
    for (m in method_names) {
      acc <- .g08_accepted(get(m, envir = table_s3), ns, m)
      if (is.null(acc)) {
        return(NULL)
      }
      out <- c(out, acc)
    }
    return(unique(out))
  }
  rebuilds_call <- any(c("match.call", "sys.call") %in% body_calls) &&
    any(c("eval", "eval.parent", "do.call", "parent.frame") %in% body_calls)
  if (rebuilds_call || any(c("...length", "...elt", "..1") %in% body_calls)) {
    return(NULL)
  }
  callees <- character(0)
  open_ended <- FALSE
  walk <- function(e) {
    if (is.call(e)) {
      parts <- as.list(e)
      has_dots <- any(vapply(parts[-1L], function(z) identical(z, quote(...)), logical(1L)))
      if (has_dots) {
        if (is.name(parts[[1L]])) {
          callees[[length(callees) + 1L]] <<- as.character(parts[[1L]])
        } else {
          open_ended <<- TRUE
        }
      }
      for (i in seq_along(parts)[-1L]) {
        if (!identical(parts[[i]], quote(expr = ))) walk(parts[[i]])
      }
    }
  }
  walk(body(f))
  if (open_ended) {
    return(NULL)
  }
  out <- setdiff(formal_names, "...")
  for (callee in unique(callees)) {
    g <- get0(callee, envir = ns, mode = "function")
    if (is.null(g) || is.primitive(g)) {
      return(NULL)
    }
    g_names <- names(formals(g))
    if ("..." %in% g_names) {
      return(NULL)
    }
    out <- c(out, g_names)
  }
  unique(out)
}

# Every call to a package function in one snippet that cannot be made. A
# snippet that is not R code (prose in \code{}, a path, a shell line) is
# skipped: only code that parses is a claim about the package's calls.
.g08_bad_calls <- function(snippet, kind, ns, exports) {
  exprs <- tryCatch(parse(text = snippet, keep.source = FALSE),
    error = function(e) NULL
  )
  if (is.null(exprs)) {
    return(character(0))
  }
  bad <- character(0)
  defined <- character(0)
  check <- function(e) {
    if (!is.call(e)) {
      return(invisible())
    }
    head <- e[[1L]]
    fname <- NA_character_
    internal_ok <- FALSE
    if (is.name(head)) {
      fname <- as.character(head)
    } else if (is.call(head) && as.character(head[[1L]]) %in% c("::", ":::") &&
      identical(as.character(head[[2L]]), "flexyBayes")) {
      fname <- as.character(head[[3L]])
      internal_ok <- identical(as.character(head[[1L]]), ":::")
      if (!internal_ok && !fname %in% exports) {
        bad[[length(bad) + 1L]] <<- paste0("flexyBayes::", fname, " is not exported")
      }
    }
    if (!is.na(fname) && fname %in% c("<-", "=") && is.name(e[[2L]])) {
      defined[[length(defined) + 1L]] <<- as.character(e[[2L]])
    }
    known <- !is.na(fname) && !fname %in% defined
    if (known && grepl("^(fb_|flexybayes)", fname) && !exists(fname, envir = ns, inherits = FALSE)) {
      bad[[length(bad) + 1L]] <<- paste0(fname, "() does not exist")
    }
    if (known && kind == "example" && is.name(head) && !fname %in% exports &&
      exists(fname, envir = ns, inherits = FALSE) &&
      !exists(fname, envir = baseenv()) && is.function(get(fname, envir = ns))) {
      bad[[length(bad) + 1L]] <<- paste0(fname, "() is internal, so an example cannot call it by name")
    }
    if (known && (fname %in% exports || internal_ok) && exists(fname, envir = ns, inherits = FALSE)) {
      f <- get(fname, envir = ns)
      if (is.function(f) && !is.primitive(f)) {
        accepted <- .g08_accepted(f, ns, fname)
        arg_names <- names(as.list(e))[-1L]
        arg_names <- arg_names[nzchar(arg_names)]
        if (!is.null(accepted)) {
          for (a in arg_names) {
            if (is.na(pmatch(a, accepted))) {
              bad[[length(bad) + 1L]] <<- paste0(fname, "(", a, " = ) reaches no argument")
            }
          }
        }
      }
    }
    parts <- as.list(e)[-1L]
    for (i in seq_along(parts)) {
      if (!identical(parts[[i]], quote(expr = ))) check(parts[[i]])
    }
  }
  for (e in exprs) check(e)
  unique(bad)
}

test_that("DC-08 guard: every documented call to a package function can be made", {
  root <- .g08_root()
  skip_if(is.na(root), "package source tree not available")
  ns <- asNamespace("flexyBayes")
  exports <- getNamespaceExports("flexyBayes")

  rd_files <- list.files(file.path(root, "man"), pattern = "\\.Rd$", full.names = TRUE)
  md_files <- c(
    file.path(root, c("README.md", "API_STABILITY.md", "inst/KNOWN_ISSUES.md")),
    list.files(file.path(root, "vignettes"), pattern = "\\.Rmd(\\.orig)?$", full.names = TRUE)
  )
  md_files <- md_files[file.exists(md_files)]
  expect_gt(length(rd_files), 0L)

  findings <- character(0)
  n_snippets <- 0L
  for (f in c(rd_files, md_files)) {
    snippets <- if (grepl("\\.Rd$", f)) .g08_rd_snippets(f) else .g08_md_snippets(f)
    n_snippets <- n_snippets + nrow(snippets)
    for (i in seq_len(nrow(snippets))) {
      hits <- .g08_bad_calls(snippets$text[i], snippets$kind[i], ns, exports)
      if (length(hits)) findings <- c(findings, paste0(basename(f), ": ", hits))
    }
  }
  findings <- unique(findings)

  # The sweep read the documentation it claims to read.
  expect_gt(n_snippets, 100L)
  expect_identical(findings, character(0))
})
