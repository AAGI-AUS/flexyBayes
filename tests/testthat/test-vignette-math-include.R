# Every vignette must reference the shared MathJax sizing include ------------
#
# rmarkdown::html_vignette loads MathJax v2 at the foot of the page, and
# v2's TeX web fonts render visibly larger than the vignette body text --
# inline math towers over its sentence. The one-place fix is
# vignettes/_mathjax-scale.html (a MathJax.Hub.Config block that must run
# BEFORE the loader, hence in_header), referenced from every vignette's
# YAML. These tests keep that binding: the asset must exist, and every
# shipped page AND its .orig source must reference it, so a new vignette
# cannot silently reintroduce the oversized-math regression.

# The vignette sources live at vignettes/ in a source checkout and are not
# installed. Resolve them relative to the test directory; return NULL when
# the layout does not carry them (e.g. running from an installed package),
# so the tests skip rather than assert on an absent tree.
.vignette_source_dir <- function() {
  cand <- testthat::test_path("..", "..", "vignettes")
  if (dir.exists(cand)) {
    return(normalizePath(cand))
  }
  NULL
}

test_that("the shared MathJax sizing asset exists and configures a scale", {
  vig_dir <- .vignette_source_dir()
  skip_if(is.null(vig_dir), "vignette sources not present in this layout")
  asset <- file.path(vig_dir, "_mathjax-scale.html")
  expect_true(file.exists(asset))
  txt <- paste(readLines(asset, warn = FALSE), collapse = "\n")
  expect_match(txt, "x-mathjax-config", fixed = TRUE)
  expect_match(txt, "scale", fixed = TRUE)
})

test_that("every vignette page and .orig source references the include", {
  vig_dir <- .vignette_source_dir()
  skip_if(is.null(vig_dir), "vignette sources not present in this layout")
  pages <- list.files(
    vig_dir,
    pattern = "^flexyBayes-.*\\.Rmd(\\.orig)?$",
    full.names = TRUE
  )
  expect_gt(length(pages), 0L)
  missing <- pages[!vapply(
    pages,
    function(p) {
      any(grepl(
        "in_header: _mathjax-scale.html",
        readLines(p, warn = FALSE),
        fixed = TRUE
      ))
    },
    logical(1L)
  )]
  expect_identical(basename(missing), character(0))
})
