# =============================================================================
# Every stochastic example carries its own seed (review of 2026-09-11).
#
# The package's mechanical invariant checker marked seed discipline present
# because it found `set.seed()` calls somewhere in the sources. Five example
# blocks generated random data without one: fb_from_asreml, fb_complete_grid,
# fb_cov, fb_from_brms and fb_inla. Two clean runs of those examples did not
# produce identical inputs, so the stated reproducibility invariant was not
# met, an invariant checked in aggregate is not checked.
#
# This gate reads each `@examples` block on its own and asks whether a seed is
# set *before* the first draw in that same block. A seed in a neighbouring
# block does not count, which is the whole point.
# =============================================================================

.es_root <- function() {
  root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    mustWork = FALSE
  )
  if (dir.exists(file.path(root, "R"))) root else NA_character_
}

# The generators that make an example non-reproducible without a seed.
.es_rng_pattern <- function() {
  paste0(
    "\\b(rnorm|runif|rbinom|rpois|rgamma|rbeta|rexp|sample|rmultinom|",
    "rcauchy|rchisq|rlnorm|rweibull|rnbinom|mvrnorm)\\s*\\("
  )
}

# Split one file's roxygen into its @examples blocks.
.es_example_blocks <- function(lines) {
  starts <- grep("^#'\\s*@examples", lines)
  lapply(starts, function(i) {
    j <- i + 1L
    block <- character(0)
    while (j <= length(lines) && startsWith(lines[[j]], "#'")) {
      body <- sub("^#'\\s*", "", lines[[j]])
      if (grepl("^@[a-zA-Z]", body)) {
        break
      }
      block <- c(block, body)
      j <- j + 1L
    }
    list(line = i, body = block)
  })
}

test_that("each example block that draws random numbers seeds itself first", {
  root <- .es_root()
  if (is.na(root)) {
    testthat::skip("source tree not reachable (installed-package test run)")
  }
  files <- list.files(
    file.path(root, "R"),
    pattern = "\\.R$",
    full.names = TRUE
  )
  rng <- .es_rng_pattern()

  offenders <- character(0)
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    for (blk in .es_example_blocks(lines)) {
      text <- paste(blk$body, collapse = "\n")
      draw <- regexpr(rng, text, perl = TRUE)
      if (draw[[1]] == -1L) {
        next
      }
      seed <- regexpr("\\bset\\.seed\\s*\\(", text, perl = TRUE)
      seeded_first <- seed[[1]] != -1L && seed[[1]] < draw[[1]]
      if (!seeded_first) {
        offenders <- c(
          offenders,
          paste0(basename(f), ":", blk$line)
        )
      }
    }
  }
  expect_equal(offenders, character(0))
})
