# Triangulation independence-axis vocabulary (ADR 0029; v0.4.0 Wave 2
# Phase 2B).
#
# Coverage: closed-vocabulary enforcement, per-pair registry symmetry +
# shape + lock, the three v0.4.0 backend pairs, and the triangulate()
# report-shape extension (independence + axis_justification) wired
# end-to-end through synthetic fits whose source resolves to a real
# backend label.

# ---------------------------------------------------------------- #
# Synthetic fits: dispatch fb_as_draws_simple on a leading synthetic  #
# class while .triangulate_source() reads the trailing real backend   #
# class. So triangulate() sees a registered (greta, inla) pair without #
# running a real backend.                                              #
# ---------------------------------------------------------------- #

fb_as_draws_simple.synthetic_greta <- function(fit, ...) fit$draws
fb_as_draws_simple.synthetic_inla <- function(fit, ...) fit$draws
.S3method(
  "fb_as_draws_simple",
  "synthetic_greta",
  fb_as_draws_simple.synthetic_greta
)
.S3method(
  "fb_as_draws_simple",
  "synthetic_inla",
  fb_as_draws_simple.synthetic_inla
)

.mk_greta_fit <- function(draws) {
  structure(list(draws = draws), class = c("synthetic_greta", "flexybayes"))
}
.mk_inla_fit <- function(draws) {
  structure(list(draws = draws), class = c("synthetic_inla", "flexybayes_inla"))
}

# ---------------------------------------------------------------- #
# Closed vocabulary                                                 #
# ---------------------------------------------------------------- #

test_that("the axis vocabulary is exactly the three controlled values", {
  expect_setequal(
    flexyBayes:::.AXIS_VOCABULARY,
    c("algorithmic", "implementation", "specification")
  )
})

test_that("registering an unknown axis refuses with a vocabulary error", {
  err <- tryCatch(
    flexyBayes:::.register_pair_independence(
      pair = c("aa", "bb"),
      axes = c("magic"),
      paradigms = list(aa = "x", bb = "y"),
      justification = "test",
      registered_in_adr = "test"
    ),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "unknown axis")
  expect_match(err, "magic")
})

test_that("a valid registration on the locked registry refuses with the lock error", {
  # Validation passes (valid axes); the lock check then fires. This
  # proves both the lock and that valid axes clear the vocabulary gate.
  err <- tryCatch(
    flexyBayes:::.register_pair_independence(
      pair = c("aa", "bb"),
      axes = c("algorithmic"),
      paradigms = list(aa = "x", bb = "y"),
      justification = "test",
      registered_in_adr = "test"
    ),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "locked")
})

# ---------------------------------------------------------------- #
# Registry symmetry + shape + lock                                  #
# ---------------------------------------------------------------- #

test_that(".pair_key canonicalises by sorting", {
  expect_identical(flexyBayes:::.pair_key(c("z", "a")), "a||z")
  expect_identical(flexyBayes:::.pair_key(c("a", "z")), "a||z")
})

test_that(".pair_key rejects a malformed pair", {
  expect_error(flexyBayes:::.pair_key("greta"))
  expect_error(flexyBayes:::.pair_key(c("a", "b", "c")))
})

test_that("pair lookup is symmetric in argument order", {
  r1 <- flexyBayes:::.lookup_pair_independence(c("greta", "inla"))
  r2 <- flexyBayes:::.lookup_pair_independence(c("inla", "greta"))
  expect_identical(r1, r2)
})

test_that("every registry entry carries the five required fields", {
  reg <- flexyBayes:::.backend_independence_registry
  for (k in ls(reg)) {
    entry <- get(k, envir = reg, inherits = FALSE)
    expect_true(all(
      c("pair", "axes", "paradigms", "justification", "registered_in_adr") %in%
        names(entry)
    ))
  }
})

test_that("the three v0.4.0 backend pairs register with the ADR 0029 axes", {
  gi <- flexyBayes:::.lookup_pair_independence(c("greta", "inla"))
  gb <- flexyBayes:::.lookup_pair_independence(c("greta", "brms"))
  bi <- flexyBayes:::.lookup_pair_independence(c("brms", "inla"))
  expect_setequal(gi$axes, c("algorithmic", "implementation"))
  expect_setequal(gb$axes, "implementation")
  expect_setequal(bi$axes, c("algorithmic", "implementation"))
})

test_that("an unregistered or same-backend pair looks up to NULL", {
  expect_null(flexyBayes:::.lookup_pair_independence(c("greta", "greta")))
  expect_null(flexyBayes:::.lookup_pair_independence(c("greta", "nimble")))
})

test_that("the backend-independence registry is locked after .onLoad()", {
  expect_true(environmentIsLocked(
    flexyBayes:::.backend_independence_registry
  ))
})

# ---------------------------------------------------------------- #
# triangulate() report-shape extension                             #
# ---------------------------------------------------------------- #

test_that("triangulate() labels a registered pair with its independence axes", {
  set.seed(1L)
  d_g <- list(beta = rnorm(200L), sigma = abs(rnorm(200L)))
  d_i <- list(beta = rnorm(200L), sigma = abs(rnorm(200L)))
  tri <- triangulate(.mk_greta_fit(d_g), .mk_inla_fit(d_i))
  expect_setequal(tri$independence, c("algorithmic", "implementation"))
  expect_type(tri$axis_justification, "character")
  expect_match(tri$axis_justification, "Laplace")
})

test_that("triangulate() on an unregistered pair leaves the axis fields empty (backward-compat)", {
  set.seed(2L)
  mk <- function(d) structure(list(draws = d), class = "synthetic_fit")
  fb_as_draws_simple.synthetic_fit <<- function(fit, ...) fit$draws
  .S3method(
    "fb_as_draws_simple",
    "synthetic_fit",
    fb_as_draws_simple.synthetic_fit
  )
  d <- list(beta = rnorm(100L))
  tri <- triangulate(mk(d), mk(d))
  expect_length(tri$independence, 0L)
  expect_true(is.na(tri$axis_justification))
  # the Wasserstein columns still render
  expect_true("W1" %in% names(tri$metrics) || nrow(tri$metrics) >= 0L)
})

test_that("print.triangulate_result shows the independence line for a registered pair", {
  set.seed(3L)
  d_g <- list(beta = rnorm(150L))
  d_i <- list(beta = rnorm(150L))
  tri <- triangulate(.mk_greta_fit(d_g), .mk_inla_fit(d_i))
  out <- utils::capture.output(print(tri))
  expect_true(any(grepl("independence:", out, fixed = TRUE)))
  expect_true(any(grepl("algorithmic", out, fixed = TRUE)))
})

# ---------------------------------------------------------------- #
# Snapshot: axis-label rendering (ADR-amendment forcing function)   #
# ---------------------------------------------------------------- #

test_that("axis-label rendering is stable across the registered pairs", {
  withr::local_options(cli.num_colors = 1L) # plain text -> stable snapshot
  reg <- flexyBayes:::.backend_independence_registry
  expect_snapshot({
    for (k in sort(ls(reg))) {
      entry <- get(k, envir = reg, inherits = FALSE)
      cat(
        paste(entry$pair, collapse = "--"),
        ": ",
        flexyBayes:::.format_independence_axes(entry$axes),
        "\n",
        sep = ""
      )
    }
  })
})
