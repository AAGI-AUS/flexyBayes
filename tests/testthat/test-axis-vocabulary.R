# Backend-independence axis vocabulary (ADR 0029; v0.4.0 Wave 2
# Phase 2B).
#
# Coverage: closed-vocabulary enforcement, per-pair registry symmetry +
# shape + lock, and the one remaining backend pair (0.9.3 withdraws two
# of the original three v0.4.0 pairs along with the engine they named --
# see NEWS.md).

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
  expect_error(flexyBayes:::.pair_key("solo"))
  expect_error(flexyBayes:::.pair_key(c("a", "b", "c")))
})

test_that("pair lookup is symmetric in argument order", {
  r1 <- flexyBayes:::.lookup_pair_independence(c("brms", "inla"))
  r2 <- flexyBayes:::.lookup_pair_independence(c("inla", "brms"))
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

test_that("the one remaining backend pair registers with the ADR 0029 axes", {
  # v0.4.0 registered three pairs: the withdrawn native engine paired
  # with inla, that same engine paired with brms, and (brms, inla).
  # 0.9.3 withdraws that engine entirely (see NEWS.md); only
  # (brms, inla) survives.
  bi <- flexyBayes:::.lookup_pair_independence(c("brms", "inla"))
  expect_setequal(bi$axes, c("algorithmic", "implementation"))
  expect_length(ls(flexyBayes:::.backend_independence_registry), 1L)
})

test_that("a same-backend pair looks up to NULL", {
  expect_null(flexyBayes:::.lookup_pair_independence(c("inla", "inla")))
})

test_that("the backend-independence registry is locked after .onLoad()", {
  expect_true(environmentIsLocked(
    flexyBayes:::.backend_independence_registry
  ))
})

# ---------------------------------------------------------------- #
# Snapshot: registered pairs and axes (ADR-amendment forcing function) #
# ---------------------------------------------------------------- #

test_that("the registered pairs and their axes are stable", {
  reg <- flexyBayes:::.backend_independence_registry
  expect_snapshot({
    for (k in sort(ls(reg))) {
      entry <- get(k, envir = reg, inherits = FALSE)
      cat(
        paste(entry$pair, collapse = "--"),
        ": ",
        paste(entry$axes, collapse = " + "),
        "\n",
        sep = ""
      )
    }
  })
})
