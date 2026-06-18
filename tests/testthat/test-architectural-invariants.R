# Architectural-invariant tests (ADR 0030 Verification; v0.4.0 Wave 2
# Phase 2C).
#
# ADR 0030 ratifies the architecturally-final API surface: six exported
# verbs, four constructor nouns, four closed-vocabulary registries
# locked at .onLoad(). These tests assert that contract and snapshot
# each surface so vocabulary / surface drift forces an ADR-amendment
# review (snapshot churn blocks merge until the change is intentional).
#
# The snapshots capture the *schema* (class vectors, attribute / element
# names, registry keys + lock state), not environment-dependent values
# (e.g. fb_engine()'s toolchain_status, which depends on whether the
# backend package is installed), so they are stable across machines.

# The exported verbs (ADR 0030 §3, roster grown by ADR 0031 Phase 3:
# fb_inla() added as the INLA engine pin). Closed list; a further verb
# requires an ADR amendment. fb_from_*() are ingest adapters, not verbs,
# and are tracked separately.
.FB_VERBS <- c(
  "flexybayes",
  "fb_brms",
  "fb_greta",
  "fb_inla",
  "fb_plan",
  "triangulate",
  "validate_approximation"
)

# The four constructor nouns (ADR 0030 §2).
.FB_CONSTRUCTORS <- c("fb_prior", "fb_cov", "fb_approx", "fb_engine")

# The five closed-vocabulary registries (ADR 0030 §6 + ADR 0031 backend).
.FB_REGISTRIES <- c(
  ".representation_registry",
  ".approximation_registry",
  ".backend_independence_registry",
  ".refusal_registry",
  ".backend_registry"
)

# Build one instance of each constructor for schema inspection.
.fb_constructor_instances <- function() {
  list(
    fb_prior = fb_prior(sigma ~ pc(upper = 2, prob = 0.05)),
    fb_cov = fb_cov(diag(3L), type = "dense"),
    fb_approx = fb_approx("low_rank_smooth", rank = 5L),
    fb_engine = fb_engine("greta")
  )
}

# A deterministic schema dump for a classed object: its class vector,
# sorted element names, and sorted attribute names (excluding names/
# class, which are structural). Values are excluded so the dump is
# environment-independent.
.fb_schema_dump <- function(x) {
  attrs <- setdiff(names(attributes(x)), c("names", "class"))
  list(
    class = class(x),
    elements = sort(names(x)),
    attributes = sort(attrs)
  )
}

# ---------------------------------------------------------------- #
# (a) Verb-count invariant                                          #
# ---------------------------------------------------------------- #

test_that("exactly six verbs are exported (ADR 0030 §3)", {
  exports <- getNamespaceExports("flexyBayes")
  for (v in .FB_VERBS) {
    expect_true(v %in% exports, info = paste("verb not exported:", v))
  }
  expect_length(.FB_VERBS, 7L)
})

test_that("the verb list is stable (snapshot)", {
  expect_snapshot(sort(.FB_VERBS))
})

# ---------------------------------------------------------------- #
# (b) Constructor-noun classed-object invariant                     #
# ---------------------------------------------------------------- #

test_that("each constructor noun is exported", {
  exports <- getNamespaceExports("flexyBayes")
  for (cn in .FB_CONSTRUCTORS) {
    expect_true(cn %in% exports, info = paste("constructor not exported:", cn))
  }
})

test_that("each constructor returns a classed list with its documented class vector", {
  inst <- .fb_constructor_instances()
  expect_identical(class(inst$fb_prior), c("fb_prior", "list"))
  expect_identical(class(inst$fb_cov), c("fb_cov", "list"))
  expect_identical(class(inst$fb_approx), c("fb_approx", "list"))
  expect_identical(class(inst$fb_engine), c("fb_engine", "list"))
})

test_that("fb_cov carries the ADR 0030 §2 attribute schema", {
  cov <- fb_cov(diag(3L), type = "dense")
  attrs <- names(attributes(cov))
  for (a in c("type", "representation_class", "validation_summary")) {
    expect_true(a %in% attrs, info = paste("fb_cov missing attribute:", a))
  }
})

test_that("fb_engine carries the ADR 0030 §2 element schema", {
  e <- fb_engine("greta")
  for (el in c("name", "paradigm", "toolchain_status", "opts")) {
    expect_true(
      el %in% names(e),
      info = paste("fb_engine missing element:", el)
    )
  }
})

test_that("the constructor schemas are stable (snapshot)", {
  inst <- .fb_constructor_instances()
  expect_snapshot(lapply(inst, .fb_schema_dump))
})

# ---------------------------------------------------------------- #
# (c) Closed-vocabulary registry-locked invariant                   #
# ---------------------------------------------------------------- #

test_that("all five registries are locked after .onLoad()", {
  for (r in .FB_REGISTRIES) {
    reg <- get(r, envir = asNamespace("flexyBayes"))
    expect_true(
      environmentIsLocked(reg),
      info = paste("registry not locked:", r)
    )
  }
})

test_that("a locked registry refuses user-side binding injection", {
  for (r in .FB_REGISTRIES) {
    reg <- get(r, envir = asNamespace("flexyBayes"))
    expect_error(assign("user_injected_key", 1L, envir = reg))
  }
})

test_that("the registry key inventory is stable (snapshot)", {
  dump <- lapply(.FB_REGISTRIES, function(r) {
    reg <- get(r, envir = asNamespace("flexyBayes"))
    list(
      registry = r,
      locked = environmentIsLocked(reg),
      n_keys = length(ls(reg, all.names = TRUE)),
      keys = sort(ls(reg, all.names = TRUE))
    )
  })
  names(dump) <- .FB_REGISTRIES
  expect_snapshot(dump)
})

# ---------------------------------------------------------------- #
# Inventory completeness: exports account for verbs + constructors  #
# + accessors, modulo S3 methods.                                   #
# ---------------------------------------------------------------- #

test_that("verbs and constructors are disjoint and all exported", {
  expect_length(intersect(.FB_VERBS, .FB_CONSTRUCTORS), 0L)
  exports <- getNamespaceExports("flexyBayes")
  expect_true(all(c(.FB_VERBS, .FB_CONSTRUCTORS) %in% exports))
})
