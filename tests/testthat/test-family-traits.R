# =============================================================================
# Per-engine family traits: one table, grounded on the engine.
#
# Two hand-maintained lists in this package disagreed about which families
# carry a residual `sigma`. `.brms_family_has_sigma()` listed gamma and
# beta; the heteroscedastic-residual gate in emit_brms.R did not. The wrong
# list governed the prior emit, so every prior route on gamma and beta --
# including the plain default with no prior argument at all -- sent brms a
# prior for a parameter the model does not have, and brms refused the fit.
# Two of the six advertised response families were unreachable on the Stan
# backend, and the capability matrix said `fits`.
#
# Three assertions hold the fix. The table agrees with brms's own
# declaration of each family's distributional parameters (the engine is the
# oracle, not a second list in our tree). Both call sites read the table.
# And no second copy of the roster exists in the package sources -- the
# duplicate is the defect, so the test looks for duplicates rather than for
# the symptom.
# =============================================================================

# The response families the entry allowlist admits (R/utils.R
# .resolve_family). The trait table has to answer for each of them.
.ft_admitted_families <- c(
  "gaussian",
  "binomial",
  "binary",
  "poisson",
  "negative_binomial",
  "negbinom",
  "gamma",
  "beta"
)

.ft_source_dir <- function() {
  candidate <- testthat::test_path("..", "..", "R")
  if (dir.exists(candidate)) {
    return(normalizePath(candidate))
  }
  NA_character_
}

# The code tokens of one R file, comments removed, so a roster written in a
# comment is not mistaken for a second implementation of the fact.
.ft_code_text <- function(file) {
  pd <- utils::getParseData(parse(file, keep.source = TRUE))
  pd <- pd[pd$terminal & pd$token != "COMMENT", , drop = FALSE]
  pd <- pd[order(pd$line1, pd$col1), , drop = FALSE]
  paste(pd$text, collapse = " ")
}

# --- the engine is the oracle -------------------------------------------- #

test_that("the residual-sigma roster matches brms's own dpars declaration", {
  skip_if_not_installed("brms")
  roster <- flexyBayes:::.fb_brms_families_with_sigma()
  for (fam in .ft_admitted_families) {
    brms_family <- flexyBayes:::.fb_family_to_brms(fam, link = NULL)
    dpars <- brms::brmsfamily(brms_family$family)$dpars
    expect_identical(
      fam %in% roster,
      "sigma" %in% dpars,
      info = paste0(
        fam,
        " -> ",
        brms_family$family,
        " (dpars: ",
        paste(dpars, collapse = ", "),
        ")"
      )
    )
  }
})

test_that("lognormal, which the emit supports, is on the roster for the same reason", {
  skip_if_not_installed("brms")
  expect_true("lognormal" %in% flexyBayes:::.fb_brms_families_with_sigma())
  expect_true("sigma" %in% brms::brmsfamily("lognormal")$dpars)
})

# --- both call sites read the table -------------------------------------- #

test_that("the prior emit predicate is the table, not a copy of it", {
  roster <- flexyBayes:::.fb_brms_families_with_sigma()
  for (fam in c(.ft_admitted_families, "lognormal", "Gamma", "Beta")) {
    expect_identical(
      flexyBayes:::.brms_family_has_sigma(fam),
      tolower(fam) %in% roster,
      info = fam
    )
  }
  # A NULL family is the gaussian default and does carry a sigma.
  expect_true(flexyBayes:::.brms_family_has_sigma(NULL))
})

test_that("a sigma prior is dropped exactly for the families with no sigma", {
  fb <- list(
    response = "y",
    intercept = TRUE,
    fixed_terms = list(list(type = "simple", var = "x")),
    random_terms = list(list(type = "simple", var = "g"))
  )
  prior <- fb_prior(
    sigma ~ half_normal(scale = 2),
    sd(group = "g") ~ half_normal(scale = 1.5)
  )
  for (fam in .ft_admitted_families) {
    fb$family <- fam
    specs <- flexyBayes:::.priors_to_brms_specs(prior, fb)
    classes <- vapply(specs, function(s) s$class, character(1))
    expect_identical(
      "sigma" %in% classes,
      fam %in% flexyBayes:::.fb_brms_families_with_sigma(),
      info = fam
    )
    # The variance-component row survives on every family.
    expect_true("sd" %in% classes, info = fam)
  }
})

test_that("the legacy scalar route is family-aware on the same table", {
  fb <- list(
    response = "y",
    intercept = TRUE,
    fixed_terms = list(list(type = "simple", var = "x")),
    random_terms = list(list(type = "simple", var = "g"))
  )
  for (fam in .ft_admitted_families) {
    fb$family <- fam
    specs <- flexyBayes:::.priors_to_brms_specs(
      NULL,
      fb,
      prior_fixed_sd = 100,
      prior_vc_sd = 3
    )
    classes <- vapply(specs, function(s) s$class, character(1))
    expect_identical(
      "sigma" %in% classes,
      fam %in% flexyBayes:::.fb_brms_families_with_sigma(),
      info = fam
    )
  }
})

# --- no second hand-maintained copy -------------------------------------- #

test_that("the residual-sigma roster is written down in exactly one file", {
  src_dir <- .ft_source_dir()
  skip_if(is.na(src_dir), "not running from the source tree")
  roster <- flexyBayes:::.fb_brms_families_with_sigma()
  files <- list.files(src_dir, pattern = "[.][Rr]$", full.names = TRUE)
  carriers <- character(0)
  for (f in files) {
    text <- .ft_code_text(f)
    vectors <- regmatches(text, gregexpr("c \\([^()]*\\)", text))[[1L]]
    for (v in vectors) {
      elements <- regmatches(v, gregexpr("\"[^\"]*\"", v))[[1L]]
      elements <- gsub("\"", "", elements, fixed = TRUE)
      if (length(elements) && setequal(elements, roster)) {
        carriers <- c(carriers, basename(f))
        break
      }
    }
  }
  expect_identical(
    carriers,
    "family_traits.R",
    label = paste0(
      "files carrying the residual-sigma roster as a literal: ",
      paste(carriers, collapse = ", ")
    )
  )
})

# The two files that used to carry the disagreeing copies both spelled
# `skew_normal`, a family the entry allowlist does not admit and which
# appears in no other context in this package. Its reappearance in either
# file is the signature of a third hand-maintained roster.
test_that("neither former carrier lists response families for the sigma decision", {
  src_dir <- .ft_source_dir()
  skip_if(is.na(src_dir), "not running from the source tree")
  for (f in c("priors_to_brms.R", "emit_brms.R")) {
    path <- file.path(src_dir, f)
    skip_if(!file.exists(path), paste(f, "not present"))
    expect_false(
      grepl("skew_normal", .ft_code_text(path), fixed = TRUE),
      label = paste0(f, " carries a residual-scale family roster again")
    )
  }
})
