# =============================================================================
# .fb_select_gate() reads the parser's type names, not a second hand-written
# list of them.
#
# Contract: review/units/U-SEL-GATE-1_2026-09-16/README.md (U-FOLD-1 finding
# 4). The gate's `us` / `fa` check used to name the strings "us", "fa",
# "unstructured" and "factor_analytic" -- and .walk() / .classify_interaction()
# (R/parse_formula.R) have never produced the last two, and pair a bare
# us()/fa() with a `:` inner factor into "us_gxe" / "fa_gxe", which the old
# check did not name at all. A base carrying `us(env):site` or
# `fa(env, k):site` -- the only forms of either structure that ever reach a
# fit on either engine -- passed the gate unnoticed.
#
# No engine fit runs anywhere in this file: .fb_select_gate() decides from
# base$fb$family, base$fb$link and base$fb$random_terms alone, so every case
# below is a real fb_from_asreml() parse (or, for the brms-grammar case, a
# real fb_from_brms() parse) wrapped in a stand-in `flexybayes`-classed list
# -- the same stub shape "a correlated random slope cannot reach selection
# through an INLA base" (test-fb-select.R) already uses for the gate's other
# engine-specific check.
# =============================================================================

fb_sel_gate_types_data <- function() {
  set.seed(1L)
  n <- 24L
  data.frame(
    y = stats::rnorm(n),
    x1 = stats::rnorm(n),
    env = factor(rep(paste0("E", seq_len(4L)), each = 6L)),
    site = factor(rep(paste0("S", seq_len(6L)), times = 4L)),
    geno = factor(rep(paste0("G", seq_len(6L)), times = 4L)),
    t = factor(rep(seq_len(6L), times = 4L))
  )
}

fb_sel_gate_asreml_base <- function(random_formula, data) {
  fb <- fb_from_asreml(
    fixed = y ~ x1, random = random_formula, residual = NULL, data = data,
    family = "gaussian", link = "identity"
  )
  structure(list(fb = fb, data = data), class = "flexybayes")
}


# --- the enumeration ------------------------------------------------------ #

test_that("the gate's verdict on every ASReml-grammar random-term type matches the parser's own name for it", {
  d <- fb_sel_gate_types_data()

  # One representative formula per type .walk() / .classify_interaction()
  # can produce (28 distinct type strings). "admit" / "refuse" is the
  # observed base behaviour this unit fixes or, for every type outside its
  # scope, leaves unchanged -- both are asserted so a future change to any
  # of them is visible here, not only to the two this unit targets.
  cases <- list(
    simple               = list(f = ~site,                 refuse = FALSE),
    units                = list(f = ~units,                refuse = FALSE),
    vm                   = list(f = ~ vm(site),             refuse = FALSE),
    ide                  = list(f = ~ ide(site),            refuse = FALSE),
    ped                  = list(f = ~ ped(site),            refuse = FALSE),
    fa                   = list(f = ~ fa(env, 1),           refuse = TRUE),
    us                   = list(f = ~ us(env),              refuse = TRUE),
    at                   = list(f = ~ at(env),              refuse = FALSE),
    corh                 = list(f = ~ corh(env),            refuse = FALSE),
    ar1                  = list(f = ~ ar1(t),               refuse = FALSE),
    ar2                  = list(f = ~ ar2(t),               refuse = FALSE),
    cor                  = list(f = ~ cor(t),               refuse = FALSE),
    str                  = list(f = ~ str(~x1),             refuse = FALSE),
    spline               = list(f = ~ spl(x1),              refuse = FALSE),
    smooth_mgcv          = list(f = ~ s(x1),                refuse = FALSE),
    at_units             = list(f = ~ at(env):units,        refuse = FALSE),
    continuous           = list(f = ~ lin(x1),              refuse = FALSE),
    polynomial           = list(f = ~ pol(x1, 2),           refuse = FALSE),
    fa_gxe               = list(f = ~ fa(env, 1):id(geno),  refuse = TRUE),
    us_gxe               = list(f = ~ us(env):id(geno),     refuse = TRUE),
    corh_gxe             = list(f = ~ corh(env):id(geno),   refuse = FALSE),
    at_simple            = list(f = ~ diag(env):id(geno),   refuse = FALSE),
    ar1_spatial          = list(f = ~ ar1(t):ar1(site),     refuse = FALSE),
    at_ar1_partial       = list(f = ~ at(env):ar1(t),       refuse = FALSE),
    vm_gxe               = list(f = ~ vm(site):id(env),     refuse = FALSE),
    nested               = list(f = ~ site:geno,            refuse = FALSE),
    combo                = list(f = ~ site:geno:env,        refuse = FALSE),
    interaction_generic  = list(f = ~ spl(x1):site,         refuse = FALSE)
  )
  expect_length(cases, 28L)

  for (nm in names(cases)) {
    base <- fb_sel_gate_asreml_base(cases[[nm]]$f, d)
    parsed_types <- vapply(
      base$fb$random_terms, function(t) t$type %||% NA_character_, character(1L)
    )
    if (isTRUE(cases[[nm]]$refuse)) {
      expect_error(
        flexyBayes:::.fb_select_gate(base),
        class = "flexybayes_refusal_selection_base_not_allowlisted",
        info = nm
      )
    } else {
      expect_true(flexyBayes:::.fb_select_gate(base), info = nm)
    }
    # The type the gate acted on is the type the real parser wrote for this
    # formula, named by hand only in this test's own table above -- never a
    # second copy inside the gate.
    expect_true(nm %in% parsed_types, info = nm)
  }
})


# --- the named defect, isolated ------------------------------------------- #

test_that("a base carrying us(env):site refuses before any fold runs", {
  # The brief's own reproduction (U-FOLD-1 finding 4, probe_01): the gate
  # returned TRUE on this base before this unit.
  d <- fb_sel_gate_types_data()
  base <- fb_sel_gate_asreml_base(~ us(env):id(site), d)
  expect_error(
    flexyBayes:::.fb_select_gate(base),
    class = "flexybayes_refusal_selection_base_not_allowlisted"
  )
})

test_that("a base carrying fa(env, k):site refuses before any fold runs", {
  d <- fb_sel_gate_types_data()
  base <- fb_sel_gate_asreml_base(~ fa(env, 2):id(site), d)
  expect_error(
    flexyBayes:::.fb_select_gate(base),
    class = "flexybayes_refusal_selection_base_not_allowlisted"
  )
})

test_that("a base mixing us(env):site with an admitted term still refuses", {
  # The bad type does not have to be the base's only random term to be
  # caught.
  d <- fb_sel_gate_types_data()
  base <- fb_sel_gate_asreml_base(~ us(env):id(site) + geno, d)
  expect_error(
    flexyBayes:::.fb_select_gate(base),
    class = "flexybayes_refusal_selection_base_not_allowlisted"
  )
})


# --- the refusal names the measurement, per family ------------------------ #

test_that("the us() refusal names the engine it fits on and the measurement that keeps it off the allowlist", {
  d <- fb_sel_gate_types_data()
  base <- fb_sel_gate_asreml_base(~ us(env):id(site), d)
  expect_error(
    flexyBayes:::.fb_select_gate(base), "brms only", fixed = TRUE
  )
  expect_error(
    flexyBayes:::.fb_select_gate(base), "0.142", fixed = TRUE
  )
  expect_error(
    flexyBayes:::.fb_select_gate(base), "0.930", fixed = TRUE
  )
})

test_that("the fa() refusal names the engine it fits on and that it is experimental", {
  d <- fb_sel_gate_types_data()
  base <- fb_sel_gate_asreml_base(~ fa(env, 1):id(site), d)
  expect_error(
    flexyBayes:::.fb_select_gate(base), "INLA only", fixed = TRUE
  )
  expect_error(
    flexyBayes:::.fb_select_gate(base), "experimental", fixed = TRUE
  )
})


# --- the brms-grammar namespace never carries either type ----------------- #

test_that("a brms-grammar base never carries us_gxe or fa_gxe, so it is not caught by name here", {
  # fb_from_brms() lowers `(x | g)` to simple_slope_cor and has no route to
  # "us_gxe" / "fa_gxe" at all (R/fb_from_brms.R); the ASReml spellings
  # us()/fa() are not part of the brms grammar. Recorded so a future change
  # to fb_from_brms() that DID start writing one of these two strings would
  # be caught by the enumeration block above, not silently admitted here.
  fb <- fb_from_brms(y ~ x1 + (x1 | site), data = fb_sel_gate_types_data(),
                      family = "gaussian")
  types <- vapply(fb$random_terms, function(t) t$type %||% NA_character_, character(1L))
  expect_false(any(types %in% c("us_gxe", "fa_gxe")))
})
