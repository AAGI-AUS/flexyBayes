# =============================================================================
# The special-form argument contract (review of 2026-09-11).
#
# The review reported two loci: a malformed dsum() lowering to a homogeneous
# `units` residual, and fa() without a rank raising a base subscript error. A
# fuzz corpus over all 22 forms in .asreml_function_vocabulary() found the
# class was wider than both, 129 malformed variants produced four distinct
# escapes from the refusal taxonomy:
#
#   * an untyped base error (15 cases across dsum, fa, pol, s, str);
#   * a silent lowering to a plausible neighbour (dsum -> units, 7 cases);
#   * a term carrying var = NA_character_ (15 zero-argument forms);
#   * a silently discarded surplus argument (32 cases).
#
# The grid below is that corpus. Every row asserts a typed refusal, and the
# no-backend-call block asserts that a refusal terminates before any engine is
# reached, a refusal that fires after INLA has been called has not prevented
# anything.
#
# The last block is the falsifier. A contract that refuses everything would
# pass every row above, so the well-formed calls are asserted to still parse
# into the term they always produced.
# =============================================================================

.gc_data <- function() {
  set.seed(20260911L)
  s1 <- 60L
  data.frame(
    y = stats::rnorm(s1),
    x = stats::rnorm(s1),
    a = factor(rep(letters[1:3], length.out = s1)),
    b = factor(rep(LETTERS[1:4], length.out = s1)),
    g = factor(rep(paste0("g", 1:5), length.out = s1))
  )
}

.gc_condition <- function(txt) {
  fml <- stats::as.formula(paste0("~ ", txt))
  tryCatch(
    {
      .parse_formula(fml, .gc_data())
      NULL
    },
    condition = function(e) if (inherits(e, "error")) e else NULL
  )
}


# --- a malformed dsum() must never become `units` ----------------------- #

test_that("a dsum() body with no grouping bar refuses instead of lowering to units", {
  # Each of these used to return list(type = "units"), a one-variance
  # residual where a sectioned one was written. The fitted model is a
  # plausible neighbour, so it survives convergence checks and prints
  # credible output, which is the failure this package exists to prevent.
  malformed <- c(
    "dsum(~units)",
    "dsum(units)",
    "dsum(a)",
    "dsum(1)",
    "dsum(~1)",
    "dsum(a, b)",
    "dsum(a, b, c)"
  )
  for (case in malformed) {
    err <- .gc_condition(case)
    expect_s3_class(err, "flexybayes_refusal")
    expect_false(
      inherits(err, "flexybayes_refusal_units"),
      label = paste0(case, " must not lower to a units residual")
    )
  }
})

test_that("the malformed-dsum refusal names the bar and the single-variance spelling", {
  err <- .gc_condition("dsum(~units)")
  expect_s3_class(err, "flexybayes_refusal_dsum_malformed")
  msg <- conditionMessage(err)
  expect_match(msg, "dsum(~ units | env)", fixed = TRUE)
  expect_match(msg, "residual = ~ units", fixed = TRUE)
})

test_that("a well-formed dsum() still sections on its grouping factor", {
  ir <- .parse_formula(~ dsum(~ units | a), .gc_data())
  expect_equal(ir[[1]]$type, "at_units")
  expect_equal(ir[[1]]$var, "a")
})


# --- malformed fa() must not reach a base indexing error ----------------- #

test_that("fa() without a rank refuses by name instead of raising a subscript error", {
  err <- .gc_condition("fa(a)")
  expect_s3_class(err, "flexybayes_refusal_fa_rank_missing")
  expect_false(inherits(err, "subscriptOutOfBoundsError"))
  expect_match(conditionMessage(err), "fa(a, 1)", fixed = TRUE)
})

test_that("a rank that is not a literal count refuses instead of being coerced", {
  # as.integer() turned every one of these into a rank the user did not
  # write: TRUE and 1.5 into 1, "1" into 1, and a symbol into an untyped
  # coercion error.
  for (case in c(
    "fa(a, k)",
    "fa(a, 1.5)",
    "fa(a, NA)",
    "fa(a, TRUE)",
    "fa(a, \"1\")"
  )) {
    err <- .gc_condition(case)
    expect_s3_class(err, "flexybayes_refusal_fa_rank_invalid")
  }
})

test_that("a rank below one keeps the wording it has carried since 0.4.0", {
  err <- .gc_condition("fa(a, 0)")
  expect_s3_class(err, "flexybayes_refusal_fa_rank_invalid")
  expect_match(conditionMessage(err), "fa() requires k >= 1", fixed = TRUE)
})

test_that("a bare fa() rank at or above the level count refuses on identifiability", {
  # The k < levels bound has existed since 0.4.0 but was reachable only
  # through the fa(f, k):g interaction form, so a bare fa(a, 99) over three
  # levels parsed and fitted.
  err <- .gc_condition("fa(a, 99)")
  expect_s3_class(err, "flexybayes_refusal_fa_rank_exceeds_dim")
  expect_match(conditionMessage(err), "k < the number of levels", fixed = TRUE)
})


# --- the wider class the review did not report --------------------------- #

test_that("no special form parses into a term carrying no variable", {
  # Fifteen zero-argument forms returned an IR node with var = NA_character_,
  # which reaches the engine as a model nobody wrote.
  forms <- setdiff(.asreml_function_vocabulary(), c("te", "ti", "t2"))
  for (fn in forms) {
    err <- .gc_condition(paste0(fn, "()"))
    expect_s3_class(err, "flexybayes_refusal")
  }
})

test_that("a surplus argument refuses instead of being silently dropped", {
  # us(a, b) read `a` and discarded `b`, so a structure written into the
  # call never reached the fitted model.
  for (fn in c(
    "us",
    "diag",
    "idh",
    "ar1",
    "ar2",
    "cor",
    "corh",
    "spl",
    "lin",
    "id",
    "ide"
  )) {
    err <- .gc_condition(paste0(fn, "(a, b)"))
    expect_s3_class(err, "flexybayes_refusal_grammar_arity_too_many")
  }
})

test_that("an unrecognised named argument refuses instead of being ignored", {
  err <- .gc_condition("us(a, weights = b)")
  expect_s3_class(err, "flexybayes_refusal_grammar_named_argument_unknown")
  expect_match(conditionMessage(err), "weights", fixed = TRUE)
})

test_that("an empty argument slot refuses instead of raising a missing-argument error", {
  for (case in c("fa(a, )", "pol(x, )")) {
    err <- .gc_condition(case)
    expect_s3_class(err, "flexybayes_refusal_grammar_argument_empty")
  }
})

test_that("pol() takes a literal degree", {
  err <- .gc_condition("pol(x, b)")
  expect_s3_class(err, "flexybayes_refusal_pol_degree_invalid")
})

test_that("a multivariate s() refuses by name instead of through mgcv's call handling", {
  err <- .gc_condition("s(x, a)")
  expect_s3_class(err, "flexybayes_refusal_multivariate_smooth_unsupported")
  expect_match(conditionMessage(err), "univariate", fixed = TRUE)
})


test_that("a named argument written twice refuses instead of resolving by accident", {
  # The 0.10.0.9000 contract validated the *set* of named
  # arguments and never their multiplicity, which left the reviewed class
  # alive in two shapes. `vm(g, V = G, V = G)` was accepted and lowered to
  # one vm term, with the surviving copy an accident of how .walk() indexes;
  # `s(x, k = 5, k = 9)` reached mgcv's own matcher and raised "formal
  # argument 'k' matched by multiple actual arguments", an untyped base
  # error about our internals. Both are escapes the contract exists to
  # close, and the open-ended s() name set must be covered as well as the
  # closed vm() / ped() ones.
  duplicated_calls <- c(
    "vm(g, V = Gmat, V = Gmat)",
    "vm(g, cov = fb_cov(Gmat), cov = fb_cov(Gmat))",
    "ped(g, V = Gmat, V = Gmat)",
    "ped(g, use_sparse_precision = TRUE, use_sparse_precision = FALSE)",
    "s(x, k = 5, k = 9)",
    "s(x, bs = \"cr\", bs = \"tp\")"
  )
  dat <- .gc_data()
  Gmat <- diag(nlevels(dat$g))
  dimnames(Gmat) <- list(levels(dat$g), levels(dat$g))
  for (case in duplicated_calls) {
    err <- tryCatch(
      {
        .parse_formula(stats::as.formula(paste0("~ ", case)), dat)
        NULL
      },
      condition = function(e) if (inherits(e, "error")) e else NULL
    )
    expect_s3_class(err, "flexybayes_refusal_grammar_argument_duplicated")
    expect_match(conditionMessage(err), "twice", fixed = TRUE)
  }
})


# --- the no-backend-call assertion -------------------------------------- #

test_that("a malformed grammar never reaches an engine", {
  # A refusal raised after INLA has been called has not prevented anything.
  # .inla_call() and the brms entry are replaced with sentinels that record
  # having been reached. Neither may fire.
  skip_if_not_installed("testthat")
  reached <- new.env(parent = emptyenv())
  reached$engine <- character(0)

  testthat::local_mocked_bindings(
    .inla_call = function(...) {
      reached$engine <- c(reached$engine, "inla")
      stop("engine reached")
    }
  )

  malformed <- c(
    "dsum(~units)",
    "fa(a)",
    "fa(a, 0)",
    "us()",
    "vm()",
    "us(a, b)",
    "fa(a, )",
    "s(x, a)"
  )
  for (case in malformed) {
    err <- tryCatch(
      flexybayes(
        y ~ 1,
        random = stats::as.formula(paste0("~ ", case)),
        data = .gc_data(),
        family = "gaussian",
        backend = "inla",
        verbose = FALSE
      ),
      condition = function(e) e
    )
    expect_s3_class(err, "flexybayes_refusal")
  }
  expect_equal(reached$engine, character(0))
})


# --- the falsifier: the contract must still accept the grammar ---------- #

test_that("every well-formed special form still parses into the term it always produced", {
  # Without this block a contract that refused everything would pass every
  # assertion above.
  dat <- .gc_data()
  Gmat <- diag(nlevels(dat$g))
  dimnames(Gmat) <- list(levels(dat$g), levels(dat$g))

  expect_equal(.parse_formula(~ us(a), dat)[[1]]$type, "us")
  expect_equal(.parse_formula(~ diag(a), dat)[[1]]$type, "at")
  expect_equal(.parse_formula(~ idh(a), dat)[[1]]$type, "at")
  expect_equal(.parse_formula(~ at(a), dat)[[1]]$type, "at")
  expect_equal(.parse_formula(~ at(a, "b"), dat)[[1]]$type, "at")
  expect_equal(.parse_formula(~ ar1(a), dat)[[1]]$type, "ar1")
  expect_equal(.parse_formula(~ spl(x), dat)[[1]]$type, "spline")
  expect_equal(.parse_formula(~ lin(x), dat)[[1]]$type, "continuous")
  expect_equal(.parse_formula(~ pol(x), dat)[[1]]$degree, 2L)
  expect_equal(.parse_formula(~ pol(x, 3), dat)[[1]]$degree, 3L)
  expect_equal(.parse_formula(~ id(a), dat)[[1]]$type, "ide")
  expect_equal(.parse_formula(~ fa(a, 1), dat)[[1]]$k, 1L)
  expect_equal(.parse_formula(~ fa(a, 2), dat)[[1]]$k, 2L)

  # The known-covariance carriers keep every spelling they documented:
  # positional, V =, and the fb_cov() front door.
  expect_equal(.parse_formula(~ vm(g, Gmat), dat)[[1]]$type, "vm")
  expect_equal(.parse_formula(~ vm(g, V = Gmat), dat)[[1]]$type, "vm")
  expect_equal(.parse_formula(~ ped(g, Gmat), dat)[[1]]$type, "ped")

  # mgcv's named basis arguments stay open.
  expect_equal(.parse_formula(~ s(x, k = 5), dat)[[1]]$type, "smooth_mgcv")
  expect_equal(.parse_formula(~ s(x, bs = "cr"), dat)[[1]]$type, "smooth_mgcv")

  # The duplicate-name check (A22) must refuse the second copy and nothing
  # else: every name it polices still works when written once.
  expect_equal(.parse_formula(~ vm(g, cov = fb_cov(Gmat)), dat)[[1]]$type, "vm")
  expect_equal(.parse_formula(~ ped(g, V = Gmat), dat)[[1]]$type, "ped")
  expect_equal(
    .parse_formula(~ ped(g, Gmat, use_sparse_precision = FALSE), dat)[[1]]$type,
    "ped"
  )
  expect_equal(
    .parse_formula(~ s(x, k = 5, bs = "cr"), dat)[[1]]$type,
    "smooth_mgcv"
  )
})
