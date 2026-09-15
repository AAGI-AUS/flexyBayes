# =============================================================================
# Selection weights, offsets and missing responses: oracle 8.
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md section 5.5, gated as
# ledger condition 8 in section 12 -- "weighted raw GLS equals unweighted
# whitening, plus a double-weight negative control".
#
# Why the negative control is the load-bearing half. A doubled weight does
# not throw, does not diverge and does not look wrong: the provider
# converges, the coefficient path is plausible, and because the doubling is
# common to every fold, resampling reports it as signal rather than as
# error. Nothing downstream of the absorption point can see it. The only
# place it is visible is here, at the boundary where the weight is applied,
# which is why this file cares less about the identity holding than about
# the check being able to see it broken.
#
# The gates:
#
#   W1a  metamorphic. Direct weighted GLS in the original basis, solved
#        through an explicit inverse and the normal equations, equals
#        unweighted least squares on data whitened by the same weighted V,
#        solved through a Cholesky forward solve and a QR. Two numerical
#        paths that share no code.
#   W1b  independent oracle. On a residual-only V the same whitened solve
#        equals stats::lm.wfit(), an implementation this package did not
#        author.
#   W2   negative control. Applying the weights a second time -- in the
#        loss after they are already in V, or twice into V -- must change
#        the answer materially AND must be caught by the audit. Failing to
#        apply them at all must also be caught.
#   W3   the offset round-trips: subtracted once before whitening, added
#        back once to predictions, and a second add-back refuses.
#   W4   a missing response never enters the objective, occupies a
#        covariance block only under the named estimand `augment`, and the
#        estimand is recorded.
#   W5   the learned transformation reproduces from the training rows it
#        records, and a leak turns it red.
#   W6   the layer is not a second covariance. At unit weights it
#        reproduces R/selection_covariance.R's V, its whitening and its
#        coefficients exactly; with weights it must depart from them; and
#        it refuses on exactly the bases that builder refuses on.
#
# TOLERANCE, MEASURED BEFORE IT WAS WRITTEN. This is deterministic dense
# linear algebra on a well-conditioned small case: kappa(V_w) = 21.4,
# kappa(X_star) = 2.32. The spread of every quantity below across 50
# identical calls is exactly 0 -- there is no stochastic step to average
# over. The observed relative error of the W1a identity is 4.3e-16, of W1b
# against lm.wfit 3.0e-16, and the audit's own worst weight error on a
# correct object is 2.2e-16. The contract's bar of 1e-8 therefore sits
# about seven and a half orders of magnitude above the measured error and
# is not tuned to any example. For contrast, the double-weighted
# coefficients differ from the correct ones by 1.1e-1 relative: the
# separation the negative control has to work with is thirteen orders of
# magnitude wider than the bar.
# =============================================================================

fb_selw_tol <- 1e-8


# --- refusal registration -------------------------------------------- #

# The module's codes are registered by
# .populate_refusal_registry_selection_weights() in
# R/selection_weights_refusals.R, which .onLoad() calls alongside the other
# populate functions. Until that wiring lands, the codes are absent and
# .fb_refusal_condition() raises a CONTRACT VIOLATION instead of the
# refusal -- which a bare expect_error() would happily pass on. This
# installs them for the test session and is a no-op once .onLoad() does it,
# so the file behaves identically before and after the wiring.
#
# The live registry is locked, deliberately, so the codes cannot be added
# to it in place. A copy is taken, extended, re-locked and rebound; the
# rebound environment is locked exactly as the original was, so
# test-refusal-registry-scaffold.R's "the registry is locked" assertion
# still holds in the same session.
fb_selw_codes <- function() {
  c(
    "selection_weights_not_positive",
    "selection_weights_absorbed_twice",
    "selection_weights_loss_route_requires_residual_only",
    "selection_offset_absorbed_twice",
    "selection_missing_response_estimand_unnamed",
    "selection_missing_response_in_objective",
    "selection_preprocessing_leak"
  )
}

# The seven codes below are registered at load through
# .populate_refusal_registry_selection_weights(), wired into .onLoad().
# An in-test registration helper stood here while that wiring did not
# exist; it mutated and re-locked the namespace registry, which is not a
# thing a test should do once the real path works.



# --- fixtures ---------------------------------------------------------- #

# Deliberately hand-built rather than fitted. Every quantity this file
# gates is deterministic linear algebra over a known V, so a model fit
# would add minutes and Monte Carlo noise and would supply nothing the
# oracles need. Three groups of eight, one variance component, unequal
# weights spanning a factor of eight.
#
# The weights must be non-constant, and that is not cosmetic. With w == 1
# everywhere, w and w^2 coincide, so a doubling is invisible to any
# arithmetic and the negative control would pass vacuously -- the exact
# failure mode the ledger's "shown to fail" rule exists to prevent. The
# spread here gives the check a separation of 3 in relative terms to work
# with, reported as `power` by the audit.
fb_selw_fixture <- function(seed = 11L, n_group = 3L, per = 8L) {
  set.seed(seed)
  n <- n_group * per
  g <- rep(seq_len(n_group), each = per)
  Z <- stats::model.matrix(~ 0 + factor(g))
  X <- cbind(
    `(Intercept)` = 1,
    x1 = stats::rnorm(n),
    x2 = stats::runif(n, -1, 1)
  )
  list(
    n = n,
    g = g,
    X = X,
    y = as.numeric(X %*% c(2, -1.5, 0.8)) + stats::rnorm(n, sd = 0.4),
    w = rep(c(0.5, 1, 2, 4), length.out = n),
    V_random = 0.7 * tcrossprod(Z),
    sigma2 = 1.3
  )
}

fb_selw_absorbed <- function(f = fb_selw_fixture(), weights = f$w, ...) {
  .fb_sel_absorb(
    y = f$y,
    X = f$X,
    V_random = f$V_random,
    residual_variance = f$sigma2,
    weights = weights,
    ...
  )
}

# Relative on quantities above 1, absolute below, so a coefficient that
# happens to sit near zero cannot pass on a vanishing denominator.
fb_selw_rel <- function(a, b) max(abs(a - b) / pmax(1, abs(b)))


# --- W1a: the metamorphic identity ------------------------------------ #

test_that("W1a: weighted raw GLS equals unweighted whitening", {
  f <- fb_selw_fixture()

  # Route (a): weighted GLS in the ORIGINAL basis. Explicit inverse and
  # normal equations -- the naive path, chosen because it shares no code
  # with the module's Cholesky forward solve. An oracle assembled from the
  # same routine as the thing it checks proves self-consistency only.
  V_w <- f$V_random + diag(f$sigma2 / f$w)
  V_inv <- solve(V_w)
  beta_gls <- as.numeric(
    solve(crossprod(f$X, V_inv %*% f$X), crossprod(f$X, V_inv %*% f$y))
  )

  # Route (b): whiten by the same weighted V, then UNWEIGHTED least
  # squares. The weights are in V and nowhere else.
  a <- fb_selw_absorbed(f)
  s <- .fb_sel_solve(a)

  expect_null(a$provider_weights)
  expect_lt(fb_selw_rel(as.numeric(s$coefficients), beta_gls), fb_selw_tol)
  expect_lt(
    fb_selw_rel(s$fitted_working, as.numeric(f$X %*% beta_gls)),
    fb_selw_tol
  )
})


test_that("W1a is deterministic: identical calls give identical numbers", {
  # The tolerance above is only meaningful if there is no spread to hide
  # in. Measured at zero over 50 calls before the tolerance was chosen;
  # this keeps it that way, and a stochastic step introduced later fails
  # here rather than being absorbed by a widened bar.
  f <- fb_selw_fixture()
  first <- as.numeric(.fb_sel_solve(fb_selw_absorbed(f))$coefficients)
  for (i in seq_len(5L)) {
    expect_identical(
      as.numeric(.fb_sel_solve(fb_selw_absorbed(f))$coefficients),
      first
    )
  }
})


# --- W1b: the independent oracle -------------------------------------- #

test_that("W1b: on a residual-only V the whitened solve is stats::lm.wfit", {
  # The strongest oracle available for this identity is base R's own
  # weighted least squares. flexyBayes did not author it, it takes the
  # weights in the same Var(y_i) = sigma^2 / w_i sense the package lowers
  # them in, and it reaches the answer by a completely different route.
  f <- fb_selw_fixture()
  a <- .fb_sel_absorb(
    y = f$y, X = f$X, V_random = NULL,
    residual_variance = f$sigma2, weights = f$w
  )
  beta_module <- as.numeric(.fb_sel_solve(a)$coefficients)
  beta_base <- as.numeric(stats::lm.wfit(f$X, f$y, f$w)$coefficients)

  expect_lt(fb_selw_rel(beta_module, beta_base), fb_selw_tol)
})


test_that("W1b: the two absorption routes agree where both are admitted", {
  # Weights in R, or weights in the provider loss. On a residual-only V
  # those are the same estimator and the contract's "or" is a genuine
  # choice; the refusal below is what stops it being offered where it is
  # not.
  f <- fb_selw_fixture()
  in_cov <- .fb_sel_absorb(
    y = f$y, X = f$X, V_random = NULL, residual_variance = f$sigma2,
    weights = f$w, absorption = "covariance"
  )
  in_loss <- .fb_sel_absorb(
    y = f$y, X = f$X, V_random = NULL, residual_variance = f$sigma2,
    weights = f$w, absorption = "provider_loss"
  )

  expect_null(in_cov$provider_weights)
  expect_equal(in_loss$provider_weights, f$w)
  # The loss route builds V unweighted: the residual block is flat.
  expect_equal(
    diag(in_loss$V),
    rep(f$sigma2, f$n),
    tolerance = fb_selw_tol
  )

  expect_lt(
    fb_selw_rel(
      as.numeric(.fb_sel_solve(in_loss)$coefficients),
      as.numeric(.fb_sel_solve(in_cov)$coefficients)
    ),
    fb_selw_tol
  )
})


# --- W2: the double-weight negative control --------------------------- #

test_that("W2: double weighting changes the answer materially", {
  # Before asking whether the check can see the defect, establish that
  # there is a defect to see. If the double-weighted fit were numerically
  # indistinguishable from the correct one, no check could separate them
  # and this file would have to say so rather than assert a pass.
  f <- fb_selw_fixture()
  a <- fb_selw_absorbed(f)
  beta_ok <- as.numeric(.fb_sel_solve(a)$coefficients)

  # The realistic defect: the weights are already in V, and the same
  # vector is handed to the provider's `weights =` argument as well.
  rw <- sqrt(f$w)
  beta_double <- as.numeric(qr.solve(a$X_star * rw, a$y_star * rw))

  separation <- fb_selw_rel(beta_double, beta_ok)
  expect_gt(separation, 1e-3)
  # Thirteen orders of magnitude of headroom over the bar the identity is
  # judged at, so the gate is not straddling its own tolerance.
  expect_gt(separation / fb_selw_tol, 1e6)
})


test_that("W2: weights in V and in the loss are caught", {
  f <- fb_selw_fixture()
  a <- fb_selw_absorbed(f)

  # A downstream consumer that whitens by a weighted V and then also
  # weights the loss. The object is corrupted AFTER construction because
  # that is how the defect actually arises: the second application is
  # always somebody else's.
  doubled <- a
  doubled$provider_weights <- f$w

  err <- tryCatch(.fb_sel_weight_audit(doubled), error = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_selection_weights_absorbed_twice")
  expect_match(conditionMessage(err), "not exactly once")
  # It names the count, not merely the mismatch, so the reader learns
  # which defect they have.
  expect_match(conditionMessage(err), "2 time")
  expect_equal(err$times_absorbed, 2, tolerance = 1e-9)

  # And the solver refuses to touch it, so the corrupted object cannot
  # reach a provider through the module's own door.
  expect_error(.fb_sel_solve(doubled), "not exactly once")
})


test_that("W2: weights applied twice into V are caught", {
  f <- fb_selw_fixture()
  a <- fb_selw_absorbed(f)

  # The other doubling route: R built as sigma^2 / w^2. Note that the
  # `absorption` field still reads "covariance" and `provider_weights` is
  # still NULL -- every flag on the object says the weights went in once.
  # The audit reads the matrix instead, which is the whole point of
  # deriving the weighting from V rather than from the record.
  a$V <- f$V_random + diag(f$sigma2 / f$w^2)

  err <- tryCatch(.fb_sel_weight_audit(a), error = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_selection_weights_absorbed_twice")
  expect_match(conditionMessage(err), "not exactly once")
  expect_equal(err$times_absorbed, 2, tolerance = 1e-9)
})


test_that("W2: weights applied nowhere are caught too", {
  # Under-absorption is the same failure of accounting and the same
  # comparison finds it: the effective weighting is 1 where w was
  # declared.
  f <- fb_selw_fixture()
  a <- fb_selw_absorbed(f)
  a$V <- f$V_random + diag(rep(f$sigma2, f$n))

  err <- tryCatch(.fb_sel_weight_audit(a), error = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_selection_weights_absorbed_twice")
  expect_equal(err$times_absorbed, 0, tolerance = 1e-9)
})


test_that("W2: a weight folded off the residual diagonal is caught", {
  # A weight that reaches a covariance rather than R leaves the diagonal
  # correct and the matrix wrong. Without the off-diagonal half of the
  # check this passes.
  f <- fb_selw_fixture()
  a <- fb_selw_absorbed(f)
  a$V[1L, 2L] <- a$V[1L, 2L] + 0.25
  a$V[2L, 1L] <- a$V[2L, 1L] + 0.25

  err <- tryCatch(.fb_sel_weight_audit(a), error = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_selection_weights_absorbed_twice")
  expect_match(conditionMessage(err), "off the residual diagonal")
})


test_that("W2: the audit passes a correct object under either route", {
  # A check that fires on everything is not a check either.
  f <- fb_selw_fixture()
  v <- .fb_sel_weight_audit(fb_selw_absorbed(f))

  expect_lt(v$worst_relative_weight_error, fb_selw_tol)
  expect_equal(v$effective_weights, f$w, tolerance = 1e-12)
  expect_equal(v$weights_in_loss, rep(1, f$n))
  # `power` is the separation between w and w^2 the declared weights
  # admit. At zero the negative control would pass vacuously, so the gate
  # asserts its own discriminating power rather than assuming it.
  expect_gt(v$power, 1)

  v_loss <- .fb_sel_weight_audit(.fb_sel_absorb(
    y = f$y, X = f$X, V_random = NULL, residual_variance = f$sigma2,
    weights = f$w, absorption = "provider_loss"
  ))
  expect_lt(v_loss$worst_relative_weight_error, fb_selw_tol)
  expect_equal(v_loss$weights_in_covariance, rep(1, f$n))
  expect_equal(v_loss$weights_in_loss, f$w)
})


test_that("W2: unit weights are declared powerless rather than passing", {
  # With w == 1 everywhere, w and w^2 coincide and no arithmetic can see a
  # doubling -- correctly, because there is nothing to double. The audit
  # says so through `power` instead of reporting a clean bill.
  f <- fb_selw_fixture()
  v <- .fb_sel_absorb(
    y = f$y, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = rep(1, f$n)
  )
  expect_identical(.fb_sel_weight_audit(v)$power, 0)
})


# --- the absorption record -------------------------------------------- #

test_that("the object records where the weights went", {
  f <- fb_selw_fixture()
  a <- fb_selw_absorbed(f)

  log <- a$absorption_log
  expect_s3_class(log, "data.frame")
  expect_setequal(
    log$quantity,
    c("weights", "offset", "missing_response", "preprocessing")
  )
  expect_true(log$absorbed[log$quantity == "weights"])
  # The site, not merely the fact: "applied" is what a wrong object would
  # also say.
  expect_match(
    log$site[log$quantity == "weights"],
    "R, and therefore V",
    fixed = TRUE
  )
  expect_match(
    log$site[log$quantity == "weights"],
    "provider loss is unweighted",
    fixed = TRUE
  )

  b <- .fb_sel_absorb(
    y = f$y, X = f$X, V_random = NULL, residual_variance = f$sigma2,
    weights = f$w, absorption = "provider_loss"
  )
  expect_match(
    b$absorption_log$site[b$absorption_log$quantity == "weights"],
    "the provider loss",
    fixed = TRUE
  )
})


test_that("the print method names the absorption site", {
  f <- fb_selw_fixture()
  out <- utils::capture.output(print.fb_sel_absorbed(fb_selw_absorbed(f)))
  expect_true(any(grepl("weights in", out, fixed = TRUE)))
  expect_true(any(grepl("covariance", out, fixed = TRUE)))
})


# --- W3: offsets ------------------------------------------------------- #

test_that("W3: an offset is subtracted once, before whitening", {
  f <- fb_selw_fixture()
  off <- seq(-1, 1, length.out = f$n)

  with_offset <- .fb_sel_absorb(
    y = f$y, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = f$w, offset = off
  )
  pre_subtracted <- .fb_sel_absorb(
    y = f$y - off, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = f$w
  )

  expect_equal(with_offset$y_working, f$y - off, tolerance = 1e-14)
  expect_lt(
    fb_selw_rel(
      as.numeric(.fb_sel_solve(with_offset)$coefficients),
      as.numeric(.fb_sel_solve(pre_subtracted)$coefficients)
    ),
    fb_selw_tol
  )
})


test_that("W3: an offset is added back once, to predictions", {
  f <- fb_selw_fixture()
  off <- seq(-1, 1, length.out = f$n)
  a <- .fb_sel_absorb(
    y = f$y, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = f$w, offset = off
  )
  s <- .fb_sel_solve(a)

  expect_lt(
    fb_selw_rel(as.numeric(s$fitted_response), s$fitted_working + off),
    fb_selw_tol
  )
  # Round trip: response-scale prediction minus offset is the working
  # prediction, exactly once over.
  expect_equal(
    as.numeric(s$fitted_response) - off,
    s$fitted_working,
    tolerance = 1e-14
  )
})


test_that("W3: an offset subtracted twice is caught", {
  f <- fb_selw_fixture()
  off <- seq(-1, 1, length.out = f$n)
  a <- .fb_sel_absorb(
    y = f$y, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = f$w, offset = off
  )
  a$y_working <- f$y - 2 * off

  err <- tryCatch(.fb_sel_weight_audit(a), error = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_selection_offset_absorbed_twice")
  expect_match(conditionMessage(err), "does not round-trip")
  expect_equal(err$times_absorbed, 2, tolerance = 1e-9)
})


test_that("W3: an offset added back twice refuses", {
  # The prediction-side twin. A doubled offset shifts every fitted value
  # and changes no coefficient, so no coefficient-level check sees it;
  # the mark on the returned vector is what makes the second call
  # visible.
  f <- fb_selw_fixture()
  off <- seq(-1, 1, length.out = f$n)
  a <- .fb_sel_absorb(
    y = f$y, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = f$w, offset = off
  )
  once <- .fb_sel_add_offset(a, rep(0, f$n))
  expect_true(attr(once, "fb_sel_offset_applied"))
  expect_error(.fb_sel_add_offset(a, once), "already carry it")
})


test_that("W3: no offset is a no-op, not a silent zero-length branch", {
  f <- fb_selw_fixture()
  a <- fb_selw_absorbed(f)
  expect_equal(a$offset_declared, rep(0, f$n))
  expect_equal(a$y_working, f$y, tolerance = 1e-14)
  expect_false(a$absorption_log$absorbed[a$absorption_log$quantity == "offset"])
})


# --- W4: missing responses --------------------------------------------- #

test_that("W4: a missing response never enters the objective", {
  f <- fb_selw_fixture()
  y_na <- f$y
  y_na[c(4L, 19L)] <- NA_real_

  for (estimand in c("omit", "augment")) {
    a <- .fb_sel_absorb(
      y = y_na, X = f$X, V_random = f$V_random,
      residual_variance = f$sigma2, weights = f$w, na_action = estimand
    )
    expect_identical(a$missing_response_rows, c(4L, 19L))
    expect_false(any(c(4L, 19L) %in% a$objective_rows))
    expect_identical(length(a$objective_rows), f$n - 2L)
    expect_true(all(is.finite(a$y_star)))
  }
})


test_that("W4: a missing row occupies a covariance block only under augment", {
  # This is the whole content of "only under a named estimand". The
  # estimand does not change the objective -- only observed outcomes ever
  # enter that -- it changes what the covariance carries, which is what
  # section 5.4's prediction blocks need.
  f <- fb_selw_fixture()
  y_na <- f$y
  y_na[c(4L, 19L)] <- NA_real_

  omit <- .fb_sel_absorb(
    y = y_na, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = f$w, na_action = "omit"
  )
  augment <- .fb_sel_absorb(
    y = y_na, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = f$w, na_action = "augment"
  )

  expect_identical(omit$covariance_rows, omit$objective_rows)
  expect_identical(augment$covariance_rows, seq_len(f$n))
  expect_identical(nrow(augment$V), f$n)
  expect_identical(nrow(omit$V), f$n - 2L)

  # And the objective is the same one under both, so the estimand buys a
  # covariance block and nothing else.
  expect_identical(augment$objective_rows, omit$objective_rows)
  expect_lt(
    fb_selw_rel(
      as.numeric(.fb_sel_solve(augment)$coefficients),
      as.numeric(.fb_sel_solve(omit)$coefficients)
    ),
    fb_selw_tol
  )
})


test_that("W4: the estimand is recorded", {
  f <- fb_selw_fixture()
  y_na <- f$y
  y_na[4L] <- NA_real_
  a <- .fb_sel_absorb(
    y = y_na, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = f$w, na_action = "augment"
  )

  expect_identical(a$na_action, "augment")
  expect_true(
    a$absorption_log$absorbed[a$absorption_log$quantity == "missing_response"]
  )
  expect_match(
    a$absorption_log$site[a$absorption_log$quantity == "missing_response"],
    "estimand `augment`",
    fixed = TRUE
  )
})


test_that("W4: the estimand vocabulary is the package's, not a second one", {
  # R/na_action.R owns augment / omit / fail. Reading the vocabulary from
  # there rather than restating it is what stops a fourth policy added
  # over there becoming an unnamed estimand over here.
  expect_setequal(.fb_na_action_choices(), c("augment", "omit", "fail"))
  f <- fb_selw_fixture()
  for (estimand in c("augment", "omit")) {
    expect_s3_class(
      .fb_sel_absorb(
        y = f$y, X = f$X, V_random = f$V_random,
        residual_variance = f$sigma2, weights = f$w, na_action = estimand
      ),
      "fb_sel_absorbed"
    )
  }
})


test_that("W4: an unobserved row reaching the objective is caught", {
  f <- fb_selw_fixture()
  y_na <- f$y
  y_na[4L] <- NA_real_
  a <- .fb_sel_absorb(
    y = y_na, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = f$w, na_action = "omit"
  )
  a$objective_rows <- sort(c(a$objective_rows, 4L))

  err <- tryCatch(.fb_sel_weight_audit(a), error = function(e) e)
  expect_s3_class(
    err,
    "flexybayes_refusal_selection_missing_response_in_objective"
  )
  expect_match(conditionMessage(err), "reached the selection objective")
  expect_identical(err$rows, 4L)
})


test_that("W4: an unobserved row in a covariance block under omit is caught", {
  f <- fb_selw_fixture()
  y_na <- f$y
  y_na[4L] <- NA_real_
  a <- .fb_sel_absorb(
    y = y_na, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = f$w, na_action = "omit"
  )
  a$covariance_rows <- sort(c(a$covariance_rows, 4L))

  err <- tryCatch(.fb_sel_weight_audit(a), error = function(e) e)
  expect_s3_class(
    err,
    "flexybayes_refusal_selection_missing_response_estimand_unnamed"
  )
  expect_match(conditionMessage(err), "occupy a covariance block under")
})


# --- W5: transformations learned on training rows only ----------------- #

test_that("W5: the learned transformation sees the training rows only", {
  f <- fb_selw_fixture()
  train <- seq_len(16L)
  a <- .fb_sel_absorb(
    y = f$y, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = f$w, train_rows = train
  )

  expect_identical(a$preprocessing$rows_used, train)
  expect_identical(a$preprocessing$n_rows_used, 16L)

  # Pinned against stats::weighted.mean(), an implementation this package
  # did not author. Comparing the recorded numbers only against
  # .fb_sel_learn_scaling() run again would prove that the function agrees
  # with itself: a defect INSIDE the learner would move the record and the
  # recomputation together and neither this gate nor the audit's
  # reproduction check would see it. The independent oracle is what closes
  # that.
  for (col in c("x1", "x2")) {
    expect_equal(
      unname(a$preprocessing$centre[[col]]),
      stats::weighted.mean(f$X[train, col], f$w[train]),
      tolerance = 1e-14
    )
    expect_equal(
      unname(a$preprocessing$scale[[col]]),
      sqrt(stats::weighted.mean(
        (f$X[train, col] - stats::weighted.mean(f$X[train, col], f$w[train]))^2,
        f$w[train]
      )),
      tolerance = 1e-14
    )
  }

  # And it does NOT equal the transformation computed on every row, so the
  # record is a fact about which rows were seen, not a label.
  on_train <- .fb_sel_learn_scaling(f$X, train, f$w)
  on_all <- .fb_sel_learn_scaling(f$X, seq_len(f$n), f$w)
  expect_equal(a$preprocessing$centre, on_train$centre, tolerance = 1e-14)
  expect_gt(max(abs(on_train$centre - on_all$centre)), 1e-3)

  # Weighted, because a weighted problem's centre is not its unweighted
  # one. If this coincided, the transformation would be ignoring w.
  unweighted <- .fb_sel_learn_scaling(f$X, train, rep(1, f$n))
  expect_gt(max(abs(on_train$centre - unweighted$centre)), 1e-3)
})


test_that("W5: a transformation learned on all rows is caught", {
  # The injected leak of section 4.1: compute the scaling over every row
  # and file it under the training rows. Nothing about the object's shape
  # changes; only the numbers do, and only a recomputation from
  # `rows_used` can see it.
  f <- fb_selw_fixture()
  train <- seq_len(16L)
  a <- .fb_sel_absorb(
    y = f$y, X = f$X, V_random = f$V_random,
    residual_variance = f$sigma2, weights = f$w, train_rows = train
  )
  on_all <- .fb_sel_learn_scaling(f$X, seq_len(f$n), f$w)
  a$preprocessing$centre <- on_all$centre
  a$preprocessing$scale <- on_all$scale

  err <- tryCatch(.fb_sel_weight_audit(a), error = function(e) e)
  expect_s3_class(err, "flexybayes_refusal_selection_preprocessing_leak")
  expect_match(conditionMessage(err), "does not reproduce from")
  expect_setequal(err$columns, c("x1", "x2"))
})


test_that("W5: an intercept column keeps scale 1 rather than dividing by 0", {
  f <- fb_selw_fixture()
  pp <- .fb_sel_learn_scaling(f$X, seq_len(f$n), f$w)
  expect_identical(pp$constant_columns, "(Intercept)")
  expect_identical(unname(pp$scale[["(Intercept)"]]), 1)
  expect_true(all(is.finite(pp$scale)))
})


# --- W6: agreement with the real covariance builder -------------------- #

# The one place in this file that needs a model fit, and it needs one for a
# reason no hand-built covariance can supply: it is the gate that stops the
# absorption layer becoming a second implementation of V.
#
# The layer re-assembles the residual block because the weights belong in R
# and `.fb_sel_whitener()` bakes a flat R into its per-block Cholesky
# factors. Re-assembly is a fork unless something holds the two together,
# and at unit weights they must be the same object: the same V, the same
# whitening, the same coefficients. If they ever drift, this is where it
# surfaces, rather than in a downstream number nobody can attribute.
#
# john.alpha is the fixture test-selection-covariance.R uses, so the two
# files gate the same base. The fixed part is `gen`, a real fixed effect --
# an intercept-only model routes to the aggregated shape, which carries no
# `$data` and no `$fb` for either module to read.
fb_selw_fit <- function() {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  d <- agridat::john.alpha
  d$gen <- factor(d$gen)
  d$rep <- factor(d$rep)
  d$block <- factor(d$block)
  flexybayes(
    fixed = yield ~ gen, random = ~ rep + rep:block,
    data = d, backend = "inla", verbose = FALSE,
    prior = .fb_test_john_alpha_prior(d)
  )
}

test_that("W6: at unit weights the layer reproduces .fb_sel_whitener()", {
  fit <- fb_selw_fit()
  X <- stats::model.matrix(~gen, fit$data)

  a <- .fb_sel_absorb_fit(fit, X = X, weights = NULL)

  # (i) the assembled covariance is the covariance module's own dense
  #     reference, elementwise -- not merely close in some summary.
  expect_equal(
    a$V,
    .fb_sel_V_dense(fit, "plugin_mean_covariance"),
    tolerance = 1e-10
  )
  expect_identical(a$covariance_estimand, "plugin_mean_covariance")

  # (ii) the whitening is the builder's whitening. Compare the whitened
  #      response, which is what every downstream number is computed from.
  w_real <- .fb_sel_whitener(fit, "plugin_mean_covariance")
  y_real <- as.numeric(.fb_sel_whiten(w_real, fit$data$yield))
  expect_lt(fb_selw_rel(a$y_star, y_real), fb_selw_tol)

  # (iii) and so the coefficients agree. A Cholesky factor is unique only
  #       up to the block ordering the builder chooses, so the estimate --
  #       which is invariant to that choice -- is the quantity to compare.
  X_real <- .fb_sel_whiten(w_real, X)
  expect_lt(
    fb_selw_rel(
      as.numeric(.fb_sel_solve(a)$coefficients),
      as.numeric(qr.solve(X_real, y_real))
    ),
    fb_selw_tol
  )
})


test_that("W6: with weights the layer departs from the unweighted builder", {
  # The converse, and the reason the layer exists. If the weighted result
  # coincided with the builder's unweighted one, the weights would not be
  # reaching R and W6's first half would be passing vacuously.
  fit <- fb_selw_fit()
  X <- stats::model.matrix(~gen, fit$data)
  n <- nrow(fit$data)
  w <- rep(c(0.5, 1, 2, 4), length.out = n)

  a <- .fb_sel_absorb_fit(fit, X = X, weights = w)
  w_real <- .fb_sel_whitener(fit, "plugin_mean_covariance")

  expect_gt(
    max(abs(diag(a$V) - diag(.fb_sel_V_dense(fit, "plugin_mean_covariance")))),
    1e-6
  )
  expect_lt(
    fb_selw_rel(
      .fb_sel_weight_audit(a)$effective_weights,
      w
    ),
    fb_selw_tol
  )
  expect_gt(
    fb_selw_rel(
      as.numeric(.fb_sel_solve(a)$coefficients),
      as.numeric(qr.solve(
        .fb_sel_whiten(w_real, X),
        as.numeric(.fb_sel_whiten(w_real, fit$data$yield))
      ))
    ),
    1e-6
  )
})


test_that("W6: the layer refuses exactly where the builder refuses", {
  # Applicability is asked of .fb_sel_block_structure(), not re-derived. A
  # layer that admitted a base the builder rejects would be answering a
  # question the substrate cannot pose.
  d <- expand.grid(site = factor(1:4), yr = factor(1:3))
  d$y <- stats::rnorm(nrow(d))
  crossed <- list(
    data = d,
    fb = list(
      response = "y",
      random_terms = list(
        list(type = "simple", var = "site"),
        list(type = "simple", var = "yr")
      ),
      residual_terms = list(list(type = "units"))
    )
  )
  expect_false(.fb_sel_block_structure(crossed)$applicable)

  err <- tryCatch(
    .fb_sel_absorb_fit(crossed, X = stats::model.matrix(~site, d)),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_selection_covariance_not_buildable")
  expect_match(conditionMessage(err), "crossed")
})


# --- refusals ---------------------------------------------------------- #

test_that("a non-positive or non-finite weight refuses by name", {
  f <- fb_selw_fixture()
  for (bad in list(
    replace(f$w, 3L, 0),
    replace(f$w, 3L, -1),
    replace(f$w, 3L, NA_real_),
    replace(f$w, 3L, Inf)
  )) {
    err <- tryCatch(
      fb_selw_absorbed(f, weights = bad),
      error = function(e) e
    )
    expect_s3_class(err, "flexybayes_refusal_selection_weights_not_positive")
    expect_match(conditionMessage(err), "are not finite and positive")
    expect_match(conditionMessage(err), "row(s) 3", fixed = TRUE)
  }
})


test_that("the provider-loss route refuses where V has a random part", {
  # Section 5.5 offers the loss route as an alternative to R. It is only
  # the same estimator on a residual-only V: weighting the loss after
  # whitening weights observation i by w_i / (v_rand_i + sigma^2), while
  # weighted GLS weights it by 1 / (v_rand_i + sigma^2 / w_i). Refusing is
  # the alternative to quietly answering a different question.
  f <- fb_selw_fixture()
  err <- tryCatch(
    fb_selw_absorbed(f, absorption = "provider_loss"),
    error = function(e) e
  )
  expect_s3_class(
    err,
    "flexybayes_refusal_selection_weights_loss_route_requires_residual_only"
  )
  expect_match(conditionMessage(err), "residual-only V")
  expect_gt(err$largest_random_entry, 0)

  # And the two really do differ, so the refusal is protecting a real
  # difference rather than a hypothetical one.
  a <- fb_selw_absorbed(f)
  rw <- sqrt(f$w)
  V_unweighted <- f$V_random + diag(rep(f$sigma2, f$n))
  Xs <- forwardsolve(t(chol(V_unweighted)), f$X)
  ys <- as.numeric(forwardsolve(t(chol(V_unweighted)), f$y))
  beta_loss <- as.numeric(qr.solve(Xs * rw, ys * rw))
  expect_gt(
    fb_selw_rel(beta_loss, as.numeric(.fb_sel_solve(a)$coefficients)),
    1e-3
  )
})


test_that("an unnamed missing-response estimand refuses, naming the named", {
  f <- fb_selw_fixture()
  err <- tryCatch(
    fb_selw_absorbed(f, na_action = "keep"),
    error = function(e) e
  )
  expect_s3_class(
    err,
    "flexybayes_refusal_selection_missing_response_estimand_unnamed"
  )
  expect_match(conditionMessage(err), "is not a named missing-response")
  expect_match(conditionMessage(err), "augment, omit, fail", fixed = TRUE)
})


test_that("na_action = 'fail' reuses the package's own refusal", {
  # A second code for a fact the package already names would be a second
  # vocabulary, which is exactly what the module is not supposed to grow.
  f <- fb_selw_fixture()
  y_na <- f$y
  y_na[4L] <- NA_real_
  err <- tryCatch(
    .fb_sel_absorb(
      y = y_na, X = f$X, V_random = f$V_random,
      residual_variance = f$sigma2, weights = f$w, na_action = "fail"
    ),
    error = function(e) e
  )
  expect_s3_class(err, "flexybayes_refusal_missing_response_refused")
  expect_match(conditionMessage(err), "have a missing response", fixed = TRUE)
})


test_that("argument-shape errors are named, not silently coerced", {
  f <- fb_selw_fixture()
  expect_error(
    .fb_sel_absorb(
      y = f$y, X = f$X[-1L, , drop = FALSE], V_random = f$V_random,
      residual_variance = f$sigma2, weights = f$w
    ),
    "rows and `y` has",
    fixed = TRUE
  )
  expect_error(
    .fb_sel_absorb(
      y = f$y, X = f$X, V_random = f$V_random, residual_variance = -1,
      weights = f$w
    ),
    "single finite positive number",
    fixed = TRUE
  )
  expect_error(
    .fb_sel_absorb(
      y = f$y, X = f$X, V_random = f$V_random,
      residual_variance = f$sigma2, weights = f$w, offset = c(1, 2)
    ),
    "finite numeric vector of length",
    fixed = TRUE
  )
  expect_error(
    .fb_sel_weight_audit(list(a = 1)),
    "must be an <fb_sel_absorbed>",
    fixed = TRUE
  )
})


# --- registration ------------------------------------------------------ #

test_that("every refusal this module can throw is registered", {
  # The package raises a contract violation when an unregistered code is
  # thrown, which means an unregistered code surfaces as the WRONG error
  # and a bare expect_error() passes on it. This is the direct check.
  registered <- fb_refusals()$reason_code
  expect_true(all(fb_selw_codes() %in% registered))
  # The two codes this file throws that belong to other modules, reused
  # rather than duplicated: the missing-response refusal from
  # R/na_action.R and the non-PD covariance refusal from
  # R/selection_refusals.R.
  expect_true("missing_response_refused" %in% registered)
  expect_true("selection_covariance_not_positive_definite" %in% registered)
  expect_true("selection_covariance_not_buildable" %in% registered)
})


test_that("the registry is still locked after the module's codes land", {
  # The shim above rebinds the registry, and it must leave it exactly as
  # locked as it found it -- otherwise this file quietly weakens an
  # invariant another file asserts.
  expect_true(
    environmentIsLocked(get(".refusal_registry", envir = asNamespace(
      "flexyBayes"
    )))
  )
})
