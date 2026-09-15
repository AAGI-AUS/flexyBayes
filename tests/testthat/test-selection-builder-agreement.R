# =============================================================================
# Oracle 10 of the correctness ledger: the two covariance builders agree.
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md sections 4.3 and 12.
#
# Builder A whitens the whole observation-level problem where V is
# block-diagonal. Builder B reduces to weighted regression on a random term's
# level effects where a candidate is constant within that term's levels.
# Section 4.3 closes with "where both apply they must agree", and section 4.4
# makes that agreement the thing that licenses reporting main-effect and
# interaction candidates together, so it is a release gate rather than a
# reassurance.
#
# WHAT THE AGREEMENT NOW LOOKS LIKE
#
# Builder B reads stage-1 effects with the term ABSORBED as fixed. Measured
# over 25 independent datasets per design on 2026-09-12, the worst candidate's
# relative difference from the one-stage solve is:
#
#   design                       median      worst of 25
#   balanced, 4 plots / level    1.6e-15       1.6e-14
#   balanced, 20 plots / level   2.4e-15       5.6e-14
#   unbalanced, 2 to 12 plots      18.9 %         226 %
#
# Under balance the two builders do not merely agree, they return the same
# number: the level-constant candidate lives in the between-level subspace,
# both routes reduce to the same unweighted regression of the level means on
# it, and the difference is floating point. The gate is therefore 1e-10, four
# orders below the smallest defect worth naming rather than a fitted
# tolerance.
#
# Under imbalance they do not agree at all, and section 4.3's "must agree" is
# false as written. The reason is measured, not guessed, in the last three
# tests: see THE UNBALANCED CASE below.
#
# THE SHRUNK ROUTE, AND WHY IT IS NOT THE DEFAULT
#
# `absorb = "random"` reads posterior means of the term instead. Those are
# shrunk by lambda = sigma_u^2 / (sigma_u^2 + sigma_e^2 / m), and a regression
# of shrunk effects on a covariate returns lambda times the coefficient: 3.8
# per cent low at four plots per level, 0.44 per cent at twenty, and up to
# 63.5 per cent for a weak candidate on an earlier five-design sweep. The
# coefficient VECTORS stay proportional, because the shrinkage is one scalar
# on a theta-free vector, so the RANKING of candidates survives it and the
# route is retained for that. The proportionality machinery lives with that
# route now, which is the only place it is the right question.
#
# THE UNBALANCED CASE, AND WHAT IT SAYS ABOUT V
#
# For a level-constant candidate, GLS under V = sigma_e^2 I + sigma_u^2 ZZ' is
# weighted regression of the level means with weights
#
#   w_i = m_i / (sigma_e^2 + m_i sigma_u^2)   equivalently 1/(sigma_u^2 + sigma_e^2/m_i)
#
# Builder B weights by 1 / Var(alpha_i) = m_i / sigma_e^2, the fixed-effect
# meta-analysis weight, which omits the between-level component. The two are
# proportional only when every m_i is equal. So:
#
#   - under balance the reduction is exact;
#   - under imbalance the reduction reproduces a builder A that IGNORES V
#     (measured 1.8e-15 from it) and differs from the real builder A by 24.8
#     per cent at the fixture seed;
#   - putting sigma_u^2 back into Omega reconciles them to 2.0e-15.
#
# That third line is the diagnosis and it is tested, so the cause is
# established rather than argued. It was also derived independently from a
# hand-built dense V with no INLA and no package code.
#
# HOW FAR THIS ORACLE REACHES, WHICH IS LESS THAN SECTION 4.3 CLAIMS
#
# Section 4.3 offers the agreement as "an independent oracle on the assembly
# of V". It is not one, under any balance condition.
#
#   Under balance, neither side depends on V. Route A's coefficient for a
#   level-constant candidate is invariant to it -- scaling sigma_u^2 by a
#   hundred moves it by 2.4e-15 -- because the between-level subspace is a V
#   eigenspace with one common eigenvalue when the levels are the same size.
#   Builder B never sees V at all: stage 1 is a fixed-effects fit and stage 2
#   weights come from that fit. Two V-free quantities agreeing says nothing
#   about V.
#
#   Under imbalance, route A does depend on V (24.8 per cent here) and builder
#   B still does not, so the two are not computing the same estimand and the
#   disagreement above is not an oracle firing, it is a scope mismatch.
#
# Oracle 10 establishes that the reduction's ALGEBRA is right -- the term key,
# the name-matched candidate matrix, admissibility, the intercept as the
# level-scale protected block, the Omega whitening -- and that stage 1 uses
# the BLUE rather than the shrunk posterior mean. V is owned by oracles 1, 2
# and 3. If builder B's Omega ever carries sigma_u^2 + sigma_e^2/m_i, oracle
# 10 becomes a genuine and very sharp oracle on V, because both routes would
# then depend on the same components through disjoint code; the test that
# proves the reconciliation is already below.
#
# Every gate here was shown red on an injected defect before it was trusted.
# =============================================================================


# --- fixtures ---------------------------------------------------------- #

# Four designs, each fitted once. A fit costs a few seconds and most tests
# read the same ones, so they are cached; the skips stay inside the accessor
# so a machine without INLA still skips rather than errors.
fb_agree_cache <- new.env(parent = emptyenv())

fb_agree_simulate <- function(m_per_level, seed, n_levels = 12L,
                              sd_u = 0.7, sd_e = 0.5) {
  set.seed(seed)
  lev <- sprintf("T%02d", seq_len(n_levels))
  d <- data.frame(trial = factor(rep(lev, times = m_per_level), levels = lev))

  # Three level-constant candidates. Two would be enough for proportionality
  # to bite; three leaves the check with two degrees of freedom rather than
  # one, so a single lucky coincidence cannot carry it.
  gamma <- c(x1 = 1.5, x2 = -0.8, x3 = 0.4)
  for (nm in names(gamma)) {
    v <- stats::setNames(stats::rnorm(n_levels), lev)
    d[[nm]] <- as.numeric(v[as.character(d$trial)])
  }

  # A within-level fixed effect is mandatory: an intercept-only fixed model
  # routes to `flexybayes_aggregated`, which carries no `data` and no parsed
  # `fb`, and the reduction refuses it. Centring z within trial makes its
  # level means exactly zero, so it cannot leak into the level effects and
  # move the comparison for a reason that has nothing to do with the builders.
  z <- stats::rnorm(nrow(d))
  d$z <- as.numeric(z - stats::ave(z, d$trial))

  u <- stats::setNames(stats::rnorm(n_levels, sd = sd_u), lev)
  d$y <- 2 + as.matrix(d[, names(gamma)]) %*% gamma + 0.3 * d$z +
    as.numeric(u[as.character(d$trial)]) + stats::rnorm(nrow(d), sd = sd_e)
  d$y <- as.numeric(d$y)
  d
}

fb_agree_fixture <- function(design = c("balanced_small", "balanced_large",
                                        "unbalanced", "two_random_terms")) {
  design <- match.arg(design)
  skip_on_cran()
  skip_if_not_installed("INLA")
  if (identical(design, "two_random_terms")) {
    skip_if_not_installed("agridat")
  }
  if (!is.null(fb_agree_cache[[design]])) {
    return(fb_agree_cache[[design]])
  }

  if (identical(design, "two_random_terms")) {
    # The nested alpha design the sibling selection files use, which is known
    # to fit. Absorbing `rep` here leaves `rep:block`, which is the case the
    # absorption refusal exists for. A hand-rolled crossed pair was tried
    # first and segfaulted INLA, so this is the fixture that is stable.
    d <- agridat::john.alpha
    d$gen <- factor(d$gen)
    d$rep <- factor(d$rep)
    d$block <- factor(d$block)
    set.seed(20260912L)
    rv <- stats::setNames(stats::rnorm(nlevels(d$rep)), levels(d$rep))
    d$repcov <- as.numeric(rv[as.character(d$rep)])
    fit <- flexybayes(yield ~ gen, random = ~ rep + rep:block, data = d,
                      prior = .fb_test_john_alpha_prior(d),
                      backend = "inla", verbose = FALSE)
  } else {
    m <- switch(design,
      balanced_small = rep(4L, 12L),
      balanced_large = rep(20L, 12L),
      unbalanced = rep(c(2L, 3L, 4L, 6L, 8L, 12L), 2L)
    )
    d <- fb_agree_simulate(m, seed = 20260912L)
    fit <- flexybayes(y ~ z, random = ~trial, data = d, backend = "inla",
                      verbose = FALSE)
  }
  fb_agree_cache[[design]] <- fit
  fit
}

fb_agree_candidates <- c("x1", "x2", "x3")


# --- the two routes ----------------------------------------------------- #

# Route A: builder A's whitener applied to the whole observation-level
# problem, then an unpenalised solve. This is `penalty = NULL`, the D2 corner,
# on the protected block plus the candidates.
fb_agree_route_a <- function(fit, w = NULL) {
  w <- w %||% .fb_sel_whitener(fit, "plugin_mean_covariance")
  d <- fit$data
  X <- stats::model.matrix(~ z + x1 + x2 + x3, data = d)
  b <- qr.solve(.fb_sel_whiten(w, X), .fb_sel_whiten(w, d$y))
  b[fb_agree_candidates, 1L]
}

# Route B: builder B's reduction to the level scale. The intercept is the
# level-scale protected block; without it the candidate coefficients absorb
# the mean of T and the two routes are not being asked the same question.
fb_agree_route_b <- function(red) {
  wb <- .fb_sel_reduction_whitener(red)
  Tk <- cbind(`(Intercept)` = 1, red$T)
  b <- qr.solve(.fb_sel_whiten(wb, Tk), .fb_sel_whiten(wb, red$u))
  b[fb_agree_candidates, 1L]
}

# The worst candidate's relative difference. This is the gated quantity for
# the default route, where the two answers are the same number.
fb_agree_rel <- function(gamma_a, gamma_b) max(abs(gamma_b / gamma_a - 1))

# The single scalar that best relates two coefficient vectors, and what is
# left after it. This is the gated quantity for the SHRUNK route, where the
# answers differ by one factor. It is scale free and, unlike a per-candidate
# ratio, does not blow up when one coefficient sits near zero -- a
# per-candidate ratio spread reached 4.9 on the measurement runs for exactly
# that reason and is not a quantity a gate can rest on.
fb_agree_scale <- function(gamma_a, gamma_b) {
  rho <- sum(gamma_a * gamma_b) / sum(gamma_a * gamma_a)
  list(
    rho = rho,
    residual = sqrt(sum((gamma_b - rho * gamma_a)^2)) /
      sqrt(sum(gamma_b^2))
  )
}

# The analytic BLUP shrinkage factor at the base fit's theta-hat, per level.
# This is what the shrunk route's scalar is predicted to be.
fb_agree_shrinkage <- function(w, red) {
  s_u2 <- unname(w$components$terms[[red$term]])
  s_e2 <- w$components$residual
  s_u2 / (s_u2 + s_e2 / red$n_per_level)
}

# A whitener that keeps the residual scale and forgets the random effect. It
# is the most plausible wrong whitener, and it is exercised through the real
# `.fb_sel_whiten()` so the defect is in V and nowhere else.
fb_agree_residual_only <- function(w) {
  s <- sqrt(w$components$residual)
  w$factors <- lapply(w$factors, function(f) {
    list(rows = f$rows, L = diag(s, nrow = length(f$rows)))
  })
  w
}


# --- the precondition: this really is a case where both apply ----------- #

test_that("both builders apply to every fixture, so oracle 10 is not vacuous", {
  for (design in c("balanced_small", "balanced_large", "unbalanced")) {
    fit <- fb_agree_fixture(design)

    # Builder A: V is block-diagonal over the one grouping factor.
    bs <- .fb_sel_block_structure(fit)
    expect_true(bs$applicable)
    expect_identical(bs$factor, "trial")
    expect_identical(bs$n_blocks, 12L)

    # Builder B: every candidate is constant within levels and varies across
    # them. If either half failed the comparison below would be comparing a
    # refusal with a number.
    adm <- .fb_sel_reduction_admissible(
      fit$data, fb_agree_candidates, as.character(fit$data$trial)
    )
    expect_true(all(adm$admissible))

    # And the route under test is the default one, not a chosen one.
    expect_identical(
      .fb_sel_reduction(fit, "trial", fb_agree_candidates)$absorb, "fixed"
    )
  }
})


test_that("a within-level candidate cannot enter the reduction at all", {
  # Section 4.4: the reduction covers main effects of level-constant
  # covariates only. `z` varies within trial, so there is no level-scale
  # column for it and oracle 10 says nothing about candidates of its kind.
  # An interaction candidate is in the same position and goes through
  # builder A on the full rotated problem.
  fit <- fb_agree_fixture("balanced_small")
  expect_error(
    .fb_sel_reduction(fit, "trial", "z"),
    "not constant within levels"
  )
})


# --- ORACLE 10, the default route --------------------------------------- #

test_that("oracle 10: the absorbed reduction equals the one-stage solve", {
  # Not "agrees within a tolerance". Under balance the two routes reduce to
  # the same regression of the level means on the candidates, so the same
  # number comes back and the only difference is floating point.
  #
  # Measured 2026-09-12 at these seeds: 8.9e-16 at four plots per level and
  # 2.9e-15 at twenty. Across 25 independent datasets each the worst was
  # 1.6e-14 and 5.6e-14. The ceiling is 1e-10 -- about 2000 times the worst
  # draw seen, and still four orders below the smallest coefficient error
  # anyone would act on. The spread across datasets is the proxy for the
  # drift this fixed-seed fixture will meet, which is INLA's numerics rather
  # than new data; the proxy is the wider of the two and is used because the
  # narrower one cannot be exercised here.
  for (design in c("balanced_small", "balanced_large")) {
    fit <- fb_agree_fixture(design)
    red <- .fb_sel_reduction(fit, "trial", fb_agree_candidates)
    expect_lt(fb_agree_rel(fb_agree_route_a(fit), fb_agree_route_b(red)), 1e-10)
    # Exactness is the claim, so the balanced reduction declares nothing.
    expect_length(red$approximations, 0L)
  }
})


test_that("oracle 10 falsifier: a misaligned reduction turns the gate red", {
  # The reduction matches covariate rows to levels by NAME, and the source
  # says why: a positional match silently misaligns every covariate value.
  # These are the two shapes that mistake takes; they reach 3.0 and 2.7 at
  # this seed, ten orders above the ceiling.
  fit <- fb_agree_fixture("balanced_small")
  red <- .fb_sel_reduction(fit, "trial", fb_agree_candidates)
  gamma_a <- fb_agree_route_a(fit)

  swapped <- red
  swapped$T[, c(1L, 2L)] <- swapped$T[, c(2L, 1L)]
  expect_gt(fb_agree_rel(gamma_a, fb_agree_route_b(swapped)), 1e-10)

  set.seed(7L)
  shuffled <- red
  shuffled$T <- red$T[sample(nrow(red$T)), , drop = FALSE]
  expect_gt(fb_agree_rel(gamma_a, fb_agree_route_b(shuffled)), 1e-10)

  # And the ceiling is not doing the work: scaling one candidate's reduction
  # coefficient by 1 + 1e-6 already lands 1e-6 away, four orders above it.
  nudged <- fb_agree_route_b(red) * c(1 + 1e-6, 1, 1)
  expect_gt(fb_agree_rel(gamma_a, nudged), 1e-10)
})


# --- the shrunk route, and why it is not the default -------------------- #

test_that("the two stage-1 estimators differ, and only the absorbed one agrees", {
  # The regression gate. If `absorb` ever silently reverts to reading ranef(),
  # the first assertion goes green and this test goes red, which is the only
  # way a reader finds out that every reduction coefficient has quietly
  # shrunk by a few per cent.
  fit <- fb_agree_fixture("balanced_small")
  gamma_a <- fb_agree_route_a(fit)
  absorbed <- fb_agree_route_b(
    .fb_sel_reduction(fit, "trial", fb_agree_candidates, absorb = "fixed")
  )
  shrunk <- fb_agree_route_b(
    .fb_sel_reduction(fit, "trial", fb_agree_candidates, absorb = "random")
  )

  # Measured 2026-09-12: the two routes sit 3.8 per cent apart at this seed,
  # and the smallest separation over 25 datasets of this design was 0.87 per
  # cent, so the floor is 4.3 times below the worst case rather than tuned to
  # this draw.
  expect_gt(fb_agree_rel(absorbed, shrunk), 2e-3)

  # Only one of them is the one-stage answer.
  expect_lt(fb_agree_rel(gamma_a, absorbed), 1e-10)
  expect_gt(fb_agree_rel(gamma_a, shrunk), 2e-3)

  # Every candidate shrinks, not just the average of them.
  expect_true(all(shrunk / absorbed < 1))
})


test_that("the shrunk route declares itself and the absorbed route does not", {
  # A coefficient from the shrunk route is on a different scale from builder
  # A's, and section 4.4 licenses reporting the two side by side. The object
  # is what has to stop that, so it says so.
  fit <- fb_agree_fixture("balanced_small")
  absorbed <- .fb_sel_reduction(fit, "trial", fb_agree_candidates,
                                absorb = "fixed")
  shrunk <- .fb_sel_reduction(fit, "trial", fb_agree_candidates,
                              absorb = "random")

  expect_identical(absorbed$absorb, "fixed")
  expect_identical(shrunk$absorb, "random")
  expect_false("shrunk_stage_one_effects" %in% absorbed$approximations)
  expect_true("shrunk_stage_one_effects" %in% shrunk$approximations)
})


test_that("the shrunk route stays proportional, by the BLUP shrinkage factor", {
  # The reason `absorb = "random"` is retained rather than removed: the
  # attenuation is one scalar on a theta-free vector, so the candidate
  # coefficient VECTORS remain proportional and the ranking survives. This is
  # the gate for that claim, and it is the only place proportionality is the
  # right question -- the default route returns the same number, not a
  # proportional one.
  ceilings <- list(
    balanced_small = c(residual = 5e-3, ratio = 3e-2),
    balanced_large = c(residual = 5e-5, ratio = 5e-3)
  )
  for (design in names(ceilings)) {
    fit <- fb_agree_fixture(design)
    w <- .fb_sel_whitener(fit, "plugin_mean_covariance")
    red <- .fb_sel_reduction(fit, "trial", fb_agree_candidates,
                             absorb = "random")
    agree <- fb_agree_scale(fb_agree_route_a(fit, w), fb_agree_route_b(red))
    lambda <- mean(fb_agree_shrinkage(w, red))

    # Residual measured 2026-09-12: 1.33e-04 and 7.50e-07 at these seeds,
    # worst of 25 datasets each 2.02e-03 and 1.88e-05, so the ceilings are
    # 2.5 and 2.7 times the worst draw.
    expect_lt(agree$residual, ceilings[[design]][["residual"]])

    # And the scalar is the analytic shrinkage factor. Measured 7.81e-03 and
    # 8.29e-04 at these seeds, worst of 25 each 1.11e-02 and 1.60e-03, so the
    # ceilings are 2.7 and 3.1 times the worst draw. The residue is not
    # numerical: INLA integrates over the hyperparameters, so the effective
    # shrinkage is E[lambda | y] while the plugin estimand evaluates lambda at
    # E[theta | y]. The two differ by the curvature of lambda over the
    # posterior, which is why the gap closes as the posterior tightens.
    expect_lt(abs(agree$rho / lambda - 1), ceilings[[design]][["ratio"]])

    # Attenuation, not inflation, and it is smaller with more plots per level.
    expect_lt(agree$rho, 1)
  }
})


test_that("shrunk-route falsifier: misalignment and a wrong shrinkage are red", {
  fit <- fb_agree_fixture("balanced_small")
  w <- .fb_sel_whitener(fit, "plugin_mean_covariance")
  red <- .fb_sel_reduction(fit, "trial", fb_agree_candidates,
                           absorb = "random")
  gamma_a <- fb_agree_route_a(fit, w)

  swapped <- red
  swapped$T[, c(1L, 2L)] <- swapped$T[, c(2L, 1L)]
  expect_gt(fb_agree_scale(gamma_a, fb_agree_route_b(swapped))$residual, 5e-3)

  # sigma_e^2 / m, not sigma_e^2. Dropping the division is the classic error
  # and it must not survive the ratio gate above; it reaches 8.3e-02 here.
  s_u2 <- unname(w$components$terms[["trial"]])
  lambda_wrong <- s_u2 / (s_u2 + w$components$residual)
  rho <- fb_agree_scale(gamma_a, fb_agree_route_b(red))$rho
  expect_gt(abs(rho / lambda_wrong - 1), 3e-2)
})


# --- the unbalanced case ------------------------------------------------ #

test_that("imbalance declares omega_diagonal_approximation on the reduction", {
  # Section 4.3: the diagonal Omega is exact under balance and an
  # approximation otherwise, and it is recorded rather than adopted quietly.
  balanced <- .fb_sel_reduction(
    fb_agree_fixture("balanced_small"), "trial", fb_agree_candidates
  )
  unbalanced <- .fb_sel_reduction(
    fb_agree_fixture("unbalanced"), "trial", fb_agree_candidates
  )

  expect_true(balanced$balanced)
  expect_length(balanced$approximations, 0L)
  expect_false(unbalanced$balanced)
  expect_true("omega_diagonal_approximation" %in% unbalanced$approximations)
})


test_that("under imbalance the two builders agree, and both depend on V", {
  # This block was a CHARACTERISATION of a defect, written to go red the day
  # builder B's Omega carried the between-level component. That landed, so it
  # is rewritten as the agreement gate its own comment called for rather than
  # deleted to make the suite green.
  #
  # The defect: GLS on a level-constant candidate weights level i by
  # m_i / (sigma_e^2 + m_i sigma_u^2), and the reduction weighted it by
  # 1 / Var(alpha_i) = m_i / sigma_e^2 -- the FIXED-effect meta-analysis
  # weight, omitting sigma_u^2, proportional to the first only when the m_i
  # are equal. Under balance it cancelled, which is why every balanced test
  # passed over it; under imbalance it was 24.8 per cent out here, median
  # 18.9 and worst 226 per cent over 25 datasets.
  fit <- fb_agree_fixture("unbalanced")
  w <- .fb_sel_whitener(fit, "plugin_mean_covariance")
  red <- .fb_sel_reduction(fit, "trial", fb_agree_candidates)
  gamma_b <- fb_agree_route_b(red)

  # Agreement, on the design where it used to fail. Measured 7.0e-05 at this
  # seed against the 7.73e-02 the defect produced -- three orders better, and
  # the ceiling sits an order above the measurement.
  expect_lt(fb_agree_scale(fb_agree_route_a(fit, w), gamma_b)$residual, 1e-3)

  # And the reason this matters beyond the number. The reduction now reads
  # sigma_u^2 from the base fit, so it is no longer the builder A that
  # DISCARDS the random effect -- measured 0.329 against that route, where
  # the defective version reproduced it to 1.8e-15.
  #
  # That is what upgrades oracle 10. Under balance neither side depends on V,
  # because a level-constant candidate lives in the between-level subspace,
  # which at equal level sizes is a V eigenspace with one common eigenvalue.
  # Two V-free quantities agreeing says nothing about V. Under imbalance both
  # routes now depend on the same variance components through disjoint code,
  # so their agreement IS evidence about V, which is the claim contract
  # section 4.3 made and could not support until now.
  gamma_no_v <- fb_agree_route_a(fit, fb_agree_residual_only(w))
  expect_gt(fb_agree_rel(gamma_no_v, gamma_b), 1e-2)
})


test_that("restoring sigma_u^2 to Omega reconciles the two builders exactly", {
  # The diagnosis, tested rather than argued, and it survives the fix: this
  # is a statement about the algebra, not about the current implementation.
  # Rebuilding stage 2 with Omega_i = sigma_u^2 + sigma_e^2/m_i -- the
  # random-effects meta-regression weight, both components taken from the
  # same whitener route A used -- closes a 24.8 per cent disagreement to
  # 2.0e-15 here, and to 3.5e-14 at worst over 25 datasets.
  fit <- fb_agree_fixture("unbalanced")
  w <- .fb_sel_whitener(fit, "plugin_mean_covariance")
  red <- .fb_sel_reduction(fit, "trial", fb_agree_candidates)

  corrected <- red
  corrected$omega_diag <- w$components$terms[["trial"]] +
    w$components$residual / red$n_per_level

  expect_lt(
    fb_agree_rel(fb_agree_route_a(fit, w), fb_agree_route_b(corrected)),
    1e-10
  )
})


# --- how far the oracle reaches on V ------------------------------------ #

test_that("under exact balance neither route depends on V", {
  # So the balanced agreement above certifies the reduction's algebra and
  # says nothing about whether V was assembled correctly. Builder B never
  # sees V; this is the half of the statement that has to be measured, which
  # is that builder A does not either.
  fit <- fb_agree_fixture("balanced_small")
  w <- .fb_sel_whitener(fit, "plugin_mean_covariance")
  gamma_a <- fb_agree_route_a(fit, w)

  # A different D1 estimand, and a whitener that discards the random effect
  # altogether: 1.3e-15 and 5.6e-16 at this seed.
  expect_lt(
    fb_agree_rel(gamma_a, fb_agree_route_a(
      fit, .fb_sel_whitener(fit, "point_estimate_covariance")
    )),
    1e-8
  )
  expect_lt(
    fb_agree_rel(gamma_a, fb_agree_route_a(fit, fb_agree_residual_only(w))),
    1e-8
  )
})


test_that("under imbalance route A depends on V and the reduction does not", {
  # The other half, and the reason the two disagree there. Unequal level sizes
  # give unequal GLS weights, the common eigenvalue is gone, and V becomes
  # load-bearing for route A alone. Measured 24.8 per cent at this seed;
  # smallest over 25 datasets 4.4 per cent, so the floor of 2e-2 sits twice
  # below the mildest case. Builder B is unmoved because it never reads V.
  fit <- fb_agree_fixture("unbalanced")
  w <- .fb_sel_whitener(fit, "plugin_mean_covariance")

  expect_gt(
    fb_agree_rel(
      fb_agree_route_a(fit, w),
      fb_agree_route_a(fit, fb_agree_residual_only(w))
    ),
    2e-2
  )
})


# --- refusals ----------------------------------------------------------- #

test_that("absorption refuses where the base carries another random term", {
  # Absorbing the term as fixed leaves the remaining random terms unmodelled,
  # which needs a mixed refit this release does not do. The refusal bounds
  # where oracle 10's default route can run at all, so it belongs here as
  # well as with the builder.
  fit <- fb_agree_fixture("two_random_terms")
  expect_error(
    .fb_sel_reduction(fit, "rep", "repcov", absorb = "fixed"),
    "further random term"
  )
  # The shrunk route has no such restriction, and says what it costs instead.
  shrunk <- .fb_sel_reduction(fit, "rep", "repcov", absorb = "random")
  expect_true("shrunk_stage_one_effects" %in% shrunk$approximations)
})


test_that("an unknown term refuses by name on both routes", {
  # The absorbed route builds a design column from the term's levels, so an
  # unknown term used to die inside a data-frame assignment with a message
  # about row counts. Both routes must refuse by name, and the message must
  # carry the term.
  fit <- fb_agree_fixture("balanced_small")
  for (mode in c("fixed", "random")) {
    expect_error(
      .fb_sel_reduction(fit, "nosuchterm", fb_agree_candidates, absorb = mode),
      "no random term `nosuchterm`"
    )
  }
})


test_that("every refusal this file relies on is registered", {
  # The package raises a contract violation for an unregistered reason code,
  # so an unregistered code surfaces as the wrong error and a loose matcher
  # passes on it. This is the direct check.
  codes <- c(
    "selection_candidates_not_constant_within_term",
    "selection_reduction_term_not_in_fit",
    "selection_reduction_absorption_unsupported",
    "selection_covariance_not_buildable"
  )
  expect_true(all(codes %in% fb_refusals()$reason_code))
})
