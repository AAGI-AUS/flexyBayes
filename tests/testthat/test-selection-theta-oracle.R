# =============================================================================
# Selection theta: oracle 3 of the correctness ledger.
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md section 12, oracle 3.
#
# The whitening operator is only as right as the theta-hat it is built from,
# and oracles 1 and 2 cannot see a wrong theta -- they check that V is
# assembled correctly FROM theta and that W inverts it. A theta that is simply
# wrong passes both. This file is the independent-engine check that closes
# that gap.
#
# Audit A29 says the load-bearing nested component is not externally
# validated: the ASReml scenario fits `rep + rep:block` and requests only
# `rep`, whose coefficient of variation is 1.118, above the comparator's own
# 0.50 cutoff, so it is not asserted either. `rep:block` is requested nowhere.
# The audit's remedy needs a converged ASReml refit and the licence server is
# down.
#
# It does not need one. Three things close it without a licence:
#
#   3a  every REML point estimate lies inside the flexyBayes credible
#       interval. This is the right cross-summary comparison, because a
#       posterior mean and a REML estimate are different summaries of the
#       same theta and are not required to coincide.
#   3b  where the posterior is tight, the two agree to within the posterior's
#       own width. No tuned tolerance: the criterion scales itself.
#   3c  lme4 reproduces the BANKED ASReml values. That is what licences lme4
#       as the stand-in while the licence server is unavailable, rather than
#       assuming it.
#
# Scope, stated so it is not overread. This validates the nested covariance
# BUILDER on a replicated nested design with both engines converged. It is
# not a general ASReml-equivalence claim, and the barley arm's spatial term
# is a separate structure with its own gate.
# =============================================================================


fb_theta_fixture <- function() {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  skip_if_not_installed("lme4")
  d <- agridat::john.alpha
  d$gen <- factor(d$gen)
  d$rep <- factor(d$rep)
  d$block <- factor(d$block)
  # lme4 needs the nested stratum as one factor; flexyBayes reads rep:block
  # from the formula. Same partition, spelled for each engine.
  d$rb <- factor(paste0(as.character(d$rep), ":", as.character(d$block)))

  fit <- flexybayes(
    fixed = yield ~ gen, random = ~ rep + rep:block,
    data = d, backend = "inla", verbose = FALSE,
    prior = .fb_test_john_alpha_prior(d)
  )
  reml <- lme4::lmer(yield ~ gen + (1 | rep) + (1 | rb), data = d, REML = TRUE)

  vc <- summary(fit)$varcomp
  lv <- as.data.frame(lme4::VarCorr(reml))
  # Map each flexyBayes component to its lme4 group. A positional match here
  # would silently compare the wrong pair.
  grp <- c(sigma = "Residual", "sd_rep" = "rep", "sd_rep:block" = "rb")

  cmp <- data.frame(
    component = vc$component,
    fb_sd = vc$estimate,
    fb_sd_se = vc$std.error,
    fb_lo = vc$conf.low,
    fb_hi = vc$conf.high,
    stringsAsFactors = FALSE
  )
  cmp$reml_sd <- vapply(
    cmp$component, function(k) lv$sdcor[lv$grp == grp[[k]]], numeric(1L)
  )
  cmp$post_cv <- cmp$fb_sd_se / cmp$fb_sd
  cmp$rel_gap <- abs(cmp$fb_sd - cmp$reml_sd) / cmp$reml_sd

  list(fit = fit, reml = reml, cmp = cmp)
}


# --- ORACLE 3a: the interval contains the independent estimate --------- #

test_that("oracle 3a: every REML estimate lies inside the posterior interval", {
  f <- fb_theta_fixture()
  cmp <- f$cmp

  inside <- cmp$reml_sd >= cmp$fb_lo & cmp$reml_sd <= cmp$fb_hi
  expect_true(
    all(inside),
    info = paste0(
      "components outside their credible interval: ",
      paste(cmp$component[!inside], collapse = ", ")
    )
  )

  # Including the one the audit says is unvalidated. Named explicitly so the
  # gap cannot silently reopen by the component being dropped from the frame.
  expect_true("sd_rep:block" %in% cmp$component)
})


test_that("oracle 3a falsifier: a displaced reference falls outside", {
  f <- fb_theta_fixture()
  cmp <- f$cmp
  # Double every reference. The residual's interval is tight, so at least one
  # component must now fall outside. A check that survives this is not
  # checking containment.
  outside_any <- any(2 * cmp$reml_sd < cmp$fb_lo | 2 * cmp$reml_sd > cmp$fb_hi)
  expect_true(outside_any)
})


# --- ORACLE 3b: agreement, scaled by the posterior's own width --------- #

test_that("oracle 3b: tight components agree to within the posterior width", {
  f <- fb_theta_fixture()
  cmp <- f$cmp

  # The criterion is |fb - reml| / reml < posterior CV. It needs no tuned
  # tolerance: a component the data pins down must agree closely, and one it
  # does not is allowed to differ by as much as the posterior is wide.
  agrees <- cmp$rel_gap < cmp$post_cv

  expect_true(agrees[cmp$component == "sigma"])
  expect_true(agrees[cmp$component == "sd_rep:block"])

  # A snapshot of which components meet the criterion. This is the
  # identification state of the design, and it should not drift unnoticed in
  # either direction: a component that starts agreeing, or stops, is a change
  # worth reading.
  expect_setequal(cmp$component[agrees], c("sigma", "sd_rep:block"))
})


test_that("oracle 3b: the nested component agrees closely in absolute terms", {
  f <- fb_theta_fixture()
  cmp <- f$cmp
  gap <- cmp$rel_gap[cmp$component == "sd_rep:block"]
  # Measured at 0.072 on 2026-09-12. The bound is deliberately loose enough
  # to be a regression gate rather than a restatement of today's number.
  expect_lt(gap, 0.15)
})


# --- ORACLE 3c: lme4 is licensed as the ASReml stand-in ---------------- #

test_that("oracle 3c: lme4 reproduces the banked ASReml components", {
  f <- fb_theta_fixture()

  path <- system.file(
    "validation/asreml-4.2/02_nested_stratum.rds", package = "flexyBayes"
  )
  skip_if(identical(path, ""), "banked ASReml fixture is not installed")
  rec <- readRDS(path)

  # A26's contract: a nonconverged fit is not admissible as a reference. This
  # scenario converged; 06 and 07 did not, which is why they are not used here.
  expect_true(isTRUE(rec$asreml_converged))

  asr_var <- rec$raw$asreml_varcomp
  asr_sd <- sqrt(asr_var)
  key <- c("units!R" = "Residual", "rep" = "rep", "rep:block" = "rb")

  lv <- as.data.frame(lme4::VarCorr(f$reml))
  reml_sd <- vapply(
    names(key), function(k) lv$sdcor[lv$grp == key[[k]]], numeric(1L)
  )

  # Two REML implementations on the same likelihood should agree to far more
  # than the tolerance any downstream claim needs. Measured agreement is to
  # roughly five significant figures.
  expect_equal(unname(reml_sd), unname(asr_sd[names(key)]), tolerance = 1e-4)
})


test_that("oracle 3c falsifier: a wrong component mapping fails", {
  f <- fb_theta_fixture()
  path <- system.file(
    "validation/asreml-4.2/02_nested_stratum.rds", package = "flexyBayes"
  )
  skip_if(identical(path, ""), "banked ASReml fixture is not installed")
  asr_sd <- sqrt(readRDS(path)$raw$asreml_varcomp)

  lv <- as.data.frame(lme4::VarCorr(f$reml))
  # Swap rep and rep:block. If the comparison cannot see a transposed mapping
  # it is comparing multisets of numbers, not components.
  swapped <- c(
    lv$sdcor[lv$grp == "rb"], lv$sdcor[lv$grp == "rep"],
    lv$sdcor[lv$grp == "Residual"]
  )
  expect_false(isTRUE(all.equal(
    swapped, unname(asr_sd[c("rep", "rep:block", "units!R")]),
    tolerance = 1e-4
  )))
})


# --- what this does and does not establish ----------------------------- #

test_that("the weakly identified component is recorded, not asserted away", {
  f <- fb_theta_fixture()
  cmp <- f$cmp
  rep_row <- cmp[cmp$component == "sd_rep", ]

  # `rep` rests on three levels. Its posterior CV exceeds one, the point
  # summaries differ by more than a posterior width, and REML still lies
  # inside the interval. All three are true at once and none of them is a
  # defect: that combination is what weak identification looks like. It is
  # recorded here so that a future change in any direction is visible.
  expect_gt(rep_row$post_cv, 1)
  expect_gt(rep_row$rel_gap, rep_row$post_cv)
  expect_true(rep_row$reml_sd >= rep_row$fb_lo && rep_row$reml_sd <= rep_row$fb_hi)
})
