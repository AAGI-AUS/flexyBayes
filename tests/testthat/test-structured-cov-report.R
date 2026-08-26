# fb_structured_cov(): identified covariance reconstruction for fa terms.
#
# 0.9.3: the only engine that ever fit an fa() term (see NEWS.md) was
# withdrawn entirely. Neither active engine represents a factor-analytic
# structured-covariance term (INLA and brms both refuse it before a fit
# object exists -- see R/structured_cov_report.R), so no surviving fit
# ever carries a fa_gxe entry in parse_info$random. fb_structured_cov()
# therefore always messages and returns an empty list, regardless of
# what the fit's parse_info$random carries -- these tests were rewritten
# from a reconstruction-correctness suite (pinned arithmetic on synthetic
# posterior draws) to a message-and-refuse suite, because there is no
# longer a code path that reconstructs anything to check.

test_that("fb_structured_cov() always messages and returns empty, even with a fa_gxe term present", {
  # A fa_gxe entry in parse_info$random can only ever appear on a
  # hand-built fixture now (no active engine emits one), but the function
  # does not special-case it -- confirms the message-and-empty path is
  # unconditional, not just the "no fa term at all" case below.
  fit <- list(
    extras = list(parse_info = list(random = list(
      list(type = "fa_gxe", outer = "env", inner = "geno", k = 2L,
           n_outer = 3L, n_inner = 10L)
    )))
  )
  class(fit) <- "flexybayes"
  expect_message(out <- fb_structured_cov(fit), "no factor-analytic")
  expect_length(out, 0L)
})

test_that("fb_structured_cov messages and returns empty without an fa term", {
  fit <- list(
    extras = list(parse_info = list(random = list()))
  )
  class(fit) <- "flexybayes"
  expect_message(out <- fb_structured_cov(fit), "no factor-analytic")
  expect_length(out, 0L)
})

test_that("fb_structured_cov rejects non-flexybayes input", {
  expect_error(fb_structured_cov(list(a = 1)), "must be a flexybayes")
})
