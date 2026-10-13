# =============================================================================
# The convergence contract: one tri-state verdict, and it can say "unknown".
#
# The verdict is TRUE, FALSE or NA, and NA means no evidence. Read as a
# logical, NA is easy to turn into a constant: `isTRUE(x %||% TRUE)` is
# always TRUE for a missing field, and isTRUE(NA) reports absent evidence as
# a failed fit. The tests show each state can be produced and reached from
# the others.
#
# Mocks, not fits: they give combinations a real fit will not produce on
# demand, such as a finite mode status with a non-finite marginal
# likelihood.
# =============================================================================


mock_inla_fit <- function(mode_status = 0, mlik = -143.3, confirm = TRUE,
                          reasons = character(0)) {
  fit <- list(
    inla = list(
      mode = if (is.null(mode_status)) NULL else list(mode.status = mode_status),
      mlik = if (is.null(mlik)) NULL else matrix(mlik, 1L, 1L)
    ),
    num_check = if (is.null(confirm)) NULL else list(pass = confirm, reasons = reasons)
  )
  class(fit) <- c("flexybayes_inla", "flexybayes", "list")
  fit
}


# --- the three states are all reachable ------------------------------- #

test_that("a clean fit returns TRUE", {
  v <- .fb_convergence_verdict(mock_inla_fit())
  expect_true(v$ok)
  expect_identical(v$engine, "inla")
  expect_length(v$reasons, 0L)
})


test_that("a diagnostic-free fit returns NA, not TRUE", {
  # The audit's definition of done, stated directly. This is the case the
  # old code turned into TRUE.
  v <- .fb_convergence_verdict(
    mock_inla_fit(mode_status = NULL, mlik = NULL, confirm = NULL)
  )
  expect_true(is.na(v$ok))
  expect_false(isTRUE(v$ok))
})


test_that("an object of no known engine returns NA", {
  v <- .fb_convergence_verdict(structure(list(), class = "list"))
  expect_true(is.na(v$ok))
  expect_match(v$reasons, "no convergence diagnostics")
})


# --- each failure mode is detected, with a named reason --------------- #

test_that("a nonzero mode status is FALSE and says so", {
  v <- .fb_convergence_verdict(mock_inla_fit(mode_status = 1))
  expect_false(v$ok)
  expect_match(paste(v$reasons, collapse = " "), "mode status is 1")
})


test_that("a non-finite marginal likelihood is FALSE", {
  v <- .fb_convergence_verdict(mock_inla_fit(mlik = -Inf))
  expect_false(v$ok)
  expect_match(paste(v$reasons, collapse = " "), "not finite")
})


test_that("a failed numerical confirmation is FALSE and carries its reasons", {
  v <- .fb_convergence_verdict(
    mock_inla_fit(confirm = FALSE, reasons = "kld above tolerance")
  )
  expect_false(v$ok)
  expect_match(paste(v$reasons, collapse = " "), "kld above tolerance")
})


test_that("failure beats missing evidence", {
  # A fit with a nonzero status AND no likelihood is FALSE, not NA. Evidence
  # of failure is still evidence, and reporting "unknown" here would hide a
  # known bad fit behind an incomplete record.
  v <- .fb_convergence_verdict(mock_inla_fit(mode_status = 2, mlik = NULL))
  expect_false(v$ok)
})


# --- absent is not failed, and not passed ----------------------------- #

test_that("a missing mode status is NA even when everything present passes", {
  # Without the engine's own convergence flag there is no basis for TRUE,
  # however good the rest looks.
  v <- .fb_convergence_verdict(mock_inla_fit(mode_status = NULL))
  expect_true(is.na(v$ok))
  expect_match(paste(v$reasons, collapse = " "), "no evidence for")
})


test_that("an uncached confirm is recomputed, not called absent", {
  # The aggregated route caches no post-fit confirm. Where the raw INLA
  # object still carries the diagnostics, the verdict recomputes them.
  v <- .fb_convergence_verdict(mock_inla_fit(confirm = NULL))
  expect_true(isTRUE(v$ok))
  expect_true(isTRUE(v$evidence$numerical_confirm))
})


test_that("a confirm that cannot be recomputed is absent, and absent is not TRUE", {
  # Nothing to recompute from: no mode status, no likelihood, no cache, so
  # .lgm_check_numerical() returns NA, not TRUE.
  bare <- .lgm_check_numerical(list())
  expect_true(is.na(bare$pass))
  expect_match(paste(bare$reasons, collapse = " "), "no numerical diagnostic")

  v <- .fb_convergence_verdict(
    mock_inla_fit(mode_status = NULL, mlik = NULL, confirm = NULL)
  )
  expect_true(is.na(v$ok))
})


test_that("an NA confirm is treated as absent, not as failure", {
  v <- .fb_convergence_verdict(mock_inla_fit(confirm = NA))
  expect_true(isTRUE(v$ok))
})


# --- the evidence travels with the verdict ---------------------------- #

test_that("the verdict carries the evidence it was formed from", {
  v <- .fb_convergence_verdict(mock_inla_fit(mode_status = 0, mlik = -12.5))
  expect_identical(v$evidence$mode_status, 0)
  expect_identical(v$evidence$mlik, -12.5)
  # A consumer that needs a stricter rule can apply it to the evidence
  # without reimplementing the verdict.
  expect_true(all(
    c("mode_status", "mlik", "numerical_confirm",
      "max_rhat", "min_ess_bulk", "n_divergent") %in% names(v$evidence)
  ))
})


# --- the falsifier for the old defects -------------------------------- #

test_that("clean, failed and undiagnosed fits give three distinct verdicts", {
  clean <- .fb_convergence_verdict(mock_inla_fit())
  broken <- .fb_convergence_verdict(mock_inla_fit(mode_status = 1))
  absent <- .fb_convergence_verdict(
    mock_inla_fit(mode_status = NULL, mlik = NULL, confirm = NULL)
  )

  states <- c(isTRUE(clean$ok), isTRUE(broken$ok), isTRUE(absent$ok))
  expect_identical(states, c(TRUE, FALSE, FALSE))
  na_states <- c(is.na(clean$ok), is.na(broken$ok), is.na(absent$ok))
  expect_identical(na_states, c(FALSE, FALSE, TRUE))
})
