# =============================================================================
# The projection structure gate, executed cell by cell.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md section 6.3.
#
# The rule this file exists to enforce: a verdict in the structure matrix
# is a claim about what projpred does, and a claim about an external
# authority is checked by running it. So every green cell here really goes
# through `projpred::get_refmodel()` and comes back accepted, and every red
# cell really goes through it and comes back refused, with the refusal
# caught and mapped onto the package's own code.
#
# The rows that have NOT been executed -- a us() or fa() translation, a
# heterogeneous residual -- refuse as well. Unrun is not a permission, and
# the only way to turn one of those rows green is to write the test that
# runs it.
# =============================================================================


fb_sel_gate_data <- function() {
  d <- fb_sel_bayes_data(seed = 4L, J = 8L, m = 6L)
  d$t <- factor(rep(seq_len(6L), 8L))
  d
}

fb_sel_gate_fit <- function(key, expr) {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("projpred")
  fb_sel_bayes_memo(key, suppressWarnings(expr))
}

# Did projpred accept this model as a reference model? The answer, and the
# message when it did not, both come from the call.
fb_sel_refmodel_verdict <- function(brmsfit) {
  out <- tryCatch(
    {
      projpred::get_refmodel(brmsfit)
      list(accepted = TRUE, message = NA_character_)
    },
    error = function(e) list(accepted = FALSE, message = conditionMessage(e))
  )
  out
}

fb_sel_gate_row <- function(gate, feature) {
  gate[gate$feature == feature, , drop = FALSE]
}


# --- green cells --------------------------------------------------------- #

test_that("a Gaussian model with a grouped intercept is projectable, and projpred agrees", {
  skip_if_not_long_tier()
  d <- fb_sel_gate_data()
  fit <- fb_sel_gate_fit("gate_intercept", fb_brms(
    y ~ p1 + x1 + (1 | site), data = d, chains = 1L, n_samples = 300L,
    warmup = 300L, seed = 4L, verbose = FALSE, mcmc_verbose = FALSE
  ))

  gate <- .fb_sel_structure_gate(fit)
  expect_true(all(gate$verdict == "supported"))
  expect_identical(
    fb_sel_gate_row(gate, "response family")$verdict, "supported"
  )
  expect_identical(fb_sel_gate_row(gate, "(1 | site)")$verdict, "supported")
  expect_silent(.fb_sel_require_projectable(fit))

  verdict <- fb_sel_refmodel_verdict(fit$brms)
  expect_true(verdict$accepted)
})


test_that("a correlated grouped slope is projectable, and projpred agrees", {
  skip_if_not_long_tier()
  d <- fb_sel_gate_data()
  fit <- fb_sel_gate_fit("gate_slope_cor", fb_brms(
    y ~ p1 + x1 + (x1 | site), data = d, chains = 1L, n_samples = 300L,
    warmup = 300L, seed = 4L, verbose = FALSE, mcmc_verbose = FALSE
  ))

  gate <- .fb_sel_structure_gate(fit)
  expect_identical(fb_sel_gate_row(gate, "(x1 | site)")$verdict, "supported")
  expect_true(all(gate$verdict == "supported"))
  expect_true(fb_sel_refmodel_verdict(fit$brms)$accepted)
})


test_that("an uncorrelated grouped slope is refused, and the reason names the spelling that works", {
  skip_if_not_long_tier()
  # This cell came back RED against the contract's own row, which grouped
  # the three grouped-term spellings together. projpred 2.10.0 refuses
  # `(x || g)` on a flexyBayes fit and on a native brms fit alike; the
  # mathematically identical two-term spelling is accepted. Both halves are
  # executed here, because a reason that names an alternative is a claim
  # about the provider like any other.
  d <- fb_sel_gate_data()
  fit <- fb_sel_gate_fit("gate_slope_uncor", fb_brms(
    y ~ p1 + x1 + (x1 || site), data = d, chains = 1L, n_samples = 300L,
    warmup = 300L, seed = 4L, verbose = FALSE, mcmc_verbose = FALSE
  ))

  gate <- .fb_sel_structure_gate(fit)
  expect_identical(fb_sel_gate_row(gate, "(x1 || site)")$verdict, "refused")
  expect_match(fb_sel_gate_row(gate, "(x1 || site)")$reason, "\\|\\|")
  expect_error(
    .fb_sel_require_projectable(fit),
    class = "flexybayes_refusal_projection_structure_unsupported"
  )

  refused <- fb_sel_refmodel_verdict(fit$brms)
  expect_false(refused$accepted)
  expect_match(refused$message, "Unexpected number of `\\|` characters")

  two_term <- fb_sel_gate_fit("gate_slope_two_term", brms::brm(
    y ~ p1 + x1 + (1 | site) + (0 + x1 | site), data = d,
    chains = 1L, iter = 400L, refresh = 0L, seed = 4L, silent = 2L
  ))
  expect_true(fb_sel_refmodel_verdict(two_term)$accepted)
})


# --- red cells ------------------------------------------------------------ #

test_that("a known-covariance carrier is refused, and projpred refuses it too", {
  skip_if_not_long_tier()
  d <- fb_sel_gate_data()
  sites <- levels(d$site)
  K <- diag(length(sites))
  dimnames(K) <- list(sites, sites)
  fit <- fb_sel_gate_fit("gate_vm", fb_brms(
    y ~ p1 + x1, random = ~ vm(site, K), data = d,
    known_matrices = list(K = K), chains = 1L, n_samples = 300L,
    warmup = 300L, seed = 4L, verbose = FALSE, mcmc_verbose = FALSE
  ))

  gate <- .fb_sel_structure_gate(fit)
  carrier <- gate[grepl("gr(", gate$feature, fixed = TRUE), , drop = FALSE]
  expect_identical(nrow(carrier), 1L)
  expect_identical(carrier$verdict[[1L]], "refused")
  expect_match(carrier$reason[[1L]], "function call")

  expect_error(
    .fb_sel_require_projectable(fit),
    class = "flexybayes_refusal_projection_structure_unsupported"
  )

  # The cell, executed: projpred is asked, and its own refusal is what the
  # verdict above is standing in for.
  verdict <- fb_sel_refmodel_verdict(fit$brms)
  expect_false(verdict$accepted)
  expect_match(verdict$message, "right-hand side of a group-term")
})


test_that("a non-linear parameter partition is refused, and projpred refuses it too", {
  skip_if_not_long_tier()
  # flexyBayes does not emit this shape, so the cell is executed against a
  # brms fit built directly -- which is also why the gate reads a bare
  # brmsfit.
  d <- fb_sel_gate_data()
  fit <- fb_sel_gate_fit("gate_nl", brms::brm(
    brms::bf(
      y ~ prot + cand,
      prot ~ 0 + Intercept + p1 + (1 | site),
      cand ~ 0 + x1 + x2,
      nl = TRUE
    ),
    data = d,
    prior = c(
      brms::set_prior("normal(0, 5)", class = "b", nlpar = "prot"),
      brms::set_prior("normal(0, 5)", class = "b", nlpar = "cand")
    ),
    chains = 1L, iter = 400L, refresh = 0L, seed = 4L, silent = 2L
  ))

  gate <- .fb_sel_structure_gate(fit)
  row <- fb_sel_gate_row(gate, "non-linear parameter partition")
  expect_identical(nrow(row), 1L)
  expect_identical(row$verdict[[1L]], "refused")

  expect_error(
    .fb_sel_require_projectable(fit),
    class = "flexybayes_refusal_projection_structure_unsupported"
  )

  verdict <- fb_sel_refmodel_verdict(fit)
  expect_false(verdict$accepted)
  expect_match(verdict$message, "non-linear")
})


test_that("a non-Gaussian response refuses under its own code", {
  # No fit needed: the family is a property of the parsed model, and the
  # two refusals are different answers a caller acts on differently.
  d <- fb_sel_gate_data()
  d$count <- stats::rpois(nrow(d), lambda = 3)
  ir <- .build_ir_polymorphic(
    fixed = count ~ p1 + x1 + (1 | site), random = NULL, residual = NULL,
    data = d, family = "poisson", link = NULL, weights = NULL,
    known_matrices = list(), prior = NULL, prior_fixed_sd = 100,
    prior_vc_sd = 1, syntax = "auto"
  )
  gate <- .fb_sel_structure_gate(ir)
  expect_identical(fb_sel_gate_row(gate, "response family")$verdict, "refused")
  expect_error(
    .fb_sel_require_projectable(ir),
    class = "flexybayes_refusal_projection_family_unsupported"
  )
})


# --- unrun cells ---------------------------------------------------------- #

test_that("a us() translation refuses as unrun rather than being assumed", {
  d <- fb_sel_gate_data()
  ir <- .build_ir_polymorphic(
    fixed = y ~ p1 + x1, random = ~ us(t):id(site), residual = NULL,
    data = d, family = "gaussian", link = NULL, weights = NULL,
    known_matrices = list(), prior = NULL, prior_fixed_sd = 100,
    prior_vc_sd = 1, syntax = "auto"
  )
  gate <- .fb_sel_structure_gate(ir)
  unrun <- gate[gate$verdict == "unrun", , drop = FALSE]
  expect_gte(nrow(unrun), 1L)
  expect_match(paste(unrun$reason, collapse = " "), "no gate has been executed")
  expect_error(
    .fb_sel_require_projectable(ir),
    class = "flexybayes_refusal_projection_structure_unsupported"
  )
})


test_that("a heterogeneous residual refuses as unrun", {
  d <- fb_sel_gate_data()
  ir <- .build_ir_polymorphic(
    fixed = y ~ p1 + x1, random = ~site, residual = ~ at(t):units,
    data = d, family = "gaussian", link = NULL, weights = NULL,
    known_matrices = list(), prior = NULL, prior_fixed_sd = 100,
    prior_vc_sd = 1, syntax = "auto"
  )
  gate <- .fb_sel_structure_gate(ir)
  residual_rows <- gate[grepl("^residual: ", gate$feature), , drop = FALSE]
  expect_identical(nrow(residual_rows), 1L)
  expect_identical(residual_rows$verdict[[1L]], "unrun")
  expect_error(
    .fb_sel_require_projectable(ir),
    class = "flexybayes_refusal_projection_structure_unsupported"
  )
})


# --- the injected defect --------------------------------------------------- #

test_that("the gate turns red when an unsupported structure is spliced into a passing model", {
  # The same IR twice. The first passes; the second carries one further
  # random term whose cell has never been run, and nothing else changes.
  d <- fb_sel_gate_data()
  ir <- .build_ir_polymorphic(
    fixed = y ~ p1 + x1 + (1 | site), random = NULL, residual = NULL,
    data = d, family = "gaussian", link = NULL, weights = NULL,
    known_matrices = list(), prior = NULL, prior_fixed_sd = 100,
    prior_vc_sd = 1, syntax = "auto"
  )
  expect_silent(.fb_sel_require_projectable(ir))

  spliced <- ir
  spliced$random_terms <- c(
    spliced$random_terms,
    list(list(type = "ar1_spatial", var = "t", row_var = "t", col_var = "site"))
  )
  gate <- .fb_sel_structure_gate(spliced)
  expect_true(any(gate$verdict == "unrun"))
  expect_error(
    .fb_sel_require_projectable(spliced),
    class = "flexybayes_refusal_projection_structure_unsupported"
  )
})


test_that("the gate reports one row per feature and never an empty verdict", {
  d <- fb_sel_gate_data()
  ir <- .build_ir_polymorphic(
    fixed = y ~ p1 + x1 + (1 | site), random = NULL, residual = NULL,
    data = d, family = "gaussian", link = NULL, weights = NULL,
    known_matrices = list(), prior = NULL, prior_fixed_sd = 100,
    prior_vc_sd = 1, syntax = "auto"
  )
  gate <- .fb_sel_structure_gate(ir)
  expect_true(is.data.frame(gate))
  expect_identical(names(gate), c("feature", "verdict", "reason"))
  expect_true(all(gate$verdict %in% c("supported", "refused", "unrun")))
  expect_true(all(nzchar(gate$reason)))
  expect_identical(anyDuplicated(gate$feature), 0L)
})
