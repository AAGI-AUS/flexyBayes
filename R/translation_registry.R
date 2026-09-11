# translation_registry.R -- what each ASReml structure becomes, and what
# grounds the claim.
#
# The capability matrix (R/capability_matrix.R) answers one question: does an
# engine fit this structure, or refuse it? The 2026-09-11 translation audit
# found that a second question was being answered implicitly and inconsistently
# across the public surfaces -- is the model that gets fitted the same model
# ASReml would fit, and how do we know? "ASReml syntax" was used as a blanket
# claim over lowerings of three quite different kinds.
#
# This table separates them. Every accepted structure carries:
#
#   vintage    which grammar the spelling belongs to. ASReml 4 removed
#              residual-side at() in favour of dsum(), and ped() in favour of
#              vm(); both spellings still work here, and both are now labelled
#              as the ASReml 3 forms they are rather than presented as current
#              syntax.
#
#   class      what the lowering promises:
#
#              exact             same response likelihood, same fixed and
#                                random design, same covariance family, same
#                                residual covariance, same missing-row
#                                estimand. Estimation and prior still differ:
#                                ASReml returns AI-REML estimates and
#                                flexyBayes returns a posterior, so equality
#                                of point estimates is not claimed by any row
#                                in this table.
#              covariance-family the fitted covariance set is the same after a
#                                documented parameter map, while
#                                identification, priors or the inferential
#                                target differ.
#              conditional       covariance-family equivalent only where a
#                                stated precondition holds. The precondition
#                                is in the `condition` column and is not
#                                currently enforced.
#              analogue          scientifically related, with a different
#                                basis, penalty, nugget or prediction rule.
#                                Must not be described as the same model.
#              refused           cannot be represented faithfully, and
#                                terminates with a typed condition.
#
#   oracle     what independently grounds the row, with the version it was
#              taken at. `none at this anchor` is a truthful entry and several
#              rows carry it: the 0.10 development line changed or added the
#              INLA factor-analytic, us()/iidKd, correlated-slope, sectioned
#              residual and combined MET paths, and the banked ASReml
#              comparisons predate all of them. Cross-engine agreement is not
#              an oracle for these rows -- both engines consume the same
#              parsed representation, so a common parser error triangulates
#              perfectly.
#
# Rendered by tools/generate_translation_registry.R to a Markdown block and to
# inst/translation/translation_registry.csv, and pinned by
# tests/testthat/test-translation-registry.R, which fails when a committed
# surface drifts from a fresh render.

.FB_TRANSLATION_CLASSES <- c(
  "exact",
  "covariance-family",
  "conditional",
  "analogue",
  "refused"
)

.FB_TRANSLATION_VINTAGES <- c(
  "ASReml 4",
  "ASReml 3 alias",
  "lme4 / brms",
  "flexyBayes extension"
)

#' The ASReml-to-flexyBayes translation registry
#'
#' @returns A data.frame with one row per accepted or refused structure and
#'   the columns `structure`, `spelling`, `vintage`, `engine`, `class`,
#'   `condition`, `oracle`, and `anchor`.
#'
#' @noRd
#' @keywords internal
.fb_translation_registry <- function() {
  rows <- list(
    list(
      structure = "Fixed factors, covariates, means model",
      spelling = "`fixed = ~ trt + x`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "exact",
      condition = "Ordinary R contrast rules; rank-deficiency policy is R's.",
      oracle = "lm() / lme4 unit tests at this anchor",
      anchor = "test-backend-conformance-open-oracle.R"
    ),
    list(
      structure = "Simple iid random intercept",
      spelling = "`~ g`, `id(g)`, `ide(g)`, `(1 | g)`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "exact",
      condition = "",
      oracle = "ASReml 4.2.0.392 at flexyBayes 0.9.0; lme4 REML at this anchor",
      anchor = "test-backend-conformance-open-oracle.R"
    ),
    list(
      structure = "Nested and interacted iid terms",
      spelling = "`~ rep:block`, `~ gen:env`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "exact",
      condition = "Each component term remains iid.",
      oracle = "ASReml 4.2.0.392 at flexyBayes 0.9.0 (RCBD / IBD / MET)",
      anchor = "test-grammar-nested-surface-parity.R"
    ),
    list(
      structure = "Diagonal heterogeneous variance over a factor",
      spelling = "`diag(f):g`, `idh(f):g`, `at(f):g`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "covariance-family",
      condition = paste(
        "Repeated identical fits on a low-signal trial do not agree: three",
        "runs on agridat::besag.met span a factor of 2.3 in a component.",
        "Mode adoption (R/inla_restart_probe.R) removed the collapse-to-zero",
        "failure this row used to show, and did not remove the remaining",
        "spread. ASReml is stable on the same rows."
      ),
      oracle = paste(
        "ASReml 4.2.0.392 at this anchor on agridat::besag.met. Structural",
        "layer exact: 7 variance parameters, 6 estimable fixed effects, 1,152",
        "observations, identical on both sides. Numerical layer FAILS its",
        "declared tolerance on the intercept -- 148.60 against 149.14, a",
        "difference of 0.55 against a tolerance of 0.35. ASReml pins two of",
        "the six genotype variances at exactly zero and the posterior does",
        "not, so the two systems weight the fixed effects differently. The",
        "discrepancy is recorded rather than excused, and the tolerance is",
        "not moved"
      ),
      anchor = "test-heterogeneous-variance-inla.R"
    ),
    list(
      structure = "Sectioned residual, one variance per section",
      spelling = "`residual = ~ dsum(~ units | env)`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "exact",
      condition = "Inner term is plain `units`; a structured inner is refused.",
      oracle = paste(
        "ASReml 4.2.0.392 at this anchor on agridat::yan.winterwheat.",
        "Structural layer exact: 10 variance parameters, 9 estimable fixed",
        "effects, 162 observations, identical on both sides. Genotype",
        "variance ratio 0.0843 against 0.0854, a difference of 0.0012",
        "against a declared tolerance of 0.05; intercept identical to four",
        "decimal places; BLUP agreement r = 1.0000. Repeatability screen",
        "stable. This row is validated at the current anchor"
      ),
      anchor = "test-sectioned-residual-inla.R"
    ),
    list(
      structure = "Sectioned residual, ASReml 3 spelling",
      spelling = "`residual = ~ at(env):units`",
      vintage = "ASReml 3 alias",
      engine = "inla, brms",
      class = "exact",
      condition = paste(
        "ASReml 4 removed residual-side at() in favour of dsum(), so current",
        "ASReml may reject this spelling. Accepted here for compatibility;",
        "dsum() is the canonical form."
      ),
      oracle = "Shares the dsum() route and its evidence",
      anchor = "test-sectioned-residual-inla.R"
    ),
    list(
      structure = "Unstructured covariance over a factor, on brms",
      spelling = "`us(f):g`",
      vintage = "ASReml 4",
      engine = "brms",
      class = "covariance-family",
      condition = "",
      oracle = "Live emission and parameter-count tests; no current ASReml comparison",
      anchor = "test-unstructured-gxe-inla.R"
    ),
    list(
      structure = "Unstructured covariance over a factor, on INLA",
      spelling = "`us(f):g`",
      vintage = "ASReml 4",
      engine = "inla",
      class = "covariance-family",
      condition = paste(
        "Lowered to an iidKd field carrying a joint Wishart prior scaled to",
        "the data. The prior geometry and the parameter coordinates are not",
        "ASReml's; a supplied SD prior is refused rather than dropped."
      ),
      oracle = "Simulated recovery and Wishart-scale tests; none at this anchor against ASReml",
      anchor = "test-wishart-scale.R"
    ),
    list(
      structure = "Factor-analytic covariance",
      spelling = "`fa(f, k):g`",
      vintage = "ASReml 4",
      engine = "inla",
      class = "conditional",
      condition = paste(
        "Loadings are anchored to one for factor j at level j, taken from",
        "factor-level order. This spans the factor-analytic family only where",
        "the leading k-by-k loading block supports the anchoring: a covariance",
        "in which the anchor level has zero common-factor loading cannot be",
        "represented with finite values. There is no pivot selection and the",
        "precondition is not enforced."
      ),
      oracle = paste(
        "Simulated recovery with favourable anchors; one stored ASReml",
        "illustration at k = 1. Not compared at this anchor: three identical",
        "fits span a factor of 15 in a component, so a single paired",
        "comparison would measure the draw rather than the translation"
      ),
      anchor = "test-factor-analytic-inla.R"
    ),
    list(
      structure = "Known covariance or precision carrier",
      spelling = "`vm(g, K)`, `vm(g, cov = fb_cov(...))`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "exact",
      condition = "Representation metadata (dense / chol / precision / blocks) must be correct.",
      oracle = "Dense and sparse carrier tests; no current ASReml fixture for the inverse form",
      anchor = "test-known-covariance-inputs.R"
    ),
    list(
      structure = "Pedigree relationship carrier",
      spelling = "`ped(g, A)`",
      vintage = "ASReml 3 alias",
      engine = "inla, brms",
      class = "exact",
      condition = paste(
        "ASReml 4 removed ped() in favour of vm(). Accepted here for",
        "compatibility; vm(g, A) is the current spelling."
      ),
      oracle = "Shares the vm() carrier machinery and its evidence",
      anchor = "test-known-covariance-inputs.R"
    ),
    list(
      structure = "Correlated random slope",
      spelling = "`(x | g)`",
      vintage = "lme4 / brms",
      engine = "inla, brms",
      class = "covariance-family",
      condition = paste(
        "brms fits this natively. On INLA it is an iid2d block whose Wishart",
        "scale is set from the data; repeated identical fits have returned",
        "materially different modes where the correlation is weakly",
        "identified."
      ),
      oracle = "Live INLA recovery tests against lme4; none against ASReml us(2):g",
      anchor = "test-correlated-random-slope.R"
    ),
    list(
      structure = "Uncorrelated random slope",
      spelling = "`(x || g)`",
      vintage = "lme4 / brms",
      engine = "inla, brms",
      class = "covariance-family",
      condition = "",
      oracle = "lme4 REML at this anchor",
      anchor = "test-random-slopes-uncor.R"
    ),
    list(
      structure = "Temporal AR1 on the random side",
      spelling = "`~ ar1(t)`, grouped AR1",
      vintage = "ASReml 4",
      engine = "inla",
      class = "exact",
      condition = paste(
        "Exact for the latent-field model only: an AR1 field plus an",
        "observation residual. This is not a translation of ASReml's",
        "residual-side AR1, which carries no nugget."
      ),
      oracle = "Independent dense-V REML and simulation at this anchor",
      anchor = "test-inla-spatial-ar1.R"
    ),
    list(
      structure = "Separable residual AR1 x AR1",
      spelling = "`residual = ~ ar1(row):ar1(col)`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "refused",
      condition = paste(
        "The package's latent field plus nugget is a four-parameter",
        "neighbour, not ASReml's nugget-free residual, so it refuses rather",
        "than substituting."
      ),
      oracle = "Refusal is the claim; asserted by name",
      anchor = "test-inla-spatial-ar1.R"
    ),
    list(
      structure = "Penalised spline",
      spelling = "`spl(x)`, `s(x)`",
      vintage = "flexyBayes extension",
      engine = "inla",
      class = "analogue",
      condition = paste(
        "Lowered to INLA's second-order random walk. ASReml's spline carries",
        "its own knot count, duplicate-point handling, prediction points and",
        "spline.scale, and no comparison of basis, null space, penalty,",
        "effective degrees of freedom or prediction has been made. Related,",
        "not the same model."
      ),
      oracle = "Emission and fit tests only",
      anchor = "test-smooth.R"
    ),
    list(
      structure = "Gaussian relative weights",
      spelling = "`weights = w`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "exact",
      condition = paste(
        "Relative weights only, giving variance proportional to sigma^2 / w.",
        "ASReml's absolute-weight / fixed-dispersion mode is a different",
        "estimand and is not carried."
      ),
      oracle = "lme4 at this anchor",
      anchor = "test-weights-lowering.R"
    ),
    list(
      structure = "Missing response retained as a prediction target",
      spelling = "`na_action = \"augment\"`",
      vintage = "flexyBayes extension",
      engine = "inla, brms",
      class = "covariance-family",
      condition = paste(
        "A representation-preserving policy rather than a universal ASReml",
        "equivalence: the design row is kept as a latent target. The",
        "estimator assumptions are the user's to check."
      ),
      oracle = "Large missing-plot simulation study and independent dense-V REML",
      anchor = "test-na-action-augment.R"
    ),
    list(
      structure = "Heterogeneous variance with one common correlation",
      spelling = "`corh(f):g`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "refused",
      condition = "Neither active engine represents it; it is not the diagonal or the unstructured form.",
      oracle = "Refusal is the claim; asserted by name",
      anchor = "test-grammar-surface-refusals.R"
    ),
    list(
      structure = "Tensor-product and multivariate smooths",
      spelling = "`te(x, z)`, `ti()`, `t2()`, `s(x, z)`",
      vintage = "lme4 / brms",
      engine = "inla, brms",
      class = "refused",
      condition = "Univariate smooths only.",
      oracle = "Refusal is the claim; asserted by name",
      anchor = "test-grammar-contract.R"
    ),
    list(
      structure = "Structured inner term inside a sectioned residual",
      spelling = "`dsum(~ ar1(a) | b)`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "refused",
      condition = "Reducing it to a per-section variance would drop the structure written.",
      oracle = "Refusal is the claim; asserted by name",
      anchor = "test-grammar-contract.R"
    )
  )

  out <- do.call(
    rbind,
    lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE))
  )
  out$condition <- .fb_squish(out$condition)
  out$note <- NULL
  rownames(out) <- NULL
  out
}


.FB_TRANSLATION_MARKER_BEGIN <- "<!-- translation-registry:begin -->"
.FB_TRANSLATION_MARKER_END <- "<!-- translation-registry:end -->"


# Render the registry as a Markdown table.
#
# @noRd
.fb_translation_markdown <- function(markers = TRUE) {
  reg <- .fb_translation_registry()
  cell <- function(x) .fb_md_escape(.fb_squish(x))
  header <- c(
    "| Requested structure | flexyBayes spelling | Syntax vintage | Engine | Translation class | Condition | Independent oracle |",
    "|---|---|---|---|---|---|---|"
  )
  body <- vapply(
    seq_len(nrow(reg)),
    function(i) {
      sprintf(
        "| %s | %s | %s | %s | **%s** | %s | %s |",
        cell(reg$structure[[i]]),
        cell(reg$spelling[[i]]),
        cell(reg$vintage[[i]]),
        cell(reg$engine[[i]]),
        cell(reg$class[[i]]),
        if (nzchar(reg$condition[[i]])) cell(reg$condition[[i]]) else "--",
        cell(reg$oracle[[i]])
      )
    },
    character(1)
  )
  out <- c(header, body)
  if (markers) {
    out <- c(.FB_TRANSLATION_MARKER_BEGIN, out, .FB_TRANSLATION_MARKER_END)
  }
  out
}


# Splice the rendered registry into a surface carrying the markers.
#
# @noRd
.fb_translation_splice <- function(path, write = TRUE) {
  lines <- readLines(path, warn = FALSE)
  b <- which(lines == .FB_TRANSLATION_MARKER_BEGIN)
  e <- which(lines == .FB_TRANSLATION_MARKER_END)
  if (length(b) != 1L || length(e) != 1L || e <= b) {
    stop(
      "Translation-registry markers not found (or malformed) in ",
      path,
      call. = FALSE
    )
  }
  out <- c(
    lines[seq_len(b - 1L)],
    .fb_translation_markdown(markers = TRUE),
    lines[seq.int(e + 1L, length(lines))]
  )
  if (write) {
    writeLines(out, path)
  }
  invisible(out)
}
