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
#              parsed representation, so a common parser error agrees
#              perfectly across them.
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

# How many engine slots a prior on each structure reaches. The head word
# before a colon is one of .FB_TRANSLATION_CARDINALITIES, and the model
# manifest records the same count for each fit (R/model_manifest.R).
.FB_TRANSLATION_CARDINALITIES <- c(
  "one-to-one",
  "one-to-many",
  "engine default",
  "refused"
)

.FB_TRANSLATION_CARDINALITY <- c(
  "Fixed factors, covariates, means model" =
    "one-to-one: a b() prior reaches its own coefficient",
  "Simple iid random intercept" = "one-to-one",
  "Nested and interacted iid terms" =
    "one-to-one: keyed by the colon-joined group",
  "Diagonal heterogeneous variance over a factor" =
    "one-to-many: one field per level of the outer factor",
  "Sectioned residual, one variance per section" =
    "one-to-many: one likelihood per section",
  "Sectioned residual, ASReml 3 spelling" =
    "one-to-many: one likelihood per section",
  "Per-trial genotype and residual variances" =
    "one-to-many: one field and one likelihood per trial",
  "Unstructured covariance over a factor, on brms" =
    "one-to-one: one standard-deviation row for the group",
  "Unstructured covariance over a factor, on INLA" = "refused",
  "Factor-analytic covariance" =
    "one-to-many: one field per common factor and one per level",
  "Known covariance or precision carrier" =
    "one-to-one: one field; a blocks carrier is refused on both engines",
  "Pedigree relationship carrier" = "one-to-one",
  "Correlated random slope, on brms" =
    "one-to-one: the intercept and slope standard deviations each take a row",
  "Correlated random slope, on INLA" = "refused",
  "Uncorrelated random slope" =
    "one-to-one: the intercept and slope standard deviations each take a field",
  "Temporal AR1 on the random side" =
    "engine default: the field's precision and correlation keep INLA's priors",
  "Separable residual AR1 x AR1" =
    "engine default: the field's precision and correlations keep INLA's priors",
  "Random term on selected levels of a factor" =
    "one-to-many: one field per named level",
  "Separable AR1 field per level of a factor" =
    "one-to-many: one field per level, each with its own hyperparameters",
  "Penalised spline" = "one-to-one: smooth() on the random walk's precision",
  "Gaussian relative weights" = "one-to-one: the residual prior is unchanged",
  "Missing response retained as a prediction target" =
    "one-to-one: the priors are unchanged",
  "Prediction over a random classify factor" =
    "one-to-one: the priors are unchanged",
  "Prediction at a level the fit did not see" =
    "one-to-one: the priors are unchanged",
  "Heterogeneous variance with one common correlation" = "refused",
  "Tensor-product and multivariate smooths" = "refused",
  "Structured inner term inside a sectioned residual" = "refused",
  "Separable residual AR1 x AR1 per section" =
    "one-to-many: one field per section, each with its own hyperparameters"
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
#'   `condition`, `oracle`, `anchor`, and `cardinality`.
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
      structure = "Per-trial genotype and residual variances",
      spelling = paste(
        "`random = ~ diag(trial):g`,",
        "`residual = ~ dsum(~ units | trial)`"
      ),
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "covariance-family",
      condition = paste(
        "Past 15 variance parameters INLA fits the model by empirical Bayes",
        "(int.strategy = \"eb\"): the fixed and random effects are",
        "conditional on the posterior mode of the variance parameters, and",
        "each variance keeps its own posterior. Other models past 15 refuse",
        "with hyperparam_budget."
      ),
      oracle = paste(
        "ASReml 4.2.0.482 at this anchor on 12 simulated trials, 26 variance",
        "parameters (data-raw/model_c_reference.R): every per-trial variance",
        "ASReml identifies within 0.17 of its ASReml standard error on INLA",
        "and 0.28 on brms"
      ),
      anchor = "test-model-c-per-trial-variances.R"
    ),
    list(
      structure = "Unstructured covariance over a factor, on brms",
      spelling = "`us(f):g`",
      vintage = "ASReml 4",
      engine = "brms",
      class = "covariance-family",
      condition = "",
      oracle =
        "Live emission and parameter-count tests; no current ASReml comparison",
      anchor = "test-unstructured-gxe-inla.R"
    ),
    list(
      structure = "Unstructured covariance over a factor, on INLA",
      spelling = "`us(f):g`",
      vintage = "ASReml 4",
      engine = "inla",
      class = "refused",
      condition = paste(
        "Estimated correlation matrices are fitted on brms only. The INLA",
        "lowering, an iidKd field under a joint Wishart prior, refuses with",
        "correlated_covariance_requires_brms and backend = \"auto\" routes to",
        "brms. The row was covariance-family until 2026-09-11, analogue until",
        "2026-09-13, and refused from then."
      ),
      oracle = paste(
        "Refusal is the claim; asserted by name. The refusal rests on ASReml",
        "4.2.0.392 on agridat::besag.met, oracle scenarios 06 and 07 (finding",
        "O04): a genetic correlation of 0.142 against 0.930 over two",
        "environments and 0.198 / 0.231 / 0.320 against 0.898 / 0.915 / 0.955",
        "over three, with eight identical calls on the response rescaled by",
        "ten none converged and correlations from -0.894 to +0.686. Under",
        "every identical Wishart prior tried INLA and Stan agreed, and no",
        "Wishart scale gave REML-sized variances (INLA scale diagnosis,",
        "2026-09-13), so the joint prior rather than the integration drives",
        "the disagreement"
      ),
      anchor = "test-inla-containment.R"
    ),
    list(
      structure = "Factor-analytic covariance",
      spelling = "`fa(f, k):g`",
      vintage = "ASReml 4",
      engine = "inla",
      class = "refused",
      condition = paste(
        "Refused from 0.13.0 with factor_analytic_not_repeatable, on both",
        "engines, and flexyBayes has no factor-analytic lowering for brms.",
        "The loadings are anchored to one for factor j at an anchor level,",
        "which represents the factor-analytic family only where that level",
        "carries a common-factor loading. From 0.17.0 the anchors are the",
        "levels with the largest common genetic variance in the data; until",
        "0.16.0 they were the first levels in factor order. The row was",
        "covariance-family until",
        "2026-09-11, conditional from then, and refused from 2026-09-16."
      ),
      oracle = paste(
        "Refusal is the claim; asserted by name. The refusal rests on INLA",
        "fits of 20 simulated trials with a known covariance, and on",
        "agridat::omer.sorghum against ASReml 4.2; the numbers are under",
        "Reading a row"
      ),
      anchor = "test-factor-analytic-inla.R"
    ),
    list(
      structure = "Known covariance or precision carrier",
      spelling = "`vm(g, K)`, `vm(g, cov = fb_cov(...))`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "exact",
      condition = paste0(
        "Representation metadata (dense / chol / precision / blocks) must be ",
        "correct."
      ),
      oracle = paste0(
        "Dense and sparse carrier tests; no current ASReml fixture for the ",
        "inverse form"
      ),
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
      structure = "Correlated random slope, on brms",
      spelling = "`(x | g)`",
      vintage = "lme4 / brms",
      engine = "brms",
      class = "covariance-family",
      condition = "brms writes the term verbatim and fits it natively.",
      oracle = paste(
        "Live emission and recovery tests; no current ASReml comparison on",
        "this engine"
      ),
      anchor = "test-correlated-random-slope.R"
    ),
    list(
      structure = "Correlated random slope, on INLA",
      spelling = "`(x | g)`",
      vintage = "lme4 / brms",
      engine = "inla",
      class = "refused",
      condition = paste(
        "Estimated correlation matrices are fitted on brms only. The INLA",
        "lowering, one iid2d field over a two-block index under a joint",
        "Wishart prior, refuses with correlated_covariance_requires_brms and",
        "backend = \"auto\" routes to brms. The row was covariance-family",
        "until 2026-09-13 and is refused from then."
      ),
      oracle = paste(
        "Refusal is the claim; asserted by name. The refusal rests on ASReml",
        "4.2.0.392 and lme4 on lme4::sleepstudy, oracle row 19 (finding S01,",
        "defect D1): a slope variance of 46.37 at the posterior mode against",
        "35.076 and 35.072, and a tenfold rescaling of the response that",
        "moved the intercept-slope correlation from 0.02 to 0.969 while REML",
        "held at 0.066. Under an identical Wishart prior INLA matched Stan",
        "within 1.5 to 6.5 per cent, and a simulated response on a kilogram",
        "scale returned a residual standard deviation of 0.13 against a truth",
        "of 250 (INLA scale diagnosis, 2026-09-13)"
      ),
      anchor = "test-inla-containment.R"
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
      engine = "inla",
      class = "covariance-family",
      condition = paste(
        "Gaussian response with the identity link, one plot per grid cell.",
        "The observation error is fixed at a variance of exp(-12) on the",
        "standardised response; the posterior medians move by less than",
        "1e-4 between exp(-10) and exp(-16). brms refuses."
      ),
      oracle = paste(
        "ASReml 4.2.0.482 at this anchor on a simulated 16 x 10 trial",
        "(data-raw/at_ar1_reference.R): every posterior median within 0.3",
        "posterior SDs of the REML estimate"
      ),
      anchor = "test-at-fields.R"
    ),
    list(
      structure = "Separable residual AR1 x AR1 per section",
      spelling = "`residual = ~ dsum(~ ar1(row):ar1(col) | trial)`",
      vintage = "ASReml 4",
      engine = "inla",
      class = "covariance-family",
      condition = paste(
        "As for the single residual field, on each section's own grid.",
        "brms refuses."
      ),
      oracle = paste(
        "ASReml 4.2.0.482 at this anchor on three simulated trials",
        "(data-raw/at_ar1_reference.R): every posterior median within 0.3",
        "posterior SDs of the REML estimate"
      ),
      anchor = "test-at-fields.R"
    ),
    list(
      structure = "Random term on selected levels of a factor",
      spelling = "`at(f, level):g`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "exact",
      condition = paste(
        "g is a plain factor. Levels are named by value or by position in",
        "levels(f); a level the data lack refuses."
      ),
      oracle = paste(
        "ASReml 4.2.0.482 at this anchor on six simulated trials, three with",
        "main plots (data-raw/at_ar1_reference.R), on both engines: every",
        "posterior median within 0.3 posterior SDs of the REML estimate"
      ),
      anchor = "test-at-fields.R"
    ),
    list(
      structure = "Separable AR1 field per level of a factor",
      spelling = paste(
        "`at(trial):ar1(row):ar1(col)`,",
        "`at(trial, level):ar1(row):ar1(col)`"
      ),
      vintage = "ASReml 4",
      engine = "inla",
      class = "exact",
      condition = paste(
        "An AR1 field with its own SD and correlations on each named level,",
        "plus the observation residual; each level's plots fill its grid.",
        "brms refuses."
      ),
      oracle = paste(
        "ASReml 4.2.0.482 at this anchor on three simulated trials",
        "(data-raw/at_ar1_reference.R): every posterior median within 0.3",
        "posterior SDs of the REML estimate"
      ),
      anchor = "test-at-fields.R"
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
      oracle =
        "Large missing-plot simulation study and independent dense-V REML",
      anchor = "test-na-action-augment.R"
    ),
    list(
      structure = "Prediction over a random classify factor",
      spelling = "`predict(fit, classify = \"gen\")` with `gen` random",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "exact",
      condition = paste(
        "ASReml's default prediction rule: the fixed part averaged over the",
        "other fixed factors, plus the random effects of the classify terms;",
        "other random terms are ignored. present, average, ignore, sed and",
        "vcov refuse with classify_argument_not_supported."
      ),
      oracle = paste(
        "ASReml 4.2.0.482 at this anchor on agridat::john.alpha and 16",
        "genotypes of agridat::besag.met (data-raw/asreml_classify_random.R):",
        "within 0.11 SE for genotypes and 0.4 SE for genotype-by-county",
        "cells, and within 0.006 SE with the variances held at ASReml's"
      ),
      anchor = "test-predict-classify-random.R"
    ),
    list(
      structure = "Prediction at a level the fit did not see",
      spelling = paste(
        "`predict(fit, newdata, allow_new_levels = TRUE,",
        "mean_by = ~ trial)`"
      ),
      vintage = "flexyBayes extension",
      engine = "inla, brms",
      class = "covariance-family",
      condition = paste(
        "A posterior predictive, with new levels drawn from their fitted",
        "variances and the observation noise included. Random slopes,",
        "splines, known covariances, factor-analytic and unstructured",
        "terms, observation weights, a known variance and non-Gaussian",
        "responses refuse by name."
      ),
      oracle = paste(
        "ASReml 4.2.0.482 at this anchor",
        "(data-raw/asreml_predict_new_level.R):",
        "at ASReml's REML variances the variance of a trial, cell or",
        "locality mean agrees with ASReml's plug-in value to 2e-15"
      ),
      anchor = "test-predict-new-levels-asreml.R"
    ),
    list(
      structure = "Heterogeneous variance with one common correlation",
      spelling = "`corh(f):g`",
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "refused",
      condition = paste(
        "flexyBayes has no lowering for it on either engine; it is not the",
        "diagonal or the unstructured form. INLA has an exchangeable group",
        "model, with one variance and one shared correlation, which",
        "flexyBayes does not use."
      ),
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
      spelling = paste(
        "`dsum(~ ar1(a) | b)`, any inner term but `units` and",
        "`ar1(a):ar1(b)`"
      ),
      vintage = "ASReml 4",
      engine = "inla, brms",
      class = "refused",
      condition = paste0(
        "Reducing it to a per-section variance would drop the structure ",
        "written."
      ),
      oracle = "Refusal is the claim; asserted by name",
      anchor = "test-grammar-contract.R"
    )
  )

  out <- do.call(
    rbind,
    lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE))
  )
  out$condition <- .fb_squish(out$condition)
  out$cardinality <- unname(.FB_TRANSLATION_CARDINALITY[out$structure])
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
    paste0(
      "| Requested structure | flexyBayes spelling | Syntax vintage | Engine ",
      "| Translation class | Prior cardinality | Condition | Independent ",
      "oracle |"
    ),
    "|---|---|---|---|---|---|---|---|"
  )
  body <- vapply(
    seq_len(nrow(reg)),
    function(i) {
      sprintf(
        "| %s | %s | %s | %s | **%s** | %s | %s | %s |",
        cell(reg$structure[[i]]),
        cell(reg$spelling[[i]]),
        cell(reg$vintage[[i]]),
        cell(reg$engine[[i]]),
        cell(reg$class[[i]]),
        cell(reg$cardinality[[i]]),
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
