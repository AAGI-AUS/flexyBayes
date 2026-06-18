#' flexyBayes: Bayesian Mixed Models with ASReml Syntax via greta, INLA, and brms
#'
#' `flexyBayes` lets you fit Bayesian mixed models using the formula syntax
#' you already know from ASReml or `lme4`/`brms`, dispatched through one of
#' three inference engines: greta (Hamiltonian Monte Carlo via TensorFlow),
#' INLA (integrated nested Laplace approximation, for the latent Gaussian
#' model class), or brms (a Stan passthrough). All current exports are at the
#' experimental `lifecycle` stage. The same fitted object supports
#' `summary()`, `predict()`, `emmeans::emmeans()`,
#' `marginaleffects::predictions()`, and the `bayesplot::*` family.
#'
#' The package's signature feature is `triangulate()`, a cross-engine
#' posterior comparison that quantifies disagreement between two fits of
#' the same model on the same data.
#'
#' @section Entry points:
#' Two ingest paths share a single internal model representation
#' (`fb_terms`):
#'
#' * [flexybayes()] — asreml-format entry: `fixed` / `random` / `rcov`
#'   formulas, `known_matrices` for kinship / pedigree, `weights` for
#'   pre-aggregated observations.
#' * [fb()] / [flexybayes()] — the universal entry. Accepts an ASReml
#'   (`fixed` / `random` / `rcov`) or brms-style (`y ~ x + (1 | g)`)
#'   formula, or a native `greta::model()`, and any `backend`
#'   (`"greta"`, `"inla"`, `"brms"`, or `"auto"`).
#' * [fb_greta()] / [fb_inla()] / [fb_brms()] — single-engine pins.
#' * [fb_prior()] — penalised-complexity-canonical prior DSL.
#' * [triangulate()] — cross-engine posterior comparison.
#'
#' @section Vignettes:
#' Sixteen vignettes ship with the package, covering: getting started,
#' asreml-shaped formulas reference, foundational regression,
#' hierarchical models, structured covariance, priors and
#' regularisation, multi-environment trials and genomics, downstream
#' analysis, spatio-temporal models, cross-engine triangulation, LGM
#' feasibility, backend internals, LGM feasibility and memory, choosing
#' an engine (the universal entry and the engine pins), extending the
#' backend registry, and big-data streaming (exact aggregation).
#'
#' @section References:
#' Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H.
#' (2017). Penalising model component complexity: A principled,
#' practical approach to constructing priors. *Statistical Science*,
#' 32(1), 1–28.
#'
#' Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
#' inference for latent Gaussian models by using integrated nested
#' Laplace approximations. *Journal of the Royal Statistical Society:
#' Series B*, 71(2), 319–392.
#'
#' Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter,
#' B., Yao, Y., Kennedy, L., Gabry, J., Bürkner, P.-C., & Modrák, M.
#' (2020). Bayesian workflow. *arXiv* 2011.01808.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end
NULL
