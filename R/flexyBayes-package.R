#' flexyBayes: Bayesian Mixed Models with ASReml Syntax via INLA and brms
#'
#' `flexyBayes` lets you fit Bayesian mixed models using the formula syntax
#' you already know from ASReml or `lme4`/`brms`, dispatched through one of
#' two active inference engines: INLA (integrated nested Laplace
#' approximation, over the latent Gaussian model class) or brms (a Stan
#' passthrough, full Hamiltonian Monte Carlo). greta and gretaR are
#' quarantined as fitting engines --- their emit code and registry
#' descriptors are retained as re-entry candidates, `backend = "auto"`
#' never selects them, and an explicit request refuses with
#' `backend_quarantined`. All current exports are at the experimental
#' `lifecycle` stage. The same fitted object supports `summary()`,
#' `predict()`, `emmeans::emmeans()`, `marginaleffects::predictions()`, and
#' the `bayesplot::*` family.
#'
#' The package's signature feature is `triangulate()`, a cross-engine
#' posterior comparison that quantifies disagreement between two fits of
#' the same model on the same data.
#'
#' @section Entry points:
#' Two ingest paths share a single internal model representation
#' (`fb_terms`):
#'
#' * [flexybayes()] — asreml-format entry: `fixed` / `random` / `residual`
#'   formulas and `known_matrices` for kinship / pedigree. Observation
#'   `weights` are parsed and then refused (`weights_not_supported`) until
#'   an active emitter consumes them.
#' * [fb()] / [flexybayes()] — the universal entry. Accepts an ASReml
#'   (`fixed` / `random` / `residual`) or brms-style (`y ~ x + (1 | g)`)
#'   formula, and any `backend` (`"inla"`, `"brms"`, or `"auto"`).
#' * [fb_inla()] / [fb_brms()] — single-engine pins. [fb_greta()] is
#'   retained and refuses.
#' * [fb_prior()] — penalised-complexity-canonical prior DSL.
#' * [triangulate()] — cross-engine posterior comparison.
#'
#' @section Vignettes:
#' Eleven vignettes ship with the package, covering: getting started,
#' the formula surface (the asreml term catalogue together with
#' structured covariance), foundational regression, hierarchical models,
#' priors and regularisation, multi-environment trials and genomics,
#' downstream analysis, spatio-temporal models, cross-engine
#' triangulation, dispatch and refusals with the backend registry, and
#' big-data streaming (exact aggregation). The dispatch-and-refusals
#' page is the technical / internals reference; the rest target a
#' general audience.
#'
#' @section Capability:
#' What each active engine fits, emits, or refuses by model class is
#' generated from one R-level table and shown in `README.md` and
#' `system.file("KNOWN_ISSUES.md", package = "flexyBayes")`. A request
#' outside that set raises a typed refusal naming the nearest implemented
#' alternative rather than fitting a neighbouring model under the
#' requested model's name.
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
