#' flexyBayes: Bayesian Mixed Models with ASReml Syntax via INLA and brms
#'
#' `flexyBayes` lets you fit Bayesian mixed models using the formula syntax
#' you already know from ASReml or `lme4`/`brms`, dispatched through one of
#' two active inference engines: INLA (integrated nested Laplace
#' approximation, over the latent Gaussian model class) or brms (a Stan
#' passthrough, full Hamiltonian Monte Carlo). A third native engine was
#' withdrawn entirely in 0.9.3 (see `NEWS.md`): a `backend` request
#' naming that withdrawn engine, or any other unrecognised name, now
#' refuses with an ordinary unknown-backend error naming the two active
#' engines; re-entry, should it ever be proposed, would be a fresh
#' implementation, not a repair of retained code. All current exports
#' are at the experimental `lifecycle` stage. The same fitted object
#' supports `summary()`, `predict()`, `emmeans::emmeans()`,
#' `marginaleffects::predictions()`, and the `bayesplot::*` family.
#'
#' @section Entry points:
#' The entry points share a single internal model representation
#' (`fb_terms`):
#'
#' * [flexybayes()] -- ASReml-format entry: `fixed` / `random` / `residual`
#'   formulas and `known_matrices` for kinship / pedigree. Observation
#'   `weights` are lowered for the Gaussian family (identity link) and
#'   refused by name for any other family -- see `?flexybayes`'s
#'   `@param weights` for the exact semantics.
#' * [fb()] / [flexybayes()] -- the universal entry. Accepts an ASReml
#'   (`fixed` / `random` / `residual`) or brms-style (`y ~ x + (1 | g)`)
#'   formula, and any `backend` (`"inla"`, `"brms"`, or `"auto"`).
#' * [fb_inla()] / [fb_brms()] -- the two active single-engine pins.
#' * [fb_prior()] -- penalised-complexity-canonical prior DSL.
#' * [vpredict()] -- ASReml's own verb for a posterior functional of the
#'   variance components: heritabilities, ratios and genetic correlations,
#'   computed on the draws rather than by the delta method.
#'
#' @section Vignettes:
#' Six vignettes ship with the package: *Mixed models in brief*, *Bayesian
#' estimation in brief*, *The ASReml syntax flexyBayes covers*, *ASReml and
#' flexyBayes on the same trials*, *A two-stage analysis of a series of
#' trials*, and *Moving from ASReml-R*. The first two set up the model class
#' and what changes under a posterior; the third is the formula surface
#' together with dispatch and the inventory of structures that are refused;
#' the fourth fits the same trials through both and shows where the answers
#' agree and where they do not; the fifth fits each trial of a series on its
#' own, carries the first stage's variances into a second stage, and fits a
#' penalised regression of trial-level covariates on the design whitened by
#' the fitted covariance; the sixth gives the flexyBayes spelling of each
#' ASReml construct, beside stored ASReml fits. The count here is checked
#' against `vignettes/` by a parity test, because a number written into prose
#' is the first thing to drift.
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
#' 32(1), 1--28.
#'
#' Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
#' inference for latent Gaussian models by using integrated nested
#' Laplace approximations. *Journal of the Royal Statistical Society:
#' Series B*, 71(2), 319--392.
#'
#' Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter,
#' B., Yao, Y., Kennedy, L., Gabry, J., Bürkner, P.-C., & Modrák, M.
#' (2020). Bayesian workflow. *arXiv* 2011.01808.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
