# flexyBayes: Bayesian Mixed Models with ASReml Syntax via INLA and brms

`flexyBayes` lets you fit Bayesian mixed models using the formula syntax
you already know from ASReml or `lme4`/`brms`, dispatched through one of
two active inference engines: INLA (integrated nested Laplace
approximation, over the latent Gaussian model class) or brms (a Stan
passthrough, full Hamiltonian Monte Carlo). A `backend` request naming
any other engine refuses with an unknown-backend error naming the two
active engines. All current exports are at the experimental `lifecycle`
stage. The same fitted object supports
[`summary()`](https://rdrr.io/r/base/summary.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html),
`marginaleffects::predictions()`, and the `bayesplot::*` family.

## Entry points

The entry points share a single internal model representation
(`fb_terms`):

- [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md),
  ASReml-format entry: `fixed` / `random` / `residual` formulas and
  `known_matrices` for kinship / pedigree. Observation `weights` are
  lowered for the Gaussian family (identity link) and refused by name
  for any other family, see
  [`?flexybayes`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)'s
  `@param weights` for the exact semantics.

- [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  /
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md),
  the universal entry. Accepts an ASReml (`fixed` / `random` /
  `residual`) or brms-style (`y ~ x + (1 | g)`) formula, and any
  `backend` (`"inla"`, `"brms"`, or `"auto"`).

- [`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
  /
  [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md),
  the two active single-engine pins.

- [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md),
  penalised-complexity-canonical prior DSL.

- [`vpredict()`](https://aagi-aus.github.io/flexyBayes/reference/vpredict.md),
  ASReml's own verb for a posterior functional of the variance
  components: heritabilities, ratios and genetic correlations, computed
  on the draws, not by the delta method.

## Vignettes

Six vignettes ship with the package: *Mixed models in brief*, *Bayesian
estimation in brief*, *The ASReml syntax flexyBayes covers*, *ASReml and
flexyBayes on the same trials*, *A two-stage analysis of a series of
trials*, and *Moving from ASReml-R*. The first two set up the model
class and what changes under a posterior. The third is the formula
surface together with dispatch and the inventory of structures that are
refused. The fourth fits the same trials through both and shows where
the answers agree and where they do not. The fifth fits each trial of a
series on its own, carries the first stage's variances into a second
stage, and fits a penalised regression of trial-level covariates on the
design whitened by the fitted covariance. The sixth gives the flexyBayes
spelling of each ASReml construct, beside stored ASReml fits. The count
here is checked against `vignettes/` by a parity test, because a number
written into prose is the first thing to drift.

## Capability

What each active engine fits, emits, or refuses by model class is
generated from one R-level table and shown in `README.md` and
`system.file("KNOWN_ISSUES.md", package = "flexyBayes")`. A request
outside that set raises a typed refusal naming the nearest implemented
alternative instead of fitting a neighbouring model under the requested
model's name.

## References

Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H.
(2017). Penalising model component complexity: A principled, practical
approach to constructing priors. *Statistical Science*, 32(1), 1–28.

Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
inference for latent Gaussian models by using integrated nested Laplace
approximations. *Journal of the Royal Statistical Society: Series B*,
71(2), 319–392.

Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B.,
Yao, Y., Kennedy, L., Gabry, J., Bürkner, P.-C., & Modrák, M. (2020).
Bayesian workflow. *arXiv* 2011.01808.

## See also

Useful links:

- <https://aagi-aus.github.io/flexyBayes/>

- <https://github.com/AAGI-AUS/flexyBayes>

- Report bugs at <https://github.com/AAGI-AUS/flexyBayes/issues>

## Author

**Maintainer**: Max Moldovan <max.moldovan@adelaide.edu.au>
([ORCID](https://orcid.org/0000-0001-9680-8474))

Authors:

- Max Moldovan <max.moldovan@adelaide.edu.au>
  ([ORCID](https://orcid.org/0000-0001-9680-8474))

Other contributors:

- Emi Tanaka ([ORCID](https://orcid.org/0000-0002-1455-259X))
  \[contributor\]

- Francis K.C. Hui ([ORCID](https://orcid.org/0000-0003-0765-3533))
  \[contributor\]

- Anabel Forte Deltell ([ORCID](https://orcid.org/0000-0001-9534-1817))
  \[contributor\]

- Grains Research and Development Corporation
  ([ROR](https://ror.org/02xwr1996)) (GRDC Project CUR2210-005OPX
  (AAGI-CU), UOA2301-005OPX (AAGI-AU)) \[funder, copyright holder\]
