# flexyBayes: Bayesian Mixed Models with ASReml Syntax via INLA and brms

`flexyBayes` lets you fit Bayesian mixed models using the formula syntax
you already know from ASReml or `lme4`/`brms`, dispatched through one of
two active inference engines: INLA (integrated nested Laplace
approximation, over the latent Gaussian model class) or brms (a Stan
passthrough, full Hamiltonian Monte Carlo). A third native engine was
withdrawn entirely in 0.9.3 (see `NEWS.md`): a `backend` request naming
that withdrawn engine, or any other unrecognised name, now refuses with
an ordinary unknown-backend error naming the two active engines;
re-entry, should it ever be proposed, would be a fresh implementation,
not a repair of retained code. All current exports are at the
experimental `lifecycle` stage. The same fitted object supports
[`summary()`](https://rdrr.io/r/base/summary.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
`emmeans::emmeans()`, `marginaleffects::predictions()`, and the
`bayesplot::*` family.

## Details

The package's signature feature is
[`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md),
a cross-engine posterior comparison that quantifies disagreement between
two fits of the same model on the same data.

## Entry points

The entry points share a single internal model representation
(`fb_terms`):

- [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  – ASReml-format entry: `fixed` / `random` / `residual` formulas and
  `known_matrices` for kinship / pedigree. Observation `weights` are
  lowered for the Gaussian family (identity link) and refused by name
  for any other family – see
  [`?flexybayes`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)'s
  `@param weights` for the exact semantics.

- [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  /
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  – the universal entry. Accepts an ASReml (`fixed` / `random` /
  `residual`) or brms-style (`y ~ x + (1 | g)`) formula, and any
  `backend` (`"inla"`, `"brms"`, or `"auto"`).

- [`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
  /
  [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md)
  – the two active single-engine pins.

- [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
  – penalised-complexity-canonical prior DSL.

- [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
  – cross-engine posterior comparison.

## Vignettes

Eight vignettes ship with the package, covering: getting started and
what changes under a posterior, the formula surface together with
dispatch and the inventory of unsupported structures, regression and
hierarchical models, priors, multi-environment trials and genomics,
spatial and temporal structure, reading and comparing a fit including
cross-engine triangulation, and fitting from sufficient statistics when
the data do not fit in memory. Each opens with a panel giving the ASReml
line, the flexyBayes line, and what the posterior adds.

## Capability

What each active engine fits, emits, or refuses by model class is
generated from one R-level table and shown in `README.md` and
`system.file("KNOWN_ISSUES.md", package = "flexyBayes")`. A request
outside that set raises a typed refusal naming the nearest implemented
alternative rather than fitting a neighbouring model under the requested
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

- Emi Tanaka ([ORCID](https://orcid.org/0000-0002-1455-259X))

- Francis K.C. Hui ([ORCID](https://orcid.org/0000-0003-0765-3533))

- Anabel Forte Deltell ([ORCID](https://orcid.org/0000-0001-9534-1817))
