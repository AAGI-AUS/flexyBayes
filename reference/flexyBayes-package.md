# flexyBayes: Bayesian Mixed Models with ASReml Syntax via greta, INLA, brms

`flexyBayes` lets you fit Bayesian mixed models using the formula syntax
you already know from ASReml or `lme4`/`brms`, dispatched through one of
three inference engines: greta (Hamiltonian Monte Carlo via TensorFlow),
INLA (integrated nested Laplace approximation, for the latent Gaussian
model class), or brms (a Stan passthrough). All current exports are at
the experimental `lifecycle` stage. The same fitted object supports
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

Two ingest paths share a single internal model representation
(`fb_terms`):

- [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  — asreml-format entry: `fixed` / `random` / `rcov` formulas,
  `known_matrices` for kinship / pedigree, `weights` for pre-aggregated
  observations.

- [`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  /
  [`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
  — the universal entry. Accepts an ASReml (`fixed` / `random` / `rcov`)
  or brms-style (`y ~ x + (1 | g)`) formula, or a native
  `greta::model()`, and any `backend` (`"greta"`, `"inla"`, `"brms"`, or
  `"auto"`).

- [`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md)
  /
  [`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
  /
  [`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md)
  — single-engine pins.

- [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
  — penalised-complexity-canonical prior DSL.

- [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
  — cross-engine posterior comparison.

## Vignettes

Sixteen vignettes ship with the package, covering: getting started,
asreml-shaped formulas reference, foundational regression, hierarchical
models, structured covariance, priors and regularisation,
multi-environment trials and genomics, downstream analysis,
spatio-temporal models, cross-engine triangulation, LGM feasibility,
backend internals, LGM feasibility and memory, choosing an engine (the
universal entry and the engine pins), extending the backend registry,
and big-data streaming (exact aggregation).

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

- <https://github.com/AAGI-AUS/flexyBayes>

- Report bugs at <https://github.com/AAGI-AUS/flexyBayes/issues>

## Author

**Maintainer**: Max Moldovan <max.moldovan@adelaide.edu.au>
([ORCID](https://orcid.org/0000-0001-9680-8474))

Authors:

- Emi Tanaka ([ORCID](https://orcid.org/0000-0002-1455-259X))

- Francis K.C. Hui ([ORCID](https://orcid.org/0000-0003-0765-3533))

- Anabel Forte Deltell ([ORCID](https://orcid.org/0000-0001-9534-1817))
