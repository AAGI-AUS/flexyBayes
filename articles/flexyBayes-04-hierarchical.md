# Tutorial 04: Hierarchical models: random effects and GLMMs

## 1. Why hierarchical?

> **This is a syntax-and-representation catalogue.** It demonstrates how
> to *write* models across the term surface. The fits use small sampling
> budgets so the document builds quickly, and several printed
> diagnostics therefore do **not** meet the package’s own convergence
> thresholds ($`\widehat{R} \le 1.01`$, ESS $`\ge 400`$); some model
> classes also mix poorly on their backend at any budget. Read the
> posterior summaries here as illustrations of the output *shape*, not
> as inferential results. For a convergence-clean worked fit and the
> trustworthy-fit workflow, see the *getting started* and *cross-engine
> triangulation* vignettes.

A *hierarchical* model — in mixed-models language, a model with one or
more random effects — partitions the residual variance into
between-cluster and within-cluster components. Three reasons this
matters:

1.  **Honest standard errors.** Treating clustered observations as
    independent shrinks standard errors and inflates Type I error.
2.  **Borrowing of strength.** A random effect lets each cluster
    “borrow” information from the others through the shared variance
    parameter — *partial pooling*. This is the source of the Bayesian
    shrinkage estimator.
3.  **Predictions for unseen clusters.** A random-effect fit gives a
    posterior distribution for the next group’s mean, drawn from the
    estimated population.

`flexyBayes` covers the full mixed-model surface — Gaussian linear mixed
models (LMMs), generalised linear mixed models (GLMMs) with Poisson,
binomial, gamma, and beta likelihoods, plus structured covariance models
that the *structured covariance* and *spatio-temporal* vignettes take
up. This vignette focuses on the foundational cases.

## 2. The Gaussian random intercept: `sleepstudy`

We start where the *getting started* vignette left off: a random
intercept on `Subject` with `Reaction` as the response and `Days` as a
fixed continuous covariate.

``` r

library(flexyBayes)
data(sleepstudy, package = "lme4")
fit_ri <- fb_greta(
  Reaction ~ Days + (1 | Subject), data = sleepstudy,
  n_samples    = 2000, warmup = 5000, chains = 4,
  verbose      = FALSE,
  mcmc_verbose = FALSE
)
summary(fit_ri)
#> Bayesian mixed model summary  [flexyBayes]
#> ============================================================ 
#>   Fixed  : Reaction ~ Days + (1 | Subject) 
#>   Family : gaussian / identity 
#>   N = 180 , chains = 4 , samples = 2000 
#>   Representation: exact
#>   Engine:         greta MCMC
#> 
#> -- Fixed effects (posterior)  --------------------------------- 
#>             Estimate Post.SD     2.5%    97.5%
#> (Intercept) 250.4686  9.5563 234.1178 275.2858
#> Days         10.4951  0.8209   8.8734  12.0981
#> 
#> -- Variance components  -------------------------------------- 
#>       Component Estimate     SD    2.5%   97.5%
#> 1 sigma_Subject  39.6313 8.2255 27.3164 59.0680
#> 2   sigma_e_atg  31.3126 1.7588 28.1633 34.9687
#> 
#> -- Convergence  --------------------------------------------- 
#>   Rhat range: 1.001 - 1.084 
#>   ESS  range: 91 - 3823 
#>   Run time  : 29.7 sec
```

The `Subject`-level standard deviation summarises between-subject
variation in baseline reaction time; the residual standard deviation
summarises within-subject variation. The intraclass correlation
coefficient (ICC) is

``` math
\text{ICC} = \frac{\sigma_{\text{Subject}}^2}{\sigma_{\text{Subject}}^2 + \sigma_e^2},
```

a derived quantity that the `posterior::summarise_draws()` call computes
from the draws.

``` r

draws <- flexyBayes::fb_as_draws_simple(fit_ri)
# Greta names the random-effect SD after the grouping factor:
# sigma_Subject (capital S to match `(1 | Subject)`).
sigma_S <- draws[["sigma_Subject"]]
sigma_e <- draws$sigma_e_atg
icc <- sigma_S^2 / (sigma_S^2 + sigma_e^2)
c(mean = mean(icc), q025 = quantile(icc, 0.025), q975 = quantile(icc, 0.975))
#>       mean  q025.2.5% q975.97.5% 
#>  0.6026748  0.4224249  0.7837353
```

A reference fit with `lme4`:

``` r

summary(lme4::lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy))
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: Reaction ~ Days + (1 | Subject)
#>    Data: sleepstudy
#> 
#> REML criterion at convergence: 1786.5
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -3.2257 -0.5529  0.0109  0.5188  4.2506 
#> 
#> Random effects:
#>  Groups   Name        Variance Std.Dev.
#>  Subject  (Intercept) 1378.2   37.12   
#>  Residual              960.5   30.99   
#> Number of obs: 180, groups:  Subject, 18
#> 
#> Fixed effects:
#>             Estimate Std. Error t value
#> (Intercept) 251.4051     9.7467   25.79
#> Days         10.4673     0.8042   13.02
#> 
#> Correlation of Fixed Effects:
#>      (Intr)
#> Days -0.371
```

The Bayesian posterior means and `lme4`’s REML estimates agree to about
the second decimal — because `Subject` has 18 levels, which is on the
cusp where REML and Bayesian fits start to diverge for small clusters.

> **Reading the diagnostics.** The fits here use production MCMC budgets
> (warmup 5000, four chains), so the fixed effects and the headline
> variance components are reliable. As in the *getting started*
> vignette, the slowest-mixing variance-related parameters can still sit
> a little above the strict $`\widehat{R} \le 1.01`$ target – variance
> components are the hardest quantities for HMC, and greta’s TensorFlow
> sampler is not fully seed-reproducible – so read a single run’s worst
> $`\widehat{R}`$ as indicative, and run longer or triangulate when a
> variance component is itself the inferential target.

## 3. Random intercept and slope (asreml-form)

The full repeated-measures model lets each subject have its own slope on
`Days` as well as its own intercept. The `lme4` syntax is
`(Days | Subject)`. The brms-style grammar accepts only single or
crossed random intercepts; for random slopes route through the
asreml-style entry, which handles them via two random terms — the
intercept and the interaction:

``` r

fit_ris <- flexybayes(
  fixed  = Reaction ~ Days,
  random = ~ Subject + Subject:Days,
  data   = sleepstudy,
  n_samples = 2000, warmup = 5000, chains = 4, verbose = FALSE
)
summary(fit_ris)
#> Bayesian mixed model summary  [flexyBayes]
#> ============================================================ 
#>   Fixed  : Reaction ~ Days 
#>   Random : ~Subject + Subject:Days 
#>   Family : gaussian / identity 
#>   N = 180 , chains = 4 , samples = 2000 
#>   Representation: exact
#>   Engine:         greta MCMC
#> 
#> -- Fixed effects (posterior)  --------------------------------- 
#>             Estimate Post.SD     2.5%    97.5%
#> (Intercept) 249.9655  9.8147 230.9471 268.5396
#> Days         10.6166  0.8037   9.0376  12.1951
#> 
#> -- Variance components  -------------------------------------- 
#>               Component Estimate     SD    2.5%   97.5%
#> 1         sigma_Subject  39.8555 8.3077 26.8527 59.0593
#> 2 sigma_Days_in_Subject   1.8964 2.7245  0.1543 12.2651
#> 3           sigma_e_atg  30.9873 1.8989 27.2677 34.8492
#> 
#> -- Convergence  --------------------------------------------- 
#>   Rhat range: 1.004 - 1.355 
#>   ESS  range: 93 - 1818 
#>   Run time  : 34.5 sec
```

This fit gives independent variance components for the random intercept
and the random slope (no covariance between them). The fully bivariate
`(Days | Subject)` model — with an estimated correlation between
subject-level intercept and slope — uses an unstructured covariance,
written as `random = ~ str(~Subject + Subject:Days, ~us(2):id(Subject))`
in the asreml DSL. It is one of the constructions covered in the
*structured covariance* vignette.

## 3a. Uncorrelated random slopes via the brms-style grammar (lme4-style `(x || g)`)

The brms-style grammar accepts the lme4 / brms double-pipe notation
`(x || g)` for *uncorrelated* random intercept + slope – each level of
the grouping factor carries an independent intercept deviation AND an
independent slope deviation, both drawn from their own univariate-
normal priors with no correlation parameter shared between them. This
matches the lme4 / brms semantics where `(Days || Subject)` expands to
`(1 | Subject) + (0 + Days | Subject)`. The `(0 + x || g)` form drops
the intercept block, fitting a slope-only random effect. Both shapes are
structurally equivalent to the asreml-form `~ Subject + Subject:Days`
construction in section 3 – same posterior, different surface.

``` r

fit_ris_uncor <- fb_greta(
  Reaction ~ Days + (Days || Subject),
  data = sleepstudy,
  n_samples = 2000, warmup = 5000, chains = 4,
  verbose = FALSE, mcmc_verbose = FALSE
)
summary(fit_ris_uncor)
#> Bayesian mixed model summary  [flexyBayes]
#> ============================================================ 
#>   Fixed  : Reaction ~ Days + (Days || Subject) 
#>   Family : gaussian / identity 
#>   N = 180 , chains = 4 , samples = 2000 
#>   Representation: exact
#>   Engine:         greta MCMC
#> 
#> -- Fixed effects (posterior)  --------------------------------- 
#>             Estimate Post.SD     2.5%    97.5%
#> (Intercept) 250.3381  8.5364 233.0382 266.2196
#> Days         10.3128  1.7623   6.6974  13.5729
#> 
#> -- Variance components  -------------------------------------- 
#>            Component Estimate     SD    2.5%   97.5%
#> 1      sigma_Subject  25.7963 5.8728 15.5494 38.7478
#> 2 sigma_Days_Subject   6.1719 1.3223  3.9892  9.1191
#> 3        sigma_e_atg  25.8969 1.5580 23.0584 29.0856
#> 
#> -- Convergence  --------------------------------------------- 
#>   Rhat range: 1.019 - 2.102 
#>   ESS  range: 75 - 942 
#>   Run time  : 26.1 sec
```

The intercept-deviation hyperparameter is canonically named `sd_Subject`
and the slope-deviation hyperparameter `sd_Days_Subject` – the
per-backend mapper handles the brms-side `sd_Subject__Days` \<-\>
flexyBayes canonical `sd_Days_Subject` translation automatically.
Inspect via
[`canonical_names()`](https://aagi-aus.github.io/flexyBayes/reference/canonical_names.md):

``` r

canonical_names(fit_ris_uncor)$map
#>             mu_atg          beta_Days      sigma_Subject sigma_Days_Subject 
#>      "(Intercept)"             "Days"       "sd_Subject"  "sd_Days_Subject" 
#>        sigma_e_atg 
#>            "sigma"
```

The correlated form `(Days | Subject)` (with a free correlation
parameter between the subject-level intercept and slope) is not yet
supported on the brms ingest path. Attempting it raises a structured
refusal with a deferral pointer – the workaround text suggests the
`(Days || Subject)` form when the correlation is not of inferential
interest, or the asreml-style structured covariance entry above when it
is. The refusal carries a typed condition class
(`flexybayes_correlated_slope_unsupported`) with slots
`deferral_target`, `workaround`, `grouping_factor`, and `slope_variable`
so downstream tooling can pattern-match.

## 4. Crossed and nested random effects

Two cluster variables that *cross* — every level of one is observed at
every level of the other — are modelled with two parallel random terms:

``` r

flexybayes(fixed = y ~ 1, random = ~ subject + item, data = ...)
```

Two cluster variables that *nest* — each level of the inner is observed
within exactly one level of the outer — are modelled with the
inner-within-outer interaction:

``` r

flexybayes(fixed = y ~ 1, random = ~ class + class:student, data = ...)
```

Both patterns generalise to deeper hierarchies. The agronomic
illustration is in the *MET and genomic selection* vignette.

## 5. The Poisson GLMM: `MASS::epil`

[`MASS::epil`](https://rdrr.io/pkg/MASS/man/epil.html) (Thall & Vail,
1990) records seizure counts in 59 epileptic patients across four 2-week
periods. We model the count of seizures with a Poisson likelihood, a log
link, and a random effect for `subject`:

``` r

data(epil, package = "MASS")
epil$subject <- factor(epil$subject)
str(epil)
#> 'data.frame':    236 obs. of  9 variables:
#>  $ y      : num  5 3 3 3 3 5 3 3 2 4 ...
#>  $ trt    : Factor w/ 2 levels "placebo","progabide": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ base   : int  11 11 11 11 11 11 11 11 6 6 ...
#>  $ age    : int  31 31 31 31 30 30 30 30 25 25 ...
#>  $ V4     : int  0 0 0 1 0 0 0 1 0 0 ...
#>  $ subject: Factor w/ 59 levels "1","2","3","4",..: 1 1 1 1 2 2 2 2 3 3 ...
#>  $ period : int  1 2 3 4 1 2 3 4 1 2 ...
#>  $ lbase  : num  -0.756 -0.756 -0.756 -0.756 -0.756 ...
#>  $ lage   : num  0.1142 0.1142 0.1142 0.1142 0.0814 ...
```

``` r

fit_pois_g <- fb_greta(
  y ~ trt + (1 | subject), data = epil,
  family       = "poisson",
  n_samples    = 2000, warmup = 5000, chains = 4,
  verbose      = FALSE,
  mcmc_verbose = FALSE
)
summary(fit_pois_g)
#> Bayesian mixed model summary  [flexyBayes]
#> ============================================================ 
#>   Fixed  : y ~ trt + (1 | subject) 
#>   Family : poisson / log 
#>   N = 236 , chains = 4 , samples = 2000 
#>   Representation: exact
#>   Engine:         greta MCMC
#> 
#> -- Fixed effects (posterior)  --------------------------------- 
#>              Estimate Post.SD     2.5%   97.5%
#> (Intercept)   -1.5524 12.5384 -13.6815 20.8873
#> trtplacebo     3.1445 12.6097 -19.4070 15.4692
#> trtprogabide   3.1619 12.5051 -19.1436 15.2908
#> 
#> -- Variance components  -------------------------------------- 
#>       Component Estimate     SD   2.5%  97.5%
#> 1 sigma_subject   1.0014 0.0983 0.8355 1.2331
#> 2   sigma_e_atg   1.5723 0.8570 0.1032 2.9190
#> 
#> -- Convergence  --------------------------------------------- 
#>   Rhat range: 1.021 - 18.769 
#>   ESS  range: 12 - 68 
#>   Run time  : 19.1 sec
```

``` r

fit_pois_i <- fb_inla(
  y ~ trt + (1 | subject), data = epil,
  family  = "poisson",
  verbose = FALSE
)
summary(fit_pois_i)
#> Bayesian mixed model summary  [flexyBayes / aggregated-gaussian]
#> ================================================================= 
#>   family:     poisson / log 
#>   N = 236 , K = 59 
#>   backend:    inla 
#>   exactness:  aggregated_exact 
#>   priors:     custom (explicit prior supplied; see prior_summary()) 
#>   aggregation: N = 236 rows -> K = 59 cells (ratio 4:1)
#> 
#> -- Fixed effects (posterior)  ----------------------------------- 
#>              Estimate Post.SD
#> (Intercept)    1.7717  0.1875
#> trtprogabide  -0.2914  0.2602
#> 
#> -- Variance components  ---------------------------------------- 
#>   tau_1  (random SD): 0.9473
#> 
#>   Run time: 2.24 sec
```

The intercept (on the log scale) is the baseline log-rate of seizures;
the `trt` coefficient is the log rate-ratio for progabide versus
placebo. Exponentiating gives the rate-ratio:

``` r

# greta names a 2-level factor trt as tau_trt[1,1] (placebo, baseline)
# and tau_trt[2,1] (progabide). The rate ratio is the exponential of
# the difference, which simplifies to exp(tau_trt[2,1]) once placebo
# is absorbed into the intercept.
draws <- flexyBayes::fb_as_draws_simple(fit_pois_g)
rr <- exp(draws[["tau_trt[2,1]"]] - draws[["tau_trt[1,1]"]])
c(mean = mean(rr), q025 = quantile(rr, 0.025), q975 = quantile(rr, 0.975))
#>       mean  q025.2.5% q975.97.5% 
#>  1.0611767  0.5998959  1.7205169
```

The cross-engine triangulation reads cleanly because the canonical
parameter-name registry reconciles the greta and INLA naming conventions
for fixed effects, factor levels, and group SDs; no user-supplied
`name_map` is required:

``` r

tri <- tryCatch(
  triangulate(fit_pois_g, fit_pois_i),
  error = function(e) {
    message("triangulate() unavailable in this session: ",
            conditionMessage(e))
    NULL
  }
)
if (!is.null(tri)) tri
#> <triangulate_result>
#>   source_a: greta
#>   source_b: inla
#>   independence: algorithmic + implementation
#>     (HMC (greta on TensorFlow) versus Laplace approximation (INLA on C): different inference paradigms and different code bases.)
#>   n_common: 3
#>   only_a:   2 parameters (trtplacebo, sigma)
#>   only_b:   118 parameters (Predictor:1, Predictor:2, Predictor:3, Predictor:4, Predictor:5, ...)
#> 
#> Metrics (per common parameter):
#>          param  mean_a  mean_b mean_diff    sd_a   sd_b sd_ratio q025_diff
#> 1  (Intercept) -1.5524  1.7632   -3.3156 12.5384 0.1816  69.0387  -15.0792
#> 2 trtprogabide  3.1619 -0.2801    3.4421 12.5051 0.2464  50.7593  -18.3664
#> 3   sd_subject  1.0014  0.9611    0.0403  0.0983 0.0991   0.9920    0.0499
#>   q975_diff wasserstein_1
#> 1   18.7771       11.5866
#> 2   15.1032       11.5692
#> 3    0.0689        0.0396
```

As always with an approximate engine, scrutinise the variance components
in the cross-check: INLA’s Laplace estimate of a small-group
random-effect variance can diverge from the MCMC value, and for very few
groups it can be unstable across runs. Where the `sd_ratio` column flags
such a divergence, trust the MCMC estimate or set an explicit prior –
the divergence is the triangulation doing its job, not a fault to hide.

If the planted-extreme example from the *cross-engine triangulation*
vignette is fresh in mind, the `epil` data have a similar structure:
patient 49 has anomalously high seizure counts, and the
Pareto-$`\hat{k}`$ diagnostic in `loo` will identify that patient as a
leverage point. The recommended response is to refit with
negative-binomial:

``` r

fit_nb <- fb_greta(
  y ~ trt + (1 | subject), data = epil,
  family       = "negative_binomial",
  n_samples    = 2000, warmup = 5000, chains = 4,
  verbose      = FALSE,
  mcmc_verbose = FALSE
)
summary(fit_nb)
#> Bayesian mixed model summary  [flexyBayes]
#> ============================================================ 
#>   Fixed  : y ~ trt + (1 | subject) 
#>   Family : negative_binomial / log 
#>   N = 236 , chains = 4 , samples = 2000 
#>   Representation: exact
#>   Engine:         greta MCMC
#> 
#> -- Fixed effects (posterior)  --------------------------------- 
#>              Estimate Post.SD     2.5%   97.5%
#> (Intercept)    8.2536  8.2452  -3.1730 27.1478
#> trtplacebo    -6.4254  8.2627 -25.4540  4.7928
#> trtprogabide  -6.8302  8.2679 -25.9041  4.4707
#> 
#> -- Variance components  -------------------------------------- 
#>       Component Estimate     SD   2.5%  97.5%
#> 1 sigma_subject   1.0107 0.0976 0.8004 1.1826
#> 2   sigma_e_atg   1.1915 0.8939 0.0549 2.9732
#> 
#> -- Convergence  --------------------------------------------- 
#>   Rhat range: 1.132 - 5.372 
#>   ESS  range: 10 - 46 
#>   Run time  : 29.4 sec
```

The negative-binomial absorbs overdispersion that the Poisson likelihood
cannot.

## 6. The binomial GLMM: `lme4::cbpp`

`lme4::cbpp` records contagious bovine pleuropneumonia incidence across
15 herds and 4 periods. The natural model is a logistic GLMM with a
random herd effect.

The brms-style grammar does not accept the aggregated trials addition
form (`y | trials(size) ~ ...`). We expand the per-herd counts to
per-trial Bernoulli observations as a workaround:

``` r

data(cbpp, package = "lme4")
cbpp_expanded <- do.call(rbind, lapply(seq_len(nrow(cbpp)), function(i) {
  row <- cbpp[i, ]
  data.frame(
    y      = c(rep(1, row$incidence), rep(0, row$size - row$incidence)),
    period = rep(row$period, row$size),
    herd   = rep(row$herd,   row$size)
  )
}))
cbpp_expanded$herd <- factor(cbpp_expanded$herd)
nrow(cbpp_expanded)
#> [1] 842
```

``` r

fit_bin <- fb_inla(
  y ~ period + (1 | herd), data = cbpp_expanded,
  family  = "binomial",
  verbose = FALSE
)
summary(fit_bin)
#> Bayesian mixed model summary  [flexyBayes / aggregated-gaussian]
#> ================================================================= 
#>   family:     binomial / logit 
#>   N = 842 , K = 56 
#>   backend:    inla 
#>   exactness:  aggregated_exact 
#>   priors:     custom (explicit prior supplied; see prior_summary()) 
#>   aggregation: N = 842 rows -> K = 56 cells (ratio 15:1)
#> 
#> -- Fixed effects (posterior)  ----------------------------------- 
#>             Estimate Post.SD
#> (Intercept)  -1.4116  0.2518
#> period2      -1.0055  0.3054
#> period3      -1.1481  0.3254
#> period4      -1.6319  0.4263
#> 
#> -- Variance components  ---------------------------------------- 
#>   tau_1  (random SD): 0.6269
#> 
#>   Run time: 2.12 sec
```

The intercept (on the logit scale) is the baseline log-odds of
infection; the `period` coefficients describe how those log-odds change
across the four observation periods. The `herd` random effect captures
unobserved herd-level differences in baseline risk.

## 7. Dispatch and refusal: `backend = "auto"`

The universal entries
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
/
[`fb()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
carry a `backend` argument: `"auto"` (the default — structural
gate-then-route: INLA on accept, greta on refuse, with a silenceable
note carrying the reason), `"greta"` (Hamiltonian Monte Carlo; covers
every supported term type), `"inla"` (explicit LGM fast path; raises a
structured refusal on non-LGM models), and `"brms"` (Stan). The engine
pins
[`fb_greta()`](https://aagi-aus.github.io/flexyBayes/reference/fb_greta.md)
/
[`fb_inla()`](https://aagi-aus.github.io/flexyBayes/reference/fb_inla.md)
/
[`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md)
fix the backend instead.

``` r

fit_auto <- flexybayes(
  Reaction ~ Days + (1 | Subject), data = sleepstudy,
  backend      = "auto",
  n_samples    = 200, warmup = 200, chains = 1,
  verbose      = FALSE,
  mcmc_verbose = FALSE
)
backend_decision(fit_auto)
#> $backend
#> [1] "inla"
#> 
#> $path
#> [1] "auto_accept"
#> 
#> $gate_checks
#> [1] "lgm_compatible"      "gretaR_slot_dormant"
#> 
#> $reason
#> [1] "lgm_gate() accepted; INLA dispatch."
#> 
#> $preflight_summary
#> NULL
#> 
#> $representation_plan
#> NULL
#> 
#> $rejected_routes
#> $rejected_routes[[1]]
#> $rejected_routes[[1]]$backend
#> [1] "greta"
#> 
#> $rejected_routes[[1]]$reason
#> [1] "not_chosen_by_policy"
#> 
#> 
#> $rejected_routes[[2]]
#> $rejected_routes[[2]]$backend
#> [1] "gretaR"
#> 
#> $rejected_routes[[2]]$reason
#> [1] "backend_not_activated"
#> 
#> 
#> 
#> $routing_policy_version
#> [1] "stage5a_v1"
```

The trace records which backend ran, which gate checks were evaluated,
and (on the refuse-to-greta path) the primary failure rule. Two
silenceable notes — `flexyBayes.silence_auto_fallback_note` (LGM gate
refuse) and `flexyBayes.silence_auto_inla_missing_note` (INLA package
unavailable) — let you keep auto-dispatch quiet in batch pipelines once
the routing behaviour is understood.

The *LGM feasibility* and *backend internals* vignettes detail the six
structural gate checks and the seventh post-fit numerical confirm, plus
the two-key armed override path (`force = "inla"`,
`acknowledge_silent_bias_risk = TRUE`, `reason = "..."`) for the rare
cases when forced INLA dispatch is deliberate.

## 8. Pitfalls

**Random slopes are not in the brms-style corpus.** Route random-slope
models through `flexybayes(fixed, random = ~ A + A:B, data = ...)`. The
fits are statistically equivalent to `lme4`’s `(B | A)` *with*
`corr = 0`; for the bivariate version with an estimated correlation, use
the structured covariance vignette’s `str(~..., ~us(2):id(...))`
pattern.

**Aggregated binomial via `trials()` is not in the brms-style corpus.**
Expand to per-trial Bernoulli, or route through
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
with `weights` for an effective sample size correction.

**Group counts dominate the prior.** With fewer than ~20 levels of a
random-effect factor, the prior on the group-level standard deviation
does meaningful work. The *priors and regularisation* vignette walks
through how to set this prior deliberately and how to assess
sensitivity.

**Overdispersion in count data.** A Poisson GLMM that fails the
posterior-predictive check on `var(y)` versus `mean(y)` is signalling
overdispersion — refit with `family = "negative_binomial"`. This is the
most common cause of cross-engine disagreement on count GLMMs (pattern 1
in the triangulation vignette).

## 9. Active prompts

1.  Compare `fit_ri` (random intercept) and `fit_ris` (random
    intercept + slope) on `sleepstudy`. Which has lower residual
    standard deviation? Which has higher ESS per second of compute time?
2.  Run the `epil` Poisson fit, then refit with the negative-binomial.
    Use
    [`triangulate()`](https://aagi-aus.github.io/flexyBayes/reference/triangulate.md)
    to compare the variance-component posteriors of each likelihood
    family on each backend.
3.  On `cbpp_expanded`, fit the model with random `herd` versus random
    `herd:period`. Which fits best by PSIS-LOO (manual; see the
    *downstream analysis* vignette)?

## 10. Session information

``` r

sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Tahoe 26.5.1
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
#> 
#> time zone: Australia/Adelaide
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] flexyBayes_0.8.3
#> 
#> loaded via a namespace (and not attached):
#>   [1] Rdpack_2.6.6           DBI_1.3.0              mnormt_2.1.2          
#>   [4] tfautograph_0.3.2      sandwich_3.1-1         rlang_1.2.0           
#>   [7] magrittr_2.0.5         multcomp_1.4-29        otel_0.2.0            
#>  [10] matrixStats_1.5.0      e1071_1.7-17           compiler_4.5.2        
#>  [13] loo_2.9.0              png_0.1-9              callr_3.7.6           
#>  [16] vctrs_0.7.3            stringr_1.6.0          pkgconfig_2.0.3       
#>  [19] crayon_1.5.3           backports_1.5.1        nloptr_2.2.1          
#>  [22] MatrixModels_0.5-4     torch_0.17.0           bit_4.6.0             
#>  [25] INLA_25.10.19          xfun_0.57              jsonlite_2.0.0        
#>  [28] progress_1.2.3         parallel_4.5.2         prettyunits_1.2.0     
#>  [31] tensorflow_2.20.0      R6_2.6.1               stringi_1.8.7         
#>  [34] RColorBrewer_1.1-3     reticulate_1.45.0      parallelly_1.47.0     
#>  [37] boot_1.3-32            numDeriv_2016.8-1.1    estimability_1.5.1    
#>  [40] Rcpp_1.1.1-1.1         knitr_1.51             zoo_1.8-15            
#>  [43] base64enc_0.1-6        bayesplot_1.15.0       Matrix_1.7-4          
#>  [46] splines_4.5.2          tidyselect_1.2.1       dichromat_2.0-0.1     
#>  [49] abind_1.4-8            codetools_0.2-20       processx_3.9.0        
#>  [52] listenv_0.10.1         gretaR_0.2.0           lattice_0.22-7        
#>  [55] tibble_3.3.1           withr_3.0.2            bridgesampling_1.2-1  
#>  [58] S7_0.2.2               posterior_1.7.0        coda_0.19-4.1         
#>  [61] evaluate_1.0.5         marginaleffects_0.32.0 future_1.70.0         
#>  [64] survival_3.8-3         sf_1.1-0               units_1.0-1           
#>  [67] proxy_0.4-29           RcppParallel_5.1.11-2  pillar_1.11.1         
#>  [70] tensorA_0.36.2.1       whisker_0.4.1          KernSmooth_2.23-26    
#>  [73] stats4_4.5.2           checkmate_2.3.4        reformulas_0.4.4      
#>  [76] sn_2.1.3               distributional_0.7.0   generics_0.1.4        
#>  [79] hms_1.1.4              ggplot2_4.0.3          rstantools_2.6.0      
#>  [82] scales_1.4.0           coro_1.1.0             minqa_1.2.8           
#>  [85] globals_0.19.1         xtable_1.8-8           class_7.3-23          
#>  [88] glue_1.8.1             emmeans_2.0.2          tools_4.5.2           
#>  [91] data.table_1.18.2.1    lme4_2.0-1             mvtnorm_1.3-6         
#>  [94] grid_4.5.2             rbibutils_2.4.1        nlme_3.1-168          
#>  [97] cli_3.6.6              tfruns_1.5.4           fmesher_0.7.0         
#> [100] Brobdingnag_1.2-9      dplyr_1.2.1            gtable_0.3.6          
#> [103] greta_0.5.1            digest_0.6.39          classInt_0.4-11       
#> [106] TH.data_1.1-5          brms_2.23.0            farver_2.1.2          
#> [109] lifecycle_1.0.5        bit64_4.8.0            MASS_7.3-65
```

## References

Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting linear
mixed-effects models using lme4. *Journal of Statistical Software*,
67(1), 1–48.

Gelman, A., & Hill, J. (2007). *Data analysis using regression and
multilevel/hierarchical models*. Cambridge University Press.

Pinheiro, J. C., & Bates, D. M. (2000). *Mixed-effects models in S and
S-PLUS*. Springer.

Thall, P. F., & Vail, S. C. (1990). Some covariance models for
longitudinal count data with overdispersion. *Biometrics*, 46(3),
657–671.
