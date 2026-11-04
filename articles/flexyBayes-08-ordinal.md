# Ordinal scores in a field trial

Growth stage, disease severity and lodging are scored on ordered scales.
A score of 3 is above a score of 2, but the step from 2 to 3 need not
equal the step from 3 to 4, so a Gaussian model of the score is the
wrong model. This vignette fits two ordinal models to a simulated
split-plot trial with a spatial field, the design of Cao et al. (2026).

## The two families

Both families model a score $`Y`$ with categories $`1, \dots, K`$
through a linear predictor $`\eta`$ and $`K - 1`$ thresholds
$`\tau_1 < \dots < \tau_{K-1}`$, with $`F`$ the logistic distribution
function.

- `family = "cumulative"`, the proportional-odds model:
  $`P(Y \le k) = F(\tau_k - \eta)`$.
- `family = "cratio"`, the continuation-ratio or sequential model:
  $`P(Y > k \mid Y \ge k) = F(\eta - \tau_k)`$. A plant reaches stage
  $`k + 1`$ only after stage $`k`$, which suits growth stages.

On both families a larger $`\eta`$ means a higher score. The thresholds
follow brms’s convention on both engines. Cao et al. write
$`\alpha_k = -\tau_k`$.

## The trial

Five rows by 18 columns of plots, in three replicates of six columns.
Sowing depth is the main-plot factor, applied to strips two columns
wide, so a replicate holds the three depths once. Ten plants are scored
per plot on six growth stages, `Z1` to `Z6`.

``` r

set.seed(20261103)

plots <- expand.grid(row = 1:5, col = 1:18)
plots$rep <- (plots$col - 1L) %/% 6L + 1L
strip <- (plots$col - 1L) %/% 2L + 1L
depth_of_strip <- unlist(lapply(1:3, function(r) sample(c(30L, 80L, 120L))))
plots$depth <- factor(depth_of_strip[strip], levels = c(30L, 80L, 120L))
plots$plot <- seq_len(nrow(plots))

# A field correlated along rows (0.3) and columns (0.6).
ar1_chol <- function(n, rho) t(chol(rho^abs(outer(1:n, 1:n, "-"))))
field <- 0.7 * ar1_chol(5L, 0.3) %*% matrix(rnorm(90L), 5L, 18L) %*%
  t(ar1_chol(18L, 0.6))

eta <- 0.8 * (plots$depth == "80") + 1.6 * (plots$depth == "120") +
  rnorm(3L, 0, 0.3)[plots$rep] +
  rnorm(9L, 0, 0.35)[3L * (plots$rep - 1L) + as.integer(plots$depth)] +
  rnorm(90L, 0, 0.5) +
  field[cbind(plots$row, plots$col)]

tau <- c(-2.0, -0.8, 0.2, 1.2, 2.3)
trial <- plots[rep(1:90, each = 10L), ]
p_le <- plogis(outer(eta[trial$plot], tau, function(e, t) t - e))
stage <- 1L + rowSums(runif(nrow(trial)) > p_le)
trial$stage <- factor(paste0("Z", stage), levels = paste0("Z", 1:6),
  ordered = TRUE)
for (v in c("row", "col", "rep", "plot")) trial[[v]] <- factor(trial[[v]])

table(trial$depth, trial$stage)
#>      
#>        Z1  Z2  Z3  Z4  Z5  Z6
#>   30   30  47  56  60  52  55
#>   80   18  25  57  63  63  74
#>   120  11  15  31  62  73 108
```

## A cumulative model with the field

The random part follows the design: replicate, main plot (`rep:depth`),
the field over rows and columns, and a plot effect for the variation the
field does not carry.

``` r

fit_cum <- flexybayes(
  stage ~ depth,
  random = ~ rep + rep:depth + ar1(row):ar1(col) + plot,
  data = trial,
  family = "cumulative",
  backend = "inla",
  verbose = FALSE
)

rnd(tidy(fit_cum))
#>               term estimate std.error conf.low conf.high
#> 1          depth80    0.646     0.513   -0.368     1.715
#> 2         depth120    1.288     0.512    0.272     2.353
#> 3 threshold[Z1|Z2]   -2.335     0.383   -3.081    -1.578
#> 4 threshold[Z2|Z3]   -1.215     0.373   -1.942    -0.483
#> 5 threshold[Z3|Z4]   -0.182     0.371   -0.905     0.558
#> 6 threshold[Z4|Z5]    0.827     0.371    0.093     1.571
#> 7 threshold[Z5|Z6]    1.901     0.375    1.170     2.639
```

The depth effects were simulated at 0.8 and 1.6 and the thresholds at
-2.0, -0.8, 0.2, 1.2 and 2.3.

``` r

rnd(summary(fit_cum)$varcomp)
#>      component estimate std.error conf.low conf.high                prior
#> 1       sd_rep    0.567     0.676    0.061     2.381 half_normal(scale=5)
#> 2 sd_rep:depth    0.567     0.256    0.222     1.204 half_normal(scale=5)
#> 3   sd_spatial    0.615     0.153    0.372     0.970 half_normal(scale=5)
#> 4      rho_row    0.091     0.425   -0.737     0.799       engine default
#> 5      rho_col    0.548     0.199    0.076     0.843       engine default
#> 6      sd_plot    0.348     0.206    0.089     0.868 half_normal(scale=5)
#>                note
#> 1 weakly identified
#> 2                  
#> 3                  
#> 4                  
#> 5                  
#> 6
```

## Category probabilities

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
probability of each category for each row of the data.
[`predict()`](https://rdrr.io/r/stats/predict.html) does the same for
new rows. With `marginal = TRUE` it integrates over the random effects,
so the probabilities describe an unseen plot at that depth and not a
plot whose random effects are zero.

``` r

new <- data.frame(depth = factor(c(30, 80, 120), levels = c(30, 80, 120)))
round(predict(fit_cum, newdata = new, marginal = TRUE, seed = 1L), 3)
#>         Z1    Z2    Z3    Z4    Z5    Z6
#> [1,] 0.138 0.145 0.187 0.189 0.162 0.180
#> [2,] 0.089 0.107 0.158 0.187 0.189 0.270
#> [3,] 0.053 0.073 0.123 0.169 0.202 0.379
```

## A sequential model, and comparing the two

The continuation-ratio model takes the same call with
`family = "cratio"`.

``` r

fit_seq <- flexybayes(
  stage ~ depth,
  random = ~ rep + rep:depth + ar1(row):ar1(col) + plot,
  data = trial,
  family = "cratio",
  backend = "inla",
  verbose = FALSE
)

rnd(tidy(fit_seq))
#>               term estimate std.error conf.low conf.high
#> 1          depth80    0.450     0.453   -0.476     1.382
#> 2         depth120    0.923     0.450    0.000     1.849
#> 3 threshold[Z1|Z2]   -2.414     0.652   -3.737    -1.089
#> 4 threshold[Z2|Z3]   -1.850     0.648   -3.165    -0.530
#> 5 threshold[Z3|Z4]   -1.024     0.645   -2.332     0.294
#> 6 threshold[Z4|Z5]   -0.273     0.645   -1.578     1.048
#> 7 threshold[Z5|Z6]    0.500     0.648   -0.806     1.829
```

The two families put different meanings on $`\eta`$, so their
coefficients are not compared directly. Their predictions are.
[`fb_kfold()`](https://aagi-aus.github.io/flexyBayes/reference/fb_kfold.md)
refits each model with one fold of plots held out and scores the
held-out plots. The ten plants of a plot share its random effects, so
the plot is scored as a whole with `group = "plot"`.

``` r

folds <- rep(rep(1:5, length.out = 90L), each = 10L)
cv_cum <- fb_kfold(fit_cum, folds, group = "plot", seed = 1L)
cv_seq <- fb_kfold(fit_seq, folds, group = "plot", seed = 1L)
rbind(cumulative = cv_cum$estimates, cratio = cv_seq$estimates)
#>                 elpd       se
#> cumulative -1482.489 21.32669
#> cratio     -1489.586 21.81711
```

A larger `elpd` predicts the held-out plots better. The data were
simulated from the cumulative model.

## Engines and priors

INLA fits both families with the random structures it fits for a
binomial response, the field among them. brms fits both families without
the field.

The thresholds carry a prior on each engine. On INLA the cumulative
family’s thresholds are hyperparameters of the likelihood, with a
uniform Dirichlet prior on the category probabilities at $`\eta = 0`$.
The continuation-ratio family’s thresholds are fixed effects with INLA’s
fixed-effect prior. On brms the cumulative thresholds carry the same
Dirichlet prior, and the continuation-ratio thresholds the prior of
class `Intercept`. The printed fit names the threshold prior.

``` r

fit_cum
#> Bayesian mixed model  [flexyBayes / inla]
#> -------------------------------------------------------------- 
#>   Fixed    : stage ~ depth 
#>   Random   : ~rep + rep:depth + ar1(row):ar1(col) + plot 
#>   Model    : G: rep iid; depth within rep; ar1(row):ar1(col) field + nugget (4 parameters); plot iid; R: units 
#>              this is not ASReml residual = ~ ar1:ar1 (3 parameters)
#>   Family   : cumulative / logit 
#>   Representation: exact (model, not inference)
#>   Engine:         INLA Laplace
#>   N        : 900 design rows, 900 observed responses
#>   na_action: augment; 0 missing response(s); 0 cell(s) completed
#>   Runtime  : 9.3 sec
#>   Versions : R version 4.6.1 (2026-06-24); flexyBayes 0.19.0; INLA 26.8.7
#>   Threads  : INLA num.threads 1:1; OMP_NUM_THREADS 8
#>   Seed     : none supplied; INLA takes no seed
#>   INLA formula: fb_ordinal_y ~ 1 + depth + f(rep, model = "iid", <prior: half-normal-SD(5)>) + f(rep_depth_id, model = "iid", <prior: half-normal-SD(5)>) + f(row_id, model = "ar1", group = col_id, control.group = list(model = "ar1"), <prior: half-normal-SD(5)>) + f(plot, model = "iid", <prior: half-normal-SD(5)>)
#>   Params   : 3 fixed, 4 random, 11 hyperparameter(s)
#>   numerical confirm: PASS
#> 
#> Separable AR1 field (posterior median, 95% interval):
#>                        parameter median    2.5%  97.5%
#>  correlation along row (rho_row) 0.1189 -0.7366 0.7988
#>  correlation along col (rho_col) 0.5793  0.0760 0.8428
#>           field SD (sigma_field) 0.5962  0.3714 0.9720
#>   The field and the nugget are separate parameters: this is a latent
#>   autoregressive field plus independent observation noise, not ASReml's
#>   nugget-free separable residual. Variances are the squares of the SDs.
#>   Field against the same model without it: not refitted: the interval of the correlation along col (rho_col) excludes zero.
#> -------------------------------------------------------------- 
#>   Object of class <flexybayes_inla>
#>   $inla: raw INLA fit (use INLA's summary, plot, etc.)
#>   $fb:   the fb_terms IR used for dispatch
#> 
#> -- Ordinal response ------------------------------------------
#>   cumulative / logit, 6 categories: Z1 < Z2 < Z3 < Z4 < Z5 < Z6
#>   Thresholds: P(Y <= k) = F(threshold_k - eta), under a Dirichlet(1) prior on the category probabilities at a linear predictor of zero, fixed by flexyBayes
#>              term estimate std.error conf.low conf.high
#>  threshold[Z1|Z2]  -2.3345    0.3835  -3.0810   -1.5783
#>  threshold[Z2|Z3]  -1.2146    0.3727  -1.9416   -0.4827
#>  threshold[Z3|Z4]  -0.1816    0.3706  -0.9053    0.5578
#>  threshold[Z4|Z5]   0.8268    0.3709   0.0934    1.5706
#>  threshold[Z5|Z6]   1.9010    0.3755   1.1696    2.6395
```

## What is refused

Each model the package does not fit is refused by name, with the reason.

- On INLA, the logit link only. A probit or complementary log-log link
  goes to brms, which does not fit the field.
- On INLA, at most 15 hyperparameters. The cumulative family’s
  thresholds count among them, five here. The continuation-ratio family
  fits its thresholds as fixed effects, so a larger random structure can
  fit there.
- Adjacent-category and stopping-ratio families, an offset, and a model
  without an intercept.
- Methods that read a mean response or a residual, such as
  [`residuals()`](https://rdrr.io/r/stats/residuals.html).

``` r

residuals(fit_cum)
#> Error:
#> ! residuals() reads a mean response, a residual or a per-row likelihood, which an ordinal fit does not have. fitted() returns each row's category probabilities.
```

## Reference

Cao Z, Megirian R, Nguyen M, Sparks AH (2026). Bayesian ordinal analysis
of cereal growth stage data from oat trials in Western Australia.
arXiv:2610.37384.
