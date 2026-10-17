# Emit the generated backend code for a deferred review object

Writes the backend code carried by a `<flexybayes_review>` object (Stan
code via
[`brms::stancode()`](https://paulbuerkner.com/brms/reference/stancode.html)
for
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
/
[`fb_brms()`](https://aagi-aus.github.io/flexyBayes/reference/fb_brms.md))
to a connection.

## Usage

``` r
cat_code(x, ...)

# S3 method for class 'flexybayes_review'
cat_code(x, file = stdout(), ...)
```

## Arguments

- x:

  A `<flexybayes_review>` object carrying the generated backend code.

- ...:

  Method-specific arguments. The `flexybayes_review` method accepts
  `file`, a connection defaulting to
  [`stdout()`](https://rdrr.io/r/base/showConnections.html).

- file:

  A connection to write the code to, defaulting to
  [`stdout()`](https://rdrr.io/r/base/showConnections.html).

## Value

Invisibly, the code as a character vector. Called for the code it writes
to `file`.

## Deprecated

`cat_code()` is deprecated in flexyBayes 0.13.0 and will be removed in a
later release. The generated code is the `code` element of the review
object, and `flexybayes(..., return_code = TRUE)` returns it directly. A
call warns with a condition of class `flexybayes_deprecated`.

## Examples

``` r
# Generating the code goes through brms::make_stancode(), so brms must be
# present.
if (requireNamespace("brms", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
  code <- flexybayes(y ~ x + (1 | g), data = d, backend = "brms",
                     return_code = TRUE)
  cat(code)
}
#> // generated with brms 2.23.0
#> functions {
#> }
#> data {
#>   int<lower=1> N;  // total number of observations
#>   vector[N] Y;  // response variable
#>   int<lower=1> K;  // number of population-level effects
#>   matrix[N, K] X;  // population-level design matrix
#>   int<lower=1> Kc;  // number of population-level effects after centering
#>   // data for group-level effects of ID 1
#>   int<lower=1> N_1;  // number of grouping levels
#>   int<lower=1> M_1;  // number of coefficients per level
#>   array[N] int<lower=1> J_1;  // grouping indicator per observation
#>   // group-level predictor values
#>   vector[N] Z_1_1;
#>   int prior_only;  // should the likelihood be ignored?
#> }
#> transformed data {
#>   matrix[N, Kc] Xc;  // centered version of X without an intercept
#>   vector[Kc] means_X;  // column means of X before centering
#>   for (i in 2:K) {
#>     means_X[i - 1] = mean(X[, i]);
#>     Xc[, i - 1] = X[, i] - means_X[i - 1];
#>   }
#> }
#> parameters {
#>   vector[Kc] b;  // regression coefficients
#>   real Intercept;  // temporary intercept for centered predictors
#>   real<lower=0> sigma;  // dispersion parameter
#>   vector<lower=0>[M_1] sd_1;  // group-level standard deviations
#>   array[M_1] vector[N_1] z_1;  // standardized group-level effects
#> }
#> transformed parameters {
#>   vector[N_1] r_1_1;  // actual group-level effects
#>   // prior contributions to the log posterior
#>   real lprior = 0;
#>   r_1_1 = (sd_1[1] * (z_1[1]));
#>   lprior += student_t_lpdf(Intercept | 3, 0.1, 2.5);
#>   lprior += normal_lpdf(sigma | 0, 4.27584420984108)
#>     - 1 * normal_lccdf(0 | 0, 4.27584420984108);
#>   lprior += normal_lpdf(sd_1 | 0, 4.27584420984108)
#>     - 1 * normal_lccdf(0 | 0, 4.27584420984108);
#> }
#> model {
#>   // likelihood including constants
#>   if (!prior_only) {
#>     // initialize linear predictor term
#>     vector[N] mu = rep_vector(0.0, N);
#>     mu += Intercept;
#>     for (n in 1:N) {
#>       // add more terms to the linear predictor
#>       mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
#>     }
#>     target += normal_id_glm_lpdf(Y | Xc, mu, b, sigma);
#>   }
#>   // priors including constants
#>   target += lprior;
#>   target += std_normal_lpdf(z_1[1]);
#> }
#> generated quantities {
#>   // actual population-level intercept
#>   real b_Intercept = Intercept - dot_product(means_X, b);
#> }
```
