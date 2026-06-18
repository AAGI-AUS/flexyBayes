# Resolved-prior summary for a flexyBayes fit

Returns the resolved priors used to fit the model – either the
user-supplied
[`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)
object, the auto-default bounded uniform on the standard-deviation scale
(Gelman, 2006), or the legacy scalar bridge (`prior_fixed_sd` +
`prior_vc_sd`). The return value is an S3 object with a
[`print()`](https://rdrr.io/r/base/print.html) method; the underlying
`fb_prior` (when applicable) is exposed under `$fb_prior` for
programmatic access.

## Usage

``` r
prior_summary(object, ...)

# S3 method for class 'flexybayes'
prior_summary(object, ...)

# S3 method for class 'flexybayes_inla'
prior_summary(object, ...)

# S3 method for class 'flexybayes_brms'
prior_summary(object, ...)

# S3 method for class 'flexybayes_direct_greta'
prior_summary(object, ...)

# Default S3 method
prior_summary(object, ...)
```

## Arguments

- object:

  A `flexybayes`, `flexybayes_inla`, or `flexybayes_direct_greta`
  object.

- ...:

  Ignored by current methods (reserved for future per-component
  selection).

## Value

A `prior_summary_flexybayes` object (list). Components:

- `kind`:

  One of `"fb_prior"`, `"legacy_scalar"`, `"no_prior_recorded"`.

- `backend`:

  The backend the fit ran on: `"greta"`, `"inla"`, or `"greta-direct"`.

- `fb_prior`:

  The `fb_prior` object (when `kind == "fb_prior"`).

- `default_origin`:

  `"auto"` when the prior was constructed by the bounded-uniform
  auto-default; `"user"` when supplied via the `prior` argument; `NA`
  for legacy / no-prior cases.

- `default_scale`, `default_basis`:

  Attributes carried by the auto-default prior naming the response-scale
  upper bound and its basis (response-scale `sd(y)`, logit-scale
  constant, log-scale constant). `NULL` when the prior was
  user-supplied.

- `fixed_sd`, `vc_sd`:

  Legacy scalar values (when `kind == "legacy_scalar"`).

- `declaration_only`:

  `TRUE` on `flexybayes_direct_greta` fits – the prior is a declaration
  of the user's greta-built model, not an enforcement.

## Details

For `flexybayes_direct_greta` fits the priors are a *declaration* of
what the user-built model graph encodes; the summary flags this with
`declaration_only = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
# live brms (Stan) fit -- needs a working Stan toolchain
data(sleepstudy, package = "lme4")
fit <- fb_brms(Reaction ~ Days + (1 | Subject),
               data = sleepstudy,
               n_samples = 100, warmup = 100, chains = 1,
               verbose = FALSE, mcmc_verbose = FALSE)
prior_summary(fit)
} # }
```
