# Emit a flexyBayes posterior as a log-density producer

Turns a fitted flexyBayes object into a log-posterior producer: a
vectorised, domain-safe, unnormalised log-posterior callable that
`proxymix::from_fb_posterior()` compresses into a closed-form
Gaussian-mixture proxy. It is the single inference-result outflow from
flexyBayes; the contract is the *log-density*, not the draws, so the
returned object is addressed purely through its callable.

## Usage

``` r
fb_log_posterior(fit, ...)

# Default S3 method
fb_log_posterior(fit, ...)

# S3 method for class 'flexybayes_brms'
fb_log_posterior(fit, ...)

# S3 method for class 'flexybayes_inla'
fb_log_posterior(fit, ...)

# S3 method for class 'flexybayes'
fb_log_posterior(fit, ...)
```

## Arguments

- fit:

  A fitted flexyBayes object. The greta classes (`flexybayes` /
  `flexybayes_direct_greta`) produce a real callable; the brms and INLA
  classes abstain.

- ...:

  Reserved for future producer options; currently unused.

## Value

For a greta fit, a bare callable `function(theta_matrix)` with the
attributes described above, ready to pass to
`proxymix::from_fb_posterior()`. For a brms or INLA fit, the function
does not return: it raises a classed `fb_c4_unavailable` condition.

## Details

The returned value is a **bare callable**
`function(theta_matrix) -> numeric`. Its input is a numeric matrix whose
rows index independent parameter draws and whose columns index
parameters, in `attr(., "parameter_names")` order, on the natural
(constrained) scale. Its output is a length-`nrow(theta_matrix)` numeric
vector of `log p(theta | data) + const` (unnormalised). The callable is
vectorised, side-effect free, and domain-safe: a row outside the
parameters' support returns `-Inf` rather than raising an error (the
consumer probes it at construction). It carries, as attributes:

- `parameter_names`:

  Character vector naming the parameters; its length fixes the proxy's
  ambient dimension. Vector-valued targets are flattened in column-major
  order with index suffixes (e.g. `beta[1,1]`, `beta[2,1]`).

- `log_normalizer`:

  The additive correction that would normalise the density, i.e.
  `-log Z`. For a posterior the marginal likelihood is generally
  unknown, so this is `NA_real_` – honest, and the consumer reports a
  shifted (not absolute) divergence.

- `support_lower`, `support_upper`:

  Length-`n_dim` numeric support bounds taken from the model's parameter
  constraints (`NA` for an unbounded coordinate). A variance / scale
  parameter is bounded below by zero, for instance. Used only to centre
  and scale the consumer's default importance proposal.

- `draws`:

  An `n` by `n_dim` numeric matrix of the fit's posterior draws on the
  natural scale, column-aligned to `parameter_names`. Used only to seed
  the consumer's default proposal; never required.

Backend support. The **greta** backend is the canonical real producer:
it evaluates the model graph's unadjusted joint density at the
free-state image of the supplied natural-scale parameters, which is the
unnormalised natural-scale log-posterior exactly. The **brms** and
**INLA** backends abstain with an informative condition – brms's
log-density lives on the Stan unconstrained scale with a version-fragile
name mapping, and INLA's posterior is a deterministic approximation, not
a sampling log-density; an honest abstain is preferred to a
plausible-but-wrong log-density.

Acyclic note. A consumer such as proxymix uses this callable without
depending on flexyBayes; flexyBayes does not list proxymix in `Imports`
or `Suggests`. The cross-package demonstration lives in a separate
integration harness, not in this package, preserving the acyclic
dependency graph.

## See also

`proxymix::from_fb_posterior()` for the consumer (compresses the
returned callable into a Gaussian-mixture proxy).

## Examples

``` r
if (FALSE) { # \dontrun{
library(greta)
n <- 30
y <- rnorm(n, 1.5, 2)
mu <- normal(0, 5)
sigma <- normal(0, 5, truncation = c(0, Inf))
yd <- as_data(y)
distribution(yd) <- normal(mu, sigma)
m <- model(mu, sigma)
fit <- fb_greta(fb_from_greta(m), n_samples = 500, warmup = 500,
                chains = 2, verbose = FALSE, mcmc_verbose = FALSE)
producer <- fb_log_posterior(fit)
attr(producer, "parameter_names")
producer(matrix(c(1.5, 2.0), nrow = 1)) # natural-scale log-posterior
## Compress with proxymix (in a separate integration harness):
## proxymix::from_fb_posterior(producer, N = 2)
} # }
```
