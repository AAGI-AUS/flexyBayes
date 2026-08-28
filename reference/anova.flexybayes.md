# Compare flexyBayes models on a plug-in information criterion

Ranks two or more fits by a DIC-shaped criterion built from the plug-in
conditional log-likelihood of
[`logLik.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/logLik.flexybayes.md)
and the recorded parameter count. The criterion is approximate and the
method says so on every print: it penalises the nominal parameter count
rather than an effective one, so it is a coarse ordering, not a
model-selection procedure.

## Usage

``` r
# S3 method for class 'flexybayes'
anova(object, ...)
```

## Arguments

- object:

  A `flexybayes` fit supplying a conditional log-likelihood and a
  recorded parameter count.

- ...:

  Further `flexybayes` fits to compare against `object`.

## Value

Invisibly, a data frame with one row per model carrying the
log-likelihood, parameter count, criterion value, and the difference
from the best model. Printed as a side effect.

## Details

Every fit compared must supply both ingredients. A fit whose
log-likelihood is `NA` – an INLA fit reports a marginal likelihood, not
a conditional one – or which records no parameter count is refused by
name, because ranking on an `NA` produces an ordering that looks
computed and is not.
