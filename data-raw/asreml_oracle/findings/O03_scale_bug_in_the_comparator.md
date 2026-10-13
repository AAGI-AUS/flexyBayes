# O03: the comparator read sectioned residuals on the wrong scale

**Found** 2026-09-11, in the oracle suite's own code, while reading scenario
05's output. Recorded because it changed a published-looking result and
because the way it was found is the same lesson the suite keeps producing.

## What it was

`fb_varcomp()` in `data-raw/asreml_oracle/compare.R` converts flexyBayes
components to the variance scale before comparing them against ASReml, which
reports variances. It tested for a standard deviation with

~~~r
v$component == "sigma" | startsWith(v$component, "sd_")
~~~

A sectioned residual has no scalar `sigma`. It reports `sigma_C1`,
`sigma_C2`, one per section, and those are standard deviations too. They
matched neither branch, so they were passed through unsquared and compared
against ASReml's variances.

## What it did

Two things, and the second is the one that matters.

**It faked a catastrophe.** Scenario 05's residual sections read 17.88
against ASReml's 299.2, which looks like a system returning answers an order
of magnitude wrong. Squared, 17.88^2 = 319.7 against 299.2, close
agreement, and a good result for the sectioned-residual translation.

**It corrupted a tolerance check.** A variance ratio is a component over the
total, and the total sums every component. In a sectioned model the six
residual sections dominate that sum, so leaving them unsquared moved the
denominator. Scenario 04's genotype variance ratio was reported as
`asr=0.0854 fb=0.0396`, a difference of -0.0458 against a tolerance of 0.05:
a pass, by a hair, and the kind of pass that invites a second look at the
translation. With the scale corrected the same comparison reads
`asr=0.0854 fb=0.0843`, a difference of **-0.0012**. The translation was
never marginal. The comparator was wrong.

## How it was found

By reading the component names. The first reading of scenario 05 compared the
two vectors positionally, concluded flexyBayes was returning residual
variances of 0.016 where ASReml returned 299, and started drafting that as a
serious defect. Printing the names showed flexyBayes orders sectioned
residuals first and genotype components second, and that the two vectors were
not aligned at all.

That is the third time in this suite's short life that a proxy was read in
place of the fact, after the prior column in O01, and after a hypothesis
about which hyperparameter rows a sectioned fit carries, which was also
wrong and also settled by printing them. A vector of numbers without its
names is a proxy.

## The guard

`fb_varcomp()` now tests `sigma`, the `sigma_` prefix and the `sd_` prefix,
and excludes correlations explicitly and not by assuming they cannot
match. The scenario results banked in `inst/validation/asreml-4.2/` are from
after the correction; the earlier ones were discarded and not kept for
comparison, because a fixture produced by a broken comparator is not
evidence of anything.
