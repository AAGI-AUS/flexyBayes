# Re-fit a flexyBayes model with modified arguments

Rebuilds the original
[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md)
call from the argument record the fit carries, applies the overrides
supplied in `...`, and re-fits.

## Usage

``` r
# S3 method for class 'flexybayes'
update(object, ...)
```

## Arguments

- object:

  A `flexybayes` fit carrying a complete argument record in
  `$extras$call_info`.

- ...:

  Named arguments overriding the recorded call, for example
  `n_samples = 2000L` or a replacement `data` frame.

## Value

A new fitted `flexybayes` object of the same engine as the original,
unless an override routes it elsewhere.

## Details

The record has to be complete for this to be safe. An argument the fit
did not record would be re-supplied as its default, so a re-fit could
quietly drop a relationship matrix, a prior scale, or a missing-response
policy, and return a different model under the same name. When any
recorded argument is absent, the method refuses and names what is
missing rather than re-fitting a model the user did not ask for. A fit
produced before its engine recorded the full set therefore refuses,
which is the intended behaviour: there is no pre-0.9.1 object to
special-case into a silent default.

## The prior a re-fit runs under

A re-fit uses the prior the original fit resolved to, taken from the
record on the fit itself rather than re-derived from the arguments. This
matters because the auto-default bounded uniform on the
standard-deviation scale fires only when the call supplies neither a
`prior` nor a `prior_vc_sd`, and a re-fit that re-issues the recorded
`prior_vc_sd` supplies one by construction. Before 0.9.1 the
auto-default therefore never fired on a re-fit: an identity
[`update()`](https://rdrr.io/r/stats/update.html) fell through to the
engine's own hyperprior and could report a variance component half the
size of the one it re-fitted, with nothing on the object to say so.

Precedence, in order:

- A `prior` in `...` is used as given. The recorded prior is not
  re-issued alongside it.

- A `prior_vc_sd` in `...` is used as given, and the recorded prior is
  *not* re-issued – the two forms describe the same variance-component
  scale, and passing both would let the recorded prior silently override
  the scalar the caller just wrote.

- Otherwise the fit's own resolved prior decides, and *how* it decides
  depends on where that prior came from – see below.

- A fit carrying the legacy scalar bridge, or no prior record at all,
  re-fits on the recorded scalars exactly as it did before.

**A policy re-fires; a bespoke prior carries.** The two things
`$extras$fb_terms$priors` can hold are not the same kind of thing, and a
re-fit treats them differently.

- **The auto-default** bounded uniform is the output of a policy – one
  bounded uniform per variance component the model has, with the bound
  read off the data – rather than a statement about particular terms. A
  re-fit passes neither `prior` nor either scalar, so the policy runs
  again over the *updated* model and data. On an identity
  [`update()`](https://rdrr.io/r/stats/update.html) the same data
  rebuild the same bound, so nothing moves. On
  `update(fit, random = ~ g + h)` the added term gets the same bounded
  uniform as its siblings.

- **A user-supplied
  [`fb_prior()`](https://aagi-aus.github.io/flexyBayes/reference/fb_prior.md)**
  is a statement about the terms it names, and is carried verbatim. A
  term the user did not name – including one an override has just added
  – keeps following the engine's own hyperprior, exactly as it did on
  the first fit. `prior_fixed_sd` and `prior_vc_sd` are re-issued
  alongside it, so a route that reads them for the terms the prior
  object does not name keeps the values the original ran under.

The asymmetry is deliberate. Re-applying a policy's old output to a new
model leaves a mixed state – the original terms priored, the added one
falling to the engine while its siblings keep the uniform – which is not
what either the policy or the user asked for. Re-deciding a *user's*
explicit prior, on the other hand, would put words in their mouth. Read
[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
on the re-fit to see what it ran under, and
[`summary.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes.md)'s
`varcomp` table to see it per component.

## The engine and representation a re-fit runs on

`backend` and `aggregate` are re-issued from the record alongside
everything else. Before 0.9.1 neither was recorded, so a re-fit took the
formal defaults: an identity
[`update()`](https://rdrr.io/r/stats/update.html) of a Stan fit came
back as an aggregated INLA fit – a different inference engine and a
different representation, under the same name and with nothing said.

What is recorded is the value the *call* asked for, not the engine the
call resolved to. A fit made under `backend = "auto"` records `"auto"`,
and a re-fit routes again from there. A fit that named `"brms"` records
`"brms"` and comes back on Stan. The distinction is deliberate: the
recorded value is the user's policy, and a re-fit whose model has
changed – `update(fit, random = ~ g + h)` – has to be free to route the
changed model rather than be pinned to the engine that suited the old
one. `aggregate` behaves the same way. A recorded `"auto"` lets the
aggregation gate re-decide, a recorded `FALSE` keeps the re-fit per-row.

An override in `...` wins, as for every other recorded argument, so
`update(fit, backend = "brms")` moves the re-fit to Stan.

What [`update()`](https://rdrr.io/r/stats/update.html) reproduces is the
model and the policy behind it, not the display settings of the session
that first ran it. `verbose` is recorded on the fit for completeness and
is deliberately *not* re-issued. A re-fit follows the current call's
display default, and `update(fit, verbose = FALSE)` quietens one on the
spot. The distinction is in what the argument can get wrong. A silently
switched engine is a different answer under the same name. A banner
printed where the first fit printed none is a banner.

## See also

[`prior_summary()`](https://aagi-aus.github.io/flexyBayes/reference/prior_summary.md)
for the resolved prior a fit carries, and
[`summary.flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/summary.flexybayes.md),
whose `varcomp` table names the prior each component ran under.
