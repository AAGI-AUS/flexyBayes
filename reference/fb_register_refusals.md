# Register a package's refusal reasons

A package that extends flexyBayes calls `fb_register_refusals()` from
its own `.onLoad()` to add the reasons it refuses with. Its users then
catch those refusals with the same `flexybayes_refusal` condition class
as the core's, find them in
[`fb_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refusals.md),
and the package raises them with
[`fb_refuse()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refuse.md).

## Usage

``` r
fb_register_refusals(refusals, package = NULL)
```

## Arguments

- refusals:

  A data frame with one row per reason and these character columns:

  - `reason_code` (required): lower-case letters, digits and
    underscores, starting with a letter. Prefix codes with something
    particular to the package so they do not collide with others.

  - `description` (required): one sentence on what the refusal means,
    shown by
    [`fb_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refusals.md)
    and used by
    [`fb_refuse()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refuse.md)
    when no message is given.

  - `since_version` (required): the package version that introduced the
    reason.

  - `message_template` (optional, default `"%s"`): the message shape,
    for documentation.

  - `registered_in_adr` (optional, default `NA`): the design record the
    reason traces to.

  - `plan_field` (optional, default `NA`): the
    [`fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/fb_plan.md)
    slot that reports the refusal before a fit, if any.

- package:

  The package that owns the reasons. `NULL` (default) takes the package
  whose code makes the call, which is the right value inside
  `.onLoad()`. Outside a package it must be given.

## Value

The registered reasons as an `fb_refusals_table`, invisibly.

## Details

Each package's reasons go into a partition named after the package. A
reason code is unique across all partitions: a code that flexyBayes or
another package already registered is refused, and so is any attempt to
register under the name `"flexyBayes"`. Registering the same code again
from the same package with identical fields does nothing, so a package
can be loaded twice in one session. Registering it with different fields
is refused. A table that fails any check registers none of its rows.

## See also

[`fb_refuse()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refuse.md)
to raise a registered reason,
[`fb_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refusals.md)
to list them.

## Examples

``` r
if (FALSE) { # \dontrun{
# In the plug-in's R/zzz.R:
.onLoad <- function(libname, pkgname) {
  flexyBayes::fb_register_refusals(
    data.frame(
      reason_code = "mypkg_threshold_negative",
      description = "The threshold must be zero or positive.",
      since_version = "0.1.0"
    )
  )
}
} # }
```
