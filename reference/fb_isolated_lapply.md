# Run jobs, each in its own R process

`fb_isolated_lapply()` is
[`lapply()`](https://rdrr.io/r/base/lapply.html) with each call of `fn`
made in a separate R process started through callr, so that a process
ended by the operating system or by a crash in compiled code is reported
as a refusal and the session goes on. It is the mechanism behind
`options(flexyBayes.isolation = "process")`, offered to packages that
extend flexyBayes for their own loops of model fits.

## Usage

``` r
fb_isolated_lapply(jobs, fn, plan = list())
```

## Arguments

- jobs:

  A list or vector, `fn` is called on each element.

- fn:

  A function of one element of `jobs`.

- plan:

  A list:

  - `workers`: how many processes run at once, a whole number of one or
    more (default 1).

  - `threads_per_worker`: the OpenMP threads each process may use, or
    `NA` (default) to leave them as the session has them. A process also
    runs one chain at a time when this is set.

  - `warm_first`: `TRUE` runs the first job alone before the others
    start, so that a Stan program compiled by it is cached for the rest
    (default `FALSE`).

  - `engine`: the engine the jobs fit with, named in a refusal (default
    `"none"`).

  - `entry_point`: the function the user called, named in a refusal
    (default `"fb_isolated_lapply()"`).

## Value

A list of the values `fn` returned, in the order of `jobs`.

## Details

Each process loads flexyBayes, and the package whose namespace `fn` is
defined in, from where this session loaded them, and starts from this
session's flexyBayes, brms, cmdstanr and INLA settings, its environment
variables and its random-number state. Messages, warnings and printed
output come back and are replayed in job order, never in completion
order. A job that fails stops new jobs from starting. The jobs before it
are waited for, and the lowest failing job's error is raised. With one
worker, each process starts from the random-number state the previous
one left, as [`lapply()`](https://rdrr.io/r/base/lapply.html) would,
with several, each starts from the session's state and the session's
state is left where it was.

A process that ends without returning raises the refusal
`isolated_process_died`, which names the engine, the entry point, the
job and the signal or exit status.

## See also

[`fb_refit()`](https://aagi-aus.github.io/flexyBayes/reference/fb_refit.md),
the refit that such jobs usually run.

## Examples

``` r
if (FALSE) { # interactive() && requireNamespace("callr", quietly = TRUE)
fb_isolated_lapply(1:2, function(i) i^2, list(workers = 2L))
}
```
