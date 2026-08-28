# emmeans support: estimation basis (default method)

The bare `flexybayes` method: reached by any fit whose engine has no
more specific override, which today is the brms backend.

## Usage

``` r
# S3 method for class 'flexybayes'
emm_basis(object, trms, xlev, grid, ...)

# S3 method for class 'flexybayes_inla'
emm_basis(object, trms, xlev, grid, ...)
```

## Arguments

- object:

  A `flexybayes_inla` fit.

- trms:

  Fixed-effect terms supplied by emmeans.

- xlev:

  Factor levels supplied by emmeans.

- grid:

  Reference grid supplied by emmeans.

- ...:

  Ignored. Present for compatibility with the generic.

## Value

A list with `X`, `bhat`, `nbasis`, `V`, `dffun`, `dfargs`.
