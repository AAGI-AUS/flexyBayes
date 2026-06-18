# emmeans support: estimation basis (greta backend)

emmeans support: estimation basis (greta backend)

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

  Ignored.

## Value

A list with `X`, `bhat`, `nbasis`, `V`, `dffun`, `dfargs`.
