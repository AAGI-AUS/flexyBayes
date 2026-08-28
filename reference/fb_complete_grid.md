# Complete a design grid before fitting

Reinstates the design cells that are absent from a data frame, returning
the Cartesian product of the index variables with the response set to
`NA` on every cell that was added. A field book already records a lost
plot this way – one row per sown plot, with the design columns filled in
and the yield missing – so this helper is for the trials where the plot
is missing from the file altogether.

## Usage

``` r
fb_complete_grid(data, index = ~row * col, response, unused_level = NULL)
```

## Arguments

- data:

  A data frame holding one row per recorded design cell.

- index:

  A one-sided formula naming the design index variables, for example
  `~ row * col`. The operators are ignored; what is read is the set of
  variables, whose observed levels are crossed.

- response:

  Name of the response column, as a single string. The reinstated cells
  carry `NA` there.

- unused_level:

  A single string, or `NULL` (the default). When supplied, design
  factors that vary across the trial are filled with this level on the
  reinstated cells and a warning names every column so filled. `NULL`
  refuses instead.

## Value

A data frame: `data` with one row appended per absent combination of the
index variables, the response `NA` on each, and the original row order
preserved ahead of them. Returned unchanged when the grid is already
complete.

## Details

The completed frame is what `na_action = "augment"` (the default) then
carries through the fit, keeping the index set a structured covariance
is built over intact. A separable `ar1(row):ar1(col)` field needs one
observation per node of the array, which is the same requirement ASReml
states for `residual = ~ ar1:ar1`.

A cell can only be reinstated when every other model column at it is
determined: it is one of the index variables, or it takes one value
across the whole trial. A column that varies – a variety label, a block
– is a real quantity nobody recorded at the absent cell, and the default
is to refuse rather than write one in. `unused_level` opens that door
explicitly, filling the varying design factors with a level of your
naming and warning with every column it wrote to. It is ASReml's nin89
LANCER coding, made a decision with a name on it rather than a default.
A varying *numeric* column is refused whatever `unused_level` says: a
measured covariate has no level to name.

## See also

[`flexybayes()`](https://aagi-aus.github.io/flexyBayes/reference/flexybayes.md),
whose `na_action` argument carries the completed grid into the fit.

## Examples

``` r
trial <- expand.grid(row = factor(1:4), col = factor(1:3))
trial$yield <- rnorm(nrow(trial))
holed <- trial[-c(2L, 7L), ]
padded <- fb_complete_grid(holed, ~ row * col, response = "yield")
nrow(padded)
#> [1] 12
sum(is.na(padded$yield))
#> [1] 2
```
