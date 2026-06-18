# List flexyBayes refusal reasons

`fb_refusals()` exposes the locked refusal-reason registry as a
browsable table: the canonical vocabulary of conditions under which
flexyBayes declines to fit, route, or validate a model — each with a
one-line description and the release it was introduced in. It is the
discovery surface for the structured refusals the package raises. Every
such refusal carries a condition class
`flexybayes_refusal_<reason_code>`, so a reason listed here can be
caught precisely, for example with
`tryCatch(fit, flexybayes_refusal_precision_not_symmetric = handler)`.

## Usage

``` r
fb_refusals(reason_code = NULL, since_version = NULL)
```

## Arguments

- reason_code:

  Optional character vector of exact reason codes to filter to. `NULL`
  (default) returns all registered reasons.

- since_version:

  Optional single version-string prefix to filter to. `NULL` (default)
  returns all.

## Value

A data frame of subclass `fb_refusals_table`, one row per matching
refusal reason, with columns `reason_code`, `description`,
`since_version`, and `plan_field`. The print method renders it as a
compact checklist.

## Details

Two optional filters narrow the listing. `reason_code` selects rows by
exact reason-code match (a single code or a vector). The `since_version`
filter selects reasons introduced in a matching release by
version-string prefix — `since_version = "0.4"` returns every reason
added in the 0.4 series.

Routing-decision reasons (surfaced by
[`fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/fb_plan.md)
and
[`backend_decision()`](https://aagi-aus.github.io/flexyBayes/reference/backend_decision.md))
and internal control-flow signals are deliberately excluded: this table
lists only refusals a user can actually encounter.

## Examples

``` r
fb_refusals()
#> flexyBayes refusal registry: 54 reasons
#> 
#>   [since 0.3.10] approximate_route_not_yet_registered
#>       Approximate covariance / dispatch carriers refuse until an
#>       approximation scheme is registered (validate_approximation()).
#> 
#>   [since 0.4.0] approximation_absent
#>       validate_approximation(): the low_rank_smooth scheme is
#>       registered but no smooth term was routed through the truncation
#>       path on this fit.
#> 
#>   [since 0.4.0] approximation_no_smooth_path
#>       s(..., representation = ): the named scheme is registered but has
#>       no smooth-basis emit path at this release.
#> 
#>   [since 0.4.0] approximation_scheme_unknown
#>       validate_approximation(): the fit carries no recognised
#>       approximation to validate (it is exact).
#> 
#>   [since 0.4.0] approximation_spec_invalid
#>       s(..., representation = ): the representation spec is not a list
#>       / fb_approx() carrying a single-string scheme.
#> 
#>   [since 0.3.10] block_not_positive_definite
#>       Block-diagonal vm/ped: at least one V_k failed the
#>       positive-definite probe.
#> 
#>   [since 0.3.10] block_partition_incomplete
#>       Block-diagonal vm/ped: the block sizes do not partition the
#>       grouping factor's level count.
#> 
#>   [since 0.4.0] blocks_empty_list
#>       vm(..., blocks = ): the block list is empty.
#> 
#>   [since 0.4.0] blocks_not_a_list
#>       vm(..., blocks = ): the block carrier is not a base-R list of
#>       covariance matrices.
#> 
#>   [since 0.4.0] blocks_not_in_known_matrices
#>       vm(..., blocks = ): the named block list is absent from
#>       known_matrices.
#> 
#>   [since 0.4.0] chol_not_in_known_matrices
#>       vm(..., chol = ): the named Cholesky factor is absent from
#>       known_matrices.
#> 
#>   [since 0.4.0] chol_not_square
#>       vm(..., chol = ): the Cholesky factor is not square.
#> 
#>   [since 0.4.0] chol_not_triangular
#>       vm(..., chol = ): the Cholesky factor is not lower-triangular.
#> 
#>   [since 0.4.0] code_flags_mutually_exclusive
#>       return_code = TRUE and review_code = TRUE were both supplied; the
#>       two code-return modes are mutually exclusive.
#> 
#>   [since 0.4.0] cov_arg_not_fb_cov
#>       vm() / ped(): the `cov` argument must be written inline as an
#>       fb_cov() carrier.
#> 
#>   [since 0.4.0] design_memory_exceeds_ceiling
#>       Preflight: the design matrix is estimated to exceed the active
#>       memory ceiling; dispatch is short-circuited before any backend
#>       code runs.
#> 
#>   [since 0.5.0] engine_pin_backend_conflict
#>       An engine pin (fb_greta / fb_inla / fb_brms) was given a
#>       `backend` argument that conflicts with the engine it pins.
#> 
#>   [since 0.7.0] fa_rank_exceeds_dim
#>       A factor-analytic term fa(x, k) was given a rank k that is not
#>       strictly below the number of levels of the outer factor. A
#>       factor-analytic covariance is identifiable only for k < n_outer:
#>       at k = n_outer the loadings and specific variances form an
#>       over-parameterised reparameterisation of the unstructured form,
#>       and at k > n_outer the lower-triangular loadings carry empty
#>       columns. This is a data-aware preflight (n_outer is known only
#>       after the term is matched against the data), complementing the
#>       data-free fa_rank_invalid (k < 1) check.
#> 
#>   [since 0.4.0] fa_rank_invalid
#>       A factor-analytic term fa(x, k) was given a rank k below 1; the
#>       factor-analytic rank must be a positive integer.
#> 
#>   [since 0.4.0] fb_cov_missing_matrix
#>       fb_cov(): the carrier matrix `M` (the first argument) was not
#>       supplied.
#> 
#>   [since 0.4.0] fb_cov_type_unknown
#>       fb_cov(): the requested carrier `type` is not one of the five
#>       known types (dense / chol / precision / blocks / low_rank).
#> 
#>   [since 0.4.0] formula_not_two_sided
#>       The model formula must be two-sided (response ~ predictors); a
#>       formula carrying no left-hand-side response was supplied.
#> 
#>   [since 0.4.1] grammar_brms_known_matrices_unsupported
#>       `known_matrices` was supplied with brms-grammar ingest via the
#>       universal entry, which has no known-matrix carrier.
#> 
#>   [since 0.4.1] grammar_brms_with_asreml_terms
#>       A brms-style bar-grouped formula was combined with ASReml
#>       `random` / `rcov` arguments on the universal entry.
#> 
#>   [since 0.6.0.9000] gretaR_below_version_floor
#>       gretaR backend: the installed gretaR is older than the activation
#>       floor.
#> 
#>   [since 0.6.0.9000] gretaR_cannot_represent_structured_cov
#>       gretaR backend: structured covariance (vm/ped/fa/us/ar1)
#>       unsupported.
#> 
#>   [since 0.6.0.9000] gretaR_family_unsupported
#>       gretaR backend: family outside gaussian/binomial/poisson.
#> 
#>   [since 0.6.0.9000] gretaR_not_installed
#>       gretaR backend: gretaR not installed and no source home set.
#> 
#>   [since 0.6.0.9000] gretaR_random_group_not_in_data
#>       gretaR backend: random-intercept grouping factor absent from
#>       data.
#> 
#>   [since 0.6.0.9000] gretaR_random_term_type_unsupported
#>       gretaR backend: only random-intercept-class random terms are
#>       supported.
#> 
#>   [since 0.4.0] heterogeneous_residual_factor_not_in_cell_key
#>       Aggregated Gaussian emit: an at(f):units heterogeneous residual
#>       factor is not in the cell key, so the cell-constant sigma
#>       property does not hold.
#> 
#>   [since 0.4.0] known_matrices_data_name_collision
#>       INLA emit: a known-matrices / blocks carrier name collides with a
#>       data column name.
#> 
#>   [since 0.4.0] known_matrix_dim_mismatch
#>       vm(): the known matrix dimension does not match the grouping
#>       factor's level count.
#> 
#>   [since 0.4.0] known_matrix_dimnames_mismatch
#>       vm(): the known matrix has differing row and column names.
#> 
#>   [since 0.4.0] known_matrix_level_mismatch
#>       vm(): the known matrix dimnames do not match (or are mis-ordered
#>       relative to) the grouping factor levels.
#> 
#>   [since 0.4.0] low_rank_rank_exceeds_basis
#>       low_rank_smooth: the requested rank meets or exceeds the
#>       truncation ceiling min(basis dimension k, n) and so is not an
#>       approximation.
#> 
#>   [since 0.4.0] low_rank_rank_invalid
#>       low_rank_smooth: the requested rank is not a single positive
#>       integer.
#> 
#>   [since 0.4.0] low_rank_requires_greta
#>       A smooth requesting the low_rank_smooth approximation was routed
#>       to a non-greta backend that cannot honour it.
#> 
#>   [since 0.4.0] low_rank_scheme_required
#>       vm()/ped(): low_rank_factor supplied without an explicit
#>       low_rank_scheme naming a registered approximation.
#> 
#>   [since 0.5.0] native_greta_requires_greta_backend
#>       A native greta model graph was passed to the universal entry /
#>       the greta pin with a non-greta backend. A native graph is
#>       greta-only by construction.
#> 
#>   [since 0.4.0] precision_not_in_known_matrices
#>       vm(..., precision = ): the named precision matrix is absent from
#>       known_matrices.
#> 
#>   [since 0.4.0] precision_not_positive_definite
#>       vm(..., precision = ): the precision matrix failed the
#>       positive-definite probe.
#> 
#>   [since 0.4.0] precision_not_square
#>       vm(..., precision = ): the precision matrix is not square.
#> 
#>   [since 0.4.0] precision_not_symmetric
#>       vm(..., precision = ): the precision matrix is not symmetric.
#> 
#>   [since 0.4.0] predict_kernel_invalid_include
#>       predict(): `include` is empty or carries values outside the
#>       prediction-kernel vocabulary.
#> 
#>   [since 0.4.0] rcov_type_unsupported_for_aggregation
#>       Aggregated Gaussian emit: the rcov term type is outside the
#>       supported aggregation scope.
#> 
#>   [since 0.4.0] representation_unknown_for_preflight
#>       Preflight: the design representation is not characterised by the
#>       preflight memory estimator.
#> 
#>   [since 0.4.0] response_not_in_data
#>       The response variable named on the formula's left-hand side is
#>       not a column of `data`.
#> 
#>   [since 0.4.0] review_code_backend_unsupported
#>       review_code = TRUE was requested with a backend other than greta;
#>       the inspect-then-fit token is currently greta-only.
#> 
#>   [since 0.4.0] smooth_variable_not_in_data
#>       The variable inside a smooth term s(x) is not a column of `data`.
#> 
#>   [since 0.5.0] stan_cannot_represent_structured_cov
#>       backend = "brms" (Stan) cannot represent an asreml
#>       structured-covariance term (vm/ped/fa/us/ar1).
#> 
#>   [since 0.4.0] tensor_smooth_unsupported
#>       A tensor-product or multivariate smooth (te(), ti(), t2()) was
#>       supplied. flexyBayes fits univariate penalised splines (s(),
#>       spl()) only.
#> 
#>   [since 0.4.0] unsupported_family
#>       The requested family is outside the set flexyBayes can emit.
#>       Refused at the family gate (.resolve_family) before any backend
#>       code runs.
#> 
#>   [since 0.4.0] vm_redundant_specification
#>       vm()/ped(): more than one covariance carrier supplied; exactly
#>       one of V / chol / precision / blocks / low_rank_factor is
#>       allowed.
#> 
fb_refusals(reason_code = "precision_not_symmetric")
#> flexyBayes refusal registry: 1 reason  (filter: reason_code in {precision_not_symmetric})
#> 
#>   [since 0.4.0] precision_not_symmetric
#>       vm(..., precision = ): the precision matrix is not symmetric.
#> 
fb_refusals(since_version = "0.4")
#> flexyBayes refusal registry: 41 reasons  (filter: since_version ~ '0.4')
#> 
#>   [since 0.4.0] approximation_absent
#>       validate_approximation(): the low_rank_smooth scheme is
#>       registered but no smooth term was routed through the truncation
#>       path on this fit.
#> 
#>   [since 0.4.0] approximation_no_smooth_path
#>       s(..., representation = ): the named scheme is registered but has
#>       no smooth-basis emit path at this release.
#> 
#>   [since 0.4.0] approximation_scheme_unknown
#>       validate_approximation(): the fit carries no recognised
#>       approximation to validate (it is exact).
#> 
#>   [since 0.4.0] approximation_spec_invalid
#>       s(..., representation = ): the representation spec is not a list
#>       / fb_approx() carrying a single-string scheme.
#> 
#>   [since 0.4.0] blocks_empty_list
#>       vm(..., blocks = ): the block list is empty.
#> 
#>   [since 0.4.0] blocks_not_a_list
#>       vm(..., blocks = ): the block carrier is not a base-R list of
#>       covariance matrices.
#> 
#>   [since 0.4.0] blocks_not_in_known_matrices
#>       vm(..., blocks = ): the named block list is absent from
#>       known_matrices.
#> 
#>   [since 0.4.0] chol_not_in_known_matrices
#>       vm(..., chol = ): the named Cholesky factor is absent from
#>       known_matrices.
#> 
#>   [since 0.4.0] chol_not_square
#>       vm(..., chol = ): the Cholesky factor is not square.
#> 
#>   [since 0.4.0] chol_not_triangular
#>       vm(..., chol = ): the Cholesky factor is not lower-triangular.
#> 
#>   [since 0.4.0] code_flags_mutually_exclusive
#>       return_code = TRUE and review_code = TRUE were both supplied; the
#>       two code-return modes are mutually exclusive.
#> 
#>   [since 0.4.0] cov_arg_not_fb_cov
#>       vm() / ped(): the `cov` argument must be written inline as an
#>       fb_cov() carrier.
#> 
#>   [since 0.4.0] design_memory_exceeds_ceiling
#>       Preflight: the design matrix is estimated to exceed the active
#>       memory ceiling; dispatch is short-circuited before any backend
#>       code runs.
#> 
#>   [since 0.4.0] fa_rank_invalid
#>       A factor-analytic term fa(x, k) was given a rank k below 1; the
#>       factor-analytic rank must be a positive integer.
#> 
#>   [since 0.4.0] fb_cov_missing_matrix
#>       fb_cov(): the carrier matrix `M` (the first argument) was not
#>       supplied.
#> 
#>   [since 0.4.0] fb_cov_type_unknown
#>       fb_cov(): the requested carrier `type` is not one of the five
#>       known types (dense / chol / precision / blocks / low_rank).
#> 
#>   [since 0.4.0] formula_not_two_sided
#>       The model formula must be two-sided (response ~ predictors); a
#>       formula carrying no left-hand-side response was supplied.
#> 
#>   [since 0.4.1] grammar_brms_known_matrices_unsupported
#>       `known_matrices` was supplied with brms-grammar ingest via the
#>       universal entry, which has no known-matrix carrier.
#> 
#>   [since 0.4.1] grammar_brms_with_asreml_terms
#>       A brms-style bar-grouped formula was combined with ASReml
#>       `random` / `rcov` arguments on the universal entry.
#> 
#>   [since 0.4.0] heterogeneous_residual_factor_not_in_cell_key
#>       Aggregated Gaussian emit: an at(f):units heterogeneous residual
#>       factor is not in the cell key, so the cell-constant sigma
#>       property does not hold.
#> 
#>   [since 0.4.0] known_matrices_data_name_collision
#>       INLA emit: a known-matrices / blocks carrier name collides with a
#>       data column name.
#> 
#>   [since 0.4.0] known_matrix_dim_mismatch
#>       vm(): the known matrix dimension does not match the grouping
#>       factor's level count.
#> 
#>   [since 0.4.0] known_matrix_dimnames_mismatch
#>       vm(): the known matrix has differing row and column names.
#> 
#>   [since 0.4.0] known_matrix_level_mismatch
#>       vm(): the known matrix dimnames do not match (or are mis-ordered
#>       relative to) the grouping factor levels.
#> 
#>   [since 0.4.0] low_rank_rank_exceeds_basis
#>       low_rank_smooth: the requested rank meets or exceeds the
#>       truncation ceiling min(basis dimension k, n) and so is not an
#>       approximation.
#> 
#>   [since 0.4.0] low_rank_rank_invalid
#>       low_rank_smooth: the requested rank is not a single positive
#>       integer.
#> 
#>   [since 0.4.0] low_rank_requires_greta
#>       A smooth requesting the low_rank_smooth approximation was routed
#>       to a non-greta backend that cannot honour it.
#> 
#>   [since 0.4.0] low_rank_scheme_required
#>       vm()/ped(): low_rank_factor supplied without an explicit
#>       low_rank_scheme naming a registered approximation.
#> 
#>   [since 0.4.0] precision_not_in_known_matrices
#>       vm(..., precision = ): the named precision matrix is absent from
#>       known_matrices.
#> 
#>   [since 0.4.0] precision_not_positive_definite
#>       vm(..., precision = ): the precision matrix failed the
#>       positive-definite probe.
#> 
#>   [since 0.4.0] precision_not_square
#>       vm(..., precision = ): the precision matrix is not square.
#> 
#>   [since 0.4.0] precision_not_symmetric
#>       vm(..., precision = ): the precision matrix is not symmetric.
#> 
#>   [since 0.4.0] predict_kernel_invalid_include
#>       predict(): `include` is empty or carries values outside the
#>       prediction-kernel vocabulary.
#> 
#>   [since 0.4.0] rcov_type_unsupported_for_aggregation
#>       Aggregated Gaussian emit: the rcov term type is outside the
#>       supported aggregation scope.
#> 
#>   [since 0.4.0] representation_unknown_for_preflight
#>       Preflight: the design representation is not characterised by the
#>       preflight memory estimator.
#> 
#>   [since 0.4.0] response_not_in_data
#>       The response variable named on the formula's left-hand side is
#>       not a column of `data`.
#> 
#>   [since 0.4.0] review_code_backend_unsupported
#>       review_code = TRUE was requested with a backend other than greta;
#>       the inspect-then-fit token is currently greta-only.
#> 
#>   [since 0.4.0] smooth_variable_not_in_data
#>       The variable inside a smooth term s(x) is not a column of `data`.
#> 
#>   [since 0.4.0] tensor_smooth_unsupported
#>       A tensor-product or multivariate smooth (te(), ti(), t2()) was
#>       supplied. flexyBayes fits univariate penalised splines (s(),
#>       spl()) only.
#> 
#>   [since 0.4.0] unsupported_family
#>       The requested family is outside the set flexyBayes can emit.
#>       Refused at the family gate (.resolve_family) before any backend
#>       code runs.
#> 
#>   [since 0.4.0] vm_redundant_specification
#>       vm()/ped(): more than one covariance carrier supplied; exactly
#>       one of V / chol / precision / blocks / low_rank_factor is
#>       allowed.
#> 
```
