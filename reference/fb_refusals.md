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
#> flexyBayes refusal registry: 114 reasons
#> 
#>   [since 0.9.2] aggregation_response_incomplete
#>       aggregate = TRUE with missing values in the response. The
#>       aggregated path compresses rows into per-cell sufficient
#>       statistics, which are not defined for a missing response, so the
#>       request has no aggregated reading. Untyped before 0.9.2.
#> 
#>   [since 0.9.2] aggregation_route_unavailable
#>       aggregate = TRUE with no aggregated emit route available --
#>       either the named backend carries no aggregated path at all, or no
#>       aggregated backend resolves for the model. Raised at two sites
#>       with site-specific messages. Untyped before 0.9.2.
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
#>   [since 0.9.0] ar1_residual_not_representable
#>       A separable AR1 process written on the residual side. INLA
#>       represents it as a latent random field plus the observation-level
#>       Gaussian noise -- four parameters -- and cannot fit ASReml's
#>       three-parameter nugget-free residual, so the residual spelling is
#>       refused rather than fitted under a name it does not have.
#> 
#>   [since 0.9.0] ar1_spatial_requires_complete_grid
#>       A separable ar1(row):ar1(col) spatial term on INLA is emitted as
#>       an AR1xAR1 latent field, faithful only with one observation per
#>       grid node. The data are an incomplete or replicated grid.
#> 
#>   [since 0.9.0] ar2_not_representable
#>       ar2() parses for the formula catalogue and has no emit path on
#>       either active engine.
#> 
#>   [since 0.9.0] asreml_function_not_recognised
#>       A call in a random or residual formula whose function is outside
#>       the recognised ASReml vocabulary. Previously such a call was read
#>       as a plain variable named after its own source text.
#> 
#>   [since 0.9.3] at_field_per_level_hyper_not_representable
#>       at(trial, level):ar1(row):ar1(col) names a single conditioned
#>       level of the grouping factor for a per-trial separable AR1 field,
#>       or (the same lowering gap read the other way) asks for field
#>       hyperparameters that vary by level rather than being shared. The
#>       supported spelling, at(trial):ar1(row):ar1(col) with no level
#>       argument, fits one field per level of `trial` on INLA via the
#>       replicate = mechanism, sharing rho_row, rho_col and the field
#>       variance across every level; there is no lowering for a
#>       level-conditioned field or for per-level (unshared)
#>       hyperparameters.
#> 
#>   [since 0.9.0] at_level_conditioning_unsupported
#>       at(f, level):g conditions a random effect on selected levels of
#>       f, which no active emitter represents. It is NOT the
#>       heterogeneous variance structure diag(f):g, and is refused rather
#>       than fitted as one.
#> 
#>   [since 0.9.0] augment_cell_not_determinable
#>       A design cell is absent from the data and the model depends on a
#>       variable whose value at that cell was never recorded.
#> 
#>   [since 0.9.0] auto_no_active_route
#>       backend = "auto" found no active backend able to faithfully fit
#>       the model: INLA refused it and brms cannot represent it. No
#>       silent fallback -- auto refuses.
#> 
#>   [since 0.9.3] binomial_response_not_binary
#>       family = "binomial" on the main (non-streaming) entry requires
#>       the response to be a numeric vector of only 0 and 1: there is no
#>       cbind(success, failure) or trials = spelling on flexybayes() /
#>       fb_inla() / fb_brms() (that exists only on flexybayes_stream()).
#>       Before this refusal, a response outside {0, 1} reached the engine
#>       unrefused and failed there: INLA with a raw subprocess error
#>       ('Binomial data ... is void', wrapped into inla_program_failed
#>       with an unrelated design-memory diagnosis), brms with an untyped
#>       condition ('Family ... requires responses to contain only two
#>       different values'), and the aggregated route with an untyped
#>       stop(). This refusal fires once, before any of the three, naming
#>       the column, the offending values, and the remedy: trials = on
#>       flexybayes_stream() for pre-aggregated counts, or recode to a 0/1
#>       numeric vector for a per-row binary outcome.
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
#>   [since 0.9.0] brms_cannot_augment_nongaussian
#>       backend = "brms" cannot carry a missing response for a
#>       non-Gaussian family: its mi() addition term is Gaussian-only, and
#>       without it brms drops the rows.
#> 
#>   [since 0.9.0] brms_cannot_represent_term
#>       A random-effect term type reached the brms formula reconstruction
#>       with no branch to lower it. It is refused by name rather than
#>       dropped from the emitted model.
#> 
#>   [since 0.9.2] brms_factor_random_slope_unsupported
#>       A bar-grammar random slope on a factor -- (f | g), (0 + f | g),
#>       (f || g) -- is outside the brms ingest. The ASReml surface of the
#>       same model fits, so the refusal names it verbatim; before 0.9.2
#>       the message listed other bar spellings, none of which expresses
#>       the model.
#> 
#>   [since 0.9.2] brms_ingest_feature_unsupported
#>       A brms formula feature outside the ingest layer's scope (an
#>       extended term the walker collected but cannot lower). Untyped
#>       before 0.9.2.
#> 
#>   [since 0.9.2] brms_random_effect_form_unsupported
#>       A bar-grammar random-effect specification outside the supported
#>       set -- multi-variable slopes and group-of-coefficients shorthand.
#>       Untyped before 0.9.2.
#> 
#>   [since 0.9.3] cell_count_exceeds_integer
#>       The aggregation planner's estimated cell count (the product of
#>       the cell key's factor level counts) exceeds R's integer limit.
#>       Recording it would coerce it to NA, so the plan's cell count and
#>       compression estimate would both silently disappear. The row-count
#>       sibling (row_count_exceeds_integer) fixed the identical defect
#>       for N; this is the cell-count member of the same class (FS-24).
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
#>   [since 0.9.1] classify_random_factor_not_supported
#>       predict(classify = ) was asked for a factor that enters the model
#>       only as a random-effects grouping term. Marginal means are built
#>       over the population-level reference grid, which such a factor is
#>       not part of.
#> 
#>   [since 0.9.1] classify_requires_emmeans
#>       predict(classify = ) builds its marginal-means table through the
#>       emmeans estimation seam, and emmeans is a suggested package
#>       rather than a hard dependency.
#> 
#>   [since 0.4.0] code_flags_mutually_exclusive
#>       return_code = TRUE and review_code = TRUE were both supplied; the
#>       two code-return modes are mutually exclusive.
#> 
#>   [since 0.9.0] conditional_loglik_not_available
#>       No plug-in conditional log-likelihood exists for this fit -- the
#>       family has no implemented form, the response or family record is
#>       absent, or the engine reports a marginal likelihood instead.
#>       Returning NA would let AIC() and anova() present a comparison
#>       that was never computed.
#> 
#>   [since 0.9.0] corh_no_equicorrelation_representation
#>       corh(f):g requests heterogeneous variances with a single shared
#>       correlation. No active backend represents an equicorrelation
#>       group-level structure; brms offers uncorrelated or fully
#>       unstructured and nothing between.
#> 
#>   [since 0.4.0] cov_arg_not_fb_cov
#>       vm() / ped(): the `cov` argument must be written inline as an
#>       fb_cov() carrier.
#> 
#>   [since 0.9.1] covariate_zero_fill_not_supported
#>       A predictor has missing values under the covariate policy
#>       `include`, which in ASReml means treating the missing covariate
#>       as zero (Reference Manual 4.2, section 3.11).
#> 
#>   [since 0.4.0] design_memory_exceeds_ceiling
#>       Preflight: the design matrix is estimated to exceed the active
#>       memory ceiling; dispatch is short-circuited before any backend
#>       code runs.
#> 
#>   [since 0.9.0] dsum_structured_inner_unsupported
#>       Residual dsum() with a structured inner term (ar1, us, ...).
#>       flexyBayes represents dsum() only as a per-region heteroscedastic
#>       variance, so a structured inner is refused rather than silently
#>       reduced to a heteroscedastic residual.
#> 
#>   [since 0.5.0] engine_pin_backend_conflict
#>       An engine pin (fb_inla / fb_brms) was given a `backend` argument
#>       that conflicts with the engine it pins.
#> 
#>   [since 0.9.0] fa_not_representable
#>       A factor-analytic covariance, fa(f, k) or fa(f, k):g. It parses
#>       for the formula catalogue and no active engine emits a
#>       factor-analytic covariance.
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
#>   [since 0.9.0] family_argument_not_recognised
#>       The `family` argument was neither a single family name nor a
#>       usable `stats` family object, or a family object and an explicit
#>       `link =` named different links. A family object supplies its own
#>       link, so resolving the disagreement either way would fit a model
#>       the call did not ask for.
#> 
#>   [since 0.4.0] fb_cov_missing_matrix
#>       fb_cov(): the carrier matrix `M` (the first argument) was not
#>       supplied.
#> 
#>   [since 0.4.0] fb_cov_type_unknown
#>       fb_cov(): the requested carrier `type` is not one of the five
#>       known types (dense / chol / precision / blocks / low_rank).
#> 
#>   [since 0.9.0] fit_lacks_posterior_draws
#>       A method needing the posterior itself -- draws on a sampling
#>       engine, marginal densities on INLA -- was called on a fit that
#>       carries neither. The refusal replaces an empty result that would
#>       read as an answer. Its sibling is loo_requires_sampler_draws,
#>       raised when the fit does carry a posterior but not the pointwise
#>       log-likelihood draws that PSIS-LOO reads.
#> 
#>   [since 0.9.2] fixed_smoother_not_supported
#>       An mgcv-style smoother (s(), te(), ti(), t2(), gp(), sos())
#>       appears in the fixed part of the formula. flexyBayes spells a
#>       univariate smooth `random = ~ spl(x)` on the INLA backend; before
#>       0.9.2 the term reached the engine and surfaced as `could not find
#>       function "s"`, on brms only after a completed sampling run.
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
#>       `random` / `residual` arguments on the universal entry.
#> 
#>   [since 0.4.0] heterogeneous_residual_factor_not_in_cell_key
#>       Aggregated Gaussian emit: an at(f):units heterogeneous residual
#>       factor is not in the cell key, so the cell-constant sigma
#>       property does not hold.
#> 
#>   [since 0.9.0] heterogeneous_residual_family_has_no_sigma
#>       A heterogeneous residual was requested for a family with no
#>       residual scale parameter (Poisson, binomial and relatives), whose
#>       dispersion is determined by the mean.
#> 
#>   [since 0.9.0] heterogeneous_residual_multiple_factors
#>       More than one heterogeneous-residual term was supplied. A
#>       residual is sectioned by one factor; nesting several would
#>       require their interaction, which the user must state explicitly.
#> 
#>   [since 0.9.0] inla_gate_refused
#>       An explicit backend = "inla" request that lgm_gate() refused. The
#>       gate's structural check-list is the message and rides on the
#>       condition in its `lgm_refusal` slot; the failing rule supplies a
#>       second condition class so one rule can be caught on its own.
#> 
#>   [since 0.9.3] inla_program_failed
#>       The INLA engine died -- a non-zero exit from the inla-program
#>       subprocess, an inla.core.safe failure, or an empty/NULL result
#>       with no R-level error at all -- past the design-memory preflight,
#>       which models per-term storage but not the sparse-Cholesky solver
#>       cost of a large latent field. Before 0.9.3 this reached the user
#>       as a raw pass-through of the engine's own message after tens of
#>       minutes of runtime, with no reason code and no remedy (FS-25).
#>       Carries the design size, the largest design this package has
#>       verified an INLA fit to complete (read from
#>       inst/validation/benchmark_scaling.csv), the binding random-effect
#>       term where determinable, and the remedies.
#> 
#>   [since 0.9.2] inla_variable_used_twice
#>       The same variable indexes both a fixed-effect term and a latent
#>       f() term in the emitted INLA formula. INLA refuses a duplicated
#>       key; before 0.9.2 the refusal surfaced as the generic
#>       inla-program-exited message.
#> 
#>   [since 0.9.0] interaction_not_representable
#>       An interaction of two structured terms that matches none of the
#>       recognised interaction patterns -- us(trait):vm(gen), say. It
#>       parsed into the catalogue's generic-interaction node, which no
#>       engine emits.
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
#>   [since 0.9.1] loo_requires_sampler_draws
#>       loo() was called on a fit whose engine does not store the
#>       pointwise log-likelihood draws PSIS-LOO reads. Sibling of
#>       fit_lacks_posterior_draws, which is raised when the fit carries
#>       no posterior at all.
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
#>   [since 0.4.0] low_rank_scheme_required
#>       vm()/ped(): low_rank_factor supplied without an explicit
#>       low_rank_scheme naming a registered approximation.
#> 
#>   [since 0.4.0] low_rank_smooth_unsupported
#>       A smooth requesting the low_rank_smooth approximation was routed
#>       to a backend that cannot honour it. No active engine represents
#>       the rank-K basis truncation, so the request refuses under `auto`
#>       as well as under an explicit `inla` / `brms`.
#> 
#>   [since 0.9.0] met_summary_not_available
#>       fb_met_summary() was called on a fit it cannot summarise: the
#>       breeder summary reads realised factor-analytic effects, and no
#>       active engine emits an fa() term, so this refusal is
#>       unconditional.
#> 
#>   [since 0.9.0] missing_covariate_not_supported
#>       A predictor has missing values under the covariate policy `fail`,
#>       which is both this package's default and ASReml's.
#> 
#>   [since 0.9.0] missing_response_refused
#>       The response has missing values and na_action = "fail" was
#>       requested.
#> 
#>   [since 0.9.0] model_matrix_not_recoverable
#>       model.matrix() was called on a fit carrying no fixed-effect
#>       formula or no fitted data, so the population-level design cannot
#>       be rebuilt. Base R would silently return an intercept-only matrix
#>       built from the calling frame.
#> 
#>   [since 0.9.0] numeric_variable_in_random_interaction
#>       A numeric variable inside a nested or multi-way random
#>       interaction. The spelling reads as a random regression to an
#>       ASReml user and as one independent deviation per cell to the
#>       emitted model; flexyBayes refuses rather than choosing one of the
#>       two readings.
#> 
#>   [since 0.9.1] pp_check_requires_predictive_draws
#>       pp_check() was called on a fit whose engine returns no posterior
#>       predictive draws, so there are no replicated datasets to overlay
#>       on the observed response.
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
#>   [since 0.9.2] prior_argument_duplicated
#>       A prior distribution call binds the same hyperparameter twice, so
#>       which of the two values applies is not decidable from the call.
#> 
#>   [since 0.9.2] prior_argument_missing
#>       A prior distribution call omits a hyperparameter the family
#>       requires. A half-specified PC prior is the sharp case: half a
#>       probability statement is not a probability statement, and the
#>       completed value reached the engine unannounced.
#> 
#>   [since 0.9.2] prior_argument_unknown
#>       A prior distribution call carries an argument name the family
#>       does not have. Before 0.9.2 the call fell through unmatched and
#>       every hyperparameter then took its emit-time default, so a
#>       misspelt argument became a silently substituted prior.
#> 
#>   [since 0.9.2] prior_distribution_not_a_call
#>       The right-hand side of a prior formula is not a distribution
#>       call, so it carries no hyperparameters.
#> 
#>   [since 0.9.2] prior_distribution_unknown
#>       The distribution named on the right-hand side of a prior formula
#>       is outside the DSL's supported set.
#> 
#>   [since 0.9.2] prior_hyperparameter_not_scalar
#>       A prior hyperparameter is not a finite numeric scalar. A
#>       zero-length or length-two value used to surface as one of R's own
#>       interpreter messages several layers downstream.
#> 
#>   [since 0.9.2] prior_hyperparameter_out_of_domain
#>       A prior hyperparameter lies outside the domain its family defines
#>       -- a non-positive scale, rate, shape or degrees of freedom, a
#>       tail probability outside (0, 1), or reversed uniform bounds. A
#>       negative standard deviation used to reach the sampler and surface
#>       as a missing-draws error.
#> 
#>   [since 0.9.2] prior_not_translatable_for_backend
#>       The selected backend has no representation for a supplied prior
#>       specification. Before 0.9.2 the INLA route discarded such a row
#>       silently while prior_summary() printed it as applied, and the
#>       brms route refused it with an untyped error.
#> 
#>   [since 0.9.2] prior_spec_empty
#>       fb_prior() was called with no specification. A prior object with
#>       no targets would silently leave every parameter on its engine
#>       default.
#> 
#>   [since 0.9.2] prior_spec_not_formula
#>       An fb_prior() argument is not a formula. The DSL's only argument
#>       shape is `target ~ distribution(...)`.
#> 
#>   [since 0.9.2] prior_spec_not_two_sided
#>       An fb_prior() formula is one-sided, so it names a distribution
#>       without naming the parameter it applies to.
#> 
#>   [since 0.9.2] prior_target_argument_missing
#>       A prior target call is missing the argument that identifies which
#>       parameter it targets -- sd() and cor() need `group`, b() and
#>       smooth() need a name.
#> 
#>   [since 0.9.2] prior_target_not_in_model
#>       A prior names a coefficient or a variance component the model
#>       does not have. The check reads the engine's own parameter list,
#>       so the refusal names the model's actual terms rather than a
#>       synthesised Stan parameter name.
#> 
#>   [since 0.9.2] prior_target_unsupported
#>       The left-hand side of a prior formula is not one of the supported
#>       targets. A bare name that is not `sigma` used to be carried as an
#>       opaque target and dropped at emit time.
#> 
#>   [since 0.4.0] representation_unknown_for_preflight
#>       Preflight: the design representation is not characterised by the
#>       preflight memory estimator.
#> 
#>   [since 0.4.0] residual_type_unsupported_for_aggregation
#>       Aggregated Gaussian emit: the residual term type is outside the
#>       supported aggregation scope.
#> 
#>   [since 0.4.0] response_not_in_data
#>       The response variable named on the formula's left-hand side is
#>       not a column of `data`.
#> 
#>   [since 0.4.0] review_code_backend_unsupported
#>       review_code = TRUE was requested on a backend with no code to
#>       review. The inspect-then-fit token is scoped to the code-emitting
#>       engines, which among the active pair is brms.
#> 
#>   [since 0.9.0] row_count_exceeds_integer
#>       A single-file or in-memory source reports more rows than R's
#>       integer type holds. Recording the count would coerce it to NA, so
#>       the aggregation plan and the fitted object would carry no row
#>       count. Partitioned .fst input keeps its total as a double and has
#>       no such ceiling.
#> 
#>   [since 0.4.0] smooth_variable_not_in_data
#>       The variable inside a smooth term s(x) is not a column of `data`.
#> 
#>   [since 0.9.0] stan_cannot_represent_ar1_field
#>       A random-side AR1 or separable AR1xAR1 latent field on brms. The
#>       field is a Kronecker autoregressive precision that only the INLA
#>       emit builds; brms has no lowering for it.
#> 
#>   [since 0.9.0] stan_cannot_represent_ar1_residual
#>       backend = "brms" (Stan) cannot represent an AR1 or separable
#>       AR1xAR1 residual covariance.
#> 
#>   [since 0.5.0] stan_cannot_represent_structured_cov
#>       backend = "brms" (Stan) cannot represent an asreml
#>       structured-covariance term (vm/ped/fa/us/ar1).
#> 
#>   [since 0.9.0] stan_cannot_represent_structured_residual
#>       backend = "brms" (Stan) cannot represent a structured residual
#>       covariance outside the iid `units` form.
#> 
#>   [since 0.9.0] str_not_representable
#>       str() declares a covariance structure over a group of terms. It
#>       parses for the formula catalogue and has no emit path on either
#>       active engine.
#> 
#>   [since 0.4.0] tensor_smooth_unsupported
#>       A tensor-product or multivariate smooth (te(), ti(), t2()) was
#>       supplied. flexyBayes fits univariate penalised splines (s(),
#>       spl()) only.
#> 
#>   [since 0.9.1] term_in_fixed_and_random
#>       The same factor appears as a fixed main effect and as a random
#>       main effect, which aliases the random deviations against the
#>       fixed level means: they are identified only by their own prior.
#> 
#>   [since 0.9.0] triangulate_incomparable_fits
#>       Two fits handed to triangulate() whose model fingerprints
#>       disagree -- a different formula, family, link, dataset or
#>       recorded prior. The distance between their posteriors would
#>       measure the difference between two questions rather than between
#>       two engines.
#> 
#>   [since 0.9.3] unknown_backend
#>       An explicit `backend` request named a value that is not a
#>       registered engine -- including a name this release withdrew
#>       entirely (see NEWS.md). The active engines are inla and brms.
#> 
#>   [since 0.4.0] unsupported_family
#>       The requested family is outside the set flexyBayes can emit.
#>       Refused at the family gate (.resolve_family) before any backend
#>       code runs.
#> 
#>   [since 0.9.0] update_call_not_reconstructable
#>       update() was called on a fit whose recorded argument set is
#>       incomplete. Re-fitting would substitute package defaults for the
#>       unrecorded arguments -- a relationship matrix or a prior scale
#>       among them -- and return a different model under the same name.
#> 
#>   [since 0.9.3] update_unnamed_argument_not_supported
#>       update() re-issues every recorded flexybayes() argument by name,
#>       matching an override in ... to it by name too; an unnamed
#>       argument -- stats::update()'s classic replacement-formula idiom,
#>       update(fit, . ~ . + z) -- has no name to match and was, before
#>       this refusal, silently discarded: the re-fit ran with the
#>       unchanged recorded call, a valid object answering a question
#>       nobody asked rather than a re-fit or a refusal. Name the slot
#>       instead: update(fit, fixed = y ~ x + z) (brms-style), or random =
#>       ~ ... / residual = ~ ... on the ASReml grammar.
#> 
#>   [since 0.9.1] variogram_requires_design_index
#>       plot(type = "variogram") needs the row / column array a
#>       structured covariance is indexed by, and the fit's terms name
#>       none.
#> 
#>   [since 0.4.0] vm_redundant_specification
#>       vm()/ped(): more than one covariance carrier supplied; exactly
#>       one of V / chol / precision / blocks / low_rank_factor is
#>       allowed.
#> 
#>   [since 0.9.3] weights_not_aggregatable
#>       aggregate = TRUE with non-constant observation weights. Weights
#>       are lowered for the per-row Gaussian route only; the aggregated
#>       route's sufficient statistics (per-cell sums and counts) do not
#>       fold in a per-observation weight, so a weighted aggregated fit is
#>       refused rather than silently ignoring the weight or answering a
#>       different question under the same name.
#> 
#>   [since 0.9.0] weights_not_supported
#>       Observation weights are parsed but not consumed by any active
#>       emitter, so a weighted fit is refused rather than silently
#>       returning the unweighted posterior.
#> 
#>   [since 0.9.3] weights_requires_gaussian
#>       `weights` is lowered for family = "gaussian" with the identity
#>       link only (brms's y | weights(w) addition term; INLA's scale = w
#>       per-observation precision multiplier -- both give Var(y_i) =
#>       sigma^2 / w_i on the Gaussian likelihood). Superseded the
#>       all-families weights_not_supported refusal for this one
#>       well-defined mapping (C6); every other family, and a non-identity
#>       link on gaussian, still refuse -- weights_not_supported remains
#>       registered for the historical all-families refusal text.
#> 
fb_refusals(reason_code = "precision_not_symmetric")
#> flexyBayes refusal registry: 1 reason  (filter: reason_code in {precision_not_symmetric})
#> 
#>   [since 0.4.0] precision_not_symmetric
#>       vm(..., precision = ): the precision matrix is not symmetric.
#> 
fb_refusals(since_version = "0.4")
#> flexyBayes refusal registry: 40 reasons  (filter: since_version ~ '0.4')
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
#>       `random` / `residual` arguments on the universal entry.
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
#>   [since 0.4.0] low_rank_scheme_required
#>       vm()/ped(): low_rank_factor supplied without an explicit
#>       low_rank_scheme naming a registered approximation.
#> 
#>   [since 0.4.0] low_rank_smooth_unsupported
#>       A smooth requesting the low_rank_smooth approximation was routed
#>       to a backend that cannot honour it. No active engine represents
#>       the rank-K basis truncation, so the request refuses under `auto`
#>       as well as under an explicit `inla` / `brms`.
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
#>   [since 0.4.0] representation_unknown_for_preflight
#>       Preflight: the design representation is not characterised by the
#>       preflight memory estimator.
#> 
#>   [since 0.4.0] residual_type_unsupported_for_aggregation
#>       Aggregated Gaussian emit: the residual term type is outside the
#>       supported aggregation scope.
#> 
#>   [since 0.4.0] response_not_in_data
#>       The response variable named on the formula's left-hand side is
#>       not a column of `data`.
#> 
#>   [since 0.4.0] review_code_backend_unsupported
#>       review_code = TRUE was requested on a backend with no code to
#>       review. The inspect-then-fit token is scoped to the code-emitting
#>       engines, which among the active pair is brms.
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
