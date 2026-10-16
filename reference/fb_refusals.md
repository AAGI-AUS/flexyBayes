# List flexyBayes refusal reasons

`fb_refusals()` exposes the locked refusal-reason registry as a
browsable table: the canonical vocabulary of conditions under which
flexyBayes declines to fit, route, or validate a model, each with a
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
refusal reason, with columns `reason_code`, `package`, `description`,
`since_version`, and `plan_field`. The print method renders it as a
compact checklist.

## Details

Two optional filters narrow the listing. `reason_code` selects rows by
exact reason-code match (a single code or a vector). The `since_version`
filter selects reasons introduced in a matching release by
version-string prefix: `since_version = "0.4"` returns every reason
added in the 0.4 series.

Routing-decision reasons (surfaced by
[`fb_plan()`](https://aagi-aus.github.io/flexyBayes/reference/fb_plan.md)
and
[`backend_decision()`](https://aagi-aus.github.io/flexyBayes/reference/backend_decision.md))
and internal control-flow signals are deliberately excluded: this table
lists only refusals a user can actually encounter.

Packages that extend flexyBayes add their own reasons with
[`fb_register_refusals()`](https://aagi-aus.github.io/flexyBayes/reference/fb_register_refusals.md).
The `package` column names the package each reason belongs to:
`"flexyBayes"` for the core, and the plug-in's name for the rest.

## Examples

``` r
fb_refusals()
#> flexyBayes refusal registry: 180 reasons
#> 
#>   [since 0.16.0] aggregate_argument_invalid
#>       `aggregate` is not TRUE, FALSE or "auto".
#> 
#>   [since 0.11.0.9000] aggregated_evidence_unavailable
#>       The marginal likelihood of a fit on the aggregated
#>       representation. The engine's value omits constants that depend on
#>       the data and, under response standardisation, scales by the
#>       number of cells instead of the number of rows, so it is not the
#>       marginal likelihood of the rows. summary() reports it as not
#>       available, naming this reason, instead of printing a number that
#>       cannot be compared.
#> 
#>   [since 0.16.0] aggregation_out_of_scope
#>       The aggregated Gaussian route or flexybayes_stream() was given a
#>       model outside the aggregation scope: a family or link it does not
#>       carry, a smooth, a random slope, a structured random term or a
#>       continuous fixed term. The gate's reason is carried in
#>       `scope_reason`.
#> 
#>   [since 0.16.0] aggregation_plan_ineligible
#>       `aggregate = TRUE` on a model the aggregation plan cannot
#>       compress. The plan's reasons are carried in `reason_codes`.
#> 
#>   [since 0.11.0.9000] aggregation_residual_prior_not_composed
#>       A Gaussian model sent to the aggregated representation, by
#>       aggregate = TRUE or by flexybayes_stream(). The aggregated
#>       Gaussian route composes the within-cell sum of squares with a
#>       fixed gamma prior on the residual precision, so a stated or
#>       default residual prior would not be the prior fitted while
#>       prior_summary() reported it. Gaussian models are fitted per row.
#>       Binomial and Poisson models still aggregate.
#> 
#>   [since 0.9.2] aggregation_response_incomplete
#>       aggregate = TRUE with missing values in the response. The
#>       aggregated path compresses rows into per-cell sufficient
#>       statistics, which are not defined for a missing response, so the
#>       request has no aggregated reading. Untyped before 0.9.2.
#> 
#>   [since 0.9.2] aggregation_route_unavailable
#>       aggregate = TRUE with no aggregated emit route available, either
#>       the named backend carries no aggregated path at all, or no
#>       aggregated backend resolves for the model, including a family or
#>       link flexybayes_stream() does not aggregate. Raised with
#>       site-specific messages. Untyped before 0.9.2.
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
#>       An AR1 process written on the residual side that flexyBayes does
#>       not fit: the 1-D and grouped forms, and the separable form for a
#>       response other than Gaussian with the identity link. The
#>       separable residual on a Gaussian response fits on INLA without a
#>       nugget.
#> 
#>   [since 0.16.0] ar1_residual_with_random_field
#>       The nugget-free autoregressive residual together with an
#>       autoregressive field on the random side. flexyBayes fits one of
#>       the two per model.
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
#>   [since 0.9.0] at_level_conditioning_unsupported
#>       at(f, level):units gives selected levels of f a residual variance
#>       of their own, and flexyBayes has no lowering for it on either
#>       engine. It is a different model from dsum(~ units | f), and is
#>       refused by name.
#> 
#>   [since 0.16.0] at_level_not_in_data
#>       at(f, level) names a level that is not a level of f in the data.
#>       Levels are literal values or positions in levels(f), as ASReml
#>       reads them.
#> 
#>   [since 0.9.0] augment_cell_not_determinable
#>       A design cell is absent from the data and the model depends on a
#>       variable whose value at that cell was never recorded.
#> 
#>   [since 0.11.0.9000] augment_not_available_inla
#>       augment() on an INLA fit, which keeps no per-row model frame for
#>       the fitted-value and residual columns. tidy() and predict() are
#>       the alternatives.
#> 
#>   [since 0.9.0] auto_no_active_route
#>       backend = "auto" found neither engine able to fit the model: the
#>       INLA route refused it and flexyBayes has no brms lowering for it.
#>       There is no silent fallback, so auto refuses.
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
#>       A missing response on brms for a family other than gaussian and
#>       lognormal, the two flexyBayes has checked with brms's mi() term.
#>       Without mi() brms drops the rows.
#> 
#>   [since 0.9.0] brms_cannot_represent_term
#>       A random-effect term type reached the brms formula reconstruction
#>       with no branch to lower it. It is refused by name instead of
#>       dropped from the emitted model.
#> 
#>   [since 0.9.2] brms_factor_random_slope_unsupported
#>       A bar-grammar random slope on a factor ((f | g), (0 + f | g), (f
#>       || g)) is outside the brms ingest. The ASReml surface of the same
#>       model fits, so the refusal names it verbatim. Before 0.9.2 the
#>       message listed other bar spellings, none of which expresses the
#>       model.
#> 
#>   [since 0.9.2] brms_ingest_feature_unsupported
#>       A brms formula feature outside the ingest layer's scope: an
#>       extended term the walker collected but cannot lower, or a part of
#>       a brms::bf() object the ingest does not carry to an engine (a
#>       non-linear declaration, a formula or fixed value for a
#>       distributional parameter, or a family that differs from the
#>       `family` argument). The extended terms were untyped before 0.9.2,
#>       and the brms::bf() parts were dropped with no message before
#>       0.14.0.
#> 
#>   [since 0.9.2] brms_random_effect_form_unsupported
#>       A bar-grammar random-effect specification outside the supported
#>       set, multi-variable slopes and group-of-coefficients shorthand.
#>       Untyped before 0.9.2.
#> 
#>   [since 0.10.0.9000] brms_stan_toolchain_compile_failure
#>       The brms path compiles a Stan program, and a broken C++ toolchain
#>       reaches R as a downstream symptom, `invalid connection`, not as
#>       the compiler's own message. On a failed fit the toolchain is
#>       probed directly, and this code is raised only when the probe
#>       confirms it cannot compile.
#> 
#>   [since 0.9.3] cell_count_exceeds_integer
#>       The aggregation planner's estimated cell count (the product of
#>       the cell key's factor level counts) exceeds R's integer limit.
#>       Recording it would coerce it to NA, so the plan's cell count and
#>       compression estimate would both silently disappear. The row-count
#>       sibling (row_count_exceeds_integer) fixed the identical defect
#>       for N, this is the cell-count member of the same class.
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
#>   [since 0.16.0] classify_argument_not_supported
#>       predict(classify = ) was given an argument of ASReml's predict()
#>       that changes the default prediction rules or asks for output the
#>       table does not carry.
#> 
#>   [since 0.9.1] classify_random_factor_not_supported
#>       predict(classify = ) was asked for a factor that enters the model
#>       as a random term the predicted-value table cannot use.
#> 
#>   [since 0.9.1] classify_requires_emmeans
#>       predict(classify = ) builds its marginal-means table through the
#>       emmeans estimation seam, and emmeans is a suggested package, not
#>       a hard dependency.
#> 
#>   [since 0.4.0] code_flags_mutually_exclusive
#>       return_code = TRUE and review_code = TRUE were both supplied. The
#>       two code-return modes are mutually exclusive.
#> 
#>   [since 0.9.0] conditional_loglik_not_available
#>       No plug-in conditional log-likelihood exists for this fit, the
#>       family has no implemented form, the response or family record is
#>       absent, or the engine reports a marginal likelihood instead.
#>       Returning NA would let AIC() and anova() present a comparison
#>       that was never computed.
#> 
#>   [since 0.10.0.9000] cor_slope_prior_not_translatable_inla
#>       (x | g) on INLA is an iid2d field whose hyperparameters are a
#>       single joint Wishart over the 2-by-2 covariance of the intercept
#>       and the slope. It has no standard-deviation slot, so a supplied
#>       SD prior cannot be translated. Refused instead of dropped.
#> 
#>   [since 0.9.0] corh_no_equicorrelation_representation
#>       corh(f):g, or cor(f), asks for one correlation shared by every
#>       pair of levels. flexyBayes has no lowering for it on either
#>       engine.
#> 
#>   [since 0.11.0.9000] correlated_covariance_requires_brms
#>       An estimated correlation matrix (a correlated random slope (x |
#>       g) or an unstructured us(f):g) was sent to INLA. INLA carries it
#>       as an iid2d or iidKd field under a joint Wishart prior that
#>       couples the variances and the correlations, and on that route the
#>       fitted covariance did not match REML's and was not invariant to
#>       rescaling the response. flexyBayes fits these structures on brms.
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
#>       memory ceiling. Dispatch is short-circuited before any backend
#>       code runs.
#> 
#>   [since 0.10.0.9000] dsum_malformed
#>       A dsum() body carrying no grouping bar. dsum(~ units | env)
#>       requests one residual variance per section. Written without the
#>       bar there is nothing to section on, and the parser used to fall
#>       through to a homogeneous `units` residual, a plausible neighbour
#>       that survives convergence checks and prints credible output.
#> 
#>   [since 0.9.0] dsum_structured_inner_unsupported
#>       Residual dsum() with a structured inner term (ar1, us, ...).
#>       flexyBayes represents dsum() only as a per-region heteroscedastic
#>       variance, so a structured inner is refused instead of silently
#>       reduced to a heteroscedastic residual.
#> 
#>   [since 0.5.0] engine_pin_backend_conflict
#>       An engine pin (fb_inla / fb_brms) was given a `backend` argument
#>       that conflicts with the engine it pins.
#> 
#>   [since 0.11.0.9000] evidence_not_computed
#>       The marginal likelihood of a fit whose engine computed none: a
#>       brms fit, where Stan sampling runs no bridge sampling. summary()
#>       and the model manifest report the evidence as not available under
#>       this reason.
#> 
#>   [since 0.9.0] fa_not_representable
#>       A factor-analytic covariance written without a grouping factor,
#>       fa(f, k). It parses for the formula catalogue, and flexyBayes
#>       does not fit a factor-analytic covariance in any form.
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
#>       A factor-analytic term fa(x, k) was given a rank that is not a
#>       positive whole number written in the formula. Registered at 0.4.0
#>       as the k < 1 floor and widened by the 0.10.0.9000 grammar
#>       contract to every non-count rank, a symbol, a string, a logical,
#>       a fraction, each of which as.integer() used to coerce into a rank
#>       the user did not write.
#> 
#>   [since 0.10.0.9000] fa_rank_missing
#>       fa() written without its rank. The rank is a number the user
#>       chooses, so it is named separately from the general arity
#>       refusal.
#> 
#>   [since 0.13.0] factor_analytic_not_repeatable
#>       A factor-analytic genotype-by-environment covariance, fa(f, k):g,
#>       sent to a fit. INLA represents it by entering each common factor
#>       as scaled copies of a latent field, but on simulated trials with
#>       a known covariance the genetic-correlation intervals of INLA fits
#>       covered the true value less often than required. flexyBayes has
#>       no factor-analytic lowering for brms, so neither engine fits the
#>       structure here. The diagonal structure diag(f):g fits on both
#>       engines, and us(f):g fits on brms.
#> 
#>   [since 0.9.0] family_argument_not_recognised
#>       The `family` argument was neither a single family name nor a
#>       usable `stats` family object, or a family object and an explicit
#>       `link =` named different links. A family object supplies its own
#>       link, so resolving the disagreement either way would fit a model
#>       the call did not ask for.
#> 
#>   [since 0.13.0] family_response_domain_violation
#>       The requested family's likelihood is defined over a response
#>       domain (gamma / lognormal strictly positive, hurdle_gamma
#>       non-negative, beta an open unit interval) that the observed
#>       response violates. Refused before any backend code runs. A
#>       missing (NA) response is exempt, it is na_action's concern, not a
#>       domain violation. Untyped before 0.13.0.
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
#>       A method needing the posterior itself (draws on a sampling
#>       engine, marginal densities on INLA) was called on a fit that
#>       carries neither. The refusal replaces an empty result that would
#>       read as an answer. Its sibling is loo_requires_sampler_draws,
#>       raised when the fit does carry a posterior but not the pointwise
#>       log-likelihood draws that PSIS-LOO reads.
#> 
#>   [since 0.14.1] fixed_effects_aliased
#>       flexybayes(): fixed-effect columns are linear combinations of the
#>       columns before them on the rows with a response, as in a
#>       factorial with empty cells, so the data give no estimate for them
#>       and an engine would report them at their prior. The refusal names
#>       the columns, fixed_alias = "drop" drops them as asreml does.
#> 
#>   [since 0.11.0.9000] fixed_offset_invalid
#>       An offset() term that does not evaluate against the data to a
#>       finite number on every row: a name that is not a column, a factor
#>       or character value, a length that is not the row count, or a
#>       missing or infinite value. An offset enters the linear predictor
#>       with a fixed coefficient of one, so an NA would be read as zero
#>       by one engine and drop the row on the other.
#> 
#>   [since 0.11.0.9000] fixed_offset_unsupported
#>       An offset() term in the fixed formula on a route that cannot
#>       carry it. The parser records every offset on the model
#>       description, and a route with no exact lowering for it refuses by
#>       name instead of fitting the model without the offset, which is
#>       what happened before offsets were recorded.
#> 
#>   [since 0.9.2] fixed_smoother_not_supported
#>       An mgcv-style smoother (s(), te(), ti(), t2(), gp(), sos())
#>       appears in the fixed part of the formula. flexyBayes spells a
#>       univariate smooth `random = ~ spl(x)` on the INLA backend. Before
#>       0.9.2 the term reached the engine and surfaced as `could not find
#>       function "s"`, on brms only after a completed sampling run.
#> 
#>   [since 0.4.0] formula_not_two_sided
#>       The model formula must be two-sided (response ~ predictors). A
#>       formula carrying no left-hand-side response was supplied.
#> 
#>   [since 0.10.0.9000] grammar_argument_duplicated
#>       A special form was given the same named argument twice. The
#>       contract validated the set of names and not their multiplicity,
#>       so a doubled name was either silently resolved to whichever copy
#>       the branch happened to index, or escaped as an untyped matcher
#>       error.
#> 
#>   [since 0.10.0.9000] grammar_argument_empty
#>       A special form written with a trailing or doubled comma, leaving
#>       R's empty symbol where a value belongs. Touching it raised
#>       "argument is missing, with no default" from wherever it was first
#>       read, a base error about flexyBayes internals instead of a
#>       statement about the formula.
#> 
#>   [since 0.10.0.9000] grammar_arity_too_few
#>       An ASReml special form was written with fewer positional
#>       arguments than it takes. Before the contract these parsed into a
#>       term whose variable was NA, which reaches the engine as a model
#>       nobody wrote.
#> 
#>   [since 0.10.0.9000] grammar_arity_too_many
#>       An ASReml special form was written with more positional arguments
#>       than it reads. The surplus used to be dropped in silence, so a
#>       structure written into the call never reached the fitted model.
#> 
#>   [since 0.4.1] grammar_brms_known_matrices_unsupported
#>       `known_matrices` was supplied with brms-grammar ingest via the
#>       universal entry, which has no known-matrix carrier.
#> 
#>   [since 0.4.1] grammar_brms_with_asreml_terms
#>       A brms-style bar-grouped formula was combined with ASReml
#>       `random` / `residual` arguments on the universal entry.
#> 
#>   [since 0.10.0.9000] grammar_named_argument_unknown
#>       A special form was given a named argument it does not answer to.
#>       An unrecognised name used to be ignored, which reads as accepted.
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
#>       residual is sectioned by one factor. Nesting several would
#>       require their interaction, which the user must state explicitly.
#> 
#>   [since 0.17.0] inla_best_of_invalid
#>       Option flexyBayes.inla_best_of, the number of INLA fits each
#>       model is run with, is not a whole number of one or more. The
#>       request is refused before anything is fitted.
#> 
#>   [since 0.16.0] inla_config_precision_not_positive_definite
#>       The latent-field precision INLA stored for a hyperparameter
#>       configuration is not positive definite, and the fit's linear
#>       constraints do not make it so: the posterior of the latent field
#>       is improper in some direction.
#> 
#>   [since 0.9.0] inla_gate_refused
#>       An explicit backend = "inla" request that lgm_gate() refused. The
#>       gate's structural check-list is the message and rides on the
#>       condition in its `lgm_refusal` slot. The failing rule supplies a
#>       second condition class so one rule can be caught on its own.
#> 
#>   [since 0.9.3] inla_program_failed
#>       The INLA engine died, a non-zero exit from the inla-program
#>       subprocess, an inla.core.safe failure, or an empty/NULL result
#>       with no R-level error at all, past the design-memory preflight,
#>       which models per-term storage but not the sparse-Cholesky solver
#>       cost of a large latent field. Before 0.9.3 this reached the user
#>       as a raw pass-through of the engine's own message after tens of
#>       minutes of runtime, with no reason code and no remedy. Carries
#>       the design size, the largest design this package has verified an
#>       INLA fit to complete (read from
#>       inst/validation/benchmark_scaling.csv), the binding random-effect
#>       term where determinable, and the remedies.
#> 
#>   [since 0.9.2] inla_variable_used_twice
#>       The same variable indexes both a fixed-effect term and a latent
#>       f() term in the emitted INLA formula. INLA refuses a duplicated
#>       key, before 0.9.2 the refusal surfaced as the generic
#>       inla-program-exited message.
#> 
#>   [since 0.9.0] interaction_not_representable
#>       An interaction of two structured terms that matches none of the
#>       recognised interaction patterns, us(trait):vm(gen), say. It
#>       parsed into the catalogue's generic-interaction node, which
#>       flexyBayes has no lowering for on either engine.
#> 
#>   [since 0.14.0] isolated_process_died
#>       A separate R process running an engine fit or a job of
#>       fb_isolated_lapply() ended before it returned: a signal ended it
#>       (SIGKILL is sent from outside the process, SIGSEGV or SIGABRT
#>       comes from a crash in compiled code), or it exited with a
#>       non-zero status. The refusal names the engine, the entry point,
#>       the fold where there is one, the exit status and the signal. The
#>       session that asked for the fit is unaffected, and nothing from
#>       the process is returned.
#> 
#>   [since 0.14.0] isolation_requires_callr
#>       options(flexyBayes.isolation = "process") runs each engine fit,
#>       fb_isolated_lapply() each job, and option flexyBayes.inla_best_of
#>       above 1 the further INLA fits, in a separate R process started
#>       through callr, which is in Suggests and not installed. The
#>       request is refused before anything is fitted, so no fit runs
#>       without the protection that was asked for.
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
#>   [since 0.14.0] known_variance_invalid
#>       `known_variance` is not one finite, non-negative variance for
#>       each row of `data`: it is not numeric, its length differs from
#>       the row count, or a value is missing, infinite or negative. Or it
#>       is a matrix that is not square, does not match the row count, is
#>       not symmetric, or is not positive definite.
#> 
#>   [since 0.14.0] known_variance_not_supported
#>       A known variance per row was asked of a model or a route that
#>       cannot carry it: a family other than the Gaussian on the identity
#>       link, the aggregated route, a known covariance matrix on brms, or
#>       a route whose refit, covariance or predictive score would leave
#>       it out. The message names the route that carries it.
#> 
#>   [since 0.11.0.9000] log_density_producer_unavailable
#>       fb_log_posterior() on a fit from either engine: flexyBayes has no
#>       log-density producer for them. INLA's posterior is a Laplace
#>       approximation, and brms's log-density is defined on the Stan
#>       unconstrained scale. The condition also carries the class
#>       fb_c4_unavailable.
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
#>       A smooth requesting the low_rank_smooth approximation. flexyBayes
#>       has no route for the rank-K basis truncation on either engine, so
#>       the request refuses under `auto` as well as under an explicit
#>       `inla` / `brms`.
#> 
#>   [since 0.16.0] marginal_cov_not_buildable
#>       fb_marginal_cov(): the model carries a random or residual
#>       structure the marginal covariance is not built for, its
#>       variance-component table lacks a variance the structure needs, or
#>       as.matrix() was asked of a structure with no variances. The
#>       refusal names the terms.
#> 
#>   [since 0.11.0.9000] marginaleffects_link_not_supported
#>       A marginaleffects prediction on a fit whose link is not the
#>       identity. Posterior draws from fb_as_draws_simple() are the
#>       alternative.
#> 
#>   [since 0.9.0] met_summary_not_available
#>       fb_met_summary() was called on a fit it cannot summarise: the
#>       breeder summary reads realised factor-analytic effects, and
#>       flexyBayes fits no fa() term on either engine, so this refusal is
#>       unconditional.
#> 
#>   [since 0.16.0] method_argument_not_supported
#>       A method was given an argument that the same method on another
#>       fit class takes and this one does not implement.
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
#>   [since 0.16.0] model_variable_not_in_data
#>       A fixed or random term names a variable that is not a column of
#>       the data.
#> 
#>   [since 0.10.0.9000] multivariate_smooth_unsupported
#>       s() written over more than one variable. flexyBayes fits
#>       univariate penalised splines only. This used to raise an untyped
#>       names-attribute error from mgcv's call handling.
#> 
#>   [since 0.10.0.9000] nested_index_name_collision
#>       Fitting an interaction random effect on INLA generates an index
#>       column, and the data already carries a column of that name.
#>       Overwriting it would index the field by whatever that column
#>       holds with nothing downstream saying so, which is why this
#>       refuses instead of renaming silently.
#> 
#>   [since 0.10.0.9000] nested_variable_not_in_data
#>       An interaction random effect names a variable that is not a
#>       column of the data, so the interaction index cannot be built.
#> 
#>   [since 0.11.0.9000] newdata_offset_missing
#>       predict() on new data from a model with an offset, where the new
#>       data lack a column the offset reads. The offset is part of every
#>       prediction and is never taken to be zero.
#> 
#>   [since 0.9.0] numeric_variable_in_random_interaction
#>       A numeric variable inside a nested or multi-way random
#>       interaction. The spelling reads as a random regression to an
#>       ASReml user and as one independent deviation per cell to the
#>       emitted model, flexyBayes refuses instead of choosing one of the
#>       two readings.
#> 
#>   [since 0.15.0.9001] plan_argument_unknown
#>       fb_plan() was given an argument it does not take, such as a
#>       misspelt name, which would otherwise be dropped and the plan made
#>       for a different model.
#> 
#>   [since 0.10.0.9000] pol_degree_invalid
#>       pol()'s degree is not a positive whole number written in the
#>       formula. Same coercion defect as fa()'s rank.
#> 
#>   [since 0.17.0] pooled_prior_requires_brms
#>       A pooled() prior, one log-normal prior with an estimated mean and
#>       SD shared by a set of per-level SDs, was sent to INLA. INLA
#>       cannot put a shared prior on separate variance parameters, so
#>       flexyBayes fits pooled() on brms only.
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
#>   [since 0.11.0.9000] predict_link_not_supported
#>       predict() without classify on an INLA fit whose link is not the
#>       identity. Posterior draws from fb_as_draws_simple() are the
#>       alternative.
#> 
#>   [since 0.16.0] predict_new_levels_family_not_supported
#>       predict(allow_new_levels = TRUE) on a fit whose response is not
#>       Gaussian on the identity link.
#> 
#>   [since 0.16.0] predict_new_levels_no_variance
#>       predict(allow_new_levels = TRUE) at a new level the fit has no
#>       variance for: a new level of f in diag(f):g, a new section of
#>       dsum(~ units | f), or a row or column off the lattice of an AR1 x
#>       AR1 field.
#> 
#>   [since 0.16.0] predict_new_levels_structure_not_supported
#>       predict(allow_new_levels = TRUE) on a model carrying a random or
#>       residual structure, observation weights or a known variance for
#>       which the draw at a new level is not built. The refusal names
#>       them.
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
#>   [since 0.11.0.9000] prior_family_has_no_residual_scale
#>       A residual prior, stated or the package default, on a family
#>       whose likelihood has no residual scale parameter: Poisson and
#>       binomial, and on brms every family whose dispersion is a function
#>       of the mean. No engine slot reads it, and the model manifest
#>       records the input as not consumed under this reason.
#> 
#>   [since 0.9.2] prior_hyperparameter_not_scalar
#>       A prior hyperparameter is not a finite numeric scalar. A
#>       zero-length or length-two value used to surface as one of R's own
#>       interpreter messages several layers downstream.
#> 
#>   [since 0.9.2] prior_hyperparameter_out_of_domain
#>       A prior hyperparameter lies outside the domain its family
#>       defines: a non-positive scale, rate, shape or degrees of freedom,
#>       a tail probability outside (0, 1), or reversed uniform bounds. A
#>       negative standard deviation used to reach the sampler and surface
#>       as a missing-draws error.
#> 
#>   [since 0.16.0] prior_not_fb_prior
#>       `prior` is neither an fb_prior() object nor NULL.
#> 
#>   [since 0.9.2] prior_not_translatable_for_backend
#>       flexyBayes does not translate a supplied prior specification for
#>       the selected backend. Before 0.9.2 the INLA route discarded such
#>       a row silently while prior_summary() printed it as applied, and
#>       the brms route refused it with an untyped error.
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
#>   [since 0.12.1] prior_support_violated
#>       An INLA fit whose mode, or whose reported 2.5 or 97.5 per cent
#>       quantile, for a hyperparameter under a bounded prior lies outside
#>       that bound. The uniform prior on a standard deviation reaches
#>       INLA as a density of -1e10 outside its support, and INLA can
#>       still stop there. Every INLA fit is checked against INLA's own
#>       record of the priors it received, and a fit outside its prior's
#>       support is refused.
#> 
#>   [since 0.9.2] prior_target_argument_missing
#>       A prior target call is missing the argument that identifies which
#>       parameter it targets, sd() and cor() need `group`, b() and
#>       smooth() need a name.
#> 
#>   [since 0.10.0.9000] prior_target_not_a_single_name
#>       A b() prior target was given something other than one coefficient
#>       name. as.character() on the unevaluated argument deparses a call
#>       head-first, so b(c("x", "z")) became the three names c("c", "x",
#>       "z") and survived fb_plan() before failing inside the emitter
#>       with an untyped base error.
#> 
#>   [since 0.9.2] prior_target_not_in_model
#>       A prior names a coefficient or a variance component the model
#>       does not have. The check reads the engine's own parameter list,
#>       so the refusal names the model's actual terms instead of a
#>       synthesised Stan parameter name.
#> 
#>   [since 0.9.2] prior_target_unsupported
#>       The left-hand side of a prior formula is not one of the supported
#>       targets. A bare name that is not `sigma` used to be carried as an
#>       opaque target and dropped at emit time.
#> 
#>   [since 0.11.0.9000] prior_vc_sd_superseded_by_prior
#>       `prior_vc_sd` supplied beside an fb_prior(). No route reads the
#>       scalar once an fb_prior() is present: INLA gives a standard
#>       deviation the prior leaves unstated the package's bounded
#>       default, and brms keeps its own. The model manifest records the
#>       scalar as not consumed under this reason.
#> 
#>   [since 0.15.0.9000] refit_hooks_require_brms
#>       fb_refit(): `extra_prior`, `stanvars` and `formula_hook` are brms
#>       objects and reach only a brms fit. The fit was not fitted with
#>       brms and no `backend = "brms"` was given.
#> 
#>   [since 0.15.0.9000] refit_not_reconstructable
#>       fb_refit(): the refit cannot be built from what the fit records.
#>       The fit keeps no fixed formula or no data, or it carries
#>       observation weights and new data arrived without `weights`, so
#>       the recorded weights cannot be matched to the new rows.
#> 
#>   [since 0.4.0] representation_unknown_for_preflight
#>       Preflight: the design representation is not characterised by the
#>       preflight memory estimator.
#> 
#>   [since 0.4.0] residual_type_unsupported_for_aggregation
#>       Aggregated Gaussian emit: the residual term type is outside the
#>       supported aggregation scope.
#> 
#>   [since 0.11.0.9000] response_not_finite
#>       A response with NaN or infinite values, usually from a
#>       transformation such as log(0) or from a data-entry error. They
#>       are not missing observations, so na_action does not apply to
#>       them, and an engine handed them either fitted the model without
#>       saying so or failed without naming the row. The fit refuses
#>       before any engine is called.
#> 
#>   [since 0.4.0] response_not_in_data
#>       The response variable named on the formula's left-hand side is
#>       not a column of `data`.
#> 
#>   [since 0.16.0] response_not_numeric
#>       The response is a factor or character vector, which no family
#>       flexyBayes fits takes.
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
#>   [since 0.14.0] spatial_field_requires_ar1_field
#>       plot(type = "spatial_field") draws the posterior mean of the
#>       separable autoregressive field random = ~ ar1(row):ar1(col), and
#>       the fit's random terms name no such field.
#> 
#>   [since 0.9.0] stan_cannot_represent_ar1_field
#>       A random-side AR1 or separable AR1xAR1 latent field on brms.
#>       flexyBayes builds the field on INLA only and has no brms lowering
#>       for it.
#> 
#>   [since 0.9.0] stan_cannot_represent_ar1_residual
#>       flexyBayes has no brms lowering for an AR1 or separable AR1xAR1
#>       residual covariance.
#> 
#>   [since 0.5.0] stan_cannot_represent_structured_cov
#>       flexyBayes has no brms lowering for a structured-covariance term
#>       in the model: a vm() / ped() term with a blocks or low-rank
#>       carrier, or a structured term written without a grouping factor.
#> 
#>   [since 0.9.0] stan_cannot_represent_structured_residual
#>       flexyBayes has no brms lowering for a residual structure other
#>       than the iid `units` residual and a residual per section.
#> 
#>   [since 0.9.0] str_not_representable
#>       str() declares a covariance structure over a group of terms. It
#>       parses for the formula catalogue and has no emit path on either
#>       active engine.
#> 
#>   [since 0.13.0] stream_model_column_not_in_source
#>       flexybayes_stream() reads named columns off the source chunk by
#>       chunk, so every column the model names (the response, a fixed or
#>       random term, `trials`, or `exposure`) has to be present in the
#>       source's head sample. Untyped before 0.13.0.
#> 
#>   [since 0.10.0.9000] stream_response_must_be_a_column
#>       Streaming aggregation reads named columns from the source chunk
#>       by chunk, so a transformed response, log(yield), has no column to
#>       read. The per-row route evaluates the expression against the
#>       data. The streaming route cannot, and says so instead of failing
#>       inside the reader.
#> 
#>   [since 0.16.0] stream_source_invalid
#>       flexybayes_stream() was given a `source` that is not an existing
#>       file, a vector of .fst shards, a data frame or a chunk generator.
#> 
#>   [since 0.13.0] stream_source_no_rows
#>       flexybayes_stream() found no rows to work with: either the head
#>       sample used to build the model, or the full read across every
#>       chunk of the source, returned zero rows. There is nothing to
#>       aggregate. Untyped before 0.13.0.
#> 
#>   [since 0.11.0.9000] streamed_rows_not_retained
#>       A row-level accessor (fitted(), residuals(), model.matrix(),
#>       insight::get_data(), emmeans::recover_data()) on a
#>       flexybayes_stream() fit that keeps fewer rows than it was fitted
#>       to. The stream keeps only the first rows it read, so values
#>       computed from them are not the values of the data.
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
#>   [since 0.9.3] unknown_backend
#>       An explicit `backend` request named a value that is not a
#>       registered engine, including a name this release withdrew
#>       entirely (see NEWS.md). The active engines are inla and brms.
#> 
#>   [since 0.4.0] unsupported_family
#>       The requested family, or the link requested with it, is outside
#>       the set flexyBayes can emit. Refused at the family gate
#>       (.resolve_family) before any backend code runs.
#> 
#>   [since 0.9.0] update_call_not_reconstructable
#>       update() was called on a fit whose recorded argument set is
#>       incomplete. Re-fitting would substitute package defaults for the
#>       unrecorded arguments, a relationship matrix or a prior scale
#>       among them, and return a different model under the same name.
#> 
#>   [since 0.11.0.9000] update_streamed_fit_not_supported
#>       update() on a flexybayes_stream() fit. The fit keeps neither the
#>       source it read nor every row, so a re-fit through flexybayes()
#>       would run on different data.
#> 
#>   [since 0.9.3] update_unnamed_argument_not_supported
#>       update() re-issues every recorded flexybayes() argument by name,
#>       matching an override in ... to it by name too. An unnamed
#>       argument, stats::update()'s classic replacement-formula idiom,
#>       update(fit, . ~ . + z), has no name to match and was, before this
#>       refusal, silently discarded: the re-fit ran with the unchanged
#>       recorded call, a valid object answering a question nobody asked
#>       instead of a re-fit or a refusal. Name the slot instead:
#>       update(fit, fixed = y ~ x + z) (brms-style), or random = ~ ... /
#>       residual = ~ ... on the ASReml grammar.
#> 
#>   [since 0.10.0.9000] us_gxe_prior_not_translatable_inla
#>       us(f):g on INLA is an iidKd field carrying a joint Wishart prior
#>       over the whole covariance. It has no standard-deviation slot, so
#>       a supplied SD prior cannot be translated. Refused instead of
#>       dropped.
#> 
#>   [since 0.9.1] variogram_requires_design_index
#>       plot(type = "variogram") needs the row / column array a
#>       structured covariance is indexed by, and the fit's terms name
#>       none.
#> 
#>   [since 0.4.0] vm_redundant_specification
#>       vm()/ped(): more than one covariance carrier supplied. Exactly
#>       one of V / chol / precision / blocks / low_rank_factor is
#>       allowed.
#> 
#>   [since 0.10.0.9000] vpredict_component_not_in_posterior
#>       A functional names a variance component the fit does not carry.
#>       Resolving the name against the caller's workspace instead would
#>       return a number computed from something other than this
#>       posterior, so it refuses and lists the components that are
#>       available.
#> 
#>   [since 0.10.0.9000] vpredict_functional_not_one_value_per_draw
#>       A functional returned something other than one number per
#>       posterior draw, so it is not a transform of the posterior and has
#>       no posterior of its own to summarise.
#> 
#>   [since 0.10.0.9000] vpredict_no_functional
#>       vpredict() was called with no functional, or with an argument
#>       that is not a two-sided `name ~ expression` formula.
#> 
#>   [since 0.16.0] weights_invalid
#>       `weights` is not a numeric vector with one value per row of the
#>       data.
#> 
#>   [since 0.9.3] weights_not_aggregatable
#>       aggregate = TRUE with non-constant observation weights. Weights
#>       are lowered for the per-row Gaussian route only. The aggregated
#>       route's sufficient statistics (per-cell sums and counts) do not
#>       fold in a per-observation weight, so a weighted aggregated fit is
#>       refused instead of silently ignoring the weight or answering a
#>       different question under the same name.
#> 
#>   [since 0.9.0] weights_not_supported
#>       Observation weights are parsed but not consumed by any active
#>       emitter, so a weighted fit is refused instead of silently
#>       returning the unweighted posterior.
#> 
#>   [since 0.9.3] weights_requires_gaussian
#>       `weights` is lowered for family = "gaussian" with the identity
#>       link only (brms's y | weights(w) addition term; INLA's scale = w
#>       per-observation precision multiplier, both give Var(y_i) =
#>       sigma^2 / w_i on the Gaussian likelihood). Superseded the
#>       all-families weights_not_supported refusal for this one
#>       well-defined mapping. Every other family, and a non-identity link
#>       on gaussian, still refuse, weights_not_supported remains
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
#>       return_code = TRUE and review_code = TRUE were both supplied. The
#>       two code-return modes are mutually exclusive.
#> 
#>   [since 0.4.0] cov_arg_not_fb_cov
#>       vm() / ped(): the `cov` argument must be written inline as an
#>       fb_cov() carrier.
#> 
#>   [since 0.4.0] design_memory_exceeds_ceiling
#>       Preflight: the design matrix is estimated to exceed the active
#>       memory ceiling. Dispatch is short-circuited before any backend
#>       code runs.
#> 
#>   [since 0.4.0] fa_rank_invalid
#>       A factor-analytic term fa(x, k) was given a rank that is not a
#>       positive whole number written in the formula. Registered at 0.4.0
#>       as the k < 1 floor and widened by the 0.10.0.9000 grammar
#>       contract to every non-count rank, a symbol, a string, a logical,
#>       a fraction, each of which as.integer() used to coerce into a rank
#>       the user did not write.
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
#>       The model formula must be two-sided (response ~ predictors). A
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
#>       A smooth requesting the low_rank_smooth approximation. flexyBayes
#>       has no route for the rank-K basis truncation on either engine, so
#>       the request refuses under `auto` as well as under an explicit
#>       `inla` / `brms`.
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
#>       The requested family, or the link requested with it, is outside
#>       the set flexyBayes can emit. Refused at the family gate
#>       (.resolve_family) before any backend code runs.
#> 
#>   [since 0.4.0] vm_redundant_specification
#>       vm()/ped(): more than one covariance carrier supplied. Exactly
#>       one of V / chol / precision / blocks / low_rank_factor is
#>       allowed.
#> 
```
