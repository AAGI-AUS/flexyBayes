# flexyBayes 0.10.0

The first release prepared for CRAN. Its criteria are held as a runnable
contract in the development repository rather than as prose. The checks
that read source answer immediately; the ones that need a build, a check
run or a person are listed as unmet until their artefact exists, and the
contract exits non-zero while any of them is outstanding. That file is
not shipped in the tarball.

## Breaking changes

* **The declared validation tier is now V2, not V3.** 0.9.2 declared V3
  and waived four of its floors (F20, F21, F24, F25); 0.9.3 shipped
  declaring V3 with that waiver already expired, and it is not re-issued.
  V2 is what the registered evidence supports. The four floors are
  recorded in `inst/validation/README.md` as the next validation arc, and
  a second consecutive waiver on the same floors is the pattern a waiver
  exists to prevent.

* **Eleven functions are withdrawn from the public API**, taking the
  export count from 44 to 33. The extreme-value and Dirichlet
  maximum-likelihood fitters (`fb_gev()`, `fb_family_gev()`, `rgev()`,
  `fb_dirichlet()`, `fb_family_dirichlet()`, `rdirichlet()`), the EMMAX
  genome-wide association pair (`fb_gwas()`, `triangulate_gwas()`), and
  the three exports that abstained unconditionally (`fb_met_summary()`,
  `fb_log_posterior()`, `fb_structured_cov()`) are no longer exported.
  They remain in the package, documented as internal and covered by their
  existing tests, and their S3 print and tidy methods stay registered so
  the objects still print.

  The reason is identity rather than arithmetic. This package's premise is
  a Bayesian posterior where REML returns a point estimate; four of the
  withdrawn fitters are maximum-likelihood or frequentist, none of them
  has a cell in the execution grid, and together they were a quarter of
  the public surface sitting outside the package's own execution oracle.
  `fb_structured_cov()` joined them late: no active engine emits an
  `fa()` term, so its only reachable outcome is a message and an empty
  list, which is the same ground the other two abstainers were withdrawn
  on. `fb_gblup_cv()`, `genomic_summary()` and `triangulate_genomic()`
  were considered for the same treatment and kept: the first is taught
  and executed in the multi-environment-trials and genomics vignette, and
  the other two operate on posterior draws.

* **`effectsize` is no longer listed as a supported downstream package.**
  It was named in the `DESCRIPTION` Description, in `Suggests`, in
  `README.md`, and in two vignettes -- one of which also said a
  `flexybayes`-class method dispatched through
  `parameters::model_parameters()`, which was never true. Removed from
  all five rather than left as a promise. The standardised-effect section
  of the downstream-analysis vignette keeps its worked example, which
  needs nothing but the draws.

* **`coda` moves from `Imports` to `Suggests`, and `splines` leaves
  `Imports` altogether.** The only trace of `coda` in `R/` was a roxygen
  `@importFrom` tag, and neither `effectiveSize()` nor `gelman.diag()`
  was called there; the test suite does use it, so it is demoted rather
  than dropped. `splines` was the same defect one release later:
  `importFrom(splines, bs)` with `bs()` called nowhere in `R/`. It moves
  to `Suggests`, where the one test that calls `splines::bs()` needs it.

## Diagnostics

* **A variance component pinned at the boundary now warns, on any term
  type on the per-row emit paths.** The existing detector covered one
  case, `sd_spatial` running to its floor against the nugget. Any other
  component could reach the same state reporting only a
  `note = "collapsed"` cell in `summary(fit)$varcomp`, which is easy to
  read past, while the convergence block reported a converged mode. Both
  engines now raise a warning naming the component, its upper credible
  bound and the residual scale, with three routes out. Silence with
  `options(flexyBayes.silence_boundary_collapse_warning = TRUE)`.

  The threshold is 0.005 of the residual SD, and it is calibrated against
  a measured floor. A sweep of 112 INLA fits with the group SD set to
  exactly zero, crossing n in {30, 60, 120, 240, 480} with 5, 10 and 20
  groups, puts the upper credible bound of a genuinely null component at
  0.0228 of the residual SD at the lowest, and the ratio is flat in both
  n and the group count because it is a floor set by the prior. The
  degenerate mode sits two orders of magnitude below that: 2.5e-04 on the
  `besag.met` fit this was built from. The threshold sits in the gap.

  Two limits on its reach, both structural. It compares against a
  residual SD, so it reads Gaussian-scale fits and is silent on families
  carrying no `sigma` row. And it needs credible bounds, so it covers the
  per-row emit paths only: the aggregated streaming emitters record
  variance components as posterior means with no quantiles, and there is
  nothing there to read an upper bound from.

* **`na_action` documentation now separates the identity from the
  arithmetic.** Under ignorability the posterior is the same whether a
  missing cell is augmented or omitted. That is a statement about the
  posterior and not a promise about what an optimiser returns: the two
  settings hand the engine an intact design and a ragged one, and on a
  weakly identified mode they can land in different places. `augment`
  remains the default and the recommendation. The caveat now stands
  wherever the identity is stated, `README.md` and the brms emit path
  included, rather than in one file.

* **Every export carries a runnable `\examples{}` block.** Nine had none.
  Five of the new blocks are wrapped in `\donttest{}` and guarded on the
  engine, so they do not fire where it is absent. The three tidier
  re-exports (`augment`, `glance`, `tidy`) document their generics on
  `man/reexports.Rd`, which carries no example; their methods do.

* **`Config/testthat/parallel` is now `false`.** It was `true`, and the
  parallel runner stalls: workers spawn, orphan, and sit at 0 per cent CPU
  indefinitely once the suite reaches the brms and INLA tests. Every
  release run so far has worked around it by forcing sequential execution
  from the outside, and `TESTTHAT_PARALLEL=false` does not reliably
  prevent the workers from spawning. A declared setting that hangs the
  suite is a defect for anyone who runs `devtools::test()`, so it is
  fixed rather than worked around.

## Corrections to the package's own claims

Each of these was a statement the code did not support. They are recorded
as fixes rather than dropped silently, and each now has a guard that fails
the lint workflow if it recurs.

* The package-level help said **twelve** vignettes ship, including one on
  arriving from an ASReml call. Eleven ship, and that page was folded into
  the getting-started vignette at 0.9.3.
* `API_STABILITY.md` said every public export carries a
  `lifecycle::badge("experimental")`. None did. The document now states
  the lifecycle stage without claiming a badge that is not there.
* `README.md` described the release as "a stable release of the supported
  capability set" six lines above a block stating that every export is
  experimental.
* `CITATION.cff` carried a `date-released` eight days earlier than the
  tree it described.
* An internal flag word reached a shipped surface in
  `inst/validation/README.md`.
* `cran-comments.md`, `SECURITY.md`, `SUPPORT.md` and `.zenodo.json` all
  still named 0.9.3. The coherence test that exists to catch exactly this
  could not: its "superseded version" pattern was written as `0\.[0-8]\.`
  and so was blind to the 0.9.x line at the moment the tree reached
  0.10.0. It now compares parsed versions and needs no edit at the next
  bump.
* Withdrawing eleven functions left their `\examples{}` blocks calling
  them by bare name, and `R CMD check` runs the examples of internal help
  pages too. Ten pages were affected. They now reach the functions through
  the namespace, and a new guard fails when any example calls something
  the package does not export.
* The boundary-collapse threshold shipped in this release at 0.05 of the
  residual SD on the strength of a single simulation, with a claim that a
  null component returns an upper bound near 0.43 of the residual SD. The
  package's own suite falsified it: seven fixtures whose components are
  null by construction warned. The claim is withdrawn and the threshold
  re-grounded on the sweep described above.
* The release contract reported its machine-evidence criteria as skipped
  and exited zero regardless, so the suite, both check runs, the grid,
  lint, pkgdown, the clean room and CI could each be red while it
  reported success. Those criteria are now unmet until evidenced, and the
  full tier exits non-zero while one is outstanding. The criterion
  asserting that a boundary-pinned component warns was satisfied by the
  existence of its test file; it runs the tests now.
* The guard asserting that every declared import is used counted mentions
  in comments and the package name inside message strings, which is how
  `splines` stayed in `Imports` unused. It now requires an imported
  symbol to be called.

# flexyBayes 0.9.3

This is the first public release of the 0.9 line. Versions 0.9.0 to 0.9.2
were staged locally and never published; their sections below stand as
the record of what changed between 0.8.3 and this release.

## Breaking changes

* **The third native fitting engine, quarantined since 0.9.0, is withdrawn
  from the package entirely.** No code path, export, S3 method, backend
  registry row, refusal code, worker script or `Suggests` entry for it
  remains. `fb_greta()`, `fb_from_greta()` and `gretaR_status()` are
  removed, and a call to any of them now fails with "could not find
  function" rather than a typed refusal. 0.9.2 was the last version whose
  code could read draws from an object fitted by that engine. Re-entry,
  should it ever be proposed, is a fresh implementation, not a repair of
  retained code. `greta` and `tensorflow` leave `Suggests`.
* **`backend = "greta"`, `backend = "gretaR"`, or any other unrecognised
  backend name now raises an ordinary `unknown_backend` refusal** naming
  the two active engines, `"inla"` and `"brms"`, before `match.arg()`
  could report a less specific error.
* **`fb_backend_status()` reports the two active engines only** (INLA,
  brms). The rows for the withdrawn engine and its dormant sibling are
  gone, with the Python/TensorFlow discovery probe that populated them.
* **The refusal-code registry drops eleven codes and adds nine** (116 in
  0.9.2, 114 now). Removed: `backend_quarantined`, the six codes of the dormant
  sibling engine's refusal family, the two codes naming the withdrawn
  engine directly, `low_rank_requires_greta` (renamed, below) and the
  stale `predict_kernel_invalid_include`, which had no live raise site.
  Added: `unknown_backend`, `low_rank_smooth_unsupported`,
  `cell_count_exceeds_integer`, `inla_program_failed`,
  `weights_requires_gaussian`, `weights_not_aggregatable`,
  `binomial_response_not_binary`, `update_unnamed_argument_not_supported`
  and `at_field_per_level_hyper_not_representable`; `fb_refusals()` is
  the authority.
* **`fb_met_summary()` abstains unconditionally** (`met_summary_not_
  available`). Its breeder summaries were derived from the withdrawn
  engine's realised factor-analytic effects, and neither active engine
  emits an `fa(env, k):gen` term. `summary()` and, on a `diag()` or
  `us()` covariance, `brms::VarCorr()` remain the way to read a
  multi-environment trial's genotype-by-environment structure.
* **`fb_log_posterior()` abstains unconditionally** (`fb_c4_unavailable`):
  neither active engine exposes an equivalent producer.
* **`update()` refuses the unnamed-formula idiom by name**
  (`update_unnamed_argument_not_supported`). `update(fit, . ~ . + z)`
  previously discarded the formula silently and re-fitted the unchanged
  model; named arguments (`fixed =`, `random =`, `residual =`) are
  unaffected.
* **GBLUP and pedigree-BLUP cross-engine comparison is two-engine**
  (INLA, brms); `triangulate_genomic()` and `genomic_summary()` keep
  their shape.
* The unreachable streaming `predict(..., output_file = )` path is
  removed; the engine-specific `predict()` methods had shadowed it since
  before this release.

## New features

* **A per-trial separable spatial field.** `random = ~
  at(trial):ar1(row):ar1(col)` fits one AR1 x AR1 field per level of
  `trial` on INLA, sharing the row correlation, column correlation and
  field variance across levels (INLA `replicate`), with the
  complete-lattice check run per level. brms has no lowering for it, as
  for the single-field spelling; a level-conditioned `at()` on this
  spelling refuses by name (`at_field_per_level_hyper_not_representable`).
* **`weights =` is lowered for the Gaussian family** (identity link) in
  the precision sense `Var(y_i) = sigma^2 / w_i` on both engines: INLA's
  per-observation `scale`, and on brms a known offset on the sigma
  distributional parameter that reproduces the same precision-weighted
  likelihood (brms's own `weights()` addition term is a likelihood-power
  quantity and was not used). Both engines match `lme4::lmer(weights = )`
  closely; doubling every weight leaves the fixed effects unchanged and
  doubles the fitted residual variance. Any other family, a non-identity
  Gaussian link, or `aggregate = TRUE` with weights refuses by name
  (`weights_requires_gaussian`, `weights_not_aggregatable`).
* **`fb_plan()` accepts the ASReml `fixed` / `random` / `residual`
  grammar**, auto-detected as `flexybayes()` does. Before, `random` and
  `residual` were silently dropped, so the plan could name a backend for a
  fixed-effects-only reading of the model (FS-22). `print.fb_plan()` and
  `print.fb_aggregation_plan()` print the row count, cell count, rows per
  cell and a plain verdict on whether aggregation will pay, with the
  threshold stated.
* **A non-syntactic factor level** (a level containing a space) no longer
  kills the INLA emit untyped: levels are legalised inside the emit and
  the user's own labels are printed back by `summary()`, `coef()`,
  `ranef()`, `predict(classify = )` and `confint()`; `backend = "auto"`
  now routes such models to INLA where it fell through to brms (FS-26).
  Purely numeric levels are left alone.
* **`glance()` returns its one-row summary on INLA and aggregated fits**
  as it did on brms fits, with sampler-specific columns `NA`.
* **A `family = "binomial"` response outside {0, 1}** refuses before any
  engine runs (`binomial_response_not_binary`), naming the column, the
  offending values and the remedy, on all three routes.
* `fb_plan()`, `backend_decision()` and `validate_approximation()` carry
  runnable examples.

## Bug fixes

* **An INLA engine death surfaces as a typed refusal**
  (`inla_program_failed`) naming the design size, the largest design this
  package has verified an INLA per-row fit to complete, the binding
  random-effect term and the remedies, instead of a raw engine message
  after tens of minutes (FS-25).
* **The aggregation planner's cell count is carried as a double.** Past
  R's integer limit the plan records `cell_count_exceeds_integer` instead
  of an `NA` count (FS-24): `aggregate = "auto"` takes the per-row route
  and says why, `aggregate = TRUE` refuses by name. A grep gate keeps the
  cast class closed.
* A streamed fit's `print()` banner names its engine and route instead of
  `(unknown engine)`.
* Gate labels spell `GxE` in ASCII; a gate message that named a
  pre-renumber vignette filename is corrected, and a test now requires
  every vignette filename named in `R/` to exist.
* `triangulate()`'s `transform_b` is exercised at a non-`NULL` value.
* The `low_rank_smooth` refusal is renamed `low_rank_smooth_unsupported`
  and says what is true: no active engine consumes a rank-truncated
  smooth basis.

## Documentation

* **The tutorial deck is eleven contiguous pages.** The page *From an
  ASReml call* is removed; its accessor tour (variance components, random
  effects, predicted means, observation counts, missing plots, the
  empirical variogram) lives in *Getting started*, section 6, with no
  `asreml()` call anywhere in the deck. Pages 06 to 11 and 16 are
  renumbered 05 to 11; the pkgdown site redirects every old address.
* **Every figure carries a caption and a reading sentence.** `html_
  vignette` leaves captions off by default, so a `fig.cap` reached the
  page only as alt text; every page now sets `fig_caption: true`.
* **Every page is copy-paste runnable**: `library(flexyBayes)` is the
  first visible line, seeds are visible, and no code depends on a hidden
  chunk. *Streaming exact aggregation* runs end to end on temporary
  `.fst` files (single file, shards, `fit = FALSE` inspection, binomial
  with `trials =`, poisson with `exposure =`), states how the Gaussian
  case stays exact as the emit implements it, and cites its methods.
* *Getting started* lists the four verbs one per row. *The formula
  surface* shows `weights =` fitting on Gaussian and refusing elsewhere;
  *Spatio-temporal models* shows the per-trial field. The kinship comment
  in *Multi-environment trials* attributes the marker-count normalisation
  to Astle and Balding (2009), which is what the code does.
* `@param random` separates terms that parse and fit from terms that
  parse and refuse by name (`fa()` on both engines). README carries a
  short paragraph on the accessor surface downstream tools read.

## Validation record

* `inst/validation/benchmark_scaling.{md,csv}` carry the 2026-08-22
  ceilings study on a realistic multi-term MET design. The flexyBayes/INLA
  ceiling on that design family is bracketed between 911,808 rows (the
  preflight refuses) and 1,823,616 rows (the engine dies after 41.6
  minutes) and has not been measured or bisected. A previously recorded
  340.9-second success at 911,808 rows was found unsupported by any run
  artefact and was corrected rather than carried forward.
* The execution grid gains cells for the space-level factor, Gaussian
  weights on both engines, and the per-trial field.
* `cran-comments.md` names the check artefacts this release produces and
  states the dormant log-posterior contract's oracle gap.

# flexyBayes 0.9.2

## Bug fixes

* **The prior mini-language refuses malformed specifications instead of
  substituting a default.** `fb_prior()` matched a distribution call
  against the family's parameter list and, when the match failed --
  which is exactly what a wrong argument name produces -- carried the
  call on unnamed, after which every emit-side read fell back to its own
  default. `half_normal(sd = 1.5)` fitted under `half_normal(scale =
  1)`, `pc(upper = 1)` fitted under a tail probability the caller never
  wrote, and a negative standard deviation reached the sampler and
  surfaced four layers downstream as a missing-draws error. Argument
  matching is now strict: an unknown or duplicated argument name, a
  missing required argument, a non-scalar or non-finite value, and a
  value outside its domain are each refused at construction with a
  `flexybayes_refusal_*` class, and both halves of a PC prior are
  required. The whole constructor surface is typed, so a caller can
  match a refusal by class rather than by message text.

* **Gamma and beta fit on the Stan backend again, on every prior route.**
  Which families carry a residual `sigma` was written down twice, and the
  two lists disagreed: the prior emit believed gamma and beta had one,
  the heteroscedastic-residual gate knew they did not. The prior emit
  won, so every route -- the plain default with no prior argument
  included -- sent brms a `sigma` prior for a parameter the model does
  not have and brms refused the fit. Two of the six advertised response
  families were unreachable on that backend. The fact now has one table
  (`R/family_traits.R`), both call sites read it, a test grounds it
  against brms's own declaration of each family's distributional
  parameters, and a second test fails if a duplicate roster reappears.

* **`prior_fixed_sd` and `prior_vc_sd` reach the engine, and
  `prior_summary()` tells one story about them.** Both arguments were
  documented as applied and neither was. `prior_fixed_sd` was absent
  from the condition that decides whether the auto-default prior fires,
  so passing it alone left the default in charge and the legacy bridge
  -- the only consumer of the scalar -- was never reached: the prior
  table was byte-identical to the same call without the argument.
  `prior_vc_sd` was consumed by the brms route and by nothing on the
  INLA route, so an INLA fit ran under INLA's own log-gamma precision
  default while `prior_summary()` printed `Lognormal(0, 3)` in its
  header and, four lines further down, said no prior had been supplied
  and the engine had used its own. Supplying either argument now applies
  it -- `prior_fixed_sd` through brms prior rows and INLA's
  `control.fixed`, `prior_vc_sd` through a brms `lognormal` row and the
  INLA expression prior that writes the same density on the
  standard-deviation scale -- on the per-row and aggregated routes
  alike, so the two engines carry the same variance-component prior and
  `triangulate()` can compare the components. Leaving either unsupplied
  keeps each engine's own default, which `prior_summary()` now names
  instead of asserting a prior that is not there. A residual prior
  declared for a family with no residual scale parameter is reported as
  declared-but-not-applied rather than silently dropped.

* **Beta, negative-binomial and binary responses fit on INLA.** The
  family gate compared a flexyBayes spelling against INLA's own
  likelihood roster, and the reconciler that maps between the two
  vocabularies runs after the gate -- so `negative_binomial`, `negbinom`
  and `binary` were refused for a naming reason on an engine that
  carries `nbinomial` and `binomial` perfectly well, three of the eight
  supported spellings. The gate now resolves the spelling first, and a
  refusal names both. Separately, the residual hyperparameter for the
  beta likelihood was emitted as `prec` where INLA declares `phi`, so a
  beta fit died inside the engine with a raw error whose own text quoted
  the correct alternatives. The keyword now comes from a table read out
  of `INLA::inla.models()`, and a test re-reads that declaration from
  the installed engine. All sixteen family x backend cells the
  capability matrix advertises now fit.

* **`(1 | g:h)` refuses on INLA exactly as `random = ~ g:h` does.** The
  two spellings are one model, and the README says the grammar is
  detected from the call shape. On INLA they diverged: the ASReml
  spelling was stopped by `lgm_gate()` with a typed refusal naming the
  engine that fits it, while the bar spelling was read as a simple
  random intercept on a group called `g:h`, passed the gate, reached the
  INLA binary, and failed there with a message describing the emitted
  formula rather than the user's model. Both surfaces are now classified
  by the same walker, so the term descriptors are identical and every
  downstream gate sees one model -- including the guard against a
  numeric variable inside a random interaction, which now fires on the
  bar surface too.

* **A prior the backend cannot carry is refused by name, instead of
  being dropped or raising a bare error.** `priors_to_inla()` translated
  three distributions on three targets and returned nothing at all for
  every other pair, while `prior_summary()` printed each of them as
  applied: `b()` and `cor()` targets were inert on INLA, as were five of
  the ten distributions on the variance components, so a
  prior-sensitivity analysis could vary a `student_t` scale and compare
  the engine's own default with itself. The brms side refused the same
  class of specification for a real reason but in an untyped condition,
  so no caller could tell "outside the translation table" from "R fell
  over". Which pairs each engine carries is now one table
  (`R/prior_translation.R`, the sibling of the family-traits table), and
  a pair outside it refuses at the fit, naming the row and both remedies
  -- re-express it in a distribution that engine does carry, or switch
  backend. The table also widened: brms takes all ten distributions on
  `sigma` and `sd(group = )`, and INLA gained `exponential` on the
  variance components (it is the PC prior in distributional form) and
  `normal` on `b()`, which arrives through `control.fixed` per
  coefficient. Every translatable pair carries an arrival test that
  reads the hyperparameter string back out of the fitted engine object
  rather than out of the translator.

* **A prior naming a term the model does not have refuses at fit entry,
  in the caller's vocabulary.** `b("nonexistent_term")` was passed
  through to brms's parser, which answered with a synthesised Stan
  parameter name (`b_nonexistent_term`) the user never wrote, untyped;
  on the INLA route the same mistake was ignored altogether. The check
  now runs against the engine's own parameter list -- `default_prior()`
  on brms, the design-matrix columns and f()-term keys on INLA -- and
  names the offending target and the model's actual terms.

* **A smoother written in the fixed part refuses instead of reaching the
  engine.** `y ~ s(x)` raised a bare `could not find function "s"` --
  from the INLA subprocess, and on brms only after a completed sampling
  run, so the compute was spent before the failure. Every mgcv-style
  smoother spelling (`s`, `te`, `ti`, `t2`, `sos`, `gp`) now refuses at
  parse time naming the flexyBayes spelling, `random = ~ spl(x)` on the
  INLA backend. Ordinary transformations in the fixed part are
  untouched.

* **A spline whose variable is also a fixed term refuses typed on
  INLA.** The emit placed `x` in the fixed part and as `f(x, model =
  "rw2")`, and INLA refuses a key used twice -- reaching the caller as
  the generic "the inla-program exited with an error" wrapper, which
  names neither the variable nor the remedy. The refusal now names both
  of INLA's own remedies and the one specific to an rw2 smooth, whose
  null space already carries the linear trend.

* **A bar-grammar random slope on a factor names the ASReml surface that
  fits.** `(trt | g)` raised an untyped error listing other bar
  spellings, none of which expresses the model, and never mentioned that
  the identical model written `random = ~ us(trt):g` completes on brms.
  The refusal is typed and names that surface verbatim -- `us(f):g` for
  the correlated form, `diag(f):g` for the uncorrelated one -- and a test
  fits both, so the alternative named is one that works. The remaining
  ingest refusals were typed alongside it.

* **Every aggregated-route refusal can be caught by class.** Three
  refusals on the way into the aggregated path raised a bare error --
  `aggregate = TRUE` with a missing response, `aggregate = TRUE` on a
  backend that carries no aggregated emit, and `aggregate = TRUE` with
  no aggregated backend resolving for the model -- while the two
  refusals further along the same path already carried a condition
  class. A caller could match three of the five by class and had to
  match the other two by message text. All five are typed now, on two
  reason codes that separate the boundaries: a missing response is a
  property of the data (`aggregation_response_incomplete`), an
  unreachable emit is a property of the engine roster
  (`aggregation_route_unavailable`). Message wording is unchanged, and
  each refusal still carries `flexybayes_aggregate_refusal` beneath its
  per-code class, so existing handlers keep working.

## New features

* **`family = "hurdle_gamma"` fits on the brms backend.** The entry
  allowlist refused it, and the engine behind it carries it natively
  (`brms::brmsfamily("hurdle_gamma")` declares `mu`, `shape` and `hu`) --
  a restriction over the engine with no documented boundary and no named
  alternative. It is admitted end to end: entry, family traits, emit,
  the capability matrix and a live fit. INLA's likelihood roster has no
  counterpart, so the family gate refuses it there and `backend = "auto"`
  routes to brms.

* **The remaining family boundaries say whose they are.** Three
  neighbouring families a field engagement reached for stay refused, and
  the refusal now records which layer the boundary sits at, checked
  against the installed engines rather than recalled: neither engine
  carries `zero_inflated_gamma` or `compound_poisson`, and INLA's
  likelihood roster does carry `tweedie` -- so that one is a flexyBayes
  boundary rather than an engine one, and is tracked as a feature
  request in `inst/KNOWN_ISSUES.md`. Each of the three names
  `hurdle_gamma` as the nearest implemented alternative.

## Documentation

* **The fixed-effect prior contract is corrected in the two vignettes
  that still stated the old one.** *Getting started* and *Default priors
  and the `fb_prior()` DSL* both described a blanket
  `normal(0, 100)` on every fixed effect as the working default. That
  was never what the engines received, and it is not what this release
  applies: `prior_fixed_sd` is honoured when it is supplied, and left
  unsupplied each engine keeps its own fixed-effect prior, which
  `prior_summary()` names. The reference documentation was corrected
  when the behaviour was fixed and the vignettes are now consistent
  with it.

* **The priors vignette carries the per-backend translation contract.**
  Section 6 previously described the cross-engine translation in a
  parenthesis that no longer matched the code. It now carries the
  target-by-distribution table for each engine, the argument each
  translating pair arrives through, and the reason `cor()` is a refusal
  on both engines rather than a gap.

* **The scaling claim now points at its measurements, and the claims
  that had none are gone.** `DESCRIPTION` said the package "scales to
  large agricultural datasets" and the README called a route "scalable"
  and the planner the "fastest way" to explore the package, none of
  which named evidence. The scaling claim is real and now says what
  produces it -- exact per-cell sufficient statistics on replicated
  designs -- and the measurements behind it ship at
  `inst/validation/benchmark_scaling.md`: the size at which the per-row
  path runs out of memory, what the streamed path costs at the same
  size, five billion rows through a flat memory envelope, and the model
  scope outside which the route refuses rather than approximates. The
  three claims with nothing behind them were reworded to say what is
  actually true (the planner needs no backend installed), because a
  quieter claim beats a strong one with no measurement under it.

## Internal

* **The capability claims are verified by execution, not by the
  registry that generates them.** `tools/execution_grid.R` runs the
  claim-derived cross product -- every capability-matrix cell, every
  prior route by family and backend, the malformed-prior corpus, the
  family allowlist against the installed engines' own rosters, the
  structure and grammar surfaces, and the prior-translation table --
  against an installed build, one isolated process per cell, and banks
  the per-cell ledger under `inst/validation/execution_grid/`.

* **`inst/validation/scenarios.yaml` registers the numerical-validation
  studies** in the schema the validation ladder reads, and `DESCRIPTION`
  declares `Config/rpkg/validationTier: V3`. Every registered row is an
  existing study. What the tier asks for and the studies do not yet
  supply is written down in `inst/validation/README.md`.

# flexyBayes 0.9.1

Work aimed at the reader who arrives from ASReml with last season's
`asreml()` call in hand, and who should be able to read the answer
without opening `$inla` or `$brms`.

## New features

* **`summary()` returns one object on every engine and every
  representation.** The two active engines used to return two
  incomparable things -- INLA's own four-slot list of its tables, brms's
  `brmssummary` -- and a fit run on the aggregated representation
  returned a third, its own list of the aggregated posterior's raw
  pieces. None of the three carried a variance-component table, so
  `summary(fit)$varcomp`, which is the first subset an ASReml user
  types, was `NULL` on every fit. It is now a data frame of `component`,
  `estimate`, `std.error`, `conf.low`, `conf.high`, `prior`, `note`,
  alongside `$fixed`, `$converge`, `$n_design`, `$n_observed` and the
  rest of an eleven-slot `summary.flexybayes` object, with a twelfth
  `$spatial_field` on a fit carrying an autoregressive latent field.
  Posterior mean is not a REML component.

* **The aggregated representation answers the same questions as the
  per-row one.** Aggregation is the default route for an ordinary
  Gaussian call with a random term, so the missing variance-component
  table was the ordinary case rather than a corner of it. An aggregated
  fit now returns the same object, with the components rebuilt from the
  engine's own posterior -- on the INLA route through the same precision
  marginal the per-row route uses, not a reciprocal square root of a
  tabulated mean. The printed banner names the representation
  (`aggregated-gaussian`, `aggregated-binomial`, `aggregated-poisson`)
  where a per-row fit names its engine, and the row-to-cell compression
  appears both on the header and in the `$model` slot. `plot(fit, type =
  "variance")` reads the same table and stopped erroring on an
  aggregated fit with it. The raw aggregated pieces the old summary
  returned -- `beta_means`, `beta_vcov`, `sigma_means`, `tau_means` --
  are unchanged at `fit$extras$summary`.

* **No convergence claim is made where no check ran.** The aggregated
  route records no numerical confirm, and the summary reported it as a
  failed one, because an absent record and a failed check tested the
  same. The line is now omitted where nothing was checked, and the
  `$converge$numerical_confirm` slot is `NA` rather than `FALSE`.

* **The variance components carry the prior each one ran under.** The
  `prior` column is a projection of `prior_summary()`, so the two
  surfaces cannot disagree about which prior reached the engine.

* **INLA's components reach the standard-deviation scale through the
  precision marginal** -- `inla.tmarginal()` then `inla.emarginal()` and
  `inla.qmarginal()` -- rather than by transforming a tabulated point
  estimate. The transform is nonlinear and decreasing, so the two
  answers differ, and they differ most where the posterior is widest.

* **`nobs(fit, type = "observed")`** reports how many responses were
  actually observed, next to `nobs(fit)`, which stays the design row
  count the engine saw. On an augmented fit those are different numbers
  and both are now visible.

* **`na_action` accepts the policy an ASReml user already writes**:
  `asreml::na.method(y = , x = )`, or the bare `list(y = , x = )` for
  readers without an asreml licence, alongside the native strings.
  Detection is by shape, never by class. The covariate half becomes
  settable in ASReml's own words: `x = "fail"` (the default, as in
  ASReml) refuses; `x = "omit"` drops the affected rows with a warning
  naming the count and the columns; `x = "include"`, which in ASReml
  means treating a missing covariate as zero, is refused by name.

* **`coef(fit, what = )` reaches the random effects and the unobserved
  cells**, and `ranef(fit)` is the same table under the name an
  `nlme` or `lme4` reader types. `what = "random"` returns one data
  frame per grouping factor with `group`, `level`, `estimate`,
  `std.error`, `conf.low` and `conf.high`, `what = "missing"` returns
  the unobserved design cells, and `what = "all"` returns all three.
  The default is the historical named numeric vector of fixed effects,
  byte for byte, so every existing caller keeps working. Where `nlme`
  or `lme4` is attached in the same session, `flexyBayes::ranef(fit)`
  is the collision-proof spelling.

* **`predict(fit, classify = )`** builds ASReml's marginal-means table
  on top of the existing `emmeans` seam, with a `level` argument
  defaulting to 0.95. It carries means and credible intervals and says
  so on the table. There is no pairwise standard-error block. On an
  INLA fit the interval comes from the Gaussian approximation of the
  fixed effects, and the banner states that rather than leaving it to
  be assumed. Without `emmeans` installed the argument refuses by name
  (`classify_requires_emmeans`) instead of failing inside a namespace
  call.

* **`fb_complete_grid()`** reinstates design cells that are absent from
  the data frame altogether, returning the crossing of the index
  variables with the response `NA` on every cell it added. It shares
  one implementation with the completion that `na_action = "augment"`
  performs, so the helper and the fit cannot drift apart. A design
  factor that varies across the trial refuses rather than being
  invented, and `unused_level =` opens that door explicitly, warning
  with every column it fills. That is ASReml's nin89 LANCER coding
  turned into a decision with a name on it.

* **`plot(fit, type = "variogram")`** draws the empirical residual
  semivariance over the design array, computed on the observed rows,
  because a residual is `NA` on an augmented row by construction. A fit
  with no design index refuses by name
  (`variogram_requires_design_index`). Nothing is fitted to the
  surface, so there is no fitted-variogram overlay to over-read.

* **The Bayesian door.** The same fit object now answers to the
  generics a Bayesian reaches for. `posterior::as_draws_df()`,
  `as_draws()` and `as_draws_matrix()` return canonical parameter names
  with the variance components on the standard-deviation scale, the
  same names `triangulate()` compares, on either engine. `loo()` passes
  through to `brms::loo()` on a sampled fit. `pp_check()` passes
  through to `brms::pp_check()` on a sampled fit, replicated datasets
  and all. On a Laplace fit each refuses by name and the message reads
  the fit rather than a memory of it, naming the WAIC and DIC that fit
  does carry and where they live. `prior_summary()` and
  `summary(fit)$converge` complete the door, and neither idiom is a
  wrapper around the other.

* **A term written on both the fixed and the random side is refused by
  name.** `flexybayes(yield ~ Variety, random = ~ Block + Variety)` is
  aliased with itself: the fixed part already estimates one mean per
  level, so the random copy's deviations are held up by their prior
  alone and the variance component reported for them reads that prior
  rather than the data. ASReml accepts the spelling and fits it, so a
  translated script arrives carrying it, and every engine here answered
  it badly in a different way -- an intermittently singular marginal
  solve on INLA, an empty variance-component table, and at benchmark
  scale two segmentation faults in the engine subprocess followed by an
  error with no reason code. The new `term_in_fixed_and_random` refusal
  fires at plan time on every engine and in `fb_plan()`, names the term,
  and gives both repairs: drop it from `random` for population-level
  means, or drop it from the fixed part for shrunk level effects.

* **The autoregressive field says so when it does not identify.** On an
  incomplete grid the latent field can collapse to a boundary -- field
  standard deviation at a floor, both correlations spanning almost the
  whole of `[-1, 1]`, its variance absorbed into the nugget -- while the
  convergence block reports a converged mode and a passing numerical
  confirm, because the optimiser did converge, to a solution with no
  field in it. A fit in that state now raises a warning naming the
  parameters that lost their identification and the three routes out
  (complete the grid, read the model on brms as well, or give the field
  an informative prior), and `summary(fit)$varcomp` carries `collapsed`
  in the `note` cell of the `sd_spatial` row. The warning is scoped to
  the field: a grouping factor whose variance component collapses is a
  different fact and keeps its quieter cell. Silence with
  `options(flexyBayes.silence_spatial_collapse_warning = TRUE)`.
  `inst/KNOWN_ISSUES.md` carries the measured behaviour.

* **The sectioned residual reaches the object, not only the printout.**
  `summary(fit)$varcomp` on a `dsum(~ units | f)` fit returned the
  grouping factors and no residual at all, while `print(summary(fit))`
  rendered a full per-level block -- so the entire point of fitting
  `dsum()` was readable and not subsettable, and reaching it meant
  opening `$brms` and transforming `b_sigma_*` draws by hand. The table
  now carries one `sigma_<level>` row per level of the sectioning
  factor, on the standard-deviation scale, taken from the same builder
  the printed block reads, with the `prior` cell naming the retargeted
  prior that actually reached the sampler.

* **`predict(classify = )` for a random-effects grouping factor refuses
  by name.** Asking a multi-environment fit for genotype means is the
  first thing a breeder does after fitting. The reference grid is built
  from the population-level design, so a factor entering only as a
  grouping term is not in it and the call died inside `emmeans` with
  "No variable named Geno in the reference grid" -- a third-party
  message with no reason code and no route onward. The new
  `classify_random_factor_not_supported` refusal names the factor and
  points at `coef(fit, what = "random")` and `ranef(fit)`, which carry
  the level effects with their intervals. Population-level marginal
  means for a random factor are planned.

## Breaking changes

* **`plot(fit, type = "pp_check")` no longer draws an
  observed-versus-fitted panel on a fit that carries no predictive
  draws.** That panel was never a posterior predictive check, and the
  documentation claimed it was. It is removed rather than retitled, and
  both entry points now raise a catchable refusal
  (`pp_check_requires_predictive_draws`) naming the diagnostics the fit
  does answer, `plot(fit, type = "residuals")` and, where there is a
  design index, `plot(fit, type = "variogram")`. On a brms-engine fit
  the same type now runs a real posterior predictive check through
  `brms::pp_check()`. A script that called it on an INLA fit for the
  picture will need one of the named alternatives.

## Minor improvements and fixes

* **A sectioned-residual fit answers the mean-model accessors.** A
  `dsum()` residual makes the emitted model distributional, and two
  readers were working on the wrong object: the formula reader indexed
  the resulting `brmsformula` by position, which reaches its parameter
  slots rather than its formula, and the coefficient reader swept every
  `b_` column, which collects the log-sigma coefficients alongside the
  mean effects. The result was `coef(fit)` reporting `sigma_EnvE1` as
  though it were an effect on the response, a fixed-effect design matrix
  that could not be reconciled with it, `$glm$y` and `$glm$residuals`
  all-`NA`, and `predict(fit, classify = )` failing in the estimability
  seam even for a plain fixed factor. All four are fixed at the two
  readers, so `coef()`, `vcov()`, `confint()` and the marginal-means
  table now describe the same parameter set.

* **An ordinary aggregated fit no longer reports a prior nobody
  supplied.** The `prior_parametrization` label was set by a class test
  on the recorded prior, which was right while the only `fb_prior` on a
  fit was one the caller wrote. The auto-default became an `fb_prior`
  object, so the commonest fit the package produces printed `custom
  (explicit prior supplied; see prior_summary())`. The label is now
  routed off the default-basis attribute `prior_summary()` reads, and
  there are three values rather than two: `per_row_equivalent` for the
  legacy scalar bridge, whose matched-prior guarantee is what that word
  claims; `package_default` for the automatic bounded uniform on the
  standard-deviation scale, which claims no such equivalence on the
  aggregated route; and `custom` for a prior from the caller.

* **`tidy(fit, effects = "random")` answers on an INLA fit.**
  `tidy.flexybayes_inla()` had no `effects` argument, so the request was
  absorbed by `...` and the fixed-effect table came back under a call
  that asked for variance components -- a wrong answer rather than an
  error, on the package's default engine. Both tidiers now take the same
  argument, route `"random"` through one builder, and stop on an
  unrecognised value. The shared builder reads through the normaliser
  the summary uses, so the request also works on an aggregated fit,
  where the raw field is a list and `nrow()` of it is `NULL`.

* **The two counts named `K` say which is which.** A fit reports the
  planner's cell count -- the product of the declared factor levels, the
  complete grid -- in `backend_decision(fit)$reason`, and the realised
  count of cell keys the data contains in the print and summary
  compression lines. On an incomplete grid these differ, and both
  printed as a bare `K`. The planner's figure is now labelled
  `K = <n> (estimated)`, in the routing reason and in the `fb_plan()`
  print; the realised count is unqualified.

* **The three prints share one header** and name the engine for what it
  did. An INLA fit no longer reports `MCMC : n chain(s) x n samples`
  over a deterministic Laplace approximation that ran no chains. The
  model line is derived from the model representation, so it describes
  the model that was written rather than the program that was emitted.

* **`update()` works on an INLA fit.** The INLA emit recorded six of the
  arguments a re-fit needs against brms's fifteen, and that gap -- not
  anything about the engine -- was the whole of the
  `update_call_not_reconstructable` refusal. Both engines now also record
  the requested `na_action`.

* **`update()` works on the aggregated routes too.** The Gaussian and
  count streaming-aggregation emits recorded eleven of the sixteen
  fields, so an ordinary Gaussian call with a random term, which routes
  there by default, refused a re-fit. All four emits now record the same
  names in the same order. `na_action` travels with the re-fit
  rather than being re-defaulted, so an augmented fit re-fits as an
  augmented fit. A fit made before this release still refuses cleanly on
  the shorter record, which is the intended behaviour and not a case to
  special-case into a silent default.

* **A re-fit runs under the prior the fit ran under.** The auto-default
  bounded uniform on the standard-deviation scale fires only when a call
  supplies neither a `prior` nor a `prior_vc_sd`, and `update()` re-issued
  the recorded `prior_vc_sd`, so the default never fired on a re-fit and
  the model fell through to the engine's own hyperprior. An identity
  `update()` -- nothing changed -- could halve a reported variance
  component and report `engine default` where the first fit reported the
  resolved uniform. `update()` now resolves the prior from the fit's own
  record, on every engine, and **a policy re-fires while a bespoke prior
  carries**. The auto-default is a policy -- one bounded uniform per
  variance component, with the bound read off the data -- so a re-fit
  passes neither `prior` nor either scalar and lets the default rebuild
  itself over the *updated* model. On an identity re-fit the same data
  give the same bound, so nothing moves. On `update(fit, random = ~ g +
  h)` the added term gets the same uniform as its siblings instead of
  falling to the engine while they keep theirs. A user-supplied
  `fb_prior()` is carried verbatim, and a term it never named still
  follows the engine's own hyperprior, exactly as on the first fit. A
  `prior` or a `prior_vc_sd` written in the `update()` call still
  outranks the record, and only one of the two forms is passed.

* **A variance component the package priored nothing for says so.** A
  partial `fb_prior()` -- one naming `sigma` but not a random term the
  model has -- left that component's `prior` cell blank in
  `summary(fit)$varcomp`, where the two words `engine default` belong.
  Reachable on a first fit, not only through `update()`.

* **A re-fit comes back on the engine and in the representation the fit
  ran on.** `backend` and `aggregate` were absent from the record, so a
  re-fit took the formal defaults for both: an identity `update()` of a
  Stan fit returned an aggregated INLA fit -- a different inference engine
  and a different representation, under the same name and with nothing
  said -- and a re-fit of a per-row fit could come back aggregated. Both
  now travel with the record, and what
  is recorded is the value the call asked for rather than the engine it
  resolved to: a fit made under `backend = "auto"` records `"auto"` and a
  re-fit routes again from there, because the recorded value is a policy
  and a re-fit whose model has changed has to be free to route the changed
  model. Both join the fields a re-fit requires, so a fit made before this
  release refuses cleanly rather than re-fitting elsewhere. An override in
  the `update()` call still wins. `verbose` is recorded alongside them and
  deliberately not carried -- `update()` reproduces the model and the
  policy behind it, not the display settings of the session that first ran
  it, and `update(fit, verbose = FALSE)` quietens a re-fit on the spot.

* **An auto fallback no longer forwards sampler settings that are not
  there.** An INLA fit records no `n_samples`, `warmup` or `chains` --
  the nested Laplace approximation runs none -- so an auto re-fit whose
  INLA attempt failed reached brms with all three absent, and brms
  refused the `iter` it builds from them (`Cannot coerce 'iter' to a
  single numeric value`). The user read a brms argument error where an
  INLA failure had happened. Absent sampler settings are now omitted and
  the receiving engine applies its own defaults.

* **A variance component priored through the legacy scalar bridge names
  that prior.** The cell read `engine default` on every row of such a
  fit, which on the Stan route is false -- the bridge writes
  `lognormal(0, prior_vc_sd)` onto the residual scale and onto every
  named group. `engine default` is now reserved for a component genuinely
  left to the engine, which on the INLA route is what the bridge leaves
  behind and where the words still belong.

* **`plot()`'s default type follows the engine, and the default draws on
  a brms fit.** The predicate deciding whether a fit carried draws read
  a greta-only slot, so a brms fit fell through to a message saying it
  had no draws to plot, and the diagnostics path had no brms branch
  underneath to reach anyway. The default is now `"residuals"` where the
  fit carries no sampler draws and `"diagnostics"` where it does, and
  the brms branch forwards to brms's own trace and density panels. Every
  explicit `type =` behaves as it did.

* **`na_action = "fail"` records itself.** With no missing response it
  used to fall through into the augment branch, completing the design
  grid and writing `augment` into the fit's record.

* **`missing_fraction` is present on every missingness path**, not only
  under `augment`, and means the same thing on each: the fraction of the
  fitted rows carrying no observed response.

* **The refusal messages for missing plots are rewritten** in the
  reader's language, each naming what ASReml does in the same situation
  and what to type instead.

* Effective sample size is recorded for the tails as well as the bulk of
  each marginal, and reported by `summary()`. A fit can mix well in the
  middle of a marginal and badly at the 2.5% bound, which is the number
  a credible interval is read off.

# flexyBayes 0.9.0

A stable release of the supported capability set, ready for testing --
published as a pre-release while the first collaborator round completes.
Entries appear newest-first within each section, so where a later change
supersedes an earlier one on the same line -- the greta quarantine over
the greta-only MET route, for instance -- the superseding entry is the
one above.

## Licence

* **The package is now licensed under MIT** (previously GPL >= 3), by the
  authors' decision of 2026-08-17. `LICENSE`, `LICENSE.md` and every
  metadata surface (`CITATION.cff`, `codemeta.json`, `.zenodo.json`)
  changed together.

## Vignette rendering

* **Inline math no longer towers over the sentence it sits in.** The
  vignette engine loads MathJax v2, whose TeX web fonts render visibly
  larger than the vignette body text. Every vignette now includes a
  shared header asset (`vignettes/_mathjax-scale.html`) that pins the
  MathJax output scale to the body font before the loader runs, and a
  test fails any vignette -- current or future -- that drops the
  reference.

## The plan, the build and the fit agree

* **`plan = TRUE` no longer names a route the fit would not take.** On a
  `residual = ~ dsum(~ units | env)` model the plan reported
  `Backend chosen: inla` and `Path: aggregated_inla` on the same screen as
  `Gate outcome: refuse_structural (residual_term_type_inla)`, while the
  live call routed to brms. Two defects combined: `.fb_aggregation_plan()`
  enumerated the structured residuals it could not aggregate and so let a
  sectioned residual through, and the plan's aggregation override treated
  an eligible aggregation plan as enough to choose the aggregated INLA
  route without asking whether the gate had accepted the model. The
  aggregation scope now names the one residual that IS cell-constant --
  plain `units` -- and the override requires an accepting gate.

* **An untracked file no longer decides what INLA fits.** The INLA gate
  for a factor-by-numeric interaction read an `.rds` under
  `inst/extdata/inla-verification/`. That directory is untracked but was
  not build-excluded, so a tarball built in a working tree carried a file
  a clean clone of the same commit did not -- and the two fitted different
  models. The directory is excluded from the build, and the artefacts are
  now a developer rehearsal hook read only under
  `options(flexyBayes.dev_inla_verification_artefacts = TRUE)`. Both
  affected term classes -- the factor-by-numeric fixed interaction and the
  uncorrelated random slope `(x || g)` -- refuse on every host, and both
  refusals now say why: their three-arbitrator verification named greta as
  one arbitrator, and greta is quarantined. The factor-by-numeric refusal
  recommended `backend = "greta"`; it names brms.

* **`tools/check_build_parity.R`** builds a tarball from the working tree
  and one from a clean `git archive` of HEAD and fails on any difference
  between their file lists. It also checks that every `.Rbuildignore` line
  compiles as a regular expression, because the file has no comment syntax
  and a `#` line containing an unbalanced parenthesis aborts `R CMD build`.

* **`backend = "auto"` on a host without INLA routes to brms.** It
  announced "routing to greta" and built a greta routing decision -- a
  branch that predates the quarantine, so every `auto` call on such a host
  advertised a quarantined engine and then refused. It now resolves
  through the same fallback the other auto paths use.

* **`family = binomial()` works.** A `stats` family object reached
  `tolower()` and a scalar `if`, and the user saw base R's "the condition
  has length > 1". Family objects and unevaluated family generators are
  accepted, and a family object supplies its own link, so `Gamma()` means
  the inverse link it names rather than the log link `family = "gamma"`
  defaults to. A family object handed a contradicting `link` refuses
  (`family_argument_not_recognised`) instead of choosing one.

* **The aggregation cell key is a set of variables, not a list of terms.**
  `y ~ a * b` produces three terms over two variables, and the planner
  multiplied the term-level counts, so a replicated 4-by-4 factorial that
  compresses 320 rows into 16 cells was sized at 256 and refused as
  unproductive. The runtime aggregator had always keyed on the union, so
  the plan and the emitter disagreed about the quantity the decision turns
  on.

* **A factor-interaction design column refuses the aggregated route.**
  Probing what the corrected planner then admits found the reason the
  first defect had been hiding: both aggregated emitters name
  model-matrix columns in the INLA formula, and INLA does not resolve a
  column called `a2:b2` even backticked. The aggregated route refuses by
  name and points at the per-row route, which uses INLA's own `a:b`
  notation and fits the model.

* **The aggregated header names the fit's own family.**
  `print()` and `summary()` on an aggregated fit said
  `aggregated-gaussian` from a string literal, so a binomial or Poisson
  fit from the count aggregator carried the wrong family on the surface
  the package uses to signal exactness.

## Fidelity of names

* **The separable AR1 spatial structure is now written on the random side.**
  This is a breaking change within the unreleased 0.9.0 line.

  ``` r
  # before
  flexybayes(yield ~ 1, residual = ~ ar1(row):ar1(col), data = trial)

  # after
  flexybayes(yield ~ 1, random = ~ ar1(row):ar1(col), data = trial)
  ```

  The INLA emit builds a latent autoregressive field plus the Gaussian
  observation nugget: four parameters on a field grid -- the row
  correlation, the column correlation, the field SD and the nugget SD.
  `residual = ~ ar1(row):ar1(col)` is ASReml's name for a *different* model,
  the nugget-free separable residual, which has three parameters and which
  the four-parameter model contains as the nugget goes to zero. On data with
  real plot-to-plot noise the two return different correlations, because a
  nugget-free model must absorb independent noise into the correlated
  structure. Emitting one under the other's name was the one place the
  package kept the original spelling for a model it was not fitting.

  A random-side field is a random effect, so the new spelling states the
  model truthfully. The residual spelling now refuses by name
  (`ar1_residual_not_representable`), naming both routes: the random-side
  field, or ASReml itself for the residual formulation. The emit, the grid
  conditions and the recovered values are unchanged -- only the spelling
  moves. brms refuses the field (`stan_cannot_represent_ar1_field`); it has
  no lowering for a Kronecker autoregressive precision.

* **`summary()` and `print()` report all four field parameters**, labelled
  and on the correlation and standard-deviation scales, rather than leaving
  INLA's precision-scale names (`Rho for row_id`, `GroupRho for row_id`) to
  be translated by the reader. The fourth parameter is exactly what
  separates this model from ASReml's, so a user who acts on the
  correlations can see what they are acting on.

* **`lgm_gate()`'s hyperparameter budget counts the field's own
  parameters** -- two for `ar1(t)`, three for the separable form, with the
  observation nugget counted on the family side. It previously counted the
  1-D field as one, and the budget is what decides whether INLA's numerical
  integration stays tractable.

* **The random / residual grammar is a closed vocabulary.** An unrecognised
  call was read as a plain variable whose name was its own source text, so
  `ar2(row)` became a lookup for a column called `"ar2(row)"` and
  `corgh(site):gen` became a crossing of two invented factors. Any call
  outside the recognised set now refuses (`asreml_function_not_recognised`),
  naming the token and the nearest spellings the parser does read.

* **The structures that parse for the formula catalogue but have no emit
  path are refused by name**: `ar2()` points at `ar1()`, `cor()` reuses the
  equicorrelation refusal, `str()` names the per-term structures to write
  instead, `fa()` names `us()` and `diag()`, and an interaction matching
  none of the recognised patterns -- `us(trait):vm(gen)`, the multi-trait
  idiom -- says so. All of these previously reached an `lgm_gate()` print on
  INLA and an untyped `stop()` on brms, two vocabularies for one fact.

* **The brms formula reconstruction refuses by name instead of stopping
  untyped** (`brms_cannot_represent_term`), naming the term type and the
  route that does fit it. `spl()` on brms is the common case and now points
  at INLA.

* **A numeric variable inside a nested or multi-way random interaction is
  refused** (`numeric_variable_in_random_interaction`). `~ Subject +
  Subject:Days` with numeric `Days` reads as a random regression to an
  ASReml user and emits as one independent deviation per Subject-by-Day
  cell -- 180 of them on `sleepstudy` -- and both readings are defensible.
  The refusal names them and the route to each: `(Days || Subject)` for the
  slope, an explicit factor conversion for the per-cell effect. A crossing
  of two factors is untouched.

## Triangulation is a gated diagnostic

* **`triangulate()` now gates the comparison before reporting it.** It
  compared any two objects that produced draws, on any parameters whose
  names happened to coincide, whatever model, data, prior or convergence
  state each carried. A distance between two such posteriors measures the
  difference between two questions, and there was nothing in the output
  that said so.

  Three gates now run first.

  - **Comparability.** Every fit records a model fingerprint at fit time:
    the canonical formula triple, family and link, the data dimensions,
    column names and a content digest, and the prior recorded on each
    variance component. Two fits whose fingerprints disagree raise
    `triangulate_incomparable_fits`, naming the first element that
    differs. There is deliberately no override argument -- comparing two
    different models is a different operation. A fit built outside the
    package carries no fingerprint, so comparability cannot be verified
    and the status is `inconclusive` rather than silently assumed.
  - **Diagnostics.** A brms fit passes on R-hat at most 1.01, bulk ESS at
    least 400 and no divergent transitions; an INLA fit passes its own
    numerical-confirm gate. A failing fit yields `inconclusive` and no
    parameter verdict at all, because disagreement between an unconverged
    fit and a converged one is not a finding.
  - **Matched priors.** A variance component is compared only when both
    fits record the same prior for it. The rest come back `not_compared`
    with the reason -- most often that the term is outside the
    default-prior walker, so each engine chose its own hyperprior.

* **The result carries a status**, `concordant` / `discordant` /
  `inconclusive`, plus a per-parameter verdict, the mean shift and the
  Wasserstein-1 distance in posterior-SD units. The thresholds behind the
  verdicts (0.1 SD, 0.1 SD, an SD ratio inside [0.9, 1.1]) are documented
  as heuristics for directing attention, in the sense Gelman et al. (2020)
  use for cross-implementation checks. They carry no error rate, and the
  man page states plainly that the package ships no simulation-based
  calibration of its own fits, so `concordant` is not a calibration claim.

* **The man page and the README stop presenting triangulation as the
  package's validation contract.** It is a diagnostic on a gated overlap.
  Agreement between two engines reading the same parsed model cannot
  detect a mistranslation of that model, and both surfaces now say so
  where the function is introduced rather than in a footnote.

## Priors reach more of the model

* **The default uniform-on-SD prior now reaches nested, `diag`/`idh`/`at`
  and `us` group terms** on the brms path. It previously walked simple,
  uncorrelated-slope, `vm` and `ped` terms only, which meant the terms
  most likely to be triangulated -- multi-stratum and
  genotype-by-environment variances -- were the ones on which the two
  engines were never asked the same question. Fits of those models with no
  explicit `prior` will shift: the variance components now carry
  `uniform(0, U)` where they previously carried whatever the engine chose.

* **Terms still outside the walker are recorded as engine defaults**
  rather than passed over in silence: the `us` correlation block (which
  keeps brms's LKJ), multi-way random interactions, and the AR1 latent
  field (which keeps INLA's own hyperprior). The record is what lets
  `triangulate()` exclude them by name instead of comparing them.

* **The heterogeneous-residual prior retarget now uses the multiplier that
  actually produced its bound.** With `sigma ~ 0 + f` the scalar residual
  prior has no parameter to attach to, so it is retargeted onto the
  log-sigma coefficients as `normal(log(sd(y)), 1)`. The code recovered
  `sd(y)` by dividing the uniform bound by 2.5 -- the superseded
  penalised-complexity default's multiplier -- while the uniform default
  had moved to `U = 5 * sd(y)`. The prior's median residual SD was
  therefore twice the data SD rather than the data SD its own rationale
  states. The bound `U` is unchanged; the divisor now reads from the same
  constant that sets it, and heterogeneous-residual fits under the default
  prior will shift slightly.

## The combined multi-environment-trial model fits

* **Interaction random effects and a sectioned residual now fit together**,
  on brms, and the capability table says so. The two halves had been tested
  apart, the combination only at the emit level, and every public surface
  said as much -- the matrix row read `emits`, meaning the Stan program
  carried both blocks and nobody had run it. It runs:
  `tests/testthat/test-met-combined.R` asserts both blocks in the emitted
  formula and in the Stan program, then fits the model on simulated
  multi-environment data and asserts every R-hat below 1.05, no divergent
  transitions, and per-environment residual variances recovered in the order
  they were simulated in. `auto` reaches brms for this class; INLA still
  refuses it structurally.

  The verified fit is 120 rows. `flexybayes()` forwards no `control` list to
  Stan, so a trial series large enough to need a raised `adapt_delta` has no
  route to one, and that is the next thing this capability needs.

* **`summary()` and `print()` report the per-level residual variances** on a
  fit with a sectioned residual. `sigma ~ 0 + f` puts the predictor on log
  sigma, so a `b_sigma_*` coefficient is neither a variance nor a standard
  deviation, and the sectioning leaves no scalar `sigma` for the
  variance-component table to report -- the residual was simply absent from
  the printed output of a model whose whole point was the residual. Both
  quantities are now summarised from the draws of `exp(2 * b)` and
  `exp(b)`, as posterior medians with a credible interval, labelled by the
  factor's own levels. They are not point transforms of a posterior mean:
  `exp()` is convex, so those are different numbers.

* **`fb_plan()` no longer reports that a model will not fit and then fits
  it.** Under `backend = "auto"` the plan read the latent-Gaussian gate's
  verdict on INLA as the fit's verdict, so every model the router hands to
  brms -- multi-stratum designs, heterogeneous variances, the combined MET
  model above -- printed `Will fit: no (preflight refused)` while also
  printing `Backend chosen: brms`. Both halves of that line were wrong: the
  preflight had sized every term, and the chosen route fits. `will_fit` now
  follows the resolved route, and where there genuinely is none the printed
  reason names which dead end it was.

* **The plan's representation table includes the residual structure.** A
  sectioned residual is half of a multi-environment-trial model and did not
  appear in the plan at all.

* **The convergence warning stopped misdiagnosing `us()` terms.** It blamed
  non-identified factor-analytic loadings, which does not apply: brms
  parameterises `us(f):g` by standard deviations and a Cholesky correlation
  factor, and those converge. What fails at one observation per cell is the
  split between the covariance diagonal and the residual, which are
  confounded there -- and the old text sent the reader to
  `fb_structured_cov()`, which abstains for `us()` terms. The note now names
  the confounding and the two remedies that work on it, replication within
  cell or an informative prior on the residual.

## Heterogeneous variances

* **`diag(f):g`, `idh(f):g` and `at(f):g` now fit**, on the brms backend: a
  separate variance for every level of the outer factor, with no covariance
  between them. This is the workhorse structure of a multi-environment trial,
  where the genotype variance differs by site, and it previously refused on
  both active backends -- it had been a greta capability, and greta is
  quarantined. ASReml treats `diag()` and `idh()` as the same structure, and
  so does flexyBayes: all three spellings emit identical code.

* **`us(f):g`** emits the correlated sibling, one character apart in the brms
  formula and the whole difference between *k* parameters and *k(k+1)/2*.

* **Heterogeneous residuals**: `dsum(~ units | f)`, and the `at(f):units`
  spelling that parses to the same node, now lower to distributional
  regression on the residual scale. Refused for families that have no
  residual scale parameter to vary -- a Poisson's dispersion is a function of
  its mean, so there is nothing there to make heterogeneous.

* **`at(f, level):g` is refused rather than treated as diagonal.**
  Conditioning a random effect on selected levels is a different model from
  varying its variance across all of them; the two share a spelling, which is
  exactly how one gets fitted under the other's name.

  The mappings are validated against ASReml before use, not after, in
  `oracle_heterogeneous.R`: ASReml's `diag(site):gen`, an explicit `lme4`
  expansion and the brms emit are fitted on the same data and checked on the
  PARAMETER COUNT before the values, because estimates alone cannot tell a
  diagonal structure from an unstructured one when the true cross-level
  correlations are zero. The two REML arms agree to 3.5e-05; all three fit
  exactly *k* variances and no covariances.

  The residual side has its own probe in the same script, with four arms so
  that neither Bayesian arm is its own oracle. Fitted through `flexybayes()`
  on 360 simulated plots, the per-site posterior-mean residual variances are
  0.1146 / 1.1516 / 4.6981 where ASReml's `dsum` gives
  0.1093 / 1.1248 / 4.6079 and `nlme::varIdent` reproduces ASReml to
  1.3e-07. The largest gap is 4.8 per cent, on the smallest of the three
  variances, which is where a posterior mean and a REML point estimate are
  least alike. Every arm fits three free residual scales, no covariance, and
  no scalar `sigma` -- the structure is asserted before any value is
  compared. Run of 2026-08-15, `results/oracle_het.log`.

## Missing responses preserve the design

* **New `na_action` argument on `flexybayes()` / `fb()`**, defaulting to
  `"augment"`. A missing response removes an observation, not a design node.
  Where the model carries a covariance indexed by the design -- a separable
  AR1(row) x AR1(col) residual over a field trial -- deleting the row of a
  lost plot changes the index set the covariance is built over, so the model
  that gets fitted is no longer the model that was written down. Under
  `"augment"` such rows are retained and carried as latent quantities, and
  design cells absent from the data frame are reinstated where their model
  variables are determinable. This is ASReml's `na.method(y = "include")`
  behaviour. `"omit"` drops them (complete-case); `"fail"` refuses.

  The layer adds no inference machinery: INLA already treats an NA response
  as a latent prediction target and marginalises it. What is new is keeping
  the design intact so the engine is handed the right index set. Previously
  a single lost plot made an `ar1(row):ar1(col)` model unfittable, since the
  emit refuses an incomplete grid.

  Under ignorability the parameter posterior is the same whether the cell is
  imputed or omitted -- augmentation preserves the *representation*, not
  information. Where missingness depends on the unobserved response itself,
  both are biased and neither this argument nor any diagnostic here will say
  so. Missing *covariates* are refused under every setting: a filled-in
  predictor is a fabricated observation, which is a deliberate departure
  from ASReml's drop-or-zero-fill and a stricter one.

  Validated against oracles flexyBayes did not author, in
  `design-preserving-missingness/`: ASReml's own device, an independently
  written observed-data REML (so ASReml is not its own oracle), and
  complete-case `lme4`. What that sweep establishes is a **REML device
  identity** -- the ASReml augment device and the independent observed-data
  REML agree to about 1e-4 relative on every design -- and not an identity
  between a flexyBayes posterior and a REML estimate, which is a different
  claim and is not made here. Complete-case deletion reproduces the device
  on iid designs and departs on the correlated-residual one, which is the
  whole subject.

* The fitted object records what the layer did, under
  `$extras$na_action`. Whether a posterior was computed on the design as
  laid out or on the plots that survived was not otherwise recoverable.

## Model-fidelity fixes

* **The vignette deck merges from thirteen pages to eleven.** Two pairs
  overlapped enough that a reader had to hold both open. *Structured
  covariance* (05) is now part of the formula-surface page (02), which
  carries the whole term catalogue plus heterogeneous residual variance,
  design variance components, the `vm()` / `ped()` carriers and the
  genotype-by-environment family in one place. *Extending backends* (15)
  is now Section 9 of *dispatch and refusals* (11), next to the gate and
  the capability trail it describes, with its live registry dump and its
  availability-versus-capability split intact. No scientific content was
  dropped, the two shipped filenames are unchanged, and the numbering
  gaps at 05 and 12--15 record the merges. One stale claim was corrected
  on the way: `ar1(env):geno` has an INLA route, which the retired page
  denied.

* **The teaching pages pass `seed`.** Tutorial 01 told the reader that
  Stan's sampler was seeded while the call above it passed none, so the
  printed diagnostics could not be reproduced from the code shown. The
  brms fits in *getting started*, *multi-environment trials and genomic
  selection* and *cross-engine triangulation* now pass a seed, and two
  knits of each return identical output. The INLA half of *getting
  started* still does not reproduce -- the collapsed `Subject` precision
  lands orders of magnitude apart between runs, single-threaded as well
  as threaded -- and the page now says so, since instability across
  re-runs of a deterministic engine carries the same warning as a failed
  sampler diagnostic. Baked numbers on those three pages moved with the
  seeds.

* **The triangulation page's posterior-predictive check draws from the
  posterior predictive.** Section 5.3 simulated from the posterior-mean
  linear predictor with the group intercepts dropped, the same error
  Tutorial 08's `loo` repair had just fixed one page over. It now uses
  `brms::posterior_predict()` and keeps the mean surface beside it for
  contrast: the two are five-fold apart at the upper percentile. The
  page's lead also stopped saying that INLA silently approximates
  outside the latent-Gaussian class -- `lgm_gate()` refuses the model
  before INLA is called -- and the walk-back on what cross-engine
  agreement buys now leads rather than trails.

* **`flexybayes()` gains `seed` and `control`, so a brms fit is
  reproducible and its sampler is tunable.** Neither was reachable: the
  entry point took no `...`, so `set.seed()` before a call did not pin
  Stan's stream and `adapt_delta` could only be raised by dropping to
  `brms::brm()` directly. The gap was load-bearing rather than cosmetic --
  the package's own heterogeneous-variance oracle measured 0.27 percent
  run-to-run drift on the per-environment variances it quotes.

  ``` r
  flexybayes(
    yield ~ site, random = ~ gen, data = trial,
    seed = 20260815L, control = list(adapt_delta = 0.95)
  )
  ```

  Both are forwarded to `brms::brm()` and recorded on the fit, so
  `update()` repeats them. `NULL` maps onto brms's own defaults (`seed =
  NA`, `control = NULL`) rather than being passed through. The INLA path
  is a deterministic Laplace approximation with no random stream and no
  adaptation phase, so both are no-ops there and say so once per fit;
  silence the note with
  `options(flexyBayes.silence_sampler_arg_note = TRUE)`. The
  convergence warning now names `adapt_delta` as a route where divergent
  transitions are the symptom, and says plainly that on an unstructured
  term at one observation per cell it lowers the divergence count without
  identifying the split.

* **Tutorial 08's `loo` now includes the subject offsets, and shows why
  that matters.** The ELPD was computed from a pointwise log-likelihood
  built out of population-level coefficients only, which answers the
  leave-one-out question for a *new* subject rather than for a new
  observation on a subject already in the data. The tell was a `p_loo` of
  57 on a three-parameter mean structure. The vignette now takes the
  array from `brms::log_lik(fit$brms)` (`p_loo` 19.5, against 18 subjects
  plus three parameters), keeps the population-level version alongside as
  the mistake a reader would otherwise make, and tabulates the two.

* **Tutorial 03 stops promising mgcv smooths in its title, converts
  INLA's precision output, and checks its own recovery claim.** The title
  named a capability no active engine has -- `s(x)` parses and is refused
  -- and the text said the hyperparameter block reported a standard
  deviation where INLA had printed a precision. The residual scale is now
  converted to $\sigma = \tau^{-1/2}$, with the interval reversed as that
  monotone-decreasing map requires, and the vignette tests rather than
  asserts that the interval covers the simulation's 0.3.

* **Tutorials 15 and 16 each run something.** Both shipped entirely
  `eval = FALSE`. Tutorial 15 now knits its backend-registry dump live,
  so the roster a reader sees is read out of the installed package.
  Tutorial 16 gains a worked example small enough to execute at build
  time -- the same 720-row model fitted with `aggregate = TRUE` and
  `FALSE`, whose fixed-effect posteriors agree to six decimal places --
  is retitled to match what a reader can run, states the streaming
  envelope once (exact cell-likelihood aggregation for iid
  exponential-family models with a small cell count, from partitioned
  `.fst`, on INLA), and labels its billion-row tables as an external
  benchmark rather than as vignette output.

* **Augmenting past roughly 30 percent missing responses now says what
  that does to the estimand.** The augmentation identity is algebraic and
  holds at any missing fraction. What weakens is the restricted
  likelihood it targets: variance components reach the boundary, the
  surface flattens, and two correct implementations legitimately stop at
  different points. A one-time warning states that, so a user comparing
  the posterior against ASReml or lme4 at a high missing fraction reads a
  disagreement as a property of the design rather than as evidence
  against either fit. Silence with
  `options(flexyBayes.silence_high_missingness_warning = TRUE)`. The
  fraction is also recorded on the fit's missing-data metadata.

* **A streaming source past `2^31 - 1` rows refuses instead of recording
  `NA`.** Four sites cast the row count with `as.integer()`, which
  returns `NA` past the integer limit, so a five-billion-row single file
  would have produced an aggregation plan and a fitted object carrying no
  row count at all rather than an error. They now raise
  `row_count_exceeds_integer` and name the route that has no such
  ceiling -- partitioned `.fst` input, whose reader keeps its total as a
  double.

* **`flexybayes_stream()` no longer offers greta as a backend.** INLA is
  the only engine with an aggregated emit and is now the default and the
  only active choice. `"greta"` stays in the recognised vocabulary so
  passing it refuses by name (`backend_quarantined`) rather than failing
  on argument matching.

* **The greta S3 shim registers whichever order the packages load in.**
  The registration was guarded on `isNamespaceLoaded("greta")`, which
  answers only for the instant `.onLoad()` runs -- so a user who attached
  flexyBayes first and greta second, the usual order, never received it.
  A `setHook(packageEvent("greta", "onLoad"), ...)` registration covers
  that direction. The registrar carries its own loadedness guard, since
  `asNamespace()` would otherwise load greta and pay exactly the cost the
  guard exists to avoid.

* **`prior_summary()` reports the whole prior, not the half this package
  chose.** On a fit with a distributional residual (`dsum(~ units | f)`)
  it printed `sigma ~ uniform(0, U)` -- a parameter the model does not
  contain, under a distribution Stan did not use, because the declared
  uniform on the SD scale is retargeted onto the per-level log-sigma
  coefficients. It also said nothing about terms the default-prior walker
  does not reach, so a reader concluded the shared default had been
  applied to both variance components of a combined GxE model when it had
  been applied to one.

  The summary now names the retarget explicitly, lists every parameter
  carrying the engine's own default with the reason, and on a brms fit
  carries `brms::prior_summary()` of the fitted model -- brms is the
  authority on what reached Stan, so the check on the declaration is the
  engine's own record rather than a second reconstruction.

* **`fb_met_summary()` abstains by name.** It raised a bare `stop()` on
  both active engines, the last export doing so. It now refuses with
  `met_summary_not_available`, and points at `summary()` and
  `brms::VarCorr()` rather than at `fb_structured_cov()`, which abstains
  for every structure an active engine can fit. `fb_met_summary()` and
  `fb_log_posterior()` carry a lifecycle note recording that they are
  superseded while greta is quarantined, and both, with
  `prior_summary()`, gain live tests on an active engine -- all three
  previously had none that ran.

* **The LGM gate's refusals carry a class.** `lgm_gate()` returns an
  `<lgm_refusal>` value with a rule-by-rule check-list, and under an
  explicit `backend = "inla"` request dispatch turned it into a bare
  `stop()`. That left one refusal a caller could not catch by name while
  every other refusal in the package carried
  `flexybayes_refusal_<code>`. The message and the check-list are
  unchanged; the throw site now raises `inla_gate_refused`, carries the
  `<lgm_refusal>` object on the condition, and adds a
  `flexybayes_lgm_<rule_id>` class so one failing rule can be caught
  without catching the rest.

* **INLA fits join the shared class graph.** An INLA fit was
  `c("flexybayes_inla", "list")` while a brms fit was
  `c("flexybayes_brms", "flexybayes", "list")`, so the two engines had no
  common parent and every generic needed a parallel method. An INLA fit is
  now `c("flexybayes_inla", "flexybayes", "list")`.

  Five parent methods had no INLA sibling and so reached
  `stats::*.default` on an INLA fit. Each now resolves its inputs from the
  slots the object carries. `nobs()` and `model.matrix()` answer correctly.
  `confint()` reads INLA's own posterior marginals through the new
  `confint.flexybayes_inla()` method, where `stats::confint.default`
  previously returned a normal approximation built from `vcov()` without
  saying so. `update()` and `anova()` refuse by name
  (`update_call_not_reconstructable`,
  `conditional_loglik_not_available`): an INLA fit records six of the
  thirteen arguments a re-fit needs, and INLA reports a marginal
  log-likelihood rather than the conditional one a criterion would need.

* **`logLik()` refuses instead of returning a silent `NA`.** The
  `flexybayes` method wrapped its evaluation in `tryCatch()` and returned
  `NA` with a warning on any failure, and returned `NA` without comment for
  any family outside gaussian, binomial and poisson. `AIC()` and `anova()`
  then consumed that `NA` as though a log-likelihood had been computed.
  Both paths now raise a typed refusal naming the family or the missing
  slot.

* **The correlated-random-slope refusal now names the model to fit
  instead.** `(x | g)` refused with a generic pointer at "(x || g)", and
  Tutorials 01 and 04 filled the gap by teaching the ASReml crossing
  `~ Subject + Subject:Days` as the random-slope route. With a numeric
  `Days` that crossing is one independent deviation per Subject-by-Day
  cell -- 180 of them on `sleepstudy` -- not a per-subject slope. The
  refusal now spells out `(Days || Subject)` for the model at hand, says
  what that model is, and warns against the crossing by name. The
  condition's `workaround` slot carries the concrete spelling rather than
  a placeholder.

* **`library(flexyBayes)` no longer loads greta.** The load hook registers
  a handful of greta S3 methods for the quarantined aggregated-greta emit,
  and guarded that on `requireNamespace("greta")` -- which loads greta,
  reticulate and TensorFlow on every session, for a path dispatch cannot
  reach. The guard is now `isNamespaceLoaded("greta")`, so the shim
  registers when the user has loaded greta themselves and costs nothing
  otherwise. greta stops appearing in every `sessionInfo()`.

* **A requested structured residual can no longer be dropped from the fitted
  model.** `backend = "auto"` could fall back from INLA to brms on a model
  carrying an `ar1(row):ar1(col)` residual, and brms has no
  residual-covariance lowering at all -- so the emitted Stan program was an
  intercept-only independent Gaussian, with neither row/column indices nor
  correlation parameters, and no error was raised. Three routes reached that
  outcome: the automatic fallback after an INLA runtime failure, the
  code-inspection modes (`return_code` / `review_code`), which resolved
  `"auto"` to brms unconditionally, and an explicit `backend = "brms"`
  request, whose capability gate matched a single reason code by name and so
  failed open on any capability added later.

  The brms capability predicate now inspects `residual_terms` as well as
  `random_terms`, against a positive allowlist (only the iid `units`
  residual, which brms carries as the family scale parameter), so a residual
  form the parser learns later defaults to a refusal rather than to silent
  omission. Two refusal codes are registered:
  `stan_cannot_represent_ar1_residual` and
  `stan_cannot_represent_structured_residual`. The code-inspection modes now
  gate their `"auto"` to brms resolution on the same predicate, and the
  explicit-brms gate refuses on any capability failure, re-raising the
  predicate's own reason code. A complete grid still emits the separable
  AR1xAR1 field on INLA, unchanged.

  The regression test that covered this path asserted on the returned
  object's class and recorded backend, both of which were correct while the
  model was wrong. It now asserts on the emitted code, and a new
  `test-residual-structure-fidelity.R` holds that contract for both
  directions: a structured residual either appears in the emitted model or
  the call refuses.

* Three refusal messages still directed users to `backend = "greta"` as a
  workaround, which the quarantine had made unreachable.

## Reshape -- two-engine faithful core (brms + INLA)

* **greta and gretaR are quarantined as fitting engines** (reshape R1). The
  active backends are now brms and INLA. `backend = "greta"` / `"gretaR"`
  refuse with a structured `backend_quarantined` reason; `backend = "auto"`
  routes to brms (when it can represent the model) or INLA, and refuses with
  `auto_no_active_route` when neither can -- there is no silent greta fallback.
  A native greta model graph can no longer be fit
  (`native_greta_fit_quarantined`). The greta / gretaR registry descriptors,
  the emit code, and the legacy greta *import* grammar (`fb_from_greta()`,
  reading an already-fitted object's draws) are retained; re-entry is repair +
  conform to the backend conformance battery, never a bare re-add. The backend
  registry gains a `quarantined` lifecycle state alongside `active` / `dormant`.
* **Closed a gap in the "no silent greta fallback" guarantee above.** A
  DATA-dependent INLA runtime failure (for example
  `ar1_spatial_requires_complete_grid` on an incomplete field-trial grid) took
  a separate `backend = "auto"` code path from a structural `lgm_gate()`
  refusal, and that path still hard-coded a greta fallback -- unreachable via
  an explicit `backend = "greta"` request, but reachable through `"auto"` on
  a model INLA's gate accepts but its emit refuses at fit time. It now
  resolves through the same brms-or-refuse logic as every other auto
  fallback, with a new regression test covering it.
* **Swept the remaining stale "re-route to greta" advice** left over from the
  quarantine, found via a `--as-cran` release rehearsal and a full-suite test
  failure. `emit_brms()`'s correlated-random-slope refusal, `lgm_gate()`'s
  generic refusal print method, and the `emit_inla()` uncorrelated-random-
  slope verification deferral all pointed users at `backend = "greta"` as a
  workaround; the last of these no longer even matched reality, since
  `backend = "auto"` already falls back to brms on that exact refusal. All
  three now name the working route (brms) instead, with the corresponding
  tests updated to match. Also fixed the package startup banner and two
  `stop()` messages (`fb_gwas()`, `fb_met_summary()`) that still listed
  greta as an active engine.
* **All sixteen vignettes reshaped for the two-engine core**, restructured
  so eleven target a general audience and two (*dispatch and refusals*,
  *extending backends*) are the technical/internals reference; the old
  `flexyBayes-11-lgm-feasibility` / `-12-backend-internals` /
  `-13-lgm-feasibility-memory` / `-14-engine-selection` vignettes are merged
  into one *dispatch and refusals* vignette. README and `_pkgdown.yml`
  rewritten to match.

## New features

* **Separable AR1(row) x AR1(col) spatial models on INLA (WP16).** A designed
  field trial's spatial structure -- `ar1(row):ar1(col)`, or a 1D `ar1(t)`,
  written as a random effect or a residual structure -- now fits on INLA as the
  grouped-AR1 latent field `f(row, model = "ar1", group = col, control.group =
  list(model = "ar1"))` plus the Gaussian nugget, with canonical `rho_row` /
  `rho_col` / `sd_spatial` hyperparameters. `backend = "auto"` routes a designed
  spatial trial to INLA. The latent-field representation is faithful under one
  observation per grid node (validated against an independent GLS/REML oracle);
  an incomplete or replicated grid refuses
  (`ar1_spatial_requires_complete_grid`) rather than silently approximating.
  This recovers -- faithfully -- the spatial capability the greta quarantine
  removed (greta emitted `ar1(row):ar1(col)` as additive iid row/col effects).
  The *spatio-temporal models* vignette is rewritten for the new INLA-native
  path, including the complete-grid refusal on a real field trial with
  missing plots and the graceful `"auto"` fallback to brms when it fires.

* *(Superseded by the quarantine above.)* **An explicit `backend = "greta"`
  request now fits crossed interaction
  random effects and a heteroscedastic per-environment residual** -- the full
  ASReml-style multi-environment-trial (MET) shape (`random = ~ gen + gen:env`,
  `residual = ~ dsum(~ units | env)`). The greta code generator already gathered
  these terms; the dispatch preflight now sizes `nested` / `combo` interaction
  random intercepts (an index gather into a per-combination latent vector) so
  the plan clears and the fit proceeds. `backend = "greta"` is a deliberate
  opt-in for the greta path; see the next entry for how `backend = "auto"` now
  routes these models.
* **`backend = "auto"` routes a multi-stratum designed experiment to brms --
  the faithful full-HMC backend.** A model with interaction random effects
  (`gen:env`, `env:rep`, `env:rep:block`) is refused structurally by INLA,
  which collapses the finest variance components to zero. Rather than falling
  back to greta -- which under-mixes them -- `auto` now routes such a model to
  brms when brms is installed and can fit it, and brms recovers every variance
  component against the ASReml / lme4 REML reference (validated on
  `agridat::besag.met`). `emit_brms()` fits these interaction random effects
  natively as `(1 | A:B)`; INLA continues to refuse them by name. The
  fallback-to-greta clause this entry originally carried was removed by the
  quarantine above: when brms cannot represent the model, `auto` now refuses
  with `auto_no_active_route`.
* *(Superseded by the quarantine above.)* **The greta backend warm-starts the
  intercept from the response mean**
  (identity / log / logit / probit link scale), shortening the initial sampler
  transient while leaving the variance components and random effects at their
  prior-draw starts so convergence diagnostics stay informative.

## Bug fixes

* **Interaction random terms no longer crash the dispatch plan.** A three-way
  interaction random term (for example `gen:loc:yearf`) previously aborted the
  preflight with a base-R "the condition has length > 1" error, then a
  "no such index at level 1" error. The level-count helper now counts distinct
  observed combinations and the term label collapses to a single string, so an
  interaction random term reaches the normal (structural-refusal or fit) path.
* **A structured `dsum()` residual no longer silently drops its spatial
  structure.** `residual = ~ dsum(~ ar1(col):ar1(row) | env)` previously parsed
  to a per-region heteroscedastic variance, discarding the separable
  `ar1():ar1()` autocorrelation with no warning and fitting a different model
  than the one written. flexyBayes represents `dsum()` only as a per-region
  heteroscedastic variance (inner `units`), so a structured inner now raises a
  clear refusal naming the dropped structure, rather than silently reducing the
  model. Per-region structured residuals are planned for a future release.

## Breaking changes

* **The residual-structure argument `rcov` is renamed to `residual`, matching
  ASReml 4.** ASReml-R renamed this argument from `rcov` (ASReml 3) to
  `residual` (ASReml 4); flexyBayes now follows the ASReml 4 name across every
  entry point (`flexybayes()` / `fb()`, `fb_from_asreml()`, and the engine pins
  `fb_greta()` / `fb_inla()` / `fb_brms()`). Supplying the old `rcov =` argument
  now raises a guiding error pointing to `residual =`; update calls such as
  `rcov = ~ at(env):units` to `residual = ~ at(env):units`. The residual grammar
  itself is unchanged. The print method now labels the residual structure
  `Residual` rather than `Rcov`, and the two associated refusal identifiers are
  renamed to match (`rcov_type_unsupported_for_aggregation` ->
  `residual_type_unsupported_for_aggregation`, `rcov_term_type_inla` ->
  `residual_term_type_inla`).

# flexyBayes 0.8.3

A documentation-accuracy and ergonomics release on the 0.8.x line. There are no
modelling-behaviour changes; the additions are accessor coverage, clearer
refusals, and metadata / vignette consistency for the first public development
release.

## New features

* **`glance()` / `augment()` gain explicit INLA methods.** Calling `glance()`
  or `augment()` on an INLA fit (`flexybayes_inla`) previously raised a bare
  "no applicable method" error, because that class does not inherit
  `flexybayes`. They now dispatch to informative refusals pointing to `tidy()`
  (coefficient-level summaries) and `summary()` / `fb_structured_cov()`
  (variance components). Net: `tidy()` covers all three backends; `glance()` /
  `augment()` cover the greta and brms classes, with a clean message on INLA.

## Minor improvements and fixes

* **`fb_met_summary()` distinguishes the wrong-backend case.** Passing an INLA
  or brms fit now returns a backend-specific message (breeder summaries need a
  greta factor-analytic fit; the INLA / brms MET path reports variance
  components via `summary()` / `fb_structured_cov()`) rather than a generic
  "not a flexybayes object" error.
* **Development-release signalling and metadata reconciliation.** `.onAttach()`
  and `README` state that all exports are experimental and the package is not
  on CRAN. The package-level documentation, `CITATION` / `codemeta` / Zenodo
  metadata, `_pkgdown.yml`, `API_STABILITY.md`, and the security / support /
  contributor docs are reconciled to the 0.8.x line and the three-backend
  (greta / INLA / brms) surface.
* **Vignette convergence disclaimers** propagated to the remaining small-budget
  reference vignettes; every vignette that prints a high R-hat now carries an
  "illustration of output shape, not inference" callout.
* `NEWS` broom-coverage wording corrected; package-level roxygen updated to
  three engines and sixteen vignettes.

# flexyBayes 0.8.2

## New features

* **`fb_log_posterior()` -- a log-posterior producer for downstream tools.**
  A new exported generic turns a fitted flexyBayes object into a vectorised,
  domain-safe, unnormalised log-posterior callable that
  `proxymix::from_fb_posterior()` compresses into a closed-form
  Gaussian-mixture proxy. This is the single inference-result outflow from
  flexyBayes, with a real
  posterior rather than a mock. The returned callable takes a numeric matrix
  (rows = parameter draws, columns = parameters on the natural / constrained
  scale) and returns one unnormalised `log p(theta | data)` per row; it
  carries `parameter_names`, an `NA` `log_normalizer` (a posterior's marginal
  likelihood is generally unknown -- reported as unknown, not fabricated), the
  parameters' `support_lower` / `support_upper` bounds, and the fit's
  posterior `draws` to seed the consumer's proposal.

  The **greta** backend is the canonical real producer: it evaluates the
  retained model graph's unadjusted joint density at the free-state image of
  the supplied natural-scale parameters, which is the unnormalised
  natural-scale log-posterior exactly (validated against an analytic
  conjugate posterior to machine precision -- correlation 1, constant
  offset). Out-of-support rows return `-Inf` rather than erroring. The
  **brms** and **INLA** backends abstain with an informative, classed
  `fb_c4_unavailable` condition: brms's log-density lives on the Stan
  unconstrained scale with a version-fragile name mapping, and INLA's
  posterior is a deterministic Laplace / grid approximation rather than a
  sampling log-density, so a plain abstain is preferred to a
  plausible-but-wrong producer (the Independent Oracle Principle). flexyBayes
  does not depend on proxymix -- the cross-package demonstration lives in a
  separate integration harness, preserving the acyclic dependency graph.

# flexyBayes 0.8.1

A hub-ergonomics release: a single broom dialect across the three backends,
two standalone maximum-likelihood fitters for distributions outside the GLM
emit path, and an INLA-led multi-environment-trial vignette.

## New features

* **`tidy()` / `glance()` / `augment()` -- one broom dialect across the
  hub.** The broom-style `tidy()` generic (re-exported from `generics`, now a
  lightweight `Imports`) returns a flat one-row-per-term `data.frame` with the
  canonical `term` / `estimate` / `std.error` / `conf.low` / `conf.high`
  columns for the greta (`flexybayes`), brms (`flexybayes_brms`), and INLA
  (`flexybayes_inla`) fit classes alike -- the INLA class gains its own
  `tidy.flexybayes_inla()` method (it does not inherit from `flexybayes`).
  `glance()` and `augment()` cover the greta and brms classes; on an INLA fit
  they error (turned into an informative refusal in 0.8.3). A cross-engine
  comparison table is now an `rbind()` of two `tidy()` outputs rather than a
  hand-built reconciliation of three different backend layouts. Mirrors the
  `tidy()` method kernR adopted, so the orchestra speaks one tidy dialect.
* **`fb_gev()` -- generalised extreme value (block-maxima) fitter.** Fits
  the location, scale, and shape of a GEV distribution to block maxima
  (annual maximum rainfall, peak yields) by dependency-free maximum
  likelihood, and reports return levels for the requested return periods.
  The family descriptor `fb_family_gev()` and the simulator `rgev()` ship
  alongside. `family = "gen_extreme_value"` in `flexybayes()` now routes to
  `fb_gev()` with an explicit pointer rather than a generic refusal -- block
  maxima have no GLM mean-link, so they do not belong on the formula emit
  path. (A scalable Bayesian GEV on INLA's native `gev` family is planned.)
* **`fb_dirichlet()` -- compositional (simplex) Dirichlet fitter.** Fits the
  concentration vector of a Dirichlet distribution to compositional rows
  (soil texture fractions, species abundance, allele frequencies) by
  maximum likelihood (the default, dependency-free) or via greta's native
  `dirichlet` distribution (`method = "greta"`), and reports the fitted mean
  composition. The family descriptor `fb_family_dirichlet()` and the
  simulator `rdirichlet()` ship alongside; `family = "dirichlet"` in
  `flexybayes()` routes here.

## Minor improvements and fixes

* The MET and genomics vignette now leads with the scalable INLA MET path
  (the diagonal genotype-by-environment model fits in seconds with
  trustworthy posteriors) as the recommended route, and presents the greta
  factor-analytic route as the slower, harder-mixing alternative for the
  stability decomposition.
* `generics` moves from `Suggests` to `Imports`, so `tidy()` / `glance()` /
  `augment()` are available without attaching `broom`; the runtime
  `registerS3method()` shim in `.onLoad()` is replaced by static `S3method()`
  registration.

# flexyBayes 0.8.0

## Genomics and MET expansion

* **`triangulate_genomic()` / `triangulate_gwas()` -- cross-engine and
  field-standard genomic triangulation.** `triangulate_genomic()` compares
  two GBLUP / pedigree analyses on heritability, the variance components,
  and the breeding values (matched by genotype) -- either two flexyBayes
  fits (does greta agree with brms and INLA?) or a flexyBayes fit against a
  generic *genomic lens* (`list(h2, var_g, var_e, gebv)`), the form a
  field-standard REML answer from sommer supplies. `triangulate_gwas()`
  compares two genome scans by the agreement that matters -- the Jaccard
  overlap of the significant marker sets, the top-marker overlap, and the
  effect correlation. Both carry the same shared-upstream caveat as
  `triangulate()` (agreement is not correspondence). flexyBayes core never
  depends on the field tools; the lens form lets the companion build the
  koine fourth opinion. Breeding-value labels are now the genotype factor
  levels on every backend (previously the greta and INLA paths used
  positional labels), so GEBVs match across engines.
* **`fb_met_summary()` -- breeder summary of a factor-analytic MET fit.**
  For a `fa(env, k):gen` factor-analytic genotype-by-environment fit it
  reports the quantities a plant breeder acts on: each genotype's overall
  performance (the across-environment mean of its realised effects) and
  stability (the across-environment spread), the genotype-by-environment
  BLUPs, the factor loadings, and the environment genetic-correlation
  matrix -- the crossover structure (negative correlations are rank
  reversals across environments). These come from the *realised* effects,
  which are identified (invariant to the loadings' rotation and sign
  ambiguity), so their posteriors are interpretable. The greta
  factor-analytic fit now monitors those realised genotype-by-environment
  effects (it previously monitored only the loadings and specific
  variances).
* **`fb_gblup_cv()` -- genomic-prediction accuracy by cross-validation.**
  The payoff layer of genomic selection: how well does a GBLUP trained on
  phenotyped genotypes predict the held-out performance of genotypes seen
  only through their markers? Each fold estimates the variance components
  by REML on the training set and predicts the held-out breeding values
  from the relationship matrix (the exact GBLUP prediction equation,
  evaluated through the spectral primitive), then reports prediction
  accuracy (predicted-observed correlation), bias (the observed-on-
  predicted slope), and the per-fold spread, with repeated-CV averaging.
  Validated against an exact full-matrix prediction; accuracy is near zero
  for an unrelated / non-heritable trait, rises with heritability, and is
  substantial when genotypes are related -- the realistic selection
  setting.
* **`fb_gwas()` -- genome-wide association scan.** A whole-genome EMMAX /
  P3D scan (Kang et al. 2010): the polygenic null mixed model is fit once
  by REML to estimate the variance components, then every marker is tested
  by exact generalised least squares under those fixed components. The
  shared spectral primitive turns the per-marker test into an `O(n)`
  weighted least squares after a single eigendecomposition, so the scan is
  feasible without a per-marker model fit, and it needs no MCMC backend --
  it is a deterministic frequentist fast path (the backends enter only at
  optional top-hit refinement). Returns marker effects, standard errors,
  chi-square statistics, p-values, Bonferroni and Benjamini-Hochberg FDR
  adjustments, the genomic-control inflation factor, and the null REML
  heritability, with `print()` and `plot()` (Manhattan / QQ) methods. The
  REML variance components are validated against an exact full-matrix
  reference and against \pkg{sommer}'s independent REML; the per-marker
  statistic is validated against exact per-marker GLS.
* **Genomic BLUP is now three-engine triangulatable.** The genomic /
  pedigree relationship random effect `vm(geno, G)` / `ped(animal, A)`
  reaches all three backends: greta and brms via the dense relationship
  matrix (brms's native `gr(geno, cov = G)` group term, which Cholesky-
  factors the covariance internally), and INLA via the precision carrier
  `vm(geno, cov = fb_cov(solve(G), type = "precision"))` (the `generic0`
  sparse-precision path -- a dense GBLUP precision is INLA-feasible but
  forgoes INLA's sparsity advantage, so it is opt-in by carrier rather
  than silent). A simulated heritability is recovered on all three
  engines with mutually-close posterior means, so `triangulate()` can
  cross-check a GBLUP fit across paradigms.
* **`genomic_summary()`** extracts the breeder-facing quantities from a
  fitted relationship model on any backend: narrow-sense heritability
  \eqn{h^2}, genomic estimated breeding values (GEBVs) with posterior
  reliability, and the genetic / residual variances. The greta GBLUP fit
  now also monitors the breeding-value vector, so GEBVs are available on
  greta, brms, and INLA alike (previously the greta path reported only
  the variance components).
* **Spectral efficiency primitive (foundations).** A shared internal
  eigendecomposition primitive now underpins the genomics / MET work: it
  decomposes a relationship matrix \eqn{K = U \Lambda U^\top} once and exposes
  the rotation that turns a structured genetic random effect into an
  independent one. This is the machinery a genome-wide scan reuses across every
  marker (the rotated model has a diagonal residual covariance, so each marker
  is an `O(n)` score test rather than a fresh mixed-model fit), and that scales
  genomic BLUP and variance-component estimation. It is positive-semidefinite
  aware -- numerical-noise negative eigenvalues are clamped and reported, while
  genuine indefiniteness is refused -- and validated against full-covariance
  generalised least squares to machine tolerance. Internal at this stage; the
  fit routes that consume it follow.
* **Genomic output contract (foundations).** A standardised genomic summary
  (narrow-sense heritability \eqn{h^2}, genomic estimated breeding values with
  posterior reliability, and -- for whole-genome marker models -- marker
  effects with posterior retention probabilities) is now computed engine-
  agnostically from posterior draws, so the greta / INLA / brms paths all feed
  the same triangulatable result. The fit-level accessor lands with the
  multi-backend GBLUP route.

## Triangulation accuracy

* `triangulate()` gained a `data_independence` argument and a
  `shared_upstream_caveat` result field (Independent Oracle Principle). It
  measures inter-fit *agreement*, and the backend-independence registry
  certifies code (not data) independence -- so if both fits share the same
  upstream data, a fabricated data fact is common-mode and their agreement
  cannot detect it. Unless the caller declares `data_independence = TRUE`, the
  result carries a caveat (surfaced prominently by `print()`) that agreement
  does not test a shared upstream data fact. Pure metadata; the metrics are
  unchanged.

## Engine reliability and diagnostics

* **Faithful default prior on the INLA backend.** The package's default
  `uniform(0, scale)` prior on every variance-component standard deviation (and
  any user-supplied `uniform()` / `half_normal()` prior) is now represented
  *exactly* on the INLA backend via an expression-prior on the log-precision,
  replacing the former PC-prior approximation. The PC approximation concentrated
  prior mass at zero and so shrank a small-group variance component more than
  the greta backend's flat uniform did, producing a cross-engine prior mismatch
  that surfaced as spurious `triangulate()` disagreement on the variance
  components for models with few groups. The two engines now carry genuinely the
  same default prior, so their variance-component posteriors agree far more
  closely. This changes the default INLA variance-component posterior for
  random-effects models fit without an explicit prior.
* **Convergence warning.** MCMC fits (`greta`, `brms`) now emit a warning when
  the sampler may not have converged (a parameter with Rhat at or above 1.1, or
  a low effective sample size), rather than surfacing it only as a print-method
  badge. Silence it for intentionally short fits with
  `options(flexyBayes.silence_convergence_warning = TRUE)`. The deterministic
  INLA path carries no such warning.
* **New `fb_structured_cov()`** reports the *identified* covariance
  \eqn{G = \Lambda\Lambda^\top + \mathrm{diag}(\psi)} for factor-analytic
  `fa()` terms, with an entrywise Rhat. The raw loadings are identified only up
  to rotation and sign, so their per-entry Rhat is meaningless; \eqn{G} is
  rotation- and sign-invariant and is the quantity whose convergence is
  interpretable. The convergence warning points to it when a factor-analytic /
  unstructured term is present.
* **Fixed** an intercept-only model on the aggregated INLA path
  (`y ~ 1 + (1 | g)`, `backend = "inla"`, `aggregate = TRUE`), which crashed in
  the fixed-effect summariser because a length-one variance vector was misread
  by `diag()` as a dimension. The one-coefficient covariance is now built
  correctly.

## Usability

* **New `fb_backend_status()`** reports which inference backends are installed
  and usable in the current session (greta additionally needs a reachable
  Python / TensorFlow stack), with an actionable install hint per backend. It
  is read-only and runs without any backend present.

## Lean-core split

* The orchestra-composition layer has been extracted to the companion package
  **flexyBayesOrchestra** (`Imports: flexyBayes`): the surrogate emulators
  (`fit_surrogate()`, `fb_surrogate_ies()`, `fb_surrogate_ppc()`), ensemble
  sources (`fb_ensemble()`, `as_fb_ensemble()`, `read_fb_ensemble()`,
  `verify_fb_ensemble()`, `register_ensemble_source()`), the PESTO
  ensemble-derived priors (`fb_pesto()`, `fb_prior_from_ensemble()`), the
  surrogate conformers (`register_pesto_surrogates()`,
  `register_kernr_surrogates()`, `register_surrogate()`), and the dormant koine
  fourth-opinion slot (`koine_status()`). Install `flexyBayesOrchestra` for
  these; `flexyBayes` alone now ships the mixed-model and cross-engine
  triangulation core.
* Dropped `kernR`, `PESTO`, `S7`, and `gretaR` from `Suggests` and removed the
  `Remotes:` field, so the package installs cleanly from CRAN-style
  repositories with no GitHub- or r-universe-only dependencies. The `greta`
  (CRAN) and `INLA` backends remain; `INLA` is served by its own
  `Additional_repositories`. The `gretaR` R-native engine is still wired as a
  dormant, opt-in backend: install `gretaR` yourself and it is detected at run
  time, but it is no longer a declared dependency of the public core.

## Methodology

* **Identified factor-analytic loadings.** `fa(x, k)` now emits the loadings
  matrix `Lambda` with the standard lower-triangular, positive-diagonal
  identification (Lopes & West 2004): the strict upper triangle is zeroed and the
  diagonal is constrained positive. Previously `Lambda` was a free `normal(0, 1)`
  matrix, which is unidentified up to rotation, sign, and column permutation for
  `k > 1`, so its posterior summaries were not interpretable. Predictions and the
  g-side covariance are unchanged -- the product `F %*% t(Lambda)` is
  rotation-invariant -- so only the loadings interpretation is affected.
* **Factor-analytic rank upper bound.** `fa(x, k)` is now refused at fit time
  when `k` is not strictly below the number of levels of `x` (new refusal
  `fa_rank_exceeds_dim`). A factor-analytic covariance is identifiable only for
  `k < n_outer`: at `k = n_outer` the loadings and specific variances are an
  over-parameterised reparameterisation of the unstructured form, and at
  `k > n_outer` the lower-triangular loadings carry empty columns. The check is
  data-aware (the number of levels is known only once the term is matched against
  the data), complementing the existing data-free `k >= 1` floor
  (`fa_rank_invalid`). The refusal message points to `us(x)` for a full
  unstructured covariance.
* **Removed the `rhat_means` triangulation metric.** `triangulate()` no longer
  reports `rhat_means`. It pooled two *different* engines' posteriors as if they
  were chains of a single sampler, which conflates genuine between-engine
  approximation bias with within-sampler non-convergence -- so it was not a valid
  convergence diagnostic, and it was mislabelled "rank-normalised R-hat" although
  it applied no rank-normalisation. Cross-engine discrepancy is reported by the
  distributional metrics: `wasserstein_1`, `sd_ratio`, `mean_diff`, and the
  quantile differences. (The per-fit, within-engine rank-normalised R-hat used as
  a convergence diagnostic is unaffected.)
* **Default-prior provenance stated in full.** The default variance-component prior
  (bounded uniform on each SD; `U = 5 * sd(y)` for Gaussian) was attributed
  "following Gelman (2006)". Gelman (2006) in fact recommends a half-t /
  half-Cauchy for variance components with few groups and cautions against a flat
  prior there; the uniform-on-SD default is a weakly-informative choice for
  *moderate* group counts, and the `5 * sd(y)` bound is a flexyBayes heuristic.
  The attribution is corrected throughout (the one-time announcement message,
  `prior_summary()`, the docstrings, the README, and the regression + priors
  vignettes). The default behaviour is unchanged; `fb_prior(half_cauchy(...))`
  remains the documented choice for small `J` (see the priors vignette).

## Documentation

* **Vignettes migrated to the engine-pin API.** Eight vignettes still showed
  the pre-0.5.0 universal `fb_brms(..., backend = "greta" / "inla" / "auto")`
  call surface, which the 0.5.0 engine-pin refactor removed (`fb_brms()` now
  pins to Stan and rejects a conflicting `backend`). Every example is updated to
  the current surface: `fb_greta()` / `fb_inla()` for a fixed engine,
  `flexybayes()` / `fb()` for the universal entry that takes a `backend`
  argument. The brms-style *grammar* is unchanged and still accepted by every
  entry; only the engine selection moved. Stale prose is corrected accordingly
  (the `backend` default is `"auto"`, not `"greta"`; the Stan/brms emit backend
  is live, not "queued"; the native-greta canonical-name map is supplied via
  `fb_from_greta(model, canonical_names = ...)`). Vignettes were re-precompiled.

# flexyBayes 0.7.0

> **Note (current state).** The orchestra-composition features described in this
> section — the ensemble-source data contract, simulator-derived priors, the
> surrogate emulators and predictive checks, and the dormant fourth-opinion slot
> — were subsequently extracted to the companion package **flexyBayesOrchestra**
> (see the lean-core split under 0.8.0 below). They are no longer part of the
> `flexyBayes` exported surface; the entries below record the 0.7.0 history.

This release reconciles three parallel development streams onto the engine, data,
and surrogate axes: gretaR activated as a fourth inference engine (the engine
axis), the ensemble-source data contract and its first consumers (the data
axis), and the surrogate carried end-to-end from a reference emulator to a
distribution-preserving predictive check (the surrogate axis).

## New features

- **gretaR activated as a fourth inference engine.** `flexybayes(..., backend =
  "gretaR")` fits supported hierarchical models through an out-of-process
  torch-NUTS worker, returning canonical parameter names and a `draws_array` that
  flows into `triangulate()` like any other engine. The integration is governed by
  a versioned backend contract with an **executable
  conformance battery** -- a new engine onboards by adding one descriptor, gated on
  `triangulate()`-agreement within an SBC-calibrated threshold. A **dormant koine
  backend slot** is provisioned (scaffolded and correct, switched on when koine's
  programmatic model builder is confirmed against the contract).
- **`fb_ensemble()` and the `ensemble_source` data contract.**
  A canonical ingest shape for calibrated parameter ensembles, with PESTO (via its
  manifest), raw `apsimx`, and file as conforming producers -- so flexyBayes is
  standalone-functional for methods work without depending on any one producer's
  internals. **`fb_pesto()`** turns a calibrated ensemble into an informative
  `fb_prior` for a downstream fit (moments method; joint and KDE methods refuse
  cleanly behind their gates).
- **`fitted()`, `residuals()`, and `logLik()` methods for INLA fits.** `fitted()`
  and `residuals()` now return values for a `flexybayes_inla` fit (previously
  `NULL`); `logLik()` returns `NA` with the right shape where INLA does
  not expose a likelihood, rather than erroring.
- **`fb_surrogate_ppc(whiten = TRUE)` whitens the outputs before the MMD.** The
  predictive draws and `observed` are Mahalanobis-whitened by their pooled
  covariance -- a label-agnostic representation that puts the outputs on a common
  footing, removing the cross-output scale heterogeneity that makes a single
  median bandwidth too coarse to resolve a distortion in a low-variance direction.
  The map depends only on the pooled cloud (not the draws-vs-observed labels), so
  the permutation test stays valid, and the same map is applied to both clouds, so
  a relative scale, shape, or mean difference is preserved. It lifts power
  substantially on localised / heterogeneous-scale distortions a raw-output check
  misses (with no blind spot found across a mean / variance / shape battery), at a
  valid level. Opt-in (a graduated ridge guards a near-singular pooled
  covariance); composes with `project` and `aggregate`.
- **`fb_surrogate_ppc(aggregate = TRUE)` runs an aggregated multi-kernel MMD.**
  An MMDAgg test (Schrab et al. 2023) over an RBF bandwidth grid around the median
  heuristic, sharing the permutations across bandwidths and aggregating by the
  weighted-quantile statistic `min_g p_g / w_g` -- it removes the bandwidth choice
  and is robust when the median heuristic is mis-scaled, at a valid level. The
  kernel weights are set by `agg_weights` (`"uniform"` default, or `"increasing"`
  / `"decreasing"` / `"centred"`, or a numeric vector): up-weighting the scale a
  distortion is expected to live at recovers power, and `"uniform"` recovers the
  min p-value aggregation. It orchestrates kernR's single-kernel `mmd_test` over
  the grid (the MMD machinery stays kernR's) and composes with `project`. Returns
  an `fb_mmd_agg` object.
- **New vignette: distribution-preserving surrogates (emulate, invert, check).**
  A self-contained walkthrough of the surrogate workflow -- ingest an
  ensemble, fit a surrogate, read its predictive distribution, invert it with
  `fb_surrogate_ies()`, and check it with `fb_surrogate_ppc()` -- runnable on the
  built-in reference surrogate with no other package.
- **`fb_surrogate_ppc(project = TRUE)` runs a factor-space MMD for large output
  counts.** For a surrogate that stores a low-rank covariance (`cov_rank`), the
  predictive-check draws and `observed` are projected onto the surrogate's `q`
  factor directions before the MMD, so the test runs in `q` dimensions, not `k`
  -- cheaper and higher-powered on the correlated cross-output structure (the MMD
  loses power as the dimension grows). It is, by design, blind to distortions
  orthogonal to the factor space (idiosyncratic per-output errors, which the
  marginal and `joint = FALSE` checks cover), so it complements the full-space
  check.
- **`fit_surrogate(cov_method = "auto")` picks the low-rank factor method.**
  In addition to `"eigen"` (default) and `"fa"`, `"auto"` selects between them
  from the heteroscedasticity of the eigen-residual diagonal -- `"fa"` when the
  residual is uneven (where eigen-truncation is biased), `"eigen"` otherwise. The
  resolved method is recorded in the fitted object.
- **`fit_surrogate(ard_refit = "cv")` now uses k-fold cross-validation.** The
  held-out predictive-density objective is summed over a k-fold partition of the
  subsample rather than a single split -- a less noisy objective, so the
  Nelder-Mead lands on sharper lengthscales.
- **`fb_surrogate_ppc()` samples a low-rank covariance factored.** When the
  surrogate stores a low-rank covariance (`cov_rank`), the joint predictive-check
  draws are sampled from the factored form directly -- a `q`-dimensional factor
  noise plus a `k`-dimensional idiosyncratic noise -- so the `k x k` covariance
  slice is never formed (the native low-rank consumer for the PPC, `O(k q)` per
  draw). Draws from the same per-row covariance as the dense path.
- **`fb_surrogate_ies()` consumes a low-rank covariance natively (Woodbury).**
  When the surrogate stores a low-rank covariance (`cov_rank`), the ES-MDA Kalman
  update folds the factored form straight in via the Woodbury identity (a
  rank-sized solve) and samples the observation perturbation factored, so the
  `d x d` data covariance is never formed -- the inversion cost scales with the
  covariance rank, not the output count. Produces the same Kalman gain as the
  dense path (a dense-joint or marginal surrogate keeps the dense path).
- **`fit_surrogate(cov_method = "fa")` learns the low-rank factor by EM.**
  In addition to the default eigen-truncation (`"eigen"`), the rank-`cov_rank`
  factor can be fitted by factor-analysis EM, which recovers the loadings more
  faithfully when the residual diagonal is heteroscedastic (the eigenvectors of
  the full covariance are otherwise contaminated by the uneven diagonal). The
  covariance diagonal stays exact either way.
- **`fit_surrogate(ard_refit = "cv")` regularises the fuller evidence-ARD.**
  The `ard_refit = TRUE` profiling maximises the in-sample evidence and can
  over-shorten lengthscales; `"cv"` instead maximises the held-out predictive
  log-density on a train/validation split, keeping the relevance sharpening
  without the in-sample over-shortening.
- **`fb_surrogate_ppc()` gains a joint draw path.** With
  `joint = TRUE` (default) the predictive-check draws are sampled from the full
  cross-output predictive covariance (`predict(cov = TRUE)` + a per-row
  Cholesky), so the emulator's cross-output correlation is exercised -- the
  second consumer of the joint covariance after `fb_surrogate_ies()`. On a
  correlated-output simulator, marginal draws (`joint = FALSE`) discard that
  structure and false-reject a faithful surrogate, so the joint path is the
  correct distribution-preserving check. A single output, or a marginals-only
  surrogate, behaves exactly as before.
- **`fit_surrogate(method = "rff")` gains `cov_rank` for large output counts.**
  Stores the cross-output noise covariance in a rank-`q` factor-analysis form
  (`U U' + diag`), `O(k q)` instead of `O(k^2)`. The per-output marginal sd and
  the covariance diagonal are reproduced exactly; only the off-diagonal
  covariance is approximated. `predict(<fb_surrogate>, newdata,
  cov = "lowrank")` returns the compact factored predictive
  (`cov_lowrank = list(u, d, infl)`), and the new `fb_cov_slice()` reconstructs
  a `k x k` covariance slice from either the dense or the factored form.
- **`fit_surrogate(method = "rff")` gains `ard_refit` (fuller evidence-ARD).**
  When `TRUE`, the prior/noise precisions `(a, t)` are re-estimated by evidence
  at each candidate lengthscale (the evidence profiled over `(a, t)`), a fuller
  -- and at equal lengthscales never lower -- marginal likelihood. It sharpens
  relevance ranking but can over-shorten lengthscales (a higher evidence does
  not guarantee better held-out recovery for finite-feature models), so it is
  opt-in (`ard_refit = FALSE` default).
- **The reference `rff` surrogate gains ARD lengthscales and a joint
  multi-output fit, both on by default.** `fit_surrogate(method = "rff")`
  now learns per-input (automatic relevance determination) lengthscales by
  maximising the model evidence (`ard = TRUE`), so irrelevant parameters are
  stretched out of the kernel; and it fits a joint matrix-normal --
  inverse-Wishart model across outputs (`joint = TRUE`), preserving the
  cross-output predictive covariance rather than only the marginals. Both
  default behaviours can be turned off (`ard = FALSE`, `joint = FALSE`) to
  recover the isotropic, per-output fit. The per-output `<o>_sd` contract is
  unchanged.
- **`predict(<fb_surrogate>, newdata, cov = TRUE)` returns the joint
  predictive distribution**: an
  `fb_surrogate_prediction` with the predictive mean matrix and a per-row
  cross-output covariance array. A marginals-only surrogate (independent
  `rff`, or a conformer reporting only `<o>_sd`) returns a diagonal
  covariance, so the call always states what it used.
- **`fb_surrogate_ies()` consumes the joint predictive covariance.** When the
  surrogate is joint, the full cross-output predictive covariance enters the
  ES-MDA data covariance (`use_surrogate_cov = TRUE`, default), so correlated
  emulator errors are assimilated coherently; a marginals-only surrogate
  falls back to the diagonal path, identical to before.
- **Sibling surrogate backends register as conformers.** PESTO's GP
  surrogate (`method = "gp"`) and kernR's conditional mean embedding
  (`method = "cme"`) are available to `fit_surrogate()` when the member
  package is installed, registered through the existing surrogate
  registry. PESTO's GP already reports a per-output predictive mean and
  variance; kernR's CME is a kernel linear smoother, so the adapter
  supplies the standard linear-smoother predictive variance
  `sigma^2 (1 + ||w||^2)` to honour the distribution-preserving `<o>_sd`
  requirement. Both register opportunistically on load, or explicitly via
  `register_pesto_surrogates()` / `register_kernr_surrogates()`. The
  surrogate methods are consumed through the contract, never copied into
  flexyBayes, and the built-in `rff` reference keeps the package
  standalone.
- **`fb_surrogate_ies()` calibrates a parameter ensemble through a
  surrogate.** An Ensemble Smoother with Multiple Data Assimilation
  (ES-MDA) drives a prior parameter ensemble toward observed data using
  the surrogate in place of the forward model; the predictive `<o>_sd` is
  folded into the data covariance, so emulator uncertainty inflates the
  assimilation noise.
- **`fb_surrogate_ppc()` checks whether a surrogate distorted the output
  distribution.** It draws from the surrogate's predictive distribution
  at held-out parameter sets and runs kernR's `mmd_ppc()` against the true
  simulator outputs there -- a distribution-level check, not a mean-only
  accuracy measure.

## Bug fixes

- **Registered surrogate conformers now predict on the standardised
  design.** `predict.fb_surrogate()` previously passed raw `newdata` to a
  registered conformer's `predict()` while training it on the standardised
  design, so train and predict disagreed on scale. The standardised design
  is now passed consistently.

---

Older release history (0.6.0 and earlier) lives in `NEWS_ARCHIVE.md`.
