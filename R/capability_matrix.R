# capability_matrix.R -- the single source of truth for the per-engine
# capability table shown on every public surface.
#
# One R-level table, rendered to Markdown by
# tools/generate_capability_matrix.R and spliced between the
# <!-- capability-matrix:begin --> / <!-- capability-matrix:end -->
# markers in README.md and inst/KNOWN_ISSUES.md. A staleness test
# (tests/testthat/test-capability-matrix.R) fails when a committed block
# drifts from a fresh render, and a behaviour anchor per row asserts the
# verdict against the emit / gate code rather than against this table.
#
# The verdict vocabulary is closed:
#
#   fits     an active engine emits the structure and a test exercises it.
#   emits    the engine generates the structure, and no live fit has yet
#            confirmed it samples acceptably. Used where the emit is
#            asserted but the inference is not.
#   refuses  the request raises rather than fitting something else.
#   n/a      the class does not apply to that engine's interface.
#
# Every row carries the anchor field naming the test that grounds it, so a
# claim on a public surface can be traced to the assertion that holds it.

# --- the table ----------------------------------------------------- #

#' Per-engine capability table
#'
#' The declared capability of each active engine by model class. This is
#' the only place the capability verdicts are written down; README.md and
#' `inst/KNOWN_ISSUES.md` carry a generated copy, and
#' `tests/testthat/test-capability-matrix.R` re-derives every verdict from
#' the gate and emit code.
#'
#' @returns A data.frame with one row per model class and the columns
#'   `model_class`, `spelling`, `inla`, `brms`, `note`, and `anchor`.
#'
#' @noRd
#' @keywords internal
.fb_capability_matrix <- function() {
  rows <- list(
    list(
      model_class = "Gaussian LMM, simple random intercept",
      spelling = "`random = ~ g` / `(1 | g)`",
      inla = "fits",
      brms = "fits",
      note = "The class both engines emit, and the one their agreement is
        checked on in the test suite.",
      anchor = "test-backend-conformance-open-oracle.R"
    ),
    list(
      model_class = paste(
        "GLMM (binomial, Poisson, negative binomial, gamma, beta),",
        "simple random effect"
      ),
      spelling = "`(1 | g)` with `family =`",
      inla = "fits",
      brms = "fits",
      note = "INLA's likelihood allowlist is read from
        `INLA::inla.models()` when INLA is installed.",
      anchor = "test-family-support.R"
    ),
    list(
      model_class = "Hurdle gamma (zero mass plus a positive gamma part)",
      spelling = "`family = \"hurdle_gamma\"`",
      inla = "refuses",
      brms = "fits",
      note = "brms-native (`dpars` mu, shape, hu); the zero-mass
        probability `hu` keeps brms's own prior. INLA's likelihood roster
        carries no counterpart, so the family gate refuses it there and
        `auto` routes to brms.",
      anchor = "test-family-support.R"
    ),
    list(
      model_class = "Uncorrelated random slope",
      spelling = "`(x || g)`",
      inla = "fits",
      brms = "fits",
      note = "Two independent iid fields on INLA, one for the intercept and
        one for the slope, checked against lme4's REML fit of the same
        model. Reported as `sd_<g>` and `sd_<x>_<g>`.",
      anchor = "test-inla-slope-and-interaction-verified.R"
    ),
    list(
      model_class = "Factor-by-numeric fixed interaction",
      spelling = "`y ~ f * x` with numeric `x`",
      inla = "fits",
      brms = "fits",
      note = "An ordinary treatment-coded design matrix on both engines.
        The INLA posterior means sit on the least-squares estimates of the
        same model under a weak prior.",
      anchor = "test-inla-slope-and-interaction-verified.R"
    ),
    list(
      model_class = "Correlated random slope",
      spelling = "`(x | g)`",
      inla = "fits",
      brms = "fits",
      note = "The uncorrelated slope with the correlation freed: one
        intercept SD, one slope SD and the correlation between them. INLA
        holds all three in a single iid2d field over a two-block index,
        with the slope block entering as a fixed-scaling copy weighted by
        `x`; brms writes `(x | g)` verbatim. Both report `sd_<g>`,
        `sd_<x>_<g>` and `cor_<g>__Intercept__<x>`. The iid2d
        hyperparameters are one joint Wishart with no standard-deviation
        slot, so on INLA an explicit SD prior refuses by name and the
        package default is reported as an engine default rather than as a
        prior that reached the fit. The joint Wishart's scale matrix is set
        from the response's own variance, and the slope block from the
        response variance per unit of the slope variable. Leaving INLA's
        default there -- a scale matrix of ones, so a prior centred on unit
        variances -- makes the mode non-unique on any response that is not
        near unit scale: on `lme4::sleepstudy` eight identical calls
        returned intercept SDs from 0.5 to 54.8 and correlations from -1.00
        to +0.75. Scaled, the same fit returns 25.7 to 26.9 against lme4's
        24.7, with the correlation between -0.07 and +0.15 around lme4's
        0.066.",
      anchor = "test-correlated-random-slope.R"
    ),
    list(
      model_class = "Nested / interaction random effects, multi-stratum",
      spelling = "`~ gen:env`, `~ env:rep:block`",
      inla = "fits",
      brms = "fits",
      note = "An iid field over the interaction index on INLA, `(1 | a:b)`
        on brms -- the same model, and the two agree with lme4's REML fit
        of it on `agridat::john.alpha`. Only observed cells get a level,
        as in lme4 and brms.",
      anchor = "test-nested-random-inla.R"
    ),
    list(
      model_class = "Heterogeneous variance by factor level",
      spelling = "`~ diag(f):g`, `~ idh(f):g`, `~ at(f):g`",
      inla = "fits",
      brms = "fits",
      note = "One variance per level of `f`, no covariance between levels.
        All three spellings emit identical code. INLA writes k independent
        iid fields, one per level; brms writes `(0 + f || g)`. Both report
        the component as `sd_<g>__<f><level>` and agree with lme4's REML
        fit of the same model. `at(f, level)` conditions on selected levels
        -- a different model -- and still refuses.",
      anchor = "test-heterogeneous-variance-inla.R"
    ),
    list(
      model_class = "Unstructured genotype-by-environment covariance",
      spelling = "`~ us(f):g`",
      inla = "fits",
      brms = "fits",
      note = "One iidKd field on INLA (k = 2..5; wider refuses and names
        fa()), `(0 + f | g)` on brms. **On INLA this route does not
        reproduce ASReml's covariance.** Paired against ASReml 4.2.0.392 on
        `agridat::besag.met` it returns a genetic correlation of 0.142 where
        ASReml returns 0.930, and 0.198, 0.231 and 0.320 where ASReml returns
        0.898, 0.915 and 0.955 over three environments, with genotype
        variances four to five times larger while the residual agrees
        closely. The gap reproduces across refits, so it is a stable wrong
        answer rather than a noisy one. Rescaling the response by ten is an
        exact transformation and must not decide whether the model fits, and
        here it does: across eight identical rescaled calls none converged and
        the returned correlation spanned -0.894 to +0.686. Read the
        correlations from the brms route, which was not exercised in that
        comparison, or from ASReml. An explicit SD prior refuses on INLA because iidKd has no
        slot for one. The correlated sibling of
        the diagonal structure --
        `k(k+1)/2` parameters against `diag()`'s `k`. At one observation
        per cell the residual variance is confounded with the diagonal
        of the covariance: the covariance block converges and `sigma`
        does not, and a longer chain does not help. Replicate within
        cell, or put an informative prior on the residual.",
      anchor = "test-heterogeneous-variance.R"
    ),
    list(
      model_class = "Heterogeneous variances with one shared correlation",
      spelling = "`~ corh(f):g`",
      inla = "refuses",
      brms = "refuses",
      note = "No active engine has an equicorrelation group-level
        structure. Use `diag(f):g` or `us(f):g`.",
      anchor = "test-heterogeneous-variance.R"
    ),
    list(
      model_class = "Heterogeneous residual by factor level",
      spelling = "`residual = ~ dsum(~ units | f)` / `~ at(f):units`",
      inla = "fits",
      brms = "fits",
      note = "One residual variance per section. INLA carries it through the
        multiple-likelihood interface, a matrix response with one Gaussian
        family per level; brms lowers it to distributional regression on
        log sigma, `sigma ~ 0 + f`. Both report `sigma_<level>` on the
        standard-deviation scale. Refused for families with no residual
        scale; `at(f, level):units` conditions on selected levels, which is
        a different model, and refuses.",
      anchor = "test-sectioned-residual-inla.R"
    ),
    list(
      model_class = paste(
        "Combined interaction random effects and heterogeneous residual",
        "(full MET)"
      ),
      spelling = "`random = ~ gen + gen:env` with the `dsum` residual",
      inla = "fits",
      brms = "fits",
      note = "The full multi-environment-trial model. Both halves refused on
        INLA before 0.10.0.9000, so the combination was unreachable there.
        It now recovers the genotype, genotype-by-environment and all three
        per-section residual standard deviations on simulated data.",
      anchor = "test-sectioned-residual-inla.R"
    ),
    list(
      model_class = "Factor-analytic genotype-by-environment covariance",
      spelling = "`~ fa(env, k):gen`",
      inla = "fits",
      brms = "refuses",
      note = "Each common factor is one latent field entered per level of `f`
        as a SCALED COPY, and INLA's copy scaling IS the loading, so the
        loadings carry a posterior rather than a point estimate. Factor j is
        anchored at level j, the zero upper triangle ASReml uses against
        rotation, so an anchored loading is absent from the summary rather
        than reported as 1. Verified at k = 1 and k = 2. brms group-level
        effects are either uncorrelated or fully unstructured, with no
        factor-analytic form between.",
      anchor = "test-capability-matrix.R"
    ),
    list(
      model_class = "Multi-trait covariance",
      spelling = "`~ us(trait):vm(gen)`",
      inla = "refuses",
      brms = "refuses",
      note = "No active engine represents a trait-by-genotype
        unstructured covariance.",
      anchor = "test-capability-matrix.R"
    ),
    list(
      model_class = "Known-covariance genomic / pedigree random effect",
      spelling = "`~ vm(g, K)`, `~ ped(a, A)`",
      inla = "fits",
      brms = "fits",
      note = "INLA takes the sparse-precision, pedigree-precision and
        block carriers, and brms additionally takes dense and Cholesky.
        Under `fb_select()`: shrinkage yes, projection no. A joint
        shrinkage prior over candidate fixed effects fits on a base
        carrying one of these terms, but the term lowers to
        `(1 | gr(g, cov = K))` and projpred refuses a function call on the
        right-hand side of a group-term `|` (probe C5), so the projection
        decision is red. `decision = \"auto\"`, the default, resolves to
        `full` there; naming `\"projection\"` refuses before the reference
        model is sampled.",
      anchor = "test-known-covariance-inputs.R"
    ),
    list(
      model_class = "Separable AR1 spatial field",
      spelling = "`random = ~ ar1(row):ar1(col)`, `random = ~ ar1(t)`",
      inla = "fits",
      brms = "refuses",
      note = "A latent AR1 field plus the Gaussian observation nugget --
        four hyperparameters, one observation per grid node. This is not
        ASReml's three-parameter nugget-free residual, so the residual
        spelling refuses and names this one.",
      anchor = "test-inla-spatial-ar1.R"
    ),
    list(
      model_class = "Per-trial separable AR1 field",
      spelling = "`random = ~ at(trial):ar1(row):ar1(col)`",
      inla = "fits",
      brms = "refuses",
      note = "One field realisation per level of `trial`, via INLA's
        `replicate =` mechanism, but the row correlation, column
        correlation and field SD are shared across every level -- not
        estimated per trial. `at(trial, level):ar1(row):ar1(col)` (a
        level argument, asking for a single conditioned trial or for
        per-trial hyperparameters) refuses by name
        (`at_field_per_level_hyper_not_representable`). brms has no
        lowering for either spelling.",
      anchor = "test-nested-ar1-field.R"
    ),
    list(
      model_class = "Grouped temporal AR1 field",
      spelling = "`random = ~ ar1(t):g`, `~ g:ar1(t)`, `~ ar1(t):id(g)`",
      inla = "fits",
      brms = "refuses",
      note = "AR1 along `t` within each level of `g` and independent across
        levels: one correlation, one field SD and the Gaussian nugget,
        three hyperparameters shared by every level. The units-by-occasions
        repeated-measures covariance of Verbyla and Cullis (1992), written
        on the random side; an occasion a unit missed stays as a design row
        with the response NA under `na_action = \"augment\"`. Either factor
        order is accepted. brms has no lowering for the grouped
        autoregressive precision and refuses.",
      anchor = "test-grouped-temporal-ar1.R"
    ),
    list(
      model_class = "Univariate P-spline",
      spelling = "`~ spl(x)`",
      inla = "fits",
      brms = "refuses",
      note = "Mapped to INLA's second-order random walk. brms has no
        lowering for the smooth basis. An analogue of ASReml's spline
        rather than a translation of it -- basis, scaling and prediction
        rule differ, and no equivalence has been established (see
        TRANSLATION.md).",
      anchor = "test-smooth.R"
    ),
    list(
      model_class = "Observation weights (Gaussian, identity link)",
      spelling = "`weights = w`",
      inla = "fits",
      brms = "fits",
      note = "Precision weighting, Var(y_i) = sigma^2 / w_i (the ASReml
        / lme4 / glm(weights=) sense): INLA's `scale = w`; on brms a
        known offset on the sigma distributional parameter, NOT brms's
        own `weights()` addition term (a different, likelihood-power
        quantity per brms's own documentation). Both engines match
        lme4::lmer(weights=) closely on a shared simulated fixture. Any
        other family, or a non-identity link on Gaussian, refuses by
        name (`weights_requires_gaussian`); `aggregate = TRUE` alongside
        weights also refuses by name (`weights_not_aggregatable`).",
      anchor = "test-weights-lowering.R"
    ),
    list(
      model_class = "Exact sufficient-statistic aggregation",
      spelling = "`aggregate = TRUE`, `flexybayes_stream()`",
      inla = "fits",
      brms = "n/a",
      note = "Exact cell-likelihood aggregation for iid exponential-family
        models with small cell count. The brms path has no aggregated
        emit.",
      anchor = "test-aggregation-equivalence-backend.R"
    ),
    list(
      model_class = "Ridge shrinkage prior on candidate fixed effects",
      spelling = "`fb_penalty(term(...) ~ ridge(scale))`",
      inla = "refuses",
      brms = "fits",
      note = "A joint shrinkage prior over a declared block of candidate
        fixed effects, compiled to one coefficient-specific `normal` row per
        candidate so the block is targeted exactly and no protected
        coefficient is touched. INLA refuses at the twelfth structural check
        (`prior_family_unsupported_inla`): a scale-mixture prior over the
        fixed effects leaves the latent field non-Gaussian given the
        hyperparameters, which is the assumption the nested Laplace
        approximation rests on. An INLA-fitted base is not blocked by that
        refusal -- `fb_select()` re-emits the reference model on brms with
        the base's structure and priors carried through, and records the
        re-emission on the certificate.",
      anchor = "test-penalty-stancode.R"
    ),
    list(
      model_class = "Laplace shrinkage prior on candidate fixed effects",
      spelling = "`fb_penalty(term(...) ~ laplace(scale))`",
      inla = "refuses",
      brms = "fits",
      note = "One coefficient-specific `double_exponential` row per
        candidate. The posterior mean is not a lasso estimate and the draws
        carry no exact zeros; the alias `lasso()` says so once per session.
        INLA refuses as for ridge.",
      anchor = "test-penalty-stancode.R"
    ),
    list(
      model_class = "Elastic-net shrinkage prior on candidate fixed effects",
      spelling = "`fb_penalty(term(...) ~ elastic_net(lambda, alpha))`",
      inla = "refuses",
      brms = "fits",
      note = "The canonical elastic-net density written as a `target +=`
        contribution over the candidate indices through `stanvars` and
        `check = FALSE`, so the block is targeted exactly while a protected
        coefficient keeps its own prior row. Fixed hyperparameters only in
        this release, because the normalising constant is then irrelevant.
        INLA refuses as for ridge.",
      anchor = "test-penalty-stancode.R"
    ),
    list(
      model_class = paste(
        "Regularised horseshoe shrinkage prior on candidate fixed effects"
      ),
      spelling = "`fb_penalty(term(...) ~ regularized_horseshoe(p0))`",
      inla = "refuses",
      brms = "fits",
      note = "brms's native `horseshoe()`, with `par_ratio` set from `p0`
        and calibrated against the number of resampling units rather than
        the number of rows. Exact block targeting ONLY when the base's
        fixed part is intercept-only: it is a class-wide prior, brms
        refuses a coefficient-specific row beside a class-wide special
        prior, so a base carrying protected fixed effects beyond the
        intercept refuses with `penalty_block_targeting_unsupported`. The
        two ways forward are `laplace()`, `ridge()` or `elastic_net()`,
        which target the block exactly, or `scope = \"all_fixed\"`, an
        opt-in under which the protected coefficients also sit under the
        global-local prior, the projection search forces them into every
        submodel, and the certificate records the named approximation
        `protected_under_global_shrinkage`. A candidate block containing a
        design column the re-emitted reference cannot standardise -- a
        factor contrast, an interaction product -- refuses with
        `penalty_candidate_columns_unscaled`. INLA refuses as for ridge.",
      anchor = "test-penalty-stancode.R"
    ),
    list(
      model_class = "R2D2 shrinkage prior on candidate fixed effects",
      spelling = "`fb_penalty(term(...) ~ r2d2(mean_R2, prec_R2, cons_D2))`",
      inla = "refuses",
      brms = "fits",
      note = "brms's native `R2D2()`, whose scales are stated on the
        explained-variance grammar rather than on a coefficient scale, so
        `p0` is refused for it. Class-wide, with the same targeting rule as
        the regularised horseshoe: exact only when the base's fixed part is
        intercept-only, otherwise `laplace()`, `ridge()` or
        `elastic_net()`, or `scope = \"all_fixed\"` with the named
        approximation. Unstandardised candidate columns refuse as for the
        horseshoe. INLA refuses as for ridge.",
      anchor = "test-penalty-stancode.R"
    )
  )

  field <- function(nm) {
    vapply(rows, function(r) .fb_squish(r[[nm]]), character(1L))
  }

  data.frame(
    model_class = field("model_class"),
    spelling = field("spelling"),
    inla = field("inla"),
    brms = field("brms"),
    note = field("note"),
    anchor = field("anchor"),
    stringsAsFactors = FALSE
  )
}

# --- rendering ----------------------------------------------------- #

# The closed verdict vocabulary. Any value outside it is a typo, and the
# renderer refuses rather than printing it onto a public surface.
.FB_CAPABILITY_VERDICTS <- c("fits", "emits", "refuses", "n/a")

# Collapse the source table's wrapped string literals to one line. The
# table is written with line breaks so it stays inside the 80-character
# budget; the rendered cell must be a single Markdown line.
.fb_squish <- function(x) {
  gsub("[[:space:]]+", " ", trimws(x))
}

# Escape the pipe characters that appear inside formula spellings such as
# `(1 | g)` so they do not close a Markdown table cell.
.fb_md_escape <- function(x) {
  gsub("|", "\\|", x, fixed = TRUE)
}

#' Render the capability table as a Markdown block
#'
#' Builds the generated block that `tools/generate_capability_matrix.R`
#' splices into `README.md` and `inst/KNOWN_ISSUES.md`. The block is
#' bounded by the begin / end markers so the generator can replace it
#' without disturbing the surrounding prose.
#'
#' @param markers Logical. When `TRUE` (the default) the returned block is
#'   wrapped in the begin / end HTML comment markers.
#'
#' @returns A single string holding the Markdown block, newline separated.
#'
#' @noRd
#' @keywords internal
.fb_capability_markdown <- function(markers = TRUE) {
  tab <- .fb_capability_matrix()

  bad <- setdiff(c(tab$inla, tab$brms), .FB_CAPABILITY_VERDICTS)
  if (length(bad) > 0L) {
    stop(
      ".fb_capability_markdown(): unknown verdict ",
      paste0("'", bad, "'", collapse = ", "),
      ". The closed vocabulary is ",
      paste(.FB_CAPABILITY_VERDICTS, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  header <- c(
    "| Model class | Spelling | INLA | brms | Notes |",
    "|---|---|:-:|:-:|---|"
  )
  body <- vapply(
    seq_len(nrow(tab)),
    function(i) {
      paste0(
        "| ",
        .fb_md_escape(tab$model_class[[i]]),
        " | ",
        .fb_md_escape(tab$spelling[[i]]),
        " | ",
        tab$inla[[i]],
        " | ",
        tab$brms[[i]],
        " | ",
        .fb_md_escape(tab$note[[i]]),
        " |"
      )
    },
    character(1L)
  )

  legend <- c(
    "",
    paste(
      "`fits` -- the engine emits the structure and a test exercises it.",
      "`emits` -- the engine generates the structure and no live fit has",
      "yet confirmed it samples acceptably. `refuses` -- the request",
      "raises rather than fitting something else. `n/a` -- the class does",
      "not apply to that engine's interface."
    ),
    "",
    paste(
      "This block is generated from `.fb_capability_matrix()` by",
      "`tools/generate_capability_matrix.R`. Edit the R table, re-run the",
      "generator, and let `tests/testthat/test-capability-matrix.R` check",
      "that every verdict still matches the gate and emit code. Do not",
      "edit the rows here by hand."
    )
  )

  block <- c(header, body, legend)
  if (isTRUE(markers)) {
    block <- c(
      .FB_CAPABILITY_MARKER_BEGIN,
      block,
      .FB_CAPABILITY_MARKER_END
    )
  }
  paste(block, collapse = "\n")
}

.FB_CAPABILITY_MARKER_BEGIN <- "<!-- capability-matrix:begin -->"
.FB_CAPABILITY_MARKER_END <- "<!-- capability-matrix:end -->"

# --- splice -------------------------------------------------------- #

#' Splice the generated capability block into a Markdown file
#'
#' Replaces whatever currently sits between the begin and end markers with
#' a fresh render. Errors when either marker is missing or duplicated,
#' rather than appending a second copy of the table.
#'
#' @param path Path to the Markdown file carrying the marker pair.
#' @param write Logical. When `TRUE` the file is rewritten in place; when
#'   `FALSE` the would-be contents are returned without touching the file
#'   (the staleness test uses this).
#'
#' @returns The full file contents as a character vector of lines,
#'   invisibly when `write` is `TRUE`.
#'
#' @noRd
#' @keywords internal
.fb_capability_splice <- function(path, write = TRUE) {
  if (!file.exists(path)) {
    stop(
      ".fb_capability_splice(): no such file '",
      path,
      "'.",
      call. = FALSE
    )
  }
  lines <- readLines(path, warn = FALSE)
  begin <- which(lines == .FB_CAPABILITY_MARKER_BEGIN)
  end <- which(lines == .FB_CAPABILITY_MARKER_END)

  if (length(begin) != 1L || length(end) != 1L || end < begin) {
    stop(
      ".fb_capability_splice(): '",
      path,
      "' must carry exactly one ",
      .FB_CAPABILITY_MARKER_BEGIN,
      " line followed by exactly one ",
      .FB_CAPABILITY_MARKER_END,
      " line. Found ",
      length(begin),
      " begin and ",
      length(end),
      " end marker(s).",
      call. = FALSE
    )
  }

  block <- strsplit(
    .fb_capability_markdown(markers = TRUE),
    "\n",
    fixed = TRUE
  )[[1L]]
  out <- c(
    if (begin > 1L) lines[seq_len(begin - 1L)] else character(0),
    block,
    if (end < length(lines)) {
      lines[seq.int(end + 1L, length(lines))]
    } else {
      character(0)
    }
  )

  if (isTRUE(write)) {
    writeLines(out, path)
    return(invisible(out))
  }
  out
}
