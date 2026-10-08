# capability_matrix.R -- the single source of truth for the per-engine
# capability table shown on every public surface.
#
# One R-level table, rendered to Markdown by
# tools/generate_capability_matrix.R and spliced between the
# <!-- capability-matrix:begin --> / <!-- capability-matrix:end -->
# markers in README.md and inst/KNOWN_ISSUES.md. README.md carries the rows
# at least one engine fits, without notes; inst/KNOWN_ISSUES.md carries every
# row with its note. A staleness test (tests/testthat/test-capability-matrix.R)
# fails when a committed block drifts from a fresh render, and a behaviour
# anchor per row asserts the verdict against the emit / gate code rather than
# against this table.
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
#
# Every `refuses` cell also carries a basis, in `inla_basis` or `brms_basis`:
#
#   "flexyBayes"   the refusal is the package's own choice or gap, and the
#                  note must not say the engine lacks the structure.
#   <description>  the description of a test_that() block in
#                  tests/testthat/test-engine-claims.R that calls the engine
#                  directly and establishes the limitation the note states.
#
# A note that says INLA or brms cannot do something, on a row that refuses,
# fails test-capability-matrix.R unless that engine's basis names such a
# block. The rule exists because the 0.10.0 review found refusals that
# described an engine limitation nobody had checked, and three were false.

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
#'   `model_class`, `spelling`, `inla`, `brms`, `note`, `anchor`,
#'   `inla_basis` and `brms_basis` (empty unless that engine's verdict is
#'   `refuses`).
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
      note = "The class both engines fit, and the one on which the tests
        check that the two agree.",
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
      note = "On INLA each family is checked against `INLA::inla.models()`
        when INLA is installed.",
      anchor = "test-family-support.R"
    ),
    list(
      model_class = "Hurdle gamma (zero mass plus a positive gamma part)",
      spelling = "`family = \"hurdle_gamma\"`",
      inla = "refuses",
      brms = "fits",
      note = "A brms family (parameters mu, shape and hu); `hu` keeps brms's
        own prior. INLA has no hurdle-gamma likelihood, so the INLA route
        refuses it and `backend = \"auto\"` sends it to brms.",
      anchor = "test-family-support.R",
      inla_basis = "INLA's likelihood list has no hurdle or zero-inflated gamma"
    ),
    list(
      model_class = "Uncorrelated random slope",
      spelling = "`(x || g)`",
      inla = "fits",
      brms = "fits",
      note = "Two independent iid fields on INLA, one for the intercept and
        one for the slope, checked against lme4's REML fit of the same
        model. Reported as `sd_<g>` and `sd_<x>_<g>` on INLA and as `sd_<g>`
        and `sd_<g>__<x>` on brms; `canonical_names()` maps both to the INLA
        names.",
      anchor = "test-inla-slope-and-interaction-verified.R"
    ),
    list(
      model_class = "Factor-by-numeric fixed interaction",
      spelling = "`y ~ f * x` with numeric `x`",
      inla = "fits",
      brms = "fits",
      note = "An ordinary treatment-coded design matrix on both engines.
        Under a weak prior the INLA posterior means sit on the least-squares
        estimates of the same model.",
      anchor = "test-inla-slope-and-interaction-verified.R"
    ),
    list(
      model_class = "Correlated random slope",
      spelling = "`(x | g)`",
      inla = "refuses",
      brms = "fits",
      note = "One intercept SD, one slope SD and the correlation between
        them. brms takes `(x | g)` as written and reports `sd_<g>`,
        `sd_<g>__<x>` and `cor_<g>__Intercept__<x>`. flexyBayes fits
        estimated correlation matrices on brms only, and `backend = \"auto\"`
        sends them there. On INLA the structure was one iid2d field under a
        joint Wishart prior. Paired against ASReml 4.2 and lme4 on
        `lme4::sleepstudy` it gave a slope variance of 46.4 at the posterior
        mode, against 35.08 and 35.07, and rescaling the response by ten
        moved the intercept-slope correlation from 0.02 to 0.969, where
        REML stays at 0.066.",
      anchor = "test-inla-containment.R",
      inla_basis = "flexyBayes"
    ),
    list(
      model_class = "Nested / interaction random effects, multi-stratum",
      spelling = "`~ gen:env`, `~ env:rep:block`",
      inla = "fits",
      brms = "fits",
      note = "An iid field over the interaction index on INLA and `(1 | a:b)`
        on brms. They are the same model, and both agree with lme4's REML
        fit of it on `agridat::john.alpha`. Only observed cells get a level,
        as in lme4 and brms.",
      anchor = "test-nested-random-inla.R"
    ),
    list(
      model_class = "Heterogeneous variance by factor level",
      spelling = "`~ diag(f):g`, `~ idh(f):g`, `~ at(f):g`",
      inla = "fits",
      brms = "fits",
      note = "One variance per level of `f`, no covariance between levels.
        The three spellings give identical code: k independent iid fields on
        INLA, `(0 + f || g)` on brms. Both report `sd_<g>__<f><level>` and
        agree with lme4's REML fit of the same model. `at(f, level):g` keeps
        the named levels of `f` only, each with its own variance, and agrees
        with ASReml on both engines.",
      anchor = "test-heterogeneous-variance-inla.R"
    ),
    list(
      model_class = "Unstructured genotype-by-environment covariance",
      spelling = "`~ us(f):g`",
      inla = "refuses",
      brms = "fits",
      note = "Written `(0 + f | g)` on brms. flexyBayes fits estimated
        correlation matrices on brms only, and `backend = \"auto\"` sends
        them there. On INLA the structure was one iidKd field under a joint
        Wishart prior. On `agridat::besag.met` it gave a genetic correlation
        of 0.142 where ASReml gives 0.930, and eight identical calls on the
        response rescaled by ten did not converge and returned correlations
        from -0.894 to +0.686. The structure has k(k+1)/2 parameters against
        the k of `diag()`. With one observation per cell the residual
        variance is confounded with the diagonal of the covariance, and a
        longer chain does not separate them: replicate within cells, or put
        an informative prior on the residual.",
      anchor = "test-heterogeneous-variance.R",
      inla_basis = "flexyBayes"
    ),
    list(
      model_class = "Heterogeneous variances with one shared correlation",
      spelling = "`~ corh(f):g`",
      inla = "refuses",
      brms = "refuses",
      note = "flexyBayes has no lowering for this structure on either
        engine. With one variance for every level and a correlation of zero
        or more, the same covariance is `random = ~ g + f:g`, which fits on
        both engines. INLA has an exchangeable group model with one variance
        and one shared correlation, which flexyBayes does not use. Use
        `diag(f):g` for separate variances without correlation, or `us(f):g`
        on brms for every correlation.",
      anchor = "test-heterogeneous-variance.R",
      inla_basis = "flexyBayes",
      brms_basis = "flexyBayes"
    ),
    list(
      model_class = "Heterogeneous residual by factor level",
      spelling = "`residual = ~ dsum(~ units | f)` / `~ at(f):units`",
      inla = "fits",
      brms = "fits",
      note = "One residual variance per section. INLA fits one Gaussian
        likelihood per level; brms fits `sigma ~ 0 + f`. Both report
        `sigma_<level>` on the standard-deviation scale. With `diag(f):g`
        beside it past 15 variance parameters, INLA fits by empirical Bayes;
        on 12 simulated trials every per-trial variance ASReml identifies is
        within 0.17 of its ASReml standard error. Refused for a family with
        no residual scale, and for `at(f, level):units`, which conditions on
        selected levels.",
      anchor = "test-sectioned-residual-inla.R"
    ),
    list(
      model_class = paste(
        "Baseline heteroscedastic MET: iid genotype, iid",
        "genotype-by-environment and a sectioned residual"
      ),
      spelling = "`random = ~ gen + gen:env` with the `dsum` residual",
      inla = "fits",
      brms = "fits",
      note = "The baseline multi-environment-trial model. On small simulated
        data both engines recover the genotype, genotype-by-environment and
        per-section residual standard deviations. `us(env):gen` with the
        sectioned residual is written out on brms but has not been checked
        against a reference fit, and `fa()` is refused on both engines.",
      anchor = "test-sectioned-residual-inla.R"
    ),
    list(
      model_class = "Factor-analytic genotype-by-environment covariance",
      spelling = "`~ fa(env, k):gen`",
      inla = "refuses",
      brms = "refuses",
      note = "flexyBayes does not fit this structure on either engine. INLA
        can represent it, one latent field per common factor entered at each
        level as a scaled copy, but repeated identical calls did not return
        the same covariance where the genetic signal is weak: on
        `agridat::besag.met` eight calls spanned a ratio of 1.23 in a
        genotype-by-environment variance and 0.12 in a genetic correlation,
        against tolerances of 1.10 and 0.05 set before the runs. Two
        simulated designs and `agridat::dasilva.maize` repeated to better
        than 1.01. The loadings are also anchored by factor-level order, and
        nothing checked that the anchor level carries a loading. flexyBayes
        has no factor-analytic lowering for brms. `diag(env):gen` fits on
        both engines, and `us(env):gen` on brms.",
      anchor = "test-capability-matrix.R",
      inla_basis = "flexyBayes",
      brms_basis = "flexyBayes"
    ),
    list(
      model_class = "Multi-trait covariance",
      spelling = "`~ us(trait):vm(gen)`",
      inla = "refuses",
      brms = "refuses",
      note = "flexyBayes has no lowering for this structure on either
        engine. brms writes it directly, as
        `(0 + trait | gr(gen, cov = K))` on data with one row per trait
        or with `mvbind()` on data with one column per trait. INLA can build
        it from a known-covariance field for each trait joined by a scaled
        `copy`. flexyBayes emits neither. Fit it with brms directly, or fit
        each trait separately.",
      anchor = "test-capability-matrix.R",
      inla_basis = "flexyBayes",
      brms_basis = "flexyBayes"
    ),
    list(
      model_class = "Known-covariance genomic / pedigree random effect",
      spelling = "`~ vm(g, K)`, `~ ped(a, A)`",
      inla = "fits",
      brms = "fits",
      note = "On INLA flexyBayes takes the sparse-precision and
        pedigree-precision carriers and refuses a dense, Cholesky or blocks
        carrier by name; on brms it takes dense, Cholesky and precision
        carriers and refuses a blocks carrier.",
      anchor = "test-known-covariance-inputs.R"
    ),
    list(
      model_class = "Separable AR1 spatial field",
      spelling = "`random = ~ ar1(row):ar1(col)`, `random = ~ ar1(t)`",
      inla = "fits",
      brms = "refuses",
      note = "A latent AR1 field plus the Gaussian observation nugget: four
        hyperparameters, one observation per grid node. ASReml's
        three-parameter residual without a nugget is written `residual = ~
        ar1(row):ar1(col)` and fits on INLA for a Gaussian response; the 1-D
        residual `~ ar1(t)` refuses. flexyBayes has no brms lowering for
        either.",
      anchor = "test-inla-spatial-ar1.R",
      brms_basis = "flexyBayes"
    ),
    list(
      model_class = "Per-trial separable AR1 field",
      spelling = "`random = ~ id(trial):ar1(row):ar1(col)`",
      inla = "fits",
      brms = "refuses",
      note = "One field realisation per level of `trial`, through INLA's
        `replicate =` mechanism. The row correlation, column correlation and
        field SD are shared by every level, not estimated per trial.
        `at(trial):ar1(row):ar1(col)` fits a field with its own SD and
        correlations on each trial, as ASReml does, and `at(trial,
        level):ar1(row):ar1(col)` on the named trials only; both agree with
        ASReml. flexyBayes has no brms lowering for the field.",
      anchor = "test-nested-ar1-field.R",
      brms_basis = "flexyBayes"
    ),
    list(
      model_class = "Grouped temporal AR1 field",
      spelling = "`random = ~ ar1(t):g`, `~ g:ar1(t)`, `~ ar1(t):id(g)`",
      inla = "fits",
      brms = "refuses",
      note = "AR1 along `t` within each level of `g`, independent across
        levels: one correlation, one field SD and the Gaussian nugget, shared
        by every level. This is the units-by-occasions repeated-measures
        covariance of Verbyla and Cullis (1992), written on the random side;
        an occasion a unit missed stays as a row with the response `NA`
        under `na_action = \"augment\"`. Either factor order is accepted.
        flexyBayes has no brms lowering for it.",
      anchor = "test-grouped-temporal-ar1.R",
      brms_basis = "flexyBayes"
    ),
    list(
      model_class = "Univariate P-spline",
      spelling = "`~ spl(x)`",
      inla = "fits",
      brms = "refuses",
      note = "Fitted on INLA as a second-order random walk. It is an analogue
        of ASReml's spline, not a translation: the basis, scaling and
        prediction rule differ, and no equivalence has been established (see
        TRANSLATION.md). flexyBayes has no brms lowering for it.",
      anchor = "test-smooth.R",
      brms_basis = "flexyBayes"
    ),
    list(
      model_class = "Observation weights (Gaussian, identity link)",
      spelling = "`weights = w`",
      inla = "fits",
      brms = "fits",
      note = "Precision weights, Var(y_i) = sigma^2 / w_i, as in ASReml,
        lme4 and `glm(weights = )`: INLA's `scale = w`, and on brms a known
        offset on sigma rather than brms's own `weights()` term, which
        multiplies each observation's log-likelihood instead. Both engines
        match `lme4::lmer(weights = )` closely on a shared simulated
        fixture. Any other family, or a non-identity link on the Gaussian,
        refuses by name, and so does `aggregate = TRUE` with weights.",
      anchor = "test-weights-lowering.R"
    ),
    list(
      model_class =
        "Exact sufficient-statistic aggregation (binomial, Poisson)",
      spelling = "`aggregate = TRUE`, `flexybayes_stream()`",
      inla = "fits",
      brms = "n/a",
      note = "Exact cell-likelihood aggregation for binomial (logit) and
        Poisson (log) models with few cells; the posterior equals the
        per-row posterior. The marginal likelihood is not reported, because
        the engine's value omits a constant that depends on the data.
        flexyBayes has no aggregated route on brms.",
      anchor = "test-aggregation-equivalence-backend.R"
    ),
    list(
      model_class = "Exact sufficient-statistic aggregation (Gaussian)",
      spelling =
        "`aggregate = TRUE`, `flexybayes_stream(family = \"gaussian\")`",
      inla = "refuses",
      brms = "n/a",
      note = "Refused by name. The aggregated Gaussian route combined the
        within-cell sums of squares with a fixed prior on the residual
        precision, so a stated or default residual prior would not be the
        prior fitted. Gaussian models are fitted per row, which is what
        `aggregate = \"auto\"` does; `flexybayes_stream(fit = FALSE)` still
        returns the sufficient statistics.",
      anchor = "test-aggregation-gaussian-containment.R",
      inla_basis = "flexyBayes"
    )
  )

  field <- function(nm) {
    vapply(rows, function(r) .fb_squish(r[[nm]] %||% ""), character(1L))
  }

  data.frame(
    model_class = field("model_class"),
    spelling = field("spelling"),
    inla = field("inla"),
    brms = field("brms"),
    note = field("note"),
    anchor = field("anchor"),
    inla_basis = field("inla_basis"),
    brms_basis = field("brms_basis"),
    stringsAsFactors = FALSE
  )
}

# --- rendering ----------------------------------------------------- #

# The closed verdict vocabulary. Any value outside it is a typo, and the
# renderer refuses rather than printing it onto a public surface.
.FB_CAPABILITY_VERDICTS <- c("fits", "emits", "refuses", "n/a")

# What each verdict means, in the words the rendered legend uses.
.FB_CAPABILITY_LEGEND <- c(
  fits = "the engine fits it, and a test checks the fit.",
  emits = "the engine writes the model, and no fit has yet been checked.",
  refuses = paste(
    "flexyBayes stops with a message saying why and what to use",
    "instead."
  ),
  `n/a` = "the request does not apply to that engine."
)

# Which surface carries which layout. README.md shows only the rows at least
# one engine fits, without notes; inst/KNOWN_ISSUES.md shows every row with
# its note. The generator and the staleness test both read this.
.fb_capability_targets <- function() {
  c("README.md" = "fits", "inst/KNOWN_ISSUES.md" = "full")
}

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
#' @param layout `"full"` (every row, with its note) or `"fits"` (the rows
#'   at least one engine fits or emits, without notes).
#'
#' @returns A single string holding the Markdown block, newline separated.
#'
#' @noRd
#' @keywords internal
.fb_capability_markdown <- function(markers = TRUE, layout = "full") {
  if (!layout %in% c("full", "fits")) {
    stop(
      ".fb_capability_markdown(): `layout` must be \"full\" or \"fits\".",
      call. = FALSE
    )
  }
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

  with_notes <- identical(layout, "full")
  if (!with_notes) {
    fitted <- tab$inla %in% c("fits", "emits") |
      tab$brms %in% c("fits", "emits")
    tab <- tab[fitted, , drop = FALSE]
  }

  header <- if (with_notes) {
    c(
      "| Model class | Spelling | INLA | brms | Notes |",
      "|---|---|:-:|:-:|---|"
    )
  } else {
    c(
      "| Model class | Spelling | INLA | brms |",
      "|---|---|:-:|:-:|"
    )
  }
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
        if (with_notes) paste0(" | ", .fb_md_escape(tab$note[[i]])) else "",
        " |"
      )
    },
    character(1L)
  )

  # The legend names only the verdicts the rendered rows use.
  used <- intersect(.FB_CAPABILITY_VERDICTS, c(tab$inla, tab$brms))
  legend <- c(
    "",
    paste0("`", used, "`: ", .FB_CAPABILITY_LEGEND[used], collapse = " ")
  )

  block <- c(header, body, legend)
  if (isTRUE(markers)) {
    block <- c(
      .FB_CAPABILITY_MARKER_BEGIN,
      .FB_CAPABILITY_GENERATED_NOTE,
      block,
      .FB_CAPABILITY_MARKER_END
    )
  }
  paste(block, collapse = "\n")
}

.FB_CAPABILITY_MARKER_BEGIN <- "<!-- capability-matrix:begin -->"
.FB_CAPABILITY_MARKER_END <- "<!-- capability-matrix:end -->"
.FB_CAPABILITY_GENERATED_NOTE <- paste(
  "<!-- Generated from .fb_capability_matrix() in R/capability_matrix.R by",
  "tools/generate_capability_matrix.R. Edit the R table and re-run the",
  "generator; tests/testthat/test-capability-matrix.R fails on a hand edit",
  "here. -->"
)

#' Render the capability table as the engine prose in `?flexybayes`
#'
#' The `backend` argument's "What each engine represents" paragraph is
#' generated from `.fb_capability_matrix()` when the documentation is built,
#' so it cannot state a capability the table does not. A test compares the
#' committed `man/flexybayes.Rd` against a fresh render.
#'
#' @returns A single string of plain prose.
#'
#' @noRd
#' @keywords internal
.fb_capability_engine_prose <- function() {
  tab <- .fb_capability_matrix()
  classes <- function(engine, verdicts) {
    paste(tab$model_class[tab[[engine]] %in% verdicts], collapse = "; ")
  }
  paste0(
    "INLA fits: ", classes("inla", c("fits", "emits")), ". ",
    "INLA refuses: ", classes("inla", "refuses"), ". ",
    "brms fits: ", classes("brms", c("fits", "emits")), ". ",
    "brms refuses: ", classes("brms", "refuses"), "."
  )
}

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
#' @param layout Passed to `.fb_capability_markdown()`; see
#'   `.fb_capability_targets()` for the layout each surface uses.
#'
#' @returns The full file contents as a character vector of lines,
#'   invisibly when `write` is `TRUE`.
#'
#' @noRd
#' @keywords internal
.fb_capability_splice <- function(path, write = TRUE, layout = "full") {
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
    .fb_capability_markdown(markers = TRUE, layout = layout),
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
