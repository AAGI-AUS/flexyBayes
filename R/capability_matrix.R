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
      note = "The certified overlap class, which both engines emit and
        `triangulate()` compares.",
      anchor = "test-backend-conformance.R"
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
      inla = "refuses",
      brms = "fits",
      note = "The INLA mapping named greta as one of its three
        verification arbitrators, so it stays deferred until the criterion
        is rebuilt around the active engines. The deferral is
        host-independent -- no local artefact lifts it. `auto` routes to
        brms.",
      anchor = "test-random-slopes-uncor.R"
    ),
    list(
      model_class = "Factor-by-numeric fixed interaction",
      spelling = "`y ~ f * x` with numeric `x`",
      inla = "refuses",
      brms = "fits",
      note = "The indexed-slope INLA mapping shares the deferred
        three-arbitrator verification with the uncorrelated random slope,
        and refuses on every host. `auto` routes to brms.",
      anchor = "test-inla-verification-artefact-policy.R"
    ),
    list(
      model_class = "Correlated random slope",
      spelling = "`(x | g)`",
      inla = "refuses",
      brms = "refuses",
      note = "Refused at ingest, before any engine is chosen. Fit
        `(x || g)` when the correlation is not of inferential interest.",
      anchor = "test-random-slopes-uncor.R"
    ),
    list(
      model_class = "Nested / interaction random effects, multi-stratum",
      spelling = "`~ gen:env`, `~ env:rep:block`",
      inla = "refuses",
      brms = "fits",
      note = "INLA collapses the finest strata, so it refuses rather than
        reporting a zero. brms emits `(1 | a:b)`.",
      anchor = "test-preflight-interaction-terms.R"
    ),
    list(
      model_class = "Heterogeneous variance by factor level",
      spelling = "`~ diag(f):g`, `~ idh(f):g`, `~ at(f):g`",
      inla = "refuses",
      brms = "fits",
      note = "One variance per level of `f`, no covariance between levels.
        All three spellings emit identical code.",
      anchor = "test-heterogeneous-variance.R"
    ),
    list(
      model_class = "Unstructured genotype-by-environment covariance",
      spelling = "`~ us(f):g`",
      inla = "refuses",
      brms = "fits",
      note = "The correlated sibling of the diagonal structure --
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
      inla = "refuses",
      brms = "fits",
      note = "Lowered to distributional regression on log sigma,
        `sigma ~ 0 + f`. Refused for families with no residual scale.",
      anchor = "test-heterogeneous-variance.R"
    ),
    list(
      model_class = paste(
        "Combined interaction random effects and heterogeneous residual",
        "(full MET)"
      ),
      spelling = "`random = ~ gen + gen:env` with the `dsum` residual",
      inla = "refuses",
      brms = "fits",
      note = "The emit carries both the group-level term and the `sigma`
        predictor, and a live fit samples cleanly on simulated
        multi-environment data. `auto` reaches brms for this class.",
      anchor = "test-met-combined.R"
    ),
    list(
      model_class = "Factor-analytic genotype-by-environment covariance",
      spelling = "`~ fa(env, k):gen`",
      inla = "refuses",
      brms = "refuses",
      note = "Parsed for the formula catalogue and refused at dispatch --
        no active engine emits a factor-analytic covariance.",
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
        block carriers, and brms additionally takes dense and Cholesky.",
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
      model_class = "Univariate P-spline",
      spelling = "`~ spl(x)`",
      inla = "fits",
      brms = "refuses",
      note = "Mapped to INLA's second-order random walk. brms has no
        lowering for the smooth basis.",
      anchor = "test-smooth.R"
    ),
    list(
      model_class = "Observation weights",
      spelling = "`weights = w`",
      inla = "refuses",
      brms = "refuses",
      note = "Parsed and recorded, consumed by no active emitter. A
        non-constant vector refuses rather than returning the
        unweighted posterior.",
      anchor = "test-weights-refusal.R"
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
