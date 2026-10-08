# =============================================================================
# The capability matrix is one R-level table, and every verdict in it is
# re-derived here from the gate and emit code.
#
# The August 2026 external audit found four public surfaces -- the primary
# man page, inst/KNOWN_ISSUES.md, the README table and two tutorials -- each
# stating a different, and each a false, account of what the package fits.
# A hand-maintained Markdown table cannot be tested, so it drifts; the
# emitted-code assertions in the rest of the suite were all green while the
# documents were wrong. Two gates close that:
#
#   (a) STALENESS. The committed block in README.md and inst/KNOWN_ISSUES.md
#       must equal a fresh render of .fb_capability_matrix(). Editing the
#       Markdown by hand fails here; the fix is to edit the R table and
#       re-run tools/generate_capability_matrix.R.
#
#   (b) BEHAVIOUR ANCHORS. Every row is probed against the live code --
#       return_code = TRUE walks the whole gate + emit path without
#       sampling -- and the recorded verdict must match what the engine
#       actually does. A `refuses` row must raise; a `fits` / `emits` row
#       must not. Where the refusal carries a registered reason code or a
#       typed condition class, the anchor asserts on that class rather than
#       on the message text.
#
#   (c) ENGINE CLAIMS. A `refuses` row whose note says that INLA or brms
#       cannot do something names, in its basis field, a block of
#       test-engine-claims.R that calls that engine directly. A refusal that
#       is flexyBayes's own choice says so, and its note makes no claim
#       about the engine. The 0.10.0 review found refusals describing engine
#       limitations nobody had checked, and three of them were false.
#
# The staleness gate reads the source tree, so it skips on an installed
# package (README.md is not installed) in the same way the withdrawn-engine
# source scan does.
# =============================================================================

.cap_pkg_file <- function(relative) {
  candidate <- testthat::test_path("..", "..", relative)
  if (file.exists(candidate)) {
    return(normalizePath(candidate))
  }
  NA_character_
}

# --- (a) staleness ---------------------------------------------------- #

test_that("the capability table renders to a well-formed Markdown block", {
  block <- flexyBayes:::.fb_capability_markdown(markers = TRUE)
  lines <- strsplit(block, "\n", fixed = TRUE)[[1L]]

  expect_identical(lines[[1L]], "<!-- capability-matrix:begin -->")
  expect_identical(lines[[length(lines)]], "<!-- capability-matrix:end -->")
  # Header, separator and one row per model class.
  tab <- flexyBayes:::.fb_capability_matrix()
  expect_true(sum(startsWith(lines, "| ")) == nrow(tab) + 1L)
  # No cell may carry an embedded newline: a wrapped source literal that
  # escaped .fb_squish() would silently break the table.
  expect_false(any(grepl("\n", tab$note, fixed = TRUE)))
})

test_that("every verdict is drawn from the closed vocabulary", {
  tab <- flexyBayes:::.fb_capability_matrix()
  vocab <- flexyBayes:::.FB_CAPABILITY_VERDICTS
  expect_true(all(tab$inla %in% vocab))
  expect_true(all(tab$brms %in% vocab))
  expect_true(all(nzchar(tab$note)))
  expect_true(all(nzchar(tab$anchor)))
})

test_that("every anchor names a test file that exists", {
  # The anchor column is the breadcrumb from a public claim back to the
  # assertion that holds it. A dangling filename would make the claim
  # unfalsifiable, which is the failure mode this whole file exists to
  # close.
  tab <- flexyBayes:::.fb_capability_matrix()
  for (anchor in unique(tab$anchor)) {
    path <- testthat::test_path(anchor)
    expect_true(file.exists(path), label = paste("anchor file", anchor))
  }
})

test_that("the committed capability blocks equal a fresh render", {
  targets <- flexyBayes:::.fb_capability_targets()
  expect_setequal(names(targets), c("README.md", "inst/KNOWN_ISSUES.md"))
  for (relative in names(targets)) {
    path <- .cap_pkg_file(relative)
    skip_if(is.na(path), paste0(relative, " is not in the installed tree"))

    committed <- readLines(path, warn = FALSE)
    fresh <- flexyBayes:::.fb_capability_splice(
      path,
      write = FALSE,
      layout = targets[[relative]]
    )
    expect_identical(
      committed,
      fresh,
      label = paste0(
        relative,
        " is stale -- re-run tools/generate_capability_matrix.R"
      )
    )
  }
})

test_that("the README block lists only what an engine fits, without notes", {
  tab <- flexyBayes:::.fb_capability_matrix()
  block <- flexyBayes:::.fb_capability_markdown(markers = FALSE, layout = "fits")
  lines <- strsplit(block, "\n", fixed = TRUE)[[1L]]
  rows <- lines[startsWith(lines, "| ")][-1L]
  fitted <- tab$inla %in% c("fits", "emits") | tab$brms %in% c("fits", "emits")
  expect_length(rows, sum(fitted))
  expect_identical(lines[[1L]], "| Model class | Spelling | INLA | brms |")
  # A row neither engine fits is absent, and so is every note.
  for (cls in tab$model_class[!fitted]) {
    expect_false(any(startsWith(rows, paste0("| ", cls, " |"))), label = cls)
  }
  expect_false(any(grepl(substr(tab$note[fitted][[1L]], 1L, 40L), rows,
    fixed = TRUE)))
  # The legend names only the verdicts the rows use.
  expect_false(grepl("`emits`", block, fixed = TRUE))
})

test_that("the engine prose in ?flexybayes equals a fresh render", {
  # The backend argument's "What each engine represents" paragraph is
  # generated from the table when the documentation is built. Before it was,
  # the hand-written paragraph said INLA refused diag() and us() and brms
  # refused (x | g), each false for several releases while the table was
  # right.
  path <- .cap_pkg_file("man/flexybayes.Rd")
  skip_if(is.na(path), "man/flexybayes.Rd is not in the installed tree")
  squish <- function(x) gsub("[[:space:]]+", " ", paste(x, collapse = " "))
  expect_true(
    grepl(
      squish(flexyBayes:::.fb_capability_engine_prose()),
      squish(readLines(path, warn = FALSE)),
      fixed = TRUE
    ),
    label = "man/flexybayes.Rd is stale -- re-run devtools::document()"
  )
})

test_that("the translation registry and the capability table agree on refusals", {
  # Two tables describe the same structures from two sides -- whether an
  # engine fits a structure, and what the fit promises -- and a refusal
  # recorded in one and not the other is the drift both exist to prevent.
  # Rows are matched on spelling; every engine a translation row names must
  # carry `refuses` exactly when the row's class is `refused`.
  cap <- flexyBayes:::.fb_capability_matrix()
  reg <- flexyBayes:::.fb_translation_registry()
  norm <- function(s) gsub("[`~[:space:]]", "", s)
  checked <- 0L
  for (i in seq_len(nrow(reg))) {
    hit <- cap[norm(cap$spelling) == norm(reg$spelling[[i]]), , drop = FALSE]
    if (nrow(hit) != 1L) {
      next
    }
    for (engine in strsplit(reg$engine[[i]], ",[[:space:]]*")[[1L]]) {
      expect_identical(
        identical(hit[[engine]], "refuses"),
        identical(reg$class[[i]], "refused"),
        label = paste0(
          reg$structure[[i]], " on ", engine, ": capability verdict ",
          hit[[engine]], ", translation class ", reg$class[[i]]
        )
      )
      checked <- checked + 1L
    }
  }
  # The four rows the 2026-09-13 containment changed, plus the corh() and
  # uncorrelated-slope rows, match on spelling today.
  expect_gte(checked, 8L)
})

# --- (c) engine claims carry a basis --------------------------------------- #

# The sentences of a note that say `engine` lacks something. "No active
# engine" and "neither engine" speak for both engines. The pattern is broad
# on purpose: a false positive costs a basis entry, a false negative lets an
# unchecked claim onto the README.
.cap_engine_claims <- function(note, engine) {
  sentences <- trimws(strsplit(note, "(?<=[.;])\\s+", perl = TRUE)[[1L]])
  lower <- tolower(sentences)
  negation <- paste0(
    "(cannot|can not|does not|do not|has no|have no|carries no|offers no|",
    "lacks|is unable)"
  )
  name <- if (identical(engine, "inla")) "inla" else "(brms|stan)"
  both <- grepl("no (active )?(engine|backend)|neither (active )?(engine|backend)",
    lower)
  named <- grepl(paste0("\\b", name, "('s)?( [a-z-]+){0,3} ", negation), lower)
  sentences[both | named]
}

# test_that() blocks of test-engine-claims.R, by description, with the
# deparsed body of each.
.cap_engine_test_blocks <- function() {
  exprs <- parse(testthat::test_path("test-engine-claims.R"),
    keep.source = FALSE)
  out <- list()
  for (e in exprs) {
    if (is.call(e) && identical(e[[1L]], as.name("test_that"))) {
      out[[as.character(e[[2L]])]] <- paste(deparse(e[[3L]]), collapse = "\n")
    }
  }
  out
}

.cap_engine_call <- c(inla = "INLA::", brms = "brms::|\\.ec_stancode\\(")

test_that("a refusing row claims an engine limitation only with an engine test behind it", {
  tab <- flexyBayes:::.fb_capability_matrix()
  blocks <- .cap_engine_test_blocks()
  claimed <- 0L
  for (i in seq_len(nrow(tab))) {
    for (engine in c("inla", "brms")) {
      basis <- tab[[paste0(engine, "_basis")]][[i]]
      label <- paste0(tab$model_class[[i]], " on ", engine)
      if (!identical(tab[[engine]][[i]], "refuses")) {
        # A basis on a cell that does not refuse would be left behind when a
        # verdict changes.
        expect_identical(basis, "", label = paste(label, "basis"))
        next
      }
      expect_true(nzchar(basis), label = paste(label, "names a basis"))
      if (!identical(basis, "flexyBayes")) {
        body <- blocks[[basis]]
        expect_false(is.null(body), label = paste0(
          label, ": test-engine-claims.R has a block named '", basis, "'"
        ))
        expect_true(
          !is.null(body) && grepl(.cap_engine_call[[engine]], body),
          label = paste0(label, ": '", basis, "' calls ", engine, " directly")
        )
      }
      for (claim in .cap_engine_claims(tab$note[[i]], engine)) {
        claimed <- claimed + 1L
        expect_false(identical(basis, "flexyBayes"), label = paste0(
          label, ": the note says \"", claim, "\", and the basis names no ",
          "engine test"
        ))
      }
    }
  }
  # The hurdle-gamma row makes one such claim, on INLA, with its test.
  expect_gte(claimed, 1L)
})

test_that("the engine-claim pattern finds the claims it exists to find", {
  # The detector itself, on the sentences that were on the README before
  # 0.14.0 and on the ones that replaced them.
  expect_length(.cap_engine_claims(
    "No active engine represents a trait-by-genotype unstructured covariance.",
    "brms"
  ), 1L)
  expect_length(.cap_engine_claims(
    "brms has no factor-analytic lowering, so no engine fits it.", "brms"
  ), 1L)
  expect_length(.cap_engine_claims(
    "INLA has no hurdle-gamma likelihood, so the INLA route refuses it.",
    "inla"
  ), 1L)
  expect_length(.cap_engine_claims(
    "INLA's likelihood roster carries no counterpart.", "inla"
  ), 1L)
  expect_length(.cap_engine_claims(
    "flexyBayes has no brms lowering for the field.", "brms"
  ), 0L)
  expect_length(.cap_engine_claims(
    "INLA has an exchangeable group model, which flexyBayes does not use.",
    "inla"
  ), 0L)
})

# --- (b) behaviour anchors -------------------------------------------- #

.cap_data <- function(seed = 11L) {
  set.seed(seed)
  d <- expand.grid(
    rep = factor(seq_len(2L)),
    gen = factor(seq_len(6L)),
    env = factor(seq_len(3L))
  )
  d$block <- factor(rep(seq_len(4L), length.out = nrow(d)))
  d$trait <- factor(rep(seq_len(2L), length.out = nrow(d)))
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rnorm(nrow(d), 10, 2)
  d$bin <- stats::rbinom(nrow(d), 1L, 0.5)
  d
}

.cap_grid <- function(n_row = 6L, n_col = 5L, seed = 12L) {
  set.seed(seed)
  g <- expand.grid(
    row = factor(seq_len(n_row)),
    col = factor(seq_len(n_col))
  )
  g$y <- stats::rnorm(nrow(g))
  g
}

.cap_kinship <- function() {
  k <- diag(6)
  dimnames(k) <- list(as.character(seq_len(6L)), as.character(seq_len(6L)))
  k
}

# Walk the gate + emit path for one engine without sampling. Returns the
# emitted object, or the condition when the request refuses.
.cap_emit <- function(engine, ..., known_matrices = list()) {
  tryCatch(
    suppressMessages(flexybayes(
      ...,
      known_matrices = known_matrices,
      backend = engine,
      return_code = TRUE,
      verbose = FALSE
    )),
    error = function(e) e
  )
}

# Assert one cell of the matrix against the live code. `verdict` is the
# recorded value; `condition_class` names the typed class the refusal is
# required to carry. Every refusal reachable from the formula surface now
# carries one, so NA is reserved for rows whose refusal is raised by a
# collaborating package rather than by flexyBayes.
.cap_expect <- function(result, verdict, label, condition_class = NA) {
  if (identical(verdict, "n/a")) {
    return(invisible(NULL))
  }
  if (identical(verdict, "refuses")) {
    expect_true(inherits(result, "error"), label = paste(label, "refuses"))
    if (!is.na(condition_class)) {
      expect_true(
        inherits(result, condition_class),
        label = paste(label, "carries", condition_class)
      )
    }
    return(invisible(NULL))
  }
  # fits / emits
  expect_false(
    inherits(result, "error"),
    label = paste0(
      label,
      " ",
      verdict,
      " (got: ",
      if (inherits(result, "error")) conditionMessage(result) else "ok",
      ")"
    )
  )
  invisible(NULL)
}

.cap_row <- function(model_class) {
  tab <- flexyBayes:::.fb_capability_matrix()
  hit <- tab[tab$model_class == model_class, , drop = FALSE]
  if (nrow(hit) != 1L) {
    stop(
      "capability anchor names a row that is not in the table (or names ",
      "it more than once): ",
      model_class,
      call. = FALSE
    )
  }
  hit
}

test_that("Gaussian LMM with a simple random intercept fits on both", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Gaussian LMM, simple random intercept")
  d <- .cap_data()
  .cap_expect(
    .cap_emit("inla", y ~ env, random = ~gen, data = d),
    r$inla,
    "INLA"
  )
  .cap_expect(
    .cap_emit("brms", y ~ env, random = ~gen, data = d),
    r$brms,
    "brms"
  )
})

test_that("a binomial GLMM with a simple random effect fits on both", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row(paste(
    "GLMM (binomial, Poisson, negative binomial, gamma, beta),",
    "simple random effect"
  ))
  d <- .cap_data()
  .cap_expect(
    .cap_emit(
      "inla",
      bin ~ env,
      random = ~gen,
      data = d,
      family = "binomial"
    ),
    r$inla,
    "INLA"
  )
  .cap_expect(
    .cap_emit(
      "brms",
      bin ~ env,
      random = ~gen,
      data = d,
      family = "binomial"
    ),
    r$brms,
    "brms"
  )
})

test_that("a hurdle gamma emits on brms and refuses on INLA", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Hurdle gamma (zero mass plus a positive gamma part)")
  d <- .cap_data()
  # Zero mass plus a positive gamma part, which is what the family is.
  set.seed(21L)
  d$hg <- ifelse(
    stats::rbinom(nrow(d), 1L, 0.7) == 1L,
    stats::rgamma(nrow(d), shape = 2, rate = 0.5),
    0
  )
  .cap_expect(
    .cap_emit(
      "inla",
      hg ~ env,
      random = ~gen,
      data = d,
      family = "hurdle_gamma"
    ),
    r$inla,
    "INLA",
    condition_class = "flexybayes_refusal_inla_gate_refused"
  )
  .cap_expect(
    .cap_emit(
      "brms",
      hg ~ env,
      random = ~gen,
      data = d,
      family = "hurdle_gamma"
    ),
    r$brms,
    "brms"
  )
})

test_that("an uncorrelated random slope fits on brms and defers on INLA", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Uncorrelated random slope")
  d <- .cap_data()
  .cap_expect(
    .cap_emit("inla", y ~ x + (x || gen), data = d),
    r$inla,
    "INLA",
    condition_class = "flexybayes_inla_simple_slope_uncor_deferred"
  )
  .cap_expect(
    .cap_emit("brms", y ~ x + (x || gen), data = d),
    r$brms,
    "brms"
  )
})

test_that("a factor-by-numeric fixed interaction fits on brms only", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Factor-by-numeric fixed interaction")
  d <- .cap_data()
  .cap_expect(
    .cap_emit("inla", y ~ env * x, data = d),
    r$inla,
    "INLA"
  )
  .cap_expect(
    .cap_emit("brms", y ~ env * x, data = d),
    r$brms,
    "brms"
  )
})

test_that("a correlated random slope refuses on INLA and emits on brms", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Correlated random slope")
  d <- .cap_data()
  .cap_expect(
    .cap_emit("inla", y ~ x + (x | gen), data = d),
    r$inla,
    "INLA",
    condition_class = "flexybayes_refusal_correlated_covariance_requires_brms"
  )
  .cap_expect(
    .cap_emit("brms", y ~ x + (x | gen), data = d),
    r$brms,
    "brms"
  )
})

test_that("interaction random effects fit on brms and refuse on INLA", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Nested / interaction random effects, multi-stratum")
  d <- .cap_data()
  .cap_expect(
    .cap_emit("inla", y ~ env, random = ~ gen + gen:env, data = d),
    r$inla,
    "INLA"
  )
  .cap_expect(
    .cap_emit("brms", y ~ env, random = ~ gen + gen:env, data = d),
    r$brms,
    "brms"
  )
})

test_that("diag / idh / at heterogeneous variances fit on brms only", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Heterogeneous variance by factor level")
  d <- .cap_data()
  .cap_expect(
    .cap_emit("inla", y ~ env, random = ~ diag(env):gen, data = d),
    r$inla,
    "INLA"
  )
  for (spelling in list(
    ~ diag(env):gen,
    ~ idh(env):gen,
    ~ at(env):gen
  )) {
    .cap_expect(
      .cap_emit("brms", y ~ env, random = spelling, data = d),
      r$brms,
      paste("brms", deparse(spelling))
    )
  }
})

test_that("us(f):g fits on brms only", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Unstructured genotype-by-environment covariance")
  d <- .cap_data()
  .cap_expect(
    .cap_emit("inla", y ~ env, random = ~ us(env):gen, data = d),
    r$inla,
    "INLA",
    condition_class = "flexybayes_refusal_correlated_covariance_requires_brms"
  )
  .cap_expect(
    .cap_emit("brms", y ~ env, random = ~ us(env):gen, data = d),
    r$brms,
    "brms"
  )
})

test_that("corh(f):g refuses on both, typed", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Heterogeneous variances with one shared correlation")
  d <- .cap_data()
  .cap_expect(
    .cap_emit("inla", y ~ env, random = ~ corh(env):gen, data = d),
    r$inla,
    "INLA",
    condition_class = "flexybayes_refusal_corh_no_equicorrelation_representation"
  )
  .cap_expect(
    .cap_emit("brms", y ~ env, random = ~ corh(env):gen, data = d),
    r$brms,
    "brms",
    condition_class = "flexybayes_refusal_corh_no_equicorrelation_representation"
  )
})

test_that("the dsum heterogeneous residual fits on brms only", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Heterogeneous residual by factor level")
  d <- .cap_data()
  .cap_expect(
    .cap_emit(
      "inla",
      y ~ env,
      random = ~gen,
      residual = ~ dsum(~ units | env),
      data = d
    ),
    r$inla,
    "INLA"
  )
  for (spelling in list(~ dsum(~ units | env), ~ at(env):units)) {
    .cap_expect(
      .cap_emit(
        "brms",
        y ~ env,
        random = ~gen,
        residual = spelling,
        data = d
      ),
      r$brms,
      paste("brms", deparse(spelling))
    )
  }
})

test_that("the combined MET model emits both blocks on brms", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row(paste(
    "Baseline heteroscedastic MET: iid genotype, iid",
    "genotype-by-environment and a sectioned residual"
  ))
  d <- .cap_data()
  .cap_expect(
    .cap_emit(
      "inla",
      y ~ env,
      random = ~ gen + gen:env,
      residual = ~ dsum(~ units | env),
      data = d
    ),
    r$inla,
    "INLA"
  )
  emitted <- .cap_emit(
    "brms",
    y ~ env,
    random = ~ gen + gen:env,
    residual = ~ dsum(~ units | env),
    data = d
  )
  .cap_expect(emitted, r$brms, "brms")

  # The emit half of the row: the reconstructed brms formula carries the
  # interaction group-level term AND the distributional predictor on
  # sigma. The fit half -- that the combination samples with acceptable
  # diagnostics, which is what moved this row from `emits` to `fits` --
  # is asserted on a live fit in test-met-combined.R.
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ gen + gen:env,
    residual = ~ dsum(~ units | env),
    data = d
  )
  bf <- flexyBayes:::.fb_to_brms_formula(fb)
  main <- paste(deparse(bf$formula), collapse = " ")
  sigma <- paste(
    vapply(
      bf$pforms,
      function(p) paste(deparse(p), collapse = " "),
      character(1L)
    ),
    collapse = " "
  )
  expect_match(main, "(1 | gen:env)", fixed = TRUE)
  expect_match(sigma, "sigma ~ 0 + env", fixed = TRUE)
})

test_that("factor-analytic GxE refuses on both engines", {
  # Refused on both from 0.13.0 under one reason code. The refusal is
  # dispatch-wide rather than per-engine because the limit is not INLA's
  # representation -- that works -- but that repeated identical calls did
  # not agree on the covariance, and brms has no factor-analytic lowering
  # to route to. A per-engine refusal would send backend = "auto" to brms
  # to collect a second, less informative message.
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Factor-analytic genotype-by-environment covariance")
  d <- .cap_data()
  .cap_expect(
    .cap_emit("inla", y ~ env, random = ~ fa(env, 1):gen, data = d),
    r$inla,
    "INLA",
    condition_class = "flexybayes_refusal_factor_analytic_not_repeatable"
  )
  .cap_expect(
    .cap_emit("brms", y ~ env, random = ~ fa(env, 1):gen, data = d),
    r$brms,
    "brms",
    condition_class = "flexybayes_refusal_factor_analytic_not_repeatable"
  )
})

test_that("a multi-trait covariance refuses on both", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Multi-trait covariance")
  d <- .cap_data()
  k <- .cap_kinship()
  .cap_expect(
    .cap_emit(
      "inla",
      y ~ env,
      random = ~ us(trait):vm(gen),
      data = d,
      known_matrices = list(K = k)
    ),
    r$inla,
    "INLA",
    condition_class = "flexybayes_refusal_interaction_not_representable"
  )
  .cap_expect(
    .cap_emit(
      "brms",
      y ~ env,
      random = ~ us(trait):vm(gen),
      data = d,
      known_matrices = list(K = k)
    ),
    r$brms,
    "brms",
    condition_class = "flexybayes_refusal_interaction_not_representable"
  )
})

test_that("known-covariance genomic terms fit on both, per carrier", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  skip_if_not_installed("Matrix")
  r <- .cap_row("Known-covariance genomic / pedigree random effect")
  d <- .cap_data()
  k <- .cap_kinship()
  q <- Matrix::Matrix(k, sparse = TRUE)
  dimnames(q) <- dimnames(k)
  # brms takes the dense carrier; INLA takes the sparse precision.
  .cap_expect(
    .cap_emit(
      "brms",
      y ~ env,
      random = ~ vm(gen, K),
      data = d,
      known_matrices = list(K = k)
    ),
    r$brms,
    "brms"
  )
  .cap_expect(
    .cap_emit(
      "inla",
      y ~ env,
      random = ~ vm(gen, cov = fb_cov(Q, type = "precision")),
      data = d,
      known_matrices = list(Q = q)
    ),
    r$inla,
    "INLA"
  )
})

test_that("the separable AR1 field fits on INLA and refuses on brms", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Separable AR1 spatial field")
  g <- .cap_grid()
  .cap_expect(
    .cap_emit("inla", y ~ 1, random = ~ ar1(row):ar1(col), data = g),
    r$inla,
    "INLA"
  )
  .cap_expect(
    .cap_emit("brms", y ~ 1, random = ~ ar1(row):ar1(col), data = g),
    r$brms,
    "brms",
    condition_class = "flexybayes_refusal_stan_cannot_represent_ar1_field"
  )
})

test_that("the per-trial separable AR1 field fits on INLA and refuses on brms", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Per-trial separable AR1 field")
  g <- .cap_grid()
  g <- rbind(
    cbind(g, trial = factor("A")),
    cbind(g, trial = factor("B"))
  )
  .cap_expect(
    .cap_emit(
      "inla", y ~ 1, random = ~ id(trial):ar1(row):ar1(col), data = g
    ),
    r$inla,
    "INLA"
  )
  .cap_expect(
    .cap_emit(
      "brms", y ~ 1, random = ~ id(trial):ar1(row):ar1(col), data = g
    ),
    r$brms,
    "brms",
    condition_class = "flexybayes_refusal_stan_cannot_represent_ar1_field"
  )
})

test_that("the grouped temporal AR1 field fits on INLA and refuses on brms", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Grouped temporal AR1 field")
  g <- .cap_grid(n_row = 8L, n_col = 5L)
  names(g)[names(g) == "row"] <- "t"
  names(g)[names(g) == "col"] <- "g"
  .cap_expect(
    .cap_emit("inla", y ~ 1, random = ~ ar1(t):g, data = g),
    r$inla,
    "INLA"
  )
  .cap_expect(
    .cap_emit("brms", y ~ 1, random = ~ ar1(t):g, data = g),
    r$brms,
    "brms",
    condition_class = "flexybayes_refusal_stan_cannot_represent_ar1_field"
  )
  # The row's spelling column lists the commuted order as the same field.
  .cap_expect(
    .cap_emit("inla", y ~ 1, random = ~ g:ar1(t), data = g),
    r$inla,
    "INLA, commuted spelling"
  )
})

test_that("a univariate P-spline fits on INLA and refuses on brms", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Univariate P-spline")
  d <- .cap_data()
  .cap_expect(
    .cap_emit("inla", y ~ 1, random = ~ spl(x), data = d),
    r$inla,
    "INLA"
  )
  .cap_expect(
    .cap_emit("brms", y ~ 1, random = ~ spl(x), data = d),
    r$brms,
    "brms",
    condition_class = "flexybayes_refusal_brms_cannot_represent_term"
  )
})

test_that("observation weights fit on both for Gaussian, typed refusal otherwise", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  r <- .cap_row("Observation weights (Gaussian, identity link)")
  d <- .cap_data()
  w <- as.numeric(seq_len(nrow(d)))
  for (engine in c("inla", "brms")) {
    .cap_expect(
      .cap_emit(engine, y ~ env, random = ~gen, data = d, weights = w),
      if (identical(engine, "inla")) r$inla else r$brms,
      engine
    )
  }
  # A family the row does not cover (r$inla / r$brms describe the
  # Gaussian-identity cell only) refuses by name on both engines.
  d_poisson <- d
  d_poisson$y <- stats::rpois(nrow(d), 3)
  for (engine in c("inla", "brms")) {
    .cap_expect(
      .cap_emit(
        engine, y ~ env, random = ~gen, data = d_poisson, weights = w,
        family = "poisson"
      ),
      "refuses",
      paste(engine, "poisson"),
      condition_class = "flexybayes_refusal_weights_requires_gaussian"
    )
  }
})

test_that("exact count aggregation is eligible on the INLA path", {
  r <- .cap_row("Exact sufficient-statistic aggregation (binomial, Poisson)")
  expect_identical(r$inla, "fits")
  expect_identical(r$brms, "n/a")
  set.seed(13)
  d <- expand.grid(
    rep = factor(seq_len(20L)),
    gen = factor(seq_len(6L)),
    env = factor(seq_len(3L))
  )
  d$y <- stats::rpois(nrow(d), 4)
  plan <- suppressMessages(
    flexybayes(
      y ~ env,
      random = ~gen,
      data = d,
      family = "poisson",
      plan = TRUE,
      verbose = FALSE
    )
  )
  expect_true(isTRUE(plan$aggregation$eligible))
})

test_that("exact Gaussian aggregation refuses by name on the INLA path", {
  skip_if_not_installed("INLA")
  r <- .cap_row("Exact sufficient-statistic aggregation (Gaussian)")
  expect_identical(r$inla, "refuses")
  expect_identical(r$brms, "n/a")
  set.seed(13)
  d <- expand.grid(
    rep = factor(seq_len(20L)),
    gen = factor(seq_len(6L)),
    env = factor(seq_len(3L))
  )
  d$y <- stats::rnorm(nrow(d), 10, 2)
  expect_error(
    suppressMessages(flexybayes(
      y ~ env,
      random = ~gen,
      data = d,
      backend = "inla",
      aggregate = TRUE,
      verbose = FALSE
    )),
    class = "flexybayes_refusal_aggregation_residual_prior_not_composed"
  )
})


# --- coverage --------------------------------------------------------- #

test_that("every capability row carries a behaviour anchor here", {
  tab <- flexyBayes:::.fb_capability_matrix()
  # Anchors either name a dedicated test file elsewhere in the suite or
  # this one. Either way the row must be probed above; the per-row tests
  # call .cap_row(), which raises on an unknown model_class, so the only
  # remaining failure mode is a row nobody probes. Guard the count.
  probed <- c(
    "Gaussian LMM, simple random intercept",
    paste(
      "GLMM (binomial, Poisson, negative binomial, gamma, beta),",
      "simple random effect"
    ),
    "Hurdle gamma (zero mass plus a positive gamma part)",
    "Uncorrelated random slope",
    "Factor-by-numeric fixed interaction",
    "Correlated random slope",
    "Nested / interaction random effects, multi-stratum",
    "Heterogeneous variance by factor level",
    "Unstructured genotype-by-environment covariance",
    "Heterogeneous variances with one shared correlation",
    "Heterogeneous residual by factor level",
    paste(
      "Baseline heteroscedastic MET: iid genotype, iid",
      "genotype-by-environment and a sectioned residual"
    ),
    "Factor-analytic genotype-by-environment covariance",
    "Multi-trait covariance",
    "Known-covariance genomic / pedigree random effect",
    "Separable AR1 spatial field",
    "Per-trial separable AR1 field",
    "Grouped temporal AR1 field",
    "Univariate P-spline",
    "Observation weights (Gaussian, identity link)",
    "Exact sufficient-statistic aggregation (binomial, Poisson)",
    "Exact sufficient-statistic aggregation (Gaussian)"
  )
  expect_setequal(tab$model_class, probed)
})
