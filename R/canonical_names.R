# canonical_names -- per-backend parameter-name resolver for flexyBayes.
#
# The `canonical_names()` S3 generic returns a list `(map, transform,
# source)` keyed by backend-native parameter names; values give the
# canonical (brms-style) parameter names; the optional transform
# slot carries value-transform functions (used by INLA's
# precision-to-SD `sqrt(1/prec)` mapping).
#
# Per-backend mapper functions register at package load via the
# internal list `.canonical_mappers`. Each mapper accepts (fit,
# fb_terms) and returns `list(map, transform)`. Users can register
# additional mappers via `register_canonical_mapper()` (internal in
# v0.2; exported in v0.3+ when third-party backend authors might use
# it).
#
# The INLA-side mapping uses summary.fixed rownames for fixed effects
# (already brms-canonical) and translates the "Precision for ..."
# hyperparameter names plus the precision-to-SD value transform. The
# brms-side mapping delegates to Stan's own canonical naming.

# ---------------------------------------------------------------- #
# Internal registry                                                  #
# ---------------------------------------------------------------- #
# Populated at the bottom of this file with the active-engine mappers
# (inla + brms). Additional mappers (e.g., for a future NIMBLE backend)
# register via `register_canonical_mapper(backend, mapper)`.

.canonical_mappers <- new.env(parent = emptyenv())


# ---------------------------------------------------------------- #
# Internal helper -- coalescing operator                            #
# ---------------------------------------------------------------- #

`%||%` <- function(x, y) if (is.null(x)) y else x


# ---------------------------------------------------------------- #
# Generic + methods                                                  #
# ---------------------------------------------------------------- #

#' Canonical parameter-name view for a flexyBayes fit
#'
#' Returns the backend-native -> canonical parameter-name map for a
#' fit, plus per-parameter value transforms where applicable (e.g.,
#' the INLA precision-to-SD `sqrt(1/prec)` transform applied to
#' hyperparameter draws before they are renamed). The canonical convention
#' follows brms (`(Intercept)`, `<term>`, `sd_<group>`, `sigma`,
#' `r_<group>[<level>]`).
#'
#' On `flexybayes_inla` and `flexybayes_brms` fits, the per-backend
#' mapper registered at package load (`inla` or `brms`) drives the
#' resolution; the returned map is cached on
#' `fit$extras$canonical_map` for fast repeated access.
#'
#' @param fit A `flexybayes_inla` or `flexybayes_brms` object.
#' @param drop Logical: if `TRUE` (default `FALSE`), drop
#'   backend-native names that are not in the registered map
#'   (e.g., INLA's `Predictor.<i>` latent-predictor draws). When
#'   `FALSE`, un-mapped names appear in the returned `$unmapped`
#'   element.
#' @param ... Additional arguments (ignored by current methods).
#' @returns A list holding the canonical-name map and its provenance.
#'
#'   \describe{
#'     \item{`map`}{Named character vector keyed by backend-native
#'       parameter name with canonical name as the value.}
#'     \item{`transform`}{Named list of `function(x) -> x'`
#'       transforms keyed by canonical name. Empty list when no
#'       transforms apply.}
#'     \item{`source`}{Character: `"registry"`, `"user"`,
#'       `"registry_fallback_verbatim"`, or `"legacy_inferred"`.}
#'     \item{`unmapped`}{Character vector of backend-native names
#'       not in the map (when `drop = FALSE`).}
#'     \item{`prior_parametrization`}{Character, present only on
#'       aggregated fits, naming which prior the fit ran under.
#'       `"per_row_equivalent"` is the legacy scalar bridge, whose
#'       precision prior makes the aggregated posterior match the
#'       per-row posterior to numerical precision.
#'       `"package_default"` is the automatic half-normal-on-SD
#'       prior, which claims no such equivalence on this route.
#'       `"custom"` is an explicit prior from the caller.}
#'   }
#' @examples
#' # The map is engine-specific, so it needs a fit to read.
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
#'   fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
#'                     verbose = FALSE)
#'   cn <- canonical_names(fit)
#'   head(cn$map)
#' }
#' }
#' @export
canonical_names <- function(fit, drop = FALSE, ...) {
  UseMethod("canonical_names")
}

#' @rdname canonical_names
#' @export
canonical_names.flexybayes_inla <- function(fit, drop = FALSE, ...) {
  if (!is.null(fit$extras$canonical_map)) {
    return(.maybe_attach_unmapped(
      fit$extras$canonical_map,
      fit,
      drop,
      backend = "inla"
    ))
  }
  fb <- fit$extras$fb_terms %||% fit$fb
  mapper <- .canonical_mappers$inla
  if (is.null(mapper)) {
    return(list(map = character(0), transform = list(), source = "no_mapper"))
  }
  res <- mapper(fit, fb)
  res$source <- res$source %||% "registry"
  .maybe_attach_unmapped(res, fit, drop, backend = "inla")
}

#' @rdname canonical_names
#' @export
canonical_names.flexybayes_brms <- function(fit, drop = FALSE, ...) {
  if (!is.null(fit$extras$canonical_map)) {
    return(.maybe_attach_unmapped(
      fit$extras$canonical_map,
      fit,
      drop,
      backend = "brms"
    ))
  }
  fb <- fit$extras$fb_terms %||% fit$extras$parse_info
  mapper <- .canonical_mappers$brms
  if (is.null(mapper)) {
    return(list(map = character(0), transform = list(), source = "no_mapper"))
  }
  res <- mapper(fit, fb)
  res$source <- res$source %||% "registry"
  .maybe_attach_unmapped(res, fit, drop, backend = "brms")
}


# ---------------------------------------------------------------- #
# Registration helper                                               #
# ---------------------------------------------------------------- #

# Internal in v0.2 -- exported in v0.3 once the third-party backend
# author surface is documented. The helper attaches a mapper to the
# .canonical_mappers registry under the given backend key.
register_canonical_mapper <- function(backend, mapper) {
  if (!is.character(backend) || length(backend) != 1L || !nzchar(backend)) {
    stop(
      "`backend` must be a non-empty length-1 character string.",
      call. = FALSE
    )
  }
  if (!is.function(mapper)) {
    stop(
      "`mapper` must be a function (fit, fb_terms) -> ",
      "list(map, transform).",
      call. = FALSE
    )
  }
  assign(backend, mapper, envir = .canonical_mappers)
  invisible(NULL)
}


# ---------------------------------------------------------------- #
# Helpers                                                            #
# ---------------------------------------------------------------- #

# Attach an `$unmapped` element to the canonical-names result when
# requested. The `unmapped` set is the difference between (the
# names present in the fit's posterior draws) and (the keys of the
# canonical map). This surfaces backend-side draws (e.g., INLA's
# Predictor.* linear-predictor samples) that the registry doesn't
# translate.
.maybe_attach_unmapped <- function(res, fit, drop, backend) {
  # Surface the aggregated-fit prior parametrization on the result so a
  # consumer can tell whether the matched-prior (per-row-equivalent)
  # guarantee holds without reaching into the fit internals. Absent on
  # non-aggregated fits (the meta slot is NULL), so the field simply does
  # not appear there. Attached on both the drop = TRUE and drop = FALSE
  # paths since it is independent of the unmapped-draws logic.
  pp <- fit$extras$aggregation_meta$prior_parametrization
  if (!is.null(pp)) {
    res$prior_parametrization <- pp
  }
  if (isTRUE(drop)) {
    return(res)
  }
  draws_names <- tryCatch(
    .backend_draw_names(fit, backend),
    error = function(e) character(0)
  )
  if (length(draws_names) == 0L) {
    res$unmapped <- character(0)
    return(res)
  }
  res$unmapped <- setdiff(draws_names, names(res$map))
  res
}

# Backend-native draw names for the unmapped set. INLA fits expose them
# via rownames(summary.fixed) + rownames(summary.hyperpar); brms fits via
# posterior::as_draws_matrix().
.backend_draw_names <- function(fit, backend) {
  if (identical(backend, "inla")) {
    if (is.null(fit$inla)) {
      return(character(0))
    }
    c(
      rownames(fit$inla$summary.fixed) %||% character(0),
      rownames(fit$inla$summary.hyperpar) %||% character(0)
    )
  } else if (identical(backend, "brms")) {
    if (is.null(fit$brms)) {
      return(character(0))
    }
    if (!requireNamespace("posterior", quietly = TRUE)) {
      return(character(0))
    }
    m <- tryCatch(
      .fb_brms_draws_matrix(fit$brms),
      error = function(e) NULL
    )
    if (is.null(m)) character(0) else colnames(m) %||% character(0)
  } else {
    character(0)
  }
}

# Parse an asreml-emit factor-coefficient name like "tau_env[3,1]"
# into its components: tag = "env", level_index = 3L. Returns NULL
# when the pattern does not match.
.parse_tau_factor <- function(nm) {
  m <- regmatches(nm, regexec("^tau_([^\\[]+)\\[(\\d+),(\\d+)\\]$", nm))[[1]]
  if (length(m) < 4L) {
    return(NULL)
  }
  list(tag = m[2], level = as.integer(m[3]), col = as.integer(m[4]))
}

# Parse an indexed-slope deviation name like
# "slope_dev_f_x_raw[2,1]" against the slope_dev_lookup keyed by
# `<factor>_<continuous>`. Returns NULL when the name does not look
# like an indexed-slope draw, or when no IR term matches the parsed
# tag. The raw-vector index `i` maps to factor level `i + 1` because
# the reference level (level 1) is fixed at zero and absent from the
# raw parameter vector.
# .fb_inla_ranef_registry() --- how ranef() reads every random-effect index
# the INLA emitter generates.
#
# emit_inla() indexes a structured random term by an INTEGER column it builds
# -- `rep_block_id` for `rep:block`, `gen_at_env_E1` for one level of a
# diag(), `geno_id` for a vm() -- and INLA names each summary.random entry
# after the column its f() term was indexed by, with `ID` running 1..n. Two
# things have to be put back before ranef() prints: the GROUP, a name the user
# never wrote and cannot look up, and the LEVELS, which no reader can match to
# a genotype, a site or an animal. They used to live in two maps that grew one
# route at a time, and the level map reached three of the generated indices,
# so a BLUP from diag(), fa(), vm(), ped() or an AR1 field could not be
# compared with ASReml's by name (review finding N5; the interaction and
# correlated-slope cases were finding S01, defect D2).
#
# One entry per random-term type .build_inla_formula() lowers. An entry is
# either a function `(term, ctx)` returning one specification per index the
# emitter writes for that term, or a string recording why the term needs
# none. The label map and the level map are both read off these entries, so
# a new index type is added here, once. test-ranef-index-specs.R reads the
# term types off the formula builder and fails when one has no entry.
#
# @noRd
# @keywords internal
.fb_inla_ranef_registry <- function() {
  own_variable <- paste0(
    "INLA is handed the grouping variable the user wrote, so summary.random ",
    "is named after it and `ID` carries its own values; .fb_random_inla() ",
    "restores a legalised factor level from `level_labels`."
  )
  list(
    simple = own_variable,
    ide = own_variable,
    id = own_variable,
    spline = paste0(
      "INLA is handed the covariate itself as the rw2 index, so `ID` carries ",
      "its distinct values, which are what the smooth is read against."
    ),
    nested = .fb_inla_spec_interaction,
    combo = .fb_inla_spec_interaction,
    # The intercept half of `(x || g)` is f(<g>) on the user's own variable;
    # only the slope half is a generated index.
    simple_slope_uncor = .fb_inla_spec_uncor_slope,
    # Refused on INLA at the gate (correlated_covariance_requires_brms), and
    # no argument of flexybayes() lifts that refusal. emit_inla() still lowers
    # the term and a direct call still fits it, so the index keeps its
    # specification rather than leaking on that route.
    simple_slope_cor = .fb_inla_spec_cor_slope,
    at_simple = .fb_inla_spec_at,
    # Refused at the gate and reachable by a direct emit_inla() call, as the
    # correlated slope is.
    us_gxe = .fb_inla_spec_us,
    fa_gxe = .fb_inla_spec_fa,
    vm = .fb_inla_spec_known_cov,
    ped = .fb_inla_spec_known_cov,
    ar1 = .fb_inla_spec_ar1,
    ar1_spatial = .fb_inla_spec_ar1,
    ar1_spatial_at = .fb_inla_spec_ar1_at
  )
}

# One specification: the group ranef() reports the table under; the levels
# in index order, NULL where there was no data to read them from; the rows of
# INLA's table those levels cover; the length of the field the index spans;
# any index the same f() term reads as `group =` or `replicate =`, which
# INLA gives no table of its own; and the index value INLA writes in `ID` on
# each row of the table, which is how a table is matched to the field by
# identity. `ID` runs 1..n over a field, and restarts at 1 in every block of
# a field read with `group =` or `replicate =` (checked on live fits of
# every type in test-ranef-index-specs.R).
#
# @noRd
# @keywords internal
.fb_inla_index_spec <- function(
  label,
  levels,
  rows = seq_along(levels),
  n_field = length(levels),
  carries = character(0),
  ids = seq_len(n_field)
) {
  list(
    label = label,
    levels = levels,
    rows = rows,
    n_field = n_field,
    carries = carries,
    ids = ids
  )
}

# .fb_inla_index_specs() --- every generated index of a model, specified.
#
# `data` is the fit's own copy, whose factor levels may have been legalised
# for INLA; `level_labels` is the map back to what the user wrote, applied per
# variable before any label is composed, so composition happens once and in
# the user's own vocabulary. `known_matrices` sizes the blocks of a
# blocks-format vm() / ped(), which the IR does not carry. Returns a named
# list keyed by the index name INLA reports.
#
# @noRd
# @keywords internal
.fb_inla_index_specs <- function(
  fb,
  data = NULL,
  level_labels = NULL,
  known_matrices = NULL
) {
  # A factor over `var`, relabelled from the legal levels back to the user's
  # own. Level relabelling is 1:1 and order-preserving, so the level
  # POSITIONS -- which are what an integer index points at -- are the same
  # either way.
  restored_factor <- function(var) {
    if (!is.data.frame(data) || !(var %in% names(data))) {
      return(NULL)
    }
    f <- factor(data[[var]])
    if (!is.null(level_labels)) {
      levels(f) <- .inla_restore_level_labels(level_labels, var, levels(f))
    }
    f
  }
  ctx <- list(
    factor = restored_factor,
    levels = function(var) {
      f <- restored_factor(var)
      if (is.null(f)) NULL else levels(f)
    },
    known_matrices = known_matrices %||% list()
  )

  registry <- .fb_inla_ranef_registry()
  specs <- list()
  for (t in fb$random_terms %||% list()) {
    tt <- t$type %||% ""
    if (tt %in% names(registry) && is.function(registry[[tt]])) {
      specs <- c(specs, registry[[tt]](t, ctx))
    }
  }
  specs
}

# `a:b` and `a:b:c` -- one iid field over the observed cells, first variable
# slowest, which is the interaction emit_inla() indexes by. The separator is
# the one the group label uses, so a level of `rep:block` reads `R1:B1`.
.fb_inla_spec_interaction <- function(t, ctx) {
  vars <- .fb_nested_vars(t)
  cells <- lapply(vars, ctx$factor)
  levs <- NULL
  if (!any(vapply(cells, is.null, logical(1L)))) {
    cells <- as.data.frame(cells)
    names(cells) <- vars
    levs <- levels(
      interaction(cells, drop = TRUE, lex.order = TRUE, sep = ":")
    )
  }
  stats::setNames(
    list(.fb_inla_index_spec(paste(vars, collapse = ":"), levs)),
    .fb_inla_interaction_index(vars)
  )
}

# `(x || g)` -- the slope half is an iid field over the levels of g, weighted
# by x.
.fb_inla_spec_uncor_slope <- function(t, ctx) {
  stats::setNames(
    list(.fb_inla_index_spec(
      paste0(t$slope_var, " | ", t$var),
      ctx$levels(t$var)
    )),
    paste0(t$var, "_for_slope_", t$slope_var)
  )
}

# `(x | g)` -- one iid2d field of length 2 * n_g, intercepts then slopes, and
# INLA returns the whole field under BOTH index names because the slope block
# enters as a copy. Each name reads its own block out of a table carrying
# both, so the two take two labels: mapping both to `<x> | <g>` once returned
# a list with two elements of the same name, and `ranef(fit)[["<x> | <g>"]]`
# reached only the first with nothing to say which block that was.
.fb_inla_spec_cor_slope <- function(t, ctx) {
  levs <- ctx$levels(t$var)
  n_g <- length(levs)
  out <- list(
    .fb_inla_index_spec(
      paste0("(Intercept) | ", t$var),
      levs,
      rows = seq_len(n_g),
      n_field = 2L * n_g
    ),
    .fb_inla_index_spec(
      paste0(t$slope_var, " | ", t$var),
      levs,
      rows = n_g + seq_len(n_g),
      n_field = 2L * n_g
    )
  )
  names(out) <- c(
    .fb_inla_cor_slope_index(t$var, t$slope_var),
    .fb_inla_cor_slope_slope_index(t$var, t$slope_var)
  )
  out
}

# diag(f):g -- one iid field per level of f, each spanning every level of g.
.fb_inla_spec_at <- function(t, ctx) {
  outer_levels <- t$outer_levels %||% character(0)
  inner <- ctx$levels(t$inner)
  out <- lapply(outer_levels, function(lev) {
    .fb_inla_index_spec(paste0(t$inner, " (", t$outer, " = ", lev, ")"), inner)
  })
  names(out) <- vapply(
    outer_levels,
    function(lev) .fb_inla_at_index(t$outer, t$inner, lev),
    character(1L),
    USE.NAMES = FALSE
  )
  out
}

# us(f):g -- one field of length n_f * n_g holding g blocked by f: position
# g + (f - 1) * n_g, so g runs fastest. A level reads `<g>:<f>`, in the order
# the group label names the two.
.fb_inla_spec_us <- function(t, ctx) {
  inner <- ctx$levels(t$inner)
  outer <- ctx$levels(t$outer)
  levs <- NULL
  if (!is.null(inner) && !is.null(outer)) {
    levs <- paste0(
      rep(inner, times = length(outer)),
      ":",
      rep(outer, each = length(inner))
    )
  }
  stats::setNames(
    list(.fb_inla_index_spec(paste0(t$inner, ":", t$outer), levs)),
    .fb_inla_us_index(t$outer, t$inner)
  )
}

# fa(f, k):g -- factor j enters from position j of .fb_fa_levels() onward, as
# its anchor field there and a scaled copy at each later level, and every level of f has
# a specific field. Each table spans the levels of g. A copy's table holds
# that level's loading-scaled contribution rather than the anchor's values,
# so it is labelled as the factor at that level.
.fb_inla_spec_fa <- function(t, ctx) {
  levs <- .fb_fa_levels(t, t$outer_levels %||% character(0))
  inner <- ctx$levels(t$inner)
  out <- list()
  for (j in seq_len(min(t$k %||% 1L, length(levs)))) {
    for (e in seq(j, length(levs))) {
      out[[.fb_inla_fa_index(t$outer, t$inner, j, levs[[e]])]] <-
        .fb_inla_index_spec(
          paste0(t$inner, " factor ", j, " (", t$outer, " = ", levs[[e]], ")"),
          inner
        )
    }
  }
  for (lev in levs) {
    out[[.fb_inla_fa_specific(t$outer, t$inner, lev)]] <- .fb_inla_index_spec(
      paste0(t$inner, " specific (", t$outer, " = ", lev, ")"),
      inner
    )
  }
  out
}

# vm() / ped() -- a generic0 field over the levels of the variable, in the
# order emit_inla() validates the known matrix against. The blocks format
# splits those levels into consecutive blocks, one field each, sized by the
# matrices in `known_matrices`; without them the block count is unknown and
# no block index is specified.
.fb_inla_spec_known_cov <- function(t, ctx) {
  levs <- ctx$levels(t$var)
  cov <- t$cov_representation
  if (is.null(cov) || !identical(cov$format, "blocks")) {
    return(stats::setNames(
      list(.fb_inla_index_spec(t$var, levs)),
      paste0(t$var, "_id")
    ))
  }
  blocks <- ctx$known_matrices[[cov$data]]
  if (!is.list(blocks) || length(blocks) == 0L) {
    return(list())
  }
  sizes <- vapply(blocks, nrow, integer(1L))
  ends <- cumsum(sizes)
  out <- lapply(seq_along(sizes), function(k) {
    block_levels <- if (is.null(levs)) {
      NULL
    } else {
      levs[seq(ends[[k]] - sizes[[k]] + 1L, ends[[k]])]
    }
    .fb_inla_index_spec(paste0(t$var, " (block ", k, ")"), block_levels)
  })
  names(out) <- paste0(t$var, "_id_block_", seq_along(sizes))
  out
}

# ar1(t) -- one AR1 field over the levels of t. ar1(row):ar1(col), ar1(t):g
# and id(trial):ar1(row):ar1(col) -- one f() term indexed by the row
# variable, reading the column variable as `group =` and the trial as
# `replicate =`. INLA lays that table out row fastest, then column, then
# trial: position row + (col - 1) * n_row + (trial - 1) * n_row * n_col. A
# level reads in the order the spelling names the variables, `TA:R1:C1`.
.fb_inla_spec_ar1 <- function(t, ctx) {
  label <- .ar1_term_spelling(t)
  if (identical(t$type, "ar1")) {
    return(stats::setNames(
      list(.fb_inla_index_spec(label, ctx$levels(t$var))),
      paste0(t$var, "_id")
    ))
  }
  per_var <- lapply(c(t$row_var, t$col_var, t$at_var), ctx$levels)
  levs <- NULL
  if (!any(vapply(per_var, is.null, logical(1L)))) {
    cells <- expand.grid(
      per_var,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    if (!is.null(t$at_var)) {
      cells <- cells[c(3L, 1L, 2L)]
    }
    levs <- do.call(paste, c(unname(as.list(cells)), sep = ":"))
  }
  # `ID` is the row index, restarting at 1 in every column and trial block.
  n_row <- length(per_var[[1L]])
  stats::setNames(
    list(.fb_inla_index_spec(
      label,
      levs,
      carries = paste0(c(t$col_var, t$at_var), "_id"),
      ids = rep(seq_len(n_row), times = length(levs) %/% max(n_row, 1L))
    )),
    paste0(t$row_var, "_id")
  )
}

# at(trial):ar1(row):ar1(col) -- one f() term per level of `trial`, indexed
# by that level's row index and reading its column index as `group =`. A
# level reads `T1:R1:C1`, row fastest, as for the shared field.
.fb_inla_spec_ar1_at <- function(t, ctx) {
  out <- lapply(t$by_levels, function(lev) {
    rows <- t$level_grids[[lev]]$rows
    cols <- t$level_grids[[lev]]$cols
    cells <- expand.grid(rows, cols, KEEP.OUT.ATTRS = FALSE,
                         stringsAsFactors = FALSE)
    levs <- paste(lev, cells[[1L]], cells[[2L]], sep = ":")
    idx <- .fb_inla_field_at_index(t, lev)
    .fb_inla_index_spec(
      paste0(.ar1_term_spelling(t), " (", lev, ")"),
      levs,
      carries = idx[["col"]],
      ids = rep(seq_along(rows), times = length(cols))
    )
  })
  names(out) <- vapply(
    t$by_levels,
    function(lev) .fb_inla_field_at_index(t, lev)[["row"]],
    character(1L)
  )
  out
}

# .fb_inla_group_labels() --- generated index name -> the group a user wrote.
#
# Read off .fb_inla_index_specs(). The variance COMPONENTS were mapped back to
# user labels when that machinery was built; the random-effect table was not,
# so `ranef()` reported a grouping factor called `rep_block_id`. Returns a
# named character vector; a group with no entry keeps its own name, which is
# right for a plain `~ g` where the index IS the variable. `known_matrices`
# is needed only to name the blocks of a blocks-format vm() / ped().
.fb_inla_group_labels <- function(fb, known_matrices = NULL) {
  specs <- .fb_inla_index_specs(fb, known_matrices = known_matrices)
  vapply(specs, function(s) s$label, character(1L))
}

# .fb_inla_group_levels() --- generated index name -> the levels a reader
# would recognise, and which rows of INLA's table they cover.
#
# The companion to .fb_inla_group_labels(), read off the same
# specifications. That map fixes the GROUP a summary.random entry is reported
# under; this one fixes the LEVELS inside it. The character-restoration path
# in .fb_random_inla() was written for the plain `~ g` case, where INLA is
# handed the factor itself, and never fires on an integer index: `rep:block`
# on john.alpha reported levels "1" to "18" against ASReml's
# `rep_R1:block_B1 ...`, and `(Days | Subject)` reported "1" to "36" with the
# subject identifiers sitting unused on the fit.
#
# Each entry is `list(levels, rows, n_field, ids)`: the labels in index
# order, the rows of the INLA table they cover, the length of the table the
# index spans, and the `ID` INLA writes on each row of that table. `rows` is
# the whole field for a term with one index, and one block of a longer field
# where two index names share it -- which is how the correlated slope's
# iid2d is laid out.
#
# @noRd
# @keywords internal
.fb_inla_group_levels <- function(
  fb,
  data,
  level_labels = NULL,
  known_matrices = NULL
) {
  if (!is.data.frame(data)) {
    return(list())
  }
  specs <- .fb_inla_index_specs(fb, data, level_labels, known_matrices)
  specs <- Filter(function(s) !is.null(s$levels), specs)
  lapply(specs, function(s) s[c("levels", "rows", "n_field", "ids")])
}

.parse_slope_dev_raw <- function(nm, lookup) {
  m <- regmatches(
    nm,
    regexec("^slope_dev_(.+)_raw\\[(\\d+),(\\d+)\\]$", nm)
  )[[1]]
  if (length(m) < 4L) {
    return(NULL)
  }
  tag <- m[2]
  raw_idx <- as.integer(m[3])
  if (is.null(lookup[[tag]])) {
    return(NULL)
  }
  lvls <- lookup[[tag]]
  # raw index i -> factor level i + 1 (treatment-coded; reference
  # level pinned at zero and absent from the raw vector).
  level_idx <- raw_idx + 1L
  level_label <- if (level_idx >= 1L && level_idx <= length(lvls)) {
    lvls[level_idx]
  } else {
    as.character(level_idx)
  }
  list(
    canonical = paste0("slope_", tag, "[", level_label, "]"),
    tag = tag,
    raw_idx = raw_idx,
    level_label = level_label
  )
}


# ---------------------------------------------------------------- #
# v0.2 mapper: INLA backend                                         #
# ---------------------------------------------------------------- #
# INLA-side names exposed via fit$inla$summary.{fixed,hyperpar}:
#   summary.fixed rownames    -- already brms-canonical:
#       "(Intercept)", "<term>" (treatment-coded factor levels)
#   summary.hyperpar rownames:
#       "Precision for the Gaussian observations" -> "sigma" + sqrt(1/prec)
#       "Precision for <group>"                   -> "sd_<group>" + sqrt(1/prec)

# The canonical names of an AR1 field's correlations, c(row = , col = ) for
# INLA's Rho and GroupRho. On an AR1 x AR1 field each is rho_<factor>, the
# factor it runs along, so ar1(col):ar1(row) and ar1(row):ar1(col) name the
# same correlation alike. A field with one AR1 factor has rho_row only.
.fb_ar1_rho_names <- function(t) {
  if (!isTRUE(t$col_ar1)) {
    return(c(row = "rho_row", col = "rho_col"))
  }
  c(row = paste0("rho_", t$row_var), col = paste0("rho_", t$col_var))
}

.mapper_inla <- function(fit, fb_terms) {
  if (is.null(fit$inla)) {
    return(list(map = character(0), transform = list()))
  }

  map <- character(0)
  transform <- list()

  # Fixed effects -- identity mapping (already brms-canonical). The
  # summary.fixed table uses the bare name ("(Intercept)") while the
  # inla.posterior.sample() latent rows use the ":1"-suffixed form
  # ("(Intercept):1"). Register both so the rename fires regardless
  # of which entry point feeds the consumer.
  fixed_names <- rownames(fit$inla$summary.fixed) %||% character(0)
  for (nm in fixed_names) {
    map[nm] <- nm
    map[paste0(nm, ":1")] <- nm
  }

  # Hyperparameters -- precision-side names; rename + sqrt(1/prec). The
  # transform is keyed by the *native* hyperpar name (not the canonical
  # target), because .fb_canonical_draws() (R/fb_posterior.R) applies
  # transforms to a fit's draws while they still carry their
  # backend-native names, before any renaming. Keying by the canonical
  # name silently no-ops, which would leave INLA precision draws
  # un-converted and hand out precisions under a standard-deviation name.
  # Hyperpar names do not carry the ":1" suffix.
  prec_to_sd <- function(prec) sqrt(1 / prec)
  # WP16: AR1 / separable AR1xAR1 spatial latent-field hyperparameters. INLA
  # names them "Rho for <idx>", "GroupRho for <idx>", "Precision for <idx>"
  # where <idx> is the f() index column (<row_var>_id). Map to canonical
  # correlation names (.fb_ar1_rho_names()) and sd_spatial; the two Rho
  # values are correlations (identity, no transform), the Precision is the
  # field precision -> SD.
  # Matched BEFORE the generic "Precision for <group>" rule so the field
  # precision resolves to sd_spatial, not sd_<idx>.
  # One field per level: its three hyperparameters carry the level in
  # their names, rho_<factor>__<by><level> and sd_spatial__<by><level>.
  field_at <- list()
  for (t in fb_terms$random_terms %||% list()) {
    if (!identical(t$type %||% "", "ar1_spatial_at")) {
      next
    }
    rho <- .fb_ar1_rho_names(t)
    for (lev in t$by_levels) {
      idx <- .fb_inla_field_at_index(t, lev)[["row"]]
      sfx <- paste0("__", t$by_var, lev)
      field_at[[paste0("Rho for ", idx)]] <- paste0(rho[["row"]], sfx)
      field_at[[paste0("GroupRho for ", idx)]] <- paste0(rho[["col"]], sfx)
      field_at[[paste0("Precision for ", idx)]] <- paste0("sd_spatial", sfx)
    }
  }
  spatial_idx <- character(0)
  spatial_rho <- list()
  for (t in c(
    fb_terms$random_terms %||% list(),
    fb_terms$residual_terms %||% list()
  )) {
    tt <- t$type %||% ""
    if (identical(tt, "ar1_spatial")) {
      idx <- paste0(t$row_var, "_id")
    } else if (identical(tt, "ar1")) {
      idx <- paste0(t$var, "_id")
    } else {
      next
    }
    spatial_idx <- c(spatial_idx, idx)
    spatial_rho[[idx]] <- .fb_ar1_rho_names(t)
  }
  # Nested / interaction random effects. The INLA emit fits these as an iid
  # field over a generated interaction index (<a>_<b>_id), so INLA names the
  # hyperparameter after that index. The generated name must not reach the
  # user's summary: map it back to the colon-joined group, which is also what
  # brms calls it, so the two engines line up on one canonical name.
  # The section levels of a dsum() residual, for naming its per-section
  # sigma rows below.
  resid_levels <- NULL
  for (t in fb_terms$residual_terms %||% list()) {
    if (identical(t$type %||% "", "at_units") && is.null(t$level)) {
      resid_levels <- t$var_levels
    }
  }

  nested_idx <- list()
  us_prec <- list()
  us_cor <- list()
  fa_load <- list()
  for (t in fb_terms$random_terms %||% list()) {
    tt <- t$type %||% ""
    # Uncorrelated random slope: the emit writes the slope half as a second
    # iid field over a generated index (<g>_for_slope_<x>). Map it back to
    # the package's own name for that component, <slope>_<group>, which is
    # also the key its prior is registered under.
    # Heterogeneous variance by factor level: one generated index per level
    # of the outer factor. Map each back to brms's name for the same
    # component, sd_<inner>__<outer><level>, so diag() reports identically
    # on both engines.
    # us(f):g -- one iidKd field. INLA names its parts "Precision for <idx>
    # (component j)" and "RhoI:J for <idx>". Map to what brms calls the same
    # quantities for `(0 + f | g)`: sd_<g>__<f><level> and
    # cor_<g>__<f><levelI>__<f><levelJ>.
    # fa(f, k):g -- three kinds of hyperparameter. The anchor field's
    # precision is the common factor's SD; each copy's Beta is a LOADING (a
    # scaling, no transform); each psi field's precision is that level's
    # specific SD, which is the same quantity diag(f):g reports, so it takes
    # the same name.
    if (identical(tt, "fa_gxe")) {
      levs <- .fb_fa_levels(t, t$outer_levels %||% character(0))
      for (j in seq_len(t$k %||% 0L)) {
        anchor <- make.names(paste0(t$inner, "_fa", j, "_", t$outer, "_",
                                    levs[[j]]))
        us_prec[[paste0("Precision for ", anchor)]] <- paste0("sd_fa", j)
        for (e in seq_along(levs)) {
          idx <- make.names(paste0(t$inner, "_fa", j, "_", t$outer, "_",
                                   levs[[e]]))
          fa_load[[paste0("Beta for ", idx)]] <-
            paste0("fa", j, "_", t$outer, levs[[e]])
        }
      }
      for (e in seq_along(levs)) {
        psi <- make.names(paste0(t$inner, "_psi_", t$outer, "_", levs[[e]]))
        us_prec[[paste0("Precision for ", psi)]] <-
          paste0("sd_", t$inner, "__", t$outer, levs[[e]])
      }
      next
    }
    if (identical(tt, "us_gxe")) {
      idx <- make.names(paste0(t$inner, "_by_", t$outer, "_id"))
      levs <- t$outer_levels %||% character(0)
      for (j in seq_along(levs)) {
        us_prec[[paste0("Precision for ", idx, " (component ", j, ")")]] <-
          paste0("sd_", t$inner, "__", t$outer, levs[[j]])
      }
      for (a in seq_along(levs)) {
        for (b in seq_along(levs)) {
          if (b <= a) next
          us_cor[[paste0("Rho", a, ":", b, " for ", idx)]] <-
            paste0("cor_", t$inner, "__", t$outer, levs[[a]],
                   "__", t$outer, levs[[b]])
        }
      }
      next
    }
    if (identical(tt, "at_simple")) {
      for (lev in t$outer_levels %||% character(0)) {
        nested_idx[[make.names(paste0(t$inner, "_at_", t$outer, "_", lev))]] <-
          paste0(t$inner, "__", t$outer, lev)
      }
      next
    }
    if (identical(tt, "simple_slope_uncor")) {
      nested_idx[[paste0(t$var, "_for_slope_", t$slope_var)]] <-
        paste0(t$slope_var, "_", t$var)
      next
    }
    # (x | g) is one iid2d field over a two-block index, so its
    # hyperparameters arrive in the same shape us() uses: two component
    # precisions and one Rho. Component 1 is the intercept block and
    # component 2 the slope block, which fixes the reading of both. The
    # standard deviations take the names the UNCORRELATED sibling already
    # reports, so the two models are read side by side; the correlation
    # takes brms's own cor_<g>__Intercept__<x>.
    if (identical(tt, "simple_slope_cor")) {
      idx <- make.names(paste0(t$var, "_cor_", t$slope_var, "_id"))
      us_prec[[paste0("Precision for ", idx, " (component 1)")]] <-
        paste0("sd_", t$var)
      us_prec[[paste0("Precision for ", idx, " (component 2)")]] <-
        paste0("sd_", t$slope_var, "_", t$var)
      us_cor[[paste0("Rho1:2 for ", idx)]] <-
        paste0("cor_", t$var, "__Intercept__", t$slope_var)
      next
    }
    if (tt %in% c("nested", "combo")) {
      vars <- if (identical(tt, "nested")) {
        c(t$outer, t$inner)
      } else {
        as.character(t$vars)
      }
      nested_idx[[paste0(paste(vars, collapse = "_"), "_id")]] <-
        paste(vars, collapse = ":")
    }
  }

  hyper_names <- rownames(fit$inla$summary.hyperpar) %||% character(0)
  for (nm in hyper_names) {
    sp_hit <- FALSE
    if (!is.null(field_at[[nm]])) {
      map[nm] <- field_at[[nm]]
      if (startsWith(nm, "Precision for ")) {
        transform[[nm]] <- prec_to_sd
      }
      next
    }
    if (!is.null(us_prec[[nm]])) {
      map[nm] <- us_prec[[nm]]
      transform[[nm]] <- prec_to_sd
      next
    }
    if (!is.null(fa_load[[nm]])) {
      # A loading, not a precision: no transform.
      map[nm] <- fa_load[[nm]]
      next
    }
    if (!is.null(us_cor[[nm]])) {
      # A correlation, not a precision: no transform.
      map[nm] <- us_cor[[nm]]
      next
    }
    if (startsWith(nm, "Precision for ")) {
      idx <- sub("^Precision for ", "", nm)
      if (!is.null(nested_idx[[idx]])) {
        map[nm] <- paste0("sd_", nested_idx[[idx]])
        transform[[nm]] <- prec_to_sd
        next
      }
    }
    for (idx in spatial_idx) {
      if (identical(nm, paste0("Rho for ", idx))) {
        map[nm] <- spatial_rho[[idx]][["row"]]
        sp_hit <- TRUE
        break
      }
      if (identical(nm, paste0("GroupRho for ", idx))) {
        map[nm] <- spatial_rho[[idx]][["col"]]
        sp_hit <- TRUE
        break
      }
      if (identical(nm, paste0("Precision for ", idx))) {
        map[nm] <- "sd_spatial"
        transform[[nm]] <- prec_to_sd
        sp_hit <- TRUE
        break
      }
    }
    if (sp_hit) {
      next
    }
    # A sectioned residual gives INLA one Gaussian likelihood per section.
    # It names the first "Precision for the Gaussian observations" and the
    # rest "...[2]", "...[3]". Map them to sigma_<level>, the same
    # vocabulary the brms path already uses for dsum(), so the two engines
    # report the sectioned residual under one set of names.
    if (length(resid_levels) > 1L &&
          startsWith(nm, "Precision for the Gaussian observations")) {
      suffix <- sub("^Precision for the Gaussian observations", "", nm)
      j <- if (nzchar(suffix)) {
        as.integer(gsub("[^0-9]", "", suffix))
      } else {
        1L
      }
      if (!is.na(j) && j >= 1L && j <= length(resid_levels)) {
        map[nm] <- paste0("sigma_", resid_levels[[j]])
        transform[[nm]] <- prec_to_sd
        next
      }
    }
    if (identical(nm, "Precision for the Gaussian observations")) {
      map[nm] <- "sigma"
      transform[[nm]] <- prec_to_sd
      next
    }
    if (startsWith(nm, "Precision for ")) {
      grp <- sub("^Precision for ", "", nm)
      map[nm] <- paste0("sd_", grp)
      transform[[nm]] <- prec_to_sd
      next
    }
    # Other hyperparameters (e.g., rho-for-ar1) left un-mapped for
    # v0.2; future ADR will register them when the structured-cov
    # work lands.
  }

  list(map = map, transform = transform)
}


# ---------------------------------------------------------------- #
# v0.2 mapper: brms-via-Stan backend (Stan passthrough)             #
# ---------------------------------------------------------------- #
# brms parameter names are the canonical names by definition (the
# brms convention IS the flexyBayes canonical convention). The
# mapper is identity over the draw-matrix column names with one
# prefix-strip rule:
#
#   b_<term>           -> <term>             (e.g. b_Days -> Days)
#   b_Intercept        -> "(Intercept)"
#   sd_<group>__Intercept -> sd_<group>      (drop the __Intercept tail)
#   sigma              -> sigma              (identity)
#
# Hyperparameters that brms exposes but flexyBayes does not yet
# canonicalise (e.g. `lp__` log-posterior, `Intercept` -- the
# scaled intercept brms emits alongside b_Intercept) are left
# un-mapped and surface via `$unmapped`.

.mapper_brms_via_stan <- function(fit, fb_terms) {
  if (is.null(fit$brms)) {
    return(list(map = character(0), transform = list()))
  }

  if (!requireNamespace("posterior", quietly = TRUE)) {
    return(list(map = character(0), transform = list()))
  }

  draw_names <- tryCatch(
    colnames(.fb_brms_draws_matrix(fit$brms)),
    error = function(e) character(0)
  )
  if (length(draw_names) == 0L) {
    return(list(map = character(0), transform = list()))
  }

  # Build a {factor, continuous, levels} lookup for every
  # factor_numeric_interaction term in the IR so we can map brms's
  # `b_<fac><level>:<con>` (or `b_<con>:<fac><level>`) parameters
  # to the canonical `slope_<fac>_<con>[<level>]` slot.
  fixed_terms_fb <- if (!is.null(fb_terms)) {
    fb_terms$fixed_terms %||% list()
  } else {
    list()
  }
  fni_terms <- Filter(
    function(t) identical(t$type, "factor_numeric_interaction"),
    fixed_terms_fb
  )

  map <- character(0)
  for (nm in draw_names) {
    if (identical(nm, "b_Intercept")) {
      map[nm] <- "(Intercept)"
      next
    }
    # Factor:continuous indexed-slope translation.
    # brms emits these as `b_<term1>:<term2>` where `<term1>` is
    # `<fac><level>` for the non-reference level. Try both
    # (factor-first, continuous-first) orientations.
    fni_canonical <- .brms_factor_numeric_interaction(nm, fni_terms)
    if (!is.null(fni_canonical)) {
      map[nm] <- fni_canonical
      next
    }
    if (startsWith(nm, "b_")) {
      bare <- sub("^b_", "", nm)
      map[nm] <- bare
      next
    }
    if (startsWith(nm, "sd_")) {
      bare <- sub("^sd_", "", nm)
      # brms emits sd_<g>__<x> for uncorrelated random slopes
      # ((x || g) / (1 + x || g)) and sd_<g>__Intercept for the
      # corresponding intercept-deviation. Canonical name is
      # sd_<x>_<g> for slopes; sd_<g> for the intercept. The
      # __Intercept rule comes first because it is a subset of the
      # __<x> rule.
      if (grepl("__Intercept$", bare)) {
        bare <- sub("__Intercept$", "", bare)
        if (nzchar(bare)) {
          map[nm] <- paste0("sd_", bare)
        }
        next
      }
      m <- regmatches(bare, regexec("^([^_]+(?:_[^_]+)*?)__(.+)$", bare))[[1]]
      if (length(m) == 3L && nzchar(m[2L]) && nzchar(m[3L])) {
        # bare matched <g>__<x>; canonicalise to sd_<x>_<g>.
        g_part <- m[2L]
        x_part <- m[3L]
        map[nm] <- paste0("sd_", x_part, "_", g_part)
        next
      }
      if (nzchar(bare)) {
        map[nm] <- paste0("sd_", bare)
      }
      next
    }
    if (identical(nm, "sigma")) {
      map[nm] <- "sigma"
      next
    }
    # Other names (lp__, Intercept-centered, r_<group>[<lvl>,Intercept],
    # cor_<group>__*, etc.) left un-mapped for v0.2; surfaced via
    # $unmapped when drop = FALSE.
  }

  list(map = map, transform = list())
}

# Helper. Translate a brms parameter name like
# `b_fb:x` or `b_x:fb` (treatment-coded factor:continuous interaction)
# to the canonical `slope_<fac>_<con>[<level>]` slot. Returns NULL
# when the name does not match any factor_numeric_interaction term
# in the IR. Reference-level rows (level == lvls[1]) are absent from
# brms's draws by design (treatment coding); they never match.
.brms_factor_numeric_interaction <- function(nm, fni_terms) {
  if (length(fni_terms) == 0L) {
    return(NULL)
  }
  if (!startsWith(nm, "b_")) {
    return(NULL)
  }
  body <- sub("^b_", "", nm)
  if (!grepl(":", body, fixed = TRUE)) {
    return(NULL)
  }
  parts <- strsplit(body, ":", fixed = TRUE)[[1L]]
  if (length(parts) != 2L) {
    return(NULL)
  }
  for (t in fni_terms) {
    fac <- t$factor
    con <- t$continuous
    lvls <- t$levels %||% NULL
    if (is.null(fac) || is.null(con) || is.null(lvls)) {
      next
    }
    # Try factor-first, then continuous-first.
    for (orient in list(c(1L, 2L), c(2L, 1L))) {
      fac_part <- parts[[orient[[1L]]]]
      con_part <- parts[[orient[[2L]]]]
      if (!identical(con_part, con)) {
        next
      }
      if (!startsWith(fac_part, fac)) {
        next
      }
      lvl <- substr(fac_part, nchar(fac) + 1L, nchar(fac_part))
      if (!nzchar(lvl)) {
        next
      }
      if (lvl %in% lvls && !identical(lvl, lvls[[1L]])) {
        return(paste0("slope_", fac, "_", con, "[", lvl, "]"))
      }
    }
  }
  NULL
}


# ---------------------------------------------------------------- #
# Package-load registration                                         #
# ---------------------------------------------------------------- #
# Run at file-source time, so the registry is populated when the
# namespace is built, before any consumer (.fb_canonical_draws() and the
# draws accessors that read it) runs.

register_canonical_mapper("inla", .mapper_inla)
register_canonical_mapper("brms", .mapper_brms_via_stan)
