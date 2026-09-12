# Penalty IR --- the joint-prior intermediate representation
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 4.2, 4.3 and
# 4.4, Unit A.
#
# `fb_prior()` expands a declaration into independent per-coefficient rows,
# which is exactly right for a coefficient that has its own prior and exactly
# wrong for a regularised horseshoe: a shared global scale with per-coefficient
# local scales is one object, and expanding it row-wise destroys the coupling
# that makes it a shrinkage prior at all. This file is the second IR the two
# paths meet in: `fb_penalty()` writes unresolved specifications, this file
# turns them into blocks of design columns, and R/penalty_compile_brms.R
# turns the blocks into one brms prior data frame plus, where a family needs
# one, one `stanvars` object.
#
# Nothing here evaluates a model or touches an engine. A block is a statement
# about columns of the candidate design, and the only external fact it
# carries is the coefficient name brms will give each of those columns ---
# which is measured against the installed brms in the compile tests rather
# than recalled here.


# --- the family table --------------------------------------------------- #

# .fb_pen_family_table() --- the one place a shrinkage family's argument
# names, its required arguments, the domain of each, its hyperparameter mode,
# its compilation route and its provider support are written down.
#
# `route` is the single most consequential column, because it is what brms
# 2.23.0 permits rather than what would be convenient (probes A1-A6):
#
#   coef_rows        one brms prior row per candidate coefficient. Exact
#                    block targeting; protected coefficients keep their own
#                    rows.
#   class_wide       one brms special prior on class "b". brms applies it to
#                    every population-level coefficient except the intercept
#                    and refuses a coefficient-specific row alongside it
#                    (probe A3), so the block must be the whole class.
#   stanvar_density  a `target +=` contribution over an index vector carried
#                    in `stanvars`. Exact targeting without a per-coefficient
#                    row, at the cost of an unnormalised density -- which is
#                    harmless only because these hyperparameters are fixed.
#
# `coordinate_wise` records whether the family can tell one column of a block
# from another. An L2 ball cannot: rotate the block and the density is
# unchanged. Everything with an L1 part or a per-coefficient local scale can,
# and its answer then depends on which contrast matrix happened to be in
# force -- which is the coordinate-semantics restriction of v3 section 5.3,
# re-applied here at `resolution = "coefficient"`.
.fb_pen_family_table <- function() {
  list(
    ridge = list(
      params = "scale",
      required = "scale",
      domain = c(scale = "positive"),
      defaults = list(),
      hyperparameter_mode = "fixed",
      route = "coef_rows",
      coordinate_wise = FALSE,
      provider_capabilities = list(brms = TRUE, inla = FALSE, glmnet = TRUE)
    ),
    laplace = list(
      params = "scale",
      required = "scale",
      domain = c(scale = "positive"),
      defaults = list(),
      hyperparameter_mode = "fixed",
      route = "coef_rows",
      coordinate_wise = TRUE,
      provider_capabilities = list(brms = TRUE, inla = FALSE, glmnet = TRUE)
    ),
    elastic_net = list(
      params = c("lambda", "alpha"),
      required = c("lambda", "alpha"),
      domain = c(lambda = "positive", alpha = "unit_open"),
      defaults = list(),
      hyperparameter_mode = "fixed",
      route = "stanvar_density",
      coordinate_wise = TRUE,
      provider_capabilities = list(brms = TRUE, inla = FALSE, glmnet = TRUE)
    ),
    regularized_horseshoe = list(
      params = c("p0", "slab_scale", "slab_df", "df", "df_global", "scope"),
      required = "p0",
      domain = c(
        p0 = "positive", slab_scale = "positive", slab_df = "positive",
        df = "positive", df_global = "positive", scope = "scope"
      ),
      defaults = list(
        slab_scale = 2, slab_df = 4, df = 1, df_global = 1,
        scope = "candidates"
      ),
      hyperparameter_mode = "learned",
      route = "class_wide",
      coordinate_wise = TRUE,
      provider_capabilities = list(brms = TRUE, inla = FALSE, glmnet = FALSE)
    ),
    r2d2 = list(
      params = c("mean_R2", "prec_R2", "cons_D2", "scope"),
      required = character(0),
      domain = c(
        mean_R2 = "unit_open", prec_R2 = "positive", cons_D2 = "positive",
        scope = "scope"
      ),
      defaults = list(
        mean_R2 = 0.5, prec_R2 = 2, cons_D2 = 0.5, scope = "candidates"
      ),
      hyperparameter_mode = "learned",
      route = "class_wide",
      coordinate_wise = TRUE,
      provider_capabilities = list(brms = TRUE, inla = FALSE, glmnet = FALSE)
    ),
    none = list(
      params = character(0),
      required = character(0),
      domain = character(0),
      defaults = list(),
      hyperparameter_mode = "fixed",
      route = "none",
      coordinate_wise = FALSE,
      provider_capabilities = list(brms = TRUE, inla = TRUE, glmnet = TRUE)
    )
  )
}


# The aliases the grammar accepts, and the canonical name each resolves to.
# The certificate records the canonical name, never the alias, so a reader
# cannot tell from the record which spelling was typed --- which is the
# point: `lasso()` and `laplace()` are one prior.
.fb_pen_family_aliases <- function() {
  c(horseshoe = "regularized_horseshoe", lasso = "laplace")
}


# The v3 grammar's group, exclusive and adaptive families. Their component
# traits stay in `.fb_sel_penalty_component_traits()` because the resolver
# still reasons about them; the public grammar refuses them by name until
# their providers ship.
.fb_pen_withdrawn_families <- function() {
  c(
    "group_lasso", "exclusive_lasso", "sparse_group_lasso",
    "group_exclusive", "adaptive"
  )
}


.fb_pen_supported_families <- function() {
  names(.fb_pen_family_table())
}


# --- brms coefficient names --------------------------------------------- #

# .fb_pen_brms_coef_names() --- what brms will call these design columns.
#
# brms rewrites a model-matrix column name before it becomes a coefficient
# name: it deletes the characters that cannot appear in a Stan identifier and
# transliterates the arithmetic ones. `x1:fb c` becomes `x1:fbc`, `fd-e`
# becomes `fdMe`, `I(x^2)` becomes `IxE2`. Selecting a prior by a name brms
# does not use is a silently unapplied prior, which is the defect
# `priors_to_brms` was hardened against at 0.9.2.
#
# The substitution table below is brms 2.23.0's. It is a claim about another
# package, so it is not trusted here: test-penalty-compile.R puts an awkward
# design through `brms::get_prior()` and requires this function to agree
# column for column. If brms changes its table, that test turns red rather
# than the priors going quietly missing.
.fb_pen_brms_coef_names <- function(x) {
  pattern <- c(
    " ", "(", ")", "[", "]", ",", "\"", "'", "?", "+", "-", "*", "/", "^",
    "=", "$"
  )
  replacement <- c(rep("", 9L), "P", "M", "MU", "D", "E", "EQ", "USD")
  out <- as.character(x)
  for (i in seq_along(pattern)) {
    out <- gsub(pattern[i], replacement[i], out, fixed = TRUE)
  }
  out
}


# --- block constructors -------------------------------------------------- #

.fb_pen_block <- function(id, family, target, selector, hyperparameters,
                          scaling) {
  entry <- .fb_pen_family_table()[[family]]
  structure(
    list(
      id = as.integer(id),
      family = family,
      target = target,
      selector = selector,
      columns = integer(0),
      coef_names = character(0),
      protected = FALSE,
      hyperparameters = hyperparameters,
      hyperparameter_mode = entry$hyperparameter_mode,
      scaling = scaling,
      centre = numeric(0),
      scale = numeric(0),
      route = entry$route,
      provider_capabilities = entry$provider_capabilities
    ),
    class = c("fb_prior_block", "list")
  )
}


.fb_pen_blocks <- function(blocks) {
  structure(blocks, class = c("fb_prior_blocks", "list"))
}


#' @export
print.fb_prior_block <- function(x, ...) {
  cat(sprintf(
    "<fb_prior_block %d> %s ~ %s\n", x$id, x$target, x$family
  ))
  if (length(x$hyperparameters) > 0L) {
    cat("  hyperparameters: ", .fb_pen_format_args(x$hyperparameters), "\n",
      sep = ""
    )
  }
  cat(sprintf(
    "  route: %s (%s)  scaling: %s\n",
    x$route, x$hyperparameter_mode, x$scaling
  ))
  if (length(x$coef_names) > 0L) {
    cat("  coefficients: ", paste(x$coef_names, collapse = ", "), "\n",
      sep = ""
    )
  } else {
    cat("  coefficients: unresolved\n")
  }
  invisible(x)
}


#' @export
print.fb_prior_blocks <- function(x, ...) {
  cat(sprintf("<fb_prior_blocks> %d block(s)\n", length(x)))
  for (b in x) {
    print(b)
  }
  invisible(x)
}


# --- resolution ---------------------------------------------------------- #

# .fb_pen_resolve() --- attach design columns and brms coefficient names to
# every specification, and settle which family owns each coordinate.
#
# `map` is a `.fb_sel_term_map()` result: the typed resolver, driven by the
# `assign` attribute `model.matrix()` writes, never by string matching on
# column names. `protected` names the terms held out of selection, in the
# same vocabulary `.fb_sel_partition_design()` uses. The intercept is
# protected under every partition.
#
# Resolution is v3's rule, unchanged: later specifications override earlier
# ones for the coordinates they name, and two equally broad selectors
# matching one coordinate refuse. Breadth is a property of the selector, not
# of how many columns it happened to hit on this data:
#
#   everything()  broadest    -- says nothing about which terms
#   matches()     narrower    -- says something about the names
#   named terms   narrowest   -- says exactly which
#
# so `term(everything()) ~ ridge(1)` followed by `term(x1) ~ laplace(0.5)` is
# the ordinary idiom and resolves without complaint, while two `matches()`
# selectors that both hit `x1` cannot be ordered and refuse.
.fb_pen_resolve <- function(penalty, map, protected = character(0)) {
  if (!inherits(penalty, "fb_penalty")) {
    stop(
      ".fb_pen_resolve(): `penalty` must be an <fb_penalty> object.",
      call. = FALSE
    )
  }
  if (!is.list(map) || is.null(map$term_labels) || is.null(map$columns)) {
    stop(
      ".fb_pen_resolve(): `map` must be a `.fb_sel_term_map()` result.",
      call. = FALSE
    )
  }

  part <- .fb_sel_partition_design(map, protected = protected)
  candidate_terms <- part$candidate_terms
  candidate_columns <- part$candidate_columns
  if (length(candidate_columns) == 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_penalty_target_not_in_model",
      message = paste0(
        "fb_penalty(): every term of the model is protected, so there is ",
        "no candidate design column for a shrinkage prior to act on."
      )
    ))
  }

  # The candidate design's own column names, and the names brms will use for
  # the SAME columns inside the extended model. The two differ, and the prior
  # rows must carry the second.
  raw_names <- map$columns$column[candidate_columns]
  brms_names <- .fb_pen_brms_coef_names(raw_names)

  resolution <- penalty$resolution
  targets <- .fb_sel_targets(map, resolution = resolution)

  specs <- penalty$specs
  claims <- vector("list", length(specs))
  for (i in seq_along(specs)) {
    claims[[i]] <- .fb_pen_claim(
      spec = specs[[i]],
      map = map,
      candidate_terms = candidate_terms,
      candidate_columns = candidate_columns,
      raw_names = raw_names,
      resolution = resolution
    )
  }

  owner <- .fb_pen_owners(claims, length(candidate_columns))

  # The coordinate-semantics gate, at the resolution the caller asked for.
  # `resolution = "term"` cannot reach inside a factor block, so it cannot
  # trip this; `resolution = "coefficient"` can, and a coordinate-wise family
  # is then selecting a contrast column rather than a term.
  .fb_pen_assert_contrast_semantics(
    specs = specs, owner = owner, targets = targets,
    raw_names = raw_names, resolution = resolution
  )

  blocks <- list()
  for (i in seq_along(specs)) {
    cols <- which(owner == i)
    if (length(cols) == 0L) {
      next
    }
    b <- .fb_pen_block(
      id = length(blocks) + 1L,
      family = specs[[i]]$family,
      target = specs[[i]]$target_label,
      selector = specs[[i]]$selector,
      hyperparameters = specs[[i]]$args,
      scaling = penalty$scaling
    )
    b$columns <- as.integer(cols)
    b$coef_names <- brms_names[cols]
    blocks[[length(blocks) + 1L]] <- b
  }

  blocks <- .fb_pen_blocks(blocks)
  .fb_pen_assert_disjoint(blocks)
  blocks <- .fb_pen_resolve_hyperparameters(
    blocks,
    p = length(candidate_columns)
  )

  # The class-wide targeting rule needs the whole `b` vector, not just the
  # candidates: whether a protected fixed effect exists is precisely what
  # decides it. Compiling would find the same thing, but the refusal belongs
  # here, before a Stan program is ever written.
  fixed_columns <- map$columns$index[map$columns$term_index > 0L]
  .fb_pen_assert_targeting(
    blocks,
    coef_names_all = .fb_pen_brms_coef_names(
      map$columns$column[fixed_columns]
    )
  )
  blocks
}


# Turn one specification into the candidate-design column indices it claims,
# together with the breadth of the selector that claimed them.
.fb_pen_claim <- function(spec, map, candidate_terms, candidate_columns,
                          raw_names, resolution) {
  sel <- spec$selector
  if (identical(sel$kind, "everything")) {
    return(list(
      columns = seq_along(candidate_columns), breadth = 3L
    ))
  }

  if (identical(sel$kind, "matches")) {
    pool <- if (identical(resolution, "coefficient")) {
      raw_names
    } else {
      map$columns$term[candidate_columns]
    }
    hit <- grepl(sel$pattern, pool)
    if (!any(hit)) {
      stop(.fb_refusal_condition(
        reason_code = "selection_penalty_target_not_in_model",
        message = paste0(
          "fb_penalty(): `", spec$target_label, "` matches no candidate ",
          if (identical(resolution, "coefficient")) {
            "coefficient. The candidate coefficients are: "
          } else {
            "term. The candidate terms are: "
          },
          paste(
            if (identical(resolution, "coefficient")) {
              raw_names
            } else {
              candidate_terms
            },
            collapse = ", "
          ),
          "."
        )
      ))
    }
    return(list(columns = which(hit), breadth = 2L))
  }

  # kind == "names"
  pool <- if (identical(resolution, "coefficient")) {
    raw_names
  } else {
    candidate_terms
  }
  unknown <- setdiff(sel$names, pool)
  if (length(unknown) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_penalty_target_not_in_model",
      message = paste0(
        "fb_penalty(): `", spec$target_label, "` names ",
        if (identical(resolution, "coefficient")) {
          "coefficient(s) "
        } else {
          "term(s) "
        },
        paste0("`", unknown, "`", collapse = ", "),
        ", which the candidate set does not contain. The candidates are: ",
        paste(pool, collapse = ", "),
        ". A protected term is not a candidate: drop it from `protected` or ",
        "from the penalty."
      )
    ))
  }
  cols <- if (identical(resolution, "coefficient")) {
    which(raw_names %in% sel$names)
  } else {
    which(map$columns$term[candidate_columns] %in% sel$names)
  }
  list(columns = cols, breadth = 1L)
}


# The later-overrides-earlier resolution, and the equal-breadth refusal.
# Returns one owning specification index per candidate column, 0 where no
# specification claimed it.
.fb_pen_owners <- function(claims, p) {
  owner <- integer(p)
  breadth <- integer(p)
  for (i in seq_along(claims)) {
    cols <- claims[[i]]$columns
    b <- claims[[i]]$breadth
    if (length(cols) == 0L) {
      next
    }
    clash <- cols[owner[cols] > 0L & breadth[cols] == b]
    if (length(clash) > 0L) {
      stop(.fb_refusal_condition(
        reason_code = "selection_penalty_target_ambiguous",
        message = paste0(
          "fb_penalty(): specification ", i, " and specification ",
          owner[clash[1L]], " are equally broad selectors and both name ",
          "candidate column(s) ", paste(clash, collapse = ", "),
          ". Later specifications override earlier ones only when they are ",
          "narrower; two selectors of the same breadth leave no order. ",
          "Write one specification per coordinate, or widen one of them."
        ),
        columns = clash
      ))
    }
    owner[cols] <- i
    breadth[cols] <- b
  }
  owner
}


# .fb_pen_assert_disjoint() --- the post-resolution invariant.
#
# Every resolved candidate design column has exactly one role. The
# later-overrides rule is what makes that true, and this is the assertion
# that says so rather than assuming it: a column reachable from two blocks is
# two densities on one coefficient, and neither the compiler nor brms would
# notice.
.fb_pen_assert_disjoint <- function(blocks) {
  cols <- unlist(lapply(blocks, function(b) b$columns), use.names = FALSE)
  if (length(cols) == 0L) {
    return(invisible(blocks))
  }
  dup <- unique(cols[duplicated(cols)])
  if (length(dup) > 0L) {
    owners <- vapply(
      dup,
      function(k) {
        paste(
          vapply(
            blocks[vapply(blocks, function(b) k %in% b$columns, logical(1L))],
            function(b) paste0(b$family, " (block ", b$id, ")"),
            character(1L)
          ),
          collapse = " and "
        )
      },
      character(1L)
    )
    stop(.fb_refusal_condition(
      reason_code = "penalty_block_overlap",
      message = paste0(
        "fb_penalty(): candidate design column(s) ",
        paste(dup, collapse = ", "),
        " are in more than one prior block -- ",
        paste(owners, collapse = "; "),
        ". One coefficient carries one prior."
      ),
      columns = dup
    ))
  }
  invisible(blocks)
}


# The coordinate-semantics restriction of v3 section 5.3, applied to the
# joint-prior layer. A coordinate-wise family aimed at one contrast column of
# a multi-column term selects a basis, not a term, and this release does not
# promise contrast-invariant within-factor selection.
#
# Which columns are basis dependent is not decided here: `.fb_sel_targets()`
# already types every target, and re-deriving the same judgement from the
# term table would be a second opinion that could drift from the first.
.fb_pen_assert_contrast_semantics <- function(specs, owner, targets,
                                              raw_names, resolution) {
  if (!identical(resolution, "coefficient")) {
    return(invisible(NULL))
  }
  tbl <- .fb_pen_family_table()
  semantics <- targets$coordinate_semantics[match(raw_names, targets$target)]
  offending <- character(0)
  families <- character(0)
  for (i in seq_along(specs)) {
    if (!isTRUE(tbl[[specs[[i]]$family]]$coordinate_wise)) {
      next
    }
    cols <- which(owner == i & semantics == "basis_dependent")
    if (length(cols) > 0L) {
      offending <- c(offending, raw_names[cols])
      families <- c(families, specs[[i]]$family)
    }
  }
  if (length(offending) == 0L) {
    return(invisible(NULL))
  }
  stop(.fb_refusal_condition(
    reason_code = "penalty_factor_contrast_ambiguous",
    message = paste0(
      "fb_penalty(): the coordinate-wise family/families ",
      paste(unique(families), collapse = ", "),
      " reach single contrast columns of multi-column terms at ",
      "resolution = \"coefficient\": ",
      paste(unique(offending), collapse = ", "),
      ". Those columns depend on the factor's reference level and contrast ",
      "matrix, so the support would move under a relevel that changes no ",
      "model. Penalise the term as a whole with resolution = \"term\", or ",
      "choose a family that acts on the block rather than on its ",
      "coordinates."
    ),
    coefficients = unique(offending)
  ))
}


# The hyperparameters that cannot be checked until `p` is known, plus the
# class-wide targeting rule. `p` is the number of candidate DESIGN COLUMNS,
# never the number of formula terms (SPEC v5 section 4.3).
.fb_pen_resolve_hyperparameters <- function(blocks, p) {
  for (i in seq_along(blocks)) {
    b <- blocks[[i]]
    if (identical(b$family, "regularized_horseshoe")) {
      p0 <- b$hyperparameters$p0
      if (p0 < 1 || p0 > p - 1) {
        stop(.fb_refusal_condition(
          reason_code = "penalty_hyperparameter_invalid",
          message = paste0(
            "fb_penalty(): `p0 = ", format(p0), "` in ", b$target,
            " ~ regularized_horseshoe() is outside [1, p - 1] for p = ", p,
            " candidate design columns. `p0` is the prior guess at the ",
            "number of non-zero coefficients, so it has to leave at least ",
            "one of each kind."
          ),
          p0 = p0, p = p
        ))
      }
      b$hyperparameters$par_ratio <- p0 / (p - p0)
      b$hyperparameters$p <- p
      blocks[[i]] <- b
    }
  }
  blocks
}


# .fb_pen_assert_targeting() --- what brms will actually do with these blocks.
#
# A class-wide special prior is not aimable. brms puts it on every
# population-level coefficient except the intercept, and refuses a
# coefficient-specific row alongside it (probe A3). So a horseshoe or R2D2
# block is admissible only when the class IS the block:
#
#   * no protected fixed effect beyond the intercept, and
#   * every candidate coefficient inside the block,
#
# unless the caller opts in with `scope = "all_fixed"`, under which the
# protected coefficients sit under the global-local prior too and the
# certificate records the approximation.
#
# `coef_names_all` is the brms `b` vector in order, intercept excluded: it is
# the only thing that knows whether a protected coefficient exists.
.fb_pen_assert_targeting <- function(blocks, coef_names_all) {
  routes <- vapply(blocks, function(b) b$route, character(1L))
  wide <- which(routes == "class_wide")
  if (length(wide) == 0L) {
    return(invisible(blocks))
  }
  if (length(wide) > 1L) {
    stop(.fb_refusal_condition(
      reason_code = "penalty_block_targeting_unsupported",
      message = paste0(
        "fb_penalty(): ", length(wide), " class-wide families were written (",
        paste(
          vapply(blocks[wide], function(b) b$family, character(1L)),
          collapse = ", "
        ),
        "). brms carries one special prior on class \"b\", and it covers ",
        "every population-level coefficient, so two of them cannot be ",
        "aimed at two blocks. Use ridge(), laplace() or elastic_net(), ",
        "which target blocks exactly."
      )
    ))
  }

  b <- blocks[[wide]]
  scope <- b$hyperparameters$scope %||% "candidates"
  covered <- unlist(
    lapply(blocks, function(x) x$coef_names),
    use.names = FALSE
  )
  outside <- setdiff(coef_names_all, b$coef_names)

  if (identical(scope, "all_fixed")) {
    extra <- setdiff(coef_names_all, covered)
    if (length(setdiff(covered, b$coef_names)) > 0L) {
      stop(.fb_refusal_condition(
        reason_code = "penalty_block_targeting_unsupported",
        message = paste0(
          "fb_penalty(): ", b$family, "(scope = \"all_fixed\") puts the ",
          "global-local prior on every population-level coefficient, so no ",
          "other block can hold one. These coefficients are also claimed ",
          "elsewhere: ",
          paste(setdiff(covered, b$coef_names), collapse = ", "), "."
        )
      ))
    }
    return(invisible(blocks))
  }

  if (length(outside) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "penalty_block_targeting_unsupported",
      message = paste0(
        "fb_penalty(): ", b$family, "() is a class-wide prior -- brms ",
        "applies it to every population-level coefficient except the ",
        "intercept -- and these coefficients are outside the block: ",
        paste(outside, collapse = ", "),
        ". Two ways forward: use ridge(), laplace() or elastic_net(), which ",
        "target blocks exactly; or write ", b$family,
        "(scope = \"all_fixed\"), an explicit opt-in under which those ",
        "coefficients also sit under the global-local prior, the projection ",
        "search forces them into every submodel, and the certificate ",
        "records approximations = \"protected_under_global_shrinkage\"."
      ),
      coefficients = outside
    ))
  }
  invisible(blocks)
}


# --- train-fold scaling --------------------------------------------------- #

# .fb_pen_scaling() --- learn the centre and scale of every block column
# inside the training rows, and store them on the block.
#
# The constants are fold-local by construction. A scale learned on the whole
# data and applied inside a fold is a leak: the assessment rows have entered
# the training design through the standardisation, and the out-of-fold score
# is then optimistic by an amount nobody measures. `train_rows` is therefore
# required, not defaulted to every row.
#
# `X` is the CANDIDATE design --- the same matrix `blocks$columns` indexes.
.fb_pen_scaling <- function(blocks, X, train_rows) {
  if (!is.matrix(X)) {
    stop(".fb_pen_scaling(): `X` must be a matrix.", call. = FALSE)
  }
  train_rows <- as.integer(train_rows)
  if (length(train_rows) < 2L) {
    stop(
      ".fb_pen_scaling(): `train_rows` must name at least two rows.",
      call. = FALSE
    )
  }
  if (any(train_rows < 1L) || any(train_rows > nrow(X))) {
    stop(
      ".fb_pen_scaling(): `train_rows` falls outside the rows of `X`.",
      call. = FALSE
    )
  }

  for (i in seq_along(blocks)) {
    b <- blocks[[i]]
    cols <- b$columns
    # An explicitly unpenalised block is left on its own scale. There is no
    # prior scale to state it against, the flat prior brms puts on the
    # coefficient is unchanged by an affine reparameterisation, and an
    # `fb_prior()` row a caller wrote for that coefficient is on the
    # scientific scale -- standardising underneath it would silently rescale
    # a prior somebody chose.
    if (identical(b$scaling, "none") || identical(b$route, "none")) {
      b$centre <- stats::setNames(rep(0, length(cols)), b$coef_names)
      b$scale <- stats::setNames(rep(1, length(cols)), b$coef_names)
      blocks[[i]] <- b
      next
    }
    sub <- X[train_rows, cols, drop = FALSE]
    centre <- colMeans(sub)
    scale <- apply(sub, 2L, stats::sd)
    bad <- !is.finite(scale) | scale <= 0
    if (any(bad)) {
      stop(.fb_refusal_condition(
        reason_code = "penalty_scaling_column_constant",
        message = paste0(
          "fb_penalty(): candidate coefficient(s) ",
          paste(b$coef_names[bad], collapse = ", "),
          " have no spread across the ", length(train_rows),
          " training rows, so train_fold_standardize has no scale to ",
          "divide by. Drop the column from the candidate set, or fit with ",
          "scaling = \"none\" and set the prior scales yourself."
        ),
        coefficients = b$coef_names[bad]
      ))
    }
    b$centre <- stats::setNames(as.numeric(centre), b$coef_names)
    b$scale <- stats::setNames(as.numeric(scale), b$coef_names)
    blocks[[i]] <- b
  }
  blocks
}


# .fb_pen_apply_scaling() --- put the stored constants on the design.
#
# Columns no block claims are left exactly as they were: an unpenalised
# candidate is not standardised, because the scale a prior is stated on is
# the only reason to move it.
.fb_pen_apply_scaling <- function(X, blocks) {
  if (!is.matrix(X)) {
    stop(".fb_pen_apply_scaling(): `X` must be a matrix.", call. = FALSE)
  }
  out <- X
  for (b in blocks) {
    if (length(b$centre) != length(b$columns)) {
      stop(
        ".fb_pen_apply_scaling(): block ", b$id, " carries no scaling ",
        "constants. Call .fb_pen_scaling() first.",
        call. = FALSE
      )
    }
    for (k in seq_along(b$columns)) {
      j <- b$columns[k]
      out[, j] <- (X[, j] - b$centre[[k]]) / b$scale[[k]]
    }
  }
  out
}


# .fb_pen_back_transform() --- posterior draws on the scientific scale.
#
# A coefficient fitted on the standardised column answers a question about a
# one-standard-deviation move; the scientific coefficient answers one about a
# one-unit move, and the two differ by the stored scale. The intercept moves
# as well, because centring the columns moved it: leaving it alone would
# return a posterior whose fitted values are wrong by the sum of the shifts.
#
# `draws` is a matrix of draws with brms-style column names (`b_x1`) or bare
# coefficient names (`x1`). Columns the blocks do not name are returned
# untouched.
.fb_pen_back_transform <- function(draws, blocks) {
  draws <- as.matrix(draws)
  nms <- colnames(draws)
  if (is.null(nms)) {
    stop(
      ".fb_pen_back_transform(): `draws` must carry column names.",
      call. = FALSE
    )
  }
  bare <- sub("^b_", "", nms)
  out <- draws

  intercept <- which(bare %in% c("Intercept", "(Intercept)"))
  shift <- if (length(intercept) == 1L) numeric(nrow(draws)) else NULL

  for (b in blocks) {
    if (length(b$scale) != length(b$coef_names)) {
      stop(
        ".fb_pen_back_transform(): block ", b$id, " carries no scaling ",
        "constants. Call .fb_pen_scaling() first.",
        call. = FALSE
      )
    }
    for (k in seq_along(b$coef_names)) {
      j <- match(b$coef_names[k], bare)
      if (is.na(j)) {
        next
      }
      out[, j] <- draws[, j] / b$scale[[k]]
      if (!is.null(shift)) {
        shift <- shift + out[, j] * b$centre[[k]]
      }
    }
  }
  if (!is.null(shift)) {
    out[, intercept] <- draws[, intercept] - shift
  }
  out
}


# --- formatting ----------------------------------------------------------- #

.fb_pen_format_args <- function(args) {
  if (length(args) == 0L) {
    return("")
  }
  paste(
    vapply(
      names(args),
      function(nm) {
        v <- args[[nm]]
        paste0(nm, " = ", if (is.character(v)) {
          paste0("\"", v, "\"")
        } else {
          format(v)
        })
      },
      character(1L)
    ),
    collapse = ", "
  )
}


# A number written into a Stan program or a brms prior string. R's default
# seven significant digits would silently round a scale on its way into the
# prior; fifteen is a round trip through the double. Scientific notation is
# left available because Stan reads it and a small scale is unreadable
# without it.
.fb_pen_num <- function(x) {
  format(x, digits = 15L, trim = TRUE)
}
