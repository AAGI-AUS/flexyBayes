# Selection candidate space --- the design columns each fold can estimate,
# fixed from the covariates on its training rows before any fold is fitted.
#
# Why a fold needs its own. A design that is identified on every row can lose
# rank on a fold's training rows. Candidates that are constant within a trial
# identify at most one column per training trial, the intercept included, so
# a fold that holds whole locations out can leave more of them than its
# trials can carry. An engine handed that design does not report it: INLA
# fitted it under its default fixed-effect prior, and the fold's score was
# the prior's prediction for the held-out locations; the generalised-least-
# squares fit set the aliased coefficients to zero; and `glmnet` chose among
# the aliased directions by its coordinate order. Each of those is a fold
# scored on a candidate space nobody declared.
#
# The contract. The declared candidate space is the ordered list of candidate
# design columns. For every fold, before any engine or solver runs:
#
#   rank    `qr()` (LINPACK dqrdc2, limited pivoting) on the training rows'
#           protected design followed by its candidate columns, in the order
#           written, with tolerance `.FB_SEL_RANK_TOL`. Where the protected
#           design has an intercept, every other column is centred on the
#           training rows first, so a column is compared with its own
#           variation there and its location cannot create or hide an alias.
#   refuse  a protected column aliased on the training rows
#           (`selection_fold_protected_unidentified`), and a factor level
#           present only in the assessment rows (`selection_fold_level_absent`).
#   drop    a candidate column aliased on the training rows. The fold fits and
#           scores the columns it keeps; a coefficient it reports is indexed
#           by the declared space, with `NA` for a dropped column; and the
#           drop is recorded with its reason.
#
# The rule reads the formulas, the covariates, the fold allocation and the
# grouping columns, and no response value, so every structural alias is
# declared before anything is fitted.

# The rank tolerance, relative to each centred column's own norm on the
# training rows. It is the `qr()` and `lm()` default, the value
# `.fb_sel_estimable_basis()` uses, and the rule by which `lm()` returns an NA
# coefficient. An exact alias leaves a residual of the order of the row count
# times machine epsilon, and a column an engine could use sits many orders
# above it.
.FB_SEL_RANK_TOL <- 1e-7


# .fb_sel_candidate_space() --- the contract for every fold, or a refusal.
#
# `route` says how the fold builders represent a fold. The benchmark route
# (`.fb_sel_fold_substrate()`) builds a design matrix, learns formula-transform
# parameters on the training rows and drops single columns. The Bayesian
# route writes each fold's model as formula text: its designs are built over
# all rows, as `.fb_sel_bayes_fold_inputs()` builds them, and it drops whole
# terms, refusing a term whose columns are aliased only in part.
#
# `groupings` names the columns a dropped column's reason is read against:
# the random-effect groupings and the resampling unit.
.fb_sel_candidate_space <- function(protected,
                                    candidates,
                                    data,
                                    folds,
                                    groupings = character(0),
                                    route = c("benchmark", "bayes"),
                                    tol = .FB_SEL_RANK_TOL) {
  route <- match.arg(route)
  by_fold <- list()
  for (split in folds) {
    train_rows <- if (identical(route, "benchmark")) split$train_rows else NULL
    protected_design <- .fb_sel_design(protected, data, train_rows)
    candidate_design <- if (is.null(candidates)) {
      NULL
    } else {
      .fb_sel_design(candidates, data, train_rows)
    }
    if (identical(route, "benchmark")) {
      .fb_sel_require_fold_local(
        list(protected_design$terms, candidate_design$terms),
        data,
        split
      )
    }
    by_fold[[split$label]] <- .fb_sel_fold_candidate_space(
      protected_design,
      candidate_design,
      data,
      split,
      groupings = groupings,
      resolution = if (identical(route, "benchmark")) "column" else "term",
      tol = tol
    )
  }

  declared <- by_fold[[1L]]$declared
  dropped <- do.call(rbind, c(
    list(data.frame(
      fold = character(0), column = character(0), term = character(0),
      reason = character(0), stringsAsFactors = FALSE
    )),
    lapply(by_fold, function(f) f$dropped)
  ))
  rownames(dropped) <- NULL
  availability <- data.frame(
    column = declared$column,
    term = declared$term,
    n_folds = length(by_fold),
    n_available = vapply(
      declared$column,
      function(column) {
        sum(vapply(
          by_fold,
          function(f) column %in% f$kept_columns,
          logical(1L)
        ))
      },
      integer(1L)
    ),
    stringsAsFactors = FALSE
  )
  rownames(availability) <- NULL

  structure(
    list(
      route = route,
      resolution = if (identical(route, "benchmark")) "column" else "term",
      tol = tol,
      rule = paste0(
        "qr() on each fold's training rows, protected columns then candidate ",
        "columns in the order written, non-intercept columns centred on the ",
        "training rows; a candidate column aliased there is dropped for the ",
        "fold"
      ),
      dropped_coefficient = "NA",
      declared = declared,
      folds = by_fold,
      dropped = dropped,
      availability = availability
    ),
    class = "fb_sel_candidate_space"
  )
}


# .fb_sel_fold_candidate_space() --- one fold's kept and dropped columns.
.fb_sel_fold_candidate_space <- function(protected_design,
                                         candidate_design,
                                         data,
                                         split,
                                         groupings = character(0),
                                         resolution = c("column", "term"),
                                         tol = .FB_SEL_RANK_TOL) {
  resolution <- match.arg(resolution)
  .fb_sel_check_fold_levels(protected_design, split)
  if (!is.null(candidate_design)) {
    .fb_sel_check_fold_levels(candidate_design, split)
  }
  train_rows <- split$train_rows

  x_protected <- protected_design$X[train_rows, , drop = FALSE]
  if (is.null(candidate_design)) {
    x_candidate <- matrix(numeric(0), nrow = length(train_rows), ncol = 0L)
    candidate_terms <- character(0)
  } else {
    labels <- attr(candidate_design$terms, "term.labels")
    keep <- colnames(candidate_design$X) != "(Intercept)"
    x_candidate <- candidate_design$X[train_rows, keep, drop = FALSE]
    candidate_terms <- labels[candidate_design$assign[keep]]
  }
  declared <- data.frame(
    column = colnames(x_candidate),
    term = candidate_terms,
    stringsAsFactors = FALSE
  )

  m <- cbind(x_protected, x_candidate)
  if ("(Intercept)" %in% colnames(x_protected)) {
    centre <- colnames(m) != "(Intercept)"
    m[, centre] <- sweep(
      m[, centre, drop = FALSE], 2L, colMeans(m[, centre, drop = FALSE])
    )
  }
  n_protected <- ncol(x_protected)
  decomposition <- qr(m, tol = tol)
  rank <- decomposition$rank
  aliased <- if (rank < ncol(m)) {
    sort(decomposition$pivot[-seq_len(rank)])
  } else {
    integer(0)
  }
  kept_index <- setdiff(seq_len(ncol(m)), aliased)

  reasons <- .fb_sel_alias_reasons(
    m, aliased, kept_index, data, train_rows, groupings, tol
  )

  aliased_protected <- aliased[aliased <= n_protected]
  if (length(aliased_protected) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_protected_unidentified",
      message = paste0(
        "fb_select(): on fold ", split$label, " the protected design is not ",
        "identified on the training rows: rank ",
        qr(x_protected, tol = tol)$rank, " over ", n_protected,
        " protected columns. These columns are linear combinations of the ",
        "protected columns before them there: ",
        paste0(
          colnames(m)[aliased_protected], " (",
          reasons[match(aliased_protected, aliased)], ")",
          collapse = "; "
        ),
        ". The protected part is the base model and no fold reduces it. Cut ",
        "the folds on a unit whose training rows identify every protected ",
        "column, or remove the column from the base."
      ),
      fold_label = split$label,
      columns = colnames(m)[aliased_protected]
    ))
  }

  aliased_candidate <- aliased - n_protected
  dropped_columns <- declared$column[aliased_candidate]
  if (identical(resolution, "term") && length(dropped_columns) > 0L) {
    touched <- unique(declared$term[aliased_candidate])
    partial <- touched[vapply(
      touched,
      function(t) !all(declared$column[declared$term == t] %in% dropped_columns),
      logical(1L)
    )]
    if (length(partial) > 0L) {
      stop(.fb_refusal_condition(
        reason_code = "selection_fold_candidate_drop_not_representable",
        message = paste0(
          "fb_select(): on fold ", split$label, " the training rows leave ",
          "some design columns of these candidate terms unidentified and ",
          "identify the others: ",
          paste0(
            "`", partial, "` (unidentified: ",
            vapply(
              partial,
              function(t) {
                paste(
                  intersect(
                    declared$column[declared$term == t], dropped_columns
                  ),
                  collapse = ", "
                )
              },
              character(1L)
            ),
            ")",
            collapse = "; "
          ),
          ". The Bayesian route writes each fold's model as formula text, ",
          "which can leave out a whole term and not one of its columns. ",
          "Score the candidates with workflow = \"unpenalized_benchmark\", ",
          "which drops single columns, or write the term as columns the ",
          "training rows identify."
        ),
        fold_label = split$label,
        terms = partial
      ))
    }
  }

  dropped <- data.frame(
    fold = rep(split$label, length(dropped_columns)),
    column = dropped_columns,
    term = declared$term[aliased_candidate],
    reason = reasons[match(aliased_candidate + n_protected, aliased)],
    stringsAsFactors = FALSE
  )
  kept_columns <- setdiff(declared$column, dropped_columns)

  # The margin the tolerance decided on: the smallest ratio of a kept
  # column's residual to its own norm, against the columns before it, and
  # the largest ratio of a dropped column's residual against every kept one.
  norms <- sqrt(colSums(m^2))
  diag_r <- abs(diag(qr.R(decomposition))[seq_len(rank)])
  kept_ratio <- diag_r / norms[decomposition$pivot[seq_len(rank)]]
  dropped_ratio <- if (length(aliased) == 0L) {
    numeric(0)
  } else {
    residual <- qr.resid(
      qr(m[, kept_index, drop = FALSE], tol = tol),
      m[, aliased, drop = FALSE]
    )
    ifelse(norms[aliased] > 0, sqrt(colSums(residual^2)) / norms[aliased], 0)
  }

  list(
    fold = split$label,
    n_train = length(train_rows),
    rank = rank,
    n_protected = n_protected,
    declared = declared,
    kept_columns = kept_columns,
    kept_terms = unique(declared$term[declared$column %in% kept_columns]),
    dropped = dropped,
    dropped_terms = if (identical(resolution, "term")) {
      unique(declared$term[aliased_candidate])
    } else {
      character(0)
    },
    margin = list(
      kept_min = if (length(kept_ratio) == 0L) NA_real_ else min(kept_ratio),
      dropped_max = if (length(dropped_ratio) == 0L) {
        NA_real_
      } else {
        max(dropped_ratio)
      }
    )
  )
}


# .fb_sel_alias_reasons() --- why each aliased column is aliased, in words.
#
# Three reasons, tried in order: the column is constant on the training rows;
# it is constant within a grouping whose training levels are already spanned
# by as many kept columns constant within it (the trial-constant case, with
# the intercept counted); or it is a linear combination of the columns kept
# before it.
.fb_sel_alias_reasons <- function(m,
                                  aliased,
                                  kept_index,
                                  data,
                                  train_rows,
                                  groupings,
                                  tol) {
  if (length(aliased) == 0L) {
    return(character(0))
  }
  groupings <- intersect(groupings, names(data))
  # For each grouping, which columns of `m` are constant within it on the
  # training rows, and how many levels it has there: one pass per grouping.
  within <- lapply(groupings, function(g) {
    key <- factor(as.character(data[[g]][train_rows]))
    counts <- tabulate(key, nbins = nlevels(key))
    means <- rowsum(m, key, reorder = TRUE) / counts
    spread <- apply(abs(m - means[as.integer(key), , drop = FALSE]), 2L, max)
    list(
      levels = nlevels(key),
      constant = spread <= tol * pmax(1, sqrt(colSums(m^2))) |
        colnames(m) == "(Intercept)"
    )
  })
  names(within) <- groupings
  vapply(
    aliased,
    function(j) {
      v <- m[, j]
      if (all(v == v[[1L]]) || sqrt(sum(v^2)) == 0) {
        return("constant on the training rows")
      }
      for (g in groupings) {
        if (!within[[g]]$constant[[j]]) {
          next
        }
        k <- within[[g]]$levels
        spanned <- sum(within[[g]]$constant[kept_index])
        if (spanned >= k) {
          return(paste0(
            "constant within `", g, "`, whose ", k, " levels on the training ",
            "rows are spanned by ", spanned, " kept columns constant within ",
            "it, the intercept included"
          ))
        }
      }
      "a linear combination of the columns kept before it on the training rows"
    },
    character(1L)
  )
}


# .fb_sel_require_no_candidate_drop() --- the penalised reference's
# out-of-fold score compiles one prior over the declared candidate block, so
# it refuses a fold that would drop a candidate.
.fb_sel_require_no_candidate_drop <- function(space) {
  if (nrow(space$dropped) == 0L) {
    return(invisible(NULL))
  }
  stop(.fb_refusal_condition(
    reason_code = "selection_fold_candidate_drop_not_representable",
    message = paste0(
      "fb_select(): the out-of-fold score of the penalised reference ",
      "compiles its prior over the declared candidate block, and the ",
      "training rows of these folds leave candidate columns unidentified: ",
      paste0(
        "fold ", space$dropped$fold, " `", space$dropped$column, "` (",
        space$dropped$reason, ")",
        collapse = "; "
      ),
      ". A prior on a coefficient the fold does not estimate has no ",
      "parameter to act on. Nothing was fitted. Cut the folds on a unit ",
      "whose training rows identify every candidate, or remove the columns ",
      "from `candidates`."
    ),
    columns = unique(space$dropped$column)
  ))
}


# .fb_sel_candidate_space_record() --- the contract as the certificate
# stores it.
.fb_sel_candidate_space_record <- function(space) {
  if (is.null(space)) {
    return(NULL)
  }
  list(
    resolution = space$resolution,
    tol = space$tol,
    rule = space$rule,
    dropped_coefficient = space$dropped_coefficient,
    declared = space$declared$column,
    n_folds = length(space$folds),
    dropped = space$dropped,
    availability = space$availability
  )
}
