# selection_api.R -- the fold data an fb_select() run keeps, and three
# functions that read it.
#
# `fb_select(keep_fold_data = TRUE)` keeps, for the penalised optimiser and
# the unpenalised benchmark, each fold's training and held-out rows, its
# response and design after the fold's own standardisation, the Cholesky
# factor L of the training covariance V_TT, and the covariance of the
# held-out rows given the training rows (R/selection_resamples.R).
# fb_select_fold_data() returns them with the whitened (L^-1) and
# protected-removed forms; fb_select_path() solves the penalised path again
# on them, at the run's alpha or another; fb_select_refit() refits a named
# column set by generalised least squares and scores its held-out clusters.
# None of them fits a model.


# --- keeping the fold data ------------------------------------------------ #

# .fb_sel_keep_fold_data() --- `out` with the fold data of the benchmark run
# `bench` attached.
.fb_sel_keep_fold_data <- function(out, bench, base, resamples, blocks,
                                   decision) {
  substrates <- bench$cv$substrates
  names(substrates) <- vapply(substrates, function(s) s$label, character(1L))
  out$fold_data <- list(
    substrates = substrates,
    blocks = blocks,
    n_lambda = decision$n_lambda %||% 20L,
    lambda_min_ratio = decision$lambda_min_ratio %||% 1e-3,
    cluster = if (length(resamples$unit %||% character(0)) > 0L) {
      .fb_sel_key(base$data, resamples$unit)
    }
  )
  out
}

# .fb_sel_check_keep_fold_data() --- refuse `keep_fold_data = TRUE` on the
# Bayesian workflow, whose folds are model refits with no whitened data.
.fb_sel_check_keep_fold_data <- function(keep_fold_data, workflow) {
  if (!is.logical(keep_fold_data) || length(keep_fold_data) != 1L ||
        is.na(keep_fold_data)) {
    stop("`keep_fold_data` must be TRUE or FALSE.", call. = FALSE)
  }
  if (keep_fold_data && identical(workflow, "bayes")) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_data_not_kept",
      message = paste0(
        "fb_select(): keep_fold_data = TRUE keeps the whitened fold data of ",
        "workflow = \"penalized_optimizer\" or \"unpenalized_benchmark\". ",
        "The Bayesian workflow refits the model in each fold and has none."
      )
    ))
  }
  invisible(keep_fold_data)
}

# .fb_sel_fold_data_get() --- the kept fold data of `x`, restricted to
# `folds`, or a refusal naming the argument that keeps it.
.fb_sel_fold_data_get <- function(x, folds = NULL, caller) {
  if (!inherits(x, "fb_select") || is.null(x$fold_data)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_data_not_kept",
      message = paste0(
        caller, ": `x` holds no fold data. Run fb_select() with ",
        "keep_fold_data = TRUE and workflow = \"penalized_optimizer\" or ",
        "\"unpenalized_benchmark\"."
      )
    ))
  }
  fold_data <- x$fold_data
  if (!is.null(folds)) {
    folds <- as.character(folds)
    unknown <- setdiff(folds, names(fold_data$substrates))
    if (length(unknown) > 0L) {
      stop(
        "`folds` names fold(s) the run does not have: ",
        paste(unknown, collapse = ", "), ". The folds are ",
        paste(names(fold_data$substrates), collapse = ", "), ".",
        call. = FALSE
      )
    }
    fold_data$substrates <- fold_data$substrates[folds]
  }
  fold_data
}

# .fb_sel_fold_data_penalty() --- the run's penalty blocks resolved over the
# candidate columns, or NULL for a run without a penalty.
.fb_sel_fold_data_penalty <- function(fold_data, blocks = fold_data$blocks) {
  if (is.null(blocks) || length(blocks) == 0L) {
    return(NULL)
  }
  first <- fold_data$substrates[[1L]]
  candidate_names <- colnames(first$x_candidate)
  declared <- first$candidate_space$declared
  .fb_sel_glmnet_blocks(
    blocks,
    candidate_names,
    n = length(first$train_rows),
    candidate_terms = declared$term[match(candidate_names, declared$column)]
  )
}

# .fb_sel_whitener_of() --- L^-1 M over a fold's training rows.
.fb_sel_whitener_of <- function(substrate) {
  factor_tt <- substrate$blocks$factor_tt
  function(m) {
    out <- .fb_sel_whiten_train(factor_tt, as.matrix(m))
    colnames(out) <- colnames(m)
    out
  }
}


# --- fb_select_fold_data() -------------------------------------------------- #

#' Whitened data of each fold of a selection run
#'
#' @description
#' Returns, for each fold of an [fb_select()] run kept with
#' `keep_fold_data = TRUE`, the data its generalised least squares and its
#' penalised path were computed from.
#'
#' @details
#' Write `V_TT = L L'` for the fold's training covariance, built from the
#' variance components in `components`. Then `whitened` holds `L^-1 y`,
#' `L^-1 X0` and `L^-1 X` on the training rows, and ordinary least squares on
#' it is generalised least squares on the original rows. `reduced` holds the
#' whitened response and penalised candidates with the protected block
#' projected out: with `Q0` an orthonormal basis of the column space of
#' `L^-1 X0`, `y - Q0 Q0' y` and `X - Q0 Q0' X`. This is the data the
#' penalised optimiser is given. Candidate columns are centred and scaled on
#' each fold's training rows (`centre`, `scale`), and a Gaussian offset is
#' already subtracted from `response`.
#'
#' @param x An `fb_select` object from `fb_select(keep_fold_data = TRUE)`.
#' @param folds Fold labels to return, or `NULL` (the default) for every
#'   fold.
#'
#' @returns A named list with one element per fold. Each element is a list
#'   with `fold`, `train_rows` and `assess_rows` (row numbers of the base's
#'   data), `response` (every row), `x_protected` and `x_candidate` (every
#'   row), `centre` and `scale`, `whitened` (a list of `y`, `x_protected` and
#'   `x_candidate` over the training rows), `reduced` (a list of `y` and `x`
#'   over the training rows), `v_assess` (the covariance of the held-out rows
#'   given the training rows) and `components` (the residual variance, each
#'   random term's variance, and any correlation parameters).
#'
#' @seealso [fb_select_path()], [fb_select_refit()].
#' @family selection fold data
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(site = factor(rep(1:6, each = 5)), x1 = rnorm(30),
#'                   x2 = rnorm(30))
#'   d$y <- 1 + 0.5 * d$x1 + rnorm(6)[d$site] + rnorm(30)
#'   base <- flexybayes(y ~ 1, random = ~site, data = d, backend = "inla",
#'                      verbose = FALSE)
#'   sel <- fb_select(base, ~ x1 + x2, workflow = "unpenalized_benchmark",
#'                    resamples = fb_resamples(unit = ~site, v = 3L),
#'                    seed = 1L, keep_fold_data = TRUE)
#'   fold <- fb_select_fold_data(sel, folds = "1")[["1"]]
#'   dim(fold$whitened$x_candidate)
#' }
#' }
#' @export
fb_select_fold_data <- function(x, folds = NULL) {
  fold_data <- .fb_sel_fold_data_get(x, folds, "fb_select_fold_data()")
  resolved <- .fb_sel_fold_data_penalty(fold_data)
  lapply(fold_data$substrates, function(s) {
    penalised <- resolved$columns %||% seq_len(ncol(s$x_candidate))
    .fb_sel_fold_data_one(s, penalised)
  })
}

# .fb_sel_fold_data_one() --- one fold's public record.
.fb_sel_fold_data_one <- function(substrate, penalised_columns) {
  train_rows <- substrate$train_rows
  whiten <- .fb_sel_whitener_of(substrate)
  reduced <- .fb_sel_glmnet_reduction(substrate, penalised_columns)
  record <- substrate$covariance_record
  list(
    fold = substrate$label,
    train_rows = train_rows,
    assess_rows = substrate$assess_rows,
    response = substrate$response,
    x_protected = substrate$x_protected,
    x_candidate = substrate$x_candidate,
    centre = substrate$learned$preprocessing_centre,
    scale = substrate$learned$preprocessing_scale,
    whitened = list(
      y = drop(whiten(matrix(substrate$response[train_rows], ncol = 1L))),
      x_protected = whiten(substrate$x_protected[train_rows, , drop = FALSE]),
      x_candidate = whiten(substrate$x_candidate[train_rows, , drop = FALSE])
    ),
    reduced = list(y = reduced$payload$y, x = reduced$payload$x),
    v_assess = .fb_sel_conditional_rows(
      substrate$assessment$v_s_given_t,
      seq_along(substrate$assess_rows)
    ),
    components = list(
      residual = record$residual,
      terms = record$terms,
      rho = record$rho %||% list()
    )
  )
}


# --- fb_select_path() ------------------------------------------------------ #

#' Penalised path of a selection run at another alpha
#'
#' @description
#' Solves the penalised path of an [fb_select()] run again on its kept fold
#' data, at the run's own elastic-net mixing `alpha` or at another, and
#' returns each fold's support at each point of the path. No model is
#' refitted: the folds' covariances and standardisation are the run's.
#'
#' @details
#' The grid is the run's: `n_lambda` points of `s = lambda / lambda_max`
#' from 1 down to `lambda_min_ratio`, each fold's `lambda_max` read from its
#' own training rows. At the run's `alpha`, the supports equal those the run
#' reported. A different `alpha` needs every penalty block to be
#' `elastic_net()`; a grouped penalty is solved by the package's own
#' proximal solver and is not re-solved here.
#'
#' @param x An `fb_select` object from
#'   `fb_select(workflow = "penalized_optimizer", keep_fold_data = TRUE)`.
#' @param alpha The elastic-net mixing parameter, a single number in
#'   `(0, 1]`, or `NULL` (the default) for the run's own.
#' @param n_lambda,lambda_min_ratio The grid, or `NULL` (the default) for the
#'   run's own.
#' @param folds Fold labels to solve, or `NULL` (the default) for every fold.
#'
#' @returns A data frame with one row per fold and path point: `fold`,
#'   `lambda_index`, `s`, `lambda` (the canonical lambda), `n_support` and
#'   `support`, the selected candidate columns sorted and joined by `";"`.
#'   Attribute `alpha` holds the alpha solved at.
#'
#' @seealso [fb_select_fold_data()], [fb_select_refit()].
#' @family selection fold data
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE) &&
#'   requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(site = factor(rep(1:6, each = 5)), x1 = rnorm(30),
#'                   x2 = rnorm(30))
#'   d$y <- 1 + 0.5 * d$x1 + rnorm(6)[d$site] + rnorm(30)
#'   base <- flexybayes(y ~ 1, random = ~site, data = d, backend = "inla",
#'                      verbose = FALSE)
#'   sel <- fb_select(
#'     base, ~ x1 + x2,
#'     penalty = fb_penalty(term(everything()) ~
#'                            elastic_net(lambda = 1, alpha = 1)),
#'     workflow = "penalized_optimizer",
#'     resamples = fb_resamples(unit = ~site, v = 3L),
#'     seed = 1L, keep_fold_data = TRUE
#'   )
#'   head(fb_select_path(sel, alpha = 0.5))
#' }
#' }
#' @export
fb_select_path <- function(x,
                           alpha = NULL,
                           n_lambda = NULL,
                           lambda_min_ratio = NULL,
                           folds = NULL) {
  fold_data <- .fb_sel_fold_data_get(x, folds, "fb_select_path()")
  blocks <- .fb_sel_path_blocks(fold_data$blocks, alpha)
  resolved <- .fb_sel_fold_data_penalty(fold_data, blocks)
  penalty <- resolved$penalty
  if (isTRUE(penalty$grouped)) {
    .fb_sel_path_refuse(paste0(
      "the run's penalty is grouped, and a grouped penalty is solved by the ",
      "package's own proximal solver, which fb_select_path() does not call."
    ))
  }
  .fb_sel_require_glmnet()
  decision <- .fb_sel_glmnet_decision(list(
    benchmark = TRUE,
    n_lambda = n_lambda %||% fold_data$n_lambda,
    lambda_min_ratio = lambda_min_ratio %||% fold_data$lambda_min_ratio
  ))
  semantics <- .fb_sel_support_semantics(resolved$family)
  convention <- .fb_sel_glmnet_convention()
  substrates <- fold_data$substrates
  parts <- lapply(substrates, function(s) {
    .fb_sel_glmnet_reduction(s, resolved$columns)
  })
  grid <- .fb_sel_glmnet_grid(parts, penalty, decision)

  # Each fold, each point of the grid ------------------------------------
  rows <- list()
  for (label in names(substrates)) {
    s <- substrates[[label]]
    p <- parts[[label]]
    convention_fold <- .fb_sel_convention_at(convention, p$payload$y)
    for (index in seq_along(grid$lambda[[label]])) {
      lambda <- grid$lambda[[label]][[index]]
      solution <- .fb_sel_glmnet_fold_solution(
        s, p, lambda, penalty$alpha_composite, convention_fold,
        penalty = penalty,
        groups = NULL,
        selects = semantics$selects,
        at_null_anchor = isTRUE(grid$s[[index]] >= 1)
      )
      beta <- solution$beta_penalised
      support <- sort(names(beta)[abs(beta) > 0])
      rows[[length(rows) + 1L]] <- data.frame(
        fold = label,
        lambda_index = index,
        s = grid$s[[index]],
        lambda = lambda,
        n_support = length(support),
        support = paste(support, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "alpha") <- penalty$alpha_composite
  out
}

# .fb_sel_path_blocks() --- the run's penalty blocks, at `alpha` where one
# is given.
.fb_sel_path_blocks <- function(blocks, alpha) {
  if (is.null(blocks) || length(blocks) == 0L) {
    .fb_sel_path_refuse(paste0(
      "the run has no penalty: it is the unpenalised benchmark, whose ",
      "models are the protected model and the full model."
    ))
  }
  if (is.null(alpha)) {
    return(blocks)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
        alpha <= 0 || alpha > 1) {
    stop("`alpha` must be a single number in (0, 1].", call. = FALSE)
  }
  families <- vapply(blocks, function(b) b$family, character(1L))
  if (!all(families == "elastic_net")) {
    .fb_sel_path_refuse(paste0(
      "`alpha` sets the mixing of an elastic_net() penalty, and the run's ",
      "penalty is ", paste(unique(families), collapse = ", "), "."
    ))
  }
  lapply(blocks, function(b) {
    b$hyperparameters$alpha <- alpha
    b
  })
}

# .fb_sel_path_refuse() --- fb_select_path() asked for a path the run's
# penalty does not give.
.fb_sel_path_refuse <- function(detail) {
  stop(.fb_refusal_condition(
    reason_code = "selection_penalty_spec_malformed",
    message = paste0("fb_select_path(): ", detail)
  ))
}


# --- fb_select_refit() ---------------------------------------------------- #

#' Refit a named column set on the folds of a selection run
#'
#' @description
#' Refits the protected model plus a named set of candidate columns, without
#' a penalty, on each fold of an [fb_select()] run kept with
#' `keep_fold_data = TRUE`, and scores each held-out cluster.
#'
#' @details
#' On each fold the coefficients are the generalised least squares estimate
#' `(X_S' V_TT^-1 X_S)^-1 X_S' V_TT^-1 y_T`, where `X_S` holds the protected
#' columns and the named candidates on the training rows and `V_TT` is the
#' fold's training covariance. A column the training rows leave aliased with
#' the columns before it is reported with an `NA` coefficient and listed in
#' `not_identified`. The held-out prediction is the one the run used
#' (`prediction` in [fb_select()]). Each cluster's score is its mean squared
#' error, root mean squared error and negative log predictive density per
#' row under the covariance of its rows given the training rows.
#'
#' @param x An `fb_select` object from `fb_select(keep_fold_data = TRUE)`.
#' @param columns Candidate column names, as in the columns of
#'   `fb_select_fold_data(x)[[1]]$x_candidate`. `character(0)` refits the
#'   protected model.
#' @param folds Fold labels to refit, or `NULL` (the default) for every fold.
#'
#' @returns A list with three data frames: `scores` (`fold`, `repeat_id`,
#'   `cluster`, `n`, `mse`, `rmse`, `nlpd`), `coefficients` (`fold`, `role`,
#'   `term`, `estimate`, on each fold's standardised candidate scale) and
#'   `not_identified` (`fold`, `term`).
#'
#' @seealso [fb_select_fold_data()], [fb_select_path()].
#' @family selection fold data
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(site = factor(rep(1:6, each = 5)), x1 = rnorm(30),
#'                   x2 = rnorm(30))
#'   d$y <- 1 + 0.5 * d$x1 + rnorm(6)[d$site] + rnorm(30)
#'   base <- flexybayes(y ~ 1, random = ~site, data = d, backend = "inla",
#'                      verbose = FALSE)
#'   sel <- fb_select(base, ~ x1 + x2, workflow = "unpenalized_benchmark",
#'                    resamples = fb_resamples(unit = ~site, v = 3L),
#'                    seed = 1L, keep_fold_data = TRUE)
#'   fb_select_refit(sel, "x1")$scores
#' }
#' }
#' @export
fb_select_refit <- function(x, columns, folds = NULL) {
  fold_data <- .fb_sel_fold_data_get(x, folds, "fb_select_refit()")
  columns <- as.character(columns)
  declared <- colnames(fold_data$substrates[[1L]]$x_candidate)
  unknown <- setdiff(columns, declared)
  if (length(unknown) > 0L) {
    stop(
      "`columns` names column(s) that are not candidates of the run: ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  scores <- list()
  coefficients <- list()
  not_identified <- list()
  for (label in names(fold_data$substrates)) {
    refit <- .fb_sel_refit_fold(
      fold_data$substrates[[label]],
      columns,
      fold_data$cluster
    )
    scores[[label]] <- refit$scores
    coefficients[[label]] <- refit$coefficients
    not_identified[[label]] <- refit$not_identified
  }
  bind <- function(parts) {
    out <- do.call(rbind, unname(parts))
    rownames(out) <- NULL
    out
  }
  list(
    scores = bind(scores),
    coefficients = bind(coefficients),
    not_identified = bind(not_identified)
  )
}

# .fb_sel_refit_fold() --- one fold's GLS refit of the protected columns and
# `columns`, its held-out cluster scores and its coefficients.
.fb_sel_refit_fold <- function(substrate, columns, cluster) {
  train_rows <- substrate$train_rows
  assess_rows <- substrate$assess_rows
  rows <- c(train_rows, assess_rows)

  # A column not finite on the fold's rows is not identified there.
  finite <- vapply(
    columns,
    function(k) all(is.finite(substrate$x_candidate[rows, k])),
    logical(1L)
  )
  used <- columns[finite]
  x <- cbind(
    substrate$x_protected,
    substrate$x_candidate[, used, drop = FALSE]
  )
  fit <- .fb_sel_gls_fit(
    substrate$blocks$factor_tt,
    x[train_rows, , drop = FALSE],
    substrate$response[train_rows]
  )
  mu <- .fb_sel_fold_predict(
    substrate$blocks,
    fit$beta,
    x[train_rows, , drop = FALSE],
    substrate$response[train_rows],
    x[assess_rows, , drop = FALSE]
  )
  key <- if (is.null(cluster)) {
    rep(substrate$label, length(assess_rows))
  } else {
    as.character(cluster)[assess_rows]
  }
  scored <- .fb_sel_score_fold(
    substrate$assessment$y,
    mu,
    substrate$assessment$v_s_given_t,
    key,
    label = substrate$label,
    repeat_id = substrate$repeat_id,
    model = "refit"
  )
  scores <- as.data.frame(scored$clusters)
  scores$model <- NULL

  estimate <- fit$beta
  estimate[fit$aliased %||% character(0)] <- NA_real_
  protected <- colnames(substrate$x_protected)
  coefficients <- data.frame(
    fold = rep(substrate$label, length(estimate)),
    role = ifelse(names(estimate) %in% protected, "protected", "candidate"),
    term = names(estimate),
    estimate = unname(estimate),
    stringsAsFactors = FALSE
  )
  missing_terms <- c(setdiff(columns, used), fit$aliased %||% character(0))
  list(
    scores = scores,
    coefficients = coefficients,
    not_identified = data.frame(
      fold = rep(substrate$label, length(missing_terms)),
      term = missing_terms,
      stringsAsFactors = FALSE
    )
  )
}
