# Grouped ablation of the penalised optimiser --- `fb_select(ablation = TRUE)`.
#
# Under shrinkage a coefficient's size is not a column's importance: the
# penalty sets the size, and a block whose columns share their signal with
# another block's can be shrunk to nothing while the other block carries it.
# The ablation asks the held-out question instead. For each declared block of
# penalised columns, the block's columns are removed from every fold's
# candidate design, the path is walked again on the same folds over the same
# grid of s = lambda / lambda_max and tuned by the same rule, and each
# held-out unit's loss at the ablated path's chosen point is paired with its
# loss at the full path's chosen point.
#
# Four choices, and why.
#
#   Columns, not terms. A block is the set of design columns the penalty
#   declared together (`support$groups`, the map the group profiles read).
#   Removing those columns leaves every other column coded exactly as in the
#   full run; removing a main-effect term from the formula would re-code the
#   interactions that remain.
#
#   The fold's covariance is reused. It is learned from the protected design
#   on the training rows (the REML plug-in) or from a refit of the base
#   model, and neither reads a candidate column, so the ablated path is
#   whitened by exactly the full path's factor and costs no covariance
#   estimate. The fold's candidate space is recomputed, so a column aliased
#   only with the dropped block's columns is estimable again.
#
#   The path is tuned again. The question is what the procedure predicts
#   without the block, so the ablated path chooses its own point by the same
#   rule. A block whose removal leaves no penalised column leaves the
#   protected columns and any unpenalised candidates, fitted by generalised
#   least squares on the same whitening.
#
#   Not nested. Each path's tuning point is chosen on the out-of-fold scores
#   of all folds, including the fold whose units are then scored, so neither
#   loss is a nested estimate. The two share that reuse, and the difference
#   is read between them; the record says so.


# .fb_sel_ablation_check() --- the request, validated before anything runs.
.fb_sel_ablation_check <- function(ablation) {
  if (!is.logical(ablation) || length(ablation) != 1L || is.na(ablation)) {
    stop("`ablation` must be TRUE or FALSE.", call. = FALSE)
  }
  ablation
}


#' Refuse an ablation on a route that walks no penalised path
#'
#' @param ablation The validated flag.
#' @param workflow The resolved workflow.
#' @return `NULL`, invisibly, or a `selection_ablation_route_undefined`
#'   refusal.
#' @noRd
#' @keywords internal
.fb_sel_ablation_require_route <- function(ablation, workflow) {
  if (!isTRUE(ablation) || identical(workflow, "penalized_optimizer")) {
    return(invisible(NULL))
  }
  stop(.fb_refusal_condition(
    reason_code = "selection_ablation_route_undefined",
    message = paste0(
      "fb_select(): ablation = TRUE drops each declared block of the ",
      "penalised optimiser's path and walks the path again, and ",
      "workflow = \"", workflow, "\" walks no penalised path as its answer. ",
      "Nothing was fitted. Use workflow = \"penalized_optimizer\" with a ",
      "penalty."
    ),
    workflow = workflow
  ))
}


# .fb_sel_ablation_contract() --- the contract every path of one ablation
# shares: the fold of every row, the scoring unit, the target, the response
# values and the loss. The design is deliberately the one-column intercept,
# as in the value contract, because an ablation varies the predictor set.
.fb_sel_ablation_contract <- function(folds, data, response, offsets,
                                      prediction, metric, scoring_unit) {
  .fb_select_score_contract(
    folds = folds,
    data = data,
    response = response,
    prediction = prediction,
    metric = metric,
    design = matrix(1, nrow(data), 1L, dimnames = list(NULL, "(Intercept)")),
    learned_on = "analysis_rows",
    stage = "ablation",
    scoring_unit = scoring_unit,
    offsets = offsets,
    candidate_space = NULL
  )
}


# .fb_sel_ablation_folds() --- the folds an ablated path was scored on, read
# off its substrates, with the unit map of the full run's folds.
.fb_sel_ablation_folds <- function(substrates, folds) {
  out <- lapply(unname(substrates), function(s) {
    list(
      fold_id = s$fold_id,
      repeat_id = s$repeat_id,
      label = s$label,
      train_rows = s$train_rows,
      assess_rows = s$assess_rows
    )
  })
  attributes(out) <- attributes(folds)
  out
}


# .fb_sel_ablation_substrate() --- one fold's substrate without the dropped
# columns. The covariance, whitening and response are the fold's own; the
# candidate space is recomputed over the columns that remain.
.fb_sel_ablation_substrate <- function(substrate, drop, data, protected,
                                       candidates, groupings) {
  split <- list(
    fold_id = substrate$fold_id,
    repeat_id = substrate$repeat_id,
    label = substrate$label,
    train_rows = substrate$train_rows,
    assess_rows = substrate$assess_rows
  )
  protected_design <- .fb_sel_design(protected, data, split$train_rows)
  candidate_design <- .fb_sel_design(candidates, data, split$train_rows)
  keep <- !colnames(candidate_design$X) %in% drop
  candidate_design$X <- candidate_design$X[, keep, drop = FALSE]
  candidate_design$assign <- candidate_design$assign[keep]
  space <- .fb_sel_fold_candidate_space(
    protected_design, candidate_design, data, split,
    groupings = groupings, resolution = "column"
  )
  columns <- setdiff(colnames(substrate$x_candidate), drop)
  substrate$x_candidate <- substrate$x_candidate[, columns, drop = FALSE]
  substrate$candidate_kept <- space$kept_columns
  substrate$candidate_space <- space
  substrate$learned$preprocessing_centre <-
    substrate$learned$preprocessing_centre[columns]
  substrate$learned$preprocessing_scale <-
    substrate$learned$preprocessing_scale[columns]
  # The pilot fit was over every declared column; the ablated path does not
  # read it.
  substrate$learned$pilot_beta <- NULL
  substrate
}


# .fb_sel_ablation_blocks() --- the resolved penalty blocks without the
# dropped columns, each naming its columns, and a block left with none
# removed.
.fb_sel_ablation_blocks <- function(blocks, candidate_names, drop) {
  out <- lapply(blocks, function(block) {
    columns <- block$columns
    if (!is.character(columns)) {
      columns <- candidate_names[columns]
    }
    block$columns <- setdiff(columns, drop)
    block
  })
  out[vapply(out, function(b) length(b$columns) > 0L, logical(1L))]
}


# .fb_sel_ablation_protected() --- every fold's generalised-least-squares fit
# of the columns an ablation leaves when it removes every penalised one: the
# protected columns and any unpenalised candidate the fold keeps.
.fb_sel_ablation_protected <- function(substrates, cluster) {
  clusters <- list()
  fitted <- list()
  for (s in substrates) {
    kept <- which(colnames(s$x_candidate) %in% s$candidate_kept)
    x <- cbind(s$x_protected, s$x_candidate[, kept, drop = FALSE])
    fit <- .fb_sel_gls_fit(
      s$blocks$factor_tt,
      x[s$train_rows, , drop = FALSE],
      s$response[s$train_rows]
    )
    mu <- .fb_sel_fold_predict(
      s$blocks,
      fit$beta,
      x[s$train_rows, , drop = FALSE],
      s$response[s$train_rows],
      x[s$assess_rows, , drop = FALSE]
    )
    key <- if (is.null(cluster)) {
      rep(s$label, length(s$assess_rows))
    } else {
      as.character(cluster)[s$assess_rows]
    }
    scored <- .fb_sel_score_fold(
      s$assessment$y, mu, s$assessment$v_s_given_t, key,
      label = s$label, repeat_id = s$repeat_id, model = "protected"
    )
    clusters[[length(clusters) + 1L]] <- scored$clusters
    fitted[[length(fitted) + 1L]] <- data.frame(
      fold = s$label,
      repeat_id = s$repeat_id,
      row = s$assess_rows,
      unit = key,
      observed = as.numeric(s$assessment$y),
      fitted = as.numeric(mu),
      stringsAsFactors = FALSE
    )
  }
  list(
    clusters = data.table::rbindlist(clusters),
    fitted = do.call(rbind, fitted)
  )
}


# .fb_sel_ablation_unit_loss() --- each held-out unit's loss at the metric
# the path was tuned on: the unit's score in each repeat, averaged over
# repeats, in radix order of the units.
.fb_sel_ablation_unit_loss <- function(clusters, metric) {
  per_unit <- .fb_sel_aggregate_scores(
    clusters, metric = metric, unit = "cluster"
  )$per_unit
  per_unit[sort(names(per_unit), method = "radix")]
}


#' Drop each declared block in turn and walk the path again
#'
#' @param data,response,protected,candidates,random,offsets The inputs the
#'   full path's substrates were built from.
#' @param blocks The resolved penalty blocks.
#' @param path The full `fb_sel_canonical_path`: its decision settings, its
#'   chosen point, its scores and its declared groups.
#' @param substrates The full path's fold substrates.
#' @param folds The folds the full path ran on.
#' @param resamples The `fb_resamples`, for the grouping columns.
#' @param prediction The prediction target.
#' @param cluster The scoring-unit key of every row, or `NULL` for the fold.
#' @param scoring_unit The scoring unit's name.
#' @param seed The run's seed.
#' @return The ablation record: the contract, one row per declared block,
#'   the paired unit table, the ablated paths' per-row predictions and
#'   chosen supports.
#' @noRd
#' @keywords internal
.fb_sel_ablation <- function(data,
                             response,
                             protected,
                             candidates,
                             random,
                             offsets,
                             blocks,
                             path,
                             substrates,
                             folds,
                             resamples,
                             prediction,
                             cluster,
                             scoring_unit,
                             seed) {
  metric <- path$selection$metric
  groups <- path$support_by_resolution$groups
  names_vec <- unique(groups$group)
  candidate_names <- colnames(substrates[[1L]]$x_candidate)
  groupings <- unique(c(
    unlist(strsplit(random, ":", fixed = TRUE)),
    resamples$unit
  ))

  contract <- .fb_sel_ablation_contract(
    folds, data, response, offsets, prediction, metric, scoring_unit
  )
  .fb_select_verify_contract(contract, "ablation")
  index <- path$selection$selected_lambda_index
  full_model <- paste0("benchmark_", index)
  full_units <- .fb_sel_ablation_unit_loss(
    path$clusters_by_model[[full_model]], metric
  )

  one_block <- function(group) {
    drop <- groups$column[groups$group == group]
    ablated <- lapply(substrates, .fb_sel_ablation_substrate,
      drop = drop, data = data, protected = protected,
      candidates = candidates, groupings = groupings
    )
    remaining <- .fb_sel_ablation_blocks(blocks, candidate_names, drop)
    run <- if (length(remaining) == 0L) {
      protected_fit <- .fb_sel_ablation_protected(ablated, cluster)
      list(
        clusters = protected_fit$clusters,
        fitted = protected_fit$fitted,
        selected_s = NA_real_,
        selected_lambda_index = NA_integer_,
        curve = NULL,
        supports = stats::setNames(
          lapply(ablated, function(s) character(0)),
          vapply(ablated, function(s) s$label, character(1L))
        )
      )
    } else {
      ablated_path <- .fb_sel_glmnet_path(
        ablated, remaining, folds = folds, decision = path$decision,
        seed = seed, cluster = cluster
      )
      if (!identical(ablated_path$selection$scores$s,
                     path$selection$scores$s)) {
        stop(
          "internal contract violation: the path without ", group,
          " walked a grid other than the full path's.",
          call. = FALSE
        )
      }
      chosen <- ablated_path$selection$selected_lambda_index
      model <- paste0("benchmark_", chosen)
      fitted <- as.data.frame(ablated_path$fitted_values)
      fitted <- fitted[
        fitted$arm == "benchmark" & fitted$lambda_index == chosen, ,
        drop = FALSE
      ]
      list(
        clusters = ablated_path$clusters_by_model[[model]],
        fitted = fitted[, c("fold", "repeat_id", "row", "unit", "observed",
                            "fitted")],
        selected_s = ablated_path$selection$selected_s,
        selected_lambda_index = chosen,
        # The ablated path's own tuning curve, which the rule chose on.
        curve = as.data.frame(ablated_path$selection$scores)[,
          c("lambda_index", "s", "score", "se", "n_units", "complexity")
        ],
        # A family that sets no coefficient to zero has no support to
        # report; its loss is read all the same.
        supports = if (isTRUE(ablated_path$support_by_resolution$selects)) {
          .fb_sel_value_supports(
            list(
              family = ablated_path$family,
              support = ablated_path$support_by_resolution,
              selection = ablated_path$selection
            ),
            folds
          )
        } else {
          list()
        }
      )
    }
    run$contract <- .fb_sel_ablation_contract(
      .fb_sel_ablation_folds(ablated, folds), data, response, offsets,
      prediction, metric, scoring_unit
    )
    run
  }

  results <- stats::setNames(lapply(names_vec, function(group) {
    tryCatch(
      one_block(group),
      flexybayes_refusal = function(e) e
    )
  }), names_vec)

  block_rows <- list()
  unit_rows <- list()
  row_rows <- list()
  support_rows <- list()
  curve_rows <- list()
  for (group in names_vec) {
    res <- results[[group]]
    members <- groups$column[groups$group == group]
    ran <- !inherits(res, "condition")
    row <- data.frame(
      block = group,
      rule = groups$rule[match(group, groups$group)],
      n_columns = length(members),
      columns = paste(members, collapse = ", "),
      ran = ran,
      reason = if (ran) {
        NA_character_
      } else {
        paste0(res$reason_code %||% "error", ": ", conditionMessage(res))
      },
      selected_s = NA_real_,
      loss_full = mean(full_units),
      loss_ablated = NA_real_,
      loss_difference = NA_real_,
      loss_se = NA_real_,
      n_units = NA_integer_,
      contract_hash = NA_character_,
      stringsAsFactors = FALSE
    )
    if (ran) {
      .fb_select_verify_contract(res$contract, "ablation")
      if (!identical(res$contract$hash, contract$hash)) {
        stop(
          "internal contract violation: the path without ", group, " was ",
          "scored on a contract other than the full path's; they differ in ",
          paste(.fb_select_contract_differences(res$contract, contract),
                collapse = ", "),
          ".",
          call. = FALSE
        )
      }
      ablated_units <- .fb_sel_ablation_unit_loss(res$clusters, metric)
      if (!identical(names(ablated_units), names(full_units))) {
        stop(
          "internal contract violation: the path without ", group,
          " scored units other than the full path's.",
          call. = FALSE
        )
      }
      difference <- ablated_units - full_units
      row$selected_s <- res$selected_s
      row$loss_ablated <- mean(ablated_units)
      row$loss_difference <- mean(difference)
      row$loss_se <- stats::sd(difference) / sqrt(length(difference))
      row$n_units <- length(difference)
      row$contract_hash <- res$contract$hash
      unit_rows[[length(unit_rows) + 1L]] <- data.frame(
        block = group,
        unit = names(full_units),
        loss_full = unname(full_units),
        loss_ablated = unname(ablated_units),
        loss_difference = unname(difference),
        stringsAsFactors = FALSE
      )
      fitted <- res$fitted
      fitted <- fitted[order(fitted$repeat_id, fitted$row, method = "radix"),
                       , drop = FALSE]
      row_rows[[length(row_rows) + 1L]] <- data.frame(
        block = group,
        repeat_id = as.integer(fitted$repeat_id),
        fold = as.character(fitted$fold),
        row = as.integer(fitted$row),
        unit = as.character(fitted$unit),
        observed = as.numeric(fitted$observed),
        predicted = as.numeric(fitted$fitted),
        lpd = NA_real_,
        stringsAsFactors = FALSE
      )
      if (!is.null(res$curve)) {
        curve_rows[[length(curve_rows) + 1L]] <- cbind(
          block = group, res$curve, stringsAsFactors = FALSE
        )
      }
      labels <- names(res$supports)
      held <- res$supports
      support_rows[[length(support_rows) + 1L]] <- data.frame(
        block = rep(group, sum(lengths(held))),
        fold = rep(labels, lengths(held)),
        column = as.character(unlist(held, use.names = FALSE)),
        stringsAsFactors = FALSE
      )
    }
    block_rows[[length(block_rows) + 1L]] <- row
  }

  bind <- function(parts, empty) {
    if (length(parts) == 0L) empty else do.call(rbind, parts)
  }
  out <- list(
    metric = metric,
    scoring_unit = scoring_unit,
    nested = FALSE,
    tuning = paste0(
      "each path's tuning point is chosen on the out-of-fold scores of all ",
      "folds, including the fold whose units are then scored, so neither ",
      "loss is a nested estimate; the full and the ablated paths share that ",
      "reuse"
    ),
    formula = paste0(
      "loss_difference = mean over the U held-out units u of d_u, where d_u ",
      "is the unit's ", metric, " under the path without the block minus ",
      "under the full path, each averaged over repeats; loss_se = sd(d_u) / ",
      "sqrt(U). With the Bayesian refit, elpd_difference = sum over u of ",
      "elpd_u(full) - elpd_u(without the block); elpd_se = sqrt(U) * ",
      "sd(elpd_u(full) - elpd_u(without the block))"
    ),
    full_selected_s = path$selection$selected_s,
    contract = contract,
    blocks = bind(block_rows, NULL),
    units = bind(unit_rows, .fb_sel_ablation_empty_units()),
    rows = bind(row_rows, .fb_sel_ablation_empty_rows()),
    supports = bind(support_rows, data.frame(
      block = character(0), fold = character(0), column = character(0),
      stringsAsFactors = FALSE
    )),
    paths = bind(curve_rows, data.frame(
      block = character(0), lambda_index = integer(0), s = numeric(0),
      score = numeric(0), se = numeric(0), n_units = integer(0),
      complexity = numeric(0), stringsAsFactors = FALSE
    ))
  )
  out$blocks$elpd_difference <- NA_real_
  out$blocks$elpd_se <- NA_real_
  out$blocks$elpd_contract_hash <- NA_character_
  # Filled by the Bayesian refit, one row per block, repeat and held-out
  # unit of the value contract.
  out$value_units <- data.frame(
    block = character(0), repeat_id = integer(0), fold = character(0),
    unit = character(0), n_rows = integer(0), elpd_full = numeric(0),
    elpd_ablated = numeric(0), elpd_difference = numeric(0),
    stringsAsFactors = FALSE
  )
  out$value_supports <- data.frame(
    block = character(0), fold = character(0), column = character(0),
    stringsAsFactors = FALSE
  )
  rownames(out$units) <- NULL
  rownames(out$rows) <- NULL
  rownames(out$supports) <- NULL
  rownames(out$paths) <- NULL
  out
}


#' Refit each ablated path's frozen supports and pair the held-out density
#'
#' The Bayesian refit of `fb_select(refit = )`, run once more for each block
#' the loss ablation ran on, with the supports the path without the block
#' chose on each fold. The refits use the full valuation's folds, seeds and
#' value contract, so their units pair with its units one by one.
#'
#' @param ablation The loss ablation record.
#' @param value The full valuation.
#' @param ... Passed to `.fb_sel_value()`: the base, candidates, resamples,
#'   prediction, folds, request, engine and seed.
#' @return The ablation record with the density columns filled.
#' @noRd
#' @keywords internal
.fb_sel_ablation_value <- function(ablation, value, ...) {
  ran <- ablation$blocks$block[ablation$blocks$ran]
  for (group in ran) {
    part <- ablation$supports[ablation$supports$block == group, ,
                              drop = FALSE]
    labels <- unique(value$folds$fold)
    supports <- stats::setNames(lapply(labels, function(label) {
      part$column[part$fold == label]
    }), labels)
    ablated <- .fb_sel_value(supports = supports, budget = NULL, ...)
    paired <- .fb_sel_value_pair_records(value, ablated)
    at <- ablation$blocks$block == group
    ablation$blocks$elpd_difference[at] <- paired$difference
    ablation$blocks$elpd_se[at] <- paired$se
    ablation$blocks$elpd_contract_hash[at] <- ablated$contract$hash
    ablation$value_units <- rbind(ablation$value_units, data.frame(
      block = group,
      repeat_id = paired$units$repeat_id,
      fold = paired$units$fold,
      unit = paired$units$unit,
      n_rows = paired$units$n_rows,
      elpd_full = paired$units$elpd,
      elpd_ablated = paired$units$elpd_against,
      elpd_difference = paired$units$difference,
      stringsAsFactors = FALSE
    ))
    # The columns each fold's refit held, read off the refit itself.
    ablation$value_supports <- rbind(ablation$value_supports, data.frame(
      block = rep(group, NROW(ablated$supports)),
      fold = as.character(ablated$supports$fold),
      column = as.character(ablated$supports$column),
      stringsAsFactors = FALSE
    ))
    in_rows <- which(ablation$rows$block == group)
    key <- paste(ablation$rows$repeat_id[in_rows], ablation$rows$row[in_rows])
    lpd_key <- paste(ablated$rows$repeat_id, ablated$rows$row)
    ablation$rows$lpd[in_rows] <- ablated$rows$lpd[match(key, lpd_key)]
  }
  ablation
}


#' The ablation's blocks and cost, as `plan = TRUE` reports them
#'
#' Nothing is fitted. The declared blocks are resolved as the path resolves
#' them, over the candidate design of every row. Each block's path walks the
#' grid and the one canonical point on every fold; a block whose removal
#' leaves no penalised column costs one least-squares fit per fold instead.
#' The folds' covariance is reused, and with the Bayesian refit each block
#' costs one engine refit per fold.
#'
#' @param resolved The resolved penalty.
#' @param candidates The candidate formula.
#' @param data The base fit's data.
#' @param resamples The `fb_resamples`.
#' @param refit Whether the Bayesian refit was asked for.
#' @return A list of the blocks and the fit counts.
#' @noRd
#' @keywords internal
.fb_sel_ablation_plan <- function(resolved, candidates, data, resamples,
                                  refit) {
  design <- .fb_sel_design(candidates, data)
  keep <- colnames(design$X) != "(Intercept)"
  candidate_names <- colnames(design$X)[keep]
  candidate_terms <- attr(design$terms, "term.labels")[design$assign[keep]]
  penalty <- .fb_sel_glmnet_blocks(
    resolved$blocks, candidate_names, n = nrow(data),
    candidate_terms = candidate_terms
  )
  groups <- .fb_sel_declared_groups(penalty)
  names_vec <- unique(groups$group)
  leaves_path <- vapply(names_vec, function(group) {
    drop <- groups$column[groups$group == group]
    length(.fb_sel_ablation_blocks(
      resolved$blocks, candidate_names, drop
    )) > 0L
  }, logical(1L))
  n_folds <- if (identical(resamples$scheme, "explicit")) {
    length(unique(resamples$fold_id))
  } else {
    as.integer(resamples$v) * as.integer(resamples$repeats)
  }
  fits_per_path <- as.integer(.fb_sel_glmnet_decision(NULL)$n_lambda + 1L)
  list(
    blocks = names_vec,
    n_folds = as.integer(n_folds),
    fits_per_path = fits_per_path,
    path_fits = as.integer(sum(leaves_path)) * as.integer(n_folds) *
      fits_per_path,
    protected_fits = as.integer(sum(!leaves_path)) * as.integer(n_folds),
    covariance_fits = 0L,
    engine_refits = if (isTRUE(refit)) {
      length(names_vec) * as.integer(n_folds)
    } else {
      0L
    },
    cv_path_fits = as.integer(n_folds) * fits_per_path
  )
}


# .fb_sel_ablation_lines() --- what the print method writes.
.fb_sel_ablation_lines <- function(ablation, digits) {
  head <- paste0(
    "each declared block dropped in turn and the path walked again on the ",
    "same folds, tuned by the same rule; held-out ", ablation$metric,
    " paired by ", ablation$scoring_unit, " (without the block minus with ",
    "it); not nested, each path's tuning point is chosen on the scores of ",
    "all folds"
  )
  body <- vapply(seq_len(nrow(ablation$blocks)), function(i) {
    b <- ablation$blocks[i, ]
    label <- paste0(b$block, " (", b$columns, ")")
    if (!isTRUE(b$ran)) {
      return(paste0(label, ": not run -- ", .fb_squish(b$reason)))
    }
    paste0(
      label, ": loss ", format(round(b$loss_difference, digits)),
      " (se ", format(round(b$loss_se, digits)), ")",
      if (is.na(b$elpd_difference)) {
        ""
      } else {
        paste0(
          ", elpd ", format(round(b$elpd_difference, digits)),
          " (se ", format(round(b$elpd_se, digits)), ")"
        )
      }
    )
  }, character(1L))
  c(head, body)
}


# --- the accessors' empty shapes -------------------------------------- #

.fb_sel_ablation_empty_units <- function() {
  data.frame(
    block = character(0), unit = character(0), loss_full = numeric(0),
    loss_ablated = numeric(0), loss_difference = numeric(0),
    stringsAsFactors = FALSE
  )
}

.fb_sel_ablation_empty_rows <- function() {
  data.frame(
    block = character(0), repeat_id = integer(0), fold = character(0),
    row = integer(0), unit = character(0), observed = numeric(0),
    predicted = numeric(0), lpd = numeric(0), stringsAsFactors = FALSE
  )
}
