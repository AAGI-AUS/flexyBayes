# Subsample stability of the penalised optimiser's support.
#
# Meinshausen, N. and Buhlmann, P. (2010). Stability selection. Journal of
# the Royal Statistical Society B 72(4), 417-473.
# doi:10.1111/j.1467-9868.2010.00740.x
#
# B subsamples each hold half of the declared units, drawn without
# replacement. On each, the penalised path is walked with the subsample's
# rows as a fold's training rows, so its covariance, whitening and candidate
# space are learned there and nowhere else, on the same data-free grid of
# s = lambda / lambda_max the cross-validated path walked. A column's
# frequency at a point is the share of the B subsamples whose solution holds
# it; its stability is its largest frequency over the region at or above the
# chosen tuning point; it is stable at or above the cutoff. A group is held
# where any member is.


# .fb_sel_stability_settings() --- the declared B, cutoff and unit, with the
# unit resolved: the declared one, else the resampling unit, else rows.
.fb_sel_stability_settings <- function(decision, resamples) {
  declared <- decision$stability_unit %||% character(0)
  resampling <- resamples$unit %||% character(0)
  unit <- if (length(declared) > 0L) declared else resampling
  list(
    subsamples = as.integer(decision$stability_subsamples %||% 0L),
    cutoff = decision$stability_cutoff %||% 0.75,
    unit = unit,
    unit_source = if (length(declared) > 0L) {
      "declared"
    } else if (length(resampling) > 0L) {
      "resampling unit"
    } else {
      "row"
    }
  )
}


# .fb_sel_stability_units() --- the units a subsample draws half of, or a
# refusal where half of them is fewer than two.
.fb_sel_stability_units <- function(data, unit) {
  .fb_sel_require_columns(data, unit, "stability_unit")
  key <- if (length(unit) > 0L) {
    .fb_sel_key(data, unit)
  } else {
    as.character(seq_len(nrow(data)))
  }
  units <- unique(key)
  half <- as.integer(floor(length(units) / 2))
  if (half < 2L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_not_usable",
      message = paste0(
        "fb_select(): stability subsamples hold half of the ",
        if (length(unit) > 0L) paste(unit, collapse = ":") else "row",
        " units, and there are ", length(units), ", so a subsample would ",
        "hold ", half, ". A subsample of fewer than two units has no ",
        "between-unit variation to learn a covariance from."
      ),
      n_units = length(units)
    ))
  }
  list(key = key, units = units, n_units = length(units),
       subsample_units = half)
}


#' Refuse a stability request on a route that walks no penalised path
#'
#' @param decision An `fb_decision`.
#' @param workflow The resolved workflow.
#' @param penalty The `fb_penalty`, or `NULL`.
#' @param benchmark The `benchmark` flag.
#' @param seed The `seed` of `fb_select()`.
#' @param resamples The `fb_resamples`, for the default unit.
#' @param data The base fit's data, whose units are counted.
#' @return `NULL`, invisibly, or a refusal.
#' @noRd
#' @keywords internal
.fb_sel_require_stability_route <- function(decision,
                                            workflow,
                                            penalty,
                                            benchmark,
                                            seed,
                                            resamples,
                                            data) {
  if ((decision$stability_subsamples %||% 0L) == 0L) {
    return(invisible(NULL))
  }
  runs_path <- identical(workflow, "penalized_optimizer") ||
    (identical(workflow, "bayes") && !is.null(penalty) && isTRUE(benchmark))
  if (!runs_path) {
    stop(.fb_refusal_condition(
      reason_code = "selection_decision_invalid",
      message = paste0(
        "fb_select(): `stability_subsamples = ",
        decision$stability_subsamples,
        "` asks for subsample stability of the penalised optimiser's ",
        "support, and this call walks no penalised path, which does not ",
        "read it. The path runs under workflow = \"penalized_optimizer\", ",
        "or beside a Bayesian answer with a penalty and benchmark = TRUE."
      )
    ))
  }
  if (is.null(seed)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_seed_required",
      message = paste0(
        "fb_select(): `stability_subsamples = ",
        decision$stability_subsamples,
        "` draws its subsamples at random and no `seed` was given. No seed ",
        "is invented here: a default seed is an unreproducible draw that ",
        "looks reproducible."
      )
    ))
  }
  .fb_sel_stability_units(
    data, .fb_sel_stability_settings(decision, resamples)$unit
  )
  invisible(NULL)
}


# .fb_sel_require_stability_family() --- stability reads a support on the s
# coordinate, which a family that does not select, or has no finite null
# anchor, does not have.
.fb_sel_require_stability_family <- function(blocks, decision) {
  if ((decision$stability_subsamples %||% 0L) == 0L) {
    return(invisible(NULL))
  }
  families <- unique(vapply(blocks, function(b) b$family, character(1L)))
  lacking <- families[!vapply(families, function(family) {
    semantics <- .fb_sel_support_semantics(family)
    isTRUE(semantics$selects) && identical(semantics$null_anchor, "finite")
  }, logical(1L))]
  if (length(lacking) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_decision_invalid",
      message = paste0(
        "fb_select(): subsample stability reads the support along ",
        "s = lambda / lambda_max, and `",
        paste(lacking, collapse = "()`, `"),
        "()` has no such path: ridge selects nothing, and the exclusive ",
        "lasso has no finite null anchor. Set stability_subsamples = 0, or ",
        "declare a family with a finite anchor."
      ),
      families = lacking
    ))
  }
  invisible(NULL)
}


#' The stability settings and their cost, as `plan = TRUE` reports them
#'
#' Nothing is fitted. The cost is counted in penalised fits: each path walks
#' the grid and the one canonical point, and each subsample also estimates
#' its own covariance, as each fold does.
#'
#' @param decision An `fb_decision`.
#' @param resamples An `fb_resamples`.
#' @param data The base fit's data.
#' @return A list of the settings, the unit counts and the fit counts.
#' @noRd
#' @keywords internal
.fb_sel_stability_plan <- function(decision, resamples, data) {
  settings <- .fb_sel_stability_settings(decision, resamples)
  n_folds <- if (identical(resamples$scheme, "explicit")) {
    length(unique(resamples$fold_id))
  } else {
    as.integer(resamples$v) * as.integer(resamples$repeats)
  }
  fits_per_path <- .fb_sel_glmnet_decision(NULL)$n_lambda + 1L
  out <- c(settings, list(
    n_units = NA_integer_,
    subsample_units = NA_integer_,
    n_folds = as.integer(n_folds),
    fits_per_path = as.integer(fits_per_path),
    path_fits = settings$subsamples * as.integer(fits_per_path),
    covariance_fits = settings$subsamples,
    cv_path_fits = as.integer(n_folds) * as.integer(fits_per_path),
    relative_cost = settings$subsamples / n_folds
  ))
  if (settings$subsamples > 0L) {
    units <- .fb_sel_stability_units(data, settings$unit)
    out$n_units <- units$n_units
    out$subsample_units <- units$subsample_units
  }
  out
}


# .fb_sel_stability_not_run() --- the record of a run that drew no
# subsamples, so the result says that stability was not asked for rather
# than saying nothing.
.fb_sel_stability_not_run <- function(decision, resamples) {
  settings <- .fb_sel_stability_settings(decision, resamples)
  c(list(ran = FALSE), settings, list(
    reason = "not requested: stability_subsamples = 0"
  ))
}


#' Subsample stability of the penalised optimiser's support
#'
#' @param data,response,protected,candidates,random,prediction,offsets The
#'   inputs the cross-validated path's substrates were built from.
#' @param blocks The resolved penalty blocks.
#' @param path The cross-validated `fb_sel_canonical_path`: its grid, its
#'   chosen point, its decision settings and its declared groups.
#' @param decision The `fb_decision` carrying the stability settings.
#' @param resamples The `fb_resamples`, for the default unit.
#' @param seed The seed the subsamples are drawn under.
#' @return The stability record: settings, draws, the per-subsample
#'   coefficient table, frequencies at every point, the per-column and
#'   per-group summary, the stable sets, the average selected count and the
#'   bound on falsely stable columns.
#' @noRd
#' @keywords internal
.fb_sel_stability <- function(data,
                              response,
                              protected,
                              candidates,
                              random,
                              prediction,
                              offsets,
                              blocks,
                              path,
                              decision,
                              resamples,
                              seed) {
  settings <- .fb_sel_stability_settings(decision, resamples)
  b_count <- settings$subsamples
  units <- .fb_sel_stability_units(data, settings$unit)
  scores <- as.data.frame(path$selection$scores)
  scores <- scores[order(scores$lambda_index), , drop = FALSE]
  if (!any(is.finite(scores$s))) {
    .fb_sel_require_stability_family(
      list(list(family = path$family)), decision
    )
  }
  region_s_min <- path$selection$selected_s
  region <- scores$lambda_index[scores$s >= region_s_min]
  groups <- path$support_by_resolution$groups
  declared <- path$support_by_resolution$history_by_fold$column

  # The draws --------------------------------------------------------------
  draws <- .fb_sel_with_seed(seed, lapply(seq_len(b_count), function(b) {
    sample(units$units, units$subsample_units)
  }))
  groupings <- unique(c(
    unlist(strsplit(random, ":", fixed = TRUE)),
    settings$unit
  ))

  # One penalised path per subsample ----------------------------------------
  coefficients <- vector("list", b_count)
  points <- vector("list", b_count)
  for (b in seq_len(b_count)) {
    held <- units$key %in% draws[[b]]
    split <- list(
      fold_id = b,
      repeat_id = 1L,
      label = paste0("subsample_", b),
      train_rows = which(held),
      assess_rows = which(!held)
    )
    space <- .fb_sel_candidate_space(
      protected, candidates, data, list(split),
      groupings = groupings, route = "benchmark"
    )
    substrate <- .fb_sel_fold_substrate(
      data,
      split,
      response = response,
      protected = protected,
      candidates = candidates,
      random = random,
      prediction = prediction,
      offsets = offsets,
      fold_space = space$folds[[split$label]]
    )
    # The cross-validated path's settings, so the subsample walks its grid.
    # The subsample's own tuning decision is scored on the undrawn units and
    # is not read.
    sub_path <- .fb_sel_glmnet_path(
      list(substrate), blocks, folds = list(split),
      decision = path$decision, seed = seed, cluster = units$key
    )
    sub_scores <- as.data.frame(sub_path$selection$scores)
    sub_scores <- sub_scores[order(sub_scores$lambda_index), , drop = FALSE]
    if (!identical(sub_scores$s, scores$s)) {
      stop(
        "internal contract violation: stability subsample ", b, " walked a ",
        "grid other than the cross-validated path's.",
        call. = FALSE
      )
    }
    coef <- as.data.frame(sub_path$coefficients_original_scale)
    coef <- coef[coef$arm == "benchmark" & coef$role == "candidate", ,
                 drop = FALSE]
    coefficients[[b]] <- data.frame(
      subsample = b,
      lambda_index = as.integer(coef$lambda_index),
      s = scores$s[match(coef$lambda_index, scores$lambda_index)],
      term = coef$term,
      estimate = coef$estimate,
      stringsAsFactors = FALSE
    )
    history <- sub_path$support_by_resolution$history
    history$subsample <- b
    points[[b]] <- history
  }
  coefficients <- do.call(rbind, coefficients)
  points <- do.call(rbind, points)

  # Frequencies at every point, for columns and for groups ------------------
  group_points <- .fb_sel_group_points(points, groups)
  group_points$subsample <- as.integer(sub("^subsample_", "",
                                           group_points$fold))
  at_level <- function(level, name, table) {
    rows <- lapply(unique(name), function(nm) {
      part <- table[name == nm, , drop = FALSE]
      n_available <- tapply(part$available, part$lambda_index, sum)
      n_selected <- tapply(part$in_support, part$lambda_index, sum)
      index <- as.integer(names(n_available))
      data.frame(
        level = level,
        name = nm,
        lambda_index = index,
        s = scores$s[match(index, scores$lambda_index)],
        n_available = as.integer(n_available),
        n_selected = as.integer(n_selected),
        frequency = as.integer(n_selected) / b_count,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }
  frequency_path <- rbind(
    at_level("column", points$column, points),
    at_level("group", group_points$group, group_points)
  )
  rownames(frequency_path) <- NULL

  summarise <- function(level, table, name) {
    part_path <- frequency_path[frequency_path$level == level, , drop = FALSE]
    rows <- lapply(unique(part_path$name), function(nm) {
      part <- part_path[part_path$name == nm, , drop = FALSE]
      subsamples_available <- length(unique(
        table$subsample[name == nm & table$available]
      ))
      largest <- max(part$frequency[part$lambda_index %in% region])
      data.frame(
        level = level,
        name = nm,
        n_available = as.integer(subsamples_available),
        max_frequency = largest,
        selected = largest >= settings$cutoff,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }
  summary <- rbind(
    summarise("column", points, points$column),
    summarise("group", group_points, group_points$group)
  )
  rownames(summary) <- NULL

  # The Meinshausen-Buhlmann bound --------------------------------------------
  in_region <- points[points$lambda_index %in% region & points$in_support, ,
                      drop = FALSE]
  selected_count <- vapply(seq_len(b_count), function(b) {
    length(unique(in_region$column[in_region$subsample == b]))
  }, integer(1L))
  q_hat <- mean(selected_count)
  n_candidates <- length(unique(declared))
  column_summary <- summary[summary$level == "column", , drop = FALSE]
  group_summary <- summary[summary$level == "group", , drop = FALSE]

  c(list(ran = TRUE), settings, list(
    method = paste0(
      "subsample stability selection (Meinshausen and Buhlmann 2010, ",
      "doi:10.1111/j.1467-9868.2010.00740.x)"
    ),
    n_units = units$n_units,
    subsample_units = units$subsample_units,
    seed = seed,
    region = "grid points with s at or above the chosen tuning point",
    region_s_min = region_s_min,
    draws = data.frame(
      subsample = rep(seq_len(b_count), each = units$subsample_units),
      unit = unlist(draws),
      stringsAsFactors = FALSE
    ),
    coefficients = coefficients,
    path = frequency_path,
    summary = summary,
    stable = column_summary$name[column_summary$selected],
    stable_groups = group_summary$name[group_summary$selected],
    q_hat = q_hat,
    n_candidates = n_candidates,
    pfer_bound = q_hat^2 / ((2 * settings$cutoff - 1) * n_candidates),
    pfer_condition = paste0(
      "the expected number of falsely stable columns is at most ",
      "q^2 / ((2 cutoff - 1) p) under the exchangeability condition of ",
      "Theorem 1, with q the average number of columns a subsample selects ",
      "in the region and p the number of penalised columns"
    )
  ))
}


# .fb_sel_stability_record() --- the settings and the answer, without the
# tables, as the certificate carries them.
.fb_sel_stability_record <- function(stability) {
  if (is.null(stability)) {
    return(NULL)
  }
  keep <- c(
    "ran", "subsamples", "cutoff", "unit", "unit_source", "n_units",
    "subsample_units", "seed", "method", "region", "region_s_min", "stable",
    "stable_groups", "q_hat", "n_candidates", "pfer_bound", "reason"
  )
  stability[intersect(keep, names(stability))]
}


# .fb_sel_stability_line() --- the one line the print methods write.
.fb_sel_stability_line <- function(stability) {
  unit <- if (length(stability$unit) > 0L) {
    paste(stability$unit, collapse = ":")
  } else {
    "row"
  }
  paste0(
    "B = ", stability$subsamples, " subsamples of ",
    stability$subsample_units, " of ", stability$n_units, " ", unit,
    " units, cutoff ", format(stability$cutoff), ", over s >= ",
    format(stability$region_s_min, digits = 4L), "; stable: ",
    if (length(stability$stable) == 0L) {
      "none"
    } else {
      paste(stability$stable, collapse = ", ")
    },
    "; at most ", format(stability$pfer_bound, digits = 3L),
    " falsely stable columns expected under exchangeability"
  )
}
