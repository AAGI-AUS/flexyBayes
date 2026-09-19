# Entrance profiles of the penalised optimiser's path, for columns and for
# declared groups of columns, and the entrance heat map.
#
# The fold-level record is `.fb_sel_glmnet_support_history()`: for each fold
# and declared column, whether the fold could estimate the column, whether its
# solution holds it at each point of the grid on s = lambda / lambda_max, and
# where it first enters. This file summarises that record across the folds.
# `tidy()`, the entrance path and the heat map all read the summary written
# here, so the table and both displays place a column at the same point.


# .fb_sel_path_transitions() --- where one membership sequence enters, leaves
# and returns, read in path order from the null anchor downwards.
#
# The point before the first is out by construction, so an entry at the first
# point counts as an entry. `first_entry_s` is the first point held,
# `last_exit_s` the last point at which the sequence is out again having been
# held, and `n_reentries` the number of entries after the first.
.fb_sel_path_transitions <- function(in_support, s) {
  n_points <- length(in_support)
  previous <- c(FALSE, in_support[-n_points])
  entries <- which(in_support & !previous)
  exits <- which(!in_support & previous)
  list(
    first_entry_s = if (length(entries) > 0L) {
      s[[entries[[1L]]]]
    } else {
      NA_real_
    },
    last_exit_s = if (length(exits) > 0L) {
      s[[exits[[length(exits)]]]]
    } else {
      NA_real_
    },
    n_reentries = max(0L, length(entries) - 1L)
  )
}


# .fb_sel_declared_groups() --- the group each penalised column belongs to,
# as the penalty declared it.
#
# A grouped family carries its own grouping: the model term under
# `group = "term"`, the specification under `group = "block"`. Any other
# family is grouped by the penalty specification that claimed the column,
# named `block_<k>` in the order the specifications were written, which is
# the naming the grouped `block` rule already uses.
.fb_sel_declared_groups <- function(resolved) {
  if (!is.null(resolved$grouping)) {
    declared <- resolved$grouping$declared
    return(data.frame(
      column = declared$column,
      group = declared$group,
      rule = resolved$grouping$rule,
      stringsAsFactors = FALSE
    ))
  }
  map <- do.call(rbind, lapply(seq_along(resolved$blocks), function(i) {
    data.frame(
      column = resolved$blocks[[i]]$column_names,
      group = paste0("block_", i),
      rule = "block",
      stringsAsFactors = FALSE
    )
  }))
  map <- map[match(resolved$column_names, map$column), , drop = FALSE]
  rownames(map) <- NULL
  map
}


# .fb_sel_group_points() --- the column membership at every point, read at
# the group's resolution: a group is available on a fold where any member is,
# and held at a point where any member is.
.fb_sel_group_points <- function(points, groups) {
  empty <- data.frame(
    fold = character(0), lambda_index = integer(0), s = numeric(0),
    group = character(0), available = logical(0), in_support = logical(0),
    n_available = integer(0), stringsAsFactors = FALSE
  )
  if (is.null(groups) || is.null(points) || nrow(points) == 0L) {
    return(empty)
  }
  points$group <- groups$group[match(points$column, groups$column)]
  points <- points[!is.na(points$group), , drop = FALSE]
  if (nrow(points) == 0L) {
    return(empty)
  }
  key <- paste(points$fold, points$lambda_index, points$group, sep = "\r")
  first <- !duplicated(key)
  out <- points[first, c("fold", "lambda_index", "s", "group"), drop = FALSE]
  out$available <- as.logical(tapply(points$available, key, any)[key[first]])
  out$in_support <- as.logical(tapply(points$in_support, key, any)[key[first]])
  out$n_available <- as.integer(tapply(points$available, key, sum)[key[first]])
  order_groups <- unique(groups$group)
  out <- out[order(
    match(out$fold, unique(points$fold)),
    match(out$group, order_groups),
    out$lambda_index
  ), , drop = FALSE]
  rownames(out) <- NULL
  out
}


# .fb_sel_group_history_by_fold() --- one row per fold and group, carrying
# the fields `.fb_sel_glmnet_support_history()` carries per column.
.fb_sel_group_history_by_fold <- function(points, groups) {
  group_points <- .fb_sel_group_points(points, groups)
  empty <- data.frame(
    fold = character(0), group = character(0), available = logical(0),
    n_columns_available = integer(0), n_points = integer(0),
    n_in_support = integer(0), first_entry_s = numeric(0),
    last_exit_s = numeric(0), n_reentries = integer(0),
    stringsAsFactors = FALSE
  )
  if (nrow(group_points) == 0L) {
    return(empty)
  }
  key <- paste(group_points$fold, group_points$group, sep = "\r")
  rows <- lapply(unique(key), function(k) {
    part <- group_points[key == k, , drop = FALSE]
    part <- part[order(part$lambda_index), , drop = FALSE]
    moves <- .fb_sel_path_transitions(part$in_support, part$s)
    data.frame(
      fold = part$fold[[1L]],
      group = part$group[[1L]],
      available = part$available[[1L]],
      n_columns_available = part$n_available[[1L]],
      n_points = nrow(part),
      n_in_support = sum(part$in_support),
      first_entry_s = moves$first_entry_s,
      last_exit_s = moves$last_exit_s,
      n_reentries = moves$n_reentries,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# .fb_sel_entry_at() --- the point on s by which at least a fraction `p` of
# the folds that could estimate a column had admitted it.
#
# `first` holds one first entry per available fold, NA for a fold on which
# the column never entered. Read in path order, from s = 1 downwards, the
# answer is the ceiling(p k)-th largest entry over k folds, so it sits on a
# grid point. A fold that never admitted the column counts against it, so
# with fewer than ceiling(p k) entries the answer is NA. At p = 0.5 this is
# the median entry the entrance path labels.
.fb_sel_entry_at <- function(first, p) {
  k <- length(first)
  if (k == 0L) {
    return(NA_real_)
  }
  needed <- ceiling(p * k)
  entered <- sort(first[!is.na(first)], decreasing = TRUE)
  if (length(entered) < needed) {
    return(NA_real_)
  }
  entered[[needed]]
}


#' The entrance profile of every column or group across the folds
#'
#' Folds available and folds entered are counted apart, because a column a
#' fold could not estimate was never on that fold's path. The median and the
#' quartiles are read in path order among the available folds, so
#' `entry_s_q25` is the larger `s` (a quarter of the folds had admitted the
#' column by then) and `entry_s_q75` the smaller.
#'
#' @param by_fold A per-fold table with `available`, `n_in_support`,
#'   `first_entry_s`, `last_exit_s` and `n_reentries`.
#' @param key The column of `by_fold` naming the column or group.
#' @return One row per name: `name`, `folds_available`, `folds_entered`,
#'   `path_entry_s` (the earliest entry on any fold), `path_exit_s`,
#'   `path_reentries`, `entry_s_median`, `entry_s_q25`, `entry_s_q75` and
#'   `entry_s_iqr`.
#' @noRd
#' @keywords internal
.fb_select_entrance_profile <- function(by_fold, key = "column") {
  empty <- data.frame(
    name = character(0), folds_available = integer(0),
    folds_entered = integer(0), path_entry_s = numeric(0),
    path_exit_s = numeric(0), path_reentries = integer(0),
    entry_s_median = numeric(0), entry_s_q25 = numeric(0),
    entry_s_q75 = numeric(0), entry_s_iqr = numeric(0),
    stringsAsFactors = FALSE
  )
  if (!is.data.frame(by_fold) || nrow(by_fold) == 0L) {
    return(empty)
  }
  rows <- lapply(unique(by_fold[[key]]), function(name) {
    part <- by_fold[by_fold[[key]] == name, , drop = FALSE]
    entered <- part$n_in_support > 0L
    first <- part$first_entry_s[part$available]
    first[!entered[part$available]] <- NA_real_
    exits <- part$last_exit_s[entered]
    q25 <- .fb_sel_entry_at(first, 0.25)
    q75 <- .fb_sel_entry_at(first, 0.75)
    data.frame(
      name = name,
      folds_available = sum(part$available),
      folds_entered = sum(entered),
      path_entry_s = if (any(entered)) {
        max(part$first_entry_s[entered])
      } else {
        NA_real_
      },
      path_exit_s = if (any(!is.na(exits))) {
        min(exits, na.rm = TRUE)
      } else {
        NA_real_
      },
      path_reentries = if (any(entered)) {
        max(part$n_reentries[entered])
      } else {
        NA_integer_
      },
      entry_s_median = .fb_sel_entry_at(first, 0.5),
      entry_s_q25 = q25,
      entry_s_q75 = q75,
      entry_s_iqr = q25 - q75,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# The viridis sequence at eleven steps, as grDevices::hcl.colors(11,
# "viridis") gives it: ordered in lightness, so it reads in grey and to the
# common colour-vision deficiencies. Written out so the display needs no
# graphics-device package.
.FB_SEL_HEAT_COLOURS <- c(
  "#4B0055", "#44286E", "#274983", "#006791", "#008298", "#009B95",
  "#00B28A", "#00C475", "#7ED357", "#C2DE34", "#FDE333"
)


#' The entrance heat map
#'
#' One row per candidate column or declared group and one column per point of
#' the tuning grid, strong penalty on the left. With `source = "folds"` a
#' cell is the fraction of the folds that could estimate the row whose
#' solution holds it at that point, and a marker sits at the row's median
#' entry. With `source = "stability"` a cell is the fraction of the stability
#' subsamples whose solution holds it, and a stable row is starred.
#'
#' @param x An `fb_select` object.
#' @param level `"column"` or `"group"`.
#' @param source `"folds"` or `"stability"`.
#' @param max_rows The largest number of rows drawn, earliest entries first.
#' @param label_cex Character expansion of the row labels.
#' @param ... Passed to `graphics::plot.window()`.
#' @return Invisibly, `values` (the matrix drawn, all rows), `s`, `rows`
#'   (name, median entry, largest stability frequency, whether stable),
#'   `level`, `source` and `selected_s`; `NULL` with a message where there is
#'   nothing to draw.
#' @noRd
#' @keywords internal
.fb_select_plot_entrance_heatmap <- function(x,
                                             level = c("column", "group"),
                                             source = c("folds", "stability"),
                                             max_rows = 40L,
                                             label_cex = 0.7,
                                             ...) {
  level <- match.arg(level)
  source <- match.arg(source)
  arm <- x$benchmark$glmnet
  scores <- arm$selection$scores

  # What there is to draw --------------------------------------------------
  if (!is.data.frame(scores) || nrow(scores) == 0L) {
    return(.fb_select_plot_unavailable(
      "the entrance heat map",
      "no penalised optimiser path was scored out of fold"
    ))
  }
  if (!any(is.finite(scores$s))) {
    return(.fb_select_plot_unavailable(
      "the entrance heat map",
      paste0(
        "this penalty has no finite null anchor, so its path has no ",
        "s = lambda / lambda_max coordinate"
      )
    ))
  }
  points <- arm$support$history
  if (!is.data.frame(points) || nrow(points) == 0L) {
    return(.fb_select_plot_unavailable(
      "the entrance heat map",
      "this penalty has no finite coefficient-support path"
    ))
  }
  stability <- arm$stability
  if (identical(source, "stability") && !isTRUE(stability$ran)) {
    return(.fb_select_plot_unavailable(
      "the stability heat map",
      paste0(
        "no stability subsamples were drawn; ask for them with ",
        "fb_decision(stability_subsamples = )"
      )
    ))
  }
  groups <- arm$support$groups
  if (identical(level, "group") && !is.data.frame(groups)) {
    return(.fb_select_plot_unavailable(
      "the group heat map",
      "this result records no declared grouping of its columns"
    ))
  }

  # The values ------------------------------------------------------------
  s_grid <- scores$s[order(scores$lambda_index)]
  n_points <- length(s_grid)
  if (identical(level, "column")) {
    by_fold <- arm$support$history_by_fold
    profile <- .fb_select_entrance_profile(by_fold, key = "column")
    cell_points <- points
    cell_points$name <- cell_points$column
  } else {
    by_fold <- arm$support$group_history_by_fold
    profile <- .fb_select_entrance_profile(by_fold, key = "group")
    cell_points <- .fb_sel_group_points(points, groups)
    cell_points$name <- cell_points$group
  }
  names_vec <- profile$name
  values <- matrix(
    NA_real_, nrow = length(names_vec), ncol = n_points,
    dimnames = list(names_vec, NULL)
  )
  if (identical(source, "folds")) {
    for (name in names_vec) {
      part <- cell_points[cell_points$name == name, , drop = FALSE]
      n_available <- tapply(part$available, part$lambda_index, sum)
      n_held <- tapply(part$in_support, part$lambda_index, sum)
      index <- as.integer(names(n_available))
      ratio <- ifelse(n_available > 0L, n_held / n_available, NA_real_)
      values[name, index] <- as.numeric(ratio)
    }
  } else {
    path <- stability$path[stability$path$level == level, , drop = FALSE]
    for (name in names_vec) {
      part <- path[path$name == name, , drop = FALSE]
      values[name, part$lambda_index] <- part$frequency
    }
  }
  summary <- if (isTRUE(stability$ran)) {
    stability$summary[stability$summary$level == level, , drop = FALSE]
  } else {
    NULL
  }
  rows <- data.frame(
    name = names_vec,
    entry_s_median = profile$entry_s_median,
    stability_max_frequency = if (is.null(summary)) {
      rep(NA_real_, length(names_vec))
    } else {
      summary$max_frequency[match(names_vec, summary$name)]
    },
    stability_selected = if (is.null(summary)) {
      rep(NA, length(names_vec))
    } else {
      summary$selected[match(names_vec, summary$name)]
    },
    stringsAsFactors = FALSE
  )

  # Row order: the earliest typical entry at the top ----------------------
  ordering <- if (identical(source, "stability")) {
    order(-rows$stability_max_frequency, -rowMeans(values, na.rm = TRUE),
          rows$name)
  } else {
    order(-ifelse(is.na(rows$entry_s_median), -Inf, rows$entry_s_median),
          -rowMeans(values, na.rm = TRUE), rows$name)
  }
  rows <- rows[ordering, , drop = FALSE]
  values <- values[rows$name, , drop = FALSE]
  rownames(rows) <- NULL
  max_rows <- as.integer(max_rows)
  if (length(max_rows) != 1L || is.na(max_rows) || max_rows < 1L) {
    stop("`max_rows` must be a positive integer.", call. = FALSE)
  }
  shown <- seq_len(min(nrow(rows), max_rows))

  # The drawing -----------------------------------------------------------
  labels <- rows$name[shown]
  if (identical(source, "stability")) {
    labels <- ifelse(rows$stability_selected[shown] %in% TRUE,
                     paste0(labels, " *"), labels)
  }
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  label_width <- max(graphics::strwidth(labels, units = "inches",
                                        cex = label_cex))
  graphics::par(mai = c(1.0, label_width + 0.3, 0.9, 1.1))
  n_rows <- length(shown)
  graphics::plot.new()
  graphics::plot.window(
    xlim = c(0.5, n_points + 0.5), ylim = c(n_rows + 0.5, 0.5),
    xaxs = "i", yaxs = "i", ...
  )
  graphics::rect(0.5, n_rows + 0.5, n_points + 0.5, 0.5,
                 col = "#D9D9D9", border = NA)
  breaks <- seq(0, 1, length.out = length(.FB_SEL_HEAT_COLOURS) + 1L)
  for (i in seq_len(n_rows)) {
    v <- values[shown[[i]], ]
    ok <- !is.na(v)
    # A row no fold could estimate stays grey.
    if (!any(ok)) {
      next
    }
    bin <- findInterval(v[ok], breaks, rightmost.closed = TRUE,
                        all.inside = TRUE)
    graphics::rect(which(ok) - 0.5, i - 0.5, which(ok) + 0.5, i + 0.5,
                   col = .FB_SEL_HEAT_COLOURS[bin], border = NA)
  }
  selected_s <- arm$selection$selected_s
  if (length(selected_s) == 1L && is.finite(selected_s)) {
    at <- match(selected_s, s_grid)
    if (!is.na(at)) {
      graphics::abline(v = at + 0.5, lty = 2L, col = "black")
    }
  }
  if (identical(source, "folds")) {
    marker <- match(rows$entry_s_median[shown], s_grid)
    ok <- !is.na(marker)
    graphics::points(marker[ok], which(ok), pch = 21L, bg = "white",
                     col = "black", cex = 0.9)
  }
  ticks <- unique(round(seq(1, n_points, length.out = min(n_points, 7L))))
  graphics::axis(1L, at = ticks, labels = signif(s_grid[ticks], 2L),
                 cex.axis = 0.8)
  graphics::axis(2L, at = seq_len(n_rows), labels = labels, las = 1L,
                 cex.axis = label_cex, tick = FALSE)
  graphics::box()
  graphics::title(
    main = if (identical(source, "folds")) {
      paste0("Entrance by ", level, " across the folds")
    } else {
      paste0("Entrance by ", level, " across the stability subsamples")
    },
    xlab = expression("relative penalty " * s == lambda / lambda[max] *
                        " (strong to weak, log scale)")
  )
  graphics::mtext(
    if (identical(source, "folds")) {
      paste0(
        "cell: share of the folds able to estimate it that hold it; ",
        "o: median entry; dashed: right edge of the chosen point"
      )
    } else {
      paste0(
        "cell: share of ", stability$subsamples,
        " subsamples that hold it; *: stable at cutoff ",
        format(stability$cutoff), " left of the dashed line"
      )
    },
    side = 3L, line = 0.3, cex = 0.75
  )
  if (length(shown) < nrow(rows)) {
    graphics::mtext(
      paste0("showing ", length(shown), " of ", nrow(rows), " rows"),
      side = 1L, line = 4, cex = 0.7, adj = 0
    )
  }

  # The colour key, in the right margin ------------------------------------
  usr <- graphics::par("usr")
  key_x <- usr[[2L]] + 0.04 * (usr[[2L]] - usr[[1L]])
  key_w <- 0.04 * (usr[[2L]] - usr[[1L]])
  key_y <- usr[[3L]] + breaks * (usr[[4L]] - usr[[3L]])
  graphics::rect(key_x, key_y[-length(key_y)], key_x + key_w, key_y[-1L],
                 col = .FB_SEL_HEAT_COLOURS, border = NA, xpd = NA)
  graphics::rect(key_x, key_y[[1L]], key_x + key_w, key_y[[length(key_y)]],
                 border = "black", xpd = NA)
  key_at <- c(0, 0.5, 1)
  graphics::text(key_x + key_w, usr[[3L]] + key_at * (usr[[4L]] - usr[[3L]]),
                 labels = format(key_at), pos = 4L, cex = 0.7, xpd = NA)

  invisible(list(
    values = values,
    s = s_grid,
    rows = rows,
    level = level,
    source = source,
    selected_s = selected_s
  ))
}
