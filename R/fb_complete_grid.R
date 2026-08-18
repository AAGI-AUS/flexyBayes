# fb_complete_grid.R -- reinstating design cells that are absent from the
# data.
#
# One implementation, two entry points. `fb_complete_grid()` is the
# user-facing helper for padding a field book before a fit;
# `.fb_complete_design_grid()` (R/na_action.R) is what
# `na_action = "augment"` calls on the way into dispatch. Both run
# `.fb_grid_complete()`, so the helper and the automatic path cannot
# drift into completing the same trial two different ways.
#
# The rule both enforce is the same one: a design cell absent from the
# data can be reinstated when every model column at that cell is
# determined -- it is one of the index variables, or it takes a single
# value across the trial. A column that varies is a real quantity nobody
# recorded there, and writing one in fabricates an observation. ASReml's
# nin89 example codes such plots as LANCER, a variety name standing in
# for "nothing was sown here"; this package makes that an explicit,
# named choice (`unused_level =`) rather than a default.


# ---------------------------------------------------------------- #
# The shared core                                                   #
# ---------------------------------------------------------------- #

# .fb_grid_varies() --- which candidate columns take more than one value.
#
# @noRd
# @keywords internal
.fb_grid_varies <- function(data, candidate_vars) {
  if (length(candidate_vars) == 0L) {
    return(character(0))
  }
  candidate_vars[vapply(
    candidate_vars,
    function(v) {
      x <- data[[v]]
      length(unique(x[!is.na(x)])) > 1L
    },
    logical(1L)
  )]
}

# .fb_grid_fillable() --- can a column take a named unused level?
#
# The hatch is for design factors. A measured covariate has no level to
# name, so a varying numeric column refuses whatever `unused_level` says.
#
# @noRd
# @keywords internal
.fb_grid_fillable <- function(data, vars) {
  vars[vapply(
    vars,
    function(v) is.factor(data[[v]]) || is.character(data[[v]]),
    logical(1L)
  )]
}

# .fb_grid_absent_rows() --- the index-variable combinations with no row.
#
# @noRd
# @keywords internal
.fb_grid_absent_rows <- function(data, index_vars) {
  levs <- lapply(index_vars, function(v) {
    x <- data[[v]]
    if (is.factor(x)) levels(x) else sort(unique(x))
  })
  names(levs) <- index_vars
  full <- expand.grid(levs, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

  key <- function(df) {
    do.call(paste, c(
      lapply(index_vars, function(v) as.character(df[[v]])),
      list(sep = "\r")
    ))
  }
  list(full = full, absent = which(!key(full) %in% key(data)))
}

# .fb_grid_complete() --- the one completion both entry points run.
#
# Returns `list(data = <possibly extended>, added = <int>, filled =
# <character>)`. `filled` names the columns an explicit `unused_level`
# wrote into, and is empty on the automatic path, which never fills.
#
# @noRd
# @keywords internal
.fb_grid_complete <- function(
  data,
  index_vars,
  response,
  candidate_vars,
  unused_level = NULL,
  remedy
) {
  unchanged <- list(data = data, added = 0L, filled = character(0))
  if (length(index_vars) == 0L || !all(index_vars %in% names(data))) {
    return(unchanged)
  }

  grid <- .fb_grid_absent_rows(data, index_vars)
  absent <- grid$absent
  if (length(absent) == 0L) {
    return(unchanged)
  }

  # ---- What would have to be invented? ---------------------------------
  candidate_vars <- setdiff(
    intersect(candidate_vars, names(data)),
    c(index_vars, response)
  )
  undetermined <- .fb_grid_varies(data, candidate_vars)
  fillable <- .fb_grid_fillable(data, undetermined)
  unfillable <- setdiff(undetermined, fillable)

  if (length(unfillable) > 0L || (length(fillable) > 0L &&
    is.null(unused_level))) {
    stop(.fb_grid_refusal(
      n_absent = length(absent),
      index_vars = index_vars,
      undetermined = undetermined,
      unfillable = if (is.null(unused_level)) character(0) else unfillable,
      remedy = remedy
    ))
  }

  # ---- Extend the level vocabulary before the rows are bound -----------
  # A factor gains the unused level on the whole column, not only on the
  # rows being added: rbind() of two factors with different level sets
  # drops to character and silently reorders the model's contrasts.
  filled <- character(0)
  if (length(fillable) > 0L) {
    for (v in fillable) {
      if (is.factor(data[[v]])) {
        data[[v]] <- factor(
          as.character(data[[v]]),
          levels = c(levels(data[[v]]), as.character(unused_level))
        )
      }
    }
    filled <- fillable
    warning(
      "flexyBayes: filled ", length(absent), " reinstated design cell(s) ",
      "with the unused level \"", as.character(unused_level), "\" in ",
      paste(filled, collapse = ", "),
      ". Those cells now carry a design assignment nobody recorded; every ",
      "estimate that reads the named column(s) is conditional on it.",
      call. = FALSE
    )
  }

  # ---- Build the reinstated rows ---------------------------------------
  add <- data[rep(1L, length(absent)), , drop = FALSE]
  rownames(add) <- NULL
  for (v in index_vars) {
    val <- grid$full[[v]][absent]
    add[[v]] <- if (is.factor(data[[v]])) {
      factor(as.character(val), levels = levels(data[[v]]))
    } else {
      val
    }
  }
  for (v in filled) {
    add[[v]] <- if (is.factor(data[[v]])) {
      factor(as.character(unused_level), levels = levels(data[[v]]))
    } else {
      as.character(unused_level)
    }
  }
  if (!is.null(response) && response %in% names(data)) {
    add[[response]] <- NA
  }

  out <- rbind(data, add)
  rownames(out) <- NULL
  list(data = out, added = length(absent), filled = filled)
}

# .fb_grid_refusal() --- the shared condition, with the caller's remedy.
#
# @noRd
# @keywords internal
.fb_grid_refusal <- function(
  n_absent,
  index_vars,
  undetermined,
  unfillable,
  remedy
) {
  head <- paste0(
    n_absent, " cell(s) of the ", paste(index_vars, collapse = " x "),
    " grid are absent from the data altogether, and ",
    paste(undetermined, collapse = ", "),
    " is not constant across the trial -- so completing the grid would ",
    "mean inventing a level of it for a plot nobody recorded. ASReml's ",
    "own nin89 example codes such plots as LANCER, a variety name ",
    "standing in for \"nothing was sown here\"; flexyBayes will not write ",
    "a level you did not. "
  )
  numeric_note <- if (length(unfillable) > 0L) {
    paste0(
      "An unused level cannot stand in for ",
      paste(unfillable, collapse = ", "),
      ": a measured covariate has no level to name. "
    )
  } else {
    ""
  }
  .fb_refusal_condition(
    reason_code = "augment_cell_not_determinable",
    message = paste0(head, numeric_note, remedy),
    family_class = "flexybayes_augment_cell_not_determinable"
  )
}


# ---------------------------------------------------------------- #
# The exported helper                                               #
# ---------------------------------------------------------------- #

#' Complete a design grid before fitting
#'
#' Reinstates the design cells that are absent from a data frame,
#' returning the Cartesian product of the index variables with the
#' response set to `NA` on every cell that was added. A field book
#' already records a lost plot this way -- one row per sown plot, with
#' the design columns filled in and the yield missing -- so this helper
#' is for the trials where the plot is missing from the file altogether.
#'
#' The completed frame is what `na_action = "augment"` (the default) then
#' carries through the fit, keeping the index set a structured covariance
#' is built over intact. A separable `ar1(row):ar1(col)` field needs one
#' observation per node of the array, which is the same requirement
#' ASReml states for `residual = ~ ar1:ar1`.
#'
#' A cell can only be reinstated when every other model column at it is
#' determined: it is one of the index variables, or it takes one value
#' across the whole trial. A column that varies -- a variety label, a
#' block -- is a real quantity nobody recorded at the absent cell, and
#' the default is to refuse rather than write one in. `unused_level`
#' opens that door explicitly, filling the varying design factors with a
#' level of your naming and warning with every column it wrote to. It is
#' ASReml's nin89 LANCER coding, made a decision with a name on it rather
#' than a default. A varying *numeric* column is refused whatever
#' `unused_level` says: a measured covariate has no level to name.
#'
#' @param data A data frame holding one row per recorded design cell.
#' @param index A one-sided formula naming the design index variables,
#'   for example `~ row * col`. The operators are ignored; what is read
#'   is the set of variables, whose observed levels are crossed.
#' @param response Name of the response column, as a single string. The
#'   reinstated cells carry `NA` there.
#' @param unused_level A single string, or `NULL` (the default). When
#'   supplied, design factors that vary across the trial are filled with
#'   this level on the reinstated cells and a warning names every column
#'   so filled. `NULL` refuses instead.
#' @returns A data frame: `data` with one row appended per absent
#'   combination of the index variables, the response `NA` on each, and
#'   the original row order preserved ahead of them. Returned unchanged
#'   when the grid is already complete.
#' @seealso [flexybayes()], whose `na_action` argument carries the
#'   completed grid into the fit.
#' @examples
#' trial <- expand.grid(row = factor(1:4), col = factor(1:3))
#' trial$yield <- rnorm(nrow(trial))
#' holed <- trial[-c(2L, 7L), ]
#' padded <- fb_complete_grid(holed, ~ row * col, response = "yield")
#' nrow(padded)
#' sum(is.na(padded$yield))
#' @export
fb_complete_grid <- function(
  data,
  index = ~ row * col,
  response,
  unused_level = NULL
) {
  .check_complete_grid_args(data, index, response, unused_level)

  index_vars <- all.vars(index)
  missing_vars <- setdiff(c(index_vars, response), names(data))
  if (length(missing_vars) > 0L) {
    stop(
      "fb_complete_grid(): `", paste(missing_vars, collapse = "`, `"),
      "` is not a column of `data`.",
      call. = FALSE
    )
  }

  out <- .fb_grid_complete(
    data = data,
    index_vars = index_vars,
    response = response,
    candidate_vars = names(data),
    unused_level = unused_level,
    remedy = paste0(
      "Supply the field-book rows instead -- one row per sown plot, with ",
      "the design columns filled in and the response set to NA -- or name ",
      "the level the empty cells should carry, ",
      "fb_complete_grid(..., unused_level = \"UNSOWN\")."
    )
  )
  out$data
}

# .check_complete_grid_args() --- argument contract for the helper.
#
# @noRd
# @keywords internal
.check_complete_grid_args <- function(data, index, response, unused_level) {
  if (!is.data.frame(data)) {
    stop(
      "fb_complete_grid(): `data` must be a data frame, not an object of ",
      "class <", paste(class(data), collapse = ", "), ">.",
      call. = FALSE
    )
  }
  if (!inherits(index, "formula") || length(index) != 2L) {
    stop(
      "fb_complete_grid(): `index` must be a one-sided formula naming the ",
      "design index variables, for example ~ row * col.",
      call. = FALSE
    )
  }
  if (length(all.vars(index)) == 0L) {
    stop(
      "fb_complete_grid(): `index` names no variables, so there is no grid ",
      "to complete.",
      call. = FALSE
    )
  }
  if (missing(response) || !is.character(response) || length(response) != 1L) {
    stop(
      "fb_complete_grid(): `response` must be the name of one column of ",
      "`data`, as a single string.",
      call. = FALSE
    )
  }
  if (!is.null(unused_level)) {
    if (length(unused_level) != 1L || is.na(unused_level)) {
      stop(
        "fb_complete_grid(): `unused_level` must be a single non-missing ",
        "value naming the level the reinstated cells should carry, or NULL ",
        "to refuse instead.",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}
