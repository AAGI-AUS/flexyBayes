# The marginal covariance of the response of a fitted model,
#
#   V = R + sum_k sigma_k^2 Z_k G_k Z_k' + K,
#
# held block by block: R the residual (sigma^2 / w_i, plus each row's known
# variance), one block per random-term variance, and K the off-diagonal part
# of a known covariance matrix. A block is a slot: `index` is the column of
# Z each row loads on (NA where the block has no effect), `value` the
# loading (NULL for one), `key` the row-level grouping within which the
# block links rows, `columns` the data columns that key is built from, and
# `structure` G where it is not the identity.
#
#   (1 | g), ide(g)        Z the indicator of g, G = I
#   a:b                    Z the indicator of the a:b cell, G = I
#   (x || g)               the intercept block, and Z = indicator(g) * x
#   diag(f):g              one block per level l of f: the indicator of g on
#                          the rows at level l, zero elsewhere
#   ar1(row):ar1(col)      Z the indicator of the lattice node, one trial,
#                          G = R_row(rho_row) (x) R_col(rho_col)
#   at(t):ar1(r):ar1(c)    Z the indicator of the trial's lattice node,
#                          G = I_trial (x) R_row(rho_row) (x) R_col(rho_col)
#   vm(g, K), ped(g, A)    Z the indicator of g, G the known covariance; a
#                          precision carrier is inverted
#
# The residual is independent units, or dsum(~ units | f) with one variance
# per level of f. Anything else refuses with `marginal_cov_not_buildable`.


#' The marginal covariance of the response of a fitted model
#'
#' @description
#' Returns the covariance of the response over the rows of a fit's data,
#' given its fitted variance components,
#' \deqn{V = R + \sum_k \sigma_k^2 Z_k G_k Z_k^\top + K,}{V = R + sum_k
#' sigma_k^2 Z_k G_k Z_k' + K,}
#' as a structured object with one block per variance component.
#' `as.matrix()` assembles it as a dense matrix.
#'
#' @details
#' Each variance enters at its posterior mean, \eqn{E[\sigma^2] =
#' \mathrm{estimate}^2 + \mathrm{std.error}^2}{E[sigma^2] = estimate^2 +
#' std.error^2}, read from the fit's variance-component table,
#' `summary(fit)$varcomp`, whose `estimate` and `std.error` are the posterior
#' mean and standard deviation of each standard deviation. For a structure
#' without a correlation parameter this is the posterior mean of V itself.
#' The two correlations of an AR1 x AR1 field enter at their posterior means
#' (`estimate`). On brms, the rows of a sectioned residual report the
#' posterior median of each standard deviation, so there the residual
#' variances are median^2 + std.error^2.
#'
#' The structures built are random intercepts (`(1 | g)`, `ide(g)`),
#' two-way interactions (`a:b`), uncorrelated random slopes (`(x || g)`),
#' `diag(f):g`, an `ar1(row):ar1(col)` field on one trial or
#' `id(trial):ar1(row):ar1(col)` on each trial, and `vm()` or `ped()` known
#' covariances, over an independent residual or one variance per level of
#' a factor (`dsum(~ units | f)`). Observation weights divide the residual
#' variance (\eqn{\sigma^2 / w_i}{sigma^2 / w_i}), a known variance per row
#' is added to it, and a known covariance matrix is added whole.
#'
#' A model with any other structure is refused with the reason code
#' `marginal_cov_not_buildable`, naming the terms.
#'
#' @param fit A fitted `flexybayes` model, or an `fb_marginal_cov` object to
#'   take a subset of or to evaluate at other variances.
#' @param rows Integer positions of the rows to cover: rows of the fit's
#'   data, or of `fit` when it is an `fb_marginal_cov` object. `NULL` (the
#'   default) covers every row.
#' @param variances `TRUE` (the default) to use the fit's own variance
#'   components; `FALSE` for the structure alone, with every variance `NA`;
#'   or a list with `residual` (the residual variance, or one per section
#'   named by level), `terms` (a variance per block, named by block) and,
#'   for an AR1 x AR1 field, `rho` (a list of `c(row = , col = )` named by
#'   block), all on the variance scale.
#'
#' @returns An object of class `fb_marginal_cov`, a list with
#' \describe{
#'   \item{`n`, `rows`}{The number of rows covered, and their positions in
#'     the fit's data.}
#'   \item{`residual`}{A list: `variance` (one value, or one per section
#'     named by level), `by` (the sectioning column, or `NULL`), `level`
#'     (each row's section, or `NULL`), `weights` (each row's observation
#'     weight, or `NULL`) and `known_variance` (each row's known variance,
#'     or `NULL`).}
#'   \item{`blocks`}{A named list, one entry per variance component:
#'     `name`, `label`, `component` (the names the variance may carry in
#'     `summary(fit)$varcomp`), `variance`, `index`, `n_levels`, `value`,
#'     `key`, `columns` and `structure` as described above, and for an
#'     AR1 x AR1 field `rho` and `kronecker`, the factors `n_at` (trials),
#'     `row` and `col` (AR1 correlation matrices) with
#'     \eqn{G = I_{n_{at}} \otimes R_{row} \otimes R_{col}}{G = I_n_at (x)
#'     row (x) col}.}
#'   \item{`known_covariance`}{The known covariance matrix over the rows,
#'     whose diagonal is `residual$known_variance`, or `NULL`.}
#'   \item{`groups`}{Whether V is block-diagonal over a grouping column:
#'     `applicable`, `column`, `blocks` (row positions per group, the
#'     grouping with the most groups that every block nests within) and
#'     `reason`.}
#'   \item{`varcomp`}{The variance-component table the values were read
#'     from, or `NULL`.}
#' }
#'
#' @seealso [fb_variance_draws()] for posterior draws of the variance
#'   components.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(site = factor(rep(1:6, each = 5)))
#'   d$y <- rnorm(6)[d$site] + rnorm(30)
#'   fit <- flexybayes(y ~ 1, random = ~site, data = d, backend = "inla",
#'                     verbose = FALSE)
#'   v <- fb_marginal_cov(fit, rows = 1:10)
#'   v$blocks$site$variance
#'   as.matrix(v)
#' }
#' }
#' @export
fb_marginal_cov <- function(fit, rows = NULL, variances = TRUE) {
  if (inherits(fit, "fb_marginal_cov")) {
    x <- fit
    if (isFALSE(variances)) {
      x <- .fb_marginal_cov_at(x, NULL)
    }
  } else {
    x <- .fb_marginal_cov_structure(fit)
    if (isTRUE(variances)) {
      x <- .fb_marginal_cov_at(x, .fb_marginal_cov_read(fit, x))
    } else if (!is.list(variances)) {
      x <- .fb_marginal_cov_at(x, NULL)
    }
  }
  if (is.list(variances)) {
    x <- .fb_marginal_cov_at(x, variances)
  } else if (!isTRUE(variances) && !isFALSE(variances)) {
    stop(
      "`variances` must be TRUE, FALSE or a list of variances.",
      call. = FALSE
    )
  }
  if (!is.null(rows)) {
    x <- .fb_marginal_cov_subset(x, rows)
  }
  x
}


#' @rdname fb_marginal_cov
#' @param x An `fb_marginal_cov` object.
#' @param cols Integer positions of the columns, for a cross-covariance
#'   between two sets of rows; defaults to `rows`.
#' @param ... Unused.
#' @export
as.matrix.fb_marginal_cov <- function(x, rows = NULL, cols = rows, ...) {
  rows <- .fb_marginal_cov_positions(x, rows, "rows")
  cols <- .fb_marginal_cov_positions(x, cols, "cols")
  missing <- c(
    if (anyNA(x$residual$variance)) "the residual",
    names(Filter(function(b) is.na(b$variance), x$blocks))
  )
  if (length(missing) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "marginal_cov_not_buildable",
      message = paste0(
        "This marginal covariance carries no variance for ",
        paste(missing, collapse = ", "),
        ". Build it with fb_marginal_cov(variances = TRUE) or supply them as ",
        "a list."
      ),
      terms = missing
    ))
  }

  out <- if (identical(rows, cols)) {
    diag(.fb_marginal_cov_residual(x, rows), nrow = length(rows))
  } else {
    m <- matrix(0, nrow = length(rows), ncol = length(cols))
    hit <- match(rows, cols)
    at <- which(!is.na(hit))
    m[cbind(at, hit[at])] <- .fb_marginal_cov_residual(x, rows[at])
    m
  }
  for (b in x$blocks) {
    out <- out + b$variance * .fb_marginal_cov_loading(b, rows, cols)
  }
  if (!is.null(x$known_covariance)) {
    off <- x$known_covariance[rows, cols, drop = FALSE]
    off[outer(rows, cols, "==")] <- 0
    out <- out + off
  }
  dimnames(out) <- NULL
  out
}


#' @export
print.fb_marginal_cov <- function(x, ...) {
  cat("<fb_marginal_cov> ", x$n, " rows\n", sep = "")
  residual <- if (is.null(x$residual$by)) {
    "independent"
  } else {
    paste0("one variance per level of `", x$residual$by, "`")
  }
  cat("  residual : ", residual, "\n", sep = "")
  for (b in x$blocks) {
    cat("  block    : ", b$label, "  variance ", format(b$variance), "\n",
        sep = "")
  }
  if (!is.null(x$known_covariance)) {
    cat("  known covariance over the rows\n")
  }
  groups <- if (isTRUE(x$groups$applicable) && !is.na(x$groups$column)) {
    paste0(length(x$groups$blocks), " blocks over `", x$groups$column, "`")
  } else if (isTRUE(x$groups$applicable)) {
    "diagonal"
  } else {
    paste0("not block-diagonal (", x$groups$reason, ")")
  }
  cat("  groups   : ", groups, "\n", sep = "")
  invisible(x)
}


# --- the structure ------------------------------------------------------- #

# .fb_marginal_cov_structure() --- the blocks, residual and grouping of a
# fit, every variance NA.
.fb_marginal_cov_structure <- function(fit) {
  fb <- fit$fb %||% fit$extras$fb_terms
  data <- fit$data %||% fit$extras$fb_dataset$data %||% fit$brms$data
  if (is.null(fb) || !is.data.frame(data)) {
    stop(.fb_refusal_condition(
      reason_code = "marginal_cov_not_buildable",
      message = paste0(
        "fb_marginal_cov(): `fit` carries no parsed model or no data, so ",
        "the structure of its covariance cannot be read from it."
      ),
      terms = character(0)
    ))
  }
  info <- fit$extras$call_info %||% list()
  blocks <- .fb_sel_cov_slots(fb, data, info$known_matrices %||% list())
  names(blocks) <- vapply(blocks, function(b) b$name, character(1L))

  residual <- list(
    variance = NA_real_,
    by = NULL,
    level = NULL,
    weights = NULL,
    known_variance = NULL
  )
  for (t in fb$residual_terms %||% list()) {
    if (identical(t$type, "at_units")) {
      section <- .fb_sel_level_index(data[[t$var]])
      residual$by <- t$var
      residual$level <- section$index
      residual$variance <- stats::setNames(
        rep(NA_real_, length(section$levels)),
        section$levels
      )
    }
  }
  weights <- .fb_ir_weights(fb) %||% info$weights
  if (.fb_weights_nonconstant(weights)) {
    residual$weights <- as.numeric(weights)
  }
  if (!is.null(.fb_known_variance_term(fb))) {
    v <- data[[.FB_KNOWN_VARIANCE_COL]]
    if (is.null(v) || length(v) != nrow(data)) {
      .fb_known_variance_refuse(paste0(
        "fb_marginal_cov(): the model carries a known variance per row, and ",
        "its data has no `", .FB_KNOWN_VARIANCE_COL, "` column holding one ",
        "value per row."
      ))
    }
    residual$known_variance <- as.numeric(v)
  }

  columns <- unique(c(
    unlist(lapply(blocks, function(b) b$columns), use.names = FALSE),
    residual$by
  ))
  frame <- data[, columns, drop = FALSE]
  known_covariance <- .fb_known_covariance_values(fb, data)

  structure(
    list(
      n = nrow(data),
      rows = seq_len(nrow(data)),
      residual = residual,
      blocks = blocks,
      known_covariance = known_covariance,
      groups = .fb_marginal_cov_groups(blocks, frame, known_covariance),
      frame = frame,
      varcomp = NULL
    ),
    class = "fb_marginal_cov"
  )
}

# .fb_marginal_cov_groups() --- the grouping V is block-diagonal over.
.fb_marginal_cov_groups <- function(blocks, frame, known_covariance) {
  if (!is.null(known_covariance) &&
        any(known_covariance[upper.tri(known_covariance)] != 0)) {
    return(list(
      applicable = FALSE,
      column = NA_character_,
      blocks = NULL,
      reason = "the known covariance correlates rows across every grouping"
    ))
  }
  g <- .fb_sel_cov_block_structure(unname(blocks), frame)
  list(
    applicable = isTRUE(g$applicable),
    column = g$group %||% NA_character_,
    blocks = g$blocks,
    reason = g$reason
  )
}

# .fb_marginal_cov_subset() --- the object over some of its rows.
.fb_marginal_cov_subset <- function(x, rows) {
  rows <- .fb_marginal_cov_positions(x, rows, "rows")
  pick <- function(v) if (is.null(v)) NULL else v[rows]
  x$blocks <- lapply(x$blocks, function(b) {
    b$index <- b$index[rows]
    b$value <- pick(b$value)
    b$key <- pick(b$key)
    b
  })
  x$residual$level <- pick(x$residual$level)
  x$residual$weights <- pick(x$residual$weights)
  x$residual$known_variance <- pick(x$residual$known_variance)
  if (!is.null(x$known_covariance)) {
    x$known_covariance <- x$known_covariance[rows, rows, drop = FALSE]
  }
  x$frame <- x$frame[rows, , drop = FALSE]
  x$groups <- .fb_marginal_cov_groups(x$blocks, x$frame, x$known_covariance)
  x$rows <- x$rows[rows]
  x$n <- length(rows)
  x
}

# .fb_marginal_cov_positions() --- `rows` checked against the object.
.fb_marginal_cov_positions <- function(x, rows, what) {
  if (is.null(rows)) {
    return(seq_len(x$n))
  }
  if (!is.numeric(rows) || anyNA(rows) || any(rows < 1) || any(rows > x$n) ||
        any(rows != round(rows))) {
    stop(
      "`", what, "` must be integer positions between 1 and ", x$n, ".",
      call. = FALSE
    )
  }
  as.integer(rows)
}


# --- the variances ------------------------------------------------------- #

# .fb_marginal_cov_read() --- the residual and every block's variance at its
# posterior mean, and an AR1 x AR1 field's correlations at theirs.
.fb_marginal_cov_read <- function(fit, x) {
  vc <- summary(fit)$varcomp
  read <- function(names, what) {
    hit <- which(vc$component %in% names)
    if (length(hit) != 1L) {
      stop(.fb_refusal_condition(
        reason_code = "marginal_cov_not_buildable",
        message = paste0(
          "fb_marginal_cov(): the variance-component table carries ",
          length(hit), " rows named ", paste(names, collapse = " or "),
          " for ", what, ", where one was expected. It carries: ",
          paste(vc$component, collapse = ", "), "."
        ),
        terms = what
      ))
    }
    vc[hit, , drop = FALSE]
  }
  second_moment <- function(row) {
    if (!is.finite(row$estimate) || !is.finite(row$std.error)) {
      stop(.fb_refusal_condition(
        reason_code = "marginal_cov_not_buildable",
        message = paste0(
          "fb_marginal_cov(): the fit reports no posterior standard ",
          "deviation for `", row$component, "`, so the posterior mean of its ",
          "variance, estimate^2 + std.error^2, cannot be formed."
        ),
        terms = row$component
      ))
    }
    row$estimate^2 + row$std.error^2
  }

  residual <- if (is.null(x$residual$by)) {
    second_moment(read("sigma", "the residual"))
  } else {
    levels <- names(x$residual$variance)
    stats::setNames(
      vapply(levels, function(l) {
        second_moment(read(paste0("sigma_", l), paste0("the residual at ", l)))
      }, numeric(1L)),
      levels
    )
  }
  terms <- numeric(0)
  rho <- list()
  for (b in x$blocks) {
    terms[[b$name]] <- second_moment(read(b$components, b$label))
    if (identical(b$structure$type, "ar1xar1")) {
      rho[[b$name]] <- c(
        row = read(b$structure$rho[["row"]], b$label)$estimate,
        col = read(b$structure$rho[["col"]], b$label)$estimate
      )
    }
  }
  list(residual = residual, terms = terms, rho = rho, varcomp = vc)
}

# .fb_marginal_cov_at() --- the object at the given variances, or with every
# variance NA when `v` is NULL.
.fb_marginal_cov_at <- function(x, v) {
  if (is.null(v)) {
    x$residual$variance[] <- NA_real_
    x$blocks <- lapply(x$blocks, function(b) {
      b$variance <- NA_real_
      b$rho <- NULL
      b$kronecker <- NULL
      b
    })
    x$varcomp <- NULL
    return(x)
  }
  residual <- v$residual
  if (!is.null(x$residual$by) && length(residual) != 1L) {
    residual <- residual[names(x$residual$variance)]
  }
  if (!is.numeric(residual) || anyNA(residual) || any(residual <= 0) ||
        !(length(residual) == 1L ||
            length(residual) == length(x$residual$variance))) {
    stop(
      "`variances$residual` must be one positive variance, or one per ",
      "level of the residual's sections.",
      call. = FALSE
    )
  }
  x$residual$variance[] <- residual
  for (nm in names(x$blocks)) {
    b <- x$blocks[[nm]]
    variance <- v$terms[nm]
    if (anyNA(variance) || !is.finite(variance) || variance < 0) {
      stop(
        "`variances$terms` must give a non-negative variance for block `",
        nm, "`.",
        call. = FALSE
      )
    }
    b$variance <- unname(variance)
    if (identical(b$structure$type, "ar1xar1")) {
      rho <- v$rho[[nm]]
      if (is.null(rho) || !all(c("row", "col") %in% names(rho)) ||
            any(abs(rho[c("row", "col")]) >= 1)) {
        stop(
          "`variances$rho` must give c(row = , col = ) correlations in ",
          "(-1, 1) for block `", nm, "`.",
          call. = FALSE
        )
      }
      b$rho <- rho[c("row", "col")]
      b$kronecker <- list(
        n_at = as.integer(b$structure$n_at),
        row = .fb_ar1_correlation(b$structure$n_row, b$rho[["row"]]),
        col = .fb_ar1_correlation(b$structure$n_col, b$rho[["col"]])
      )
    }
    x$blocks[[nm]] <- b
  }
  x$varcomp <- v$varcomp
  x
}


# --- assembly ------------------------------------------------------------ #

# .fb_marginal_cov_residual() --- the residual's diagonal over some rows.
.fb_marginal_cov_residual <- function(x, rows) {
  r <- x$residual
  out <- if (is.null(r$by)) {
    r$variance[[1L]]
  } else {
    unname(r$variance[r$level[rows]])
  }
  if (!is.null(r$weights)) {
    out <- out / r$weights[rows]
  }
  if (!is.null(r$known_variance)) {
    out <- out + r$known_variance[rows]
  }
  out
}

# .fb_marginal_cov_loading() --- Z G Z' of one block between two row sets.
.fb_marginal_cov_loading <- function(b, rows, cols) {
  type <- b$structure$type %||% "identity"
  if (identical(type, "ar1xar1")) {
    return(.fb_sel_ar1_field_correlation(b, b$rho, rows, cols))
  }
  if (identical(type, "known")) {
    return(b$structure$k[b$index[rows], b$index[cols], drop = FALSE])
  }
  same <- outer(b$index[rows], b$index[cols], "==")
  same[is.na(same)] <- FALSE
  if (is.null(b$value)) {
    same * 1
  } else {
    same * outer(b$value[rows], b$value[cols])
  }
}

# .fb_ar1_correlation() --- the AR1 correlation matrix over n ordered levels.
.fb_ar1_correlation <- function(n, rho) {
  rho^abs(outer(seq_len(n), seq_len(n), "-"))
}

# .fb_sel_ar1_field_correlation() --- the correlation of a per-trial
# AR1 x AR1 field between two sets of rows: zero across trials.
.fb_sel_ar1_field_correlation <- function(slot, rho, rows_a, rows_b) {
  n_row <- slot$structure$n_row
  n_col <- slot$structure$n_col
  position <- function(rows) {
    i <- slot$index[rows] - 1L
    list(
      trial = i %/% (n_row * n_col),
      row = (i %/% n_col) %% n_row,
      col = i %% n_col
    )
  }
  a <- position(rows_a)
  b <- position(rows_b)
  outer(a$trial, b$trial, "==") *
    rho[["row"]]^abs(outer(a$row, b$row, "-")) *
    rho[["col"]]^abs(outer(a$col, b$col, "-"))
}


# --- the blocks, slot by slot -------------------------------------------- #

# .fb_sel_level_index() --- a grouping column as integer levels, built as the
# emitters build their indices: `factor()` of the column.
.fb_sel_level_index <- function(x) {
  f <- factor(x)
  list(index = as.integer(f), levels = levels(f))
}

# .fb_sel_cov_slot() --- one variance's structure. `components` are the names
# the variance may carry in a fit's variance-component table: the INLA and
# brms spellings differ for a slope and for a known-covariance term.
.fb_sel_cov_slot <- function(name,
                             label,
                             components,
                             index,
                             n_levels,
                             value = NULL,
                             key = NULL,
                             columns = character(0),
                             structure = NULL) {
  list(
    name = name,
    label = label,
    components = components,
    index = index,
    n_levels = n_levels,
    value = value,
    key = key,
    columns = columns,
    structure = structure
  )
}

# .fb_sel_cov_slots() --- the model's random structure as slots, or a refusal
# naming every term, random or residual, that is not built.
.fb_sel_cov_slots <- function(fb, data, known_matrices = list()) {
  residual <- fb$residual_terms %||% list()
  residual_types <- vapply(
    residual,
    function(t) t$type %||% "unknown",
    character(1L)
  )
  refused <- residual_types[!residual_types %in% c("units", "at_units")]
  if (length(refused) > 0L) {
    refused <- paste0("residual ", refused)
  }

  slots <- list()
  for (t in fb$random_terms %||% list()) {
    type <- t$type %||% "unknown"
    built <- switch(
      type,
      simple = ,
      ide = ,
      id = .fb_sel_slots_intercept(t$var, data),
      nested = .fb_sel_slots_nested(t$outer, t$inner, data),
      simple_slope_uncor = .fb_sel_slots_slope(t, data),
      at_simple = if (is.null(t$level)) {
        .fb_sel_slots_diag(t$outer, t$inner, data)
      },
      ar1_spatial = if (isTRUE(t$col_ar1)) {
        .fb_sel_slots_ar1(t$row_var, t$col_var, data, t$at_var)
      },
      vm = ,
      ped = .fb_sel_slots_known(t, data, known_matrices),
      NULL
    )
    if (is.null(built)) {
      named <- unlist(t[c("var", "outer", "inner", "row_var", "col_var",
                          "vars", "at_var", "level")])
      refused <- c(
        refused,
        paste0(type, "(", paste(named, collapse = ", "), ")")
      )
    } else {
      slots <- c(slots, built)
    }
  }

  if (length(refused) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "marginal_cov_not_buildable",
      message = paste0(
        "fb_marginal_cov(): the marginal covariance is built from random ",
        "intercepts, two-way interactions, uncorrelated random slopes, ",
        "diag(f):g, an ar1(row):ar1(col) field on one trial or ",
        "id(trial):ar1(row):ar1(col) on each trial, and vm() or ped() known ",
        "covariances, over an independent residual or one variance per ",
        "section (dsum(~ units | f)), and this model also carries: ",
        paste(unique(refused), collapse = ", "),
        "."
      ),
      terms = unique(refused)
    ))
  }
  slots
}

.fb_sel_slots_intercept <- function(var, data) {
  levels <- .fb_sel_level_index(data[[var]])
  list(.fb_sel_cov_slot(
    name = var,
    label = var,
    components = paste0("sd_", var),
    index = levels$index,
    n_levels = length(levels$levels),
    key = as.character(data[[var]]),
    columns = var
  ))
}

.fb_sel_slots_nested <- function(outer, inner, data) {
  key <- paste0(
    as.character(data[[outer]]), ":", as.character(data[[inner]])
  )
  name <- paste0(outer, ":", inner)
  levels <- .fb_sel_level_index(key)
  list(.fb_sel_cov_slot(
    name = name,
    label = name,
    components = paste0("sd_", name),
    index = levels$index,
    n_levels = length(levels$levels),
    key = key,
    columns = c(outer, inner)
  ))
}

.fb_sel_slots_slope <- function(t, data) {
  group <- t$var
  covariate <- t$slope_var
  out <- if (isTRUE(t$with_intercept)) {
    .fb_sel_slots_intercept(group, data)
  } else {
    list()
  }
  levels <- .fb_sel_level_index(data[[group]])
  c(out, list(.fb_sel_cov_slot(
    name = paste0(covariate, "_", group),
    label = paste0("(", covariate, " || ", group, ")"),
    components = c(
      paste0("sd_", covariate, "_", group),
      paste0("sd_", group, "__", covariate)
    ),
    index = levels$index,
    n_levels = length(levels$levels),
    value = as.numeric(data[[covariate]]),
    key = as.character(data[[group]]),
    columns = group
  )))
}

.fb_sel_slots_diag <- function(outer, inner, data) {
  at <- as.character(data[[outer]])
  levels <- .fb_sel_level_index(data[[inner]])
  lapply(levels(factor(data[[outer]])), function(l) {
    index <- levels$index
    index[is.na(at) | at != l] <- NA_integer_
    name <- paste0(inner, "__", outer, l)
    .fb_sel_cov_slot(
      name = name,
      label = paste0("diag(", outer, "):", inner, " at ", l),
      components = paste0("sd_", name),
      index = index,
      n_levels = length(levels$levels),
      key = paste0(at, ":", as.character(data[[inner]])),
      columns = c(outer, inner)
    )
  })
}

.fb_sel_slots_ar1 <- function(row_var, col_var, data, at_var = NULL) {
  row <- .fb_sel_level_index(data[[row_var]])
  col <- .fb_sel_level_index(data[[col_var]])
  n_row <- length(row$levels)
  n_col <- length(col$levels)
  label <- paste0("ar1(", row_var, "):ar1(", col_var, ")")
  trial <- 1L
  n_at <- 1L
  key <- NULL
  columns <- c(row_var, col_var)
  if (!is.null(at_var)) {
    at <- .fb_sel_level_index(data[[at_var]])
    trial <- at$index
    n_at <- length(at$levels)
    key <- as.character(data[[at_var]])
    columns <- at_var
    label <- paste0("id(", at_var, "):", label)
  }
  list(.fb_sel_cov_slot(
    name = "spatial",
    label = label,
    components = "sd_spatial",
    # Trial, then row, then column, as in kronecker(I, R_row, R_col).
    index = ((trial - 1L) * n_row + row$index - 1L) * n_col + col$index,
    n_levels = n_at * n_row * n_col,
    key = key,
    columns = columns,
    structure = list(
      type = "ar1xar1",
      n_row = n_row,
      n_col = n_col,
      n_at = n_at,
      rho = c(row = paste0("rho_", row_var), col = paste0("rho_", col_var))
    )
  ))
}

# The known covariance of a vm() or ped() term, as the brms emitter densifies
# it (`.fb_brms_data2()`), aligned to the levels of the grouping factor: by
# name where the matrix carries dimnames, by position where it does not, as
# INLA's generic0 index requires.
.fb_sel_slots_known <- function(t, data, known_matrices) {
  var <- t$var
  cov <- t$cov_representation
  format <- cov$format %||% "dense"
  symbol <- cov$data %||% t$mat
  levels <- .fb_sel_level_index(data[[var]])
  m <- if (is.character(symbol) && length(symbol) == 1L) {
    known_matrices[[symbol]]
  }
  if (is.null(m) ||
        !format %in% c("dense", "chol", "precision",
                       "pedigree_sparse_precision")) {
    return(NULL)
  }
  k <- switch(
    format,
    dense = as.matrix(m),
    chol = tcrossprod(as.matrix(m)),
    precision = solve(as.matrix(m)),
    pedigree_sparse_precision = as.matrix(solve(m))
  )
  names_k <- rownames(k)
  if (!is.null(names_k) && all(levels$levels %in% names_k)) {
    k <- k[levels$levels, levels$levels, drop = FALSE]
  } else if (nrow(k) != length(levels$levels)) {
    return(NULL)
  }
  dimnames(k) <- NULL
  list(.fb_sel_cov_slot(
    name = var,
    label = paste0(t$type, "(", var, ")"),
    components = c(paste0("sd_", var, "_id"), paste0("sd_", var)),
    index = levels$index,
    n_levels = length(levels$levels),
    columns = var,
    structure = list(type = "known", k = k)
  ))
}


# --- the blocks V is diagonal over ----------------------------------------- #

# .fb_sel_cov_block_structure() --- the grouping column every slot nests
# within, with the most blocks among those that qualify, or the reason there
# is none. A slot with a structure matrix and no key links rows across every
# level of its index, so it admits no grouping.
.fb_sel_cov_block_structure <- function(slots, data) {
  n <- nrow(data)
  if (length(slots) == 0L) {
    return(list(
      applicable = TRUE,
      group = NA_character_,
      blocks = as.list(seq_len(n)),
      reason = "no random terms, so V is diagonal"
    ))
  }
  linking <- Filter(function(s) !is.null(s$structure) && is.null(s$key), slots)
  if (length(linking) > 0L) {
    return(list(
      applicable = FALSE,
      reason = paste0(
        paste(vapply(linking, function(s) s$label, character(1L)),
              collapse = ", "),
        " correlates rows across every level it indexes"
      )
    ))
  }
  best <- NULL
  for (g in unique(unlist(lapply(slots, function(s) s$columns)))) {
    gv <- as.character(data[[g]])
    nests <- all(vapply(slots, function(s) {
      keep <- !is.na(s$index)
      spans <- tapply(gv[keep], s$key[keep], function(x) {
        length(unique(x))
      })
      all(spans <= 1L)
    }, logical(1L)))
    if (!nests) {
      next
    }
    blocks <- unname(split(seq_len(n), factor(gv)))
    if (is.null(best) || length(blocks) > length(best$blocks)) {
      best <- list(group = g, blocks = blocks)
    }
  }
  if (is.null(best)) {
    return(list(
      applicable = FALSE,
      reason = "no grouping column holds every random term"
    ))
  }
  list(
    applicable = TRUE,
    group = best$group,
    blocks = best$blocks,
    reason = paste0("every random term nests within `", best$group, "`")
  )
}
