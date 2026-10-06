# fixed_design.R: the design matrix of a fixed formula, its columns' terms,
# and the order in which ASReml 4.2 enters them.

# The rank tolerance for an aliased column, relative to each centred column's
# norm: the qr() and lm() default.
.FB_FIXED_RANK_TOL <- 1e-7


#' The design matrix of a fixed formula
#'
#' Builds the fixed-effect design matrix of `formula` over `data`, as
#' flexyBayes builds it for its aliasing check, with the term of each column
#' and the order in which ASReml 4.2 enters the columns.
#'
#' Factor levels and contrasts are those of all of `data`. With `train_rows`,
#' a transform that learns from the rows it is evaluated on (`scale()`,
#' `poly()`, `splines::ns()`, `splines::bs()`) learns from the training rows
#' alone and is applied with those parameters to every row, the rule of
#' [stats::predict()]. Missing values are kept as rows of `NA`.
#'
#' @param formula A one-sided formula of fixed terms.
#' @param data A data frame.
#' @param train_rows Row numbers of the training rows, or `NULL` (default)
#'   to learn transforms from every row.
#' @returns A list: `X` (the matrix), `terms`, `frame` (the model frame),
#'   `xlev` (the factor levels), `assign` (the term index of each column),
#'   `column_terms` (the term label of each column) and `asreml_order` (the
#'   column positions in ASReml's order: term by term, and within a term by
#'   descending factor levels with the first factor varying slowest).
#' @rdname fb_extension
#' @export
fb_fixed_design <- function(formula, data, train_rows = NULL) {
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  tt <- stats::terms(mf)
  if (!is.null(train_rows)) {
    train_terms <- stats::terms(stats::model.frame(
      formula,
      data = data[train_rows, , drop = FALSE],
      na.action = stats::na.pass
    ))
    if (!identical(attr(train_terms, "predvars"), attr(tt, "predvars"))) {
      tt <- train_terms
      mf <- stats::model.frame(tt, data = data, na.action = stats::na.pass)
    }
  }
  X <- stats::model.matrix(tt, mf)
  design <- list(
    X = X,
    terms = tt,
    frame = mf,
    xlev = stats::.getXlevels(tt, mf),
    assign = attr(X, "assign")
  )
  design$column_terms <- .fb_fixed_column_terms(design)
  design$asreml_order <- .fb_fixed_asreml_order(design)
  design
}


# .fb_fixed_asreml_order() --- the columns in the order ASReml 4.2 enters
# them. A term whose column names do not parse into factor levels keeps the
# design's order.
.fb_fixed_asreml_order <- function(design) {
  X <- design$X
  labels <- attr(design$terms, "term.labels")
  out <- integer(0)
  for (a in sort(unique(design$assign))) {
    columns <- which(design$assign == a)
    if (a == 0L || length(columns) < 2L) {
      out <- c(out, columns)
      next
    }
    vars <- strsplit(labels[[a]], ":", fixed = TRUE)[[1L]]
    parts <- strsplit(colnames(X)[columns], ":", fixed = TRUE)
    key <- vapply(seq_along(vars), function(k) {
      levels <- design$xlev[[vars[[k]]]]
      vapply(parts, function(p) {
        if (length(p) != length(vars)) {
          return(NA_integer_)
        }
        if (is.null(levels)) {
          return(if (identical(p[[k]], vars[[k]])) 0L else NA_integer_)
        }
        match(p[[k]], paste0(vars[[k]], levels))
      }, integer(1L))
    }, integer(length(columns)))
    key <- matrix(key, nrow = length(columns))
    if (anyNA(key)) {
      out <- c(out, columns)
      next
    }
    out <- c(out, columns[do.call(order, lapply(seq_len(ncol(key)),
                                              function(k) -key[, k]))])
  }
  out
}


# .fb_fixed_column_terms() --- the term label of each column.
.fb_fixed_column_terms <- function(design) {
  labels <- c("(Intercept)", attr(design$terms, "term.labels"))
  labels[design$assign + 1L]
}
