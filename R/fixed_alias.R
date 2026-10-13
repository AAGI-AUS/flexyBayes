# Aliased fixed-effect columns in flexybayes().
#
# A fixed-effect column that is a linear combination of the columns before it
# on the rows with a response has no estimate: the data say nothing about
# it. INLA and brms both fit such a model without a message and report the
# column at its prior, and a prediction that depends on it is the prior's.
# flexybayes() therefore refuses such a model by default
# (`fixed_effects_aliased`). With `fixed_alias = "drop"` it drops the
# columns by the rule measured against asreml 4.2 (R/fixed_design.R): terms
# in model order,
# within a term by descending factor levels with the term's first factor
# varying slowest, and a column that is a combination of the columns before
# it is dropped. A dropped column is fixed at zero, which is how asreml
# reports it, and a prediction whose row is not a combination of the fitting
# rows is `NA` with status "Aliased".

# The prior precision that fixes a dropped column at zero on INLA, on the
# coordinates INLA fits in: a prior standard deviation of 1e-5 there.
.FB_FIXED_ALIAS_PIN_PREC <- 1e10

# The fixed-term types whose design `stats::model.matrix()` builds as the
# engines do. A model with any other fixed term is not checked.
.FB_FIXED_ALIAS_TYPES <- c(
  "factor",
  "continuous",
  "interaction",
  "factor_interaction",
  "factor_numeric_interaction",
  "expression"
)


# .fb_fixed_alias(): the fixed-effect columns aliased on the rows with a
# response, in the order and by the rule asreml 4.2 drops them.
#
# Returns NULL when the check does not apply (no fixed term beside the
# intercept, or a term type outside `.FB_FIXED_ALIAS_TYPES`). Otherwise
# `list(columns, terms, n_columns, n_rows, nbasis)`, where `nbasis` has one
# row per design column and one column per dropped column and spans the
# directions the fitting rows do not estimate.
.fb_fixed_alias <- function(fb, data) {
  terms_in <- fb$fixed_terms %||% list()
  types <- vapply(terms_in, function(t) t$type %||% "", character(1L))
  if (length(terms_in) == 0L || !all(types %in% .FB_FIXED_ALIAS_TYPES)) {
    return(NULL)
  }
  labels <- vapply(
    terms_in,
    function(t) {
      t$label %||% t$var %||% paste(t$vars %||% character(0), collapse = ":")
    },
    character(1L)
  )
  rhs <- stats::as.formula(paste0(
    "~ ",
    if (isTRUE(fb$intercept)) "1 + " else "0 + ",
    paste(labels, collapse = " + ")
  ))
  design <- fb_fixed_design(rhs, data)
  X <- design$X
  response_vars <- intersect(all.vars(str2lang(fb$response)), names(data))
  rows <- which(
    stats::complete.cases(data[response_vars]) & stats::complete.cases(X)
  )
  x_fit <- X[rows, , drop = FALSE]

  m <- x_fit
  if ("(Intercept)" %in% colnames(m)) {
    centre <- colnames(m) != "(Intercept)"
    m[, centre] <- sweep(
      m[, centre, drop = FALSE], 2L, colMeans(m[, centre, drop = FALSE])
    )
  }
  entry <- design$asreml_order
  decomposition <- qr(m[, entry, drop = FALSE], tol = .FB_FIXED_RANK_TOL)
  rank <- decomposition$rank
  dropped <- if (rank < ncol(m)) {
    sort(entry[decomposition$pivot[-seq_len(rank)]])
  } else {
    integer(0)
  }
  kept <- setdiff(seq_len(ncol(X)), dropped)

  nbasis <- matrix(
    0,
    nrow = ncol(X),
    ncol = length(dropped),
    dimnames = list(colnames(X), colnames(X)[dropped])
  )
  if (length(dropped) > 0L) {
    # Each dropped column is, on the fitting rows, the kept columns times `a`.
    a <- qr.coef(
      qr(x_fit[, kept, drop = FALSE], tol = .FB_FIXED_RANK_TOL),
      x_fit[, dropped, drop = FALSE]
    )
    nbasis[kept, ] <- -a
    nbasis[cbind(dropped, seq_along(dropped))] <- 1
  }
  list(
    columns = colnames(X)[dropped],
    terms = design$column_terms[dropped],
    n_columns = ncol(X),
    n_rows = length(rows),
    nbasis = nbasis
  )
}


# .fb_fixed_alias_resolve(): refuse an aliased fixed part, or record it
# on the IR for the emit to drop.
.fb_fixed_alias_resolve <- function(fb, data, fixed_alias, backend, aggregate) {
  alias <- .fb_fixed_alias(fb, data)
  if (is.null(alias) || length(alias$columns) == 0L) {
    return(fb)
  }
  n <- length(alias$columns)
  shown <- utils::head(alias$columns, 12L)
  listed <- paste0(
    paste(shown, collapse = ", "),
    if (n > length(shown)) paste0(" and ", n - length(shown), " more") else ""
  )
  if (identical(fixed_alias, "refuse")) {
    stop(.fb_refusal_condition(
      reason_code = "fixed_effects_aliased",
      message = paste0(
        "flexybayes(): ", n, " of the ", alias$n_columns, " fixed-effect ",
        "columns are linear combinations of the columns before them on the ",
        alias$n_rows, " rows with a response, so the data give no estimate ",
        "for them: ", listed, ". Fitted as written, they would be reported ",
        "at their prior. Pass fixed_alias = \"drop\" to drop them as asreml ",
        "does, with predictions that depend on them returned as NA, or ",
        "remove the terms that create them."
      ),
      columns = alias$columns
    ))
  }
  if (identical(backend, "brms")) {
    stop(
      "flexybayes(): fixed_alias = \"drop\" is fitted on INLA, which fixes ",
      "the ", n, " aliased columns (", listed, ") at zero. Use ",
      "backend = \"inla\" or \"auto\".",
      call. = FALSE
    )
  }
  if (isTRUE(aggregate)) {
    stop(
      "flexybayes(): fixed_alias = \"drop\" is fitted on the per-row INLA ",
      "route, and aggregate = TRUE asks for the aggregated one. Use ",
      "aggregate = \"auto\" or FALSE.",
      call. = FALSE
    )
  }
  fb$fixed_alias <- alias
  fb
}


# .fb_inla_pin_aliased(): control.fixed with each dropped column fixed at
# zero, on the coordinates INLA fits in. `names` are INLA's column names.
.fb_inla_pin_aliased <- function(control_fixed, names) {
  if (length(names) == 0L) {
    return(control_fixed)
  }
  pin <- function(slot, value) {
    entries <- stats::setNames(as.list(rep(value, length(names))), names)
    if (is.null(slot)) {
      return(entries)
    }
    if (!is.list(slot)) {
      return(c(entries, list(default = slot)))
    }
    slot[names] <- entries
    slot
  }
  control_fixed$mean <- pin(control_fixed$mean, 0)
  control_fixed$prec <- pin(control_fixed$prec, .FB_FIXED_ALIAS_PIN_PREC)
  control_fixed
}


# .fb_inla_aliased_names(): INLA's names for the dropped columns: the
# same positions in the design built on the data INLA is handed, whose
# factor levels may have been made syntactic.
.fb_inla_aliased_names <- function(fb, data) {
  alias <- fb$fixed_alias
  if (is.null(alias) || length(alias$columns) == 0L) {
    return(character(0))
  }
  inla_names <- rownames(alias$nbasis)
  legal <- .fb_inla_fixed_coef_names(fb, data)
  user <- setdiff(inla_names, "(Intercept)")
  if (length(legal) == length(user)) {
    inla_names[match(user, inla_names)] <- legal
  }
  inla_names[match(alias$columns, rownames(alias$nbasis))]
}


# .fb_fixed_alias_attach(): the dropped columns recorded on an INLA fit,
# and reported as asreml reports them: estimate zero, no standard error.
.fb_fixed_alias_attach <- function(fit, fb) {
  alias <- fb$fixed_alias
  if (is.null(alias) || !is.list(fit) || is.null(fit$extras)) {
    return(fit)
  }
  fit$extras$fixed_alias <- list(
    mode = "drop",
    rule = paste0(
      "terms in model order; within a term, columns by descending factor ",
      "levels, the first factor varying slowest; a column that is a linear ",
      "combination of the columns before it on the rows with a response is ",
      "fixed at zero (asreml 4.2)"
    ),
    dropped = data.frame(
      column = alias$columns,
      term = alias$terms,
      stringsAsFactors = FALSE
    ),
    n_rows = alias$n_rows,
    nbasis = alias$nbasis
  )
  if (!is.null(fit$extras$call_info)) {
    fit$extras$call_info$fixed_alias <- "drop"
  }
  sf <- fit$inla$summary.fixed
  if (!is.null(sf)) {
    at <- match(
      alias$columns,
      .inla_restore_term_labels(fit$level_labels, rownames(sf))
    )
    sf[at, intersect(c("mean", "0.5quant", "mode"), colnames(sf))] <- 0
    unknown <- c("sd", "0.025quant", "0.975quant", "kld")
    sf[at, intersect(unknown, colnames(sf))] <- NA_real_
    fit$inla$summary.fixed <- sf
  }
  fit
}


# .fb_fixed_alias_status(): "Estimable" or "Aliased" for each row of a
# design, or NULL on a fit that dropped nothing. A row is estimable exactly
# when it is a linear combination of the fitting rows.
.fb_fixed_alias_status <- function(object, X) {
  nbasis <- object$extras$fixed_alias$nbasis
  if (is.null(nbasis) || ncol(nbasis) == 0L) {
    return(NULL)
  }
  nbasis <- nbasis[colnames(X), , drop = FALSE]
  residual <- X %*% nbasis
  scale <- 1 + abs(X) %*% abs(nbasis)
  aliased <- rowSums(abs(residual) > .FB_FIXED_RANK_TOL * scale) > 0L
  unname(ifelse(aliased, "Aliased", "Estimable"))
}
