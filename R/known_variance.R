# known_variance.R -- a known variance per row, added to the residual.
#
# `flexybayes(known_variance = v)` fits
#
#   y_i = x_i' beta + z_i' u + e_i,   e_i ~ N(0, sigma^2 + v_i),
#
# with v_i known and sigma^2 estimated. It is the second stage of a two-stage
# analysis, where y_i is a first-stage estimate (a variety mean in one trial)
# and v_i is the variance the first stage attached to it. `weights = w` is a
# different model, Var(e_i) = sigma^2 / w_i, which scales the estimated
# variance instead of adding to it, so the two do not share an argument.
#
# The lowering, one per engine.
#
#   brms   y | se(s, sigma = TRUE) with s_i = sqrt(v_i). brms writes the
#          likelihood as normal_lpdf(Y | mu, sqrt(square(sigma) + se2)), and
#          with a sigma sub-formula (a heterogeneous residual, the weights
#          offset) sigma is that row's value, so the two add as they do on
#          INLA.
#   INLA   a row effect f(i, w, model = "iid") whose precision is fixed at 1,
#          with w_i = sqrt(v_i) / s on the standardised response (s is its
#          scale). Row i's variance on the original scale gains
#          s^2 w_i^2 = v_i. The Gaussian likelihood's own precision offset
#          cannot carry this: it is one value for every row, and INLA fixes
#          it with no prior.
#
# The INLA row effect sits in the linear predictor, so INLA's fitted values
# include its posterior mean, where brms's exclude the known variance, which
# is part of the likelihood there. fitted() on an INLA fit subtracts it
# (.fb_inla_known_variance_effect()), and the effect is left out of every
# random-effect table.
#
# Where the values live. The vector is checked against the rows it was given
# with, then carried as a column of the analysis data, so it follows every
# row operation the missing-data layer performs: a row dropped under `omit`
# drops its variance, and a cell the design grid adds, which has no response
# and so no likelihood, gets zero. The IR records the column and the vector
# as supplied, which update() re-issues.
#
# A vector of zeros is the model without a known variance under another
# spelling, so nothing is attached and the engines receive the call they
# would have received without it.

.FB_KNOWN_VARIANCE_COL <- "flexybayes_known_variance"
.FB_BRMS_KNOWN_SE_COL <- "flexybayes_known_se"
.FB_INLA_KNOWN_INDEX <- "flexybayes_known_index"
.FB_INLA_KNOWN_SCALE <- "flexybayes_known_scale"

.fb_known_variance_columns <- function() {
  c(
    .FB_KNOWN_VARIANCE_COL,
    .FB_BRMS_KNOWN_SE_COL,
    .FB_INLA_KNOWN_INDEX,
    .FB_INLA_KNOWN_SCALE
  )
}

# .fb_known_variance_check() --- refuse a vector that is not one finite,
# non-negative variance per row of `data`.
.fb_known_variance_check <- function(known_variance, data) {
  n <- nrow(data)
  values <- function(k) paste0(k, if (k == 1L) " value is" else " values are")
  problem <- if (!is.numeric(known_variance) || is.matrix(known_variance)) {
    paste0("a ", class(known_variance)[[1L]], " was given")
  } else if (length(known_variance) != n) {
    paste0(
      "it has ", length(known_variance), " values and `data` has ", n,
      " rows"
    )
  } else if (any(!is.finite(known_variance))) {
    paste0(values(sum(!is.finite(known_variance))), " missing or not finite")
  } else if (any(known_variance < 0)) {
    paste0(values(sum(known_variance < 0)), " negative")
  }
  if (is.null(problem)) {
    return(invisible(NULL))
  }
  stop(.fb_refusal_condition(
    reason_code = "known_variance_invalid",
    message = paste0(
      "`known_variance` must be a numeric vector with one finite, ",
      "non-negative variance for each row of `data`, and ", problem, ". ",
      "A row whose first-stage estimate carries no error takes 0."
    )
  ))
}

# .fb_known_variance_refuse() --- the refusal for a request that cannot carry
# a known variance, with the reason and the remedy in `detail`.
.fb_known_variance_refuse <- function(detail) {
  stop(.fb_refusal_condition(
    reason_code = "known_variance_not_supported",
    message = detail
  ))
}

# .fb_known_variance_attach() --- the vector as a column of the analysis data,
# and its record on the IR. A vector of zeros attaches nothing.
.fb_known_variance_attach <- function(fb, data, known_variance, aggregate) {
  if (!any(known_variance > 0)) {
    return(list(fb = fb, data = data))
  }
  fam <- fb$family %||% ""
  link_ok <- is.null(fb$link) || identical(fb$link, "identity")
  if (!identical(fam, "gaussian") || !link_ok) {
    .fb_known_variance_refuse(paste0(
      "`known_variance` adds a known variance to the residual variance of a ",
      "Gaussian response on the identity link, Var(y_i) = sigma^2 + v_i. ",
      "This model is family = \"", fam, "\"",
      if (!link_ok) paste0(", link = \"", fb$link, "\"") else "",
      ", which has no additive residual variance for it to join."
    ))
  }
  if (isTRUE(aggregate)) {
    .fb_known_variance_refuse(paste0(
      "`aggregate = TRUE` refused: the aggregated route compresses rows into ",
      "cell statistics, and a variance that differs from row to row does ",
      "not survive the compression. Leave `aggregate` at \"auto\", which ",
      "fits a Gaussian model row by row."
    ))
  }
  taken <- intersect(.fb_known_variance_columns(), names(data))
  if (length(taken) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "nested_index_name_collision",
      message = paste0(
        "`known_variance` is carried in columns named ",
        paste0("`", .fb_known_variance_columns(), "`", collapse = ", "),
        ", and `data` already has ",
        paste0("`", taken, "`", collapse = ", "),
        ". Rename ", if (length(taken) > 1L) "them" else "it", "."
      )
    ))
  }
  data[[.FB_KNOWN_VARIANCE_COL]] <- as.numeric(known_variance)
  fb$addition_terms <- c(
    fb$addition_terms %||% list(),
    list(list(
      type = "known_variance",
      column = .FB_KNOWN_VARIANCE_COL,
      values = as.numeric(known_variance)
    ))
  )
  list(fb = fb, data = data)
}

# .fb_known_variance_term() --- the IR's known-variance record, or NULL.
.fb_known_variance_term <- function(fb) {
  for (t in fb$addition_terms %||% list()) {
    if (identical(t$type %||% "", "known_variance")) {
      return(t)
    }
  }
  NULL
}

# .fb_known_variance_fill() --- zero on the `n_added` cells the design grid
# appended to the data (R/fb_complete_grid.R), which carry row 1's value in
# every column the model does not name. Such a cell has no response, and no
# first stage estimated it.
.fb_known_variance_fill <- function(fb, data, n_added) {
  if (is.null(.fb_known_variance_term(fb)) || n_added < 1L) {
    return(data)
  }
  rows <- nrow(data) - seq_len(n_added) + 1L
  data[[.FB_KNOWN_VARIANCE_COL]][rows] <- 0
  data
}

# .fb_known_variance_values() --- the known variance of each row of `data`,
# or NULL when the model carries none.
.fb_known_variance_values <- function(fb, data) {
  if (is.null(.fb_known_variance_term(fb))) {
    return(NULL)
  }
  as.numeric(data[[.FB_KNOWN_VARIANCE_COL]])
}

# .fb_inla_known_variance_formula_term() --- the INLA row effect.
.fb_inla_known_variance_formula_term <- function() {
  paste0(
    "f(", .FB_INLA_KNOWN_INDEX, ", ", .FB_INLA_KNOWN_SCALE,
    ", model = \"iid\", ",
    "hyper = list(prec = list(initial = 0, fixed = TRUE)))"
  )
}

# .fb_inla_known_variance_effect() --- the posterior mean of the row effect on
# each row, on the response scale, or NULL. INLA's linear predictor carries
# it and flexyBayes's fitted values do not.
.fb_inla_known_variance_effect <- function(object) {
  fb <- object$extras$fb_terms %||% object$fb
  v <- .fb_known_variance_values(fb, object$data)
  effect <- object$inla$summary.random[[.FB_INLA_KNOWN_INDEX]]
  if (is.null(v) || is.null(effect)) {
    return(NULL)
  }
  scale <- object$extras$response_standardisation$scale %||% 1
  # summary.random has been mapped back to the response scale, a factor of
  # `scale`, so each row adds w_i times its mean, w_i = sqrt(v_i) / scale.
  sqrt(v) / scale * effect$mean[match(seq_along(v), effect$ID)]
}


# --- fb_select() ------------------------------------------------------- #
#
# The fold covariance of the benchmark and optimiser workflows under
# `covariance = "base_model"` carries each row's known variance on V's
# diagonal, sigma^2 + v_i, with sigma^2 read off the fold refit, which is
# itself fitted with the known variances of the rows it is handed
# (R/selection_cov_model.R). Every other route leaves it out, and refuses: the
# REML plug-in estimates one residual variance common to every row, and the
# Bayesian workflow, the refit of each fold's chosen support and the
# projection's reference fit score or project posterior draws whose
# predictive variance omits it.

# .fb_sel_known_variance() --- each row's known variance, aligned to the
# base's data, or NULL when the base carries none.
.fb_sel_known_variance <- function(base) {
  fb <- base$fb %||% base$extras$fb_terms
  if (is.null(.fb_known_variance_term(fb))) {
    return(NULL)
  }
  v <- base$data[[.FB_KNOWN_VARIANCE_COL]]
  if (is.null(v) || length(v) != nrow(base$data)) {
    .fb_known_variance_refuse(paste0(
      "fb_select(): the base carries a known variance per row, and its data ",
      "has no `", .FB_KNOWN_VARIANCE_COL, "` column holding one value per ",
      "row. Refit the base with flexybayes(known_variance = ) on the data ",
      "fb_select() is to read."
    ))
  }
  as.numeric(v)
}

# .fb_sel_known_variance_refuse_route() --- a route that would leave the
# known variance out.
.fb_sel_known_variance_refuse_route <- function(route) {
  .fb_known_variance_refuse(paste0(
    "fb_select(): the base carries a known variance per row, which ", route,
    " would leave out. The penalised optimiser and the unpenalised ",
    "benchmark carry it into the fold covariance with ",
    "covariance = \"base_model\"."
  ))
}

# .fb_sel_known_variance_refuse_reml() --- the REML plug-in's refusal.
.fb_sel_known_variance_refuse_reml <- function() {
  .fb_sel_known_variance_refuse_route(paste0(
    "the REML plug-in (covariance = \"reml_plugin\"), which estimates one ",
    "residual variance common to every row,"
  ))
}

# .fb_select_require_known_variance_route() --- refuse, before any fit, a
# route that would leave a base's known variance out.
.fb_select_require_known_variance_route <- function(base,
                                                    workflow,
                                                    covariance,
                                                    value_request,
                                                    block_scope = FALSE) {
  if (is.null(.fb_sel_known_variance(base))) {
    return(invisible(NULL))
  }
  # The reference fit of a block-scoped prior carries the known variance in
  # its likelihood, and with no out-of-fold score nothing else is fitted.
  if (identical(workflow, "bayes") && isTRUE(block_scope)) {
    return(invisible(NULL))
  }
  if (identical(workflow, "bayes")) {
    .fb_sel_known_variance_refuse_route(paste0(
      "the Bayesian workflow, which scores each fold refit's posterior ",
      "predictive draws,"
    ))
  }
  if (identical(covariance, "reml_plugin")) {
    .fb_sel_known_variance_refuse_reml()
  }
  if (!is.null(value_request)) {
    .fb_sel_known_variance_refuse_route(paste0(
      "the Bayesian refit of each fold's chosen support (`refit`), whose ",
      "held-out log predictive density reads posterior draws,"
    ))
  }
  invisible(NULL)
}

# .fb_brms_known_se_newdata() --- brms reads the se() column of `newdata`
# even for a posterior mean, which does not depend on it. A new row's known
# variance is not part of its mean, so an absent column is supplied as zero.
.fb_brms_known_se_newdata <- function(object, newdata) {
  if (
    is.null(.fb_known_variance_term(.fb_fit_ir(object))) ||
      .FB_BRMS_KNOWN_SE_COL %in% names(newdata)
  ) {
    return(newdata)
  }
  newdata[[.FB_BRMS_KNOWN_SE_COL]] <- 0
  newdata
}
