# ordinal_methods.R -- methods for ordinal fits (class flexybayes_ordinal).
#
# An ordinal fit reports its coefficients and thresholds on the latent scale
# and its category probabilities from posterior draws: fitted() at the fitted
# random effects, predict() at random effects of zero or integrated over
# their variances. The engine's own threshold rows (the fixed intercept of
# the cumulative INLA route, the threshold rows of the continuation-ratio
# expansion, brms's Intercept[k]) are removed from every coefficient
# surface. A method that reads a mean response, a residual or a per-row
# likelihood refuses by name.


# Draws ------------------------------------------------------------------------

#' Evaluate under a seed when one is given
#'
#' @param seed An integer seed or `NULL`.
#' @param expr The expression.
#' @returns The value of `expr`.
#' @noRd
#' @keywords internal
.fb_ordinal_seeded <- function(seed, expr) {
  if (is.null(seed)) {
    return(expr)
  }
  .fb_with_seed(seed, expr)
}

#' Validate a number of posterior draws
#'
#' @param n_draws The requested number.
#' @returns An integer.
#' @noRd
#' @keywords internal
.fb_ordinal_check_n_draws <- function(n_draws) {
  if (
    !is.numeric(n_draws) ||
      length(n_draws) != 1L ||
      !is.finite(n_draws) ||
      n_draws < 1
  ) {
    stop("`n_draws` must be a single positive number of draws.", call. = FALSE)
  }
  as.integer(n_draws)
}

#' Rows of a latent draw matrix, by name
#'
#' @param latent A matrix, latent rows by draws.
#' @param wanted The row names wanted.
#' @returns A matrix, draws by `wanted`.
#' @noRd
#' @keywords internal
.fb_ordinal_latent_rows <- function(latent, wanted) {
  idx <- match(wanted, rownames(latent))
  if (anyNA(idx)) {
    stop(
      "flexyBayes internal error: the INLA draws carry no row ",
      paste(wanted[is.na(idx)][seq_len(min(3L, sum(is.na(idx))))],
        collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  t(latent[idx, , drop = FALSE])
}

#' Joint posterior draws of an INLA ordinal fit
#'
#' @param object The ordinal fit.
#' @param n_draws The number of draws.
#' @param seed An integer seed or `NULL`.
#' @param rows User rows whose linear predictor is wanted.
#' @returns A list: `tau` (draws by thresholds), `beta` (draws by
#'   coefficients, named as coef() names them), `eta` (draws by `rows`,
#'   thresholds excluded, or `NULL`) and `variance` (the summed variance of
#'   the random terms in each draw).
#' @noRd
#' @keywords internal
.fb_ordinal_draws_inla <- function(object, n_draws, seed, rows = integer(0)) {
  ord <- object$extras$ordinal
  fit <- object$inla
  if (is.null(fit$misc$configs)) {
    stop(.fb_refusal_condition(
      reason_code = "ordinal_draws_unavailable",
      message = paste0(
        "This ordinal fit kept no posterior configurations to draw category ",
        "probabilities from, as a refit with readout = \"varcomp\" or ",
        "\"fitted\" does not. Refit it with fb_refit(readout = \"draws\") ",
        "or in full."
      )
    ))
  }
  samples <- .fb_ordinal_seeded(
    seed,
    .fb_inla_posterior_sample(
      object,
      n_draws,
      seed = if (is.null(seed)) 0L else as.integer(seed),
      num.threads = "1:1"
    )
  )
  latent <- do.call(cbind, lapply(samples, function(z) z$latent[, 1L]))
  rownames(latent) <- rownames(samples[[1L]]$latent)
  hyper <- do.call(rbind, lapply(samples, function(z) z$hyperpar))

  # Thresholds, and the part of the stored predictor that is not eta -------
  K <- ord$K
  if (identical(ord$route, "pom")) {
    tau <- .fb_ordinal_pom_alpha(hyper, K, scale = "user")
    shift <- .fb_ordinal_latent_rows(latent, "(Intercept):1")[, 1L]
  } else {
    cut <- .fb_ordinal_latent_rows(latent, paste0(ord$engine_rows, ":1"))
    tau <- -cut
    # The stored predictor of a user row is its first boundary row,
    # cut_1 + eta.
    shift <- cut[, 1L]
  }

  # Coefficients and linear predictor ----------------------------------
  coef_rows <- setdiff(rownames(fit$summary.fixed), ord$engine_rows)
  beta <- .fb_ordinal_latent_rows(latent, paste0(coef_rows, ":1"))
  colnames(beta) <- .inla_restore_term_labels(object$level_labels, coef_rows)
  eta <- if (length(rows)) {
    .fb_ordinal_latent_rows(
      latent,
      paste0("Predictor:", ord$user_rows[rows])
    ) - shift
  }
  # A model with no random term and no threshold hyperparameter (the
  # continuation ratio with fixed effects only) carries no hyperparameter.
  precision <- grep("^Precision for ", colnames(hyper), value = TRUE)
  variance <- if (length(precision)) {
    rowSums(1 / hyper[, precision, drop = FALSE])
  } else {
    rep(0, length(samples))
  }
  list(tau = tau, beta = beta, eta = eta, variance = variance)
}

#' Posterior draws of a brms ordinal fit
#'
#' @param object The ordinal fit.
#' @param n_draws The number of draws, at most the number stored.
#' @param seed An integer seed or `NULL`, for the choice of draws.
#' @returns A list as `.fb_ordinal_draws_inla()` returns, `eta` `NULL`.
#' @noRd
#' @keywords internal
.fb_ordinal_draws_brms <- function(object, n_draws, seed) {
  ord <- object$extras$ordinal
  m <- .fb_brms_draws_matrix(object$brms)
  keep <- if (n_draws < nrow(m)) {
    .fb_ordinal_seeded(seed, sort(sample.int(nrow(m), n_draws)))
  } else {
    seq_len(nrow(m))
  }
  m <- m[keep, , drop = FALSE]
  cn <- colnames(m)
  tau_cols <- paste0("b_Intercept[", seq_len(ord$K - 1L), "]")
  b_cols <- setdiff(.brms_mean_b_cols(cn, .brms_dpar_names(object$brms)),
    tau_cols)
  beta <- as.matrix(m[, b_cols, drop = FALSE])
  colnames(beta) <- sub("^b_", "", b_cols)
  sd_cols <- grep("^sd_.*__Intercept$", cn, value = TRUE)
  list(
    tau = unname(as.matrix(m[, tau_cols, drop = FALSE])),
    beta = beta,
    eta = NULL,
    variance = if (length(sd_cols)) {
      rowSums(as.matrix(m[, sd_cols, drop = FALSE])^2)
    } else {
      rep(0, nrow(m))
    }
  )
}

#' Posterior draws of an ordinal fit, on either engine
#'
#' @param object The ordinal fit.
#' @param n_draws,seed As for the engine readers.
#' @returns A list as `.fb_ordinal_draws_inla()` returns.
#' @noRd
#' @keywords internal
.fb_ordinal_draws <- function(object, n_draws, seed) {
  if (identical(object$extras$ordinal$engine, "brms")) {
    return(.fb_ordinal_draws_brms(object, n_draws, seed))
  }
  .fb_ordinal_draws_inla(object, n_draws, seed)
}

#' Posterior mean category probabilities over draws
#'
#' @param eta A matrix, draws by rows.
#' @param tau A matrix, draws by thresholds.
#' @param variance The random-effect variance to integrate over in each
#'   draw, or zeros.
#' @param family,link The family and link.
#' @param summary `TRUE` for the mean over draws, `FALSE` for every draw.
#' @returns A matrix, rows by K, or an array, draws by rows by K.
#' @noRd
#' @keywords internal
.fb_ordinal_probs_over_draws <- function(
  eta,
  tau,
  variance,
  family,
  link,
  summary = TRUE
) {
  n_draws <- nrow(eta)
  K <- ncol(tau) + 1L
  rule <- .fb_gauss_hermite()
  each <- if (summary) NULL else array(NA_real_, c(n_draws, ncol(eta), K))
  total <- 0
  for (s in seq_len(n_draws)) {
    p <- .fb_ordinal_probs_integrated(
      eta[s, ],
      tau[s, ],
      variance[[s]],
      family,
      link,
      rule
    )
    if (summary) {
      total <- total + p
    } else {
      each[s, , ] <- p
    }
  }
  if (summary) total / n_draws else each
}


#' Canonical posterior draws of an ordinal fit, thresholds included
#'
#' The canonical draws of the fit with the engine's threshold rows removed
#' and the thresholds added under their own names, read from the same draws.
#'
#' @param object The ordinal fit.
#' @param n_draws The number of draws (INLA).
#' @returns A named list of draw vectors.
#' @noRd
#' @keywords internal
.fb_ordinal_canonical_draws <- function(object, n_draws) {
  ord <- object$extras$ordinal
  K <- ord$K
  raw <- if (identical(ord$engine, "brms")) {
    fb_as_draws_simple(object)
  } else {
    fb_as_draws_simple(object, n_samples = n_draws)
  }
  tau <- if (identical(ord$engine, "brms")) {
    do.call(cbind, raw[paste0("b_Intercept[", seq_len(K - 1L), "]")])
  } else if (identical(ord$route, "pom")) {
    .fb_ordinal_pom_alpha(raw, K, scale = "user")
  } else {
    -do.call(cbind, raw[paste0(ord$engine_rows, ":1")])
  }
  draws <- .fb_canonical_draws(object, drop = TRUE, raw = raw)
  engine <- c(
    ord$engine_rows,
    paste0("b_Intercept[", seq_len(K - 1L), "]")
  )
  draws <- draws[!names(draws) %in% engine]
  names_tau <- .fb_ordinal_threshold_names(ord$levels)
  for (k in seq_len(K - 1L)) {
    draws[[names_tau[[k]]]] <- unname(tau[, k])
  }
  draws
}


# Design and random-term support ------------------------------------------------

#' The coefficient design of new rows
#'
#' @param object The ordinal fit.
#' @param newdata A data frame.
#' @param coef_names The coefficient names, as coef() gives them.
#' @returns A matrix, rows by `coef_names`.
#' @noRd
#' @keywords internal
.fb_ordinal_design <- function(object, newdata, coef_names) {
  trms <- .fb_fixef_terms(object)
  X <- .fb_fixef_model_matrix(
    trms,
    newdata,
    c("(Intercept)", coef_names),
    .fb_fit_data(object)
  )
  X[, coef_names, drop = FALSE]
}

#' Refuse an integrated prediction over a random term it does not cover
#'
#' A new unit's random effects are integrated as one normal deviation whose
#' variance is the sum of the fitted variances. That holds for random
#' intercepts, their interactions and an AR1 field, whose marginal variance
#' is its fitted variance.
#'
#' @param object The ordinal fit.
#' @returns Invisibly `TRUE`.
#' @noRd
#' @keywords internal
.fb_ordinal_check_marginal <- function(object) {
  types <- vapply(
    .fb_fit_ir(object)$random_terms %||% list(),
    function(t) t$type %||% "",
    character(1L)
  )
  covered <- c("simple", "nested", "combo", "ar1", "ar1_spatial")
  off <- unique(types[!types %in% covered])
  if (length(off)) {
    stop(.fb_refusal_condition(
      reason_code = "predict_new_levels_structure_not_supported",
      message = paste0(
        "predict(marginal = TRUE) on an ordinal fit integrates random ",
        "intercepts, their interactions and an AR1 field, and this model ",
        "has a random term of type ",
        paste(off, collapse = ", "),
        ". Predict with marginal = FALSE, which sets the random effects to ",
        "zero."
      )
    ))
  }
  invisible(TRUE)
}


# Fitted and predicted probabilities --------------------------------------------

#' Category probabilities of an ordinal fit
#'
#' @description
#' `fitted()` returns, for each row the model was fitted to, the posterior
#' mean probability of each category, at the fitted values of the random
#' effects. `predict()` returns the same probabilities for the rows of
#' `newdata`, with the random effects set to zero, or with
#' `marginal = TRUE` integrated over their fitted variances, which is the
#' probability for a new unit in a new block.
#'
#' @details
#' Each probability is averaged over `n_draws` joint posterior draws of the
#' thresholds, the coefficients and the variances, so it carries the
#' uncertainty of all three. With `marginal = TRUE` the random terms of a
#' new unit enter as one normal deviation whose variance is the sum of the
#' fitted variances (random intercepts, their interactions and an AR1
#' field, whose marginal variance is its fitted variance), and the
#' probability is integrated over it by 40-node Gauss-Hermite quadrature in
#' each draw. On brms, `fitted()` reads [brms::posterior_epred()].
#'
#' The probability of a set of categories is the sum of their columns, for
#' example `rowSums(p[, c("Z5", "Z6")])` for a stage of Z5 or later.
#'
#' @param object An ordinal fit, from `flexybayes()` with
#'   `family = "cumulative"` or `"cratio"`.
#' @param newdata A data frame holding the fixed-effect variables, or
#'   `NULL` (the default) for the rows the model was fitted to.
#' @param marginal `FALSE` (the default) sets every random effect to zero.
#'   `TRUE` integrates over them.
#' @param n_draws The number of posterior draws. Default `1000`. On brms it
#'   is at most the number the sampler stored.
#' @param seed An integer seed for the draws, or `NULL` (the default) to use
#'   the session's random number stream.
#' @param summary `TRUE` (the default) returns the posterior mean. `FALSE`
#'   returns every draw.
#' @param ... Ignored by `fitted()`. In `predict()`, an argument of the
#'   other `predict()` methods that has no meaning here is refused by name.
#'
#' @returns A matrix with one row per row of the data and one column per
#'   category, named by the category labels. With `summary = FALSE`, an
#'   array, draws by rows by categories.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(
#'     trt = factor(rep(c("a", "b"), each = 60)),
#'     block = factor(rep(1:6, times = 20))
#'   )
#'   eta <- 0.8 * (d$trt == "b") + rnorm(6, 0, 0.3)[d$block]
#'   d$score <- cut(eta + rlogis(120), c(-Inf, -0.5, 0.5, 1.5, Inf),
#'     labels = c("low", "mid", "high", "top"), ordered_result = TRUE)
#'   fit <- flexybayes(score ~ trt, random = ~block, data = d,
#'     family = "cumulative", backend = "inla", verbose = FALSE)
#'   head(fitted(fit, seed = 1L))
#'   predict(fit, newdata = data.frame(trt = c("a", "b")), marginal = TRUE,
#'     seed = 1L)
#' }
#' }
#' @name ordinal_probabilities
#' @export
fitted.flexybayes_ordinal <- function(
  object,
  n_draws = 1000L,
  seed = NULL,
  summary = TRUE,
  ...
) {
  ord <- object$extras$ordinal
  n_draws <- .fb_ordinal_check_n_draws(n_draws)
  if (identical(ord$engine, "brms")) {
    ep <- .fb_ordinal_seeded(seed, brms::posterior_epred(
      object$brms,
      ndraws = min(n_draws, brms::ndraws(object$brms))
    ))
    out <- if (isTRUE(summary)) apply(ep, c(2L, 3L), mean) else unname(ep)
  } else {
    dr <- .fb_ordinal_draws_inla(
      object,
      n_draws,
      seed,
      rows = seq_len(nrow(object$data))
    )
    out <- .fb_ordinal_probs_over_draws(
      dr$eta,
      dr$tau,
      rep(0, nrow(dr$eta)),
      ord$family,
      ord$link,
      summary = isTRUE(summary)
    )
  }
  .fb_ordinal_label_probs(out, ord$levels)
}

#' @rdname ordinal_probabilities
#' @export
predict.flexybayes_ordinal <- function(
  object,
  newdata = NULL,
  marginal = FALSE,
  n_draws = 1000L,
  seed = NULL,
  summary = TRUE,
  ...
) {
  .fb_ordinal_refuse_predict_args(list(...))
  ord <- object$extras$ordinal
  n_draws <- .fb_ordinal_check_n_draws(n_draws)
  if (isTRUE(marginal)) {
    .fb_ordinal_check_marginal(object)
  }
  if (is.null(newdata)) {
    newdata <- .fb_fit_data(object)
  }
  dr <- .fb_ordinal_draws(object, n_draws, seed)
  X <- .fb_ordinal_design(object, newdata, colnames(dr$beta))
  eta <- dr$beta %*% t(X)
  variance <- if (isTRUE(marginal)) dr$variance else rep(0, nrow(eta))
  out <- .fb_ordinal_probs_over_draws(
    eta,
    dr$tau,
    variance,
    ord$family,
    ord$link,
    summary = isTRUE(summary)
  )
  .fb_ordinal_label_probs(out, ord$levels)
}

#' Name the category dimension of a probability matrix or array
#'
#' @param p A matrix (rows by K) or an array (draws by rows by K).
#' @param levels The category labels.
#' @returns `p`, with the category dimension named.
#' @noRd
#' @keywords internal
.fb_ordinal_label_probs <- function(p, levels) {
  if (is.matrix(p)) {
    dimnames(p) <- list(NULL, levels)
  } else {
    dimnames(p) <- list(NULL, NULL, levels)
  }
  p
}

#' Refuse a predict() argument that has no meaning for an ordinal fit
#'
#' @param dots The `...` of predict(), as a list.
#' @returns Invisibly `TRUE`.
#' @noRd
#' @keywords internal
.fb_ordinal_refuse_predict_args <- function(dots) {
  why <- c(
    type = "it returns category probabilities.",
    se.fit = "summary = FALSE returns every posterior draw.",
    classify = "it predicts the rows of `newdata`.",
    allow_new_levels = paste0(
      "a new unit's random effects are integrated with marginal = TRUE."
    ),
    re_formula = paste0(
      "marginal = FALSE sets the random effects to zero and marginal = TRUE ",
      "integrates over them."
    ),
    mean_by = "it predicts the rows of `newdata`.",
    level = "summary = FALSE returns every posterior draw."
  )
  named <- intersect(names(dots), names(why))
  if (length(named)) {
    .fb_refuse_method_argument(
      "predict() on an ordinal fit",
      named[[1L]],
      why[[named[[1L]]]]
    )
  }
  invisible(TRUE)
}


# Coefficients and tables -------------------------------------------------------

#' Drop the engine's threshold rows from a coefficient table
#'
#' @param object The ordinal fit.
#' @param names Coefficient names.
#' @returns A logical vector, `TRUE` for the user's coefficients.
#' @noRd
#' @keywords internal
.fb_ordinal_keep <- function(object, names) {
  !names %in% object$extras$ordinal$engine_rows
}

#' @rdname coef.flexybayes
#' @export
coef.flexybayes_ordinal <- function(
  object,
  what = c("fixed", "thresholds", "random", "missing", "all"),
  ...
) {
  what <- match.arg(what)
  thresholds <- object$extras$ordinal$thresholds
  tau <- stats::setNames(thresholds$estimate, thresholds$term)
  if (identical(what, "thresholds")) {
    return(tau)
  }
  parent <- if (inherits(object, "flexybayes_inla")) {
    coef.flexybayes_inla
  } else {
    coef.flexybayes
  }
  out <- parent(object, what = what)
  if (identical(what, "fixed")) {
    return(out[.fb_ordinal_keep(object, names(out))])
  }
  if (identical(what, "all")) {
    out$fixed <- out$fixed[.fb_ordinal_keep(object, names(out$fixed))]
    out$thresholds <- tau
  }
  out
}

#' @export
vcov.flexybayes_ordinal <- function(object, ...) {
  v <- NextMethod()
  keep <- .fb_ordinal_keep(object, rownames(v))
  v[keep, keep, drop = FALSE]
}

#' @export
confint.flexybayes_ordinal <- function(object, ...) {
  ci <- NextMethod()
  ci[.fb_ordinal_keep(object, rownames(ci)), , drop = FALSE]
}

#' @export
tidy.flexybayes_ordinal <- function(
  x,
  conf.int = TRUE,
  conf.level = 0.95,
  effects = c("fixed", "random"),
  ...
) {
  effects <- match.arg(effects)
  parent <- if (inherits(x, "flexybayes_inla")) {
    tidy.flexybayes_inla
  } else {
    tidy.flexybayes
  }
  out <- parent(
    x,
    conf.int = conf.int,
    conf.level = conf.level,
    effects = effects,
    ...
  )
  if (identical(effects, "random")) {
    return(out)
  }
  out <- out[.fb_ordinal_keep(x, out$term), , drop = FALSE]
  thresholds <- x$extras$ordinal$thresholds
  if (!isTRUE(conf.int)) {
    thresholds <- thresholds[c("term", "estimate", "std.error")]
  }
  out <- rbind(out, thresholds[names(out)])
  rownames(out) <- NULL
  out
}

#' @export
summary.flexybayes_ordinal <- function(object, ...) {
  out <- NextMethod()
  ord <- object$extras$ordinal
  out$thresholds <- ord$thresholds
  out$ordinal <- list(
    family = ord$family,
    link = ord$link,
    levels = ord$levels,
    dropped = ord$dropped,
    threshold_prior = ord$threshold_prior
  )
  class(out) <- c("summary.flexybayes_ordinal", class(out))
  out
}

#' @export
print.summary.flexybayes_ordinal <- function(x, digits = 4L, ...) {
  NextMethod()
  .fb_ordinal_print_block(x$ordinal, NULL, digits)
  invisible(x)
}

#' @export
print.flexybayes_ordinal <- function(x, ...) {
  NextMethod()
  ord <- x$extras$ordinal
  .fb_ordinal_print_block(ord, ord$thresholds, 4L)
  invisible(x)
}

#' The ordinal block of print() and summary()
#'
#' @param ord The ordinal record (or its summary subset).
#' @param thresholds The threshold table to print, or `NULL`.
#' @param digits Rounding.
#' @returns Invisibly `NULL`.
#' @noRd
#' @keywords internal
.fb_ordinal_print_block <- function(ord, thresholds, digits) {
  cat(
    "\n-- Ordinal response ",
    strrep("-", 42),
    "\n",
    sep = ""
  )
  cat(
    "  ",
    ord$family,
    " / ",
    ord$link,
    ", ",
    length(ord$levels),
    " categories: ",
    paste(ord$levels, collapse = " < "),
    "\n",
    sep = ""
  )
  cat(
    "  Thresholds: ",
    if (identical(ord$family, "cumulative")) {
      "P(Y <= k) = F(threshold_k - eta)"
    } else {
      "P(Y > k | Y >= k) = F(eta - threshold_k)"
    },
    ", under ",
    ord$threshold_prior,
    "\n",
    sep = ""
  )
  if (length(ord$dropped)) {
    cat(
      "  Levels with no observation, dropped: ",
      paste(ord$dropped, collapse = ", "),
      "\n",
      sep = ""
    )
  }
  if (!is.null(thresholds)) {
    body <- thresholds
    num <- vapply(body, is.numeric, logical(1L))
    body[num] <- lapply(body[num], round, digits = digits)
    print(body, row.names = FALSE)
  }
  invisible(NULL)
}

#' @export
model.matrix.flexybayes_ordinal <- function(object, ...) {
  X <- NextMethod()
  X[, colnames(X) != "(Intercept)", drop = FALSE]
}


# Refusals ----------------------------------------------------------------------

#' Refuse a method an ordinal fit has no answer for
#'
#' @param method The method, as the user called it.
#' @param instead What answers the question on an ordinal fit.
#' @returns Does not return.
#' @noRd
#' @keywords internal
.fb_ordinal_refuse_method <- function(method, instead) {
  stop(.fb_refusal_condition(
    reason_code = "ordinal_method_not_supported",
    message = paste0(
      method,
      " reads a mean response, a residual or a per-row likelihood, which ",
      "an ordinal fit does not have. ",
      instead
    ),
    method = method
  ))
}

#' @export
residuals.flexybayes_ordinal <- function(object, ...) {
  .fb_ordinal_refuse_method(
    "residuals()",
    "fitted() returns each row's category probabilities."
  )
}

#' @export
augment.flexybayes_ordinal <- function(x, ...) {
  .fb_ordinal_refuse_method(
    "augment()",
    "fitted() returns each row's category probabilities."
  )
}

#' @export
plot.flexybayes_ordinal <- function(x, type = NULL, ...) {
  if (is.null(type)) {
    type <- if (.fb_has_sampler_draws(x)) "diagnostics" else "effects"
  }
  if (type %in% c("residuals", "variogram")) {
    .fb_ordinal_refuse_method(
      paste0("plot(type = \"", type, "\")"),
      "The displays \"effects\", \"variance\" and \"blups\" read the latent scale."
    )
  }
  plot.flexybayes(x, type = type, ...)
}

#' @exportS3Method emmeans::recover_data
recover_data.flexybayes_ordinal <- function(object, ...) {
  .fb_ordinal_refuse_method(
    "emmeans",
    paste0(
      "predict() returns category probabilities for a reference grid passed ",
      "as `newdata`, and coef() and vcov() the latent-scale coefficients."
    )
  )
}

#' @exportS3Method emmeans::emm_basis
emm_basis.flexybayes_ordinal <- function(object, ...) {
  .fb_ordinal_refuse_method(
    "emmeans",
    paste0(
      "predict() returns category probabilities for a reference grid passed ",
      "as `newdata`, and coef() and vcov() the latent-scale coefficients."
    )
  )
}

#' @exportS3Method marginaleffects::get_predict
get_predict.flexybayes_ordinal <- function(model, ...) {
  .fb_ordinal_refuse_method(
    "marginaleffects",
    "predict() returns category probabilities for the rows of `newdata`."
  )
}
