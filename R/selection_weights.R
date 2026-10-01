# Selection weights, offsets and missing responses --- the absorption point
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md section 5.5, and
# oracle 8 of the correctness ledger in section 12.
#
# Section 5.5 is four sentences and every one of them is an "exactly
# once" rule:
#
#   observation weights live in R, and therefore in V, OR in the
#     provider loss, never both;
#   Gaussian offsets are subtracted before whitening and added back to
#     predictions;
#   only observed outcomes enter the selection objective, and augmented
#     or missing rows enter prediction or covariance blocks only under a
#     named estimand;
#   every transformation learned from data is learned on training rows.
#
# "Exactly once" cannot be asserted, only audited, and that is the whole
# design of this file. There is one absorption point --- the residual
# block R is formed in one expression, `.fb_sel_residual_block()`, and
# the offset is subtracted in one expression --- and the object carries
# enough of its own arithmetic that a second, independent function can
# recover from the MATRICES how many times each quantity was applied and
# refuse if the answer is not one. `.fb_sel_weight_audit()` reads V and
# the response, never a bookkeeping flag, because a flag records the
# intention and the matrices record what happened.
#
# The failure this guards against is not hypothetical. Whitening by a V
# that already carries the weights and then handing the same weights to
# a provider's `weights =` argument gives a loss weighted by w^2. It
# converges, it produces a plausible coefficient path, and nothing
# downstream can see it: the doubling is common to every fold, so
# resampling reports it as signal.
#
# Weight semantics is the package's own, not a second dialect:
# Var(y_i) = sigma^2 / w_i, the ASReml / lme4 / glm(weights =) precision
# sense that .refuse_unsupported_weights() in R/dispatch.R lowers for
# Gaussian identity only. Missing-response estimands reuse the augment /
# omit / fail vocabulary of R/na_action.R, read from
# .fb_na_action_choices() so the two cannot drift.
#
# THIS FILE DOES NOT OWN THE COVARIANCE. R/selection_covariance.R does.
# `.fb_sel_absorb_fit()` is the door a fit comes through and it reads the
# variance components, the Z designs and the applicability rule from that
# module rather than re-deriving any of them; `.fb_sel_absorb()` takes
# numbers only so the oracles can run on a known covariance without a
# model fit. The one thing this layer necessarily re-assembles is the
# residual block, because the weights belong in R and the block builder
# bakes a flat R into its Cholesky factors -- see `.fb_sel_absorb_fit()`
# for the full note and for the gate that holds the two together at unit
# weights.


# --- the single absorption point ------------------------------------- #

# .fb_sel_residual_block() --- R's diagonal, and the ONLY place in the
# selection module where an observation weight is allowed to touch a
# covariance.
#
# Var(y_i) = sigma^2 / w_i, so the residual variance of observation i is
# the base residual variance divided by its weight. Keeping this in one
# named expression is not decoration: the defect this file exists to
# catch is the same weight being applied at two sites, and two sites
# cannot arise from one function.
.fb_sel_residual_block <- function(residual_variance, weights) {
  residual_variance / weights
}


# .fb_sel_check_weights() --- weights are finite and strictly positive.
#
# Zero is refused rather than treated as a deletion. A zero-weight row
# has infinite residual variance, which is a row removed from the
# likelihood by a route that never appears in the row count, so the
# object would report n observations and fit fewer. If a row is to be
# dropped, it is dropped by name.
.fb_sel_check_weights <- function(weights, n) {
  if (!is.numeric(weights) || length(weights) != n) {
    stop(
      "`weights` must be a numeric vector of length ", n, ".",
      call. = FALSE
    )
  }
  bad <- which(!is.finite(weights) | weights <= 0)
  if (length(bad) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_weights_not_positive",
      message = paste0(
        "fb_select(): observation weight(s) at row(s) ",
        paste(utils::head(bad, 10L), collapse = ", "),
        if (length(bad) > 10L) ", ..." else "",
        " are not finite and positive. A zero weight deletes a row ",
        "without changing the row count and a missing weight is not a ",
        "weight; neither is absorbed silently."
      ),
      rows = bad
    ))
  }
  invisible(NULL)
}


# .fb_sel_base_weights() --- the observation weights the base was fitted
# with, one per row of `base$data`, or NULL. A constant vector is the
# unweighted model.
.fb_sel_base_weights <- function(base) {
  w <- .fb_ir_weights(base$fb %||% list()) %||% base$extras$call_info$weights
  if (.fb_weights_nonconstant(w)) as.numeric(w) else NULL
}


# The column of a weighted base's data that carries each row's weight into
# the fold refits.
.FB_SEL_WEIGHTS_COL <- ".fb_sel_weight"

# .fb_sel_with_weights_column() --- `data` with the base's weights as a
# column, or unchanged for an unweighted base.
.fb_sel_with_weights_column <- function(base, data) {
  w <- .fb_sel_base_weights(base)
  if (is.null(w)) {
    return(data)
  }
  if (length(w) != nrow(data)) {
    .fb_sel_refuse_weights_not_carried(paste0(
      "a base whose ", length(w), " weights do not match the ", nrow(data),
      " rows of its data"
    ))
  }
  data[[.FB_SEL_WEIGHTS_COL]] <- w
  data
}


# .fb_sel_refuse_weights_not_carried() --- a weighted base on a route that
# does not carry the weights refuses by name.
.fb_sel_refuse_weights_not_carried <- function(route) {
  stop(.fb_refusal_condition(
    reason_code = "selection_weights_not_carried",
    message = paste0(
      "fb_select(): the base carries observation weights, and ", route,
      " does not carry them into its fold fits and scores. Refit the base ",
      "without `weights` to select on the unweighted model."
    ),
    route = route
  ))
}


# .fb_sel_require_weights_by_route() --- `fb_resamples(weights_by = )` runs
# on a weighted base, on the benchmark and optimiser folds and the Bayesian
# refit of their supports, and refuses by name everywhere else.
.fb_sel_require_weights_by_route <- function(base, resamples, workflow,
                                             decision) {
  by <- resamples$weights_by
  .fb_sel_require_columns(base$data, by, "weights_by")
  if (is.null(.fb_sel_base_weights(base))) {
    stop(.fb_refusal_condition(
      reason_code = "selection_weights_by_unweighted_base",
      message = paste0(
        "fb_select(): `weights_by` re-estimates the base's observation ",
        "weights in each fold, and the base was fitted without `weights`. ",
        "Fit the base with the weights of the whole data first."
      )
    ))
  }
  route <- if (identical(workflow, "bayes")) {
    "the Bayesian workflow"
  } else if ((decision$stability_subsamples %||% 0L) > 0L) {
    "the stability subsamples"
  } else {
    NULL
  }
  if (!is.null(route)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_weights_not_carried",
      message = paste0(
        "fb_select(): ", route, " does not re-estimate the weights in each ",
        "fold, so `weights_by = ~ ", paste(by, collapse = " + "),
        "` would not be applied. Ask for workflow = \"penalized_optimizer\" ",
        "or \"unpenalized_benchmark\" with stability_subsamples = 0, or ",
        "drop `weights_by` to keep the base's weights in every fold."
      ),
      route = route
    ))
  }
  invisible(NULL)
}


# .fb_sel_fold_weights() --- one fold's observation weights, re-estimated by
# group from its training rows (`fb_resamples(weights_by = )`).
#
# s2_g is the mean squared residual of group g's training rows in an
# unweighted refit of the base, random effects included, and
# w_g = median(s2) / s2_g over the groups with training rows. Held-out rows,
# and rows of a group with no training rows, get w = 1.
.fb_sel_fold_weights <- function(base,
                                 split,
                                 by,
                                 provider,
                                 engine = NULL,
                                 seed = NULL) {
  spec <- .fb_sel_bayes_base_spec(base)
  fold_data <- spec$data
  # A constant weight is the unweighted model.
  fold_data[[.FB_SEL_WEIGHTS_COL]] <- 1
  if (identical(provider, "inla")) {
    fold_data[[spec$response]][
      c(split$assess_rows, split$not_estimable_rows)
    ] <- NA
    fitted_rows <- seq_len(nrow(fold_data))
  } else {
    fold_data <- fold_data[split$train_rows, , drop = FALSE]
    fitted_rows <- split$train_rows
  }
  fit <- .fb_sel_bayes_refit(
    base,
    fold_data,
    .fb_sel_bayes_fixed(
      spec$response,
      spec$protected_labels,
      character(0),
      spec$offsets
    ),
    provider,
    engine = engine,
    seed = seed,
    na_action = if (identical(provider, "inla")) "augment" else NULL,
    carry_known_variance = TRUE,
    fixed_alias = "drop",
    readout = "fitted"
  )
  y <- as.numeric(.fb_response_values(list(response = spec$response),
                                      spec$data))
  mu <- rep(NA_real_, nrow(spec$data))
  mu[fitted_rows] <- stats::fitted(fit)
  train <- split$train_rows
  key <- .fb_sel_key(spec$data, by)
  s2 <- vapply(
    split((y[train] - mu[train])^2, key[train]),
    mean,
    numeric(1L)
  )
  w_group <- stats::median(s2) / s2
  weights <- rep(1, nrow(spec$data))
  hit <- key %in% names(w_group)
  weights[hit] <- w_group[key[hit]]
  weights[split$assess_rows] <- 1
  .fb_sel_check_weights(weights, nrow(spec$data))
  list(
    weights = weights,
    record = list(by = by, s2 = s2, weights = w_group)
  )
}


# --- learned transformations ----------------------------------------- #

# .fb_sel_learn_scaling() --- the one transformation this layer learns
# from data, and the rows it is entitled to see.
#
# Weighted, because a weighted problem's centre is not its unweighted
# centre, and learned on `rows` only. `rows` is recorded on the result
# so `.fb_sel_weight_audit()` can recompute the numbers from the rows
# they claim to come from; a centre computed over the assessment rows
# and then labelled as training does not survive that recomputation.
#
# A column with no weighted spread over `rows` keeps scale 1 rather than
# being divided by zero --- an intercept column is the ordinary case,
# not an error.
.fb_sel_learn_scaling <- function(X, rows, weights) {
  Xr <- X[rows, , drop = FALSE]
  wr <- weights[rows]
  sw <- sum(wr)
  centre <- as.numeric(crossprod(Xr, wr) / sw)
  dev <- sweep(Xr, 2L, centre, "-")
  scale <- sqrt(as.numeric(crossprod(dev^2, wr) / sw))
  degenerate <- !is.finite(scale) | scale <= .Machine$double.eps^0.5
  scale[degenerate] <- 1
  names(centre) <- colnames(X)
  names(scale) <- colnames(X)
  list(
    rows_used = rows,
    n_rows_used = length(rows),
    centre = centre,
    scale = scale,
    constant_columns = colnames(X)[degenerate],
    weighted = TRUE
  )
}


# --- the absorbed problem -------------------------------------------- #

# .fb_sel_absorb() --- weights, offset and missing responses absorbed
# once each, with the arithmetic left on the object so a downstream
# consumer can check the count rather than take it on trust.
#
# Arguments
#   y                  response, NA where the outcome is not observed.
#   X                  design over all rows, protected and candidate
#                      columns together; the reduction of section 5.2
#                      splits them afterwards, on the whitened matrix.
#   V_random           sum_k sigma_k^2 Z_k Z_k', the random part of V
#                      with NO residual block and NO weights in it, and
#                      positive semi-definite. NULL means a
#                      residual-only V.
#   residual_variance  sigma^2, the base residual variance the weights
#                      scale.
#   weights            relative observation weights, Var(y_i) =
#                      sigma^2 / w_i. NULL is the unweighted model.
#   offset             Gaussian offset on the response scale. NULL is
#                      no offset.
#   absorption         where the weights go. "covariance" puts them in
#                      R and hands the provider an unweighted loss;
#                      "provider_loss" builds V unweighted and hands the
#                      weights on, and is admitted only where the two
#                      coincide (see below).
#   na_action          the named missing-response estimand, from
#                      .fb_na_action_choices(): "omit" analyses the
#                      observed cells, "augment" keeps the unobserved
#                      rows in the covariance blocks so section 5.4 can
#                      predict at them, "fail" refuses.
#   train_rows         rows any learned transformation may see. NULL
#                      means every row in the objective, which is the
#                      no-fold case.
#
# The provider-loss route and the residual-only restriction. Whitening
# by V and then minimising sum_i w_i r*_i^2 weights observation i by
# w_i / (v_rand_i + sigma^2), while weighted GLS weights it by
# 1 / (v_rand_i + sigma^2 / w_i). Those agree for every w only when
# v_rand is zero. Section 5.5's "or in the provider loss" is therefore
# an equivalence on a residual-only V and a different estimator
# otherwise, so the route refuses in the presence of a random part
# instead of quietly answering the other question.
.fb_sel_absorb <- function(
  y,
  X,
  V_random = NULL,
  residual_variance,
  weights = NULL,
  offset = NULL,
  absorption = c("covariance", "provider_loss"),
  na_action = "omit",
  train_rows = NULL
) {
  absorption <- match.arg(absorption)

  y <- as.numeric(y)
  X <- as.matrix(X)
  n <- length(y)
  if (nrow(X) != n) {
    stop("`X` has ", nrow(X), " rows and `y` has ", n, ".", call. = FALSE)
  }
  if (
    !is.numeric(residual_variance) ||
      length(residual_variance) != 1L ||
      !is.finite(residual_variance) ||
      residual_variance <= 0
  ) {
    stop(
      "`residual_variance` must be a single finite positive number.",
      call. = FALSE
    )
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(ncol(X)))
  }

  # ---- the named missing-response estimand --------------------------
  #
  # Read the vocabulary from R/na_action.R rather than restating it, so
  # a fourth policy added there cannot silently become an unnamed
  # estimand here.
  named <- .fb_na_action_choices()
  if (
    !is.character(na_action) ||
      length(na_action) != 1L ||
      !(na_action %in% named)
  ) {
    stop(.fb_refusal_condition(
      reason_code = "selection_missing_response_estimand_unnamed",
      message = paste0(
        "fb_select(): `",
        if (is.character(na_action) && length(na_action) == 1L) {
          na_action
        } else {
          "<not a single string>"
        },
        "` is not a named missing-response estimand. Rows whose ",
        "response is missing enter prediction or covariance blocks ",
        "only under a named estimand, and the named ones are: ",
        paste(named, collapse = ", "), "."
      ),
      supplied = na_action
    ))
  }

  missing_rows <- which(is.na(y))
  if (length(missing_rows) > 0L && identical(na_action, "fail")) {
    # The package already owns this refusal; a second code for the same
    # fact would be a second vocabulary.
    stop(.fb_refusal_condition(
      reason_code = "missing_response_refused",
      message = paste0(
        length(missing_rows), " observation(s) have a missing response ",
        "and na_action = \"fail\". Use na_action = \"augment\" to keep ",
        "them in the covariance blocks the prediction estimand reads, ",
        "or na_action = \"omit\" to analyse the observed cells only."
      ),
      family_class = "flexybayes_missing_response_refused"
    ))
  }

  objective_rows <- setdiff(seq_len(n), missing_rows)
  if (length(objective_rows) == 0L) {
    stop(
      "no observed response: the selection objective would be empty.",
      call. = FALSE
    )
  }
  # Only "augment" lets an unobserved row occupy a covariance block. It
  # never lets one into the objective; that is a separate rule and the
  # audit checks it separately.
  covariance_rows <- if (identical(na_action, "augment")) {
    seq_len(n)
  } else {
    objective_rows
  }

  # ---- the random part ----------------------------------------------
  if (is.null(V_random)) {
    V_random <- matrix(0, nrow = n, ncol = n)
  }
  V_random <- as.matrix(V_random)
  dimnames(V_random) <- NULL
  if (nrow(V_random) != n || ncol(V_random) != n) {
    stop("`V_random` must be ", n, " by ", n, ".", call. = FALSE)
  }

  # ---- weights and the route -----------------------------------------
  weights_declared <- if (is.null(weights)) rep(1, n) else as.numeric(weights)
  .fb_sel_check_weights(weights_declared, n)

  if (identical(absorption, "provider_loss")) {
    largest_random <- max(abs(V_random))
    if (largest_random > .Machine$double.eps^0.5) {
      stop(.fb_refusal_condition(
        reason_code = "selection_weights_loss_route_requires_residual_only",
        message = paste0(
          "fb_select(): the provider-loss weight route needs a ",
          "residual-only V; this V carries a random part with largest ",
          "entry ", format(largest_random, digits = 4), ". Weighting ",
          "the loss after whitening weights observation i by ",
          "w_i / (v_rand_i + sigma^2), while weighted GLS weights it ",
          "by 1 / (v_rand_i + sigma^2 / w_i); those agree only where ",
          "v_rand is zero. Use absorption = \"covariance\"."
        ),
        largest_random_entry = largest_random
      ))
    }
  }

  # The one absorption decision: the weights go into R, or they go to
  # the provider, and the two branches are the two arms of one `if` so
  # a third site cannot appear.
  if (identical(absorption, "covariance")) {
    r_diag <- .fb_sel_residual_block(residual_variance, weights_declared)
    provider_weights <- NULL
  } else {
    r_diag <- .fb_sel_residual_block(residual_variance, rep(1, n))
    provider_weights <- weights_declared
  }

  V_random_block <- V_random[covariance_rows, covariance_rows, drop = FALSE]
  V <- V_random_block +
    diag(r_diag[covariance_rows], nrow = length(covariance_rows))

  # ---- the offset, subtracted once, before whitening -----------------
  offset_declared <- if (is.null(offset)) rep(0, n) else as.numeric(offset)
  if (length(offset_declared) != n || any(!is.finite(offset_declared))) {
    stop(
      "`offset` must be a finite numeric vector of length ", n, ".",
      call. = FALSE
    )
  }
  y_working <- y[objective_rows] - offset_declared[objective_rows]
  X_working <- X[objective_rows, , drop = FALSE]

  # ---- whitening the objective block ---------------------------------
  #
  # Under "augment" the covariance block is larger than the objective,
  # so the objective's own sub-block is the one that is factored. V^-1
  # is never formed; W is a forward solve, as in
  # R/selection_covariance.R. A V_random that is not positive
  # semi-definite fails here as a caller error rather than a typed
  # refusal: the covariance builder upstream owns V's admissibility and
  # already refuses a block with no Cholesky factor by name.
  obj_in_cov <- match(objective_rows, covariance_rows)
  ch <- tryCatch(
    chol(V[obj_in_cov, obj_in_cov, drop = FALSE]),
    error = function(e) NULL
  )
  if (is.null(ch)) {
    # The covariance module's own code for this fact, registered in
    # R/selection_refusals.R. Weighting cannot make a positive-definite V
    # indefinite -- every weight is positive -- so reaching here means the
    # random part was inadmissible before the weights touched it, which is
    # that module's refusal to own, not a second one.
    stop(.fb_refusal_condition(
      reason_code = "selection_covariance_not_positive_definite",
      message = paste0(
        "fb_select(): the weighted covariance over the objective rows is ",
        "not positive definite, so it has no Cholesky factor. No jitter ",
        "was applied: a repair here would be silent misspecification."
      )
    ))
  }
  L <- t(ch)
  y_star <- as.numeric(forwardsolve(L, y_working))
  X_star <- forwardsolve(L, X_working)
  colnames(X_star) <- colnames(X)

  # ---- the learned transformation, on training rows only -------------
  if (is.null(train_rows)) {
    train_rows <- objective_rows
  }
  train_rows <- sort(unique(as.integer(train_rows)))
  learn_rows <- intersect(train_rows, objective_rows)
  if (length(learn_rows) == 0L) {
    stop(
      "`train_rows` contains no row with an observed response.",
      call. = FALSE
    )
  }
  preprocessing <- .fb_sel_learn_scaling(X, learn_rows, weights_declared)

  out <- structure(
    list(
      n_input = n,
      objective_rows = objective_rows,
      covariance_rows = covariance_rows,
      missing_response_rows = missing_rows,
      na_action = na_action,
      absorption = absorption,
      weights_declared = weights_declared,
      provider_weights = provider_weights,
      residual_variance = residual_variance,
      residual_block = r_diag,
      V_random_block = V_random_block,
      offset_declared = offset_declared,
      y_input = y,
      y_working = y_working,
      X_working = X_working,
      V = V,
      L = L,
      y_star = y_star,
      X_star = X_star,
      preprocessing = preprocessing,
      absorption_log = .fb_sel_absorption_log(
        absorption = absorption,
        weights_declared = weights_declared,
        offset_declared = offset_declared,
        na_action = na_action,
        n_missing = length(missing_rows),
        preprocessing = preprocessing
      )
    ),
    class = "fb_sel_absorbed"
  )

  # Audit at construction, so no object leaves this function unchecked.
  # The negative controls corrupt an object AFTER it is built, which is
  # the only way a double absorption arises in practice: the second
  # application is always somebody else's.
  .fb_sel_weight_audit(out)
  out
}


# .fb_sel_absorb_fit() --- the fit-facing constructor.
#
# `.fb_sel_absorb()` above takes the random part and sigma^2 as numbers so
# the oracles can exercise it on a known covariance without a model fit.
# That is right for the algebra and wrong as an interface: a caller who
# assembles their own V has forked the covariance, which is the one thing
# this module must not let happen. This is the door a real fit comes
# through, and it reads every input from R/selection_covariance.R rather
# than accepting one:
#
#   sigma^2 and the component variances  .fb_sel_variance_components()
#   the Z designs, matched BY NAME       .fb_sel_random_designs()
#   whether the base is admissible       .fb_sel_block_structure()
#
# so the estimand label, the E[sigma^2] convention, the name-matching of
# Z to theta and the applicability rule are that module's, once.
#
# WHAT IS STILL DUPLICATED, AND WHY. The random part is re-assembled here
# as sum_k sigma_k^2 Z_k Z_k' rather than reused from
# `.fb_sel_whitener()`, because that builder bakes a FLAT residual block
# into its per-block Cholesky factors and exposes no hook for a weighted
# R. Weights must enter R -- absorbing them anywhere downstream of the
# factorisation is the provider-loss route, which is a different
# estimator whenever V has a random part. The durable fix is a `weights`
# argument on `.fb_sel_whitener()` so R is formed once, in the builder
# that owns it; until that lands, `test-selection-weights.R` gates the
# two against each other at unit weights, where they must agree exactly.
# A duplicate that is checked against its original every run is a
# different thing from a fork, but it is not as good as not duplicating.
#
# Dense, and deliberately bounded. `max_n` is the same cap and the same
# reasoning as `.fb_sel_V_dense()`: at 73,341 observations a dense V is
# about 43 GB, and the block-structured route is the builder's, not this
# layer's.
.fb_sel_absorb_fit <- function(
  object,
  X,
  weights = NULL,
  offset = NULL,
  estimand = "plugin_mean_covariance",
  max_n = 4000L,
  ...
) {
  bs <- .fb_sel_block_structure(object)
  if (!isTRUE(bs$applicable)) {
    # Refuse on exactly the structures the covariance builder refuses on,
    # by asking it rather than by re-deriving the rule. A layer that
    # admitted a base the builder rejects would be answering a question
    # the substrate cannot pose.
    stop(.fb_refusal_condition(
      reason_code = "selection_covariance_not_buildable",
      message = paste0(
        "fb_select(): the absorption layer builds on the block-Cholesky ",
        "covariance and that builder does not apply because ", bs$reason,
        ". The random-term reduction is the route for this base."
      )
    ))
  }

  vc <- .fb_sel_variance_components(object, estimand)
  Z <- .fb_sel_random_designs(object)
  n <- nrow(object$data)
  if (n > max_n) {
    stop(.fb_refusal_condition(
      reason_code = "selection_dense_reference_too_large",
      message = paste0(
        "fb_select(): the absorption layer assembles a dense ", n, " by ",
        n, " covariance (about ",
        format(round(n^2 * 8 / 2^30, 1), nsmall = 1),
        " GB), above the cap of ", max_n, "."
      )
    ))
  }

  missing_var <- setdiff(names(Z), names(vc$terms))
  if (length(missing_var) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_covariance_not_buildable",
      message = paste0(
        "fb_select(): no variance component was reported for random ",
        "term(s) ", paste(missing_var, collapse = ", "),
        ". Z and theta must match by name."
      )
    ))
  }

  V_random <- matrix(0, nrow = n, ncol = n)
  for (nm in names(Z)) {
    V_random <- V_random + vc$terms[[nm]] * tcrossprod(Z[[nm]])
  }

  # The base model's own offset is the offset of this problem. A second one
  # declared on top of it would be subtracted twice, so the two are not
  # combined.
  base_offset <- .fb_offset_values(object$fb, object$data)
  if (!is.null(base_offset) && !is.null(offset)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_offset_absorbed_twice",
      message = paste0(
        "fb_select(): the offset does not round-trip; the base model already ",
        "carries ", paste0("`", object$fb$offsets, "`", collapse = ", "),
        " and a second offset was declared. The base's offset is absorbed ",
        "from the model; declare none here."
      )
    ))
  }
  offset <- offset %||% base_offset

  out <- .fb_sel_absorb(
    y = object$data[[object$fb$response]],
    X = X,
    V_random = V_random,
    residual_variance = vc$residual,
    weights = weights,
    offset = offset,
    ...
  )
  out$covariance_estimand <- vc$estimand
  out
}


# .fb_sel_absorption_log() --- where each quantity went, in words.
#
# This is the field that makes "exactly once" auditable by a reader
# rather than only by a function. It states the site, not merely the
# fact, because "weights: applied" is what a wrong object would also
# say.
.fb_sel_absorption_log <- function(
  absorption,
  weights_declared,
  offset_declared,
  na_action,
  n_missing,
  preprocessing
) {
  weights_site <- if (identical(absorption, "covariance")) {
    paste0(
      "R, and therefore V: R_ii = sigma^2 / w_i. The provider loss is ",
      "unweighted and `provider_weights` is NULL."
    )
  } else {
    paste0(
      "the provider loss: V is built with unit residual weights and ",
      "`provider_weights` carries w."
    )
  }
  data.frame(
    quantity = c("weights", "offset", "missing_response", "preprocessing"),
    absorbed = c(
      any(weights_declared != 1),
      any(offset_declared != 0),
      n_missing > 0L,
      TRUE
    ),
    site = c(
      weights_site,
      paste0(
        "subtracted from the response once, before whitening; added ",
        "back once by .fb_sel_add_offset()."
      ),
      paste0(
        "estimand `", na_action, "`: unobserved rows are excluded from ",
        "the objective under every estimand and occupy covariance ",
        "blocks only under `augment`."
      ),
      paste0(
        "weighted centre and scale learned on ",
        preprocessing$n_rows_used, " training row(s) only."
      )
    ),
    stringsAsFactors = FALSE
  )
}


# --- the audit (oracle 8) -------------------------------------------- #

# .fb_sel_weight_audit() --- recover from the object's MATRICES how many
# times each quantity was applied, and refuse if the answer is not one.
#
# The weight arithmetic. The absorbed covariance must be the random part
# plus a DIAGONAL residual block, so subtracting the random part leaves
#
#   R_implied = V - V_random,   R_implied_ii = sigma^2 / w_absorbed_i
#
# and w_absorbed = sigma^2 / diag(R_implied). Requiring R_implied to be
# diagonal is the off-diagonal half of the same check: a weight folded
# into a covariance rather than into R shows up there and nowhere in the
# diagonal. The total weighting the provider then sees is w_absorbed
# times whatever loss weight is attached, and under the contract that
# product is w exactly once. All three ways to break it land in the same
# number:
#
#   weights put in V twice          w_absorbed = w^2, loss 1   -> w^2
#   weights in V and in the loss    w_absorbed = w,   loss w   -> w^2
#   weights in neither              w_absorbed = 1,   loss 1   -> 1
#
# so one comparison catches over-absorption on either route and
# under-absorption as well.
#
# The covariance half reads V rather than the `absorption` field on
# purpose. A flag records what the constructor meant; the matrix records
# what it did, and the defect being hunted is exactly the case where
# those differ --- a doubled R leaves `absorption` reading "covariance"
# and `provider_weights` reading NULL, so every flag on the object still
# says the weights went in once. The loss half necessarily reads
# `provider_weights`, because a loss weight lives in no matrix: it IS the
# instruction to the provider, and the object is the only place it
# exists. That asymmetry is why the two are multiplied rather than
# compared: whichever side carries the extra copy, the product moves.
#
# Where the check has no power, and why that is correct. With w_i == 1
# for every i, w and w^2 coincide and a doubling is not detectable from
# any arithmetic --- because with unit weights there is nothing to
# double. The verdict's `power` field reports the largest departure of
# the declared weights from 1, which is the separation the check has to
# work with; a fixture with unit weights cannot falsify this gate and
# the test fixtures therefore use genuinely unequal weights.
.fb_sel_weight_audit <- function(x, tolerance = 1e-8) {
  if (!inherits(x, "fb_sel_absorbed")) {
    stop("`x` must be an <fb_sel_absorbed>.", call. = FALSE)
  }

  obj <- x$objective_rows
  cov_rows <- x$covariance_rows

  # ---- only observed outcomes enter the objective --------------------
  leaked <- intersect(obj, x$missing_response_rows)
  if (length(leaked) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_missing_response_in_objective",
      message = paste0(
        "fb_select(): ", length(leaked), " row(s) with no observed ",
        "response reached the selection objective (row(s) ",
        paste(utils::head(leaked, 10L), collapse = ", "),
        "). Only observed outcomes enter the objective, under every ",
        "estimand."
      ),
      rows = leaked
    ))
  }
  if (!identical(x$na_action, "augment")) {
    stray <- intersect(cov_rows, x$missing_response_rows)
    if (length(stray) > 0L) {
      stop(.fb_refusal_condition(
        reason_code = "selection_missing_response_estimand_unnamed",
        message = paste0(
          "fb_select(): ", length(stray), " unobserved row(s) occupy a ",
          "covariance block under estimand `", x$na_action, "`. An ",
          "augmented row enters a covariance block only under the ",
          "named estimand `augment`."
        ),
        rows = stray
      ))
    }
  }

  # ---- the weighting, recovered from V --------------------------------
  r_implied <- x$V - x$V_random_block
  r_diag <- diag(r_implied)
  off_diagonal <- r_implied
  diag(off_diagonal) <- 0
  scale_ref <- max(1, max(abs(r_diag)))
  if (max(abs(off_diagonal)) > tolerance * scale_ref) {
    stop(.fb_refusal_condition(
      reason_code = "selection_weights_absorbed_twice",
      message = paste0(
        "fb_select(): observation weights are absorbed off the residual ",
        "diagonal, not exactly once into R. V minus the random part ",
        "must be diagonal, and its largest off-diagonal entry here is ",
        format(max(abs(off_diagonal)), digits = 4), "."
      )
    ))
  }
  if (any(!is.finite(r_diag)) || any(r_diag <= 0)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_weights_absorbed_twice",
      message = paste0(
        "fb_select(): observation weights are absorbed in a way that ",
        "leaves a non-positive residual variance, not exactly once. V ",
        "minus the random part must be sigma^2 / w, which is positive ",
        "for every finite positive weight."
      )
    ))
  }
  w_in_covariance <- x$residual_variance / r_diag
  names(w_in_covariance) <- NULL

  w_in_loss <- if (is.null(x$provider_weights)) {
    rep(1, length(obj))
  } else {
    x$provider_weights[obj]
  }
  effective <- w_in_covariance[match(obj, cov_rows)] * w_in_loss
  declared <- x$weights_declared[obj]

  rel <- abs(effective / declared - 1)
  worst <- max(rel)
  if (worst > tolerance) {
    # Name the count, not just the mismatch. log(effective)/log(w) is
    # how many times the weight went in, and it reads 2 for either
    # doubling route and 0 for none, which tells the reader which
    # defect they have rather than only that they have one.
    varying <- declared != 1
    times <- if (any(varying)) {
      stats::median(log(effective[varying]) / log(declared[varying]))
    } else {
      NA_real_
    }
    how <- if (is.finite(times)) {
      paste0("about ", format(round(times, 3)), " time(s) over")
    } else {
      "an unrecoverable number of times over"
    }
    stop(.fb_refusal_condition(
      reason_code = "selection_weights_absorbed_twice",
      message = paste0(
        "fb_select(): observation weights are absorbed ", how,
        " the objective rows, not exactly once. The covariance and the ",
        "provider loss together imply a weighting whose largest ",
        "relative departure from the declared weights is ",
        format(worst, digits = 4), ", against a bar of ",
        format(tolerance), ". Weights live in R, and therefore in V, ",
        "or in the provider loss -- never both."
      ),
      worst = worst,
      times_absorbed = times
    ))
  }

  # ---- the offset, recovered from the response ------------------------
  implied_offset <- x$y_input[obj] - x$y_working
  declared_offset <- x$offset_declared[obj]
  off_scale <- max(1, max(abs(declared_offset)))
  off_worst <- max(abs(implied_offset - declared_offset)) / off_scale
  if (off_worst > tolerance) {
    applied <- if (any(declared_offset != 0)) {
      stats::median((implied_offset / declared_offset)[declared_offset != 0])
    } else {
      NA_real_
    }
    stop(.fb_refusal_condition(
      reason_code = "selection_offset_absorbed_twice",
      message = paste0(
        "fb_select(): the offset does not round-trip; the working ",
        "response differs from the input response by ",
        if (is.finite(applied)) {
          paste0("about ", format(round(applied, 3)), " times ")
        } else {
          "a quantity that is not "
        },
        "the declared offset (largest relative departure ",
        format(off_worst, digits = 4), "). It is subtracted once ",
        "before whitening and added back once to predictions."
      ),
      worst = off_worst,
      times_absorbed = applied
    ))
  }

  # ---- the learned transformation saw only its own rows ---------------
  #
  # Recomputed from `rows_used` alone. A centre computed over every row
  # and then filed under the training rows does not survive this, which
  # is what makes the fold-leak falsifier turn the gate red rather than
  # merely change a label. The drift scale is relative above 1 and
  # absolute below it, so a centre that happens to sit near zero cannot
  # hide a leak behind a vanishing denominator.
  pp <- x$preprocessing
  redone <- .fb_sel_learn_scaling(
    x$X_working[match(pp$rows_used, obj), , drop = FALSE],
    seq_along(pp$rows_used),
    x$weights_declared[pp$rows_used]
  )
  drift <- pmax(
    abs(redone$centre - pp$centre) / pmax(1, abs(pp$centre)),
    abs(redone$scale - pp$scale) / pmax(1, abs(pp$scale))
  )
  if (any(drift > tolerance)) {
    bad <- names(drift)[drift > tolerance]
    stop(.fb_refusal_condition(
      reason_code = "selection_preprocessing_leak",
      message = paste0(
        "fb_select(): recorded preprocessing does not reproduce from ",
        "its own rows: ", paste(utils::head(bad, 10L), collapse = ", "),
        ". The recorded centre and scale claim ", pp$n_rows_used,
        " training row(s) and do not follow from them, so the ",
        "transformation saw rows it was not entitled to."
      ),
      columns = bad
    ))
  }

  invisible(list(
    weights_in_covariance = w_in_covariance,
    weights_in_loss = w_in_loss,
    effective_weights = effective,
    declared_weights = declared,
    worst_relative_weight_error = worst,
    worst_relative_offset_error = off_worst,
    worst_relative_preprocessing_drift = max(drift),
    power = max(abs(declared - 1)),
    tolerance = tolerance
  ))
}


# --- solving, and adding the offset back ------------------------------ #

# .fb_sel_solve() --- the unpenalised solve on the absorbed problem.
#
# Two properties matter more than the arithmetic. It audits before it
# solves, so a double-absorbed object cannot reach a solver through this
# door; and the loss it minimises is unweighted exactly when the weights
# are already in V, which is section 5.5's first sentence expressed as
# control flow rather than as a comment.
.fb_sel_solve <- function(x) {
  .fb_sel_weight_audit(x)
  if (is.null(x$provider_weights)) {
    beta <- qr.solve(x$X_star, x$y_star)
  } else {
    rw <- sqrt(x$provider_weights[x$objective_rows])
    beta <- qr.solve(x$X_star * rw, x$y_star * rw)
  }
  names(beta) <- colnames(x$X_star)
  fitted_working <- as.numeric(x$X_working %*% beta)
  list(
    coefficients = beta,
    fitted_working = fitted_working,
    fitted_response = .fb_sel_add_offset(x, fitted_working)
  )
}


# .fb_sel_add_offset() --- the add-back half of the offset rule.
#
# Marks its own output. Adding the offset to a value that already
# carries it is the prediction-side twin of weighting an already-
# whitened loss, and the mark is what lets the second call refuse
# instead of returning a plausible number: a doubled offset shifts every
# fitted value and changes no coefficient, so no coefficient-level check
# can see it.
.fb_sel_add_offset <- function(x, eta) {
  if (!inherits(x, "fb_sel_absorbed")) {
    stop("`x` must be an <fb_sel_absorbed>.", call. = FALSE)
  }
  if (isTRUE(attr(eta, "fb_sel_offset_applied"))) {
    stop(.fb_refusal_condition(
      reason_code = "selection_offset_absorbed_twice",
      message = paste0(
        "fb_select(): the offset does not round-trip; these values ",
        "already carry it. A Gaussian offset is subtracted once before ",
        "whitening and added back once to predictions."
      )
    ))
  }
  if (length(eta) != length(x$objective_rows)) {
    stop(
      "`eta` must have one value per objective row (",
      length(x$objective_rows), ").",
      call. = FALSE
    )
  }
  out <- as.numeric(eta) + x$offset_declared[x$objective_rows]
  attr(out, "fb_sel_offset_applied") <- TRUE
  out
}


#' @export
print.fb_sel_absorbed <- function(x, ...) {
  cat("<fb_sel_absorbed>\n")
  cat("  rows in       : ", x$n_input, "\n", sep = "")
  cat(
    "  objective     : ", length(x$objective_rows),
    " (observed outcomes only)\n",
    sep = ""
  )
  cat(
    "  covariance    : ", length(x$covariance_rows),
    " under estimand `", x$na_action, "`\n",
    sep = ""
  )
  cat("  weights in    : ", x$absorption, "\n", sep = "")
  cat(
    "  offset        : ",
    if (any(x$offset_declared != 0)) {
      "subtracted once, before whitening"
    } else {
      "none"
    },
    "\n",
    sep = ""
  )
  cat(
    "  preprocessing : weighted centre and scale on ",
    x$preprocessing$n_rows_used, " training row(s)\n",
    sep = ""
  )
  invisible(x)
}
