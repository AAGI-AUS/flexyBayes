# Selection covariance from the base model --- the fold covariance the
# benchmark and optimiser workflows whiten by, read off a refit of the base
# model on each fold's training rows.
#
# Why this file exists beside R/selection_resamples.R. The default fold
# covariance comes from `.fb_sel_reml_components()`, an independent REML
# estimator that knows random intercepts and two-way interactions and forms V
# densely. `fb_select(covariance = "base_model")` reads V from the model the
# analyst fitted instead: in each training fold the base model is refitted
# through flexyBayes with its own formula, random and residual structure,
# priors and engine, and V is assembled from that refit's variance
# components.
#
# The point estimate. Each variance enters V at its posterior mean,
# E[sigma^2] = E[sigma]^2 + Var(sigma), read off the refit's variance-component
# table, whose `estimate` and `std.error` are the posterior mean and standard
# deviation of the standard deviation. V is linear in the variances, so for a
# structure without a correlation parameter this plug-in is the posterior
# mean of V itself, the estimand `.fb_sel_variance_components()` names
# `plugin_mean_covariance`. The square of the posterior mean of sigma would be
# smaller by the posterior variance of sigma, which is largest where a
# component is weakly identified, and that is where a wrong V does the most
# damage. The two correlations of an autoregressive field enter at their
# posterior means. V is not linear in them, so there the plug-in is the
# posterior mean of each parameter rather than of V, and the record says so.
#
# The fold refit. INLA is handed every row of the data with the assessment
# responses set to missing (`na_action = "augment"`): a row without a
# response adds nothing to the likelihood, so the hyperparameter posterior is
# the training rows', while every level, lattice node and known-matrix row
# stays where the base put it. brms is handed the training rows.
#
# Structures. V and its blocks are read through fb_marginal_cov()
# (R/marginal_cov.R), which names the structures it builds. This route also
# refuses a sectioned residual, and anything fb_marginal_cov() refuses, with
# `selection_covariance_not_buildable`.
#
# Blocking. V is formed densely only while one n x n double matrix over the
# data is within `options(flexyBayes.selection_dense_gb)`, default 0.5 GiB
# (8,192 rows). The dense route's peak memory is about five times that
# matrix (review/units/U-COV-MODEL_2026-09-19, probe 02), and its time grows
# with the cube of the rows. Above the bound V is factored block by block
# over the grouping every slot nests within, choosing the one with the most
# blocks, and no n x n matrix is formed; a V with no such grouping refuses
# with `selection_dense_covariance_too_large`. A base with one field per trial
# is factored by trial at any size. Below the bound both sources
# take the dense route, so one set of components gives one answer, bit for
# bit, whichever source supplied it.


# --- the dense bound ------------------------------------------------------ #

.FB_SEL_DENSE_GB <- 0.5

# .fb_sel_dense_gb() --- the bound, in GiB, on one dense n x n matrix.
.fb_sel_dense_gb <- function() {
  gb <- getOption("flexyBayes.selection_dense_gb", .FB_SEL_DENSE_GB)
  if (!is.numeric(gb) || length(gb) != 1L || !is.finite(gb) || gb <= 0) {
    stop(
      "`options(flexyBayes.selection_dense_gb = )` must be a single positive ",
      "number of GiB.",
      call. = FALSE
    )
  }
  gb
}

# .fb_sel_dense_admitted() --- is a dense n x n double matrix within it.
.fb_sel_dense_admitted <- function(n) {
  8 * as.numeric(n)^2 <= .fb_sel_dense_gb() * 2^30
}

# .fb_sel_require_dense() --- refuse, before any fold runs, a covariance that
# would be formed densely above the bound. `why` says what forms it densely
# and what else there is.
.fb_sel_require_dense <- function(n, why) {
  if (.fb_sel_dense_admitted(n)) {
    return(invisible(NULL))
  }
  bound <- .fb_sel_dense_gb()
  stop(.fb_refusal_condition(
    reason_code = "selection_dense_covariance_too_large",
    message = paste0(
      "fb_select(): the fold covariance over these ", n, " rows would be ",
      "formed as a dense ", n, " x ", n, " matrix of ",
      format(round(8 * as.numeric(n)^2 / 2^30, 2), nsmall = 2),
      " GiB, above the bound of ", bound, " GiB set by ",
      "options(flexyBayes.selection_dense_gb = ). ", why, " A larger bound ",
      "admits the dense route at its cost: about five times that matrix in ",
      "memory, and time that grows with the cube of the rows."
    ),
    n_rows = n,
    bound_gb = bound
  ))
}


# --- the model, before any fold ------------------------------------------ #

# .fb_sel_marginal_cov() --- fb_marginal_cov(), its refusals restated as this
# route's.
.fb_sel_marginal_cov <- function(fit, ...) {
  tryCatch(
    fb_marginal_cov(fit, ...),
    flexybayes_refusal_marginal_cov_not_buildable = function(e) {
      stop(.fb_refusal_condition(
        reason_code = "selection_covariance_not_buildable",
        message = paste0(
          "fb_select(): ",
          sub("^fb_marginal_cov\\(\\): ", "", conditionMessage(e))
        ),
        terms = e$terms
      ))
    }
  )
}

# .fb_sel_refuse_cov_terms() --- the refusal naming the terms this route does
# not build although fb_marginal_cov() does.
.fb_sel_refuse_cov_terms <- function(terms, detail = NULL) {
  stop(.fb_refusal_condition(
    reason_code = "selection_covariance_not_buildable",
    message = paste0(
      "fb_select(): covariance = \"base_model\" builds the fold ",
      "covariance from random intercepts, two-way interactions, ",
      "uncorrelated random slopes, diag(f):g, an ar1(row):ar1(col) field ",
      "on one trial or id(trial):ar1(row):ar1(col) on each trial, and vm() ",
      "or ped() known covariances, over ",
      "an independent residual, and this base also carries: ",
      paste(unique(terms), collapse = ", "),
      ".",
      if (!is.null(detail)) paste0(" ", detail)
    ),
    terms = unique(terms)
  ))
}

# .fb_sel_cov_model() --- everything the base-model covariance needs that no
# fold changes: the structure of V, the engine, and how every fold will be
# factored. The dense bound is checked here, before any refit.
.fb_sel_cov_model <- function(base, engine = NULL, seed = NULL) {
  base <- .fb_sel_bayes_normalise_base(base)
  provider <- .fb_sel_bayes_provider(base)
  n <- nrow(base$data)
  marginal <- .fb_sel_marginal_cov(base, variances = FALSE)
  if (!is.null(marginal$residual$by)) {
    .fb_sel_refuse_cov_terms("residual at_units")
  }
  slots <- marginal$blocks
  blocking <- marginal$groups

  # At 48 trials of 72 plots the dense route's peak was 2.1 GB against 0.9 GB
  # by trial (review/units/U-AR_2026-09-24).
  field_per_trial <- isTRUE(blocking$applicable) && any(vapply(
    slots,
    function(s) !is.null(s$structure) && !is.null(s$key),
    logical(1L)
  ))
  # A known covariance with off-diagonal entries links rows across any
  # grouping, so V is formed densely (R/known_variance.R).
  known_covariance <- marginal$known_covariance
  if (!is.null(known_covariance) &&
        any(known_covariance[upper.tri(known_covariance)] != 0)) {
    .fb_sel_require_dense(n, paste0(
      "The known covariance in `known_variance` correlates rows across ",
      "every grouping, so the covariance cannot be factored block by block."
    ))
    field_per_trial <- FALSE
    blocking <- list(applicable = FALSE)
  }
  whitening <- if (.fb_sel_dense_admitted(n) && !field_per_trial) {
    list(
      mode = "dense",
      group = NA_character_,
      n_blocks = NA_integer_,
      largest_block = NA_integer_
    )
  } else if (isTRUE(blocking$applicable)) {
    list(
      mode = "by_block",
      group = blocking$column,
      n_blocks = length(blocking$blocks),
      largest_block = max(lengths(blocking$blocks)),
      blocks = blocking$blocks
    )
  } else {
    .fb_sel_require_dense(n, paste0(
      "The covariance of this base is not block-diagonal by any grouping of ",
      "its random terms (", blocking$reason, "), so it cannot be factored ",
      "block by block."
    ))
  }
  whitening$bound_gb <- .fb_sel_dense_gb()

  structure(
    list(
      base = base,
      provider = provider,
      engine = engine,
      seed = seed,
      n = n,
      slots = slots,
      # The structure of V with no variances, evaluated at each fold's.
      marginal = marginal,
      # Each row's known variance, added to the residual on V's diagonal in
      # every fold, or NULL (R/known_variance.R).
      known_variance = marginal$residual$known_variance,
      # The known covariance over the base's rows, whose off-diagonal part
      # V carries beside the diagonal above, or NULL.
      known_covariance = known_covariance,
      # Each row's observation weight, dividing the residual on V's diagonal,
      # or NULL.
      weights = marginal$residual$weights,
      builder = if (isTRUE(blocking$applicable)) {
        "block_cholesky"
      } else {
        "dense_cholesky"
      },
      whitening = whitening,
      engine_version = as.character(utils::packageVersion(
        if (identical(provider, "inla")) "INLA" else "brms"
      ))
    ),
    class = "fb_sel_cov_model"
  )
}


# --- one fold ---------------------------------------------------------- #

# .fb_sel_cov_fold_refit() --- the base refitted on a fold, and the posterior
# mean of each of its variances.
.fb_sel_cov_fold_refit <- function(base,
                                   split,
                                   provider,
                                   engine = NULL,
                                   seed = NULL,
                                   slots = list()) {
  spec <- .fb_sel_bayes_base_spec(base)
  fold_data <- spec$data
  if (!is.null(split$weights)) {
    fold_data[[.FB_SEL_WEIGHTS_COL]] <- split$weights
  }
  if (identical(provider, "inla")) {
    fold_data[[spec$response]][
      c(split$assess_rows, split$not_estimable_rows)
    ] <- NA
  } else {
    fold_data <- fold_data[split$train_rows, , drop = FALSE]
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
    # The candidate space has already refused, or dropped, every protected
    # column aliased on this fold's training rows.
    fixed_alias = "drop",
    readout = "varcomp"
  )
  out <- .fb_sel_cov_fold_values(fit, slots)
  out$rows_passed <- nrow(fold_data)
  out
}

# .fb_sel_cov_fold_values() --- the fold refit's residual and slot variances,
# and an AR1 x AR1 field's correlations, read through fb_marginal_cov().
.fb_sel_cov_fold_values <- function(fit, slots) {
  marginal <- .fb_sel_marginal_cov(fit)
  if (!is.null(marginal$residual$by)) {
    .fb_sel_refuse_cov_terms("residual at_units")
  }
  terms <- numeric(0)
  rho <- list()
  for (s in slots) {
    b <- marginal$blocks[[s$name]]
    if (is.null(b)) {
      .fb_sel_refuse_cov_terms(
        s$label,
        "The fold refit carries no variance for it."
      )
    }
    terms[[s$name]] <- b$variance
    if (!is.null(b$rho)) {
      rho[[s$name]] <- b$rho
    }
  }
  list(residual = marginal$residual$variance[[1L]], terms = terms, rho = rho)
}

# .fb_sel_cov_fold() --- one fold's components and the covariance they give,
# with the record the certificate carries.
.fb_sel_cov_fold <- function(model, split, f) {
  fit_seed <- as.integer(model$seed %||% 0L) + 1000L * as.integer(f)
  started <- proc.time()
  components <- .fb_sel_cov_fold_refit(
    model$base,
    split,
    provider = model$provider,
    engine = model$engine,
    seed = fit_seed,
    slots = model$slots
  )
  seconds <- unname((proc.time() - started)[["elapsed"]])
  if (!is.null(split$weights)) {
    model$weights <- split$weights
  }
  # A held-out row is scored at w = 1: its own weight is not known from the
  # training rows.
  if (!is.null(model$weights)) {
    model$weights[split$assess_rows] <- 1
  }
  list(
    spec = .fb_sel_cov_fold_spec(model, components),
    record = list(
      fold = split$label,
      residual = components$residual,
      terms = components$terms,
      rho = components$rho,
      rows_passed = components$rows_passed %||% NA_integer_,
      seed = fit_seed,
      seconds = seconds
    )
  )
}

# .fb_sel_cov_fold_spec() --- the components as a covariance specification.
#
# On the dense route every slot is a dense Z over all rows, a structure matrix
# folded in as Z L with G = L L', so `.fb_sel_fold_blocks()` reads it as it
# reads the REML plug-in's. On the block route the slots travel as they are
# and each block is assembled by `.fb_sel_block_covariance()`.
.fb_sel_cov_fold_spec <- function(model, components) {
  estimand <- "plugin_mean_covariance"
  mode <- model$whitening$mode
  if (identical(mode, "dense")) {
    z <- list()
    for (s in model$slots) {
      z[[s$name]] <- .fb_sel_slot_dense(s, model$n, components)
    }
    spec <- .fb_sel_cov_spec(
      z,
      components$residual,
      components$terms[names(z)],
      estimand,
      known_variance = model$known_variance,
      weights = model$weights,
      known_covariance = model$known_covariance
    )
  } else {
    names_slots <- vapply(model$slots, function(s) s$name, character(1L))
    spec <- .fb_sel_cov_spec(
      list(),
      components$residual,
      components$terms[names_slots],
      estimand,
      known_variance = model$known_variance,
      weights = model$weights
    )
    spec$slots <- model$slots
    spec$rho <- components$rho
    spec$blocks <- model$whitening$blocks
    spec$marginal <- fb_marginal_cov(model$marginal, variances = components)
    spec$marginal$residual$weights <- model$weights
  }
  spec$whitening <- mode
  spec
}

# .fb_sel_slot_dense() --- one slot's Z over all rows; with a structure matrix,
# Z L where G = L L'. Z is an indicator, so Z L is L's rows read by index.
.fb_sel_slot_dense <- function(slot, n, components) {
  z <- matrix(0, nrow = n, ncol = slot$n_levels)
  keep <- which(!is.na(slot$index))
  if (is.null(slot$structure)) {
    z[cbind(keep, slot$index[keep])] <- if (is.null(slot$value)) {
      1
    } else {
      slot$value[keep]
    }
    return(z)
  }
  g <- if (identical(slot$structure$type, "ar1xar1")) {
    rho <- components$rho[[slot$name]]
    kronecker(
      .fb_sel_ar1_correlation(slot$structure$n_row, rho[["row"]]),
      .fb_sel_ar1_correlation(slot$structure$n_col, rho[["col"]])
    )
  } else {
    slot$structure$k
  }
  factor_g <- tryCatch(chol(g), error = function(e) NULL)
  if (is.null(factor_g)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_covariance_not_positive_definite",
      message = paste0(
        "fb_select(): the structure matrix of ", slot$label, " is not ",
        "positive definite, so it has no Cholesky factor. No jitter was ",
        "applied."
      )
    ))
  }
  # With one field per trial, G = I (x) g and its factor is I (x) chol(g).
  n_nodes <- nrow(g)
  node <- (slot$index[keep] - 1L) %% n_nodes + 1L
  trial <- (slot$index[keep] - 1L) %/% n_nodes
  for (a in unique(trial)) {
    at <- which(trial == a)
    z[keep[at], a * n_nodes + seq_len(n_nodes)] <-
      t(factor_g)[node[at], , drop = FALSE]
  }
  z
}

# .fb_sel_ar1_correlation() --- the AR1 correlation matrix over n ordered
# levels.
.fb_sel_ar1_correlation <- function(n, rho) {
  rho^abs(outer(seq_len(n), seq_len(n), "-"))
}


# --- block by block -------------------------------------------------------- #

# .fb_sel_block_covariance() --- the covariance between two sets of rows of
# one block, assembled by fb_marginal_cov() at the fold's variances and
# weights. It sums the same terms in the same order as the dense route: an
# indicator slot's entries equal the dense V's exactly, a per-trial field's
# to rounding.
.fb_sel_block_covariance <- function(cov_spec, rows_a, rows_b) {
  as.matrix(cov_spec$marginal, rows = rows_a, cols = rows_b)
}

# .fb_sel_fold_blocks_by_group() --- `.fb_sel_fold_blocks()` block by block.
#
# V_TT, V_SS and V_ST are block-diagonal over the grouping, so each block's
# training rows are factored on their own, the gain V_ST V_TT^-1 is taken
# within each block, and the conditional covariance of the assessment rows is
# held per block. Positions index the fold's training and assessment rows in
# the order the fold lists them, as the dense route's matrices do.
.fb_sel_fold_blocks_by_group <- function(cov_spec,
                                         train_rows,
                                         assess_rows,
                                         prediction = c("new_cluster",
                                                        "known_cluster"),
                                         label = "1") {
  prediction <- match.arg(prediction)
  n <- max(train_rows, assess_rows, unlist(cov_spec$blocks, use.names = FALSE))
  train_pos <- integer(n)
  train_pos[train_rows] <- seq_along(train_rows)
  assess_pos <- integer(n)
  assess_pos[assess_rows] <- seq_along(assess_rows)

  factors <- list()
  conditionals <- list()
  gains <- list()
  logdet <- 0
  for (rows in cov_spec$blocks) {
    t_rows <- rows[train_pos[rows] > 0L]
    t_rows <- t_rows[order(train_pos[t_rows])]
    s_rows <- rows[assess_pos[rows] > 0L]
    s_rows <- s_rows[order(assess_pos[s_rows])]
    factor_b <- NULL
    if (length(t_rows) > 0L) {
      factor_b <- tryCatch(
        chol(.fb_sel_block_covariance(cov_spec, t_rows, t_rows)),
        error = function(e) NULL
      )
      if (is.null(factor_b)) {
        stop(.fb_refusal_condition(
          reason_code = "selection_fold_covariance_not_factorisable",
          message = paste0(
            "fb_select(): a block of the training covariance of fold ",
            label, " is not positive definite, so it has no Cholesky ",
            "factor. No jitter was applied: a repair here is a ",
            "misspecification carried into every score this fold contributes."
          )
        ))
      }
      logdet <- logdet + 2 * sum(log(diag(factor_b)))
      factors[[length(factors) + 1L]] <- list(
        positions = train_pos[t_rows],
        lower = t(factor_b)
      )
    }
    if (length(s_rows) > 0L) {
      v_ss <- .fb_sel_block_covariance(cov_spec, s_rows, s_rows)
      if (identical(prediction, "known_cluster") && !is.null(factor_b)) {
        v_st <- .fb_sel_block_covariance(cov_spec, s_rows, t_rows)
        gain <- t(backsolve(factor_b, forwardsolve(t(factor_b), t(v_st))))
        v_ss <- v_ss - tcrossprod(gain, v_st)
        v_ss <- (v_ss + t(v_ss)) / 2
        gains[[length(gains) + 1L]] <- list(
          assess_positions = assess_pos[s_rows],
          train_positions = train_pos[t_rows],
          gain = gain
        )
      }
      conditionals[[length(conditionals) + 1L]] <- list(
        positions = assess_pos[s_rows],
        v = v_ss
      )
    }
  }

  block_of <- integer(length(assess_rows))
  within <- integer(length(assess_rows))
  for (b in seq_along(conditionals)) {
    p <- conditionals[[b]]$positions
    block_of[p] <- b
    within[p] <- seq_along(p)
  }

  list(
    factor_tt = structure(
      list(blocks = factors, n = length(train_rows)),
      class = "fb_sel_block_factor"
    ),
    v_s_given_t = structure(
      list(
        blocks = conditionals,
        block_of = block_of,
        within = within,
        n = length(assess_rows)
      ),
      class = "fb_sel_block_conditional"
    ),
    gain = if (identical(prediction, "known_cluster")) {
      structure(list(blocks = gains), class = "fb_sel_block_gain")
    },
    prediction = prediction,
    logdet_v_tt = logdet,
    mode = "by_block"
  )
}

# .fb_sel_whiten_train() --- L^-1 M for the fold's training rows, with L the
# lower Cholesky factor of V_TT, dense or block by block.
.fb_sel_whiten_train <- function(factor_tt, m) {
  m <- as.matrix(m)
  if (!inherits(factor_tt, "fb_sel_block_factor")) {
    return(forwardsolve(t(factor_tt), m))
  }
  out <- matrix(0, nrow = nrow(m), ncol = ncol(m))
  for (b in factor_tt$blocks) {
    out[b$positions, ] <- forwardsolve(
      b$lower,
      m[b$positions, , drop = FALSE]
    )
  }
  out
}

# .fb_sel_conditional_rows() --- V_S|T over some of the assessment rows.
.fb_sel_conditional_rows <- function(v_s_given_t, index) {
  if (!inherits(v_s_given_t, "fb_sel_block_conditional")) {
    return(v_s_given_t[index, index, drop = FALSE])
  }
  out <- matrix(0, nrow = length(index), ncol = length(index))
  block <- v_s_given_t$block_of[index]
  for (b in unique(block)) {
    at <- which(block == b)
    w <- v_s_given_t$within[index[at]]
    out[at, at] <- v_s_given_t$blocks[[b]]$v[w, w, drop = FALSE]
  }
  out
}

# .fb_sel_conditional_variance() --- the diagonal of V_S|T.
.fb_sel_conditional_variance <- function(v_s_given_t) {
  if (!inherits(v_s_given_t, "fb_sel_block_conditional")) {
    return(diag(v_s_given_t))
  }
  out <- numeric(v_s_given_t$n)
  for (b in v_s_given_t$blocks) {
    out[b$positions] <- diag(b$v)
  }
  out
}


# --- the record -------------------------------------------------------- #

# .fb_sel_cv_covariance_record() --- the fold covariance as the certificate
# states it, read from the model and the fold objects.
.fb_sel_cv_covariance_record <- function(
  cov_model,
  substrates,
  weighted = FALSE,
  weights_by = NULL
) {
  folds <- unname(lapply(substrates, function(s) s$covariance_record))
  equivalence_reason <- paste0(
    "V is held fixed at each fold's plug-in, so every point on the path is ",
    "the exact posterior mode, given that V, of the base's structure under ",
    "the declared prior; the Bayesian workflow integrates V instead"
  )
  if (is.null(cov_model)) {
    return(list(
      source = "reml_plugin",
      point_estimate = paste0(
        "restricted maximum likelihood on each fold's training rows, for ",
        "random intercepts and two-way interactions, profiled over the ",
        "variance ratios with the residual variance in closed form"
      ),
      engine = NA_character_,
      engine_version = NA_character_,
      known_variance = FALSE,
      # Whether each row's residual was divided by its observation weight.
      weights = weighted,
      # The groups whose weights each fold re-estimated, or NULL.
      weights_by = weights_by,
      builder = NA_character_,
      whitening = list(mode = "dense", bound_gb = .fb_sel_dense_gb()),
      equivalence_level = "map_exact_fixed_covariance",
      equivalence_reason = equivalence_reason,
      folds = folds
    ))
  }
  list(
    source = "base_model",
    point_estimate = paste0(
      "the posterior mean of each variance, E[sigma^2] = estimate^2 + ",
      "std.error^2 from the refit's variance-component table, which for a ",
      "structure without a correlation parameter is the posterior mean of ",
      "V; an autoregressive field's correlations at their posterior means"
    ),
    engine = cov_model$provider,
    engine_version = cov_model$engine_version,
    # Whether each row's known variance was added to the residual on V's
    # diagonal (R/known_variance.R).
    known_variance = !is.null(cov_model$known_variance),
    weights = weighted,
    weights_by = weights_by,
    builder = cov_model$builder,
    whitening = cov_model$whitening[
      c("mode", "group", "n_blocks", "largest_block", "bound_gb")
    ],
    equivalence_level = "map_exact_fixed_covariance",
    equivalence_reason = equivalence_reason,
    folds = folds
  )
}
