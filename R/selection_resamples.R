# Selection resampling --- folds, the prediction estimand, and the scoring
# unit for fb_select().
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md sections 3.3, 4.1, 4.2,
# 5.4, 9 and 12. This file carries oracles 4, 7 and 9 of the correctness
# ledger, and none of them is trusted before it has been seen red on its own
# defect.
#
#   Oracle 4  no fold leakage. Folds fall on original units, before any
#             rotation, and every quantity learned from data is learned on
#             the training rows. The gate is a falsifier rather than a
#             bookkeeping check: the substrate is rebuilt with the assessment
#             rows corrupted, and anything the fold declares as learned that
#             moves was not learned on the training rows.
#   Oracle 7  the prediction estimand. V_TT is built and factored; V_ST
#             V_TT^-1 is a solve, never an inverse; and for a genuinely new
#             cluster V_ST is zero, so the marginal fixed part is what the
#             fold is scored on.
#   Oracle 9  the scoring unit. The standard-error unit is the independent
#             cluster, not the row and not the correlated fold score. One
#             standard-error selection is not reproducible without it, which
#             is why the definition sits in the contract.
#
# Why the order of section 4.2 is not an implementation detail. An eigen or
# Cholesky rotation mixes observations within a covariance block, so a row of
# the rotated design is a contrast among original rows rather than an
# observation. Splitting rotated rows therefore leaves cluster mates on both
# sides of the split under a different name, and the resulting score answers
# an easier question than the one asked. Folds are assigned first, on the
# original units, and every rotation happens inside a fold.

# --- the resampling specification ------------------------------------ #
#
# `fb_resamples()` is the public constructor and lives in R/fb_resamples.R,
# beside the other classed inputs. Everything below reads the object it
# returns; nothing below builds one.


# .fb_sel_variables() --- a one-sided formula or a character vector, resolved
# to plain column names. Anything else is an argument-shape error rather than
# a refusal: there is no user-facing condition to report, only a wrong type.
.fb_sel_variables <- function(x, arg) {
  if (is.null(x)) {
    return(character(0))
  }
  if (inherits(x, "formula")) {
    if (length(x) != 2L) {
      stop(
        "`",
        arg,
        "` must be a one-sided formula, as in ~ LOCATION.",
        call. = FALSE
      )
    }
    return(all.vars(x))
  }
  if (is.character(x)) {
    return(x)
  }
  stop(
    "`",
    arg,
    "` must be a one-sided formula or a character vector of ",
    "column names.",
    call. = FALSE
  )
}


# .fb_sel_count() --- a single whole number at or above a floor.
.fb_sel_count <- function(x, arg, minimum) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      x != as.integer(x) ||
      x < minimum
  ) {
    stop(
      "`",
      arg,
      "` must be a single whole number of at least ",
      minimum,
      ".",
      call. = FALSE
    )
  }
  as.integer(x)
}


# --- reproducibility -------------------------------------------------- #

# .fb_sel_with_seed() --- evaluate under a named seed and put the caller's
# random stream back exactly as it was found.
#
# The restoration matters more than it looks. A fold assignment that silently
# advanced the session's stream would make everything downstream of an
# fb_select() call depend on whether fb_select() had been run, which is the
# kind of dependence that shows up months later as an irreproducible figure.
.fb_sel_with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(
      suppressWarnings(rm(
        list = ".Random.seed",
        envir = globalenv()
      )),
      add = TRUE
    )
  }
  set.seed(seed)
  force(expr)
}


# --- fold assignment on original units -------------------------------- #

# .fb_sel_key() --- the composite level of one or more columns, as a
# character vector. `paste` with a separator rather than `interaction()`,
# because interaction() carries the full crossed level set and the empty
# cells would be dealt into folds as if they were units.
.fb_sel_key <- function(data, vars) {
  parts <- lapply(vars, function(v) as.character(data[[v]]))
  do.call(paste, c(parts, list(sep = ":")))
}


# .fb_sel_folds() --- resolve a specification against the data into splits of
# ORIGINAL row indices.
#
# Normative (contract section 4.2 step 1): assignment happens here, on the
# original units, before any covariance is built and before any rotation. The
# indices a split carries index the rows of `data` as the user supplied them,
# and they are the only currency the rest of the module uses.
#
# The seed rule (contract section 3.1): `repeats > 1` randomises, so it
# requires a seed and never invents one. A single unrepeated replicate is
# dealt systematically -- units in first-appearance order, round robin into
# `v` folds within each stratum -- which touches no random stream and is
# therefore reproducible without a seed. Supplying a seed switches the deal to
# a permutation, and the same seed twice gives the same permutation.
.fb_sel_folds <- function(spec, data, seed = NULL) {
  if (!inherits(spec, "fb_resamples")) {
    stop("`spec` must be an <fb_resamples>.", call. = FALSE)
  }
  n <- nrow(data)

  .fb_sel_require_columns(data, spec$unit, "unit")
  .fb_sel_require_columns(data, spec$strata, "strata")

  unit_key <- if (length(spec$unit) > 0L) {
    .fb_sel_key(data, spec$unit)
  } else {
    as.character(seq_len(n))
  }
  units <- unique(unit_key)

  stratum <- .fb_sel_unit_strata(data, spec$strata, unit_key, units)

  if (identical(spec$scheme, "explicit")) {
    assignments <- list(.fb_sel_explicit_folds(
      spec$fold_id,
      unit_key,
      units,
      n
    ))
    assignment <- "explicit"
  } else if (is.null(seed)) {
    if (spec$repeats > 1L) {
      stop(.fb_refusal_condition(
        reason_code = "selection_seed_required",
        message = paste0(
          "fb_select(): `repeats = ",
          spec$repeats,
          "` randomises the fold assignment and no `seed` was given. No seed ",
          "is invented here: a default seed is an unreproducible assignment ",
          "that looks reproducible."
        )
      ))
    }
    assignments <- list(.fb_sel_deal(units, units, stratum, spec$v))
    assignment <- "systematic"
  } else {
    assignments <- vector("list", spec$repeats)
    .fb_sel_with_seed(seed, {
      for (r in seq_len(spec$repeats)) {
        assignments[[r]] <- .fb_sel_deal(
          sample(units),
          units,
          stratum,
          spec$v
        )
      }
    })
    assignment <- "random"
  }

  splits <- list()
  for (r in seq_along(assignments)) {
    fold_of_row <- unname(assignments[[r]][unit_key])
    for (f in .fb_sel_order_labels(unique(fold_of_row))) {
      assess_rows <- which(fold_of_row == f)
      train_rows <- which(fold_of_row != f)
      label <- if (length(assignments) > 1L) paste0(r, ".", f) else f
      if (length(train_rows) == 0L || length(assess_rows) == 0L) {
        stop(.fb_refusal_condition(
          reason_code = "selection_fold_not_usable",
          message = paste0(
            "fb_select(): fold ",
            label,
            " has ",
            length(train_rows),
            " training and ",
            length(assess_rows),
            " assessment rows. A fold with nothing to learn on, or nothing ",
            "to score on, is not a fold."
          )
        ))
      }
      splits[[length(splits) + 1L]] <- list(
        fold_id = f,
        repeat_id = r,
        label = label,
        train_rows = train_rows,
        assess_rows = assess_rows
      )
    }
  }

  structure(
    splits,
    class = "fb_sel_folds",
    unit = spec$unit,
    unit_key = unit_key,
    strata = spec$strata,
    scheme = spec$scheme,
    assignment = assignment,
    seed = seed,
    n_rows = n,
    n_units = length(units)
  )
}


# .fb_sel_require_columns() --- the named columns exist, or the refusal names
# the ones that do not.
.fb_sel_require_columns <- function(data, vars, arg) {
  if (length(vars) == 0L) {
    return(invisible(NULL))
  }
  absent <- setdiff(vars, names(data))
  if (length(absent) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_variable_not_in_data",
      message = paste0(
        "fb_select(): `",
        arg,
        "` names column(s) absent from the data: ",
        paste(absent, collapse = ", "),
        "."
      )
    ))
  }
  invisible(NULL)
}


# .fb_sel_unit_strata() --- the stratum of each unit. A stratifying variable
# that varies inside a unit puts the unit in two strata at once, which can
# only be resolved by splitting the unit, so it refuses instead.
.fb_sel_unit_strata <- function(data, strata_vars, unit_key, units) {
  if (length(strata_vars) == 0L) {
    out <- rep("all", length(units))
    names(out) <- units
    return(out)
  }
  skey <- .fb_sel_key(data, strata_vars)
  spread <- tapply(skey, unit_key, function(x) length(unique(x)))
  split_units <- names(spread)[spread > 1L]
  if (length(split_units) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_strata_not_constant_within_unit",
      message = paste0(
        "fb_select(): `",
        paste(strata_vars, collapse = ", "),
        "` is not constant within unit(s): ",
        paste(utils::head(split_units, 5L), collapse = ", "),
        if (length(split_units) > 5L) ", ..." else "",
        ". A unit that belongs to two strata cannot be dealt into a fold ",
        "without being split."
      )
    ))
  }
  out <- skey[match(units, unit_key)]
  names(out) <- units
  out
}


# .fb_sel_deal() --- round-robin units into `v` folds within each stratum.
#
# Dealing rather than cutting keeps the fold sizes within one unit of each
# other in every stratum, which is what makes the per-fold scores comparable.
# `order_units` fixes the order the deal walks in: the natural order for the
# systematic path, a permutation for the seeded one.
.fb_sel_deal <- function(order_units, units, stratum, v) {
  label <- character(length(units))
  names(label) <- units
  for (s in unique(stratum)) {
    in_stratum <- order_units[stratum[order_units] == s]
    if (length(in_stratum) == 0L) {
      next
    }
    label[in_stratum] <- as.character(
      rep(seq_len(v), length.out = length(in_stratum))
    )
  }
  label
}


# .fb_sel_explicit_folds() --- validate a user-supplied partition and lift it
# from rows to units.
#
# The load-bearing check is the last one. A `fold_id` that is not constant
# within a unit draws the fold boundary through the middle of the unit, which
# is exactly the leakage the `unit` argument exists to prevent, arriving by
# the back door.
.fb_sel_explicit_folds <- function(fold_id, unit_key, units, n) {
  if (length(fold_id) != n) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_id_invalid",
      message = paste0(
        "fb_select(): `fold_id` has ",
        length(fold_id),
        " entries and the data has ",
        n,
        " rows."
      )
    ))
  }
  if (anyNA(fold_id)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_id_invalid",
      message = paste0(
        "fb_select(): `fold_id` carries ",
        sum(is.na(fold_id)),
        " missing value(s). A row with no fold is neither trained on nor ",
        "scored, and dropping it silently would change the estimand."
      )
    ))
  }
  fid <- as.character(fold_id)
  if (length(unique(fid)) < 2L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_id_invalid",
      message = paste0(
        "fb_select(): `fold_id` describes ",
        length(unique(fid)),
        " fold, so there is no holdout."
      )
    ))
  }
  spread <- tapply(fid, unit_key, function(x) length(unique(x)))
  split_units <- names(spread)[spread > 1L]
  if (length(split_units) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_id_invalid",
      message = paste0(
        "fb_select(): `fold_id` splits unit(s) ",
        paste(utils::head(split_units, 5L), collapse = ", "),
        if (length(split_units) > 5L) ", ..." else "",
        " across folds. A unit on both sides of a fold boundary is the ",
        "leakage `unit` exists to prevent."
      )
    ))
  }
  out <- fid[match(units, unit_key)]
  names(out) <- units
  out
}


# .fb_sel_order_labels() --- fold labels in a readable order: numerically
# where they parse as numbers, so fold 10 does not sort before fold 2.
.fb_sel_order_labels <- function(x) {
  num <- suppressWarnings(as.numeric(x))
  if (!anyNA(num)) x[order(num)] else x[order(x, method = "radix")]
}


# --- frozen designs --------------------------------------------------- #

# .fb_sel_design() --- the model matrix over ALL rows, under the full data's
# factor levels and contrasts.
#
# Building it once and subsetting per fold is what freezes the levels and the
# contrast map across folds (contract section 3.3). It is not a leak: the
# design is a function of the recorded covariate values and the declared
# contrast scheme, and carries nothing about the response. Rebuilding it
# inside each fold would instead let a level absent from a training fold
# silently renumber every dummy column, so that the same coefficient means a
# different contrast in different folds.
.fb_sel_design <- function(formula, data) {
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  tt <- stats::terms(mf)
  X <- stats::model.matrix(tt, mf)
  list(
    X = X,
    terms = tt,
    frame = mf,
    xlev = stats::.getXlevels(tt, mf),
    assign = attr(X, "assign")
  )
}


# .fb_sel_scored_design() --- the design an arm's folds are scored on: the
# protected block beside the candidate block without its intercept, over all
# rows and before any fold-local preprocessing, as each fold's substrate
# builds it. The score contract hashes its values, so it is recorded once on
# the out-of-fold object.
.fb_sel_scored_design <- function(protected, candidates, data) {
  x_protected <- .fb_sel_design(protected, data)$X
  if (is.null(candidates)) {
    return(x_protected)
  }
  x_candidate <- .fb_sel_design(candidates, data)$X
  cbind(
    x_protected,
    x_candidate[, colnames(x_candidate) != "(Intercept)", drop = FALSE]
  )
}


# .fb_sel_check_fold_levels() --- an assessment fold carrying a level absent
# from its training rows is a named outcome, never a silent drop.
#
# The coefficient of an unseen level is not estimable from the training rows,
# so the fold's prediction for those rows is undefined. Dropping the rows
# would quietly change the denominator of the score; dropping the level would
# quietly change the model. The refusal names both the term and the level.
.fb_sel_check_fold_levels <- function(design, split) {
  frame <- design$frame
  for (nm in names(frame)) {
    column <- frame[[nm]]
    if (!is.factor(column) && !is.character(column)) {
      next
    }
    values <- as.character(column)
    absent <- setdiff(
      unique(values[split$assess_rows]),
      unique(values[split$train_rows])
    )
    if (length(absent) > 0L) {
      stop(.fb_refusal_condition(
        reason_code = "selection_fold_level_absent",
        message = paste0(
          "fb_select(): fold ",
          split$label,
          " carries level(s) of `",
          nm,
          "` absent from its training rows: ",
          paste(utils::head(absent, 5L), collapse = ", "),
          if (length(absent) > 5L) ", ..." else "",
          ". The level's contrast is not estimable on this fold."
        )
      ))
    }
  }
  invisible(NULL)
}


# --- preprocessing, learned on the training rows ---------------------- #

# .fb_sel_preprocessing() --- centre and scale, from the training rows only.
#
# A column with no spread among the training rows keeps a scale of one: it
# carries no information on this fold, and dividing by its zero spread would
# replace a harmless constant with NaN somewhere far from here.
.fb_sel_preprocessing <- function(X, train_rows) {
  if (ncol(X) == 0L) {
    return(list(centre = numeric(0), scale = numeric(0)))
  }
  block <- X[train_rows, , drop = FALSE]
  centre <- colMeans(block)
  spread <- apply(block, 2L, stats::sd)
  spread[!is.finite(spread) | spread < sqrt(.Machine$double.eps)] <- 1
  list(centre = centre, scale = spread)
}


# .fb_sel_apply_preprocessing() --- apply a learned centring and scaling.
.fb_sel_apply_preprocessing <- function(X, preprocessing) {
  if (ncol(X) == 0L) {
    return(X)
  }
  out <- sweep(X, 2L, preprocessing$centre, "-")
  sweep(out, 2L, preprocessing$scale, "/")
}


# --- the covariance a fold works in ----------------------------------- #

# .fb_sel_benchmark_random() --- a base's random terms in the spelling
# `.fb_sel_Z_from_vars()` reads: the group of a random intercept, and
# `outer:inner` for a two-way interaction.
#
# Every other term refuses here by name. Those two are the structures the
# covariance builder admits, and a term missing from its list is missing
# from V: a random slope or a multi-way interaction used to reach it as
# nothing, and a two-way interaction as two crossed intercepts, so the
# benchmark scored a covariance the base does not have.
.fb_sel_benchmark_random <- function(fb) {
  out <- character(0)
  refused <- character(0)
  for (t in fb$random_terms %||% list()) {
    type <- t$type %||% "unknown"
    if (type %in% c("simple", "ide", "id")) {
      out <- c(out, t$var)
    } else if (identical(type, "nested")) {
      out <- c(out, paste0(t$outer, ":", t$inner))
    } else {
      refused <- c(refused, type)
    }
  }
  if (length(refused) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_covariance_not_buildable",
      message = paste0(
        "fb_select(): the generalised-least-squares benchmark builds its ",
        "fold covariance from random intercepts and two-way interactions, ",
        "and this base also carries random terms of these types: ",
        paste(unique(refused), collapse = ", "),
        ". The Bayesian workflow refits the base's own structure on every ",
        "fold; the benchmark needs a base without those terms."
      )
    ))
  }
  unique(out)
}


# .fb_sel_Z_from_vars() --- the random-effect indicators for a data frame,
# where the grouping is named directly rather than read off a fitted base.
#
# It does not build them. `.fb_sel_random_designs()` is the one place in the
# module that turns a grouping into an indicator matrix, and a second
# implementation of that here is exactly the divergence this module exists to
# avoid, so this translates the `outer:inner` spelling into the random-term
# list that function reads and hands the work over.
.fb_sel_Z_from_vars <- function(data, vars) {
  if (length(vars) == 0L) {
    return(list())
  }
  needed <- unique(unlist(strsplit(vars, ":", fixed = TRUE)))
  .fb_sel_require_columns(data, needed, "random")
  random_terms <- lapply(vars, function(v) {
    parts <- strsplit(v, ":", fixed = TRUE)[[1L]]
    if (length(parts) == 1L) {
      list(type = "simple", var = parts[[1L]])
    } else if (length(parts) == 2L) {
      list(type = "nested", outer = parts[[1L]], inner = parts[[2L]])
    } else {
      stop(.fb_refusal_condition(
        reason_code = "selection_fold_variable_not_in_data",
        message = paste0(
          "fb_select(): `random` term `",
          v,
          "` names more than two groupings. The admitted spellings are ",
          "`group` and `outer:inner`."
        )
      ))
    }
  })
  .fb_sel_random_designs(list(
    data = data,
    fb = list(random_terms = random_terms)
  ))
}


# .fb_sel_reml_components() --- variance components from the TRAINING rows.
#
# Contract section 4.2 step 2 refits the base model per fold, and the
# production route reads that refit through `.fb_sel_variance_components()`
# under the D1 estimand. Where a full engine refit per fold is not what the
# caller wants to pay for, this is the independent REML plugin of the
# sensitivity set (section 4.5, estimand 3): the profiled restricted
# likelihood in the ratios phi_k = sigma_k^2 / sigma_e^2, with sigma_e^2 in
# closed form at the optimum. Either producer feeds the same
# `.fb_sel_cov_spec()`, and the estimand label travels with it.
#
# It is deterministic --- no random starts, no jitter --- so two calls on the
# same training rows agree bit for bit, which is what lets the leakage gate
# compare with `identical()` and need no tolerance at all.
.fb_sel_reml_components <- function(y, X, Z, train_rows, label = "1") {
  y_t <- as.numeric(y)[train_rows]
  x_t <- X[train_rows, , drop = FALSE]

  # Reduce the training design to an estimable basis first. A rank-deficient
  # X'H^-1X has no Cholesky factor and the criterion would return its
  # sentinel everywhere, which optim would read as a flat surface.
  decomposition <- qr(x_t)
  rank <- decomposition$rank
  x_t <- x_t[, sort(decomposition$pivot[seq_len(rank)]), drop = FALSE]

  n_t <- length(train_rows)
  if (n_t <= rank) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_not_usable",
      message = paste0(
        "fb_select(): fold ",
        label,
        " has ",
        n_t,
        " training rows for ",
        rank,
        " estimable columns, so it has no residual degrees of freedom and no ",
        "variance component can be estimated on it."
      )
    ))
  }

  z_t <- lapply(Z, function(z) {
    block <- z[train_rows, , drop = FALSE]
    block[, colSums(block != 0) > 0, drop = FALSE]
  })
  n_components <- length(z_t)

  profile <- function(phi) {
    h <- diag(1, nrow = n_t)
    for (k in seq_len(n_components)) {
      h <- h + phi[[k]] * tcrossprod(z_t[[k]])
    }
    factor_h <- tryCatch(chol(h), error = function(e) NULL)
    if (is.null(factor_h)) {
      return(NULL)
    }
    solve_h <- function(m) backsolve(factor_h, forwardsolve(t(factor_h), m))
    xtx <- crossprod(x_t, solve_h(x_t))
    factor_x <- tryCatch(chol(xtx), error = function(e) NULL)
    if (is.null(factor_x)) {
      return(NULL)
    }
    beta <- backsolve(
      factor_x,
      forwardsolve(t(factor_x), crossprod(x_t, solve_h(y_t)))
    )
    residual <- y_t - drop(x_t %*% beta)
    weighted_ss <- sum(residual * solve_h(residual))
    list(
      criterion = (n_t - rank) *
        log(weighted_ss / (n_t - rank)) +
        2 * sum(log(diag(factor_h))) +
        2 * sum(log(diag(factor_x))),
      sigma2 = weighted_ss / (n_t - rank)
    )
  }

  criterion <- function(log_phi) {
    fit <- profile(exp(log_phi))
    if (is.null(fit)) 1e10 else fit$criterion
  }

  if (n_components == 0L) {
    phi <- numeric(0)
    converged <- TRUE
  } else {
    optimum <- stats::optim(
      par = rep(log(0.5), n_components),
      fn = criterion,
      method = "L-BFGS-B",
      lower = rep(log(1e-8), n_components),
      upper = rep(log(1e8), n_components)
    )
    phi <- exp(optimum$par)
    converged <- identical(optimum$convergence, 0L)
  }

  fit <- profile(phi)
  terms <- phi * fit$sigma2
  names(terms) <- names(Z)

  list(
    residual = fit$sigma2,
    terms = terms,
    estimand = "reml_plugin_covariance",
    converged = converged,
    n_train = n_t,
    rank = rank
  )
}


# .fb_sel_cov_spec() --- the covariance a fold works in, as its parts rather
# than as a matrix.
#
# `Z` holds one indicator per random term over ALL rows, so a fold can take
# both its training block and its cross-block from the same object;
# `.fb_sel_Z_from_vars()` builds it from a data frame and
# `.fb_sel_random_designs()` from a fitted base. `residual` and `terms` are on
# the VARIANCE scale and `terms` is matched to `Z` by name.
.fb_sel_cov_spec <- function(
  Z,
  residual,
  terms,
  estimand = "reml_plugin_covariance"
) {
  if (length(Z) > 0L && !setequal(names(Z), names(terms))) {
    stop(
      "`Z` and `terms` must carry the same names: Z has ",
      paste(names(Z), collapse = ", "),
      " and terms has ",
      paste(names(terms), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.numeric(residual) || length(residual) != 1L || residual <= 0) {
    stop("`residual` must be a single positive variance.", call. = FALSE)
  }
  structure(
    list(Z = Z, residual = residual, terms = terms, estimand = estimand),
    class = "fb_sel_cov_spec"
  )
}


# --- oracle 7: the prediction estimand -------------------------------- #

# .fb_sel_fold_blocks() --- V_TT, its factor, and the assessment blocks.
#
# Contract section 5.4. With T training and S assessment rows,
#
#   E[y_S | y_T] = X_S beta + V_ST V_TT^-1 (y_T - X_T beta)
#   V_S_given_T  = V_SS - V_ST V_TT^-1 V_TS
#
# Only V_TT is factored, and V_ST V_TT^-1 is obtained as a solve against that
# factor. V^-1 is never formed: for a deposit-sized problem it cannot be, and
# for a small one it would still be the slower and less accurate route.
#
# `prediction = "new_cluster"` is the default because a held-out location's
# own effect is not estimable from the training rows. Scoring it conditionally
# hands the model the very random effect it is being asked to do without, and
# rewards it for having been handed it. The structural cross-covariance is
# computed either way and kept, so the estimand's effect is a number on the
# certificate rather than an argument about it.
.fb_sel_fold_blocks <- function(
  cov_spec,
  train_rows,
  assess_rows,
  prediction = c("new_cluster", "known_cluster"),
  label = "1"
) {
  prediction <- match.arg(prediction)
  n_train <- length(train_rows)
  n_assess <- length(assess_rows)

  v_tt <- diag(cov_spec$residual, nrow = n_train, ncol = n_train)
  v_ss <- diag(cov_spec$residual, nrow = n_assess, ncol = n_assess)
  v_st_structural <- matrix(0, nrow = n_assess, ncol = n_train)

  for (nm in names(cov_spec$Z)) {
    z_t <- cov_spec$Z[[nm]][train_rows, , drop = FALSE]
    z_s <- cov_spec$Z[[nm]][assess_rows, , drop = FALSE]
    variance <- cov_spec$terms[[nm]]
    v_tt <- v_tt + variance * tcrossprod(z_t)
    v_ss <- v_ss + variance * tcrossprod(z_s)
    v_st_structural <- v_st_structural + variance * tcrossprod(z_s, z_t)
  }

  # The indicator designs carry the data's row names, and a crossproduct
  # carries them into the covariance. `.fb_sel_V_dense()` strips them for the
  # same reason: a covariance block is indexed by position within the fold,
  # not by an original row label, and keeping the label invites a downstream
  # join on a name that means something else here.
  dimnames(v_tt) <- NULL
  dimnames(v_ss) <- NULL
  dimnames(v_st_structural) <- NULL

  v_st <- if (identical(prediction, "new_cluster")) {
    matrix(0, nrow = n_assess, ncol = n_train)
  } else {
    v_st_structural
  }

  factor_tt <- tryCatch(chol(v_tt), error = function(e) NULL)
  if (is.null(factor_tt)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_covariance_not_factorisable",
      message = paste0(
        "fb_select(): the training covariance of fold ",
        label,
        " is not positive definite, so it has no Cholesky factor. No jitter ",
        "was applied: a repair here is a misspecification carried into every ",
        "score this fold contributes."
      )
    ))
  }

  # A = V_ST V_TT^-1, as S solves against the training factor.
  gain <- t(backsolve(factor_tt, forwardsolve(t(factor_tt), t(v_st))))
  conditional <- v_ss - tcrossprod(gain, v_st)
  # Symmetric in exact arithmetic; the two triangles differ by rounding of
  # order 1e-16 and an unsymmetrised matrix makes chol() fail for no reason.
  conditional <- (conditional + t(conditional)) / 2

  list(
    v_tt = v_tt,
    factor_tt = factor_tt,
    v_ss = v_ss,
    v_st = v_st,
    v_st_structural = v_st_structural,
    gain = gain,
    v_s_given_t = conditional,
    prediction = prediction,
    logdet_v_tt = 2 * sum(log(diag(factor_tt)))
  )
}


# .fb_sel_gls_fit() --- generalised least squares as whitened least squares.
#
# This is the lambda = 0 corner of the module's objective (contract section
# 5.1) and, under D2, the corner that proves the penalised path. Aliased
# columns come back as zero rather than NA, with the rank recorded, so a
# rank-deficient candidate block yields an estimable-space answer instead of
# an NA that propagates silently into every prediction.
.fb_sel_gls_fit <- function(factor_tt, X_train, y_train) {
  x_white <- forwardsolve(t(factor_tt), X_train)
  y_white <- forwardsolve(t(factor_tt), as.numeric(y_train))
  if (ncol(X_train) == 0L) {
    return(list(
      beta = numeric(0),
      rank = 0L,
      residual_ss = sum(y_white^2)
    ))
  }
  decomposition <- qr(x_white)
  beta <- qr.coef(decomposition, y_white)
  aliased <- is.na(beta)
  beta[aliased] <- 0
  names(beta) <- colnames(X_train)
  list(
    beta = beta,
    rank = decomposition$rank,
    aliased = colnames(X_train)[aliased],
    residual_ss = sum((y_white - drop(x_white %*% beta))^2)
  )
}


# .fb_sel_fold_predict() --- the conditional mean of the assessment rows.
.fb_sel_fold_predict <- function(blocks, beta, X_train, y_train, X_assess) {
  fixed_assess <- if (length(beta) == 0L) {
    rep(0, nrow(X_assess))
  } else {
    drop(X_assess %*% beta)
  }
  fixed_train <- if (length(beta) == 0L) {
    rep(0, nrow(X_train))
  } else {
    drop(X_train %*% beta)
  }
  fixed_assess + drop(blocks$gain %*% (as.numeric(y_train) - fixed_train))
}


# --- oracle 4: the fold substrate, and the leak that must turn it red -- #

# .fb_sel_fold_substrate() --- everything one fold learns, and the assessment
# blocks it will be scored against.
#
# The object separates `learned` from `assessment` deliberately, and the
# separation is the contract rather than a convenience. Every element of
# `learned` is a quantity the fold claims to have taken from its training rows
# alone; `.fb_sel_leakage_gate()` holds it to that claim by rebuilding the
# fold with the assessment rows corrupted and asking whether any of them
# moved.
.fb_sel_fold_substrate <- function(
  data,
  split,
  response,
  protected,
  candidates = NULL,
  random = character(0),
  prediction = c("new_cluster", "known_cluster"),
  cov_spec = NULL,
  offsets = character(0)
) {
  prediction <- match.arg(prediction)
  # The response through the resolver the fit itself uses, so a transformed
  # response (`log(y)`, `I(10 * y)`) is the modelled quantity. Read as a
  # column name it was NULL, and the fold's variance components stopped
  # inside the optimiser.
  y <- .fb_response_values(list(response = response), data)
  if (is.null(y)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_fold_variable_not_in_data",
      message = paste0(
        "fb_select(): the response `", response, "` is neither a column of ",
        "the data nor an expression that evaluates over its rows."
      )
    ))
  }
  y <- as.numeric(y)
  # A Gaussian offset is subtracted here, before anything is whitened or
  # learned, and the fold works in y - o from then on. The components, the
  # pilot fit, the predictions and the scored residuals are all of y - o, and
  # a residual of y - o against a prediction of y - o is the residual of y
  # against the prediction plus the offset, so no add-back is needed to
  # score. It is evaluated on the data this call receives, so the leakage
  # gate's perturbed rebuilds see their own offsets.
  off <- .fb_offset_total(offsets, data)
  if (!is.null(off)) {
    y <- y - off
  }

  protected_design <- .fb_sel_design(protected, data)
  candidate_design <- if (is.null(candidates)) {
    list(X = matrix(numeric(0), nrow = nrow(data), ncol = 0L))
  } else {
    .fb_sel_design(candidates, data)
  }
  .fb_sel_check_fold_levels(protected_design, split)
  if (!is.null(candidates)) {
    .fb_sel_check_fold_levels(candidate_design, split)
  }

  # Drop the candidate block's intercept: the protected block owns it, and
  # two intercepts would be an aliased pair in every fold.
  x_candidate <- candidate_design$X
  keep <- colnames(x_candidate) != "(Intercept)"
  x_candidate <- x_candidate[, keep, drop = FALSE]

  # Learned on the training rows, and only there.
  preprocessing <- .fb_sel_preprocessing(x_candidate, split$train_rows)
  x_candidate <- .fb_sel_apply_preprocessing(x_candidate, preprocessing)

  Z <- .fb_sel_Z_from_vars(data, random)
  if (is.null(cov_spec)) {
    components <- .fb_sel_reml_components(
      y,
      protected_design$X,
      Z,
      split$train_rows,
      label = split$label
    )
    cov_spec <- .fb_sel_cov_spec(
      Z,
      components$residual,
      components$terms,
      components$estimand
    )
  }

  blocks <- .fb_sel_fold_blocks(
    cov_spec,
    split$train_rows,
    split$assess_rows,
    prediction = prediction,
    label = split$label
  )

  x_full <- cbind(protected_design$X, x_candidate)
  pilot <- .fb_sel_gls_fit(
    blocks$factor_tt,
    x_full[split$train_rows, , drop = FALSE],
    y[split$train_rows]
  )

  structure(
    list(
      fold_id = split$fold_id,
      repeat_id = split$repeat_id,
      label = split$label,
      train_rows = split$train_rows,
      assess_rows = split$assess_rows,
      prediction = prediction,
      covariance_estimand = cov_spec$estimand,
      response = y,
      x_protected = protected_design$X,
      x_candidate = x_candidate,
      blocks = blocks,
      # Every entry here is claimed to be a function of the training rows
      # alone. The leakage gate holds the claim to account.
      learned = list(
        components_residual = cov_spec$residual,
        components_terms = cov_spec$terms,
        preprocessing_centre = preprocessing$centre,
        preprocessing_scale = preprocessing$scale,
        pilot_beta = pilot$beta,
        logdet_v_tt = blocks$logdet_v_tt
      ),
      assessment = list(
        y = y[split$assess_rows],
        v_s_given_t = blocks$v_s_given_t
      )
    ),
    class = "fb_sel_substrate"
  )
}


# .fb_sel_perturb_rows() --- a deterministic, level-preserving corruption of
# the named rows.
#
# Numeric and logical columns move; factors and characters are left alone,
# because changing a level would change the design rather than the values and
# the gate would then report a difference that is not a leak. The map is
# affine and its coefficients are arbitrary; all that matters is that it moves
# every numeric value and moves nothing else.
.fb_sel_perturb_rows <- function(data, rows) {
  for (j in seq_along(data)) {
    column <- data[[j]]
    if (is.numeric(column)) {
      data[[j]][rows] <- -3 * column[rows] + 7
    } else if (is.logical(column)) {
      data[[j]][rows] <- !column[rows]
    }
  }
  data
}


# .fb_sel_leakage_gate() --- oracle 4, as a falsifier rather than as
# bookkeeping.
#
# `build_fn(data, split, ...)` must return an object carrying a `learned`
# list. The gate rebuilds the fold twice:
#
#   with the ASSESSMENT rows corrupted --- anything in `learned` that moves
#   was read from the assessment rows, and the fold is leaking;
#   with the TRAINING rows corrupted --- anything in `learned` that does NOT
#   move is inert, and the gate would have passed on it whatever the builder
#   did.
#
# The second half is why the gate cannot pass vacuously. A recorded provenance
# label, or a check that only ever asks the first question, both stay green on
# a field that is a hard-coded constant. Reporting inert fields turns that
# into a visible failure instead of an invisible one.
#
# Comparison is by `identical()` and carries no tolerance, because both builds
# are deterministic: a field untouched by the perturbation is bitwise equal,
# not nearly equal.
.fb_sel_leakage_gate <- function(build_fn, data, split, ...) {
  reference <- build_fn(data, split, ...)
  learned <- reference$learned
  if (!is.list(learned) || length(learned) == 0L || is.null(names(learned))) {
    stop(
      "the builder returned no named `learned` list, so the leakage gate ",
      "would pass on an empty claim.",
      call. = FALSE
    )
  }

  moved <- function(other) {
    names(learned)[
      !vapply(
        names(learned),
        function(nm) identical(learned[[nm]], other$learned[[nm]]),
        logical(1L)
      )
    ]
  }

  assessment_moved <- moved(
    build_fn(.fb_sel_perturb_rows(data, split$assess_rows), split, ...)
  )
  training_moved <- moved(
    build_fn(.fb_sel_perturb_rows(data, split$train_rows), split, ...)
  )

  list(
    clean = length(assessment_moved) == 0L &&
      length(training_moved) == length(learned),
    leaked = assessment_moved,
    inert = setdiff(names(learned), training_moved),
    checked = names(learned)
  )
}


# .fb_sel_certify_no_leakage() --- run the gate and refuse on a leak.
.fb_sel_certify_no_leakage <- function(build_fn, data, split, ...) {
  report <- .fb_sel_leakage_gate(build_fn, data, split, ...)
  if (length(report$leaked) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_leakage_detected",
      message = paste0(
        "fb_select(): fold ",
        split$label,
        " changed these quantities when only its assessment rows changed, so ",
        "they were not learned on its training rows: ",
        paste(report$leaked, collapse = ", "),
        "."
      )
    ))
  }
  if (length(report$inert) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_leakage_detected",
      message = paste0(
        "fb_select(): fold ",
        split$label,
        " reports these quantities as learned, but they did not move when ",
        "the training rows changed either, so the gate cannot see them: ",
        paste(report$inert, collapse = ", "),
        "."
      )
    ))
  }
  invisible(report)
}


# --- oracle 9: the scoring unit --------------------------------------- #

# .fb_sel_score_fold() --- the per-row and per-cluster records of one fold.
#
# Contract section 9. Two tables come back because the difference between them
# is the point. The row table is what a row-level scorer would aggregate, and
# it treats every observation as independent evidence. The cluster table is
# the reported one: the standard-error unit is the independent cluster, and a
# row-level standard error on correlated data understates it, so a
# one-standard-error rule computed at row level selects a different model and
# does so unreproducibly.
#
# The cluster negative log predictive density uses the JOINT block of
# V_S_given_T for the cluster; the row one uses the marginal variance of each
# row. They are equal only where the block is diagonal, and where they differ
# it is the joint one that is the density of what was observed.
.fb_sel_score_fold <- function(
  y_assess,
  mu_assess,
  v_s_given_t,
  cluster,
  label = "1",
  repeat_id = 1L,
  model = "full"
) {
  residual <- as.numeric(y_assess) - as.numeric(mu_assess)
  marginal_variance <- diag(v_s_given_t)
  rows <- data.table::data.table(
    model = model,
    fold = label,
    repeat_id = repeat_id,
    cluster = as.character(cluster),
    residual = residual,
    squared_error = residual^2,
    nlpd_row = 0.5 *
      (log(2 * pi) + log(marginal_variance) + residual^2 / marginal_variance)
  )

  groups <- split(
    seq_along(residual),
    factor(as.character(cluster), levels = unique(as.character(cluster)))
  )
  clusters <- data.table::rbindlist(lapply(names(groups), function(key) {
    index <- groups[[key]]
    block <- v_s_given_t[index, index, drop = FALSE]
    factor_block <- chol(block)
    standardised <- forwardsolve(t(factor_block), residual[index])
    data.table::data.table(
      model = model,
      fold = label,
      repeat_id = repeat_id,
      cluster = key,
      n = length(index),
      mse = mean(residual[index]^2),
      rmse = sqrt(mean(residual[index]^2)),
      nlpd = (0.5 *
        (length(index) *
          log(2 * pi) +
          2 * sum(log(diag(factor_block))) +
          sum(standardised^2))) /
        length(index)
    )
  }))

  list(rows = rows, clusters = clusters)
}


# .fb_sel_aggregate_scores() --- the reported score and its standard error, at
# the declared independent unit.
#
# Averaging within the unit first is what makes the standard error correct. The
# same cluster scored under three repeats is the same cluster measured three
# times, not three independent pieces of evidence, and treating the three fold
# scores as independent divides the standard error by a factor the data never
# earned.
.fb_sel_aggregate_scores <- function(
  clusters,
  metric = c("rmse", "mse", "nlpd"),
  unit = c("cluster", "resample")
) {
  metric <- match.arg(metric)
  unit <- match.arg(unit)
  key <- if (identical(unit, "cluster")) {
    as.character(clusters$cluster)
  } else {
    as.character(clusters$repeat_id)
  }
  per_unit <- vapply(
    split(clusters[[metric]], factor(key, levels = unique(key))),
    mean,
    numeric(1L)
  )
  if (length(per_unit) < 2L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_scoring_unit_undeclared",
      message = paste0(
        "fb_select(): the scoring unit `",
        unit,
        "` yields ",
        length(per_unit),
        " independent unit(s). A standard error needs at least two, and a ",
        "one-standard-error rule without one is not defined."
      )
    ))
  }
  list(
    metric = metric,
    unit = unit,
    score = mean(per_unit),
    se = stats::sd(per_unit) / sqrt(length(per_unit)),
    n_units = length(per_unit),
    per_unit = per_unit
  )
}


# .fb_sel_select_model() --- the decision rule over a scored path.
#
# `scores` carries one row per candidate model with `model`, `score`, `se` and
# `complexity`. The one-standard-error rule takes the standard error AT THE
# MINIMUM, which is the only place the rule is defined, and then the simplest
# model whose score is inside it. With the standard error computed at the row
# level that threshold is too tight and the rule collapses to the minimum,
# which is why section 9's unit is part of the contract rather than of the
# implementation.
.fb_sel_select_model <- function(
  scores,
  rule = c("one_se", "min", "tolerance"),
  tolerance = NULL
) {
  rule <- match.arg(rule)
  best <- which.min(scores$score)
  chosen <- switch(
    rule,
    min = best,
    one_se = {
      threshold <- scores$score[[best]] + scores$se[[best]]
      inside <- which(scores$score <= threshold)
      inside[which.min(scores$complexity[inside])]
    },
    tolerance = {
      if (!is.numeric(tolerance) || length(tolerance) != 1L) {
        stop(
          "`tolerance` must be a single number for rule = \"tolerance\".",
          call. = FALSE
        )
      }
      threshold <- scores$score[[best]] * (1 + tolerance)
      inside <- which(scores$score <= threshold)
      inside[which.min(scores$complexity[inside])]
    }
  )
  list(
    rule = rule,
    selected = scores$model[[chosen]],
    index = chosen,
    minimum = scores$model[[best]],
    threshold = switch(
      rule,
      min = scores$score[[best]],
      one_se = scores$score[[best]] + scores$se[[best]],
      tolerance = scores$score[[best]] * (1 + tolerance)
    )
  )
}


# --- the driver -------------------------------------------------------- #

# .fb_sel_cv() --- folds, substrate, prediction and scoring, in the order
# section 4.2 fixes.
#
# `candidate_sets` is a named list of integer column indices into the
# candidate block. The default is the pair the contract's `increments` field
# needs: the protected model alone, and the protected model with every
# candidate. A longer nested list makes a complexity path, which is what the
# one-standard-error rule selects along before any penalty code exists.
.fb_sel_cv <- function(
  data,
  response,
  protected,
  candidates = NULL,
  random = character(0),
  resamples = fb_resamples(),
  seed = NULL,
  prediction = c("new_cluster", "known_cluster"),
  candidate_sets = NULL,
  cluster = NULL,
  metric = c("rmse", "mse", "nlpd"),
  offsets = character(0)
) {
  prediction <- match.arg(prediction)
  metric <- match.arg(metric)

  # Step 1 of section 4.2, and it comes first for the reason the contract
  # gives: a rotation mixes observations, so a split of rotated rows is not
  # a holdout of anything.
  folds <- .fb_sel_folds(resamples, data, seed = seed)

  cluster_vars <- if (is.null(cluster)) resamples$unit else cluster
  .fb_sel_require_columns(data, cluster_vars, "cluster")
  cluster_key <- if (length(cluster_vars) > 0L) {
    .fb_sel_key(data, cluster_vars)
  } else {
    as.character(seq_len(nrow(data)))
  }

  # The candidate block's width is a property of the design, so the model
  # list is fixed once, before the folds, rather than rebuilt inside each of
  # them where a fold could quietly change what "the full model" means.
  n_candidates <- if (is.null(candidates)) {
    0L
  } else {
    sum(colnames(.fb_sel_design(candidates, data)$X) != "(Intercept)")
  }
  sets <- if (is.null(candidate_sets)) {
    list(protected_only = integer(0), full = seq_len(n_candidates))
  } else {
    candidate_sets
  }

  rows <- list()
  # Kept per model rather than as one table, because flexyBayes is not a
  # data.table-aware namespace: `dt[i, ]` inside it falls back to data.frame
  # semantics, so the module never subsets a data.table and collects instead.
  clusters <- list()
  substrates <- list()

  for (split in folds) {
    substrate <- .fb_sel_fold_substrate(
      data,
      split,
      response = response,
      protected = protected,
      candidates = candidates,
      random = random,
      prediction = prediction,
      offsets = offsets
    )
    substrates[[split$label]] <- substrate

    for (name in names(sets)) {
      columns <- sets[[name]]
      x <- cbind(
        substrate$x_protected,
        substrate$x_candidate[, columns, drop = FALSE]
      )
      fit <- .fb_sel_gls_fit(
        substrate$blocks$factor_tt,
        x[split$train_rows, , drop = FALSE],
        substrate$response[split$train_rows]
      )
      mu <- .fb_sel_fold_predict(
        substrate$blocks,
        fit$beta,
        x[split$train_rows, , drop = FALSE],
        substrate$response[split$train_rows],
        x[split$assess_rows, , drop = FALSE]
      )
      scored <- .fb_sel_score_fold(
        substrate$assessment$y,
        mu,
        substrate$assessment$v_s_given_t,
        cluster_key[split$assess_rows],
        label = split$label,
        repeat_id = split$repeat_id,
        model = name
      )
      rows[[length(rows) + 1L]] <- scored$rows
      clusters[[name]] <- c(clusters[[name]], list(scored$clusters))
    }
  }

  by_model <- lapply(clusters, data.table::rbindlist)
  aggregated <- data.table::rbindlist(lapply(names(by_model), function(name) {
    summary <- .fb_sel_aggregate_scores(by_model[[name]], metric = metric)
    data.table::data.table(
      model = name,
      score = summary$score,
      se = summary$se,
      n_units = summary$n_units,
      se_unit = summary$unit,
      complexity = length(sets[[name]])
    )
  }))

  list(
    folds = folds,
    substrates = substrates,
    # The design every fold was cut from, recorded once: the values the
    # score contract identifies.
    design = .fb_sel_scored_design(protected, candidates, data),
    rows = data.table::rbindlist(rows),
    clusters = data.table::rbindlist(by_model),
    by_model = by_model,
    aggregated = aggregated,
    metric = metric,
    prediction = prediction,
    seed = seed,
    # The four axes of SPEC v5 section 2, as far as this route resolves them.
    # This is the unpenalised BENCHMARK: a closed-form GLS on a covariance
    # estimated once per fold and plugged in. It is never labelled Bayesian.
    workflow = "unpenalized_benchmark",
    provider = "gls_reml_plugin",
    # Every quantity a fold learns -- the REML components, the candidate
    # standardisation and the pilot fit -- is read off its training rows,
    # which the leakage gate holds each substrate to.
    preprocessing = list(learned_on = "analysis_rows"),
    # Read from the fold objects, never from an argument. Each substrate
    # carries the estimand of the `fb_sel_cov_spec` it actually worked in,
    # which `.fb_sel_reml_components()` set; a certificate built from this
    # field cannot disagree with what ran, which is the defect SPEC v5
    # section 9.2 closes. `unique()` rather than the first fold's value, so a
    # run in which the folds disagreed reports two estimands rather than
    # silently reporting one of them.
    covariance_estimand = unique(vapply(
      substrates,
      function(z) z$covariance_estimand,
      character(1L)
    ))
  )
}
