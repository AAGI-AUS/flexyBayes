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
# Structures. V = sigma_e^2 I + sum_k sigma_k^2 Z_k G_k Z_k', one term (a
# slot) per variance:
#
#   (1 | g), ide(g)        Z the indicator of g, G = I
#   a:b                    Z the indicator of the a:b cell, G = I
#   (x || g)               the intercept slot, and Z = indicator(g) * x
#   diag(f):g              one slot per level l of f: the indicator of g on
#                          the rows at level l, zero elsewhere
#   ar1(row):ar1(col)      Z the indicator of the lattice node, one trial,
#                          G = R_row(rho_row) (x) R_col(rho_col)
#   at(t):ar1(r):ar1(c)    Z the indicator of the trial's lattice node,
#                          G = I_trial (x) R_row(rho_row) (x) R_col(rho_col)
#   vm(g, K), ped(g, A)    Z the indicator of g, G the known covariance; a
#                          precision carrier is inverted
#
# Anything else refuses with `selection_covariance_not_buildable`.
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


# --- the base's covariance, slot by slot --------------------------------- #

# .fb_sel_level_index() --- a grouping column as integer levels, built as the
# emitters build their indices: `factor()` of the column.
.fb_sel_level_index <- function(x) {
  f <- factor(x)
  list(index = as.integer(f), levels = levels(f))
}

# .fb_sel_cov_slot() --- one variance's structure.
#
# `index` is the column of Z each row loads on (NA where the slot has no
# effect), `value` the loading (NULL for one), `key` the row-level grouping
# within which the slot links rows, and `columns` the data columns that key
# is built from. `components` are the names the variance may carry in a
# fit's variance-component table: the INLA and brms spellings differ for a
# slope and for a known-covariance term.
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

# .fb_sel_cov_slots() --- the base's random and residual structure as slots,
# or a refusal naming every term the route does not build.
.fb_sel_cov_slots <- function(fb, data, known_matrices = list()) {
  residual <- fb$residual_terms %||% list()
  residual_types <- vapply(
    residual,
    function(t) t$type %||% "unknown",
    character(1L)
  )
  refused <- residual_types[residual_types != "units"]
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
      reason_code = "selection_covariance_not_buildable",
      message = paste0(
        "fb_select(): covariance = \"base_model\" builds the fold ",
        "covariance from random intercepts, two-way interactions, ",
        "uncorrelated random slopes, diag(f):g, an ar1(row):ar1(col) field ",
        "on one trial or id(trial):ar1(row):ar1(col) on each trial, and vm() ",
        "or ped() known covariances, over ",
        "an independent residual, and this base also carries: ",
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
      rho = c(row = "rho_row", col = "rho_col")
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


# --- the model, before any fold ------------------------------------------ #

# .fb_sel_cov_model() --- everything the base-model covariance needs that no
# fold changes: the slots, the engine, and how every fold will be factored.
# The dense bound is checked here, before any refit.
.fb_sel_cov_model <- function(base, engine = NULL, seed = NULL) {
  base <- .fb_sel_bayes_normalise_base(base)
  provider <- .fb_sel_bayes_provider(base)
  info <- base$extras$call_info %||% list()
  data <- base$data
  n <- nrow(data)
  slots <- .fb_sel_cov_slots(base$fb, data, info$known_matrices %||% list())
  blocking <- .fb_sel_cov_block_structure(slots, data)

  # At 48 trials of 72 plots the dense route's peak was 2.1 GB against 0.9 GB
  # by trial (review/units/U-AR_2026-09-24).
  field_per_trial <- isTRUE(blocking$applicable) && any(vapply(
    slots,
    function(s) !is.null(s$structure) && !is.null(s$key),
    logical(1L)
  ))
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
      group = blocking$group,
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
      # Each row's known variance, added to the residual on V's diagonal in
      # every fold, or NULL (R/known_variance.R).
      known_variance = .fb_sel_known_variance(base),
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
    fixed_alias = "drop"
  )
  out <- .fb_sel_cov_components(fit, slots)
  out$rows_passed <- nrow(fold_data)
  out
}

# .fb_sel_cov_components() --- the residual and every slot's variance at its
# posterior mean, and an autoregressive field's correlations at theirs.
.fb_sel_cov_components <- function(fit, slots) {
  vc <- .fb_sel_summary(fit)$varcomp
  read <- function(names, what) {
    hit <- which(vc$component %in% names)
    if (length(hit) != 1L) {
      stop(.fb_refusal_condition(
        reason_code = "selection_covariance_not_buildable",
        message = paste0(
          "fb_select(): the fold refit's variance-component table carries ",
          length(hit), " rows named ", paste(names, collapse = " or "),
          " for ", what, ", where one was expected. It carries: ",
          paste(vc$component, collapse = ", "), "."
        )
      ))
    }
    vc[hit, , drop = FALSE]
  }
  second_moment <- function(row) {
    if (!is.finite(row$estimate) || !is.finite(row$std.error)) {
      stop(.fb_refusal_condition(
        reason_code = "selection_covariance_not_buildable",
        message = paste0(
          "fb_select(): the fold refit reports no posterior standard ",
          "deviation for `", row$component, "`, so the posterior mean of its ",
          "variance, estimate^2 + std.error^2, cannot be formed."
        )
      ))
    }
    row$estimate^2 + row$std.error^2
  }

  residual <- second_moment(read("sigma", "the residual"))
  terms <- numeric(0)
  rho <- list()
  for (s in slots) {
    terms[[s$name]] <- second_moment(read(s$components, s$label))
    if (identical(s$structure$type, "ar1xar1")) {
      rho[[s$name]] <- c(
        row = read(s$structure$rho[["row"]], s$label)$estimate,
        col = read(s$structure$rho[["col"]], s$label)$estimate
      )
    }
  }
  list(residual = residual, terms = terms, rho = rho)
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
      known_variance = model$known_variance
    )
  } else {
    names_slots <- vapply(model$slots, function(s) s$name, character(1L))
    spec <- .fb_sel_cov_spec(
      list(),
      components$residual,
      components$terms[names_slots],
      estimand,
      known_variance = model$known_variance
    )
    spec$slots <- model$slots
    spec$rho <- components$rho
    spec$blocks <- model$whitening$blocks
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
# one block, from the slots.
#
# The same sums as the dense route, in the same order: the residual on the
# diagonal, then each slot's variance times its loadings' products. An
# indicator slot's entries equal the dense V's exactly, a per-trial field's
# to rounding.
.fb_sel_block_covariance <- function(cov_spec, rows_a, rows_b) {
  out <- if (identical(rows_a, rows_b)) {
    diag(.fb_sel_residual_diagonal(cov_spec, rows_a), nrow = length(rows_a))
  } else {
    matrix(0, nrow = length(rows_a), ncol = length(rows_b))
  }
  for (s in cov_spec$slots) {
    if (identical(s$structure$type, "ar1xar1")) {
      out <- out + cov_spec$terms[[s$name]] *
        .fb_sel_ar1_field_correlation(s, cov_spec$rho[[s$name]], rows_a, rows_b)
      next
    }
    same <- outer(s$index[rows_a], s$index[rows_b], "==")
    same[is.na(same)] <- FALSE
    loading <- if (is.null(s$value)) {
      same * 1
    } else {
      same * outer(s$value[rows_a], s$value[rows_b])
    }
    out <- out + cov_spec$terms[[s$name]] * loading
  }
  dimnames(out) <- NULL
  out
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
  weighted = FALSE
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
    builder = cov_model$builder,
    whitening = cov_model$whitening[
      c("mode", "group", "n_blocks", "largest_block", "bound_gb")
    ],
    equivalence_level = "map_exact_fixed_covariance",
    equivalence_reason = equivalence_reason,
    folds = folds
  )
}
