# The Bayesian valuation of a frozen support --- `fb_select(refit = )`.
#
# The penalised optimiser chooses a support on each fold's training rows. This
# file refits that support, frozen, as a Bayesian mixed model on the base's
# own engine, structure and priors, scores the fold's assessment rows by their
# log predictive density, and reads off the refit's posterior the share of
# between-trial variance a declared block of fixed effects explains. Nothing
# is selected again: each fold's model is its protected fixed part, the
# columns its support holds and the base's random terms.
#
# Three choices, and why.
#
#   Support columns, not terms. The optimiser's support is a set of design
#   columns, and a fold may keep one column of a factor and drop another
#   (known issue 19). The formula refits write terms and cannot say that. Each
#   support column is written into the fold's data as a numeric column of its
#   own, taken from the candidate design learned on the fold's training rows,
#   so the refit holds exactly the frozen support and a fold is never refused
#   or dropped for a partial factor.
#
#   The unseen cluster in closed form. Under `new_cluster` every random term's
#   level on an assessment row is one no training row carries, so given a
#   posterior draw its effect is a Gaussian with the term's variance and mean
#   zero. Where every random term is an independent intercept on a Gaussian
#   likelihood, the draw's predictive for the row is N(f_d, s2_d + sum_k
#   s2_kd), with f_d the fixed part and offset: the unseen effects are
#   integrated exactly and only the posterior of the parameters is averaged
#   over draws. Elsewhere (a slope, a structured term, `known_cluster`) the
#   predictive is the engine's own draws of the linear predictor.
#
#   The block share. Gelman, Goodrich, Gabry and Vehtari (2019) define the
#   Bayesian R-squared of a draw as the variance of the fitted values over
#   the variance of the fitted values plus the residual variance. Here the
#   "fitted values" are the declared block's contribution to the linear
#   predictor, averaged within each level of the between-trial grouping, and
#   the "residual variance" is that grouping's random-intercept variance:
#   what the block explains of the between-trial variance, over what it
#   explains plus what is left.

.FB_SEL_VALUE_FIELDS <- c("method", "block", "group", "n_draws")

.FB_SEL_VALUE_FORMULA <- paste0(
  "elpd = sum over the U held-out units u of elpd_u, where elpd_u is the sum ",
  "over the unit's assessment rows i of log((1/D) sum_d p(y_i | theta_d)), ",
  "theta_d the D posterior draws of the refit of the fold that held u out, ",
  "averaged over repeats; se = sqrt(U) * sd(elpd_u)"
)


# --- the request ------------------------------------------------------------ #

#' Validate `fb_select(refit = )` against the base, before anything is fitted
#'
#' @param refit `NULL` or the named list the caller passed.
#' @param base The normalised base fit.
#' @param candidates The candidate formula.
#' @return `NULL`, or `list(method, n_draws, block_terms, block_label,
#'   group, group_reason)`.
#' @noRd
#' @keywords internal
.fb_sel_value_request <- function(refit, base, candidates) {
  if (is.null(refit)) {
    return(NULL)
  }
  if (!is.list(refit) || is.null(names(refit)) ||
        any(!nzchar(names(refit)))) {
    stop(
      "`refit` must be NULL or a named list; list(method = \"bayes\") ",
      "refits with the defaults.",
      call. = FALSE
    )
  }
  unknown <- setdiff(names(refit), .FB_SEL_VALUE_FIELDS)
  if (length(unknown) > 0L) {
    stop(
      "`refit` carries ", paste0("`", unknown, "`", collapse = ", "),
      ", which it does not read. Its fields are ",
      paste0("`", .FB_SEL_VALUE_FIELDS, "`", collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!identical(refit$method, "bayes")) {
    stop(
      "`refit$method` must be \"bayes\": each fold's frozen support is ",
      "refitted as a Bayesian model on the base's own engine.",
      call. = FALSE
    )
  }
  n_draws <- refit$n_draws %||% 4000L
  if (!is.numeric(n_draws) || length(n_draws) != 1L || !is.finite(n_draws) ||
        n_draws < 1 || n_draws != round(n_draws)) {
    stop(
      "`refit$n_draws` must be one whole number of posterior draws, at ",
      "least 1.",
      call. = FALSE
    )
  }

  data <- base$data
  protected_labels <- .fb_sel_protected_labels(base$fb)
  candidate_labels <- attr(stats::terms(candidates, data = data), "term.labels")
  block <- refit$block
  if (is.null(block)) {
    block_terms <- candidate_labels
    block_label <- "every candidate"
  } else {
    if (!inherits(block, "formula") || length(block) != 2L) {
      stop(
        "`refit$block` must be a one-sided formula over terms of the base's ",
        "fixed part or of `candidates`.",
        call. = FALSE
      )
    }
    block_terms <- attr(stats::terms(block, data = data), "term.labels")
    absent <- setdiff(block_terms, c(protected_labels, candidate_labels))
    if (length(block_terms) == 0L || length(absent) > 0L) {
      stop(
        "`refit$block` names ",
        if (length(absent) > 0L) {
          paste0(paste0("`", absent, "`", collapse = ", "), ", which ")
        } else {
          "no term, and a block "
        },
        "must name terms of the base's fixed part (",
        paste(protected_labels, collapse = ", "), ") or of `candidates` (",
        paste(candidate_labels, collapse = ", "), ").",
        call. = FALSE
      )
    }
    block_label <- paste(block_terms, collapse = " + ")
  }

  intercepts <- unique(unlist(lapply(
    base$fb$random_terms %||% list(),
    function(t) if (identical(t$type, "simple")) t$var else NULL
  )))
  group <- refit$group
  group_reason <- NA_character_
  if (!is.null(group)) {
    if (!is.character(group) || length(group) != 1L ||
          !group %in% intercepts) {
      stop(
        "`refit$group` must name the grouping column of a random intercept ",
        "in the base, whose variance is the between-trial variance the ",
        "block share is stated against. The base's random intercepts are ",
        "grouped by: ",
        if (length(intercepts) > 0L) {
          paste(intercepts, collapse = ", ")
        } else {
          "none"
        },
        ".",
        call. = FALSE
      )
    }
  } else if (length(intercepts) == 1L) {
    group <- intercepts
  } else {
    group <- NA_character_
    group_reason <- if (length(intercepts) == 0L) {
      paste0(
        "the base has no random intercept, so there is no between-trial ",
        "variance to state the share against"
      )
    } else {
      paste0(
        "the base has random intercepts on ",
        paste(intercepts, collapse = ", "),
        "; name the one whose variance is the between-trial variance with ",
        "refit = list(group = )"
      )
    }
  }

  list(
    method = "bayes",
    n_draws = as.integer(n_draws),
    block_terms = block_terms,
    block_label = block_label,
    group = group,
    group_reason = group_reason
  )
}


#' Refuse a refit on a route that freezes no support per fold
#'
#' @param request The validated request, or `NULL`.
#' @param workflow The resolved workflow.
#' @return `NULL`, invisibly, or a `selection_refit_support_undefined`
#'   refusal.
#' @noRd
#' @keywords internal
.fb_sel_value_require_route <- function(request, workflow) {
  if (is.null(request) || identical(workflow, "penalized_optimizer")) {
    return(invisible(NULL))
  }
  stop(.fb_refusal_condition(
    reason_code = "selection_refit_support_undefined",
    message = paste0(
      "fb_select(): refit = list(method = \"bayes\") refits the support the ",
      "penalised optimiser froze on each fold, and workflow = \"", workflow,
      "\" freezes none: ",
      if (identical(workflow, "bayes")) {
        paste0(
          "the Bayesian workflow scores every candidate on every fold, or ",
          "decides on the whole sample"
        )
      } else {
        "the unpenalised benchmark holds every candidate"
      },
      ". Nothing was fitted. Use workflow = \"penalized_optimizer\" with a ",
      "penalty that selects."
    ),
    workflow = workflow
  ))
}


# --- the frozen supports ---------------------------------------------------- #

#' Each fold's support at the optimiser's chosen point
#'
#' @param arm The optimiser arm (`benchmark$arms$glmnet`).
#' @param folds The folds the arm ran on.
#' @return A named list, one character vector of design columns per fold
#'   label, or a `selection_refit_support_undefined` refusal where the
#'   family selects nothing.
#' @noRd
#' @keywords internal
.fb_sel_value_supports <- function(arm, folds) {
  if (!isTRUE(arm$support$selects) || !isTRUE(arm$selection$ran)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_refit_support_undefined",
      message = paste0(
        "fb_select(): refit = list(method = \"bayes\") refits the support ",
        "the penalised optimiser froze on each fold, and a ", arm$family,
        " path shrinks every coefficient without setting one to zero, so ",
        "it has no support to freeze. Nothing was refitted. Use laplace(), ",
        "elastic_net() or a grouped family, which set coefficients to zero."
      ),
      penalty_family = arm$family
    ))
  }
  path <- as.data.frame(arm$support$by_lambda)
  chosen <- path[
    path$arm == "benchmark" &
      path$lambda_index == arm$selection$selected_lambda_index, ,
    drop = FALSE
  ]
  labels <- vapply(folds, function(s) as.character(s$label), character(1L))
  stats::setNames(lapply(labels, function(label) {
    held <- strsplit(
      as.character(chosen$terms[as.character(chosen$fold) == label]),
      ", ",
      fixed = TRUE
    )
    held <- unlist(held)
    held[nzchar(held)]
  }), labels)
}


# --- the value contract ----------------------------------------------------- #

#' The score contract two valuations must share to be differenced unit by unit
#'
#' The fold of every row, the scoring unit, the prediction target and the
#' response values are identified, through the same encoding and hash as
#' every other score contract. The design is deliberately not: a paired
#' comparison varies the predictor set, so every value contract carries the
#' same one-column intercept design and the hash does not depend on the
#' candidates or the support.
#'
#' @return An `fb_select_score_contract` at stage `"value"`.
#' @noRd
#' @keywords internal
.fb_sel_value_contract <- function(folds, data, response, offsets, prediction) {
  .fb_select_score_contract(
    folds = folds,
    data = data,
    response = response,
    prediction = prediction,
    metric = "elpd",
    design = matrix(1, nrow(data), 1L, dimnames = list(NULL, "(Intercept)")),
    learned_on = "analysis_rows",
    stage = "value",
    scoring_unit = NULL,
    offsets = offsets,
    candidate_space = NULL
  )
}


# --- the predictive draws --------------------------------------------- #

#' Score assessment rows from per-draw Gaussian predictive moments
#'
#' @param y The assessment rows' observed response.
#' @param mean Draws in rows, assessment rows in columns.
#' @param variance The same shape as `mean`, or one variance per draw.
#' @return `list(lpd, mean, sd)`, one entry per assessment row: the log of
#'   the draw-averaged density, and the mean and standard deviation of the
#'   mixture over draws.
#' @noRd
#' @keywords internal
.fb_sel_value_score <- function(y, mean, variance) {
  n_draws <- nrow(mean)
  if (!is.matrix(variance)) {
    variance <- matrix(variance, n_draws, length(y))
  }
  logdens <- stats::dnorm(
    matrix(y, nrow = n_draws, ncol = length(y), byrow = TRUE),
    mean = mean,
    sd = sqrt(variance),
    log = TRUE
  )
  centre <- colMeans(mean)
  centred <- mean - matrix(centre, n_draws, length(y), byrow = TRUE)
  spread <- colMeans(centred^2)
  list(
    lpd = apply(logdens, 2L, .fb_sel_logmeanexp),
    mean = centre,
    sd = sqrt(colMeans(variance) + spread)
  )
}


#' Score assessment rows a block of rows at a time
#'
#' Each block's predictive moments are built, scored and dropped before the
#' next, so no draws-by-rows matrix larger than one block is held.
#'
#' @param y The assessment rows' observed response.
#' @param moments `function(columns)` returning `list(mean, variance)` for
#'   those assessment rows, in the shapes `.fb_sel_value_score()` takes.
#' @param size Rows per block.
#' @return As `.fb_sel_value_score()`.
#' @noRd
#' @keywords internal
.fb_sel_value_score_blocks <- function(y, moments, size = 1000L) {
  n <- length(y)
  parts <- lapply(seq.int(1L, n, by = size), function(first) {
    columns <- seq.int(first, min(n, first + size - 1L))
    m <- moments(columns)
    .fb_sel_value_score(y[columns], m$mean, m$variance)
  })
  list(
    lpd = unlist(lapply(parts, function(p) p$lpd), use.names = FALSE),
    mean = unlist(lapply(parts, function(p) p$mean), use.names = FALSE),
    sd = unlist(lapply(parts, function(p) p$sd), use.names = FALSE)
  )
}


#' The independent random intercepts of an INLA refit whose levels on every
#' assessment row are unseen in training
#'
#' @return `list(ok, names, reason)`: the `f()` index names, or why the
#'   closed form does not apply.
#' @noRd
#' @keywords internal
.fb_sel_value_inla_unseen <- function(fit, split) {
  calls <- list()
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1L]], as.name("f"))) {
        calls[[length(calls) + 1L]] <<- e
      } else {
        for (part in as.list(e)[-1L]) walk(part)
      }
    }
  }
  walk(fit$extras$formula[[3L]])
  data_inla <- fit$inla$.args$data
  names_out <- character(0)
  for (e in calls) {
    args <- as.list(e)[-1L]
    arg_names <- names(args) %||% rep("", length(args))
    positional <- args[!nzchar(arg_names)]
    named <- arg_names[nzchar(arg_names)]
    shown <- paste(deparse(e[[2L]]), collapse = "")
    if (length(positional) != 1L || !is.name(positional[[1L]]) ||
          !all(named %in% c("model", "hyper")) ||
          !identical(args$model, "iid")) {
      return(list(ok = FALSE, reason = paste0(
        "the random term f(", shown, ", ...) is not an independent ",
        "intercept, so its unseen effect is read off the draws"
      )))
    }
    name <- as.character(positional[[1L]])
    index <- data_inla[[name]]
    n_engine <- length(split$train_rows) + length(split$assess_rows) +
      length(split$not_estimable_rows)
    if (is.null(index) || length(index) != n_engine) {
      return(list(ok = FALSE, reason = paste0(
        "the index of f(", name, ") is not one value per engine row"
      )))
    }
    seen <- intersect(
      stats::na.omit(index[split$assess_rows]),
      stats::na.omit(index[split$train_rows])
    )
    if (length(seen) > 0L || anyNA(index[split$assess_rows])) {
      return(list(ok = FALSE, reason = paste0(
        "some assessment rows carry a level of f(", name, ") that training ",
        "rows carry, so the effect is conditioned on and read off the draws"
      )))
    }
    names_out <- c(names_out, name)
  }
  list(ok = TRUE, names = names_out, reason = NA_character_)
}


#' The block share of each posterior draw
#'
#' @param design The block's design columns on the fold's training rows.
#' @param key The between-trial grouping on those rows.
#' @param beta The block's coefficients, one column per draw.
#' @param sigma2 The grouping's random-intercept variance per draw.
#' @return One share per draw, or `NULL` with fewer than two levels.
#' @noRd
#' @keywords internal
.fb_sel_value_share_draws <- function(design, key, beta, sigma2) {
  key <- as.character(key)
  if (length(unique(key)) < 2L) {
    return(NULL)
  }
  if (ncol(design) == 0L) {
    return(rep(0, length(sigma2)))
  }
  sums <- rowsum(design, key)
  means <- sums / as.numeric(table(key)[rownames(sums)])
  contribution <- means %*% beta
  numerator <- apply(contribution, 2L, stats::var)
  numerator / (numerator + sigma2)
}


#' One fold's predictive moments and block-share draws on INLA
#'
#' @noRd
#' @keywords internal
.fb_sel_value_fold_inla <- function(fit, inputs, split, rhs, offsets,
                                    prediction, n_draws, seed, block, group) {
  engine_data <- inputs$engine_data
  X <- .fb_sel_design(rhs, engine_data)$X
  fixed_names <- rownames(fit$inla$summary.fixed)
  names_match <- setequal(colnames(X), fixed_names)
  unseen <- if (identical(prediction, "new_cluster")) {
    .fb_sel_value_inla_unseen(fit, split)
  } else {
    list(ok = FALSE, reason = paste0(
      "prediction = \"known_cluster\" conditions on the training rows' ",
      "levels, which the engine's draws carry"
    ))
  }
  precisions <- paste0("Precision for ", unseen$names)
  named <- all(precisions %in% rownames(fit$inla$summary.hyperpar))
  analytic <- names_match && isTRUE(unseen$ok) && named
  reason <- if (analytic) {
    NA_character_
  } else if (!names_match) {
    "the engine's fixed-effect names are not the design's"
  } else if (!isTRUE(unseen$ok)) {
    unseen$reason
  } else {
    "a random term's precision is not named after its index"
  }

  sample_seed <- as.integer(abs(seed %||% 1L) %% .Machine$integer.max)
  if (sample_seed == 0L) {
    sample_seed <- 1L
  }
  draws <- if (analytic) {
    .fb_inla_posterior_sample(
      fit, as.integer(n_draws), seed = sample_seed, num.threads = "1:1",
      selection = stats::setNames(
        as.list(rep(1L, length(fixed_names))), fixed_names
      )
    )
  } else {
    .fb_inla_posterior_sample(
      fit, as.integer(n_draws), seed = sample_seed, num.threads = "1:1"
    )
  }
  hyper <- matrix(
    unlist(lapply(draws, function(z) z$hyperpar), use.names = FALSE),
    ncol = length(draws[[1L]]$hyperpar),
    byrow = TRUE,
    dimnames = list(NULL, names(draws[[1L]]$hyperpar))
  )
  sigma2 <- 1 / hyper[, "Precision for the Gaussian observations"]
  latent_names <- rownames(draws[[1L]]$latent)
  latent <- function(rows) {
    index <- match(rows, latent_names)
    vapply(draws, function(z) z$latent[index, 1L], numeric(length(index)))
  }

  assess <- split$assess_rows
  y <- inputs$assessment$y
  if (analytic) {
    beta <- matrix(latent(paste0(colnames(X), ":1")), nrow = ncol(X))
    offset <- .fb_offset_total(offsets, engine_data)
    offset <- if (is.null(offset)) numeric(length(assess)) else offset[assess]
    variance <- sigma2 + rowSums(1 / hyper[, precisions, drop = FALSE])
    scored <- .fb_sel_value_score_blocks(y, function(columns) {
      rows <- assess[columns]
      mean <- t(X[rows, , drop = FALSE] %*% beta) +
        matrix(offset[columns], nrow(hyper), length(columns), byrow = TRUE)
      list(mean = mean, variance = variance)
    })
  } else {
    predictor <- grep("^Predictor:", latent_names)
    if (length(predictor) != nrow(engine_data)) {
      stop(
        "the sampled latent field carries ", length(predictor),
        " Predictor entries for ", nrow(engine_data),
        " engine rows, so the assessment rows cannot be indexed into it.",
        call. = FALSE
      )
    }
    scored <- .fb_sel_value_score_blocks(y, function(columns) {
      index <- predictor[assess[columns]]
      mean <- matrix(
        vapply(draws, function(z) z$latent[index, 1L], numeric(length(index))),
        nrow = length(index)
      )
      list(mean = t(mean), variance = sigma2)
    })
  }

  share <- NULL
  share_reason <- NA_character_
  precision_group <- paste0("Precision for ", group)
  if (!is.na(group) && !precision_group %in% colnames(hyper)) {
    share_reason <- paste0(
      "the refit carries no precision named `", precision_group, "`"
    )
  } else if (!is.na(group)) {
    engine_names <- block$engine_names$inla
    missing <- setdiff(engine_names, sub(":1$", "", latent_names))
    if (length(missing) > 0L) {
      share_reason <- paste0(
        "the refit's draws carry no coefficient named ",
        paste(missing, collapse = ", ")
      )
    } else {
      beta_block <- if (length(engine_names) == 0L) {
        matrix(numeric(0), nrow = 0L, ncol = length(draws))
      } else {
        matrix(
          latent(paste0(engine_names, ":1")), nrow = length(engine_names)
        )
      }
      share <- .fb_sel_value_share_draws(
        block$design, block$key, beta_block, 1 / hyper[, precision_group]
      )
    }
  }

  list(
    scored = scored,
    method = if (analytic) {
      "analytic_gaussian"
    } else {
      "posterior_predictive_draws"
    },
    reason = reason,
    n_draws = length(draws),
    share = share,
    share_reason = share_reason
  )
}


#' One fold's predictive moments and block-share draws on brms
#'
#' @noRd
#' @keywords internal
.fb_sel_value_fold_brms <- function(fit, inputs, prediction, n_draws, seed,
                                    block, group) {
  newdata <- inputs$assess_data
  total <- brms::ndraws(fit$brms)
  ids <- if (total <= n_draws) {
    seq_len(total)
  } else {
    unique(as.integer(round(seq(1, total, length.out = n_draws))))
  }
  draws <- as.matrix(fit$brms)[ids, , drop = FALSE]
  sigma2 <- draws[, "sigma"]^2
  sd_names <- grep("^sd_", colnames(draws), value = TRUE)
  groups <- sub("^sd_(.*?)__.*$", "\\1", sd_names, perl = TRUE)
  coefs <- sub("^sd_.*?__", "", sd_names, perl = TRUE)
  numeric_coef <- vapply(coefs, function(c) {
    identical(c, "Intercept") || is.numeric(newdata[[c]])
  }, logical(1L))
  analytic <- identical(prediction, "new_cluster") && all(numeric_coef)
  reason <- if (analytic) {
    NA_character_
  } else if (!identical(prediction, "new_cluster")) {
    paste0(
      "prediction = \"known_cluster\" conditions on the training rows' ",
      "levels, which the engine's draws carry"
    )
  } else {
    paste0(
      "the group-level coefficients ",
      paste(unique(coefs[!numeric_coef]), collapse = ", "),
      " are not numeric columns, so the unseen effects are read off the draws"
    )
  }

  y <- inputs$assessment$y
  if (analytic) {
    scored <- .fb_sel_value_score_blocks(y, function(columns) {
      part <- newdata[columns, , drop = FALSE]
      mean <- .fb_sel_with_seed(seed %||% 1L, brms::posterior_epred(
        fit$brms, newdata = part, re_formula = NA, draw_ids = ids
      ))
      # Each group's unseen effect on a row: sum over its coefficients of
      # z^2 sd^2, plus 2 z_a z_b cor sd_a sd_b over each correlated pair.
      variance <- matrix(sigma2, nrow(mean), ncol(mean))
      for (g in unique(groups)) {
        in_group <- which(groups == g)
        z <- lapply(coefs[in_group], function(c) {
          if (identical(c, "Intercept")) rep(1, nrow(part)) else part[[c]]
        })
        s <- draws[, sd_names[in_group], drop = FALSE]
        for (a in seq_along(in_group)) {
          variance <- variance + outer(s[, a]^2, z[[a]]^2)
          for (b in seq_len(a - 1L)) {
            cor_name <- paste0(
              "cor_", g, "__", coefs[in_group][[b]], "__", coefs[in_group][[a]]
            )
            if (cor_name %in% colnames(draws)) {
              variance <- variance + 2 * outer(
                draws[, cor_name] * s[, a] * s[, b], z[[a]] * z[[b]]
              )
            }
          }
        }
      }
      list(mean = mean, variance = variance)
    })
  } else {
    predictive <- .fb_sel_predict_fold_brms(
      fit, seq_len(nrow(newdata)), newdata,
      prediction = prediction, n_draws = n_draws, seed = seed
    )
    scored <- .fb_sel_value_score_blocks(y, function(columns) {
      list(
        mean = predictive$mu[, columns, drop = FALSE],
        variance = predictive$sigma^2
      )
    })
  }

  share <- NULL
  share_reason <- NA_character_
  sd_group <- paste0("sd_", group, "__Intercept")
  if (!is.na(group) && !sd_group %in% colnames(draws)) {
    share_reason <- paste0("the refit carries no draw named `", sd_group, "`")
  } else if (!is.na(group)) {
    # brms names population-level coefficients in design order; the block's
    # columns are found by position among the non-intercept ones.
    b_names <- setdiff(
      grep("^b_", colnames(draws), value = TRUE),
      "b_Intercept"
    )
    position <- match(block$engine_names$brms, block$design_columns)
    if (length(b_names) != length(block$design_columns)) {
      share_reason <- paste0(
        "the refit carries ", length(b_names), " population-level ",
        "coefficients against ", length(block$design_columns),
        " design columns"
      )
    } else {
      beta_block <- t(draws[, b_names[position], drop = FALSE])
      share <- .fb_sel_value_share_draws(
        block$design, block$key, beta_block, draws[, sd_group]^2
      )
    }
  }

  list(
    scored = scored,
    method = if (analytic) {
      "analytic_gaussian"
    } else {
      "posterior_predictive_draws"
    },
    reason = reason,
    n_draws = length(ids),
    share = share,
    share_reason = share_reason
  )
}


# --- the driver ------------------------------------------------------------- #

#' Refit each fold's frozen support and score its assessment rows
#'
#' @param base The normalised base fit.
#' @param candidates The candidate formula.
#' @param resamples The resampling specification.
#' @param prediction The prediction target.
#' @param folds The folds the optimiser ran on.
#' @param arm The optimiser arm, or `NULL` where `supports` is given.
#' @param request The validated `refit` request.
#' @param engine The brms sampler settings, or `NULL`.
#' @param seed The run's seed.
#' @param budget Seconds, or `NULL`.
#' @param supports `NULL` to read each fold's support off `arm` at its chosen
#'   point, or a named list of design columns per fold label, which the
#'   grouped ablation passes for a path without one block.
#' @param candidate_space The optimiser's candidate space. Under
#'   `protected_alias = "drop"` each fold's refit holds the protected
#'   columns that fold kept.
#' @return The valuation record stored at `fb_select()$value`.
#' @noRd
#' @keywords internal
.fb_sel_value <- function(base,
                          candidates,
                          resamples,
                          prediction,
                          folds,
                          arm = NULL,
                          request,
                          engine = NULL,
                          seed = NULL,
                          budget = NULL,
                          supports = NULL,
                          candidate_space = NULL) {
  supports <- supports %||% .fb_sel_value_supports(arm, folds)
  drop_protected <- identical(candidate_space$protected_alias, "drop")
  provider <- .fb_sel_bayes_provider(base)
  spec <- .fb_sel_bayes_base_spec(base)
  data <- spec$data
  .fb_sel_require_relabelable(base$fb, prediction)
  .fb_sel_require_columns(
    data,
    unique(c(
      spec$response,
      all.vars(spec$protected),
      spec$random_vars,
      .fb_offset_vars(spec$offsets)
    )),
    "data"
  )
  # The protected part is written as formula text, which an INLA refit
  # evaluates over training and assessment rows together; the support
  # columns are written as values learned on the training rows.
  .fb_sel_require_fold_local_formula(list(spec$protected), data, folds)
  cost <- .fb_sel_cost_estimate(
    base,
    resamples,
    n_folds = length(folds),
    fits_per_fold = 1L,
    probe_seconds = .fb_select_base_seconds(base)
  )
  .fb_sel_bayes_check_budget(cost, budget)

  unit_key <- attr(folds, "unit_key") %||% as.character(seq_len(nrow(data)))
  contract <- .fb_sel_value_contract(
    folds, data, spec$response, spec$offsets, prediction
  )
  .fb_select_verify_contract(contract, "value")

  protected_design <- .fb_sel_design(spec$protected, data)
  protected_labels_all <- attr(protected_design$terms, "term.labels")
  protected_block <- colnames(protected_design$X)[
    protected_design$assign > 0L &
      protected_labels_all[pmax(protected_design$assign, 1L)] %in%
        request$block_terms
  ]
  existing_names <- names(data)
  group <- request$group
  base_seed <- as.integer(seed %||% 0L)

  fold_job <- function(f) {
    split <- folds[[f]]
    label <- as.character(split$label)
    held <- supports[[label]]
    candidate_design <- .fb_sel_design(candidates, data, split$train_rows)
    candidate_terms <- attr(candidate_design$terms, "term.labels")[
      candidate_design$assign[match(held, colnames(candidate_design$X))]
    ]
    column_names <- utils::tail(
      make.unique(c(existing_names, paste0("fb_support_", make.names(held)))),
      length(held)
    )
    fold_data <- data
    if (!is.null(split$weights)) {
      fold_data[[.FB_SEL_WEIGHTS_COL]] <- split$weights
    }
    for (j in seq_along(held)) {
      fold_data[[column_names[[j]]]] <- as.numeric(
        candidate_design$X[, held[[j]]]
      )
    }
    protected_labels <- spec$protected_labels
    protected_formula <- spec$protected
    protected_names <- stats::setNames(protected_block, protected_block)
    if (drop_protected) {
      # The protected columns this fold kept, written as numeric columns the
      # way the support is, so the refit leaves out the dropped ones.
      dropped <- candidate_space$folds[[label]]$protected_dropped$column
      kept <- setdiff(colnames(protected_design$X), c("(Intercept)", dropped))
      kept_names <- utils::tail(
        make.unique(c(
          names(fold_data),
          paste0("fb_protected_", make.names(kept))
        )),
        length(kept)
      )
      for (j in seq_along(kept)) {
        fold_data[[kept_names[[j]]]] <- as.numeric(
          protected_design$X[, kept[[j]]]
        )
      }
      protected_labels <- kept_names
      protected_formula <- stats::as.formula(paste0(
        "~ ", paste(c("1", kept_names), collapse = " + ")
      ))
      protected_names <- stats::setNames(
        kept_names[match(protected_block, kept)],
        protected_block
      )
      protected_names <- protected_names[!is.na(protected_names)]
      # Held out with the scored rows, so the refit does not train on them.
      fold_data[[spec$response]][split$not_estimable_rows] <- NA
    }
    inputs <- .fb_sel_bayes_fold_inputs(
      fold_data,
      split,
      response = spec$response,
      protected = protected_formula,
      candidates = NULL,
      random = spec$random_vars,
      prediction = prediction
    )
    fixed <- .fb_sel_bayes_fixed(
      spec$response, protected_labels, column_names, spec$offsets
    )
    rhs <- stats::as.formula(paste0(
      "~ ",
      paste(c("1", protected_labels, column_names), collapse = " + ")
    ))
    fit_seed <- base_seed + 1000L * f + 1L
    t_fit <- proc.time()
    fit <- .fb_sel_bayes_refit(
      base,
      if (identical(provider, "inla")) {
        inputs$engine_data
      } else {
        inputs$train_data
      },
      fixed,
      provider,
      engine = engine,
      seed = fit_seed,
      na_action = if (identical(provider, "inla")) "augment" else NULL
    )

    # The declared block on this fold: its protected columns and the support
    # columns whose term it names, on the training rows.
    in_block <- candidate_terms %in% request$block_terms
    design_columns <- colnames(.fb_sel_design(rhs, inputs$engine_data)$X)
    design_columns <- design_columns[design_columns != "(Intercept)"]
    block <- list(
      design = cbind(
        protected_design$X[
          split$train_rows, names(protected_names), drop = FALSE
        ],
        candidate_design$X[split$train_rows, held[in_block], drop = FALSE]
      ),
      key = if (is.na(group)) NULL else data[[group]][split$train_rows],
      engine_names = list(
        inla = c(unname(protected_names), column_names[in_block]),
        brms = c(unname(protected_names), column_names[in_block])
      ),
      design_columns = design_columns
    )

    out <- if (identical(provider, "inla")) {
      .fb_sel_value_fold_inla(
        fit, inputs, split, rhs, spec$offsets, prediction,
        request$n_draws, fit_seed, block, group
      )
    } else {
      .fb_sel_value_fold_brms(
        fit, inputs, prediction, request$n_draws, fit_seed, block, group
      )
    }
    seconds <- unname((proc.time() - t_fit)["elapsed"])
    rows <- data.frame(
      repeat_id = as.integer(split$repeat_id),
      fold = label,
      row = as.integer(split$assess_rows),
      row_key = paste0("repeat:", split$repeat_id, "/row:", split$assess_rows),
      unit = as.character(unit_key[split$assess_rows]),
      observed = as.numeric(inputs$assessment$y),
      predictive_mean = out$scored$mean,
      predictive_sd = out$scored$sd,
      lpd = out$scored$lpd,
      method = out$method,
      prediction = prediction,
      n_draws = as.integer(out$n_draws),
      stringsAsFactors = FALSE
    )
    list(
      label = label,
      repeat_id = as.integer(split$repeat_id),
      rows = rows,
      held = held,
      method = out$method,
      reason = out$reason,
      share = out$share,
      share_reason = out$share_reason,
      seconds = seconds
    )
  }

  plan <- .fb_sel_worker_plan(length(folds), provider = provider)
  t_start <- proc.time()
  results <- .fb_sel_fold_lapply(
    seq_along(folds),
    fold_job,
    plan,
    warm_first = identical(provider, "brms")
  )
  elapsed <- unname((proc.time() - t_start)["elapsed"])

  rows <- do.call(rbind, lapply(results, function(r) r$rows))
  rows <- rows[
    order(rows$repeat_id, rows$row, method = "radix"), ,
    drop = FALSE
  ]
  rows$contract_hash <- contract$hash
  rownames(rows) <- NULL
  units <- .fb_sel_value_units(rows)

  supports_table <- do.call(rbind, lapply(results, function(r) {
    data.frame(
      repeat_id = rep(r$repeat_id, length(r$held)),
      fold = rep(r$label, length(r$held)),
      column = r$held,
      stringsAsFactors = FALSE
    )
  }))
  rownames(supports_table) <- NULL

  list(
    method = request$method,
    engine = provider,
    prediction = prediction,
    n_draws = request$n_draws,
    seed = seed,
    optimizer = if (is.null(arm)) {
      NULL
    } else {
      list(
        provider = arm$provider,
        family = arm$family,
        selected_lambda_index = arm$selection$selected_lambda_index,
        selected_s = arm$selection$selected_s
      )
    },
    supports = supports_table,
    folds = data.frame(
      repeat_id = vapply(results, function(r) r$repeat_id, integer(1L)),
      fold = vapply(results, function(r) r$label, character(1L)),
      method = vapply(results, function(r) r$method, character(1L)),
      reason = vapply(results, function(r) r$reason, character(1L)),
      seconds = vapply(results, function(r) r$seconds, numeric(1L)),
      stringsAsFactors = FALSE
    ),
    rows = rows,
    units = units,
    # Under protected_alias = "drop", the held-out rows of each fold with no
    # estimable prediction, which are not scored.
    not_estimable = if (drop_protected) {
      data.frame(
        fold = vapply(folds, function(f) as.character(f$label), character(1L)),
        n_not_estimable = vapply(
          folds,
          function(f) length(f$not_estimable_rows),
          integer(1L)
        ),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    },
    summary = .fb_sel_value_summary(units, contract),
    block_share = .fb_sel_value_block_share(results, request),
    contract = contract,
    parallel = plan,
    timing = list(
      estimate_seconds = cost$seconds,
      elapsed_seconds = elapsed
    )
  )
}


# --- the tables and their summaries ----------------------------------- #

#' One row per repeat and held-out unit: its rows' summed log density
#'
#' @noRd
#' @keywords internal
.fb_sel_value_units <- function(rows) {
  key <- paste(rows$repeat_id, rows$unit, sep = "\r")
  first <- !duplicated(key)
  units <- data.frame(
    repeat_id = rows$repeat_id[first],
    fold = rows$fold[first],
    unit = rows$unit[first],
    stringsAsFactors = FALSE
  )
  index <- match(key, key[first])
  units$n_rows <- as.integer(tabulate(index, nbins = nrow(units)))
  units$elpd <- as.numeric(vapply(
    split(rows$lpd, factor(index, levels = seq_len(nrow(units)))),
    sum,
    numeric(1L)
  ))
  units$contract_hash <- rows$contract_hash[first]
  units <- units[
    order(units$repeat_id, units$unit, method = "radix"), ,
    drop = FALSE
  ]
  rownames(units) <- NULL
  units
}


#' The summed log density over held-out units and its standard error
#'
#' @noRd
#' @keywords internal
.fb_sel_value_summary <- function(units, contract, value = "elpd") {
  per_unit <- vapply(
    split(units[[value]], factor(units$unit, levels = unique(units$unit))),
    mean,
    numeric(1L)
  )
  unit <- contract$scoring_unit
  list(
    elpd = sum(per_unit),
    se = if (length(per_unit) > 1L) {
      sqrt(length(per_unit)) * stats::sd(per_unit)
    } else {
      NA_real_
    },
    n_units = length(per_unit),
    unit = paste(unit, collapse = ":"),
    formula = .FB_SEL_VALUE_FORMULA
  )
}


#' The block share: per fold and pooled over folds
#'
#' @noRd
#' @keywords internal
.fb_sel_value_block_share <- function(results, request) {
  numerator <- paste0(
    "V_t(m_t): the sample variance, with divisor T - 1, over the T levels t ",
    "of `", request$group, "` on the fold's training rows, of m_t, the mean ",
    "over level t's training rows of sum_j x_ij beta_j over the block's ",
    "columns j in the fold's model"
  )
  denominator <- paste0(
    "V_t(m_t) + sigma^2, with sigma^2 the variance of the `", request$group,
    "` random intercept"
  )
  empty <- data.frame(
    block = character(0), median = numeric(0), lower = numeric(0),
    upper = numeric(0), n_draws = integer(0), stringsAsFactors = FALSE
  )
  reasons <- unique(stats::na.omit(vapply(
    results, function(r) r$share_reason, character(1L)
  )))
  if (is.na(request$group) || length(reasons) > 0L ||
        any(vapply(results, function(r) is.null(r$share), logical(1L)))) {
    return(list(
      available = FALSE,
      reason = if (is.na(request$group)) {
        request$group_reason
      } else if (length(reasons) > 0L) {
        paste(reasons, collapse = "; ")
      } else {
        "a fold's training rows carry fewer than two levels of the grouping"
      },
      block = request$block_label,
      block_terms = request$block_terms,
      group = request$group,
      numerator = NA_character_,
      denominator = NA_character_,
      by_fold = cbind(fold = character(0), empty),
      pooled = empty
    ))
  }
  summarise <- function(draws) {
    q <- stats::quantile(draws, c(0.025, 0.5, 0.975), names = FALSE, type = 7)
    data.frame(
      block = request$block_label,
      median = q[[2L]],
      lower = q[[1L]],
      upper = q[[3L]],
      n_draws = length(draws),
      stringsAsFactors = FALSE
    )
  }
  by_fold <- do.call(rbind, lapply(results, function(r) {
    cbind(fold = r$label, summarise(r$share), stringsAsFactors = FALSE)
  }))
  rownames(by_fold) <- NULL
  list(
    available = TRUE,
    reason = NA_character_,
    block = request$block_label,
    block_terms = request$block_terms,
    group = request$group,
    numerator = numerator,
    denominator = denominator,
    by_fold = by_fold,
    pooled = summarise(unlist(lapply(results, function(r) r$share)))
  )
}


# --- two valuations, unit by unit ------------------------------------- #

#' Pair two valuations on one value contract
#'
#' @param x,against Two `fb_select` results.
#' @return `list(units, difference, se, n_units)`, or a
#'   `selection_value_contract_mismatch` refusal.
#' @noRd
#' @keywords internal
.fb_sel_value_pair <- function(x, against) {
  if (!inherits(against, "fb_select")) {
    stop("`against` must be an <fb_select> result.", call. = FALSE)
  }
  a <- x$value
  b <- against$value
  if (is.null(a) || is.null(b)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_value_contract_mismatch",
      message = paste0(
        "fb_select(): ",
        if (is.null(a)) "`x`" else "`against`",
        " carries no valuation, so there are no units to pair. Run it with ",
        "refit = list(method = \"bayes\")."
      ),
      fields = "value"
    ))
  }
  .fb_sel_value_pair_records(a, b)
}


#' Pair two valuation records on one value contract
#'
#' @param a,b Two valuation records (`fb_select()$value`).
#' @return As `.fb_sel_value_pair()`.
#' @noRd
#' @keywords internal
.fb_sel_value_pair_records <- function(a, b) {
  .fb_select_verify_contract(a$contract, "value")
  .fb_select_verify_contract(b$contract, "value")
  if (!identical(a$contract$hash, b$contract$hash)) {
    fields <- .fb_select_contract_differences(a$contract, b$contract)
    stop(.fb_refusal_condition(
      reason_code = "selection_value_contract_mismatch",
      message = paste0(
        "fb_select(): the two valuations do not hold out the same rows in ",
        "the same folds under the same target, so their units do not pair. ",
        "Their value contracts differ in: ",
        paste(fields, collapse = ", "),
        ". Run both with the same data, `resamples`, `seed` and `prediction`."
      ),
      fields = fields
    ))
  }
  units <- a$units[, c("repeat_id", "fold", "unit", "n_rows", "elpd")]
  units$elpd_against <- b$units$elpd
  units$difference <- units$elpd - units$elpd_against
  summary <- .fb_sel_value_summary(units, a$contract, value = "difference")
  list(
    units = units,
    difference = summary$elpd,
    se = summary$se,
    n_units = summary$n_units
  )
}


# --- the accessors' empty shapes -------------------------------------- #

.fb_sel_value_empty_rows <- function() {
  data.frame(
    repeat_id = integer(0), fold = character(0), row = integer(0),
    row_key = character(0), unit = character(0), observed = numeric(0),
    predictive_mean = numeric(0), predictive_sd = numeric(0),
    lpd = numeric(0), method = character(0), prediction = character(0),
    n_draws = integer(0), contract_hash = character(0),
    stringsAsFactors = FALSE
  )
}

.fb_sel_value_empty_units <- function() {
  data.frame(
    repeat_id = integer(0), fold = character(0), unit = character(0),
    n_rows = integer(0), elpd = numeric(0), contract_hash = character(0),
    stringsAsFactors = FALSE
  )
}
