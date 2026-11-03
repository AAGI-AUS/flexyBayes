# fb_kfold.R -- exact K-fold cross-validation of an ordinal fit.
#
# Each fold is held out once. The model is refitted without the fold's
# responses, and each held-out group is scored by the log of its posterior
# predictive probability, its rows taken jointly. The groups' scores sum to
# the expected log predictive density (elpd).

#' Exact K-fold cross-validation of an ordinal fit
#'
#' @description
#' Refits the model once per fold with that fold's responses held out, and
#' scores each held-out group by its log predictive density, the group's
#' rows taken jointly. With `group` naming the plot, a plot is scored by the
#' probability of all of its plants' categories together, which is how
#' ordinal models of sub-sampled field trials are compared when
#' leave-one-out approximations are unreliable.
#'
#' @details
#' On INLA a fold's responses are set to missing and its rows kept, so every
#' latent node, an AR1 field's included, stays in the model, and the held-out
#' linear predictors are drawn from the refit's joint posterior. On brms the
#' model is refitted to the other folds' rows, and a held-out group's random
#' effect at a level the refit did not see is drawn from its fitted variance
#' (brms's `allow_new_levels = TRUE` with `sample_new_levels = "gaussian"`).
#'
#' The score of group `g` is
#' `log(mean over draws s of prod over its rows i of P(y_i | draw s))`, the
#' `elpd` is the sum over groups, and its standard error is
#' `sqrt(G) * sd(score)` over the `G` groups. Rows with a missing response
#' are not scored.
#'
#' Every category must keep a training observation in every fold, or the
#' refit would have fewer thresholds than the model it is scored against.
#'
#' @param fit An ordinal fit, from `flexybayes()` with
#'   `family = "cumulative"` or `"cratio"`.
#' @param folds The fold of each row of the data the fit was given: a vector
#'   with one label per row and at least two distinct labels.
#' @param group The name of a column of the fit's data whose levels are
#'   scored jointly, for example the plot, or `NULL` (the default) to score
#'   each row on its own. The rows of a group must share a fold.
#' @param n_draws The number of posterior draws per fold. Default `1000`.
#' @param seed An integer seed, or `NULL` (the default) to use the session's
#'   random number stream. Fold `f` uses `seed + f`.
#' @param ... Further arguments of [fb_refit()] for each fold's refit, for
#'   example `control` on brms.
#'
#' @returns An object of class `fb_kfold`, a list with `estimates`, a data
#'   frame holding `elpd` and its standard error `se`, `pointwise`, one row
#'   per scored group with `group`, `fold`, `n_rows` and `elpd`, and `folds`,
#'   `engine` and `group`.
#'
#' @references Cao Z, Megirian R, Nguyen M, Sparks AH (2026). Bayesian
#'   ordinal analysis of cereal growth stage data from oat trials in Western
#'   Australia. arXiv:2610.37384.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(
#'     trt = factor(rep(c("a", "b"), each = 60)),
#'     plot = factor(rep(1:24, each = 5))
#'   )
#'   eta <- 0.8 * (d$trt == "b") + rnorm(24, 0, 0.5)[d$plot]
#'   d$score <- cut(eta + rlogis(120), c(-Inf, -0.5, 0.5, 1.5, Inf),
#'     labels = c("low", "mid", "high", "top"), ordered_result = TRUE)
#'   fit <- flexybayes(score ~ trt, random = ~plot, data = d,
#'     family = "cumulative", backend = "inla", verbose = FALSE)
#'   folds <- rep(rep(1:4, 6), each = 5)
#'   fb_kfold(fit, folds, group = "plot", seed = 1L)
#' }
#' }
#' @export
fb_kfold <- function(
  fit,
  folds,
  group = NULL,
  n_draws = 1000L,
  seed = NULL,
  ...
) {
  .fb_kfold_check_fit(fit)
  n_draws <- .fb_ordinal_check_n_draws(n_draws)
  ord <- fit$extras$ordinal
  ir <- .fb_fit_ir(fit)
  data <- .fb_fit_data(fit)
  groups <- .fb_kfold_groups(data, folds, group)
  if (!ir$response %in% names(data)) {
    stop(.fb_refusal_condition(
      reason_code = "kfold_response_not_a_column",
      message = paste0(
        "fb_kfold() holds a fold out by setting its responses to missing, ",
        "and the response `",
        ir$response,
        "` is an expression, not a column of the data. Add it to the data ",
        "as a column and fit the model to that column."
      )
    ))
  }
  codes <- .fb_ordinal_codes(ir, data, ord)

  # One refit per fold -----------------------------------------------------
  fold_ids <- unique(folds)
  pointwise <- vector("list", length(fold_ids))
  for (f in seq_along(fold_ids)) {
    held <- folds == fold_ids[[f]]
    .fb_kfold_check_training(codes[!held], ord, fold_ids[[f]])
    scored <- which(held & !is.na(codes))
    if (!length(scored)) {
      next
    }
    fold_seed <- if (is.null(seed)) NULL else as.integer(seed) + f
    loglik <- if (identical(ord$engine, "brms")) {
      .fb_kfold_loglik_brms(fit, data, held, scored, codes, n_draws,
        fold_seed, ...)
    } else {
      .fb_kfold_loglik_inla(fit, data, held, scored, codes, n_draws,
        fold_seed, ...)
    }
    by_group <- split(seq_along(scored), groups[scored])
    pointwise[[f]] <- data.frame(
      group = names(by_group),
      fold = rep(fold_ids[[f]], length(by_group)),
      n_rows = lengths(by_group, use.names = FALSE),
      elpd = vapply(
        by_group,
        function(i) .fb_log_mean_exp(rowSums(loglik[, i, drop = FALSE])),
        numeric(1L)
      ),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  # Sum over groups ---------------------------------------------------------
  pointwise <- do.call(rbind, pointwise)
  structure(
    list(
      estimates = data.frame(
        elpd = sum(pointwise$elpd),
        se = sqrt(nrow(pointwise)) * stats::sd(pointwise$elpd)
      ),
      pointwise = pointwise,
      folds = length(fold_ids),
      engine = ord$engine,
      group = group
    ),
    class = "fb_kfold"
  )
}

#' @export
print.fb_kfold <- function(x, ...) {
  cat(
    "Exact ",
    x$folds,
    "-fold cross-validation (",
    x$engine,
    "), ",
    nrow(x$pointwise),
    if (is.null(x$group)) {
      " rows scored"
    } else {
      paste0(" groups of `", x$group, "` scored jointly")
    },
    "\n",
    sep = ""
  )
  cat(
    "  elpd = ",
    format(round(x$estimates$elpd, 2L), nsmall = 2L),
    " (SE ",
    format(round(x$estimates$se, 2L), nsmall = 2L),
    ")\n",
    sep = ""
  )
  invisible(x)
}


# Checks ------------------------------------------------------------------------

#' Refuse a fit fb_kfold() does not score
#'
#' @param fit The fit.
#' @returns Invisibly `TRUE`.
#' @noRd
#' @keywords internal
.fb_kfold_check_fit <- function(fit) {
  if (!inherits(fit, "flexybayes_ordinal")) {
    stop(.fb_refusal_condition(
      reason_code = "kfold_family_not_supported",
      message = paste0(
        "fb_kfold() scores fits with family = \"cumulative\" or \"cratio\" ",
        "in this version. For a brms fit of another family, ",
        "brms::kfold(fit$brms) cross-validates it."
      )
    ))
  }
  invisible(TRUE)
}

#' The scoring group of each row, after checking folds and groups
#'
#' @param data The fit's data.
#' @param folds The fold labels.
#' @param group The group column name or `NULL`.
#' @returns A character vector, one group label per row.
#' @noRd
#' @keywords internal
.fb_kfold_groups <- function(data, folds, group) {
  refuse <- function(why) {
    stop(.fb_refusal_condition(
      reason_code = "kfold_folds_invalid",
      message = paste0("fb_kfold(): ", why)
    ))
  }
  if (length(folds) != nrow(data) || anyNA(folds)) {
    refuse(paste0(
      "`folds` must hold one fold label for each of the ",
      nrow(data),
      " rows of the fit's data, with none missing."
    ))
  }
  if (length(unique(folds)) < 2L) {
    refuse("`folds` must hold at least two distinct folds.")
  }
  if (is.null(group)) {
    return(as.character(seq_len(nrow(data))))
  }
  if (!is.character(group) || length(group) != 1L || !group %in% names(data)) {
    refuse("`group` must name one column of the fit's data.")
  }
  groups <- as.character(data[[group]])
  spread <- tapply(folds, groups, function(f) length(unique(f)))
  if (any(spread > 1L)) {
    refuse(paste0(
      "the rows of a group must share a fold, and ",
      sum(spread > 1L),
      " group(s) of `",
      group,
      "` span more than one (",
      paste(utils::head(names(spread)[spread > 1L], 3L), collapse = ", "),
      ")."
    ))
  }
  groups
}

#' Refuse a fold that leaves a category without a training observation
#'
#' @param train_codes The category codes of the training rows.
#' @param ord The ordinal record.
#' @param fold The fold label.
#' @returns Invisibly `TRUE`.
#' @noRd
#' @keywords internal
.fb_kfold_check_training <- function(train_codes, ord, fold) {
  absent <- setdiff(seq_len(ord$K), train_codes)
  if (length(absent)) {
    stop(.fb_refusal_condition(
      reason_code = "kfold_category_missing_in_training",
      message = paste0(
        "fb_kfold(): fold ",
        fold,
        " holds out every observation of ",
        if (length(absent) == 1L) "category " else "categories ",
        paste(ord$levels[absent], collapse = ", "),
        ", so its refit would have fewer thresholds than the model. ",
        "Assign folds that keep each category in training, or merge rare ",
        "categories."
      )
    ))
  }
  invisible(TRUE)
}


# Held-out log likelihood -------------------------------------------------------

#' The log of a mean of exponentials, computed stably
#'
#' @param x A numeric vector.
#' @returns `log(mean(exp(x)))`.
#' @noRd
#' @keywords internal
.fb_log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

#' Refit one fold quietly
#'
#' The refit runs with `verbose = FALSE` (and `mcmc_verbose = FALSE` on
#' brms) unless `...` says otherwise. INLA prints a note on standard output
#' for every fit with a missing response, about the link it uses for those
#' rows' fitted values. The scores read the linear predictor, so the
#' refit's standard output is not shown. Messages and warnings pass.
#'
#' @param fit The ordinal fit.
#' @param data The fold's data.
#' @param readout The fb_refit() readout.
#' @param ... Passed to fb_refit().
#' @returns The refit.
#' @noRd
#' @keywords internal
.fb_kfold_refit <- function(fit, data, readout, ...) {
  dots <- list(...)
  if (is.null(dots$verbose)) {
    dots$verbose <- FALSE
  }
  if (is.null(dots$mcmc_verbose) && inherits(fit, "flexybayes_brms")) {
    dots$mcmc_verbose <- FALSE
  }
  # A fold's held-out rows are missing responses by construction, so the
  # high-missingness warning says nothing about the refit.
  old <- options(flexyBayes.silence_high_missingness_warning = TRUE)
  on.exit(options(old), add = TRUE)
  refit <- NULL
  utils::capture.output(
    refit <- do.call(
      fb_refit,
      c(list(fit, data = data, readout = readout), dots)
    )
  )
  refit
}

#' Held-out log likelihood of an INLA fold
#'
#' @param fit The ordinal fit.
#' @param data The fit's data.
#' @param held Logical, the fold's rows.
#' @param scored The held-out rows with an observed response.
#' @param codes The category codes of every row.
#' @param n_draws,seed Draws and seed.
#' @param ... Passed to fb_refit().
#' @returns A matrix, draws by `scored`.
#' @noRd
#' @keywords internal
.fb_kfold_loglik_inla <- function(
  fit,
  data,
  held,
  scored,
  codes,
  n_draws,
  seed,
  ...
) {
  ord <- fit$extras$ordinal
  response <- .fb_fit_ir(fit)$response
  fold_data <- data
  fold_data[[response]][held] <- NA
  refit <- .fb_kfold_refit(fit, fold_data, "draws", ...)
  dr <- .fb_ordinal_draws_inla(refit, n_draws, seed, rows = scored)
  .fb_kfold_category_loglik(dr$eta, dr$tau, codes[scored], ord)
}

#' Held-out log likelihood of a brms fold
#'
#' @inheritParams .fb_kfold_loglik_inla
#' @returns A matrix, draws by `scored`.
#' @noRd
#' @keywords internal
.fb_kfold_loglik_brms <- function(
  fit,
  data,
  held,
  scored,
  codes,
  n_draws,
  seed,
  ...
) {
  refit <- .fb_kfold_refit(fit, data[!held, , drop = FALSE], "full", ...)
  newdata <- data[scored, , drop = FALSE]
  newdata[[.FB_ORDINAL_RESPONSE]] <- codes[scored]
  .fb_ordinal_seeded(seed, brms::log_lik(
    refit$brms,
    newdata = newdata,
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian",
    ndraws = min(n_draws, brms::ndraws(refit$brms))
  ))
}

#' Log probability of each observed category in each draw
#'
#' @param eta A matrix, draws by rows.
#' @param tau A matrix, draws by thresholds.
#' @param y The category codes of the rows.
#' @param ord The ordinal record.
#' @returns A matrix, draws by rows.
#' @noRd
#' @keywords internal
.fb_kfold_category_loglik <- function(eta, tau, y, ord) {
  out <- matrix(NA_real_, nrow(eta), ncol(eta))
  pick <- cbind(seq_along(y), y)
  for (s in seq_len(nrow(eta))) {
    p <- .fb_ordinal_probs(eta[s, ], tau[s, ], ord$family, ord$link)
    out[s, ] <- log(p[pick])
  }
  out
}
