# inla_best_of.R -- every INLA fit is run k times and the fit with the highest
# marginal likelihood is kept.
#
# Identical one-thread INLA calls in fresh R processes can reach different
# posterior modes, while repeated fits within one R session agree
# (inst/KNOWN_ISSUES.md, item 26). The k fits are therefore the fit in this
# session and k - 1 fits in fresh R processes, run one after another through
# .fb_isolated_call(). The spread of the k fits is recorded at
# `$extras$reproducibility$best_of`.

# Relative difference of a standard deviation across the k fits above which
# print() and summary() say so.
.FB_BEST_OF_SD_NOTE <- 0.05

#' The number of INLA fits each model is run with
#'
#' Called before the first fit.
#'
#' @return A whole number, one or more: option `flexyBayes.inla_best_of`,
#'   3 when unset.
#' @noRd
#' @keywords internal
.fb_inla_best_of_k <- function() {
  value <- getOption("flexyBayes.inla_best_of", 3L)
  ok <- is.numeric(value) &&
    length(value) == 1L &&
    is.finite(value) &&
    value >= 1 &&
    value == round(value)
  if (!ok) {
    stop(.fb_refusal_condition(
      reason_code = "inla_best_of_invalid",
      message = paste0(
        "option flexyBayes.inla_best_of must be a whole number, one or more ",
        "(the number of INLA fits each model is run with). It is: ",
        .fb_prior_deparse(value),
        "."
      )
    ))
  }
  k <- as.integer(value)
  if (k > 1L && !.fb_isolation_callr_available()) {
    stop(.fb_refusal_condition(
      reason_code = "isolation_requires_callr",
      message = paste0(
        "flexybayes(): option flexyBayes.inla_best_of = ",
        k,
        " runs ",
        k - 1L,
        " of the INLA fits in separate R processes through callr, and ",
        "callr is not installed. Install it with install.packages(\"callr\"), ",
        "or set options(flexyBayes.inla_best_of = 1) to fit once."
      ),
      entry_point = "flexybayes()"
    ))
  }
  k
}


#' One INLA fit in a fresh R process
#'
#' Tests replace it to return fits with chosen marginal likelihoods.
#'
#' @param args The arguments of `INLA::inla()`.
#' @return The fit the separate process returned.
#' @noRd
#' @keywords internal
.fb_inla_isolated_fit <- function(args) {
  .fb_isolated_call(
    .fb_inla_child_fit,
    list(args = args),
    engine = "inla",
    entry_point = "flexybayes()"
  )
}

# Defined at the top level so that only `args`, and no frame of the caller,
# is sent to the separate process.
.fb_inla_child_fit <- function(args) {
  .inla_call(args)
}


#' The log marginal likelihood of an INLA fit, `-Inf` when it has none
#' @noRd
.fb_best_of_mlik <- function(fit) {
  m <- fit$mlik
  value <- if (is.null(m) || !length(m)) NA_real_ else as.numeric(m[1L, 1L])
  if (is.finite(value)) value else -Inf
}


#' The posterior mean of each standard deviation of an INLA fit
#'
#' @return A named numeric vector, one entry per `Precision for` row of the
#'   hyperparameter table, `NA` where the fit kept no marginal.
#' @noRd
.fb_best_of_sds <- function(fit) {
  hp <- fit$summary.hyperpar
  nm <- if (is.null(hp)) character(0) else rownames(hp)
  nm <- nm[startsWith(nm, "Precision for ")]
  marginals <- fit$marginals.hyperpar %||% list()
  out <- vapply(
    nm,
    function(x) {
      m <- marginals[[x]]
      if (is.null(m)) {
        return(NA_real_)
      }
      INLA::inla.emarginal(function(p) 1 / sqrt(p), m)
    },
    numeric(1)
  )
  names(out) <- sub("^Precision for ", "", nm)
  out
}


#' Run the remaining k - 1 fits and keep the best of the k
#'
#' @param fit The fit made in this session.
#' @param args The arguments it was made with.
#' @param k The number of fits, from `.fb_inla_best_of_k()`, read before
#'   the first fit.
#' @param fa_term The model's fa(f, k):g term, or `NULL`. With one, the
#'   record also holds `fa_gate`, the agreement of the k fits on the
#'   genetic covariance (`.fb_fa_gate()`).
#' @return A list: `fit`, the fit with the highest marginal likelihood (the
#'   first of equals), and `record`, `NULL` when `k` is 1, otherwise
#'   `k`, `mlik`, `mode_status`, `chosen`, `mlik_range`, `sd` (one row per
#'   fit) and `sd_rel_diff`, the largest difference of each standard
#'   deviation across the fits relative to its smallest value.
#' @noRd
#' @keywords internal
.fb_inla_best_of <- function(fit, args, k, fa_term = NULL) {
  if (k == 1L) {
    return(list(fit = fit, record = NULL))
  }
  mlik <- numeric(k)
  status <- integer(k)
  sds <- vector("list", k)
  fits <- vector("list", k)
  best <- fit
  chosen <- 1L
  for (j in seq_len(k)) {
    if (j > 1L) {
      fit <- tryCatch(
        .fb_inla_isolated_fit(args),
        error = function(e) {
          e$message <- paste0(
            "INLA fit ",
            j,
            " of ",
            k,
            ", run in a separate R process, failed: ",
            conditionMessage(e)
          )
          stop(e)
        }
      )
    }
    mlik[[j]] <- .fb_best_of_mlik(fit)
    status[[j]] <- as.integer(fit$mode$mode.status %||% NA_integer_)
    sds[[j]] <- .fb_best_of_sds(fit)
    if (!is.null(fa_term)) {
      fits[[j]] <- fit
    }
    if (mlik[[j]] > mlik[[chosen]]) {
      best <- fit
      chosen <- j
    }
  }
  sd_mat <- do.call(rbind, sds)
  sd_rel_diff <- if (is.null(sd_mat) || !ncol(sd_mat)) {
    numeric(0)
  } else {
    apply(sd_mat, 2L, function(s) (max(s) - min(s)) / min(s))
  }
  mlik[!is.finite(mlik)] <- NA_real_
  list(
    fit = best,
    record = list(
      k = k,
      mlik = mlik,
      mode_status = status,
      chosen = chosen,
      mlik_range = if (all(is.na(mlik))) {
        NA_real_
      } else {
        diff(range(mlik, na.rm = TRUE))
      },
      sd = sd_mat,
      sd_rel_diff = sd_rel_diff,
      fa_gate = if (is.null(fa_term)) NULL else .fb_fa_gate(fits, fa_term)
    )
  )
}


#' Move the recorded marginal likelihoods to the original response scale
#'
#' @param record A `.fb_inla_best_of()` record, or `NULL`.
#' @param shift What the back-transform subtracts from the fit's marginal
#'   likelihood.
#' @noRd
.fb_best_of_shift <- function(record, shift) {
  if (!is.null(record)) {
    record$mlik <- record$mlik - shift
  }
  record
}


#' The line print() and summary() show when the k fits disagree
#'
#' @param record A `.fb_inla_best_of()` record.
#' @return A string, or `NULL` when no standard deviation differs by more
#'   than 5% across the fits.
#' @noRd
.fb_best_of_line <- function(record) {
  d <- record$sd_rel_diff
  if (!length(d) || !any(d > .FB_BEST_OF_SD_NOTE, na.rm = TRUE)) {
    return(NULL)
  }
  worst <- which.max(d)
  paste0(
    "fit ",
    record$chosen,
    " of ",
    record$k,
    " kept (highest marginal likelihood); the SD for ",
    names(d)[[worst]],
    " differed by ",
    format(round(100 * d[[worst]], 1), nsmall = 1),
    "% across the ",
    record$k,
    " fits"
  )
}


# Thrown by .fb_inla_best_of_k() and triggered in
# tests/testthat/test-inla-best-of.R.
.populate_refusal_registry_inla_best_of <- function() {
  .register_refusal(
    reason_code = "inla_best_of_invalid",
    description = paste0(
      "Option flexyBayes.inla_best_of, the number of INLA fits each model ",
      "is run with, is not a whole number of one or more. The request is ",
      "refused before anything is fitted."
    ),
    message_template = paste0(
      "option flexyBayes.inla_best_of must be a whole number, one or more. ",
      "It is: %s."
    ),
    registered_in_adr = "PLAN_0.17.0_repeatability_fa_speed_09Oct2026 (R1)",
    plan_field = NA_character_,
    since_version = "0.17.0"
  )
}
