# fb_select() --- the public verb
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md sections 2, 3.1, 11.
#
# The partition principle, which is why this function exists in this shape.
# The base model, the extended base model and a penalised selection are ONE
# computation under three partitions of the design:
#
#   base model M0      everything protected, no candidates, no penalty
#   extended base M1   base plus named candidates, still no penalty
#   selection          base protected, candidate block penalised
#
# Everything expensive is shared -- theta-hat, V, the whitening, the folds,
# the prediction estimand, the scoring unit and the certificate. What differs
# is where the vertical line falls in the design matrix.
#
# So `penalty = NULL` is not a degenerate case kept for tidiness. It is the
# lambda = 0 corner of the penalised problem, it exercises every part of the
# substrate except the provider adapter, and it is what makes the eventual
# penalty a small diff against machinery already proved rather than a large
# one against machinery being written at the same time.
#
# This release implements that corner only. A penalty refuses by name.


# --- the gate -------------------------------------------------------- #

# .fb_select_gate() --- the seven-axis check, answered before anything runs.
#
# A base is admitted only when every axis passes, and the refusal names the
# measurement rather than the conclusion. `us()` and `fa()` on INLA stay out
# regardless of how a particular fit behaved: the measured genetic
# correlation is 0.142 against ASReml's 0.930, and multiplying the response
# by ten -- an exact transformation under which a correlation is invariant --
# does not reproduce it.
.fb_select_gate <- function(base) {
  if (!inherits(base, "flexybayes")) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_not_allowlisted",
      message = "fb_select(): `base` must be a fitted flexyBayes model."
    ))
  }
  if (is.null(base$fb) || is.null(base$data)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_shape_unreadable",
      message = paste0(
        "fb_select(): this base fit does not carry the parsed model and the ",
        "data the selection substrate reads. Fits on the aggregated route, ",
        "which an intercept-only fixed model takes, have a different shape."
      )
    ))
  }

  fam <- base$fb$family %||% NA_character_
  link <- base$fb$link %||% NA_character_
  if (!identical(fam, "gaussian") || !identical(link, "identity")) {
    stop(.fb_refusal_condition(
      reason_code = "selection_family_not_supported",
      message = paste0(
        "fb_select(): the whitening identity is exact only for a Gaussian ",
        "response with an identity link, and this fit is ", fam, "/", link,
        ". For other families the observation variance depends on the mean, ",
        "no single V survives integrating the random effects, and an IRLS ",
        "formulation changes the weights as the coefficients change. A ",
        "penalised-likelihood route is a separate programme with its own ",
        "validation."
      )
    ))
  }

  # Structures whose covariance is not admitted, by name and with the
  # measurement that put them there.
  types <- vapply(
    base$fb$random_terms %||% list(),
    function(t) t$type %||% NA_character_, character(1L)
  )
  bad <- intersect(types, c("us", "fa", "unstructured", "factor_analytic"))
  if (length(bad) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_not_allowlisted",
      message = paste0(
        "fb_select(): the base carries ", paste(bad, collapse = ", "),
        ", which is not on the selection allowlist. INLA's us() returns a ",
        "genetic correlation of 0.142 where ASReml returns 0.930, and ",
        "rescaling the response by ten does not reproduce it, so a rotation ",
        "built on that covariance would be badly misspecified while looking ",
        "perfectly reproducible."
      )
    ))
  }

  invisible(TRUE)
}


# --- the certificate --------------------------------------------------- #

# .fb_select_certificate() --- what was done, and every place it was
# approximate.
#
# Contract section 11.1. This is the package's strongest existing idea --
# stated inability -- carried into selection. A certificate that recorded
# only the happy path would be decoration; the fields that earn it are
# `approximations` and `refusals`, which say where the answer is weaker than
# it looks.
.fb_select_certificate <- function(base, cv, resamples, candidates,
                                   covariance, penalty) {
  conv <- .fb_convergence_verdict(base)
  vc <- try(.fb_sel_variance_components(base, covariance), silent = TRUE)
  bs <- try(.fb_sel_block_structure(base), silent = TRUE)

  structure(
    list(
      base = list(
        engine = base$fb$source %||% NA_character_,
        family = base$fb$family,
        link = base$fb$link,
        n_observed = nrow(base$data),
        convergence = conv$ok,
        convergence_reasons = conv$reasons
      ),
      covariance = list(
        estimand = covariance,
        builder = if (inherits(bs, "try-error") || !isTRUE(bs$applicable)) {
          "random_term_reduction"
        } else {
          "block_cholesky"
        },
        components = if (inherits(vc, "try-error")) NULL else vc$terms,
        residual = if (inherits(vc, "try-error")) NA_real_ else vc$residual
      ),
      design = list(
        protected = attr(cv, "protected") %||% character(0),
        candidates = if (is.null(candidates)) character(0) else all.vars(candidates)
      ),
      resampling = list(
        unit = resamples$unit,
        v = resamples$v,
        repeats = resamples$repeats,
        n_folds = length(cv$folds),
        prediction = cv$prediction,
        scoring_unit = cv$aggregated$se_unit[[1L]] %||% NA_character_,
        seed = cv$seed
      ),
      penalty = if (is.null(penalty)) "none: every candidate protected" else "penalised",
      # Named because a reader must be able to see what was not exact.
      approximations = character(0),
      refusals = character(0),
      contract = "SPEC_penalised_selection_v3_2026-09-11"
    ),
    class = "fb_select_certificate"
  )
}


# --- the verb ---------------------------------------------------------- #

#' Covariance-aware selection over a fitted flexyBayes model
#'
#' Scores a candidate block against a protected base model, out of fold, on
#' data whitened by the base fit's own covariance. Folds are assigned on
#' original observational units before any rotation, every quantity learned
#' from data is learned inside the training fold, and the prediction estimand
#' is named rather than implied.
#'
#' @param base A fitted flexyBayes model. Its fixed effects are the protected
#'   block and are never penalised.
#' @param candidates A one-sided formula over columns of the base fit's data.
#' @param penalty Reserved. Only `NULL` is implemented in this release, which
#'   protects every candidate and gives the extended base model, its
#'   out-of-fold score, and the increment over the base.
#' @param resamples A resampling specification; see the internal
#'   `fb_resamples()` constructor.
#' @param prediction `"new_cluster"` (default) or `"known_cluster"`. A
#'   held-out cluster's own effect is not estimable from the training rows,
#'   so scoring it conditionally would reward the model for a random effect
#'   it was handed.
#' @param covariance The covariance summary used to whiten. The default,
#'   `"plugin_mean_covariance"`, is exact: V is linear in the variances, so
#'   `E[V]` needs only the marginals.
#' @param metric `"rmse"`, `"mse"` or `"nlpd"`.
#' @param seed Required whenever the resampling draws more than once.
#'
#' @return An object of class `fb_select` carrying the out-of-fold scores,
#'   the increment over the protected model, and the selection certificate.
#'
#' @export
fb_select <- function(base,
                      candidates,
                      penalty = NULL,
                      resamples = fb_resamples(),
                      prediction = c("new_cluster", "known_cluster"),
                      covariance = "plugin_mean_covariance",
                      metric = c("rmse", "mse", "nlpd"),
                      seed = NULL) {
  prediction <- match.arg(prediction)
  metric <- match.arg(metric)
  .fb_select_gate(base)

  if (!is.null(penalty)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_penalty_not_implemented",
      message = paste0(
        "fb_select(): a penalty is not implemented in this release. ",
        "`penalty = NULL` protects every candidate and returns the extended ",
        "base model with its out-of-fold score and increment, which is the ",
        "lambda = 0 corner of the penalised problem."
      )
    ))
  }

  fixed_labels <- vapply(
    base$fb$fixed_terms %||% list(),
    function(t) t$label %||% t$var %||% NA_character_, character(1L)
  )
  fixed_labels <- fixed_labels[!is.na(fixed_labels) & nzchar(fixed_labels)]
  protected <- stats::as.formula(paste0(
    "~ ", if (length(fixed_labels) > 0L) paste(fixed_labels, collapse = " + ") else "1"
  ))
  random_vars <- unique(unlist(lapply(
    base$fb$random_terms %||% list(),
    function(t) if (identical(t$type, "simple")) t$var else c(t$outer, t$inner)
  )))

  cv <- .fb_sel_cv(
    data = base$data,
    response = base$fb$response,
    protected = protected,
    candidates = candidates,
    random = random_vars %||% character(0),
    resamples = resamples,
    seed = seed,
    prediction = prediction,
    metric = metric
  )

  agg <- cv$aggregated
  base_row <- agg[agg$model == "protected_only", ]
  full_row <- agg[agg$model == "full", ]
  increment <- list(
    protected_only = base_row$score,
    full = full_row$score,
    # rmse and mse fall as a model improves, so the increment is signed to
    # read the same way whichever metric was asked for: positive is better.
    improvement = base_row$score - full_row$score,
    metric = metric,
    se_unit = full_row$se_unit,
    n_units = full_row$n_units
  )

  structure(
    list(
      base = base,
      scores = agg,
      clusters = cv$clusters,
      increment = increment,
      folds = cv$folds,
      certificate = .fb_select_certificate(
        base, cv, resamples, candidates, covariance, penalty
      )
    ),
    class = "fb_select"
  )
}


#' @export
print.fb_select <- function(x, digits = 4L, ...) {
  cert <- x$certificate
  cat("<fb_select> ", cert$penalty, "\n", sep = "")
  cat(
    "  base      : ", cert$base$family, "/", cert$base$link,
    ", ", cert$base$n_observed, " rows, convergence ",
    if (is.na(cert$base$convergence)) "UNVERIFIED"
    else if (cert$base$convergence) "ok" else "FAILED",
    "\n",
    sep = ""
  )
  cat("  covariance: ", cert$covariance$builder, " / ",
      cert$covariance$estimand, "\n", sep = "")
  cat("  resampling: ", cert$resampling$n_folds, " folds, prediction ",
      cert$resampling$prediction, ", scored at ",
      cert$resampling$scoring_unit, "\n", sep = "")
  cat("\n  out-of-fold ", x$increment$metric, ":\n", sep = "")
  cat("    protected only : ", round(x$increment$protected_only, digits), "\n", sep = "")
  cat("    with candidates: ", round(x$increment$full, digits), "\n", sep = "")
  cat("    improvement    : ", round(x$increment$improvement, digits),
      "  (", x$increment$n_units, " independent ",
      x$increment$se_unit, ")\n", sep = "")
  if (length(cert$approximations) > 0L) {
    cat("  approximations: ", paste(cert$approximations, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}


#' @export
glance.fb_select <- function(x, ...) {
  data.frame(
    metric = x$increment$metric,
    protected_only = x$increment$protected_only,
    full = x$increment$full,
    improvement = x$increment$improvement,
    n_units = x$increment$n_units,
    se_unit = x$increment$se_unit,
    n_folds = length(x$folds),
    prediction = x$certificate$resampling$prediction,
    covariance = x$certificate$covariance$estimand,
    builder = x$certificate$covariance$builder,
    stringsAsFactors = FALSE
  )
}


#' @export
print.fb_select_certificate <- function(x, ...) {
  cat("<fb_select_certificate> ", x$contract, "\n", sep = "")
  cat("  covariance : ", x$covariance$builder, " / ", x$covariance$estimand, "\n", sep = "")
  cat("  prediction : ", x$resampling$prediction, "\n", sep = "")
  cat("  scored at  : ", x$resampling$scoring_unit, "\n", sep = "")
  cat("  seed       : ", x$resampling$seed %||% "none", "\n", sep = "")
  cat("  approximations: ",
      if (length(x$approximations) == 0L) "none declared"
      else paste(x$approximations, collapse = ", "), "\n", sep = "")
  invisible(x)
}
