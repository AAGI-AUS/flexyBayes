# S3 methods for flexybayes_glm (overrides of glm/lm methods with Bayesian
# semantics)

#' Summary for flexybayes GLM-compatible object
#'
#' Similar to `summary.glm()` but with Bayesian posterior statistics
#' instead of p-values.
#'
#' @param object A `flexybayes_glm` object, reached as `fit$glm` on any
#'   fitted `flexybayes` object.
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, the printed coefficient table as a data.frame, or
#'   `NULL` when the object carries no fixed effects.
#' @export
summary.flexybayes_glm <- function(object, ...) {
  beta <- object$coefficients
  V <- attr(object, "posterior_vcov")

  if (length(beta) == 0 || is.null(V) || nrow(V) == 0) {
    cat("No fixed effects to summarise.\n")
    return(invisible(NULL))
  }

  se <- sqrt(diag(V))
  draws <- attr(object, "posterior_draws")
  exact <- !is.null(draws) && is.matrix(draws) && nrow(draws) > 0L &&
    identical(colnames(draws), names(beta))

  # Probability of direction: the posterior mass on the side of zero the
  # bulk of the distribution falls on, so it ranges over [0.5, 1].
  #
  # This column previously reported 2 * (1 - pd) under the heading
  # `Pr(>|z|)`, described in the footnote as the probability of
  # direction. It is neither: it is a doubled posterior sign-tail area,
  # formatted to look like a frequentist p-value. The heading is gone,
  # and the quantity is now computed from the draws where they exist.
  if (exact) {
    dm <- as.matrix(unname(draws))
    pd <- pmax(colMeans(dm > 0), colMeans(dm < 0))
    ci <- t(apply(dm, 2L, stats::quantile, probs = c(0.025, 0.975)))
    pd_label <- "pd"
    basis <- "posterior draws"
  } else {
    pd <- stats::pnorm(abs(beta / se))
    z <- stats::qnorm(0.975)
    ci <- cbind(beta - z * se, beta + z * se)
    pd_label <- "pd_normal_approx"
    basis <- "normal approximation (no draws available)"
  }

  coef_table <- cbind(
    Estimate = beta,
    `Post.SD` = se,
    `2.5%` = ci[, 1L],
    `97.5%` = ci[, 2L],
    pd
  )
  colnames(coef_table)[5L] <- pd_label
  rownames(coef_table) <- names(beta)

  cat("\nBayesian GLM summary (flexyBayes)\n")
  cat("Family:", object$family$family, "\n")
  cat("Link:  ", object$family$link, "\n\n")
  cat("Coefficients (posterior means, 95% credible intervals):\n")
  print(round(coef_table, 4L))
  cat("\nBasis:", basis, "\n")
  cat(
    "pd = posterior probability of direction, in [0.5, 1]. It is not a\n",
    "p-value and no null hypothesis is being tested.\n"
  )
  cat("Residual df:", object$df.residual, "\n")

  invisible(coef_table)
}

#' Credible intervals for flexybayes_glm
#'
#' Posterior quantile credible intervals, computed from the fixed-effect
#' draws.
#'
#' Where the draws are unavailable the method falls back to the
#' normal-approximation interval \eqn{\hat\beta \pm z_{\alpha/2}
#' \mathrm{sd}}, and marks the returned matrix with
#' `attr(, "interval_basis") == "normal_approximation"`. The two agree
#' only for a symmetric posterior; on a skewed one they can differ
#' materially, which is why the basis is reported rather than assumed.
#' Earlier versions documented this method as quantile-based while
#' always returning the approximation.
#'
#' @param object A `flexybayes_glm` object, reached as `fit$glm` on any
#'   fitted `flexybayes` object.
#' @param parm A character vector naming the parameters to report, or
#'   `NULL` (the default) for every fixed effect.
#' @param level A single numeric giving the credible level, defaulting
#'   to 0.95.
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns A two-column matrix of interval bounds, carrying an
#'   `interval_basis` attribute of `"posterior_quantile"` or
#'   `"normal_approximation"`.
#' @export
confint.flexybayes_glm <- function(object, parm = NULL, level = 0.95, ...) {
  beta <- object$coefficients
  V <- attr(object, "posterior_vcov")

  if (length(beta) == 0 || is.null(V)) {
    return(matrix(nrow = 0, ncol = 2))
  }

  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)
  draws <- attr(object, "posterior_draws")
  exact <- !is.null(draws) && is.matrix(draws) && nrow(draws) > 0L &&
    identical(colnames(draws), names(beta))

  if (exact) {
    # as.matrix() drops the draws object's own dimnames names ("draw" /
    # "variable"), which apply() otherwise carries into the result and
    # prints as a stray header above the coefficient names.
    ci_mat <- t(apply(as.matrix(unname(draws)), 2L, stats::quantile,
                      probs = probs))
    basis <- "posterior_quantile"
  } else {
    se <- sqrt(diag(V))
    z <- stats::qnorm(1 - alpha / 2)
    ci_mat <- cbind(beta - z * se, beta + z * se)
    basis <- "normal_approximation"
  }

  rownames(ci_mat) <- names(beta)
  colnames(ci_mat) <- paste0(round(probs * 100, 1), "%")

  if (!is.null(parm)) {
    ci_mat <- ci_mat[parm, , drop = FALSE]
  }

  attr(ci_mat, "interval_basis") <- basis
  ci_mat
}
