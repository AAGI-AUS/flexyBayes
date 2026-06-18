# S3 methods for flexybayes_glm (overrides of glm/lm methods with Bayesian
# semantics)

#' Summary for flexybayes GLM-compatible object
#'
#' Similar to `summary.glm()` but with Bayesian posterior statistics
#' instead of p-values.
#'
#' @param object A flexybayes_glm object (accessed via `fit$glm`)
#' @param ... Additional arguments (ignored)
#' @export
summary.flexybayes_glm <- function(object, ...) {
  beta <- object$coefficients
  V <- attr(object, "posterior_vcov")

  if (length(beta) == 0 || is.null(V) || nrow(V) == 0) {
    cat("No fixed effects to summarise.\n")
    return(invisible(NULL))
  }

  se <- sqrt(diag(V))
  # Probability of direction (proportion of posterior on same side as mean)
  # Approximate using normal assumption
  pd <- pnorm(abs(beta / se))

  coef_table <- cbind(
    Estimate = beta,
    `Post.SD` = se,
    `z value` = beta / se,
    `Pr(>|z|)` = 2 * (1 - pd)
  )

  cat("\nBayesian GLM summary (flexyBayes)\n")
  cat("Family:", object$family$family, "\n")
  cat("Link:  ", object$family$link, "\n\n")
  cat("Coefficients (posterior means):\n")
  printCoefmat(
    coef_table,
    P.values = TRUE,
    has.Pvalue = TRUE,
    signif.stars = FALSE
  )
  cat(
    "\n(Note: Pr(>|z|) is approximate posterior probability of direction,\n",
    " not a frequentist p-value.)\n"
  )
  cat("Residual df:", object$df.residual, "\n")

  invisible(coef_table)
}

#' Credible intervals for flexybayes_glm
#'
#' Returns posterior quantile-based credible intervals.
#'
#' @param object A flexybayes_glm object
#' @param parm Parameter names
#' @param level Credible level
#' @param ... Additional arguments (ignored)
#' @export
confint.flexybayes_glm <- function(object, parm = NULL, level = 0.95, ...) {
  beta <- object$coefficients
  V <- attr(object, "posterior_vcov")

  if (length(beta) == 0 || is.null(V)) {
    return(matrix(nrow = 0, ncol = 2))
  }

  se <- sqrt(diag(V))
  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2)

  ci_mat <- cbind(
    beta - z * se,
    beta + z * se
  )
  rownames(ci_mat) <- names(beta)
  colnames(ci_mat) <- paste0(round(c(alpha / 2, 1 - alpha / 2) * 100, 1), "%")

  if (!is.null(parm)) {
    ci_mat <- ci_mat[parm, , drop = FALSE]
  }

  ci_mat
}
