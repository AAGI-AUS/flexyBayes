# vpredict.R -- exact posterior functionals of the variance components.
#
# ASReml's vpredict() takes a function of variance components -- a
# heritability, a genetic correlation, a variance ratio -- and returns the
# estimate with a DELTA-METHOD standard error. The delta method linearises
# the function at the REML estimate, so it does not know where the
# parameter space ends: an ASReml heritability can come back with an
# interval running below zero, and a correlation with one running past one.
# A8 §2.3 has a worked instance, a specific variance of 0.017 with a
# standard error of 0.137.
#
# The Bayesian version is not an approximation of anything. Every variance
# component already carries draws, so a function of them is evaluated draw
# by draw and the answer is that function's exact posterior. It respects
# the boundary because every draw does, and it needs no asymptotics.
#
# The same call works on either engine because it reads canonical draws,
# not engine-native ones: .fb_canonical_draws() has already renamed and
# put every variance component on the standard-deviation scale.

#' Posterior functionals of the variance components
#'
#' The exact posterior of any function of the fitted variance components --
#' a heritability, a genetic correlation, a variance ratio, a repeatability.
#' Each functional is written as a two-sided formula whose left side names
#' the result and whose right side is an expression over the component
#' names that [summary()] reports.
#'
#' Components are on the **standard-deviation** scale, as everywhere else
#' in this package, so a variance ratio squares them. Writing a name the
#' fit does not carry is refused, and the refusal lists the names it does.
#'
#' The interval is a quantile interval of the functional's own posterior.
#' Where the delta method would put a heritability's interval below zero,
#' this one cannot: every draw of `sd_g^2 / (sd_g^2 + sigma^2)` lies in
#' \[0, 1\], so the interval does too.
#'
#' @param object A fitted `flexybayes` model.
#' @param ... One or more two-sided formulas, `name ~ expression`.
#' @param n_draws Number of posterior draws to evaluate over.
#' @param prob Width of the credible interval.
#' @param draws When `TRUE`, return the matrix of functional draws (one
#'   column per functional) instead of the summary.
#'
#' @returns A `data.frame` with one row per functional and the columns
#'   `functional`, `estimate` (posterior mean), `median`, `std.error`
#'   (posterior standard deviation), `conf.low` and `conf.high`. With
#'   `draws = TRUE`, a numeric matrix of draws instead.
#'
#' @examples
#' \dontrun{
#' fit <- flexybayes(yield ~ 1, random = ~ gen, data = trial)
#' vpredict(fit, h2 ~ sd_gen^2 / (sd_gen^2 + sigma^2))
#' }
#'
#' @export
vpredict <- function(object, ...) {
  UseMethod("vpredict")
}

#' @rdname vpredict
#' @export
vpredict.flexybayes <- function(
  object,
  ...,
  n_draws = 1000L,
  prob = 0.95,
  draws = FALSE
) {
  specs <- list(...)
  if (!length(specs)) {
    stop(.fb_refusal_condition(
      reason_code = "vpredict_no_functional",
      message = paste0(
        "vpredict() needs at least one functional, written as a two-sided ",
        "formula: vpredict(fit, h2 ~ sd_gen^2 / (sd_gen^2 + sigma^2))."
      )
    ))
  }
  if (!is.numeric(prob) || length(prob) != 1L || prob <= 0 || prob >= 1) {
    stop(.fb_refusal_condition(
      reason_code = "vpredict_no_functional",
      message = "`prob` must be a single number strictly between 0 and 1."
    ))
  }

  post <- .fb_canonical_draws(object, n_samples = n_draws)
  post <- post[vapply(post, is.numeric, logical(1L))]
  if (!length(post)) {
    stop(.fb_refusal_condition(
      reason_code = "vpredict_component_not_in_posterior",
      message = paste0(
        "This fit carries no posterior draws to build a functional from."
      )
    ))
  }
  n <- length(post[[1L]])

  # baseenv(), not the calling frame: a component name the fit does not
  # carry must refuse by name rather than silently resolve to whatever the
  # caller happens to have in their workspace under that name.
  env <- list2env(post, parent = baseenv())

  labels <- character(0)
  values <- list()
  for (i in seq_along(specs)) {
    spec <- specs[[i]]
    if (!inherits(spec, "formula") || length(spec) != 3L) {
      stop(.fb_refusal_condition(
        reason_code = "vpredict_no_functional",
        message = paste0(
          "Each functional is a two-sided formula, `name ~ expression`. ",
          "Argument ", i, " is not one."
        )
      ))
    }
    label <- paste(deparse(spec[[2L]]), collapse = "")
    expr <- spec[[3L]]

    unknown <- setdiff(all.vars(expr), names(post))
    if (length(unknown)) {
      stop(.fb_refusal_condition(
        reason_code = "vpredict_component_not_in_posterior",
        message = paste0(
          "`", label, "` names ",
          if (length(unknown) > 1L) "components " else "a component ",
          paste0("`", unknown, "`", collapse = ", "),
          ", which this fit does not carry. It carries: ",
          paste(sort(names(post)), collapse = ", "),
          ". Variance components are on the standard-deviation scale, so a ",
          "variance ratio squares them."
        )
      ))
    }

    val <- eval(expr, envir = env)
    if (!is.numeric(val) || length(val) != n) {
      stop(.fb_refusal_condition(
        reason_code = "vpredict_functional_not_one_value_per_draw",
        message = paste0(
          "`", label, "` must evaluate to one number per posterior draw. ",
          "It returned ",
          if (!is.numeric(val)) "a non-numeric value" else {
            paste0(length(val), " values against ", n, " draws")
          },
          "."
        )
      ))
    }
    labels <- c(labels, label)
    values[[label]] <- as.numeric(val)
  }

  mat <- matrix(
    unlist(values, use.names = FALSE),
    nrow = n,
    dimnames = list(NULL, labels)
  )
  if (isTRUE(draws)) {
    return(mat)
  }

  lo <- (1 - prob) / 2
  q <- function(v, p) unname(stats::quantile(v, probs = p, names = FALSE))
  data.frame(
    functional = labels,
    estimate = vapply(values, mean, numeric(1L), USE.NAMES = FALSE),
    median = vapply(values, stats::median, numeric(1L), USE.NAMES = FALSE),
    std.error = vapply(values, stats::sd, numeric(1L), USE.NAMES = FALSE),
    conf.low = vapply(values, q, numeric(1L), p = lo, USE.NAMES = FALSE),
    conf.high = vapply(values, q, numeric(1L), p = 1 - lo, USE.NAMES = FALSE),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
