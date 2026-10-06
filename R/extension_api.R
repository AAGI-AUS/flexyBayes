# extension_api.R: what a package that extends flexyBayes may read from a fit.


#' Reading a fit, for packages that extend flexyBayes
#'
#' A package built on flexyBayes reads a fitted model through these
#' functions, [fb_refit()], [fb_marginal_cov()], [fb_isolated_lapply()] and
#' the refusal registry ([fb_register_refusals()], [fb_refuse()]), and never
#' through the package's internal functions. `fb_extension_api()` gives the
#' version of that set; a package checks it in its `.onLoad()` and refuses to
#' load against a version it was not written for.
#'
#' * `fb_model_terms()` returns the model description of a fit: the parsed
#'   fixed, random, residual and addition terms, the response, the family
#'   and the link, as the engines received them.
#' * `fb_response_values()` evaluates a response, a column name or an
#'   expression such as `"log(y)"`, over `data`.
#' * `fb_offset_values()` sums the offsets of a model over the rows of
#'   `data`, refusing an offset that is not a finite number on every row.
#' * `fb_model_weights()` returns the observation weights a model was
#'   fitted with.
#' * `fb_known_variance()` describes the known variance or covariance of a
#'   model and where it is carried in the data.
#' * `fb_engine_settings()` checks and returns sampler settings.
#' * `fb_convergence()` returns the convergence verdict that
#'   `summary(fit)` reports, with its evidence.
#' * `fb_inla_posterior_sample()` draws joint posterior samples of an INLA
#'   fit on the scale of the data.
#' * `fb_check_refit()` refuses a fit that cannot be refitted as the model
#'   it was fitted as.
#' * `fb_fixed_design()` builds the design matrix of a fixed formula.
#'
#' @param fit,object A flexyBayes fit.
#' @param terms A model description from `fb_model_terms()`, or a fit.
#' @param response The response as written in the model: a column name or
#'   an R expression in the columns of `data`.
#' @param offsets The offset terms of a model as written, for example
#'   `"offset(log(area))"`; `fb_model_terms(fit)$offsets`.
#' @param newdata `TRUE` when `data` holds new rows for a prediction, so
#'   that a missing offset column is refused with a message that says so.
#' @param engine `NULL`, an `fb_engine` object, or a named list of sampler
#'   settings.
#' @param n The number of posterior samples.
#' @param ... Further arguments to [INLA::inla.posterior.sample()], for
#'   example `seed`, `num.threads` or `selection`.
#' @param caller Text placed in front of a refusal's message, naming the
#'   function that asked.
#'
#' @returns
#' * `fb_extension_api()`: an integer, the version of the set. It is 1.
#' * `fb_model_terms()`: a list, or `NULL` for an object that carries none.
#' * `fb_response_values()`: a vector with one value per row of `data`, or
#'   `NULL` when the response cannot be evaluated there.
#' * `fb_offset_values()`: a numeric vector with one value per row of `data`,
#'   or `NULL` when the model has no offset.
#' * `fb_model_weights()`: the numeric weights, or `NULL` for an unweighted
#'   model.
#' * `fb_known_variance()`: `NULL` for a model without a known variance;
#'   otherwise a list with `values` (one variance per row of the fitted data),
#'   `matrix` (the known covariance, or `NULL` for a vector), `column` (the
#'   data column that carries each row's variance) and `inla_index` (the
#'   index of the row effect that carries it in an INLA formula).
#' * `fb_engine_settings()`: a named list, empty for `NULL`.
#' * `fb_convergence()`: a list with `engine`, `ok` (`TRUE`, `FALSE` or
#'   `NA`), `reasons` and `evidence`.
#' * `fb_inla_posterior_sample()`: the list [INLA::inla.posterior.sample()]
#'   returns, with the linear predictor and the hyperparameters on the
#'   scale of the data.
#' * `fb_check_refit()`: `fit`, invisibly.
#'
#' @name fb_extension
#'
#' @examples
#' fb_extension_api()
#' d <- data.frame(y = c(1, 4, 9), area = c(1, 2, 3))
#' fb_response_values("sqrt(y)", d)
#' fb_offset_values("offset(log(area))", d)
#' fb_engine_settings(list(chains = 2L))
#' fb_fixed_design(~ area, d)$X
NULL


#' @rdname fb_extension
#' @export
fb_extension_api <- function() {
  1L
}


#' @rdname fb_extension
#' @export
fb_model_terms <- function(fit) {
  .fb_fit_ir(fit)
}


#' @rdname fb_extension
#' @export
fb_response_values <- function(response, data) {
  .fb_response_values(list(response = response), data)
}


#' @rdname fb_extension
#' @export
fb_offset_values <- function(offsets, data, newdata = FALSE) {
  .fb_offset_total(offsets, data, newdata = newdata)
}


#' @rdname fb_extension
#' @export
fb_model_weights <- function(terms) {
  .fb_ir_weights(.fb_extension_terms(terms))
}


#' @rdname fb_extension
#' @export
fb_known_variance <- function(terms) {
  term <- .fb_known_variance_term(.fb_extension_terms(terms))
  if (is.null(term)) {
    return(NULL)
  }
  list(
    values = term$values,
    matrix = term$matrix,
    column = .FB_KNOWN_VARIANCE_COL,
    inla_index = .FB_INLA_KNOWN_INDEX
  )
}


#' @rdname fb_extension
#' @export
fb_engine_settings <- function(engine) {
  if (is.null(engine)) {
    return(list())
  }
  if (inherits(engine, "fb_engine")) {
    return(engine$opts %||% list())
  }
  if (!is.list(engine) || is.null(names(engine)) ||
        any(!nzchar(names(engine)))) {
    stop(
      "`engine` must be NULL or a named list of brms sampler settings.",
      call. = FALSE
    )
  }
  unknown <- setdiff(names(engine), .fb_engine_opt_names)
  if (length(unknown)) {
    stop(
      "`engine` contains unrecognised setting",
      if (length(unknown) == 1L) " " else "s ",
      paste(shQuote(unknown), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  engine
}


#' @rdname fb_extension
#' @export
fb_convergence <- function(fit) {
  .fb_convergence_verdict(fit)
}


#' @rdname fb_extension
#' @export
fb_inla_posterior_sample <- function(fit, n, ...) {
  if (!inherits(fit, "flexybayes_inla") || is.null(fit$inla)) {
    stop("`fit` must be a flexyBayes fit on INLA.", call. = FALSE)
  }
  .fb_inla_posterior_sample(fit, n, ...)
}


#' @rdname fb_extension
#' @export
fb_check_refit <- function(fit, caller = "") {
  recorded <- fit$extras$call_info$fixed
  if (inherits(recorded, "brmsformula")) {
    fb <- fit$fb %||% .fb_fit_ir(fit)
    .fb_brms_refuse_unfitted_parts(
      recorded,
      family = fb$family,
      link = fb$link,
      caller = caller
    )
  }
  invisible(fit)
}


# .fb_extension_terms() --- a model description, from itself or from a fit.
.fb_extension_terms <- function(terms) {
  if (inherits(terms, "flexybayes")) .fb_fit_ir(terms) else terms
}
