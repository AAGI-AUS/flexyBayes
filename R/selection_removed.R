# selection_removed.R: the penalised-selection functions removed in 0.16.0.


#' Penalised selection functions removed from flexyBayes
#'
#' @description
#' These functions are deprecated. Penalised selection is no longer part of
#' flexyBayes. `fb_select()`, `fb_select_path()`,
#' `fb_select_refit()`, `fb_select_fold_data()`, `fb_penalty()`,
#' `fb_resamples()` and `fb_decision()` remain for this release so that a
#' call to them stops with a refusal that says so, of class
#' `flexybayes_refusal_selection_removed`. They will be removed in the next
#' release.
#'
#' A penalised regression over candidate fixed effects can be fitted on the
#' design whitened by the fitted marginal covariance, from
#' [fb_marginal_cov()], with any penalised-regression package. The section
#' "Penalised selection on the whitened design" of
#' `vignette("flexyBayes-06-two-stage", package = "flexyBayes")` shows it
#' with glmnet.
#'
#' @param ... Ignored.
#' @returns None; each function signals the refusal.
#' @keywords internal
#' @name fb_selection_removed
NULL


#' @rdname fb_selection_removed
#' @export
fb_select <- function(...) {
  .fb_selection_removed("fb_select()")
}

#' @rdname fb_selection_removed
#' @export
fb_select_path <- function(...) {
  .fb_selection_removed("fb_select_path()")
}

#' @rdname fb_selection_removed
#' @export
fb_select_refit <- function(...) {
  .fb_selection_removed("fb_select_refit()")
}

#' @rdname fb_selection_removed
#' @export
fb_select_fold_data <- function(...) {
  .fb_selection_removed("fb_select_fold_data()")
}

#' @rdname fb_selection_removed
#' @export
fb_penalty <- function(...) {
  .fb_selection_removed("fb_penalty()")
}

#' @rdname fb_selection_removed
#' @export
fb_resamples <- function(...) {
  .fb_selection_removed("fb_resamples()")
}

#' @rdname fb_selection_removed
#' @export
fb_decision <- function(...) {
  .fb_selection_removed("fb_decision()")
}


# .fb_selection_removed() --- the refusal every removed function raises.
.fb_selection_removed <- function(fn) {
  stop(.fb_refusal_condition(
    reason_code = "selection_removed",
    message = paste0(
      fn, " is no longer part of flexyBayes: penalised selection was ",
      "removed in version 0.16.0. To fit a penalised regression on the ",
      "fixed effects of a mixed model, whiten the design with the fit's ",
      "marginal covariance from fb_marginal_cov() and pass it to glmnet, ",
      "as shown in vignette(\"flexyBayes-06-two-stage\"), section ",
      "\"Penalised selection on the whitened design\"."
    ),
    removed = fn
  ))
}


.populate_refusal_registry_selection_removed <- function() {
  .register_refusal(
    reason_code = "selection_removed",
    description = paste0(
      "A penalised-selection function was called: fb_select(), ",
      "fb_select_path(), fb_select_refit(), fb_select_fold_data(), ",
      "fb_penalty(), fb_resamples() or fb_decision(). Penalised selection ",
      "was removed from flexyBayes in 0.16.0, and the two-stage vignette ",
      "shows a penalised regression on the whitened design with glmnet."
    ),
    message_template = "%s is no longer part of flexyBayes.",
    registered_in_adr = "PARTB_coupling_table_06Oct2026 (D6)",
    plan_field = NA_character_,
    since_version = "0.15.0.9000"
  )
}
