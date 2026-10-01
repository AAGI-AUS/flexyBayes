# selection_fold_refit.R -- what an INLA fold refit inside fb_select()
# computes.
#
# A fold refit is read for one thing: its variance components (the fold
# covariance), its fitted values (the weights of `weights_by`) or its
# posterior draws on the held-out rows (the Bayesian route). While such a
# refit runs, emit_inla() reads the plan set here and computes that readout
# only: no DIC, WAIC or latent-field marginals, posterior configurations only
# for draws, and no second fit testing an AR1 field. None of these changes
# the mode INLA finds or the readout taken from it.
#
# INLA's optimiser starts from its own default. Started from the base fit's
# mode, it stopped at a lower marginal likelihood than from the default (by
# 0.006 log units with an identified AR1 field, and by 0.38 with an
# unidentified one, where a variance moved by 7 per cent).

.fb_fold_refit_state <- new.env(parent = emptyenv())


#' The readout plan of an INLA fold refit
#'
#' @param readout `"varcomp"`, `"fitted"` or `"draws"`.
#' @returns A list: `readout`, and `config`, whether INLA keeps the posterior
#'   configurations that draws are sampled from.
#' @noRd
#' @keywords internal
.fb_fold_refit_plan <- function(readout) {
  readout <- match.arg(readout, c("varcomp", "fitted", "draws"))
  list(readout = readout, config = identical(readout, "draws"))
}


#' Evaluate a fold refit under a readout plan
#'
#' @param plan A `.fb_fold_refit_plan()` record, or `NULL` for a full fit.
#' @param expr The refit.
#' @returns The value of `expr`.
#' @noRd
#' @keywords internal
.fb_with_fold_refit <- function(plan, expr) {
  previous <- .fb_fold_refit_state$plan
  .fb_fold_refit_state$plan <- plan
  on.exit(.fb_fold_refit_state$plan <- previous, add = TRUE)
  force(expr)
}


#' The readout plan in force, or `NULL` outside a fold refit
#'
#' @returns A `.fb_fold_refit_plan()` record or `NULL`.
#' @noRd
#' @keywords internal
.fb_fold_refit_active <- function() {
  .fb_fold_refit_state$plan
}


#' INLA's `control.compute` for a fit
#'
#' @param plan A `.fb_fold_refit_plan()` record, or `NULL` for a full fit.
#' @returns The `control.compute` list.
#' @noRd
#' @keywords internal
.fb_inla_control_compute <- function(plan = NULL) {
  if (is.null(plan)) {
    return(list(
      config = TRUE,
      return.marginals = TRUE,
      dic = TRUE,
      waic = TRUE
    ))
  }
  list(
    config = isTRUE(plan$config),
    return.marginals = FALSE,
    dic = FALSE,
    waic = FALSE
  )
}
