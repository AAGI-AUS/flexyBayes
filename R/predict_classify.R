# predict_classify.R -- the marginal-means table an ASReml user takes to a
# meeting.
#
# `predict(asreml_fit, classify = "Variety")` is the object a breeder
# reads, and the flexyBayes equivalent is `predict(fit, classify =
# "Variety")`. The table is built through the emmeans seam this package
# already registers (recover_data / emm_basis in R/emmeans_support.R), so
# the marginal means come from one reference-grid implementation rather
# than a second one written here.
#
# Two things the table is not, and says so on every print. It is not a
# standard-error-of-difference table: no SED, no LSD, no critical value.
# And on an INLA fit the interval is not INLA's own marginal but the
# Gaussian approximation of the joint fixed-effect posterior, which is
# what the emm_basis seam supplies; the banner names that rather than
# letting the table read as something it is not.


# ---------------------------------------------------------------- #
# Reading the classify argument                                     #
# ---------------------------------------------------------------- #

# .fb_classify_vars() --- the factors a classify request names.
#
# Accepts ASReml's own spelling ("Variety", "Variety:env") and the
# one-sided formula an R user reaches for (~ Variety, ~ Variety * env).
# The operator carries no meaning here: what a classify names is a set of
# factors whose combinations the table has one row for.
#
# @noRd
# @keywords internal
.fb_classify_vars <- function(classify) {
  if (inherits(classify, "formula")) {
    if (length(classify) != 2L) {
      stop(
        "predict(): `classify` given as a two-sided formula. It names the ",
        "factors the means table is broken down by, so it takes the ",
        "one-sided form -- classify = ~ Variety, or ~ Variety * env.",
        call. = FALSE
      )
    }
    vars <- all.vars(classify)
  } else if (is.character(classify)) {
    vars <- unlist(strsplit(classify, "[:*+]"))
    vars <- trimws(vars)
    vars <- vars[nzchar(vars)]
  } else {
    stop(
      "predict(): `classify` must be a character value (\"Variety\", ",
      "\"Variety:env\") or a one-sided formula (~ Variety), not an object ",
      "of class <", paste(class(classify), collapse = ", "), ">.",
      call. = FALSE
    )
  }
  if (length(vars) == 0L) {
    stop(
      "predict(): `classify` names no variables, so there is no table to ",
      "build.",
      call. = FALSE
    )
  }
  unique(vars)
}

# .fb_classify_column() --- one emmeans output column under our name.
#
# emmeans names the interval bounds by the degrees of freedom it used:
# `asymp.LCL` / `asymp.UCL` on the z-intervals this package's emm_basis
# asks for (df = Inf), `lower.CL` / `upper.CL` when a finite df is in
# play. Both spellings are read so the returned column names do not
# depend on which one emmeans chose.
#
# @noRd
# @keywords internal
.fb_classify_column <- function(tab, candidates) {
  hit <- candidates[candidates %in% names(tab)]
  if (length(hit) == 0L) {
    return(rep(NA_real_, nrow(tab)))
  }
  as.numeric(tab[[hit[[1L]]]])
}


# ---------------------------------------------------------------- #
# The table                                                         #
# ---------------------------------------------------------------- #

# .fb_predict_classify() --- build the marginal-means table.
#
# @noRd
# @keywords internal
.fb_predict_classify <- function(object, classify, level = 0.95) {
  .fb_check_level(level)
  .fb_check_classify_emmeans()

  vars <- .fb_classify_vars(classify)
  .fb_check_classify_random_factor(object, vars)
  spec <- stats::as.formula(paste("~", paste(vars, collapse = " * ")))
  emm <- emmeans::emmeans(object, specs = spec, level = level)
  tab <- as.data.frame(summary(emm, level = level))

  present <- vars[vars %in% names(tab)]
  out <- tab[, present, drop = FALSE]
  # C4/FS-26: on an INLA fit whose model touched a non-syntactic factor
  # level, the reference grid emmeans built (via R/emmeans_support.R's
  # recover_data(), which reads object$data) carries the legalised
  # (make.names()) labels; restore the user's own here so the table a
  # breeder reads shows what they wrote. object$level_labels is NULL
  # for every fit this legalisation never touched (every non-INLA
  # engine, and an INLA fit with no non-syntactic level), so this is a
  # no-op there.
  for (v in present) {
    if (is.factor(out[[v]])) {
      levels(out[[v]]) <- .inla_restore_level_labels(
        object$level_labels,
        v,
        levels(out[[v]])
      )
    } else if (is.character(out[[v]])) {
      out[[v]] <- .inla_restore_level_labels(object$level_labels, v, out[[v]])
    }
  }
  out$estimate <- .fb_classify_column(tab, c("emmean", "response", "rate"))
  out$std.error <- .fb_classify_column(tab, c("SE", "std.error"))
  out$conf.low <- .fb_classify_column(tab, c("asymp.LCL", "lower.CL"))
  out$conf.high <- .fb_classify_column(tab, c("asymp.UCL", "upper.CL"))
  rownames(out) <- NULL

  structure(
    out,
    class = c("fb_predict_classify", "data.frame"),
    classify = vars,
    level = level,
    engine = .fb_fit_engine(object)
  )
}

# .fb_emmeans_available() --- is the estimation seam's package present?
#
# A named predicate rather than an inline requireNamespace(), so the
# refusal below is reachable from a test without uninstalling emmeans.
#
# @noRd
# @keywords internal
.fb_emmeans_available <- function() {
  requireNamespace("emmeans", quietly = TRUE)
}

# .fb_check_classify_emmeans() --- refuse by name, not by namespace.
#
# emmeans is a suggested package because most of flexyBayes does not need
# it. A bare requireNamespace() error carries no reason code and cannot
# be caught by class, so the guard raises the registered refusal and
# names what the fit already answers without it.
#
# @noRd
# @keywords internal
.fb_check_classify_emmeans <- function() {
  if (.fb_emmeans_available()) {
    return(invisible(TRUE))
  }
  stop(.fb_refusal_condition(
    reason_code = "classify_requires_emmeans",
    message = paste0(
      "predict(classify = ) builds the marginal-means table through ",
      "emmeans, which is not installed. emmeans is a suggested package ",
      "rather than a hard dependency because most of flexyBayes does not ",
      "need it. Install it with install.packages(\"emmeans\"), or read the ",
      "fixed effects directly from summary(fit)$fixed and the ",
      "random-effect predictions from coef(fit, what = \"random\")."
    ),
    family_class = "flexybayes_classify_requires_emmeans"
  ))
}

# .fb_check_classify_random_factor() --- refuse a marginal mean the
# reference grid cannot hold.
#
# emmeans builds its reference grid from the population-level design, so
# a factor entering the model only as a random-effects grouping term is
# not in it and the request dies inside emmeans with "No variable named
# <f> in the reference grid" -- a third-party message carrying no reason
# code, no mention of this package, and no route onward. Asking a MET fit
# for genotype means is the commonest thing a breeder does after fitting,
# and ASReml serves it, so this is a gap worth naming rather than a
# malformed call.
#
# Only a factor the fit knows as a grouping term is refused here. A name
# the model does not carry at all is left to emmeans, whose message for
# that case is already about the right thing.
#
# @noRd
# @keywords internal
.fb_check_classify_random_factor <- function(object, vars) {
  grid_vars <- tryCatch(
    all.vars(.fb_fixef_terms(object)),
    error = function(e) character(0)
  )
  missing_vars <- setdiff(vars, grid_vars)
  if (length(missing_vars) == 0L) {
    return(invisible(TRUE))
  }

  groups <- tryCatch(
    names(.fb_summary_random(object)) %||% character(0),
    error = function(e) character(0)
  )
  offending <- intersect(missing_vars, groups)
  if (length(offending) == 0L) {
    return(invisible(TRUE))
  }

  entry <- .lookup_refusal("classify_random_factor_not_supported")
  stop(.fb_refusal_condition(
    reason_code = "classify_random_factor_not_supported",
    message = sprintf(entry$message_template, offending[[1L]]),
    factor = offending[[1L]]
  ))
}

# .fb_classify_newdata_note() --- the two prediction paths do not mix.
#
# `newdata` predicts the rows a caller hands over; `classify` averages
# the fitted model over the levels it saw. Supplying both is a request
# for two different things, and only one of them is answered, so the
# other is not passed over in silence.
#
# @noRd
# @keywords internal
.fb_classify_newdata_note <- function(newdata) {
  if (is.null(newdata)) {
    return(invisible(NULL))
  }
  warning(
    "predict(): `newdata` is ignored on the classify path. The means ",
    "table averages the fitted model over the levels it was fitted to. ",
    "To predict rows you supply -- untested lines, a new site -- call ",
    "predict(fit, newdata = ..., allow_new_levels = \"sample\") without ",
    "`classify`.",
    call. = FALSE
  )
  invisible(NULL)
}

# .fb_check_level() --- the credible level is a proportion.
#
# @noRd
# @keywords internal
.fb_check_level <- function(level) {
  ok <- is.numeric(level) && length(level) == 1L && !is.na(level) &&
    level > 0 && level < 1
  if (!ok) {
    stop(
      "predict(): `level` is the credible level as a proportion, strictly ",
      "between 0 and 1 -- 0.95 for a 95% interval.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# ---------------------------------------------------------------- #
# Printing                                                          #
# ---------------------------------------------------------------- #

#' Print a classify means table
#'
#' Prints the table with a banner naming what the numbers are. The
#' estimate is a posterior mean of a marginal mean, and the table is not
#' a standard-error-of-difference table -- ASReml's `predict()` prints
#' one beside the means and this does not, so the banner says so rather
#' than leaving the absence to be noticed.
#'
#' On an INLA fit the banner carries a second line. The interval there
#' comes from the Gaussian approximation of the joint fixed-effect
#' posterior, which is what the estimation basis supplies to emmeans, and
#' not from INLA's own marginal densities.
#'
#' @param x An `fb_predict_classify` table, as returned by
#'   `predict(fit, classify = )`.
#' @param digits Number of significant digits for the printed table.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the table it prints.
#' @export
print.fb_predict_classify <- function(x, digits = 4L, ...) {
  level <- attr(x, "level") %||% 0.95
  cat(sprintf(
    "Predicted means: %s  (%g%% credible interval)\n",
    paste(attr(x, "classify") %||% character(0), collapse = " x "),
    100 * level
  ))
  cat(
    "Estimate = posterior mean of the marginal mean. Not an SED table.\n"
  )
  if (identical(attr(x, "engine"), "inla")) {
    cat(
      "INLA: intervals from the Gaussian approximation of the fixed ",
      "effects.\n",
      sep = ""
    )
  }
  body <- as.data.frame(x)
  num <- vapply(body, is.numeric, logical(1L))
  body[num] <- lapply(body[num], round, digits = digits)
  print(body, row.names = FALSE)
  invisible(x)
}
