# penalty_block_nl.R -- the reference model of a `scope = "block"` penalty.
#
# A class-wide prior (horseshoe, R2D2) with `scope = "block"` is fitted as the
# brms non-linear formula
#
#   y ~ a + b,   a ~ <intercept, protected fixed effects, random terms>,
#                b ~ 0 + <candidate block>,
#
# so the class-wide prior sits on part `b` alone and the protected effects
# keep the base model's own priors in part `a`. brms names the draws of the
# two parts `b_a_<coef>` and `b_b_<coef>`; every reader in this package sees
# them under the names a linear fit of the same model would carry.


#' Whether a penalty asks for a block-scoped class-wide prior
#'
#' @param penalty An `fb_penalty`, or `NULL`.
#' @returns `TRUE` or `FALSE`.
#' @noRd
#' @keywords internal
.fb_pen_block_scope <- function(penalty) {
  if (is.null(penalty)) {
    return(FALSE)
  }
  any(vapply(
    penalty$specs,
    function(spec) identical(spec$args$scope, "block"),
    logical(1L)
  ))
}


#' Write a linear brms formula as the two-part non-linear formula
#'
#' @param form The linear brms formula, plain or a `brmsformula` carrying a
#'   `sigma` sub-formula.
#' @param nl_block The compiler's `nl_block` record, with `candidate_terms`
#'   added by the reference fit.
#' @returns A `brmsformula` with `nl = TRUE`.
#' @noRd
#' @keywords internal
.fb_nl_block_formula <- function(form, nl_block) {
  main <- .brms_main_formula(form)
  tt <- stats::terms(main)
  labels <- attr(tt, "term.labels")
  candidates <- nl_block$candidate_terms

  absent <- setdiff(candidates, labels)
  if (length(absent) > 0L) {
    stop(
      ".fb_nl_block_formula(): the candidate terms ",
      paste(absent, collapse = ", "),
      " are not in the reference formula.",
      call. = FALSE
    )
  }
  .fb_nl_block_check_names(c(all.vars(main), labels))

  protected <- setdiff(labels, candidates)
  protected <- ifelse(
    grepl("|", protected, fixed = TRUE),
    paste0("(", protected, ")"),
    protected
  )
  offsets <- attr(tt, "offset")
  if (length(offsets) > 0L) {
    protected <- c(
      protected,
      vapply(
        offsets,
        function(i) deparse(attr(tt, "variables")[[i + 1L]]),
        character(1L)
      )
    )
  }
  intercept <- if (attr(tt, "intercept") == 1L) "1" else "0"
  lhs <- paste(deparse(main[[2L]]), collapse = " ")
  a_part <- nl_block$protected_part
  b_part <- nl_block$block_part

  parts <- list(
    stats::as.formula(paste0(lhs, " ~ ", a_part, " + ", b_part)),
    stats::as.formula(paste0(
      a_part, " ~ ", paste(c(intercept, protected), collapse = " + ")
    )),
    stats::as.formula(paste0(
      b_part, " ~ 0 + ", paste(candidates, collapse = " + ")
    ))
  )
  sigma_form <- if (inherits(form, "brmsformula")) {
    form$pforms$sigma
  } else {
    NULL
  }
  if (!is.null(sigma_form)) {
    parts <- c(parts, list(sigma_form))
  }
  do.call(brms::bf, c(parts, list(nl = TRUE)))
}

# A variable or term named as one of the two parts would be read as the part.
.fb_nl_block_check_names <- function(names_used) {
  taken <- intersect(
    c(.FB_NL_PROTECTED_PART, .FB_NL_BLOCK_PART),
    names_used
  )
  if (length(taken) == 0L) {
    return(invisible(NULL))
  }
  stop(.fb_refusal_condition(
    reason_code = "penalty_block_part_name_taken",
    message = paste0(
      "fb_penalty(): scope = \"block\" fits the model as the two non-linear ",
      "parts `", .FB_NL_PROTECTED_PART, "` and `", .FB_NL_BLOCK_PART,
      "`, and the model uses ",
      paste0("`", taken, "`", collapse = " and "),
      " as a variable name. Rename the column."
    ),
    names = taken
  ))
}


#' Move the protected prior rows onto the protected part
#'
#' A row on the mean model's coefficients, intercept or group-level standard
#' deviations goes to part `a`; brms has no `Intercept` class in a
#' non-linear part, so the intercept row becomes the `b` row of coefficient
#' `Intercept`. Rows of a distributional parameter and of `sigma` stay.
#'
#' @param prior A `brmsprior` data frame, or `NULL`.
#' @param nl_block The compiler's `nl_block` record.
#' @returns `prior` with its rows moved, attributes kept.
#' @noRd
#' @keywords internal
.fb_nl_block_prior <- function(prior, nl_block) {
  if (is.null(prior) || NROW(prior) == 0L) {
    return(prior)
  }
  dpar <- as.character(prior$dpar %||% rep("", nrow(prior)))
  nlpar <- as.character(prior$nlpar %||% rep("", nrow(prior)))
  dpar[is.na(dpar)] <- ""
  nlpar[is.na(nlpar)] <- ""
  mean_row <- !nzchar(dpar) & !nzchar(nlpar)

  intercept <- mean_row & prior$class == "Intercept"
  prior$class[intercept] <- "b"
  prior$coef[intercept] <- "Intercept"

  moved <- mean_row & prior$class %in% c("b", "sd")
  prior$nlpar[moved] <- nl_block$protected_part
  prior
}


#' The parts of a package-built two-part non-linear brms fit, or NULL
#'
#' @param brmsfit A `brmsfit`.
#' @returns The two part names, or `NULL` for any other fit.
#' @noRd
#' @keywords internal
.fb_brms_nl_parts <- function(brmsfit) {
  form <- brmsfit$formula
  if (!inherits(form, "brmsformula") || !isTRUE(attr(form$formula, "nl"))) {
    return(NULL)
  }
  parts <- c(.FB_NL_PROTECTED_PART, .FB_NL_BLOCK_PART)
  if (!all(parts %in% names(form$pforms))) {
    return(NULL)
  }
  parts
}


#' Draw names of a two-part fit under the names of the linear model
#'
#' `b_a_regionR2` and `b_b_x1` become `b_regionR2` and `b_x1`;
#' `sdb_b_x1` becomes `sdb_x1`; `hs_global_b` becomes `hs_global`;
#' `sd_g__a_Intercept` becomes `sd_g__Intercept`; `r_g__a[l,Intercept]`
#' becomes `r_g[l,Intercept]`.
#'
#' @param nm Draw names.
#' @param parts The two part names.
#' @returns `nm`, renamed.
#' @noRd
#' @keywords internal
.fb_nl_canonical_names <- function(nm, parts) {
  alt <- paste(parts, collapse = "|")
  nm <- sub(paste0("^b_(", alt, ")_"), "b_", nm)
  nm <- sub(paste0("^sdb_(", alt, ")_"), "sdb_", nm)
  nm <- sub(paste0("^(hs_[a-z]+|R2D2_[A-Za-z0-9]+)_(", alt, ")$"), "\\1", nm)
  sd_cor <- grepl("^(sd|cor)_", nm)
  nm[sd_cor] <- gsub(paste0("__(", alt, ")_"), "__", nm[sd_cor])
  nm <- sub(paste0("^r_(.+?)__(", alt, ")\\["), "r_\\1[", nm)
  nm
}


#' Posterior draws of a brms fit as a matrix, under linear-model names
#'
#' @param brmsfit A `brmsfit`.
#' @returns A numeric matrix of draws.
#' @noRd
#' @keywords internal
.fb_brms_draws_matrix <- function(brmsfit) {
  draws <- as.matrix(posterior::as_draws_matrix(brmsfit))
  parts <- .fb_brms_nl_parts(brmsfit)
  if (!is.null(parts)) {
    colnames(draws) <- .fb_nl_canonical_names(colnames(draws), parts)
  }
  draws
}
