# brms_two_part.R: reading a two-part non-linear brms fit under the names of
# the linear model.
#
# A fit refitted through fb_refit(formula_hook = ) may be the brms non-linear
# model `y ~ a + b`, with part `a` holding the intercept, the fixed effects
# and the random terms of the linear model and part `b` a block of further
# fixed effects. brms names the draws of the parts `b_a_<coef>` and
# `b_b_<coef>`; coef(), vcov(), confint(), the draws and the prior summary
# read them under the names a linear fit of the same model carries.

# The two parts: the linear model's own terms, and the added block.
.FB_NL_PROTECTED_PART <- "a"
.FB_NL_BLOCK_PART <- "b"


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
