# inla_covariate_standardisation.R, INLA fits on standardised covariates.
#
# A Gaussian identity-link fit runs on a standardised response
# (R/inla_standardise.R), which makes it one model whatever units the
# response is recorded in. Nothing did the same for the covariates, and on a
# fit that carries no response standardisation the covariate is where the
# remaining scale sensitivity sits. Measured at the base of this unit, 30
# groups of 8 with a standard-normal covariate: multiplying the covariate by
# a thousand moved the Poisson coefficient by 1.8 per cent and the Bernoulli
# coefficient by 1.7 per cent, against a run-to-run spread of 2e-7 over five
# identical fits.
#
# Each eligible numeric covariate is divided by its standard deviation before
# the fit, and every reading of it is multiplied back afterwards. Two
# properties make that exact.
#
#   * It is a division only. Centring would move the intercept by a linear
#     combination of the slopes, and INLA reports marginals instead of the
#     joint posterior such a combination would need, so the intercept could
#     not be mapped back exactly. Dividing maps every quantity by a scalar of
#     its own.
#   * The linear predictor is unchanged row by row: beta (x / d) with
#     beta* = beta d is the same number as beta x. The likelihood, the fitted
#     values, the marginal likelihood, DIC, WAIC and CPO are therefore the
#     same on both coordinate systems, and only the fixed coefficients, the
#     uncorrelated slope field and that field's precision move.
#
# What the standardisation buys is the prior. A quantity with no stated prior
# takes INLA's own default on the standardised coordinates, which is written
# for a covariate of unit spread. And a prior the package or the caller did
# state is restated on those coordinates, so the model is the one that was
# asked for whatever units the covariate is recorded in.


# --- what is standardised ------------------------------------------------- #

#' The numeric covariates an INLA fit is standardised on
#'
#' @param fb An `fb_terms` object.
#' @param data The data the fit is built from.
#' @returns `NULL` when nothing is eligible, otherwise a list with
#'   `divisor`, the standard deviation of each standardised data column.
#'   `coef_divisor`, the factor each fixed-effect coefficient is multiplied
#'   by on the standardised coordinates, keyed by design column.
#'   `slope_divisor`, the same for each uncorrelated slope field, keyed by
#'   the INLA index the field is written over. And `slope_key_divisor`, keyed
#'   by the prior DSL's name for that field's standard deviation.
#' @noRd
#' @keywords internal
.fb_inla_covariate_standardisation <- function(fb, data) {
  if (!is.data.frame(data)) {
    return(NULL)
  }
  columns <- .fb_inla_covariate_columns(fb, data)
  if (!length(columns)) {
    return(NULL)
  }

  # The rows a default is learned on: those whose working response is
  # observed, as every other covariate-scaled quantity reads them
  # (R/inla_default_priors.R).
  rows <- .fb_inla_default_rows(fb, data)
  divisor <- numeric(0)
  for (nm in columns) {
    x <- suppressWarnings(as.numeric(data[[nm]][rows]))
    x <- x[is.finite(x)]
    s <- if (length(x) >= 2L) stats::sd(x) else NA_real_
    # A covariate with no spread has nothing to divide by, and dividing by
    # one is the fit the package already made.
    if (is.finite(s) && s > 0) {
      divisor[[nm]] <- s
    }
  }
  if (!length(divisor)) {
    return(NULL)
  }

  slope <- .fb_inla_covariate_slope_divisor(fb, divisor)
  list(
    divisor = divisor,
    coef_divisor = .fb_inla_covariate_coef_divisor(fb, data, divisor),
    slope_divisor = slope$field,
    slope_key_divisor = slope$key
  )
}

#' The data columns eligible for standardisation
#'
#' A column is eligible when it enters the model only in ways a division
#' passes through exactly: as a continuous fixed term, as the numeric half of
#' a factor-by-numeric interaction, as a part of a numeric interaction, or as
#' the covariate of an uncorrelated random slope.
#'
#' Everything else protects the column instead of standardising it. A
#' formula transform reads the column through a function, and `log(x / d)` is
#' not `log(x)` divided by anything, so the coefficient could not be mapped
#' back. The response, an offset and an addition term are not covariates. A
#' grouping factor indexes a field. It does not weight one. The correlated
#' slope's covariate scales a Wishart prior built from `var(y) / var(x)`
#' (R/emit_inla.R), which a division would state twice. A protected column is
#' left exactly as the caller recorded it.
#'
#' @param fb An `fb_terms` object.
#' @param data The data the fit is built from.
#' @returns A character vector of column names.
#' @noRd
#' @keywords internal
.fb_inla_covariate_columns <- function(fb, data) {
  eligible <- character(0)
  for (term in fb$fixed_terms %||% list()) {
    type <- term$type %||% ""
    if (identical(type, "continuous")) {
      eligible <- c(eligible, term$var %||% character(0))
    } else if (identical(type, "factor_numeric_interaction")) {
      eligible <- c(eligible, term$continuous %||% character(0))
    } else if (identical(type, "interaction")) {
      eligible <- c(eligible, term$vars %||% character(0))
    }
  }
  for (term in fb$random_terms %||% list()) {
    if (identical(term$type %||% "", "simple_slope_uncor")) {
      eligible <- c(eligible, term$slope_var %||% character(0))
    }
  }
  eligible <- unique(eligible[nzchar(eligible)])
  eligible <- eligible[eligible %in% names(data)]
  if (!length(eligible)) {
    return(character(0))
  }
  numeric_col <- vapply(
    eligible,
    function(nm) is.numeric(data[[nm]]) && !is.matrix(data[[nm]]),
    logical(1)
  )
  setdiff(eligible[numeric_col], .fb_inla_covariate_protected(fb))
}

# Every variable name an expression reads, or nothing when it does not parse.
.fb_inla_covariate_expression_vars <- function(text) {
  out <- character(0)
  for (s in text[nzchar(text %||% character(0))]) {
    out <- c(out, tryCatch(all.vars(str2lang(s)), error = function(e) {
      character(0)
    }))
  }
  out
}

#' Columns no standardisation may touch
#'
#' @param fb An `fb_terms` object.
#' @returns A character vector of column names.
#' @noRd
#' @keywords internal
.fb_inla_covariate_protected <- function(fb) {
  protected <- .fb_inla_covariate_expression_vars(c(
    fb$response %||% character(0),
    fb$offsets %||% character(0)
  ))
  for (term in fb$fixed_terms %||% list()) {
    if (identical(term$type %||% "", "expression")) {
      protected <- c(
        protected,
        .fb_inla_covariate_expression_vars(term$label %||% character(0))
      )
    }
  }
  for (term in fb$random_terms %||% list()) {
    # [[ ]] because $var would partially match a combo term's `vars`.
    protected <- c(
      protected,
      term[["var"]] %||% character(0),
      term$vars %||% character(0),
      term$outer %||% character(0),
      term$inner %||% character(0)
    )
    if (!identical(term$type %||% "", "simple_slope_uncor")) {
      protected <- c(protected, term$slope_var %||% character(0))
    }
  }
  for (term in fb$addition_terms %||% list()) {
    protected <- c(
      protected,
      .fb_inla_covariate_expression_vars(term$label %||% character(0)),
      term$var %||% character(0)
    )
  }
  unique(protected)
}

#' The factor each fixed-effect coefficient carries on the standardised
#' coordinates
#'
#' A design column is the product of the variables its term names, so its
#' coefficient is multiplied by the product of their divisors. A column whose
#' term names no standardised variable is absent, which the caller reads as a
#' factor of one.
#'
#' @param fb An `fb_terms` object.
#' @param data The data the fit is built from.
#' @param divisor The per-column divisors.
#' @returns A named numeric vector keyed by design column.
#' @noRd
#' @keywords internal
.fb_inla_covariate_coef_divisor <- function(fb, data, divisor) {
  form <- .fb_inla_fixed_design_formula(fb)
  if (is.null(form)) {
    return(numeric(0))
  }
  mf <- tryCatch(stats::model.frame(form, data = data), error = function(e) {
    NULL
  })
  if (is.null(mf)) {
    return(numeric(0))
  }
  tt <- attr(mf, "terms")
  mm <- stats::model.matrix(tt, mf)
  factors <- attr(tt, "factors")
  assign <- attr(mm, "assign")

  out <- numeric(0)
  for (j in seq_len(ncol(mm))) {
    a <- assign[[j]]
    if (a == 0L) {
      next
    }
    vars <- rownames(factors)[factors[, a] > 0]
    hit <- vars[vars %in% names(divisor)]
    if (!length(hit)) {
      next
    }
    out[[colnames(mm)[[j]]]] <- prod(divisor[hit])
  }
  out
}

#' The factor each uncorrelated slope field carries
#'
#' @param fb An `fb_terms` object.
#' @param divisor The per-column divisors.
#' @returns A list with `field`, keyed by the INLA index the field is written
#'   over, and `key`, keyed by the prior DSL's name for its standard
#'   deviation.
#' @noRd
#' @keywords internal
.fb_inla_covariate_slope_divisor <- function(fb, divisor) {
  field <- numeric(0)
  key <- numeric(0)
  for (term in fb$random_terms %||% list()) {
    if (!identical(term$type %||% "", "simple_slope_uncor")) {
      next
    }
    sv <- term$slope_var %||% ""
    if (!nzchar(sv) || !sv %in% names(divisor)) {
      next
    }
    field[[paste0(term$var, "_for_slope_", sv)]] <- divisor[[sv]]
    key[[paste0(sv, "_", term$var)]] <- divisor[[sv]]
  }
  list(field = field, key = key)
}


# --- stated priors on the standardised coordinates ------------------------ #

#' Every uncorrelated slope prior restated on the standardised covariate
#'
#' The standard deviation of an uncorrelated slope is in response units per
#' covariate unit, so a covariate divided by `d` multiplies it by `d`. The
#' map is `.fb_inla_map_hyper_entry()`, which divides, so it is handed
#' `1 / d`. Every prior family that reaches it is closed under that map, and
#' one that is not refuses there, so no misstated prior reaches INLA.
#'
#' @param hyper_ctrl The `priors_to_inla()` list.
#' @param cstd The covariate standardisation record.
#' @returns `hyper_ctrl` with each slope entry restated.
#' @noRd
#' @keywords internal
.fb_inla_covariate_standardise_hyper <- function(hyper_ctrl, cstd) {
  keys <- intersect(names(hyper_ctrl), names(cstd$slope_key_divisor))
  for (k in keys) {
    hyper_ctrl[[k]] <- .fb_inla_map_hyper_entry(
      hyper_ctrl[[k]],
      k = 1 / cstd$slope_key_divisor[[k]]
    )
  }
  hyper_ctrl
}

#' The stated fixed-effect priors restated on the standardised coordinates
#'
#' A coefficient multiplied by `d` has its prior mean multiplied by `d` and
#' its precision divided by `d^2`. Only a coefficient the caller or the
#' package gave a prior is restated. One with neither keeps INLA's own
#' default on the standardised coordinates. The intercept is untouched: the
#' standardisation does not centre.
#'
#' @param control_fixed The `control.fixed` list the emit built.
#' @param cstd The covariate standardisation record.
#' @returns The restated `control.fixed` list.
#' @noRd
#' @keywords internal
.fb_inla_covariate_standardise_control_fixed <- function(control_fixed, cstd) {
  d <- cstd$coef_divisor
  if (!length(d) || !length(control_fixed)) {
    return(control_fixed)
  }
  out <- control_fixed
  out$mean <- .fb_inla_covariate_fixed_slot(
    control_fixed$mean, d, function(v, k) v * k
  )
  out$prec <- .fb_inla_covariate_fixed_slot(
    control_fixed$prec, d, function(v, k) v / k^2
  )
  out
}

# One control.fixed slot restated per coefficient. A scalar slot is a
# statement about every coefficient at once, and the standardised ones no
# longer share a value, so it is expanded into the named form INLA reads,
# with `default` carrying the scalar for the rest.
.fb_inla_covariate_fixed_slot <- function(slot, d, f) {
  if (is.null(slot)) {
    return(NULL)
  }
  if (!is.list(slot)) {
    entries <- lapply(names(d), function(nm) f(slot, d[[nm]]))
    names(entries) <- names(d)
    return(c(entries, list(default = slot)))
  }
  fallback <- slot[["default"]]
  for (nm in names(d)) {
    if (!is.null(slot[[nm]])) {
      slot[[nm]] <- f(slot[[nm]], d[[nm]])
    } else if (!is.null(fallback)) {
      slot[[nm]] <- f(fallback, d[[nm]])
    }
  }
  slot
}


# --- the result on the original scale ------------------------------------- #

# The factor a hyperparameter is multiplied by. A slope field divided by `d`
# has its standard deviation divided by `d`, so its precision multiplies by
# `d^2`. The residual, every other field and every correlation are unchanged,
# because the covariate carries no information about them.
.fb_inla_covariate_hyper_factor <- function(name, cstd) {
  if (!startsWith(name, "Precision for ")) {
    return(1)
  }
  field <- sub("^Precision for ", "", name)
  if (!field %in% names(cstd$slope_divisor)) {
    return(1)
  }
  cstd$slope_divisor[[field]]^2
}

#' An INLA result fitted on standardised covariates, on the original scale
#'
#' The fixed coefficients, the uncorrelated slope fields and those fields'
#' precisions, with their marginals. Nothing else is touched: the linear
#' predictor is the same number on both coordinate systems, so the fitted
#' values, the marginal likelihood, DIC, WAIC and CPO need no correction, and
#' `$.args`, `$mode` and `$misc` keep the engine's own coordinates, since they
#' are what INLA ran and what its samplers read.
#'
#' @param fit An `inla` result.
#' @param cstd The covariate standardisation the fit was run under.
#' @returns The result on the original scale.
#' @noRd
#' @keywords internal
.fb_inla_covariate_back_transform <- function(fit, cstd) {
  d <- cstd$coef_divisor
  s <- cstd$slope_divisor

  # Fixed effects -------------------------------------------------------

  coef_b <- function(nm) if (nm %in% names(d)) 1 / d[[nm]] else 1
  fit$summary.fixed <- .fb_inla_map_rows(
    fit$summary.fixed, function(nm) 0, coef_b
  )
  for (nm in names(fit$marginals.fixed)) {
    fit$marginals.fixed[[nm]] <- .fb_inla_map_marginal(
      fit$marginals.fixed[[nm]], 0, coef_b(nm)
    )
  }

  # The uncorrelated slope fields ---------------------------------------

  for (nm in intersect(names(fit$summary.random), names(s))) {
    fit$summary.random[[nm]] <- .fb_inla_map_table(
      fit$summary.random[[nm]], 0, 1 / s[[nm]]
    )
    if (!is.null(fit$marginals.random[[nm]])) {
      fit$marginals.random[[nm]] <- lapply(
        fit$marginals.random[[nm]],
        .fb_inla_map_marginal,
        a = 0,
        b = 1 / s[[nm]]
      )
    }
  }

  # Their precisions ----------------------------------------------------

  hyper_b <- function(nm) .fb_inla_covariate_hyper_factor(nm, cstd)
  fit$summary.hyperpar <- .fb_inla_map_rows(
    fit$summary.hyperpar, function(nm) 0, hyper_b
  )
  for (nm in names(fit$marginals.hyperpar)) {
    fit$marginals.hyperpar[[nm]] <- .fb_inla_map_marginal(
      fit$marginals.hyperpar[[nm]], 0, hyper_b(nm)
    )
  }
  # The internal hyperparameters are log precisions, which shift.
  internal_a <- function(nm) {
    if (!startsWith(nm, "Log precision for ")) {
      return(0)
    }
    log(.fb_inla_covariate_hyper_factor(
      sub("^Log precision for ", "Precision for ", nm), cstd
    ))
  }
  fit$internal.summary.hyperpar <- .fb_inla_map_rows(
    fit$internal.summary.hyperpar, internal_a, function(nm) 1
  )
  for (nm in names(fit$internal.marginals.hyperpar)) {
    fit$internal.marginals.hyperpar[[nm]] <- .fb_inla_map_marginal(
      fit$internal.marginals.hyperpar[[nm]], internal_a(nm), 1
    )
  }

  fit
}

# The factor each latent draw row is multiplied by, in the order INLA writes
# them. A row the standardisation did not touch takes one.
.fb_inla_covariate_draw_factor <- function(prefixes, cstd) {
  b <- rep(1, length(prefixes))
  d <- c(cstd$coef_divisor, cstd$slope_divisor)
  for (nm in intersect(unique(prefixes), names(d))) {
    b[prefixes == nm] <- 1 / d[[nm]]
  }
  b
}
