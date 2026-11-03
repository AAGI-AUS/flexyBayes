# ordinal.R -- ordinal responses: the cumulative and continuation-ratio
# families.
#
# Both families write the K - 1 thresholds tau_1, ..., tau_{K-1} in brms's
# convention. The cumulative family models P(Y <= k) = F(tau_k - eta) and the
# continuation-ratio family P(Y > k | Y >= k) = F(eta - tau_k), so a larger
# linear predictor eta moves an observation to a higher category in both. On
# INLA the cumulative family is the likelihood "pom" and the
# continuation-ratio family a binomial fit to the rows of its binary
# expansion. On brms both are brms's own families.
#
# This file holds the response record, the engine set-up on INLA, the
# threshold readers, the category probabilities and the refusal codes. The
# methods are in R/ordinal_methods.R and exact K-fold in R/fb_kfold.R.

# The two ordinal family spellings flexybayes() admits.
.FB_ORDINAL_FAMILIES <- c("cumulative", "cratio")

# The engine-side column names: the response coded 1..K (or its binary
# expansion), and the threshold factor of the expansion.
.FB_ORDINAL_RESPONSE <- "fb_ordinal_y"
.FB_ORDINAL_CUT <- "fb_ordinal_cut"

# Number of posterior draws used to summarise the cumulative family's
# thresholds on INLA, which INLA reports only as hyperparameters on its own
# internal scale.
.FB_ORDINAL_THRESHOLD_DRAWS <- 4000L

# The seed of that summary, so that print() and summary() show the same
# numbers on every call. The session's random stream is restored afterwards.
.FB_ORDINAL_THRESHOLD_SEED <- 20261103L

# The concentration of the Dirichlet prior INLA's "pom" likelihood puts on the
# category probabilities at a linear predictor of zero. INLA's own default is
# 3. At 1, the uniform distribution on those probabilities, the thresholds of
# the sealed recovery design (20 seeds, model with the field) are recovered
# with less bias at the same coverage, and the model without the field agrees
# with brms within 0.16 brms posterior SD, against 0.42 at 3.
.FB_ORDINAL_POM_DIRICHLET <- 1

# Nodes of the Gauss-Hermite rule that integrates a category probability over
# the random effects of a new unit.
.FB_ORDINAL_GH_NODES <- 40L


# Family predicates ------------------------------------------------------------

#' The ordinal family a family argument names, or NULL
#'
#' @param family A family name, a family object or `NULL`.
#' @returns `"cumulative"`, `"cratio"` or `NULL`.
#' @noRd
#' @keywords internal
.fb_ordinal_family <- function(family) {
  fam <- if (inherits(family, "family")) family$family else family
  if (!is.character(fam) || length(fam) != 1L) {
    return(NULL)
  }
  fam <- tolower(fam)
  if (fam %in% .FB_ORDINAL_FAMILIES) fam else NULL
}

#' Whether an IR or a fit has an ordinal response
#'
#' @param x An `fb_terms` IR or a flexyBayes fit.
#' @returns `TRUE` or `FALSE`.
#' @noRd
#' @keywords internal
.fb_is_ordinal <- function(x) {
  fb <- if (inherits(x, "flexybayes")) .fb_fit_ir(x) else x
  !is.null(.fb_ordinal_family(fb$family))
}

#' The largest number of categories INLA's cumulative likelihood takes
#'
#' INLA's "pom" likelihood declares one hyperparameter per threshold, ten at
#' INLA 26.8.7, so it fits at most eleven categories. Read from the installed
#' engine's own declaration when INLA is present.
#'
#' @returns An integer.
#' @noRd
#' @keywords internal
.fb_ordinal_inla_max_categories <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    return(11L)
  }
  length(INLA::inla.models()$likelihood$pom$hyper) + 1L
}


# Response record --------------------------------------------------------------

#' The categories of an ordinal response
#'
#' The categories are the values the data hold, in order: the levels of an
#' ordered factor that occur, or the distinct values of an integer-valued
#' numeric response. A level of an ordered factor that no row holds is
#' dropped, with a message, because no observation can inform its
#' threshold.
#'
#' @param fb An `fb_terms` IR with an ordinal family.
#' @param data The data the IR was built from.
#' @returns A list: `family`, `levels` (category labels), `values` (the
#'   response values, character or numeric, in order), `K`, `dropped` and
#'   `type` (`"ordered"` or `"integer"`).
#' @noRd
#' @keywords internal
.fb_ordinal_record <- function(fb, data) {
  fam <- .fb_ordinal_family(fb$family)
  y <- .fb_response_values(fb, data)
  if (is.null(y)) {
    stop(.fb_refusal_condition(
      reason_code = "response_not_in_data",
      message = paste0(
        "Response variable \"",
        fb$response,
        "\" not found in data."
      )
    ))
  }

  # Ordered factor or integer-valued numeric ------------------------------
  dropped <- character(0)
  if (is.ordered(y)) {
    seen <- as.character(unique(y[!is.na(y)]))
    values <- levels(y)[levels(y) %in% seen]
    dropped <- setdiff(levels(y), values)
    type <- "ordered"
  } else if (is.numeric(y)) {
    finite <- y[is.finite(y)]
    off <- which(finite != round(finite))
    if (length(off)) {
      stop(.fb_refusal_condition(
        reason_code = "ordinal_response_not_integer",
        message = paste0(
          "family = \"",
          fam,
          "\" models an ordered category, and the response `",
          fb$response,
          "` holds values that are not whole numbers (",
          paste(utils::head(format(finite[off]), 3L), collapse = ", "),
          "). Code the categories as whole numbers in their order, or ",
          "pass an ordered factor, factor(y, levels = , ordered = TRUE)."
        ),
        column = fb$response
      ))
    }
    values <- sort(unique(finite))
    type <- "integer"
  } else {
    stop(.fb_refusal_condition(
      reason_code = "ordinal_response_unordered",
      message = paste0(
        "family = \"",
        fam,
        "\" models an ordered category, and the response `",
        fb$response,
        "` is a ",
        class(y)[[1L]],
        ", which does not state an order. Pass an ordered factor, ",
        "factor(y, levels = <categories, lowest first>, ordered = TRUE), ",
        "or whole numbers in category order."
      ),
      column = fb$response
    ))
  }

  # At least two categories -----------------------------------------------
  if (length(values) < 2L) {
    stop(.fb_refusal_condition(
      reason_code = "ordinal_response_too_few_categories",
      message = paste0(
        "family = \"",
        fam,
        "\" needs at least two observed categories, and the response `",
        fb$response,
        "` holds ",
        if (length(values)) {
          paste0("one, ", format(values[[1L]]), ".")
        } else {
          "none."
        }
      ),
      column = fb$response
    ))
  }
  if (length(dropped)) {
    message(
      "flexyBayes: the ordered response `",
      fb$response,
      "` has no observation at ",
      if (length(dropped) == 1L) "level " else "levels ",
      paste(dropped, collapse = ", "),
      ", so the model has ",
      length(values),
      " categories, the levels observed."
    )
  }
  list(
    family = fam,
    levels = as.character(values),
    values = values,
    K = length(values),
    dropped = dropped,
    type = type
  )
}

#' The response coded 1..K
#'
#' @param fb An `fb_terms` IR carrying `$ordinal`.
#' @param data A data frame holding the response.
#' @param record The ordinal record, `fb$ordinal` by default.
#' @returns An integer vector, `NA` where the response is missing.
#' @noRd
#' @keywords internal
.fb_ordinal_codes <- function(fb, data, record = fb$ordinal) {
  y <- .fb_response_values(fb, data)
  key <- if (is.factor(y)) as.character(y) else y
  codes <- match(key, record$values)
  unknown <- which(!is.na(y) & is.na(codes))
  if (length(unknown)) {
    stop(.fb_refusal_condition(
      reason_code = "ordinal_response_category_unknown",
      message = paste0(
        "The response `",
        fb$response,
        "` holds ",
        if (length(unknown) == 1L) "a value" else "values",
        " outside the model's ",
        record$K,
        " categories (",
        paste(utils::head(unique(as.character(y[unknown])), 3L),
          collapse = ", "),
        "). The categories are ",
        paste(record$levels, collapse = ", "),
        "."
      ),
      column = fb$response
    ))
  }
  as.integer(codes)
}

#' Refuse an ordinal fixed formula without an intercept
#'
#' The thresholds of an ordinal model take the place of the intercept, and
#' the emits rely on the intercept to give every factor its treatment
#' contrasts.
#'
#' @param fb An `fb_terms` IR.
#' @returns Invisibly `TRUE`.
#' @noRd
#' @keywords internal
.fb_ordinal_refuse_unsupported <- function(fb) {
  fam <- .fb_ordinal_family(fb$family)
  if (is.null(fam)) {
    return(invisible(TRUE))
  }
  if (!isTRUE(fb$intercept)) {
    stop(.fb_refusal_condition(
      reason_code = "ordinal_fixed_without_intercept",
      message = paste0(
        "family = \"",
        fam,
        "\" estimates one threshold between each pair of adjacent ",
        "categories, and the thresholds take the place of the intercept. ",
        "Write the fixed formula with its intercept: drop `0 +` or `- 1`."
      )
    ))
  }
  if (length(fb$offsets)) {
    stop(.fb_refusal_condition(
      reason_code = "ordinal_offset_not_supported",
      message = paste0(
        "family = \"",
        fam,
        "\" is fitted without an offset in this version: the category ",
        "probabilities predict() and fitted() report do not carry one. ",
        "Enter the offset variable as a fixed covariate."
      )
    ))
  }
  invisible(TRUE)
}

#' Refuse an engine column that the data already hold
#'
#' @param data A data frame.
#' @param name The column name the engine needs.
#' @returns Invisibly `TRUE`.
#' @noRd
#' @keywords internal
.fb_ordinal_engine_column <- function(data, name) {
  if (name %in% names(data)) {
    stop(.fb_refusal_condition(
      reason_code = "nested_index_name_collision",
      message = paste0(
        "Fitting an ordinal response needs a column named `",
        name,
        "`, and `data` already has one. Rename that column."
      )
    ))
  }
  invisible(TRUE)
}


# INLA set-up ------------------------------------------------------------------

#' The IR and data an ordinal model is fitted from on INLA
#'
#' The cumulative family keeps one row per observation and writes the
#' response as its category code. The continuation-ratio family writes each
#' observation as one Bernoulli row per category boundary it reached,
#' `1{Y > k}` for k = 1, ..., min(Y, K - 1), with a threshold factor entered
#' first in an intercept-free formula, so each threshold is a fixed effect
#' and every other factor keeps its treatment contrasts. A missing response
#' is one row with a missing outcome.
#'
#' @param fb The `fb_terms` IR.
#' @param data The data, level-legalised.
#' @returns `NULL` for a family that is not ordinal, otherwise a list: `fb`,
#'   the IR the engine fits, `data`, the rows it fits, `record`, `route`
#'   (`"pom"` or `"expansion"`) and `user_rows`, the engine row that holds
#'   each user row (its first boundary row under the expansion).
#' @noRd
#' @keywords internal
.fb_ordinal_inla_setup <- function(fb, data) {
  fam <- .fb_ordinal_family(fb$family)
  if (is.null(fam)) {
    return(NULL)
  }
  record <- fb$ordinal %||% .fb_ordinal_record(fb, data)
  fb$ordinal <- record
  codes <- .fb_ordinal_codes(fb, data, record)
  K <- record$K
  .fb_ordinal_engine_column(data, .FB_ORDINAL_RESPONSE)
  engine_fb <- fb
  engine_fb$response <- .FB_ORDINAL_RESPONSE

  # Cumulative: one row per observation ------------------------------------
  if (identical(fam, "cumulative")) {
    data[[.FB_ORDINAL_RESPONSE]] <- codes
    return(list(
      fb = engine_fb,
      data = data,
      record = record,
      route = "pom",
      user_rows = seq_len(nrow(data))
    ))
  }

  # Continuation ratio: the binary expansion -------------------------------
  .fb_ordinal_engine_column(data, .FB_ORDINAL_CUT)
  reached <- ifelse(is.na(codes), 1L, pmin(codes, K - 1L))
  rows <- rep(seq_len(nrow(data)), reached)
  k <- sequence(reached)
  expanded <- data[rows, , drop = FALSE]
  rownames(expanded) <- NULL
  expanded[[.FB_ORDINAL_RESPONSE]] <- as.integer(codes[rows] > k)
  if (K > 2L) {
    expanded[[.FB_ORDINAL_CUT]] <- factor(k, levels = seq_len(K - 1L))
    cut_term <- list(
      type = "factor",
      var = .FB_ORDINAL_CUT,
      n_levels = K - 1L,
      levels = as.character(seq_len(K - 1L)),
      label = .FB_ORDINAL_CUT
    )
    engine_fb$intercept <- FALSE
    engine_fb$fixed_terms <- c(list(cut_term), fb$fixed_terms)
  }
  list(
    fb = engine_fb,
    data = expanded,
    record = record,
    route = "expansion",
    user_rows = match(seq_len(nrow(data)), rows)
  )
}

#' Fix the intercept of a cumulative INLA fit at zero
#'
#' The thresholds of INLA's "pom" likelihood absorb the intercept. Keeping
#' the intercept in the formula keeps every factor on treatment contrasts,
#' and fixing it at zero, the way an aliased column is fixed
#' (R/fixed_alias.R), leaves the model INLA documents.
#'
#' @param control_fixed The `control.fixed` list.
#' @returns The list with the intercept fixed at 0.
#' @noRd
#' @keywords internal
.fb_ordinal_pin_intercept <- function(control_fixed) {
  control_fixed$mean.intercept <- 0
  control_fixed$prec.intercept <- .FB_FIXED_ALIAS_PIN_PREC
  control_fixed
}

#' The fixed-effect rows of an INLA ordinal fit that are not coefficients
#'
#' The fixed intercept of the cumulative route and the threshold rows of
#' the expansion.
#'
#' @param route `"pom"` or `"expansion"`.
#' @param K The number of categories.
#' @returns A character vector of INLA row names.
#' @noRd
#' @keywords internal
.fb_ordinal_engine_rows <- function(route, K) {
  if (identical(route, "pom") || K == 2L) {
    return("(Intercept)")
  }
  paste0(.FB_ORDINAL_CUT, seq_len(K - 1L))
}


# Thresholds -------------------------------------------------------------------

#' Threshold names, one per pair of adjacent categories
#'
#' @param levels The category labels.
#' @returns `threshold[<lower>|<upper>]` for each boundary.
#' @noRd
#' @keywords internal
.fb_ordinal_threshold_names <- function(levels) {
  K <- length(levels)
  paste0("threshold[", levels[-K], "|", levels[-1L], "]")
}

#' Summarise threshold draws
#'
#' @param draws A matrix, draws by thresholds.
#' @param names The threshold names.
#' @returns A data frame: `term`, `estimate` (posterior mean), `std.error`,
#'   `conf.low` and `conf.high` (2.5% and 97.5% quantiles).
#' @noRd
#' @keywords internal
.fb_ordinal_summarise_draws <- function(draws, names) {
  q <- apply(draws, 2L, stats::quantile, probs = c(0.025, 0.975))
  data.frame(
    term = names,
    estimate = colMeans(draws),
    std.error = apply(draws, 2L, stats::sd),
    conf.low = q[1L, ],
    conf.high = q[2L, ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Thresholds from INLA's cumulative hyperparameters
#'
#' INLA's "pom" likelihood carries theta_1, ..., theta_{K-1} with
#' alpha_1 = theta_1 and alpha_k = alpha_{k-1} + exp(theta_k), and models
#' P(Y <= k) = F(alpha_k - eta), so alpha is tau in brms's convention.
#'
#' The engine reports the hyperparameters on two scales. Its summary table
#' and `inla.hyperpar.sample(intern = TRUE)` give theta itself.
#' `inla.posterior.sample()` gives each through the likelihood's declared
#' `from.theta`, which at INLA 26.8.7 is the identity for theta_1 and
#' `exp()` for the rest, so the user-scale draws are mapped back through the
#' declared `to.theta` before the formula above is applied.
#'
#' @param hyper A matrix of hyperparameter draws, or a named list of draw
#'   vectors.
#' @param K The number of categories.
#' @param scale `"internal"` for theta itself, `"user"` for draws through
#'   `from.theta`.
#' @returns A matrix, draws by K - 1 thresholds.
#' @noRd
#' @keywords internal
.fb_ordinal_pom_alpha <- function(hyper, K, scale = c("internal", "user")) {
  scale <- match.arg(scale)
  theta_names <- paste0("theta", seq_len(K - 1L), " for POM")
  theta <- if (is.matrix(hyper)) {
    hyper[, theta_names, drop = FALSE]
  } else {
    do.call(cbind, lapply(theta_names, function(nm) hyper[[nm]]))
  }
  if (identical(scale, "user")) {
    declared <- INLA::inla.models()$likelihood$pom$hyper
    for (k in seq_len(K - 1L)) {
      theta[, k] <- declared[[paste0("theta", k)]]$to.theta(theta[, k])
    }
  }
  alpha <- theta
  if (K > 2L) {
    for (k in seq(2L, K - 1L)) {
      alpha[, k] <- alpha[, k - 1L] + exp(theta[, k])
    }
  }
  unname(alpha)
}

#' The threshold table of an INLA ordinal fit
#'
#' Read from the fixed-effect table under the expansion, where
#' tau_k = -coefficient of the k-th threshold row, and from draws of the
#' hyperparameters under the cumulative likelihood.
#'
#' @param inla_fit The INLA result.
#' @param route `"pom"` or `"expansion"`.
#' @param record The ordinal record.
#' @returns A data frame as `.fb_ordinal_summarise_draws()` returns.
#' @noRd
#' @keywords internal
.fb_ordinal_thresholds_inla <- function(inla_fit, route, record) {
  K <- record$K
  nms <- .fb_ordinal_threshold_names(record$levels)
  if (identical(route, "pom")) {
    hyper <- .fb_with_seed(
      .FB_ORDINAL_THRESHOLD_SEED,
      INLA::inla.hyperpar.sample(
        .FB_ORDINAL_THRESHOLD_DRAWS,
        inla_fit,
        intern = TRUE,
        improve.marginals = TRUE
      )
    )
    return(.fb_ordinal_summarise_draws(.fb_ordinal_pom_alpha(hyper, K), nms))
  }
  sf <- inla_fit$summary.fixed
  rows <- .fb_ordinal_engine_rows(route, K)
  data.frame(
    term = nms,
    estimate = -sf[rows, "mean"],
    std.error = sf[rows, "sd"],
    conf.low = -sf[rows, "0.975quant"],
    conf.high = -sf[rows, "0.025quant"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' The threshold draws of a brms ordinal fit
#'
#' brms writes the thresholds as `b_Intercept[k]`, on the scale of the
#' uncentred predictors and in the convention this package reports.
#'
#' @param brmsfit The brms fit.
#' @param K The number of categories.
#' @returns A matrix, draws by K - 1 thresholds.
#' @noRd
#' @keywords internal
.fb_ordinal_thresholds_brms_draws <- function(brmsfit, K) {
  m <- .fb_brms_draws_matrix(brmsfit)
  cols <- paste0("b_Intercept[", seq_len(K - 1L), "]")
  unname(as.matrix(m[, cols, drop = FALSE]))
}


# Category probabilities -------------------------------------------------------

#' The inverse link of an ordinal family
#'
#' @param x A numeric vector or matrix.
#' @param link `"logit"`, `"probit"` or `"cloglog"`.
#' @returns `F(x)`, of the same shape.
#' @noRd
#' @keywords internal
.fb_ordinal_inv_link <- function(x, link) {
  switch(
    link,
    logit = stats::plogis(x),
    probit = stats::pnorm(x),
    cloglog = -expm1(-exp(x)),
    stop("flexyBayes internal error: ordinal link \"", link, "\".",
      call. = FALSE)
  )
}

#' Category probabilities given the linear predictor and the thresholds
#'
#' @param eta A numeric vector of linear predictors, thresholds excluded.
#' @param tau The K - 1 thresholds.
#' @param family `"cumulative"` or `"cratio"`.
#' @param link The link.
#' @returns A matrix, `length(eta)` by K.
#' @noRd
#' @keywords internal
.fb_ordinal_probs <- function(eta, tau, family, link) {
  n <- length(eta)
  K <- length(tau) + 1L
  if (identical(family, "cumulative")) {
    cum <- .fb_ordinal_inv_link(
      matrix(tau, n, K - 1L, byrow = TRUE) - eta,
      link
    )
    return(cbind(cum, 1) - cbind(0, cum))
  }
  # Continuation ratio: q_k = P(Y > k | Y >= k) = F(eta - tau_k)
  q <- .fb_ordinal_inv_link(eta - matrix(tau, n, K - 1L, byrow = TRUE), link)
  p <- matrix(0, n, K)
  carried <- rep(1, n)
  for (k in seq_len(K - 1L)) {
    p[, k] <- carried * (1 - q[, k])
    carried <- carried * q[, k]
  }
  p[, K] <- carried
  p
}

#' Gauss-Hermite nodes and weights
#'
#' The Golub-Welsch eigenvalue construction, for integrals of the form
#' `int f(x) exp(-x^2) dx`.
#'
#' @param n The number of nodes.
#' @returns A list: `nodes` and `weights`.
#' @noRd
#' @keywords internal
.fb_gauss_hermite <- function(n = .FB_ORDINAL_GH_NODES) {
  i <- seq_len(n - 1L)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(i, i + 1L)] <- sqrt(i / 2)
  jacobi[cbind(i + 1L, i)] <- sqrt(i / 2)
  e <- eigen(jacobi, symmetric = TRUE)
  list(nodes = e$values, weights = sqrt(pi) * e$vectors[1L, ]^2)
}

#' Category probabilities integrated over a normal random effect
#'
#' `int p_k(eta + u) N(u; 0, v) du`, by Gauss-Hermite quadrature.
#'
#' @param eta A numeric vector of linear predictors.
#' @param tau The thresholds.
#' @param v The variance of the random effect, a single number.
#' @param family,link The family and link.
#' @param rule A `.fb_gauss_hermite()` rule.
#' @returns A matrix, `length(eta)` by K.
#' @noRd
#' @keywords internal
.fb_ordinal_probs_integrated <- function(
  eta,
  tau,
  v,
  family,
  link,
  rule = .fb_gauss_hermite()
) {
  if (v <= 0) {
    return(.fb_ordinal_probs(eta, tau, family, link))
  }
  out <- 0
  for (j in seq_along(rule$nodes)) {
    shifted <- eta + sqrt(2 * v) * rule$nodes[[j]]
    out <- out + rule$weights[[j]] *
      .fb_ordinal_probs(shifted, tau, family, link)
  }
  out / sqrt(pi)
}


# Fit records ------------------------------------------------------------------

#' Finish an INLA ordinal fit
#'
#' Puts the user's IR and rows back on the fit, records the engine route,
#' reads the thresholds, and marks the class.
#'
#' @param out The `flexybayes_inla` object the emit built.
#' @param setup The `.fb_ordinal_inla_setup()` record.
#' @param user_fb The IR the user's call built.
#' @param user_data The user's rows, level-legalised.
#' @returns The fit, of class `flexybayes_ordinal` first.
#' @noRd
#' @keywords internal
.fb_ordinal_finish_inla <- function(out, setup, user_fb, user_data) {
  record <- setup$record
  user_fb$ordinal <- record
  out$fb <- user_fb
  out$data <- user_data
  out$extras$model_info$n_obs <- nrow(user_data)
  out$extras$model_info$family <- record$family
  out$extras$model_info$engine_family <- .resolve_inla_family(user_fb)
  out$extras$model_info$link <- user_fb$link
  # A refit that reads draws only keeps no hyperparameter marginals, and
  # the threshold table is not part of that readout (R/fold_refit.R).
  thresholds <- if (is.null(.fb_fold_refit_active())) {
    .fb_ordinal_thresholds_inla(out$inla, setup$route, record)
  }
  out$extras$ordinal <- c(
    record,
    list(
      link = user_fb$link,
      route = setup$route,
      engine = "inla",
      engine_rows = .fb_ordinal_engine_rows(setup$route, record$K),
      user_rows = setup$user_rows,
      thresholds = thresholds,
      threshold_prior = if (identical(setup$route, "pom")) {
        paste0(
          "a Dirichlet(",
          .FB_ORDINAL_POM_DIRICHLET,
          ") prior on the category probabilities at a linear predictor of ",
          "zero, fixed by flexyBayes"
        )
      } else {
        "the fixed-effect prior of each threshold row"
      }
    )
  )
  class(out) <- c("flexybayes_ordinal", class(out))
  out
}

#' Finish a brms ordinal fit
#'
#' @param out The `flexybayes_brms` object the emit built.
#' @param user_fb The IR the user's call built.
#' @param user_data The user's rows.
#' @returns The fit, of class `flexybayes_ordinal` first.
#' @noRd
#' @keywords internal
.fb_ordinal_finish_brms <- function(out, user_fb, user_data) {
  record <- user_fb$ordinal
  out$fb <- user_fb
  out$data <- user_data
  draws <- .fb_ordinal_thresholds_brms_draws(out$brms, record$K)
  out$extras$ordinal <- c(
    record,
    list(
      link = user_fb$link,
      route = "brms",
      engine = "brms",
      engine_rows = paste0("Intercept[", seq_len(record$K - 1L), "]"),
      user_rows = seq_len(nrow(user_data)),
      thresholds = .fb_ordinal_summarise_draws(
        draws,
        .fb_ordinal_threshold_names(record$levels)
      ),
      threshold_prior = if (!is.null(out$extras$ordinal_threshold_statement)) {
        paste0(
          "a Dirichlet(",
          .FB_ORDINAL_POM_DIRICHLET,
          ") prior on the category probabilities at a linear predictor of ",
          "zero, fixed by flexyBayes"
        )
      } else {
        "brms's prior on class Intercept"
      }
    )
  )
  class(out) <- c("flexybayes_ordinal", class(out))
  out
}

#' The brms data and IR of an ordinal model
#'
#' The response is written as its category code in an engine column, so
#' brms counts exactly the K categories the data hold.
#'
#' @param fb The `fb_terms` IR.
#' @param data The data.
#' @returns `NULL` for a family that is not ordinal, otherwise a list: `fb`
#'   and `data`.
#' @noRd
#' @keywords internal
.fb_ordinal_brms_setup <- function(fb, data) {
  if (is.null(.fb_ordinal_family(fb$family))) {
    return(NULL)
  }
  record <- fb$ordinal %||% .fb_ordinal_record(fb, data)
  fb$ordinal <- record
  .fb_ordinal_engine_column(data, .FB_ORDINAL_RESPONSE)
  data[[.FB_ORDINAL_RESPONSE]] <- .fb_ordinal_codes(fb, data, record)
  engine_fb <- fb
  engine_fb$response <- .FB_ORDINAL_RESPONSE
  list(fb = engine_fb, data = data, user_fb = fb)
}


#' The threshold prior of a cumulative brms fit
#'
#' INLA's "pom" likelihood puts a Dirichlet prior on the category
#' probabilities at a linear predictor of zero, and flexyBayes fixes its
#' concentration at 1 (`.FB_ORDINAL_POM_DIRICHLET`). At concentration 1 the
#' Dirichlet density is constant, so the prior it induces on the thresholds
#' is the product of the link's densities at each threshold, on the ordered
#' thresholds: logistic(0, 1) under the logit link, standard normal under the
#' probit link. brms gets the same prior as one Stan statement, with its own
#' threshold prior replaced by a flat row.
#'
#' brms centres the fixed-effect columns, and its `Intercept` vector holds the
#' thresholds at the column means. The thresholds at a linear predictor of
#' zero are `Intercept + dot_product(means_X, b)`, the `b_Intercept` brms
#' reports. The shift does not depend on the thresholds, so its Jacobian is 1.
#' (`bf(center = FALSE)` is not used: under brms 2.23 a cumulative model
#' without random effects then passes thresholds of 0 to its likelihood.)
#'
#' A threshold prior the caller wrote (a class `Intercept` row that did not
#' come from `prior_fixed_sd`) is kept, and no statement is added.
#'
#' @param fb The `fb_terms` IR the engine fits.
#' @param specs The brms prior specifications built so far.
#' @returns `NULL` unless the family is cumulative, otherwise a list:
#'   `specs`, the specifications with the threshold row, `stanvars`, the
#'   statement as a brms `stanvars` object or `NULL`, and `statement`, its
#'   text or `NULL`.
#' @noRd
#' @keywords internal
.fb_ordinal_brms_threshold_prior <- function(fb, specs) {
  if (!identical(.fb_ordinal_family(fb$family), "cumulative")) {
    return(NULL)
  }
  if (!identical(.FB_ORDINAL_POM_DIRICHLET, 1)) {
    stop(
      "flexyBayes internal error: the brms threshold prior is written for ",
      "a Dirichlet concentration of 1.",
      call. = FALSE
    )
  }
  is_threshold <- vapply(
    specs,
    function(s) identical(s$class, "Intercept"),
    logical(1L)
  )
  from_scalar <- vapply(
    specs,
    function(s) identical(s$input, "scalar:prior_fixed_sd"),
    logical(1L)
  )
  if (any(is_threshold & !from_scalar)) {
    return(list(specs = specs, stanvars = NULL, statement = NULL))
  }
  link <- fb$link %||% "logit"
  thresholds <- if (length(fb$fixed_terms) > 0L) {
    "Intercept + dot_product(means_X, b)"
  } else {
    "Intercept"
  }
  density <- switch(
    link,
    logit = paste0("logistic_lpdf(", thresholds, " | 0, 1)"),
    probit = paste0("std_normal_lpdf(", thresholds, ")"),
    stop(
      "flexyBayes internal error: cumulative link \"", link, "\" on brms.",
      call. = FALSE
    )
  )
  statement <- paste0("target += ", density, ";")
  specs <- c(
    specs[!(is_threshold & from_scalar)],
    list(list(
      string = "",
      class = "Intercept",
      coef = NA_character_,
      group = NA_character_,
      lb = NA_real_,
      ub = NA_real_,
      input = "threshold:Intercept"
    ))
  )
  list(
    specs = specs,
    stanvars = brms::stanvar(scode = statement, block = "model"),
    statement = statement
  )
}


# Refusal registry -------------------------------------------------------------

.populate_refusal_registry_ordinal <- function() {
  reg <- function(code, description) {
    .register_refusal(
      reason_code = code,
      description = description,
      message_template = "%s",
      registered_in_adr = "0.19.0 ordinal families",
      plan_field = NA_character_,
      since_version = "0.19.0"
    )
  }
  reg(
    "ordinal_response_unordered",
    paste0(
      "An ordinal family was given a response whose order is not stated: ",
      "a factor that is not ordered, a character or a logical vector."
    )
  )
  reg(
    "ordinal_response_not_integer",
    "An ordinal family was given a numeric response that is not whole numbers."
  )
  reg(
    "ordinal_response_too_few_categories",
    "An ordinal family was given a response with fewer than two categories."
  )
  reg(
    "ordinal_response_category_unknown",
    paste0(
      "A response value lies outside the categories of the ordinal model ",
      "it is used with."
    )
  )
  reg(
    "ordinal_fixed_without_intercept",
    paste0(
      "An ordinal model was written without an intercept, whose place its ",
      "thresholds take."
    )
  )
  reg(
    "ordinal_offset_not_supported",
    "An ordinal model was given an offset, which this version does not fit."
  )
  reg(
    "ordinal_method_not_supported",
    paste0(
      "A method that reads a mean response, a residual or a per-row ",
      "likelihood was called on an ordinal fit."
    )
  )
  reg(
    "ordinal_draws_unavailable",
    paste0(
      "Category probabilities were asked of an ordinal INLA fit that kept ",
      "no posterior configurations to draw from."
    )
  )
  reg(
    "kfold_family_not_supported",
    "fb_kfold() was called on a fit whose family it does not score."
  )
  reg(
    "kfold_folds_invalid",
    paste0(
      "fb_kfold() was given folds or groups it cannot use: the wrong ",
      "length, fewer than two folds, or a group split across folds."
    )
  )
  reg(
    "kfold_response_not_a_column",
    paste0(
      "fb_kfold() was called on a fit whose response is an expression, not ",
      "a column it can set to missing."
    )
  )
  reg(
    "kfold_category_missing_in_training",
    paste0(
      "A fold of fb_kfold() leaves a response category with no training ",
      "observation."
    )
  )
  invisible(TRUE)
}
