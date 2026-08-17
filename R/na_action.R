# ------------------------------------------------------------------ #
# Missing responses: the design-preserving augmentation layer          #
# ------------------------------------------------------------------ #
#
# A missing response removes an observation, not a design node. The
# distinction matters whenever the model carries a covariance structure
# indexed by the design -- a separable AR1(row) x AR1(col) residual over
# a field trial, say. Deleting the row of a lost plot does not merely
# discard one number: it changes the index set the covariance is built
# over, so the model that gets fitted is a different model from the one
# that was written down.
#
# The classical treatment is Yates's missing-plot technique, put on a
# likelihood footing by Houtman and Speed (1984) and Verbyla and Cullis
# (1992): carry the missing cell as a free quantity, and the analysis of
# the completed grid reproduces the observed-data likelihood exactly.
# Read in a Bayesian frame that construction is data augmentation
# (Tanner and Wong 1987): the missing response is a latent quantity,
# integrated out. Under ignorability (Rubin 1976) the posterior for the
# model parameters is then the same whether the cell is imputed or
# omitted from the likelihood -- what augmentation preserves is the
# REPRESENTATION, not information.
#
# That is why this layer adds no inference machinery. INLA already
# treats an NA response as a latent prediction target, and marginalises
# it; brms reaches the same object through `mi()`. The work here is
# keeping the design intact so the engine is handed the right index set:
#
#   na_action = "augment" -- retain rows whose response is missing, and
#     complete the design grid where the absent cells are determinable.
#     Matches ASReml's na.method(y = "include").
#   na_action = "omit"    -- drop them (complete-case). A structured
#     covariance over a broken grid then refuses downstream.
#   na_action = "fail"    -- refuse if any response is missing.
#
# Ignorability is an assumption, not a property of the device. Under
# missingness that depends on the unobserved response itself, augmenting
# and deleting are both biased, and equally so; the layer does not
# detect this and cannot.
#
# Missing COVARIATES are refused throughout. ASReml drops or zero-fills
# them; a zero-filled covariate is a fabricated observation, so this is
# a deliberate departure and a stricter one.

# .fb_na_action_choices() --- the closed vocabulary.
.fb_na_action_choices <- function() c("augment", "omit", "fail")

# .fb_response_missing() --- index of rows whose response is NA.
.fb_response_missing <- function(fb, data) {
  resp <- fb$response
  if (is.null(resp) || !nzchar(resp) || !resp %in% names(data)) {
    return(integer(0))
  }
  which(is.na(data[[resp]]))
}

# .fb_design_index_vars() --- the variables a structured covariance is
# indexed by. These are the columns whose completeness the covariance
# representation depends on: an ar1_spatial term is emitted as a field
# over the row x col grid, so both must be present for every node.
# Returns character(0) when no term carries a design index, in which
# case there is no grid to complete and augmentation reduces to
# retaining the missing-response rows.
.fb_design_index_vars <- function(fb) {
  out <- character(0)
  terms_all <- c(fb$random_terms %||% list(), fb$residual_terms %||% list())
  for (t in terms_all) {
    ty <- t$type %||% ""
    if (identical(ty, "ar1_spatial")) {
      out <- c(out, t$row_var, t$col_var)
    } else if (identical(ty, "ar1")) {
      out <- c(out, t$var)
    }
  }
  unique(out[!is.na(out) & nzchar(out)])
}

# .fb_model_vars() --- every column the model refers to, so the grid
# completion can tell which values it would have to invent.
.fb_model_vars <- function(fb) {
  out <- character(0)
  collect <- function(tl) {
    for (t in tl %||% list()) {
      for (f in c("var", "row_var", "col_var", "group", "inner", "outer")) {
        v <- t[[f]]
        if (is.character(v)) out <<- c(out, v)
      }
      if (is.character(t$vars)) out <<- c(out, t$vars)
      if (!is.null(t$expr)) out <<- c(out, all.vars(t$expr))
      if (is.character(t$label) && !grepl("[^A-Za-z0-9._]", t$label)) {
        out <<- c(out, t$label)
      }
    }
  }
  collect(fb$fixed_terms)
  collect(fb$random_terms)
  collect(fb$residual_terms)
  unique(out[!is.na(out) & nzchar(out)])
}

# .fb_complete_design_grid() --- reinsert absent design cells.
#
# Returns the data with one NA-response row added per absent combination
# of the design index variables, or refuses when a model variable's
# value at an absent cell cannot be determined. Determinable means: it
# is one of the index variables, or it takes a single value across the
# whole dataset. Anything else -- a genotype label, a block, a measured
# covariate -- is a real quantity that was never observed at that cell,
# and inventing one would fabricate an observation.
#
# The refusal is not a limitation of the device so much as a reflection
# of how the data arrive. A field book lists every plot that was sown,
# so a lost yield is a row present with an NA, which needs no invention.
# A cell absent from the data frame entirely is a cell whose design
# assignment nobody recorded.
.fb_complete_design_grid <- function(fb, data, index_vars) {
  if (length(index_vars) == 0L) {
    return(list(data = data, added = 0L))
  }
  if (!all(index_vars %in% names(data))) {
    return(list(data = data, added = 0L))
  }

  levs <- lapply(index_vars, function(v) {
    x <- data[[v]]
    if (is.factor(x)) levels(x) else sort(unique(x))
  })
  names(levs) <- index_vars
  full <- expand.grid(levs, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

  key <- function(df) {
    do.call(paste, c(lapply(index_vars, function(v) as.character(df[[v]])),
                     list(sep = "\r")))
  }
  have <- key(data)
  want <- key(full)
  absent <- which(!want %in% have)
  if (length(absent) == 0L) {
    return(list(data = data, added = 0L))
  }

  # Which columns would have to be invented for the absent rows?
  needed <- setdiff(intersect(.fb_model_vars(fb), names(data)), index_vars)
  undetermined <- needed[vapply(
    needed,
    function(v) length(unique(data[[v]][!is.na(data[[v]])])) > 1L,
    logical(1)
  )]
  if (length(undetermined) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "augment_cell_not_determinable",
      message = paste0(
        "na_action = \"augment\" cannot complete the design grid: ",
        length(absent), " cell(s) of the ",
        paste(index_vars, collapse = " x "),
        " grid are absent from the data, and the model also depends on ",
        paste(undetermined, collapse = ", "),
        ", whose value at those cells is not recorded anywhere. ",
        "Inventing it would fabricate an observation. Supply the absent ",
        "cells as rows with the design columns filled in and the ",
        "response set to NA -- which is how a field book records a lost ",
        "plot -- or use na_action = \"omit\" to analyse the observed ",
        "cells only."
      ),
      family_class = "flexybayes_augment_cell_not_determinable"
    ))
  }

  add <- data[rep(1L, length(absent)), , drop = FALSE]
  rownames(add) <- NULL
  for (v in index_vars) {
    val <- full[[v]][absent]
    add[[v]] <- if (is.factor(data[[v]])) {
      factor(as.character(val), levels = levels(data[[v]]))
    } else {
      val
    }
  }
  add[[fb$response]] <- NA
  # Constant columns keep their single value; everything else the model
  # does not refer to is carried from row 1 and never read.
  out <- rbind(data, add)
  rownames(out) <- NULL
  list(data = out, added = length(absent))
}

# .fb_apply_na_action() --- the entry point, called from flexybayes()
# once the IR exists and before dispatch.
#
# Returns list(data = <possibly modified>, meta = <record of what was
# done>). The metadata travels into the fit so a reader can see whether
# a posterior was computed on the design as laid out or on the plots
# that happened to survive -- a distinction the fit object could not
# previously express.
.fb_apply_na_action <- function(fb, data, na_action) {
  na_action <- match.arg(na_action, .fb_na_action_choices())
  miss <- .fb_response_missing(fb, data)

  # A missing covariate is refused under every na_action. There is no
  # ignorability argument that makes a fabricated predictor safe, and
  # the device does not extend to one.
  cov_vars <- setdiff(intersect(.fb_model_vars(fb), names(data)), fb$response)
  bad_cov <- cov_vars[vapply(
    cov_vars, function(v) anyNA(data[[v]]), logical(1)
  )]
  if (length(bad_cov) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "missing_covariate_not_supported",
      message = paste0(
        "Missing values in the covariate(s) ",
        paste(bad_cov, collapse = ", "),
        ". The missing-response device does not extend to predictors: ",
        "a filled-in covariate is a fabricated observation, not a ",
        "marginalised latent quantity. Supply the values, or drop the ",
        "affected rows deliberately before fitting."
      ),
      family_class = "flexybayes_missing_covariate_not_supported"
    ))
  }

  if (identical(na_action, "fail") && length(miss) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "missing_response_refused",
      message = paste0(
        length(miss), " observation(s) have a missing response and ",
        "na_action = \"fail\". Use na_action = \"augment\" to carry ",
        "them as latent quantities and keep the design intact, or ",
        "na_action = \"omit\" to analyse the observed cells only."
      ),
      family_class = "flexybayes_missing_response_refused"
    ))
  }

  if (identical(na_action, "omit")) {
    if (length(miss) > 0L) {
      data <- data[-miss, , drop = FALSE]
      rownames(data) <- NULL
    }
    return(list(
      data = data,
      meta = list(
        na_action = "omit", n_missing_response = length(miss),
        n_cells_completed = 0L, design_index_vars = character(0)
      )
    ))
  }

  # augment
  index_vars <- .fb_design_index_vars(fb)
  completed <- .fb_complete_design_grid(fb, data, index_vars)
  n_unobserved <- length(miss) + completed$added
  .fb_warn_high_missingness(n_unobserved, nrow(completed$data))
  list(
    data = completed$data,
    meta = list(
      na_action = "augment",
      n_missing_response = length(miss),
      n_cells_completed = completed$added,
      design_index_vars = index_vars,
      missing_fraction = if (nrow(completed$data) > 0L) {
        n_unobserved / nrow(completed$data)
      } else {
        NA_real_
      }
    )
  )
}


# .fb_warn_high_missingness() --- say what a high missing fraction does to
# the estimand, once per session.
#
# The augmentation identity is algebraic and does not weaken as the missing
# fraction rises. What weakens is the restricted likelihood it targets:
# past roughly a third unobserved, variance components reach the boundary,
# the surface flattens, and two correct implementations stop at different
# points. That is a property of the design, not of the device, and the fit
# proceeds -- but a user comparing this posterior against ASReml or lme4
# needs to know it before reading a disagreement as a defect in either.
#
# Once per session rather than once per fit, in the manner of the routing
# notes: a simulation study fitting the same design a thousand times
# should say this once.
.FB_HIGH_MISSING_FRACTION <- 0.30

.fb_warn_high_missingness <- function(n_unobserved, n_rows) {
  if (n_rows <= 0L || n_unobserved <= 0L) {
    return(invisible(NULL))
  }
  frac <- n_unobserved / n_rows
  if (frac <= .FB_HIGH_MISSING_FRACTION) {
    return(invisible(NULL))
  }
  if (isTRUE(getOption("flexyBayes.silence_high_missingness_warning", FALSE))) {
    return(invisible(NULL))
  }
  if (.emit_state_get("high_missingness_warning")) {
    return(invisible(NULL))
  }
  .emit_state_set("high_missingness_warning", TRUE)
  warning(
    "flexyBayes: ", format(round(100 * frac, 1)), "% of the fitted rows ",
    "carry no observed response (", n_unobserved, " of ", n_rows,
    "). Above roughly 30% the variance components can sit near a ",
    "boundary or on a flat ridge, where two correct implementations ",
    "legitimately stop at different points -- so a disagreement with ",
    "ASReml, lme4 or another engine at this missing fraction is a ",
    "property of the estimand rather than evidence against either fit. ",
    "The augmentation identity itself is unaffected. Silence via ",
    "options(flexyBayes.silence_high_missingness_warning = TRUE).",
    call. = FALSE
  )
  invisible(NULL)
}
