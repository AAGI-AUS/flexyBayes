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
# Missing COVARIATES have their own policy, written the way ASReml
# writes it (`na.method(x = )`). The default refuses: ASReml's own
# default is `x = "fail"`, and its alternative -- zero-filling the
# missing predictor (Manual section 3.11) -- fabricates an observation,
# so that setting is refused by name rather than reproduced. `x = "omit"`
# drops the affected rows and says how many and which columns.

# .fb_na_action_choices() --- the closed vocabulary.
.fb_na_action_choices <- function() c("augment", "omit", "fail")

# .fb_na_covariate_choices() --- the closed vocabulary for the covariate
# half of the policy, in ASReml's spelling (there is no flexyBayes-native
# spelling: the argument exists to accept the one users already write).
.fb_na_covariate_choices <- function() c("fail", "omit", "include")


# ------------------------------------------------------------------ #
# Accepting the policy an ASReml user already writes                   #
# ------------------------------------------------------------------ #
#
# ASReml states the same two decisions through one object,
# `na.method(y = , x = )`: what to do with a missing RESPONSE and what to
# do with a missing COVARIATE. flexyBayes accepts that object, the bare
# list a reader without an asreml licence can write by hand, and its own
# native strings, and reduces all three to one internal policy.
#
# Two recorded properties of `asreml::na.method()` drive the reduction
# (asreml 4.2.0.392, recorded 2026-08-17; fixture in
# tests/testthat/helper-asreml-shapes.R):
#
#   1. The value is a PLAIN list -- class() is "list", the only attribute
#      is `names`. There is no class to dispatch on, so detection is by
#      shape: a list whose names are drawn from `x` / `y` and nothing
#      else.
#
#   2. An unsupplied argument is NOT reduced to its default scalar. The
#      whole default vector comes back, so `na.method(y = "include")`
#      carries `x = c("fail", "include", "omit")`. The effective policy
#      is the first element, as `match.arg()` would take it. Reading the
#      slot as a scalar silently mis-reads a partially specified call.

# .fb_na_method_y_vocabulary() --- ASReml's response-side words and the
# flexyBayes word each one means. The native spellings map to themselves
# so a normalised policy can be normalised again without changing.
.fb_na_method_y_vocabulary <- function() {
  c(
    include = "augment",
    omit = "omit",
    fail = "fail",
    augment = "augment"
  )
}

# .fb_match_na_slot() --- reduce one slot of an incoming policy to a
# single word of a closed vocabulary.
#
# `value` is whatever arrived in the slot: a scalar, or the whole default
# vector of an unsupplied `na.method()` argument (property 2 above). The
# first element is the effective policy either way. An unrecognised word
# is named rather than silently defaulted.
.fb_match_na_slot <- function(value, choices, slot) {
  if (is.null(value) || length(value) == 0L) {
    return(choices[[1L]])
  }
  if (!is.character(value)) {
    stop(
      "`na_action`: the `", slot, "` slot must be a character value, not ",
      "an object of class <", paste(class(value), collapse = ", "), ">.",
      call. = FALSE
    )
  }
  word <- value[[1L]]
  if (!word %in% choices) {
    stop(
      "`na_action`: \"", word, "\" is not a recognised `", slot,
      "` policy. Use one of ", paste0("\"", choices, "\"", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  word
}

# .fb_is_na_method_shape() --- the shape test, and only the shape test.
#
# True for `asreml::na.method(...)` output and for the hand-written list
# that stands in for it. Never `inherits()` and never `class()`: the
# recorded value carries no class of its own, so there is nothing to
# match on.
.fb_is_na_method_shape <- function(v) {
  if (!is.list(v) || length(v) == 0L) {
    return(FALSE)
  }
  nms <- names(v)
  if (is.null(nms) || any(!nzchar(nms))) {
    return(FALSE)
  }
  all(nms %in% c("x", "y")) && any(c("x", "y") %in% nms)
}

# .fb_normalise_na_action() --- one internal policy from any accepted
# spelling.
#
# Returns `list(y = <"augment"|"omit"|"fail">, x = <"fail"|"omit"|
# "include">)`. Idempotent: a policy that has already been normalised
# passes through unchanged, which is what lets flexybayes() normalise
# once for the record and hand the same value on to be applied.
#
# A slot the caller left out takes ASReml's own default for it --
# `y = "include"` (so `"augment"`, which is also the flexyBayes default)
# and `x = "fail"`.
.fb_normalise_na_action <- function(na_action) {
  if (.fb_is_na_method_shape(na_action)) {
    y_word <- .fb_match_na_slot(
      na_action$y,
      names(.fb_na_method_y_vocabulary()),
      "y"
    )
    return(list(
      y = unname(.fb_na_method_y_vocabulary()[[y_word]]),
      x = .fb_match_na_slot(
        na_action$x,
        .fb_na_covariate_choices(),
        "x"
      )
    ))
  }

  if (is.list(na_action)) {
    stop(
      "`na_action` was given as a list, but its names are ",
      if (is.null(names(na_action))) {
        "absent"
      } else {
        paste0("`", paste(names(na_action), collapse = "`, `"), "`")
      },
      ". A policy list is the shape asreml::na.method() returns: named ",
      "`y` (the response policy) and optionally `x` (the covariate ",
      "policy), for example list(y = \"include\", x = \"fail\").",
      call. = FALSE
    )
  }

  list(
    y = match.arg(na_action, .fb_na_action_choices()),
    x = .fb_na_covariate_choices()[[1L]]
  )
}


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
#
# The completion itself is `.fb_grid_complete()` in R/fb_complete_grid.R,
# which the exported `fb_complete_grid()` also runs. One implementation,
# two entry points, so the helper a user calls by hand and the automatic
# path a default fit takes cannot complete the same trial two different
# ways. What differs here is the candidate set -- only the columns the
# model refers to matter to a fit, where the standalone helper has no IR
# and considers every column -- and the remedy the refusal names.
#
# Constant columns keep their single value; everything else the model
# does not refer to is carried from row 1 and never read.
.fb_complete_design_grid <- function(fb, data, index_vars) {
  out <- .fb_grid_complete(
    data = data,
    index_vars = index_vars,
    response = fb$response,
    candidate_vars = .fb_model_vars(fb),
    unused_level = NULL,
    remedy = paste0(
      "Supply the field-book rows instead -- one row per sown plot, with ",
      "the design columns filled in and the response set to NA, which is ",
      "how a field book already records a lost plot -- or pad the array ",
      "yourself with ",
      "fb_complete_grid(data, ~ row * col, response, ",
      "unused_level = \"UNSOWN\"), which names the level the empty cells ",
      "carry; or use na_action = \"omit\" to analyse the observed cells ",
      "only."
    )
  )
  list(data = out$data, added = out$added)
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
  policy <- .fb_normalise_na_action(na_action)

  # ---- Covariates first -------------------------------------------------
  # The response policy is about which design cells the engine is handed.
  # The covariate policy is about which rows exist at all, so it settles
  # before the design is counted or completed.
  cov_step <- .fb_apply_covariate_policy(fb, data, policy$x)
  data <- cov_step$data

  miss <- .fb_response_missing(fb, data)

  if (identical(policy$y, "fail")) {
    if (length(miss) > 0L) {
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
    # Nothing missing, so nothing to do. The grid is deliberately NOT
    # completed here: completing it would introduce the missing responses
    # this setting exists to refuse. Before 0.9.1 control fell through
    # into the augment branch, which completed the grid and recorded the
    # run as `augment` -- so the record did not echo the argument.
    return(list(
      data = data,
      meta = .fb_na_action_meta(
        policy = policy,
        n_missing_response = 0L,
        n_cells_completed = 0L,
        design_index_vars = character(0),
        n_design = nrow(data),
        n_unobserved = 0L,
        n_covariate_rows_dropped = cov_step$n_dropped
      )
    ))
  }

  if (identical(policy$y, "omit")) {
    if (length(miss) > 0L) {
      data <- data[-miss, , drop = FALSE]
      rownames(data) <- NULL
    }
    return(list(
      data = data,
      meta = .fb_na_action_meta(
        policy = policy,
        n_missing_response = length(miss),
        n_cells_completed = 0L,
        design_index_vars = character(0),
        n_design = nrow(data),
        n_unobserved = 0L,
        n_covariate_rows_dropped = cov_step$n_dropped
      )
    ))
  }

  # ---- augment ----------------------------------------------------------
  index_vars <- .fb_design_index_vars(fb)
  completed <- .fb_complete_design_grid(fb, data, index_vars)
  n_unobserved <- length(miss) + completed$added
  .fb_warn_high_missingness(n_unobserved, nrow(completed$data))
  list(
    data = completed$data,
    meta = .fb_na_action_meta(
      policy = policy,
      n_missing_response = length(miss),
      n_cells_completed = completed$added,
      design_index_vars = index_vars,
      n_design = nrow(completed$data),
      n_unobserved = n_unobserved,
      n_covariate_rows_dropped = cov_step$n_dropped
    )
  )
}

# .fb_na_action_meta() --- the record that travels into the fit.
#
# One builder for all three response policies so the fields cannot drift
# apart between branches, which is how `missing_fraction` came to exist
# on the augment path alone.
#
# `n_design` is the number of rows handed to the engine and `n_observed`
# the number of those rows carrying an observed response, so
# `missing_fraction` reads the same way on every path: the fraction of
# the FITTED rows that carry no observed response. On `omit` and `fail`
# that is zero by construction -- how many rows were dropped or would
# have been is `n_missing_response`, which is a different number and has
# its own field.
.fb_na_action_meta <- function(
  policy,
  n_missing_response,
  n_cells_completed,
  design_index_vars,
  n_design,
  n_unobserved,
  n_covariate_rows_dropped
) {
  list(
    na_action = policy$y,
    na_covariate = policy$x,
    n_missing_response = as.integer(n_missing_response),
    n_cells_completed = as.integer(n_cells_completed),
    design_index_vars = design_index_vars,
    n_design = as.integer(n_design),
    n_observed = as.integer(n_design - n_unobserved),
    n_covariate_rows_dropped = as.integer(n_covariate_rows_dropped),
    missing_fraction = if (n_design > 0L) {
      n_unobserved / n_design
    } else {
      NA_real_
    }
  )
}

# .fb_apply_covariate_policy() --- the covariate half of na.method().
#
# Returns list(data = <possibly reduced>, n_dropped = <int>). Called
# before anything looks at the response, because two of the three
# policies decide whether a row is in the analysis at all.
#
# The default (`x = "fail"`) is ASReml's own default, and it is the
# setting that refuses. `x = "include"` is the one ASReml behaviour this
# package will not reproduce.
.fb_apply_covariate_policy <- function(fb, data, x_policy) {
  cov_vars <- setdiff(intersect(.fb_model_vars(fb), names(data)), fb$response)
  bad_cov <- cov_vars[vapply(
    cov_vars, function(v) anyNA(data[[v]]), logical(1)
  )]
  if (length(bad_cov) == 0L) {
    return(list(data = data, n_dropped = 0L))
  }

  named <- paste(bad_cov, collapse = ", ")

  if (identical(x_policy, "fail")) {
    stop(.fb_refusal_condition(
      reason_code = "missing_covariate_not_supported",
      message = paste0(
        "The predictor(s) ", named, " have missing values, and the ",
        "covariate policy is `fail` -- which is also what ASReml does by ",
        "default (na.method(x = \"fail\")). A missing predictor is not a ",
        "missing observation: the response device carries an unobserved ",
        "yield as a latent quantity the engine integrates out, and there ",
        "is no equivalent for a predictor whose value was never ",
        "recorded. Supply the values, drop the affected rows yourself, ",
        "or say so in the call with ",
        "na_action = list(y = \"include\", x = \"omit\"), which drops ",
        "them and reports how many."
      ),
      family_class = "flexybayes_missing_covariate_not_supported"
    ))
  }

  if (identical(x_policy, "include")) {
    stop(.fb_refusal_condition(
      reason_code = "covariate_zero_fill_not_supported",
      message = paste0(
        "The predictor(s) ", named, " have missing values and the ",
        "covariate policy is `include`. ASReml treats a missing covariate ",
        "as zero under that setting (ASReml-R Reference Manual 4.2, ",
        "section 3.11). A zero is a value the plot did not have: the fit ",
        "would report a coefficient estimated partly from observations ",
        "nobody made, and nothing downstream would say so. Supply the ",
        "values, or pass ",
        "na_action = list(y = \"include\", x = \"omit\") to drop the ",
        "affected rows and see the count."
      ),
      family_class = "flexybayes_covariate_zero_fill_not_supported"
    ))
  }

  # x = "omit": drop the affected rows, and say what went.
  drop_rows <- Reduce(
    `|`,
    lapply(bad_cov, function(v) is.na(data[[v]])),
    init = rep(FALSE, nrow(data))
  )
  n_dropped <- sum(drop_rows)
  warning(
    "flexyBayes: dropped ", n_dropped, " row(s) of ", nrow(data),
    " whose predictor value was missing in ", named,
    " (na_action covariate policy x = \"omit\"). Deletion is a ",
    "complete-case analysis of the remaining rows; a covariance indexed ",
    "by the design is built over the reduced index set.",
    call. = FALSE
  )
  data <- data[!drop_rows, , drop = FALSE]
  rownames(data) <- NULL
  list(data = data, n_dropped = as.integer(n_dropped))
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
