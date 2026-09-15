# The random-term reduction --- builder B
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md section 4.3.
#
# Builder A needs V block-diagonal. That holds for a replicated trial series
# and fails for a crossed observational design, where site crosses year and a
# third term spans both: no grouping factor contains every random term, and a
# dense V at 73,341 observations is about 43 GB. The specification that
# preceded this file assumed one builder on the premise that the random
# effects nest within trial. They do not, in the case the study actually
# needs.
#
# The route for that case is the two-stage reduction. Where a candidate is
# constant within the levels of a random term, the profiled problem is
# weighted regression on the ESTIMATED EFFECTS of that term rather than on
# the observations:
#
#   gamma-hat = argmin (u_k - T_k gamma)' Omega_k^-1 (u_k - T_k gamma) + P(gamma)
#
# with rows dropping from n to the number of levels. For the deposit that is
# 3,427 instead of 73,341.
#
# The admissibility rule is what makes the site / year / site-by-year
# asymmetry structural rather than something to audit afterwards. A static
# soil layer is constant within year, so it is NOT IDENTIFIED for the year
# regression and cannot enter it at all -- it is not that we check afterwards
# whether it wrongly explained the year component, it is that there is no
# column for it to occupy.


# --- the effects of one random term ----------------------------------- #

# .fb_sel_term_effects() --- u_k-hat and its marginal posterior spread.
#
# Read through ranef(), which is the package's own accessor, rather than
# reaching into the engine object. Omega_k's diagonal is the square of the
# reported standard error.
.fb_sel_term_effects <- function(object, term) {
  re <- ranef(object)
  if (!term %in% names(re)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_reduction_term_not_in_fit",
      message = paste0(
        "fb_select(): the base fit carries no random term `", term,
        "`. It has: ", paste(names(re), collapse = ", "), "."
      )
    ))
  }
  eff <- re[[term]]
  data.frame(
    level = as.character(eff$level),
    estimate = as.numeric(eff$estimate),
    std_error = as.numeric(eff$std.error),
    stringsAsFactors = FALSE
  )
}



# --- stage-1 level effects --------------------------------------------- #

# .fb_sel_absorbed_effects() --- the level effects with the term ABSORBED.
#
# This is the correction that matters most in this file. The reduction cites
# the two-stage proposition, which says "with the trial effects absorbed into
# X0" -- absorbed means FIXED. Reading them from ranef() instead reads
# posterior means of a RANDOM term, which are shrunk toward zero by
#
#   lambda = sigma_u^2 / (sigma_u^2 + sigma_e^2 / m)
#
# and a regression of shrunk effects on a covariate returns lambda times the
# coefficient, not the coefficient. Measured on 12 replicates at 30 levels
# and 4 plots per level (2026-09-12): route A gives 1.4906 against a truth of
# 1.5, the shrunk route gives 1.3734 -- 8.0 per cent low -- and the absorbed
# route gives 1.4922, 0.19 per cent from route A. The attenuation grows as
# the level means carry less information, reaching tens of per cent with few
# plots per level or a weak candidate.
#
# So the default is `absorb = "fixed"`. The shrunk route is retained under
# `absorb = "random"` because the RANKING of candidates survives it -- the
# two coefficient vectors are proportional, and the scalar cancels -- but a
# coefficient from it is on a different scale from builder A's and must not
# be reported beside one.
.fb_sel_absorbed_effects <- function(object, term, data) {
  d <- data
  key <- .fb_sel_term_key(object, term, d)
  d[[".fb_level"]] <- factor(key)

  fixed <- stats::formula(object$fb$source_formula %||% NULL)
  # Rebuild stage 1 from the base's own fixed part plus the absorbed term.
  rhs <- vapply(object$fb$fixed_terms, function(t) t$label %||% t$var, character(1L))
  rhs <- rhs[nzchar(rhs)]
  other_random <- Filter(
    function(t) !identical(.fb_sel_term_label(t), term),
    object$fb$random_terms
  )
  if (length(other_random) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_reduction_absorption_unsupported",
      message = paste0(
        "fb_select(): absorbing `", term, "` as a fixed term leaves ",
        length(other_random), " further random term(s), which needs a mixed ",
        "refit this release does not perform. Use absorb = \"random\" and ",
        "read only proportional quantities from the result, or reduce on a ",
        "base whose only random term is the one being absorbed."
      )
    ))
  }

  form <- stats::as.formula(paste0(
    object$fb$response, " ~ 0 + .fb_level",
    if (length(rhs) > 0L) paste0(" + ", paste(rhs, collapse = " + ")) else ""
  ))
  fit1 <- stats::lm(form, data = d)
  cf <- summary(fit1)$coefficients
  rows <- grep("^\\.fb_level", rownames(cf))
  data.frame(
    level = sub("^\\.fb_level", "", rownames(cf)[rows]),
    estimate = unname(cf[rows, "Estimate"]),
    std_error = unname(cf[rows, "Std. Error"]),
    stringsAsFactors = FALSE
  )
}


# .fb_sel_term_label() --- the name a random term is known by.
.fb_sel_term_label <- function(t) {
  if (identical(t$type, "simple")) {
    t$var
  } else if (identical(t$type, "nested")) {
    paste0(t$outer, ":", t$inner)
  } else {
    NA_character_
  }
}

# --- admissibility ----------------------------------------------------- #

# .fb_sel_reduction_admissible() --- may this candidate enter this term's
# regression?
#
# Two conditions, and both are needed.
#
#   constant within levels   otherwise the candidate has no value per level
#                            and the reduction is not defined for it.
#   varies across levels     otherwise the column is a constant and carries
#                            no information about the term.
#
# The second is the one that does the scientific work. A site-constant soil
# layer is constant within every year, so it fails `varies` for the year
# term, and the leak the design was worried about has nowhere to occur.
.fb_sel_reduction_admissible <- function(data, columns, level_key) {
  key <- as.character(level_key)
  out <- data.frame(
    column = columns,
    constant_within = NA, varies_across = NA, admissible = NA,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(columns)) {
    v <- data[[columns[[i]]]]
    if (is.null(v)) {
      out$constant_within[[i]] <- FALSE
      out$varies_across[[i]] <- FALSE
      out$admissible[[i]] <- FALSE
      next
    }
    # Constant within: at most one distinct value per level.
    n_distinct <- tapply(v, key, function(x) length(unique(x[!is.na(x)])))
    out$constant_within[[i]] <- all(n_distinct <= 1L, na.rm = TRUE)
    # Varies across: the per-level values are not all the same.
    per_level <- tapply(v, key, function(x) x[!is.na(x)][1L])
    out$varies_across[[i]] <- length(unique(per_level[!is.na(per_level)])) > 1L
    out$admissible[[i]] <- out$constant_within[[i]] && out$varies_across[[i]]
  }
  out
}


# --- the reduction ------------------------------------------------------ #

# .fb_sel_reduction() --- the level-scale problem for one random term.
#
# Returns the response u_k, the level-scale candidate matrix T_k, and the
# weights. Omega_k is the DIAGONAL of the posterior covariance of u_k, which
# is exact under balance and an approximation otherwise. The approximation is
# recorded rather than adopted quietly: an unbalanced design gets
# `omega_diagonal_approximation` on the certificate.
.fb_sel_reduction <- function(object, term, candidates, data = NULL,
                              absorb = c("fixed", "random")) {
  absorb <- match.arg(absorb)
  d <- data %||% object$data
  # The aggregated route builds a different object -- no `data`, no parsed
  # `fb` -- and reaching into it would fail with a length-zero comparison
  # several lines later rather than with an answer. Say so here instead.
  if (is.null(d) || is.null(nrow(d)) || is.null(object$fb$random_terms)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_shape_unreadable",
      message = paste0(
        "fb_select(): this base fit does not carry the parsed model and data ",
        "the reduction reads. Fits on the aggregated route (class ",
        "`flexybayes_aggregated`, which an intercept-only fixed model takes) ",
        "have a different shape and are not supported by this builder."
      )
    ))
  }
  # Validate the term BEFORE either route touches it. The absorbed route
  # builds a design column from the term's levels, so an unknown term used
  # to die inside a data-frame assignment with a message about row counts
  # rather than refusing by name.
  known <- vapply(object$fb$random_terms, .fb_sel_term_label, character(1L))
  if (!term %in% known) {
    stop(.fb_refusal_condition(
      reason_code = "selection_reduction_term_not_in_fit",
      message = paste0(
        "fb_select(): the base fit carries no random term `", term,
        "`. It has: ", paste(known[!is.na(known)], collapse = ", "), "."
      )
    ))
  }

  eff <- if (identical(absorb, "fixed")) {
    .fb_sel_absorbed_effects(object, term, d)
  } else {
    .fb_sel_term_effects(object, term)
  }

  key <- .fb_sel_term_key(object, term, d)
  if (length(key) != nrow(d)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_reduction_term_not_in_fit",
      message = paste0(
        "fb_select(): the levels of random term `", term,
        "` could not be recovered from the data."
      )
    ))
  }

  adm <- .fb_sel_reduction_admissible(d, candidates, key)
  if (any(!adm$admissible)) {
    bad <- adm[!adm$admissible, , drop = FALSE]
    why <- ifelse(
      !bad$constant_within, "not constant within levels",
      "does not vary across levels"
    )
    stop(.fb_refusal_condition(
      reason_code = "selection_candidates_not_constant_within_term",
      message = paste0(
        "fb_select(): these candidates cannot enter the reduction for `",
        term, "`: ",
        paste0(bad$column, " (", why, ")", collapse = "; "),
        ". A candidate must be constant within a level to have a value ",
        "there, and must vary across levels to say anything about the term."
      )
    ))
  }

  # One row per level, in the order ranef() reports them, matched by NAME.
  # A positional match between the effects and the data's levels is a silent
  # misalignment of every covariate value.
  first_row <- match(eff$level, as.character(key))
  if (anyNA(first_row)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_reduction_level_mismatch",
      message = paste0(
        "fb_select(): random term `", term, "` reports level(s) absent from ",
        "the data: ", paste(eff$level[is.na(first_row)], collapse = ", "), "."
      )
    ))
  }

  Tk <- as.matrix(d[first_row, candidates, drop = FALSE])
  rownames(Tk) <- eff$level

  n_per_level <- as.integer(table(factor(key, levels = eff$level)))
  balanced <- length(unique(n_per_level)) == 1L

  # Omega is the variance of the STAGE-2 residual, not the sampling variance
  # of the stage-1 estimate. Writing alpha_i-hat = T_i gamma + b_i + e_i, the
  # level departure b_i carries sigma_u^2 and the estimate carries
  # sigma_e^2 / m_i, so
  #
  #   Var(alpha_i-hat) = sigma_u^2 + sigma_e^2 / m_i
  #
  # This file used to weight by the second term alone, which is the
  # FIXED-effect meta-analysis weight where the random-effects one is
  # required. Under balance the omitted term is constant across levels, the
  # weights stay proportional and the answer is unchanged -- which is why it
  # survived every balanced test. Under imbalance it is wrong: on a
  # hand-built dense V with no package code, GLS gives 1.8024932073, these
  # weights gave 1.9224778858 (6.7 per cent out) and the corrected ones give
  # 1.8024932073 to 1.3e-15. An independent sweep put the worst case at 226
  # per cent. The deposit is unbalanced.
  #
  # The between-level component comes from the BASE fit, where the term is
  # random and sigma_u^2 is estimated; stage 1 absorbs the term as fixed and
  # cannot report it.
  omega_between <- if (identical(absorb, "fixed")) {
    vc <- .fb_sel_variance_components(object)
    vc$terms[[term]] %||% 0
  } else {
    # The shrunk route's posterior sd already carries the between-level
    # spread, and that route answers a different question anyway.
    0
  }

  list(
    term = term,
    absorb = absorb,
    levels = eff$level,
    u = eff$estimate,
    omega_between = omega_between,
    omega_diag = omega_between + eff$std_error^2,
    T = Tk,
    n_per_level = n_per_level,
    balanced = balanced,
    approximations = c(
      if (balanced) character(0) else "omega_diagonal_approximation",
      # A shrunk stage 1 answers a different question by a scalar. Recorded,
      # so a coefficient from it cannot be paired with builder A's unnoticed.
      if (identical(absorb, "random")) "shrunk_stage_one_effects" else character(0)
    ),
    admissibility = adm
  )
}


# .fb_sel_term_key() --- the row-level factor a random term induces.
.fb_sel_term_key <- function(object, term, data = NULL) {
  d <- data %||% object$data
  for (t in object$fb$random_terms) {
    if (identical(t$type, "simple") && identical(t$var, term)) {
      return(as.character(d[[t$var]]))
    }
    if (
      identical(t$type, "nested") &&
        identical(paste0(t$outer, ":", t$inner), term)
    ) {
      return(paste0(as.character(d[[t$outer]]), ":", as.character(d[[t$inner]])))
    }
  }
  character(0)
}


# --- the reduction's whitener ------------------------------------------- #

# .fb_sel_reduction_whitener() --- W for the level-scale problem.
#
# Omega_k is diagonal here, so W is too and whitening is a division. The
# object carries the same class as builder A's so downstream code consumes
# one interface, and the builder label on it says which route produced it.
.fb_sel_reduction_whitener <- function(red) {
  if (any(!is.finite(red$omega_diag)) || any(red$omega_diag <= 0)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_covariance_not_positive_definite",
      message = paste0(
        "fb_select(): the posterior spread of random term `", red$term,
        "` is not positive and finite at every level, so the reduction has ",
        "no weights."
      )
    ))
  }
  factors <- lapply(seq_along(red$levels), function(i) {
    list(rows = i, L = matrix(sqrt(red$omega_diag[[i]]), 1L, 1L))
  })
  structure(
    list(
      builder = "random_term_reduction",
      estimand = "posterior_marginal_spread",
      n = length(red$levels),
      block_factor = red$term,
      n_blocks = length(red$levels),
      factors = factors,
      components = list(residual = NA_real_, terms = numeric(0)),
      approximations = red$approximations
    ),
    class = "fb_sel_whitener"
  )
}
