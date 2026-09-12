# The two decisions that are not the projection.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md section 6.2.
#
# `full` and `posterior_threshold` share the projection's return shape as far
# as the shape is meaningful, and are explicit about the fields they cannot
# fill. A decision object whose `chosen_terms` is populated by one route and
# absent from another is a shape a consumer has to branch on twice; a field
# that is present and empty, beside a `decision` label saying why, is one.
#
# Neither route runs a search, so neither carries an elpd-by-size table.
# `full` chooses nothing at all -- that is what makes it `full` -- and
# `posterior_threshold` chooses by reading the reference posterior directly.


# --- the reference posterior, unselected -------------------------------- #

# .fb_sel_decision_full() --- the reference model with no selection.
#
# Every candidate is retained: `full` is the answer for a user who wanted
# shrinkage and did not want a decision, and reporting it as "nothing was
# chosen" would misdescribe it.
.fb_sel_decision_full <- function(fit, candidate_terms, protected_terms) {
  list(
    decision = "full",
    chosen_terms = as.character(candidate_terms %||% character(0)),
    # The vocabulary the chosen set is written in, so a consumer matching
    # against a design does not have to guess. `full` chooses formula terms.
    chosen_vocabulary = "term",
    forced_terms = as.character(protected_terms %||% character(0)),
    chosen_size = length(candidate_terms %||% character(0)),
    label = paste0(
      "the reference posterior under the joint shrinkage prior, with no ",
      "selection: every candidate is retained and shrunk"
    ),
    elpd_by_size = NULL,
    ranking = NULL,
    rank_stability = NULL,
    projected_draws = NULL,
    settings = list(method = "full")
  )
}


# --- the posterior-probability rule -------------------------------------- #

# .fb_sel_posterior_threshold() --- retain candidate j when
# P(|beta_j| > effect_size) > threshold, on the SCIENTIFIC scale.
#
# Three things this is not, all of them stated on the object because each has
# been read into a rule of this shape before:
#
#   * it is not a posterior inclusion probability. A continuous shrinkage
#     prior puts no mass at zero, so there is no event "beta_j = 0" whose
#     probability could be reported. The quantity here is the posterior
#     probability of a magnitude exceeding a threshold the analyst declared;
#   * it is not calibrated. Nothing here controls a false-discovery rate or
#     any other frequentist quantity, and the threshold is not tuned;
#   * it is not scale-free. `effect_size` is compared against the coefficient
#     of a one-unit move in the candidate column, so the draws are
#     back-transformed off the fold-standardised scale before the comparison
#     and the result would be a different set at a different scaling.
#
# The draws are read off the fitted brms object and the block scaling is
# undone by `.fb_pen_back_transform()`, which is the same function the
# coefficient table uses -- one back-transformation, not two that could
# disagree.
.fb_sel_posterior_threshold <- function(fit,
                                        blocks,
                                        decision,
                                        candidate_terms = NULL,
                                        protected_terms = NULL) {
  spec <- .fb_sel_decision_spec(decision)
  effect_size <- spec$effect_size
  if (!is.numeric(effect_size) || length(effect_size) != 1L ||
      !is.finite(effect_size)) {
    stop(
      ".fb_sel_posterior_threshold(): the decision carries no `effect_size`; ",
      "fb_decision(method = \"posterior_threshold\") refuses without one.",
      call. = FALSE
    )
  }

  brmsfit <- if (inherits(fit, "brmsfit")) fit else fit$brms
  if (is.null(brmsfit)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_shape_unreadable",
      message = paste0(
        "fb_select(): the posterior-threshold decision reads the reference ",
        "model's posterior draws, and this object carries none."
      )
    ))
  }

  draws <- as.matrix(brmsfit)
  coef_names <- unique(unlist(
    lapply(blocks %||% list(), function(b) b$coef_names),
    use.names = FALSE
  ))
  if (length(coef_names) == 0L) {
    stop(
      ".fb_sel_posterior_threshold(): no candidate coefficient is named by ",
      "the compiled prior blocks, so there is nothing to threshold.",
      call. = FALSE
    )
  }
  # Only the columns the blocks name, plus the intercept the
  # back-transformation moves with them.
  keep <- colnames(draws) %in% c(paste0("b_", coef_names), "b_Intercept")
  scaled <- .fb_pen_back_transform(draws[, keep, drop = FALSE], blocks)

  present <- paste0("b_", coef_names) %in% colnames(scaled)
  probability <- vapply(
    coef_names[present],
    function(nm) mean(abs(as.numeric(scaled[, paste0("b_", nm)])) >
      effect_size),
    numeric(1L)
  )
  retained <- probability > spec$threshold

  list(
    decision = "posterior_threshold",
    chosen_terms = names(probability)[retained],
    # This is the one decision that chooses coefficient by coefficient: the
    # probability is a property of a single draw column, not of a term.
    chosen_vocabulary = "coefficient",
    forced_terms = as.character(protected_terms %||% character(0)),
    chosen_size = sum(retained),
    coefficients = data.frame(
      coefficient = names(probability),
      probability = unname(probability),
      retained = unname(retained),
      stringsAsFactors = FALSE
    ),
    label = paste0(
      "candidates whose posterior probability of exceeding a magnitude of ",
      format(effect_size),
      " on the scientific scale is above ",
      format(spec$threshold),
      "; a decision rule with no inclusion-probability interpretation"
    ),
    unreadable = coef_names[!present],
    elpd_by_size = NULL,
    ranking = NULL,
    rank_stability = NULL,
    projected_draws = NULL,
    settings = list(
      method = "posterior_threshold",
      threshold = spec$threshold,
      effect_size = effect_size,
      scale = "scientific"
    )
  )
}
