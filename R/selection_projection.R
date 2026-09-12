# The decision layer --- projection predictive selection, and its gate
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 6.1 and 6.3.
#
# Two things live here and they are deliberately separate.
#
# The STRUCTURE GATE answers, before any search runs, whether the
# reference model is one this provider can project at all. Its verdicts
# come from runs: `get_refmodel()` accepted a grouped random intercept and
# a validated grouped k-fold search ran on it (probes C1-C3); it refused a
# known-covariance carrier (probe C5) and a non-linear partition (probe
# C6). Every row that has NOT been executed refuses exactly as a row that
# was executed and came back red. Unrun is not a permission, and turning a
# row green means writing a test that runs it.
#
# The PROJECTION itself runs the search, scores it on the caller's own
# folds, and applies one of three decision rules. Two properties of the
# result matter more than the number it returns. The standard error is
# aggregated at the resampling unit, not at the row, because rows inside a
# cluster are not independent replicates and a row-level SE understates
# the spread by roughly the square root of the cluster size. And the
# projected draws are labelled `projected_parameters`: they are the
# reference posterior projected onto a submodel, not a posterior anyone
# sampled, and no consumer may read them as one.


# --- the structure gate ------------------------------------------------ #

# How projpred 2.10.0 answered for each structure, and where.
#
# The verdicts that come from runs are marked (probe Cn); the rest refuse
# as unrun.
.fb_sel_projpred_random_support <- function(type) {
  supported <- function(reason) list(verdict = "supported", reason = reason)
  refused <- function(reason) list(verdict = "refused", reason = reason)
  unrun <- function(what) {
    list(
      verdict = "unrun",
      reason = paste0(
        "no gate has been executed for ", what, " under projection, and an ",
        "unrun cell refuses"
      )
    )
  }
  switch(
    type,
    simple = supported(paste0(
      "grouped random intercept: get_refmodel() accepted it and a ",
      "validated grouped k-fold search ran (probes C1-C3)"
    )),
    simple_slope_cor = supported(paste0(
      "grouped random slope: get_refmodel() accepted the grouped-term ",
      "construct (probes C1-C3)"
    )),
    # Executed, and it corrects the row the contract carried. projpred
    # 2.10.0 refuses the `||` spelling outright -- "Unexpected number of
    # `|` characters in group terms" -- on a flexyBayes fit and on a native
    # brms fit alike, while the same model written as two terms is
    # accepted. The refusal is about the spelling, not about the model, so
    # the reason says which spelling works.
    simple_slope_uncor = refused(paste0(
      "projpred 2.10.0 refuses the `||` spelling with \"Unexpected number ",
      "of `|` characters in group terms\"; the same model written as ",
      "(1 | g) + (0 + x | g) is accepted"
    )),
    vm = refused(paste0(
      "a known-covariance carrier lowers to (1 | gr(g, cov = K)), and ",
      "projpred refuses a function call on the right-hand side of a ",
      "group-term `|` (probe C5)"
    )),
    ped = refused(paste0(
      "a pedigree carrier lowers to (1 | gr(g, cov = A)), and projpred ",
      "refuses a function call on the right-hand side of a group-term ",
      "`|` (probe C5)"
    )),
    us = unrun("a us() translation"),
    us_gxe = unrun("a us() translation"),
    fa = unrun("an fa() translation"),
    fa_gxe = unrun("an fa() translation"),
    corh = unrun("a corh() translation"),
    corh_gxe = unrun("a corh() translation"),
    at_simple = unrun("a heterogeneous variance structure"),
    ar1 = unrun("an autoregressive field"),
    ar1_spatial = unrun("an autoregressive field"),
    smooth_mgcv = unrun("a smooth term"),
    spline = unrun("a smooth term"),
    unrun(paste0("the random-term type `", type, "`"))
  )
}

.fb_sel_projpred_residual_support <- function(type) {
  if (type %in% c("units", "ide", "id", "idv")) {
    return(list(
      verdict = "supported",
      reason = paste0(
        "an independent homogeneous residual is the family's own error term"
      )
    ))
  }
  list(
    verdict = "unrun",
    reason = paste0(
      "no gate has been executed for the residual structure `", type,
      "` under projection, and an unrun cell refuses"
    )
  )
}

.fb_sel_gate_row <- function(feature, support) {
  data.frame(
    feature = feature,
    verdict = support$verdict,
    reason = support$reason,
    stringsAsFactors = FALSE
  )
}

# Term labels of a fitted brms model, with a group term's parentheses
# restored: `stats::terms()` reports `1 | g` and projpred names the same
# term `(1 | g)` in a ranking and expects that spelling in `search_terms`.
.fb_sel_reference_terms <- function(brmsfit) {
  form <- tryCatch(stats::formula(brmsfit), error = function(e) NULL)
  if (is.null(form)) {
    return(character(0))
  }
  fixed <- if (inherits(form, "brmsformula")) form$formula else form
  labels <- tryCatch(
    attr(stats::terms(fixed), "term.labels"),
    error = function(e) character(0)
  )
  vapply(
    labels,
    function(l) if (grepl("|", l, fixed = TRUE)) paste0("(", l, ")") else l,
    character(1L),
    USE.NAMES = FALSE
  )
}

# Is this brms formula a non-linear partition? Read from the formula
# object the fit carries, not from the call that produced it.
.fb_sel_brms_is_nonlinear <- function(brmsfit) {
  form <- tryCatch(stats::formula(brmsfit), error = function(e) NULL)
  if (is.null(form)) {
    return(FALSE)
  }
  isTRUE(attr(form$formula, "nl")) || length(form$pforms %||% list()) > 0L
}

# .fb_sel_structure_gate() --- the section 6.3 matrix, for THIS model.
#
# One row per structural feature the model actually carries, plus the
# family row, so the table describes the model in front of it rather than
# reciting a matrix. Accepts a flexyBayes IR, a fitted flexyBayes model, or
# a bare `brmsfit` --- the last because two cells of the matrix (a
# known-covariance carrier and a non-linear partition) can only be executed
# against models that flexyBayes either lowers to brms or cannot build at
# all, and a gate whose red cells cannot be exercised is not a gate.
.fb_sel_structure_gate <- function(fb_or_fit) {
  rows <- list()

  if (inherits(fb_or_fit, "brmsfit")) {
    fam <- tryCatch(stats::family(fb_or_fit), error = function(e) NULL)
    fam_name <- fam$family %||% NA_character_
    link_name <- fam$link %||% NA_character_
    rows[[length(rows) + 1L]] <- .fb_sel_gate_row(
      "response family",
      if (identical(fam_name, "gaussian") && identical(link_name, "identity")) {
        list(
          verdict = "supported",
          reason = "Gaussian response with an identity link"
        )
      } else {
        list(
          verdict = "refused",
          reason = paste0(
            "projection is Gaussian-only in this release, and this reference ",
            "is ", fam_name, "/", link_name
          )
        )
      }
    )
    if (.fb_sel_brms_is_nonlinear(fb_or_fit)) {
      rows[[length(rows) + 1L]] <- .fb_sel_gate_row(
        "non-linear parameter partition",
        list(
          verdict = "refused",
          reason = "projpred does not support non-linear models (probe C6)"
        )
      )
    }
    for (label in .fb_sel_reference_terms(fb_or_fit)) {
      if (!grepl("|", label, fixed = TRUE)) {
        next
      }
      group <- sub("^.*\\|+\\s*", "", sub("\\)$", "", label))
      rows[[length(rows) + 1L]] <- .fb_sel_gate_row(
        label,
        if (grepl("(", group, fixed = TRUE)) {
          list(
            verdict = "refused",
            reason = paste0(
              "projpred refuses a function call on the right-hand side of ",
              "a group-term `|` (probe C5)"
            )
          )
        } else if (grepl("||", label, fixed = TRUE)) {
          .fb_sel_projpred_random_support("simple_slope_uncor")
        } else {
          list(
            verdict = "supported",
            reason = paste0(
              "grouped random term: get_refmodel() accepted it and a ",
              "validated grouped k-fold search ran (probes C1-C3)"
            )
          )
        }
      )
    }
    return(do.call(rbind, rows))
  }

  ir <- if (inherits(fb_or_fit, "fb_terms")) {
    fb_or_fit
  } else {
    .fb_fit_ir(fb_or_fit)
  }
  if (is.null(ir)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_shape_unreadable",
      message = paste0(
        "fb_select(): this fit does not carry the parsed model the ",
        "projection structure gate reads."
      )
    ))
  }

  fam_name <- ir$family %||% NA_character_
  link_name <- ir$link %||% NA_character_
  rows[[length(rows) + 1L]] <- .fb_sel_gate_row(
    "response family",
    if (identical(fam_name, "gaussian") && identical(link_name, "identity")) {
      list(
        verdict = "supported",
        reason = "Gaussian response with an identity link"
      )
    } else {
      list(
        verdict = "refused",
        reason = paste0(
          "projection is Gaussian-only in this release, and this reference ",
          "is ", fam_name, "/", link_name
        )
      )
    }
  )
  rows[[length(rows) + 1L]] <- .fb_sel_gate_row(
    "population-level effects",
    list(
      verdict = "supported",
      reason = "the fixed part is what the search selects over (probe C2)"
    )
  )

  random_labels <- .fb_sel_random_labels(ir)
  types <- vapply(
    ir$random_terms %||% list(),
    function(t) t$type %||% "unknown",
    character(1L)
  )
  for (i in seq_along(types)) {
    rows[[length(rows) + 1L]] <- .fb_sel_gate_row(
      random_labels[[i]],
      .fb_sel_projpred_random_support(types[[i]])
    )
  }

  residual_types <- vapply(
    ir$residual_terms %||% list(),
    function(t) t$type %||% "unknown",
    character(1L)
  )
  for (type in residual_types) {
    rows[[length(rows) + 1L]] <- .fb_sel_gate_row(
      paste0("residual: ", type),
      .fb_sel_projpred_residual_support(type)
    )
  }

  do.call(rbind, rows)
}

# .fb_sel_require_projectable() --- the same table, as a gate.
#
# The family is checked first and refuses under its own code, because
# "this response is not Gaussian" and "this covariance structure has no
# executed gate" are different answers and a caller acts differently on
# each.
.fb_sel_require_projectable <- function(fb_or_fit) {
  gate <- .fb_sel_structure_gate(fb_or_fit)
  family_row <- gate[gate$feature == "response family", , drop = FALSE]
  unsupported_family <- nrow(family_row) > 0L &&
    !identical(family_row$verdict[[1L]], "supported")
  if (unsupported_family) {
    stop(.fb_refusal_condition(
      reason_code = "projection_family_unsupported",
      message = paste0("fb_select(): ", family_row$reason[[1L]], "."),
      gate = gate
    ))
  }
  bad <- gate[gate$verdict != "supported", , drop = FALSE]
  if (nrow(bad) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "projection_structure_unsupported",
      message = paste0(
        "fb_select(): projection does not support ",
        paste(
          paste0("`", bad$feature, "` (", bad$reason, ")"),
          collapse = "; "
        ),
        ". Pass decision = \"full\" to keep the shrinkage posterior with no ",
        "selection."
      ),
      gate = gate
    ))
  }
  invisible(gate)
}


# --- folds ------------------------------------------------------------- #

# .fb_sel_projpred_folds() --- an <fb_sel_folds> as projpred's fold vector.
#
# projpred takes one integer per ROW naming that row's assessment fold;
# `.fb_sel_folds()` returns the complementary shape, a list of splits. The
# unit key travels with the vector as an attribute because the scoring unit
# -- the cluster -- is not recoverable from the fold numbers: several
# clusters share a fold.
.fb_sel_projpred_folds <- function(folds) {
  if (is.numeric(folds) && is.null(dim(folds))) {
    out <- as.integer(folds)
    if (anyNA(out) || any(out < 1L)) {
      stop(
        ".fb_sel_project(): `folds` must name a fold from 1 upwards for ",
        "every row.",
        call. = FALSE
      )
    }
    return(out)
  }
  if (!inherits(folds, "fb_sel_folds")) {
    stop(
      ".fb_sel_project(): `folds` must be an <fb_sel_folds> or an integer ",
      "vector over the rows.",
      call. = FALSE
    )
  }
  repeats <- unique(vapply(folds, function(s) s$repeat_id, integer(1L)))
  if (length(repeats) > 1L) {
    stop(
      ".fb_sel_project(): projpred's k-fold cross-validation takes one ",
      "partition, and this resampling scheme carries ", length(repeats),
      " repeats. Pass a single-repeat scheme, or score the repeats ",
      "separately.",
      call. = FALSE
    )
  }
  n <- attr(folds, "n_rows")
  out <- integer(n)
  for (i in seq_along(folds)) {
    out[folds[[i]]$assess_rows] <- i
  }
  if (any(out == 0L)) {
    stop(
      ".fb_sel_project(): ", sum(out == 0L), " row(s) are in no assessment ",
      "fold.",
      call. = FALSE
    )
  }
  attr(out, "unit_key") <- attr(folds, "unit_key")
  attr(out, "unit") <- attr(folds, "unit")
  out
}


# --- prediction target --------------------------------------------------- #

# .fb_sel_predict_cluster() --- the two prediction targets of section 6.3.
#
# `known_cluster` conditions on the fitted level: the cluster's own random
# effect is in the posterior and the prediction uses it. `new_cluster`
# treats the level as unseen and draws its random effect from the estimated
# population, which is `allow_new_levels = TRUE` with
# `sample_new_levels = "gaussian"` (probe C7). The two answers differ in
# spread, not usually in location: an unseen cluster carries the
# random-effect variance that a known one has already resolved.
#
# The caller supplies `newdata`. Under `new_cluster` the grouping level in
# it is expected to be one the fit has not seen -- which is what an
# assessment fold's held-out cluster is.
.fb_sel_predict_cluster <- function(fit,
                                    newdata,
                                    prediction = c(
                                      "new_cluster", "known_cluster"
                                    ),
                                    ndraws = NULL,
                                    seed = NULL) {
  prediction <- match.arg(prediction)
  brmsfit <- if (inherits(fit, "brmsfit")) fit else fit$brms
  if (is.null(brmsfit)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_shape_unreadable",
      message = paste0(
        "fb_select(): the predictive draws are taken from the reference ",
        "model's brms fit, and this object carries none."
      )
    ))
  }
  args <- list(object = brmsfit, newdata = newdata)
  if (!is.null(ndraws)) {
    args$ndraws <- as.integer(ndraws)
  }
  if (identical(prediction, "new_cluster")) {
    args$allow_new_levels <- TRUE
    args$sample_new_levels <- "gaussian"
  }
  draw <- function() do.call(brms::posterior_epred, args)
  if (is.null(seed)) draw() else .fb_sel_with_seed(seed, draw())
}


# --- unit-level scoring -------------------------------------------------- #

# .fb_sel_unit_elpd() --- elpd and its standard error at the scoring unit.
#
# projpred keeps the pointwise log predictive density of every observation
# under every submodel size in `$summaries`, and under the reference model
# in `$summaries$ref`. Summing those within a cluster and taking the spread
# ACROSS clusters is the only aggregation that matches the resampling: the
# rows of one cluster are not independent replicates, and the row-level
# standard error projpred's own summary reports is smaller than the
# cluster-level one by roughly the square root of the cluster size.
#
# The difference to the reference is computed as a PAIRED difference per
# unit before its spread is taken. Between-unit variance cancels in the
# pairing, which is what makes the difference the better-identified
# quantity and the one the decision rules use.
.fb_sel_unit_elpd <- function(summaries, unit_key = NULL) {
  ref_lppd <- as.numeric(summaries$ref$lppd)
  n <- length(ref_lppd)
  if (is.null(unit_key) || length(unit_key) != n) {
    unit_key <- as.character(seq_len(n))
  }
  unit_key <- as.character(unit_key)
  unit <- factor(unit_key, levels = unique(unit_key))
  n_units <- nlevels(unit)

  by_unit <- function(x) as.numeric(tapply(x, unit, sum))
  se_of <- function(x) {
    if (length(x) < 2L) NA_real_ else stats::sd(x) * sqrt(length(x))
  }

  ref_units <- by_unit(ref_lppd)
  sizes <- seq_along(summaries$sub) - 1L
  n_nonfinite <- sum(!is.finite(ref_lppd))

  rows <- vector("list", length(summaries$sub))
  unit_matrix <- matrix(NA_real_, nrow = n_units, ncol = length(summaries$sub))
  for (k in seq_along(summaries$sub)) {
    lppd <- as.numeric(summaries$sub[[k]]$lppd)
    n_nonfinite <- n_nonfinite + sum(!is.finite(lppd))
    u <- by_unit(lppd)
    unit_matrix[, k] <- u
    paired <- u - ref_units
    rows[[k]] <- data.frame(
      size = sizes[[k]],
      elpd = sum(u),
      se = se_of(u),
      delta_to_ref = sum(paired),
      delta_se = se_of(paired),
      stringsAsFactors = FALSE
    )
  }

  list(
    table = do.call(rbind, rows),
    unit_matrix = unit_matrix,
    reference_units = ref_units,
    n_units = n_units,
    unit_labels = levels(unit),
    n_nonfinite = as.integer(n_nonfinite)
  )
}

# projpred's own standard error, for the same sizes, kept beside the
# unit-level one. It is derived from the pointwise densities treating every
# ROW as an independent replicate, so on clustered data it is the smaller
# of the two by roughly the square root of the cluster size. Recording both
# is what lets a reader see which ruler a decision was made with.
.fb_sel_row_level_se <- function(vs, sizes) {
  perf <- tryCatch(
    summary(vs, stats = "elpd")$perf_sub,
    error = function(e) NULL
  )
  if (!is.data.frame(perf) || !all(c("size", "elpd.se") %in% names(perf))) {
    return(rep(NA_real_, length(sizes)))
  }
  perf$elpd.se[match(sizes, perf$size)]
}


# --- the decision -------------------------------------------------------- #

.fb_sel_decision_spec <- function(decision = NULL) {
  spec <- list(
    method = "projection",
    metric = "elpd",
    rule = "one_se",
    tolerance = NULL,
    search = "forward",
    nterms_max = NULL,
    validate_search = TRUE,
    threshold = 0.5,
    # The practical effect size the posterior-threshold decision compares
    # against. NULL here rather than 0: `fb_decision()` refuses to build a
    # posterior-threshold decision without one, and a zero default would
    # make P(|beta| > 0) the quantity, which is one under any continuous
    # prior whatever the data say.
    effect_size = NULL
  )
  if (!is.null(decision)) {
    if (!is.list(decision)) {
      stop(
        ".fb_sel_project(): `decision` must be a list (an <fb_decision>).",
        call. = FALSE
      )
    }
    for (nm in intersect(names(decision), names(spec))) {
      spec[[nm]] <- decision[[nm]]
    }
  }
  spec$rule <- match.arg(spec$rule, c("one_se", "min", "tolerance"))
  spec$search <- match.arg(spec$search, c("forward", "L1"))
  spec$validate_search <- isTRUE(spec$validate_search)
  if (identical(spec$rule, "tolerance")) {
    if (
      is.null(spec$tolerance) ||
        !is.numeric(spec$tolerance) ||
        length(spec$tolerance) != 1L ||
        !is.finite(spec$tolerance) ||
        spec$tolerance < 0
    ) {
      stop(
        ".fb_sel_project(): rule = \"tolerance\" needs a single ",
        "non-negative `tolerance` in nats per unit.",
        call. = FALSE
      )
    }
  }
  spec
}

# .fb_sel_decide() --- the three rules, on the unit-level table.
#
#   min        the size with the highest elpd
#   one_se     the smallest size whose loss against that best size is
#              within one standard error OF THAT PAIRED DIFFERENCE. The
#              paired SE, not the best size's own SE: the two models are
#              scored on the same units, so the between-unit variance
#              cancels and the un-paired SE would be the wrong ruler
#   tolerance  the smallest size losing less than `tolerance` nats per
#              unit against the REFERENCE model, which is the quantity
#              section 6.1 names as the one the rule operates on
#
# `floor_size` is the smallest submodel that still contains every forced
# term. A smaller size drops a protected term and is not an admissible
# answer whatever its elpd.
.fb_sel_decide <- function(unit_elpd, spec, floor_size = 0L) {
  tab <- unit_elpd$table
  admissible <- tab$size >= floor_size
  if (!any(admissible)) {
    admissible <- rep(TRUE, nrow(tab))
  }
  idx <- which(admissible)
  best <- idx[which.max(tab$elpd[idx])]

  if (identical(spec$rule, "min")) {
    return(list(size = tab$size[[best]], row = best, rule = "min"))
  }

  if (identical(spec$rule, "tolerance")) {
    loss_per_unit <- -tab$delta_to_ref / unit_elpd$n_units
    within_tol <- is.finite(loss_per_unit[idx]) &
      loss_per_unit[idx] < spec$tolerance
    ok <- idx[within_tol]
    if (length(ok) == 0L) {
      return(list(size = tab$size[[best]], row = best, rule = "tolerance"))
    }
    pick <- ok[which.min(tab$size[ok])]
    return(list(size = tab$size[[pick]], row = pick, rule = "tolerance"))
  }

  # one_se: the paired spread of (best size - candidate size) across units.
  best_units <- unit_elpd$unit_matrix[, best]
  within <- vapply(
    idx,
    function(k) {
      paired <- best_units - unit_elpd$unit_matrix[, k]
      loss <- sum(paired)
      se <- if (length(paired) < 2L) {
        NA_real_
      } else {
        stats::sd(paired) * sqrt(length(paired))
      }
      isTRUE(loss <= (if (is.finite(se)) se else 0))
    },
    logical(1L)
  )
  ok <- idx[within]
  if (length(ok) == 0L) {
    return(list(size = tab$size[[best]], row = best, rule = "one_se"))
  }
  pick <- ok[which.min(tab$size[ok])]
  list(size = tab$size[[pick]], row = pick, rule = "one_se")
}


# --- the search ---------------------------------------------------------- #

# `search_terms` that force every protected and random term into every
# submodel. projpred combines the elements of this vector, so an element
# carrying the forced block plus one free term is the smallest step the
# search can take, and no path exists that omits the block (probe C4).
.fb_sel_search_terms <- function(forced, free) {
  forced <- unique(forced[nzchar(forced)])
  free <- setdiff(unique(free[nzchar(free)]), forced)
  if (length(forced) == 0L) {
    return(free)
  }
  block <- paste(forced, collapse = " + ")
  c(block, paste0(block, " + ", free))
}

# .fb_sel_project() --- the projection predictive decision.
#
# Everything the object records is read off the run: the settings from the
# `vsel` object rather than from the arguments, the rankings and their
# fold-wise agreement from `ranking()` and `cv_proportions()`, the elpd
# table from the pointwise summaries aggregated at the unit, and the
# projected draws with the label that keeps them from being read as the
# reference posterior.
.fb_sel_project <- function(fit,
                            folds,
                            protected_terms = NULL,
                            random_terms = NULL,
                            decision = NULL,
                            seed = NULL,
                            ndraws_project = 400L,
                            ...) {
  .check_installed(
    "projpred",
    "Package 'projpred' is required for the projection decision. ",
    "Install via:\n  install.packages('projpred')\n"
  )
  if (
    is.null(seed) ||
      length(seed) != 1L ||
      !is.numeric(seed) ||
      !is.finite(seed)
  ) {
    stop(.fb_refusal_condition(
      reason_code = "selection_seed_required",
      message = paste0(
        "fb_select(): the projection search and its projected draws are ",
        "stochastic, and no `seed` was given."
      )
    ))
  }

  brmsfit <- if (inherits(fit, "brmsfit")) fit else fit$brms
  if (is.null(brmsfit)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_shape_unreadable",
      message = paste0(
        "fb_select(): the projection reads the reference model's brms fit, ",
        "and this object carries none."
      )
    ))
  }
  .fb_sel_require_projectable(fit)

  spec <- .fb_sel_decision_spec(decision)
  if (!identical(spec$metric, "elpd")) {
    stop(
      ".fb_sel_project(): the projection decision is on elpd in this ",
      "release; `", spec$metric, "` is scored by the out-of-fold path, not ",
      "by the projection.",
      call. = FALSE
    )
  }

  fold_vector <- .fb_sel_projpred_folds(folds)
  unit_key <- attr(fold_vector, "unit_key")
  scoring_unit <- attr(fold_vector, "unit") %||% character(0)
  if (length(scoring_unit) == 0L) {
    scoring_unit <- "row"
  }

  all_terms <- .fb_sel_reference_terms(brmsfit)
  forced <- unique(c(
    protected_terms %||% attr(fit, "protected_terms") %||% character(0),
    random_terms %||% attr(fit, "random_terms") %||% character(0)
  ))
  forced <- intersect(all_terms, forced)
  free <- setdiff(all_terms, forced)
  search_terms <- .fb_sel_search_terms(forced, free)
  nterms_max <- spec$nterms_max %||% (length(forced) + length(free))

  t0 <- proc.time()
  # The seed goes to projpred AND around the whole call, and the reference
  # model is BUILT inside that seeded region. `get_refmodel()` draws from
  # R's stream itself; building it outside left two runs of the same seeded
  # call with rankings that agreed and elpd values that did not, differing
  # by up to a tenth of a nat. Measured, not assumed: with the construction
  # moved inside, three successive runs return the same table to the last
  # digit.
  vs <- .fb_sel_with_seed(seed, {
    refmodel <- projpred::get_refmodel(brmsfit)
    projpred::cv_varsel(
      refmodel,
      cv_method = "kfold",
      K = max(fold_vector),
      folds = fold_vector,
      method = spec$search,
      nterms_max = nterms_max,
      validate_search = spec$validate_search,
      search_terms = search_terms,
      verbose = FALSE,
      seed = as.integer(seed),
      ...
    )
  })
  t_search <- unname((proc.time() - t0)["elapsed"])

  rk <- projpred::ranking(vs)
  stability <- tryCatch(
    projpred::cv_proportions(rk),
    error = function(e) NULL
  )

  unit_elpd <- .fb_sel_unit_elpd(vs$summaries, unit_key)
  unit_elpd$table$se_row_level <- .fb_sel_row_level_se(vs, unit_elpd$table$size)
  # projpred measures a submodel's SIZE in search steps, and the forced
  # block is one step however many terms it contains (its ranking entry is
  # the whole block, `p1 + (1 | site)`). The smallest admissible submodel is
  # therefore the one step that carries the block, not one step per forced
  # term.
  floor_size <- if (length(forced) > 0L) 1L else 0L
  chosen <- .fb_sel_decide(unit_elpd, spec, floor_size = floor_size)
  ranking_full <- rk$fulldata
  chosen_entries <- if (chosen$size > 0L) {
    utils::head(ranking_full, chosen$size)
  } else {
    character(0)
  }
  # The forced block is one ranking entry carrying several terms; every
  # other entry is a single term. Only the block is expanded, by identity
  # against the string this function built -- splitting every entry on "+"
  # would take `I(a + b)` apart as well.
  forced_block <- if (length(forced) > 0L) {
    paste(forced, collapse = " + ")
  } else {
    ""
  }
  chosen_terms <- unique(unlist(
    lapply(
      chosen_entries,
      function(e) if (identical(e, forced_block)) forced else e
    ),
    use.names = FALSE
  ))

  t0 <- proc.time()
  prj <- .fb_sel_with_seed(seed, projpred::project(
    vs,
    nterms = chosen$size,
    ndraws = as.integer(ndraws_project),
    seed = as.integer(seed)
  ))
  projected <- as.matrix(prj)
  t_project <- unname((proc.time() - t0)["elapsed"])
  attr(projected, "fb_draw_label") <- "projected_parameters"

  list(
    ranking = ranking_full,
    foldwise_rankings = rk$foldwise,
    rank_stability = stability,
    elpd_by_size = unit_elpd$table,
    chosen_size = chosen$size,
    chosen_terms = chosen_terms,
    forced_terms = forced,
    projected_draws = list(
      draws = projected,
      label = "projected_parameters",
      nterms = chosen$size,
      ndraws = nrow(projected)
    ),
    scoring = list(
      unit = scoring_unit,
      n_units = unit_elpd$n_units,
      n_nonfinite_lppd = unit_elpd$n_nonfinite,
      reference_elpd = sum(unit_elpd$reference_units),
      se_basis = paste0(
        "pointwise log predictive densities summed within the ",
        "resampling unit"
      )
    ),
    settings = list(
      search = vs$method,
      nterms_max = vs$nterms_max,
      validate_search = vs$validate_search,
      optimistic = !isTRUE(vs$validate_search),
      cv_method = vs$cv_method,
      K = vs$K,
      rule = spec$rule,
      tolerance = spec$tolerance,
      search_terms = search_terms,
      seed = as.integer(seed),
      projpred_version = as.character(vs$projpred_version)
    ),
    timing = list(cv_varsel = t_search, project = t_project)
  )
}
