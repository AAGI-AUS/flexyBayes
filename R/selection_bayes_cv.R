# The unpenalised Bayesian reference --- a per-fold engine refit, scored on
# posterior predictive draws.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 2, 3.4, 7, 9.2,
# 10 and 14 (Unit C).
#
# Why this file exists beside R/selection_resamples.R. The B0 route
# (`.fb_sel_cv()`) fits a closed-form GLS per fold on variance components
# estimated by an independent REML plugin inside that fold. It is correct and
# it stays, but it plugs a covariance in rather than integrating over it, so
# it is a benchmark and not a Bayesian answer. Section 2 of the contract
# separates the two on four axes, and this file owns the `bayes` workflow at
# `penalty = NULL`: the base model's own engine and priors are refitted on
# every fold's training rows, twice --- protected-only and protected plus
# candidates --- and the held-out rows are scored on draws from the posterior
# predictive rather than on a plug-in mean.
#
# Three things the route must get right, and where each is enforced.
#
#   The estimand.  `covariance_estimand = "integrated"` is a property of the
#   fit, not of an argument. It is set here and read from the fold object by
#   anything that certifies the run, which is the same rule section 9.2
#   applies to the benchmark route's `reml_plugin_covariance`.
#
#   The prediction target.  `new_cluster` must not condition on the held-out
#   unit's random effect. Because a fold may split a unit (folds are on units
#   only when `fb_resamples(unit = )` says so), "absent from training" is
#   enforced rather than assumed: every random-effect grouping column is
#   relabelled on the assessment rows to a level no training row carries, so
#   the engine has nothing to condition on and the predictive picks up the
#   random-effect variance instead. Under `known_cluster` the labels are left
#   alone and the engine conditions on whatever training rows share the level.
#   This is the Bayesian counterpart of B0's `V_ST = 0` rule.
#
#   The fold boundary.  Everything the fold learns is declared in a `learned`
#   list and held to that claim by `.fb_sel_leakage_gate()` --- the same
#   falsifier the benchmark route uses, run against this route's builder.
#
# --- what was measured on this machine, not recalled -------------------- #
#
# R 4.6.1, INLA 26.8.7, brms 2.23.0 / CmdStan via cmdstanr 0.9.0.
#
# `INLA::inla.posterior.sample(n, result)` returns a list of length `n`. Each
# element is a plain list of three components:
#
#   $hyperpar  named numeric --- for a Gaussian mixed model, "Precision for
#              the Gaussian observations" and one "Precision for <group>" per
#              random term;
#   $latent    a one-column matrix whose rownames run "Predictor:1" ...
#              "Predictor:n" for the n data rows, then "<group>:1" ...
#              "<group>:q" for each random term's levels, then
#              "(Intercept):1" and one row per fixed-effect coefficient;
#   $logdens   list(hyperpar, latent, joint).
#
# `result$misc$configs$contents` carries the same blocking as `$tag`,
# `$start`, `$length`. The "Predictor:" block is in data-row order and has
# exactly `nrow(data)` entries, which is what lets an assessment row index
# into it directly. A row whose response is `NA` keeps its Predictor entry, so
# the NA-response refit is the INLA prediction route: measured on a twelve-
# site model with two sites held out, the held-out sites' random effects came
# back with mean near zero and a draw-to-draw standard deviation of 1.32 and
# 1.41 against 0.37 and 0.43 for observed sites, which is the marginal rather
# than the conditional effect.
#
# `INLA::inla.posterior.sample()` is NOT reproducible on this build. Measured:
# two calls with `seed = 11L, num.threads = "1:1"` differ by up to 0.21 in the
# latent vector, with `parallel.configs = FALSE` by 0.21, with
# `skew.corr = FALSE` by 0.28, in `inla.mode = "classic"` by 0.25, and with
# the fit itself run at `num.threads = "1:1"` by 0.21; seeding R and passing
# `seed = 0L` differs by 1.14. The underlying `INLA::inla.qsample(seed = 7L)`
# IS bit-reproducible on the same build, so the loss happens above it. The
# seed is passed anyway --- it is the documented contract --- and the INLA
# arm's reproducibility is therefore stated at the score level within Monte
# Carlo error rather than claimed at the draw level. The INLA optimiser is
# itself only near-deterministic: two identical calls agreed to 3.7e-08 on the
# fixed-effect means, not bitwise.
#
# `brms::posterior_epred(newdata, allow_new_levels = TRUE,
# sample_new_levels = "gaussian", draw_ids = )` IS bit-reproducible under
# `set.seed()`, and a brms refit under the same `seed` and `chains` returns a
# bit-identical draws matrix. A new level's predictive mean carried a
# draw-to-draw standard deviation of 1.89 against 0.18 for a fitted level.
#
# Where the time goes, and why `n_draws` defaults high. On a 120-row,
# twelve-site model a three-fold run cost 12 s of INLA fitting and 76 s of
# `inla.posterior.sample()`, and the sampling cost was flat between 400 and
# 6,000 draws --- it is a per-call overhead, not a per-draw one. Draws are
# therefore close to free, and they are worth taking: under `new_cluster` the
# held-out unit's random effect is integrated by Monte Carlo, so when the
# cluster scale exceeds the residual scale only the draws that land within a
# residual standard deviation of the observation contribute to that row's
# density and the effective sample size is a fraction of `n_draws`. The
# per-unit elpd moved by about 0.3 nats between 1,500 and 6,000 draws on that
# model. This is a property of draw-based predictive scoring rather than of
# either engine: `brms::posterior_predict()` integrates the same effect the
# same way.

# --- the base's own specification -------------------------------------- #

# .fb_sel_bayes_normalise_base() --- one shape for two engines.
#
# The two backends do not return the same object. An INLA fit carries the
# parsed model at `$fb` and the fitted data at `$data`; a brms fit carries the
# same parsed model at `$extras$fb_terms` and NO data at all --- its
# `$extras$fb_dataset` is metadata only (`origin = "metadata-only"`, `$data`
# NULL), and `$brms$data` is the model frame, which holds the model's own
# variables and nothing else. A candidate block names columns the base model
# does not, so on the brms path the data has to come from the caller.
#
# This function puts both shapes into the one the selection substrate reads,
# on a copy, so `.fb_select_gate()` and everything downstream can be reused
# rather than re-written per engine.
.fb_sel_bayes_normalise_base <- function(base, data = NULL) {
  fb <- base$fb %||% base$extras$fb_terms
  if (is.null(fb)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_shape_unreadable",
      message = paste0(
        "fb_select(): this base fit carries neither `$fb` nor ",
        "`$extras$fb_terms`, so the parsed model the selection substrate ",
        "reads is not on it. Fits on the aggregated route, which an ",
        "intercept-only fixed model takes, have that shape."
      )
    ))
  }
  resolved <- data %||%
    base$data %||%
    base$extras$fb_dataset$data %||%
    base$brms$data
  if (is.null(resolved) || !is.data.frame(resolved)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_base_shape_unreadable",
      message = paste0(
        "fb_select(): this base fit does not retain the data it was fitted ",
        "to -- the brms path records the dataset as metadata only -- and no ",
        "`data` was given. Pass the data frame the base was fitted to; the ",
        "candidate block names columns that are not in the base model, so ",
        "the engine's own model frame cannot stand in for it."
      )
    ))
  }
  base$fb <- fb
  base$data <- .fb_sel_with_weights_column(base, as.data.frame(resolved))
  base
}


# .fb_sel_grouping_columns() --- the data columns each random term's effect
# is indexed by, read off the parsed model the emitters read.
#
# A column here is one both emitters build the term's index from: the group
# of an intercept or a slope (`var`; a slope's covariate is not a grouping),
# every factor of an interaction (`outer`, `inner`, `vars`), and the inner
# factor of a structure written over an outer one (`diag(f):g`, `us(f):g`,
# `fa(f, k):g`). The outer factor there selects which variance or covariance
# block applies, and is not relabelled: given a level no training row has,
# brms refuses it as a new level of the slope factor and INLA builds a new,
# empty field for it.
#
# A smooth of a covariate (`spl(x)`, `s(x)`, `gp(x)`) and a residual term
# (`units`, `at(f):units`) index no cluster, so they carry no column and
# nothing to relabel. `unrelabelled` names every other term with no column,
# so that no structure passes through `new_cluster` with its levels
# untouched: a known-covariance term (`vm()`, `ped()`) indexes its matrix by
# level, and a new level is a row the matrix does not have; an
# autoregressive field with no id() factor indexes one lattice across all
# rows; and a type this function does not list is named in the refusal too.
.fb_sel_grouping_columns <- function(fb) {
  columns <- character(0)
  unrelabelled <- character(0)
  # RANGE and ROW of id(TRIAL):ar1(RANGE):ar1(ROW) are plot positions. In
  # TRIAL:RANGE, relabelling TRIAL already gives the held-out rows new levels.
  field_group <- character(0)
  field_index <- character(0)
  for (t in fb$random_terms %||% list()) {
    if (identical(t$type, "ar1_spatial") && length(t$at_var) > 0L) {
      field_group <- c(field_group, t$at_var)
      field_index <- c(field_index, t$row_var, t$col_var)
    }
  }
  for (t in fb$random_terms %||% list()) {
    type <- t$type %||% "unknown"
    if (!nzchar(type)) {
      type <- "unknown"
    }
    found <- switch(
      type,
      simple = ,
      ide = ,
      id = ,
      simple_slope_uncor = ,
      simple_slope_cor = t$var,
      nested = ,
      combo = {
        vars <- .fb_nested_vars(t)
        if (any(vars %in% field_group)) setdiff(vars, field_index) else vars
      },
      at_simple = ,
      us_gxe = ,
      fa_gxe = ,
      corh_gxe = t$inner,
      # A field written id(g):ar1(r):ar1(c) is one realisation per level of
      # g with shared parameters, so a held-out level is a new, independent
      # realisation. Without id(g) the field has no column to relabel.
      ar1_spatial = t$at_var,
      spline = ,
      smooth_mgcv = ,
      smoother = ,
      gp = ,
      units = ,
      at_units = character(0),
      NULL
    )
    if (is.null(found)) {
      named <- unlist(t[c("var", "outer", "inner", "row_var", "col_var")])
      unrelabelled <- c(
        unrelabelled,
        paste0(type, "(", paste(named, collapse = ", "), ")")
      )
    } else {
      columns <- c(columns, found)
    }
  }
  list(columns = unique(columns), unrelabelled = unique(unrelabelled))
}


# .fb_sel_require_relabelable() --- refuse, before anything is fitted, a
# `new_cluster` request on a base carrying a term that cannot take a new
# level. `known_cluster` keeps every level and passes.
.fb_sel_require_relabelable <- function(fb, prediction) {
  if (!identical(prediction, "new_cluster")) {
    return(invisible(NULL))
  }
  unrelabelled <- .fb_sel_grouping_columns(fb)$unrelabelled
  if (length(unrelabelled) == 0L) {
    return(invisible(NULL))
  }
  stop(.fb_refusal_condition(
    reason_code = "selection_new_cluster_not_representable",
    message = paste0(
      "fb_select(): prediction = \"new_cluster\" gives each assessment row ",
      "a grouping level no training row carries, and these random terms ",
      "cannot take a new level without changing the structure they index: ",
      paste(unrelabelled, collapse = ", "),
      ". A known-covariance term indexes its matrix by level, and an ",
      "autoregressive field indexes its lattice. Use prediction = ",
      "\"known_cluster\", which keeps the levels, or a base without the term."
    ),
    terms = unrelabelled
  ))
}


# .fb_sel_require_fold_local_formula() --- refuse, before anything is fitted,
# a formula variable whose value on a row depends on other rows, on a route
# that refits the formula text on each fold.
#
# The engines evaluate that text over the rows they are handed: an INLA fold
# refit over the training and assessment rows together, so `scale(z1)` there
# is centred with the assessment rows included. Each fold of `folds` is
# checked, protected part and candidates alike.
.fb_sel_require_fold_local_formula <- function(formulas, data, folds) {
  terms_list <- lapply(
    Filter(Negate(is.null), formulas),
    function(f) stats::terms(f, data = data)
  )
  for (split in folds) {
    .fb_sel_require_fold_local(terms_list, data, split, route = "bayes")
  }
  invisible(NULL)
}


# .fb_sel_bayes_base_spec() --- the response, the protected fixed part, the
# grouping variables and the data, read off a fitted base.
#
# The same four things `fb_select()` reads, read the same way. They are
# re-derived here rather than passed in because the fold refit has to write a
# NEW fixed formula (protected plus a candidate subset) and needs the
# protected side as terms rather than as a fitted design.
.fb_sel_bayes_base_spec <- function(base) {
  if (isTRUE(base$fb$response_is_expression)) {
    stop(
      "the Bayesian selection route sets the assessment rows' response to ",
      "NA in the fold refit, and this base's response is an expression (`",
      base$fb$response,
      "`) rather than a column, so there is no column to set.",
      call. = FALSE
    )
  }
  labels <- vapply(
    base$fb$fixed_terms %||% list(),
    function(t) t$label %||% t$var %||% NA_character_,
    character(1L)
  )
  labels <- labels[!is.na(labels) & nzchar(labels)]
  protected <- stats::as.formula(paste0(
    "~ ",
    if (length(labels) > 0L) paste(labels, collapse = " + ") else "1"
  ))
  list(
    response = base$fb$response,
    protected = protected,
    protected_labels = labels,
    random_vars = .fb_sel_grouping_columns(base$fb)$columns,
    # The base's offsets travel into every fold refit's formula.
    offsets = base$fb$offsets %||% character(0),
    data = base$data
  )
}


# .fb_sel_provider_installed() --- is a provider package on this machine.
#
# One line, and it exists as a function rather than as an inline
# `requireNamespace()` for one reason: `selection_provider_absent` has to have
# a triggering test, and the only way to test it on a machine where the
# provider IS installed is to be able to substitute this answer. A refusal
# nobody has seen fire is a string.
.fb_sel_provider_installed <- function(package) {
  requireNamespace(package, quietly = TRUE)
}


# .fb_sel_bayes_provider() --- which engine fitted the base, and is it here.
#
# The class of the fitted object is the oracle, not a recorded string:
# `call_info$backend` holds the routing POLICY the call asked for ("auto"),
# which is not the engine the request landed on.
.fb_sel_bayes_provider <- function(base) {
  provider <- if (inherits(base, "flexybayes_brms")) {
    "brms"
  } else if (inherits(base, "flexybayes_inla")) {
    "inla"
  } else {
    stop(
      "`base` must be a fitted flexybayes_brms or flexybayes_inla model; ",
      "this one is <",
      paste(class(base), collapse = ", "),
      ">.",
      call. = FALSE
    )
  }
  package <- if (identical(provider, "brms")) "brms" else "INLA"
  if (!.fb_sel_provider_installed(package)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_provider_absent",
      message = paste0(
        "fb_select(): the Bayesian workflow refits the base on every fold ",
        "with its own engine, and this base was fitted by ",
        provider,
        ", whose package `",
        package,
        "` is not installed. No substitute provider is chosen here: a ",
        "different engine is a different answer under the same name. ",
        "Install it with ",
        if (identical(package, "INLA")) {
          paste0(
            "install.packages(\"INLA\", repos = c(getOption(\"repos\"), ",
            "INLA = \"https://inla.r-inla-download.org/R/stable\"))"
          )
        } else {
          "install.packages(\"brms\")"
        },
        "."
      )
    ))
  }
  provider
}


# --- the fold's engine inputs, and the claim the gate holds them to ----- #

# .fb_sel_bayes_heldout_labels() --- the assessment rows' grouping levels,
# relabelled so no training row carries them.
#
# One fresh level per ORIGINAL level, not one per row: the assessment rows of
# a unit stay in the same new cluster, so their predictive draws keep the
# within-unit correlation the model says they have while losing the link to
# the training rows. Relabelling is a no-op in the common case where the fold
# already holds the whole unit out, and it is the only thing that makes
# `new_cluster` mean what it says when a fold splits one.
.fb_sel_bayes_heldout_labels <- function(data, assess_rows, random_vars) {
  for (v in random_vars) {
    column <- as.character(data[[v]])
    column[assess_rows] <- .fb_sel_heldout_level(column[assess_rows])
    data[[v]] <- factor(column)
  }
  data
}


# .fb_sel_heldout_level() --- the held-out label of an original level.
.fb_sel_heldout_level <- function(level) {
  paste0(level, "::fb_heldout")
}


# .fb_sel_check_heldout_labels() --- every relabelled grouping column carries,
# row by row, the level that row is meant to carry.
#
# The map is one-to-one and compared element-wise: a training row keeps its
# own level, an assessment row carries the held-out label of its own level,
# and no held-out label is a level a training row carries. A check on the
# number of levels, or on rows per level, passes a relabelling that gives one
# site's held-out rows another site's label; and a training level already
# spelled like a held-out label would score the held-out rows as a known
# cluster. Either refuses here, before any engine is called.
.fb_sel_check_heldout_labels <- function(
  data,
  engine_data,
  split,
  random_vars
) {
  assess_rows <- c(
    split$assess_rows, split$unobserved_rows, split$not_estimable_rows
  )
  for (v in random_vars) {
    original <- as.character(data[[v]])
    intended <- original
    intended[assess_rows] <- .fb_sel_heldout_level(original[assess_rows])
    relabelled <- as.character(engine_data[[v]])
    wrong <- which(
      is.na(relabelled) != is.na(intended) |
        (!is.na(relabelled) & !is.na(intended) & relabelled != intended)
    )
    train_rows <- setdiff(seq_along(original), assess_rows)
    shared <- intersect(relabelled[assess_rows], relabelled[train_rows])
    shared <- shared[!is.na(shared)]
    if (length(wrong) == 0L && length(shared) == 0L) {
      next
    }
    detail <- if (length(wrong) > 0L) {
      one <- length(wrong) == 1L
      paste0(
        if (one) "row " else "rows ",
        paste(utils::head(wrong, 5L), collapse = ", "),
        if (one) " carries " else " carry ",
        paste0("'", utils::head(relabelled[wrong], 5L), "'", collapse = ", "),
        " where the map gives ",
        paste0("'", utils::head(intended[wrong], 5L), "'", collapse = ", ")
      )
    } else {
      one <- length(shared) == 1L
      paste0(
        if (one) "the held-out label " else "the held-out labels ",
        paste0("'", utils::head(shared, 5L), "'", collapse = ", "),
        if (one) " is already a level" else " are already levels",
        " of training rows, so the held-out rows would be scored as a known ",
        "cluster. Rename the training ",
        if (one) "level" else "levels",
        ", or use prediction = \"known_cluster\""
      )
    }
    stop(.fb_refusal_condition(
      reason_code = "selection_new_cluster_not_representable",
      message = paste0(
        "fb_select(): prediction = \"new_cluster\" relabels `",
        v,
        "` on the assessment rows of fold ",
        split$label,
        ", and the relabelled column is not the one-to-one map of its ",
        "levels: ",
        detail,
        "."
      ),
      terms = v
    ))
  }
  invisible(NULL)
}


# .fb_sel_bayes_fold_inputs() --- everything one fold hands the engine, and
# the quantities it claims to have taken from its training rows alone.
#
# `learned` is the leakage gate's surface and every entry is chosen so the
# gate can see it: each moves when the TRAINING rows move and none moves when
# only the assessment rows move. `engine_response` is the strongest of them
# --- it is the response vector actually handed to the engine, carrying `NA`
# at every assessment row, so a builder that slipped an assessment response
# into training turns the gate red at the exact object the engine reads.
#
# Nothing is declared that the fit does not depend on. This route hands the
# engine a data frame, not a rotated design, so there is no fold-local
# standardisation to declare here --- the candidate block's centring and
# scaling belong to the penalised route, where a prior scale is written on
# the standardised column. Declaring them anyway would put a quantity in
# front of the gate that no answer in this file reads.
.fb_sel_bayes_fold_inputs <- function(
  data,
  split,
  response,
  protected,
  candidates = NULL,
  random = character(0),
  prediction = c("new_cluster", "known_cluster")
) {
  prediction <- match.arg(prediction)
  y <- as.numeric(data[[response]])

  # A grouping column that also carries a fixed effect cannot be relabelled:
  # the relabel would change the fixed design rather than only the grouping.
  # An argument-shape error rather than a refusal --- there is no user-facing
  # condition to report, only a model this route cannot represent.
  model_vars <- unique(c(
    all.vars(protected),
    if (is.null(candidates)) character(0) else all.vars(candidates)
  ))
  clash <- intersect(random, model_vars)
  if (identical(prediction, "new_cluster") && length(clash) > 0L) {
    stop(
      "prediction = \"new_cluster\" relabels the grouping columns on the ",
      "assessment rows, and ",
      paste(clash, collapse = ", "),
      " carries a fixed effect as well as a grouping, so relabelling would ",
      "change the fixed design. Use prediction = \"known_cluster\".",
      call. = FALSE
    )
  }

  protected_design <- .fb_sel_design(protected, data)
  candidate_design <- if (is.null(candidates)) {
    list(X = matrix(numeric(0), nrow = nrow(data), ncol = 0L))
  } else {
    .fb_sel_design(candidates, data)
  }
  .fb_sel_check_fold_levels(protected_design, split)
  if (!is.null(candidates)) {
    .fb_sel_check_fold_levels(candidate_design, split)
  }

  # The candidate block's intercept belongs to the protected block; two of
  # them would be an aliased pair in every fold.
  x_candidate <- candidate_design$X
  x_candidate <- x_candidate[,
    colnames(x_candidate) != "(Intercept)",
    drop = FALSE
  ]

  engine_data <- if (identical(prediction, "new_cluster")) {
    relabelled <- .fb_sel_bayes_heldout_labels(
      data,
      c(split$assess_rows, split$unobserved_rows, split$not_estimable_rows),
      random
    )
    .fb_sel_check_heldout_labels(data, relabelled, split, random)
    relabelled
  } else {
    data
  }
  engine_response <- y
  engine_response[split$assess_rows] <- NA_real_
  engine_data[[response]] <- engine_response

  x_full <- cbind(protected_design$X, x_candidate)
  x_train <- x_full[split$train_rows, , drop = FALSE]
  y_train <- y[split$train_rows]

  learned <- list(
    engine_response = engine_response,
    train_response = y_train,
    train_xtx = crossprod(x_train),
    train_xty = crossprod(x_train, y_train)
  )

  structure(
    list(
      fold_id = split$fold_id,
      repeat_id = split$repeat_id,
      label = split$label,
      train_rows = split$train_rows,
      assess_rows = split$assess_rows,
      prediction = prediction,
      covariance_estimand = "integrated",
      response = y,
      engine_data = engine_data,
      train_data = engine_data[split$train_rows, , drop = FALSE],
      assess_data = engine_data[split$assess_rows, , drop = FALSE],
      x_protected = protected_design$X,
      x_candidate = x_candidate,
      candidate_names = colnames(x_candidate),
      learned = learned,
      assessment = list(y = y[split$assess_rows])
    ),
    class = "fb_sel_bayes_inputs"
  )
}


# --- the per-fold refit ------------------------------------------------- #

# .fb_sel_bayes_fixed() --- the fixed formula for one candidate subset.
.fb_sel_bayes_fixed <- function(response,
                                protected_labels,
                                candidate_labels,
                                offsets = character(0)) {
  right <- c(protected_labels, candidate_labels, offsets)
  stats::as.formula(paste0(
    response,
    " ~ ",
    if (length(right) > 0L) paste(right, collapse = " + ") else "1"
  ))
}


# .fb_sel_bayes_bar_terms() --- the random terms a brms-grammar base wrote
# into its model formula, as `(lhs | group)` terms, or none.
#
# Empty when the base recorded a `random` formula (the ASReml grammar), which
# the refit re-issues as it is.
.fb_sel_bayes_bar_terms <- function(info) {
  if (!is.null(info$random)) {
    return(character(0))
  }
  model <- info$fixed
  if (inherits(model, "brmsformula")) {
    model <- model$formula
  }
  if (!inherits(model, "formula")) {
    return(character(0))
  }
  labels <- attr(stats::terms(model), "term.labels")
  bars <- labels[grepl("|", labels, fixed = TRUE)]
  if (length(bars) == 0L) {
    return(character(0))
  }
  paste0("(", bars, ")")
}


# .fb_sel_bayes_refit() --- the base model, refitted on a fold, with the same
# engine, structure and priors.
#
# The base's own `call_info` is the record of how it was fitted, and every
# field of it that describes the MODEL is re-issued unchanged; only the fixed
# formula, the data and the seed differ. `prior` is re-issued only when the
# base carried a user prior: the package's own default prior object is stamped
# with `fb_prior_default_basis`, and re-passing it as an explicit prior would
# take a different path through the emitter than the default it came from.
#
# `readout` names what the caller reads off an INLA refit ("varcomp",
# "fitted" or "draws"), and the refit computes only that
# (R/selection_fold_refit.R); NULL fits in full.
.fb_sel_bayes_refit <- function(
  base,
  data,
  fixed,
  provider,
  engine = NULL,
  seed = NULL,
  na_action = NULL,
  carry_known_variance = FALSE,
  fixed_alias = NULL,
  readout = NULL
) {
  info <- base$extras$call_info %||% list()
  # Each row's weight travels in its column of the fold's data
  # (`.fb_sel_with_weights_column()`), so it follows every subset of rows.
  weights <- NULL
  if (!is.null(.fb_sel_base_weights(base))) {
    weights <- data[[.FB_SEL_WEIGHTS_COL]]
    if (is.null(weights)) {
      .fb_sel_refuse_weights_not_carried(
        "a per-fold Bayesian refit handed data without the weights column"
      )
    }
  }
  data <- data[setdiff(names(data), .FB_SEL_WEIGHTS_COL)]
  # A known variance per row (R/known_variance.R) travels in its column of
  # the fold's data. The fold covariance carries it into the refit; every
  # other refit here feeds a posterior predictive score that leaves it out,
  # and refuses.
  known_variance <- NULL
  if (!is.null(.fb_known_variance_term(base$fb %||% base$extras$fb_terms))) {
    if (!isTRUE(carry_known_variance)) {
      .fb_sel_known_variance_refuse_route(
        "a Bayesian refit whose posterior predictive draws are scored"
      )
    }
    known_variance <- data[[.FB_KNOWN_VARIANCE_COL]]
    if (is.null(known_variance)) {
      .fb_known_variance_refuse(paste0(
        "fb_select(): the base carries a known variance per row, and the ",
        "data handed to the fold refit has no `", .FB_KNOWN_VARIANCE_COL,
        "` column to read it from."
      ))
    }
    data <- data[setdiff(names(data), .fb_known_variance_columns())]
  }
  priors <- base$fb$priors
  if (!is.null(attr(priors, "fb_prior_default_basis"))) {
    priors <- NULL
  }
  # A base written in the brms grammar, `y ~ x + (x || g)`, records its
  # random terms in the model formula and no `random`. The fold's formula is
  # rebuilt from the fixed labels, so the random terms are carried across
  # here; without them every fold refit of such a base fitted no random
  # effect.
  bars <- .fb_sel_bayes_bar_terms(info)
  if (length(bars) > 0L) {
    fixed <- .fb_sel_extend_fixed(fixed, bars)
  }

  opts <- .fb_engine_opts(engine) %||% list()
  n_samples <- opts$n_samples %||% info$n_samples
  warmup <- opts$warmup %||% info$warmup
  chains <- opts$chains %||% info$chains

  args <- list(
    fixed = fixed,
    random = info$random,
    residual = info$residual,
    data = data,
    family = info$family %||% base$fb$family,
    link = info$link,
    known_matrices = info$known_matrices %||% list(),
    prior = priors,
    control = info$control,
    # `seed` is a Stan sampler setting. The INLA backend says so, once per
    # call, and a fold loop would say it once per fold; the nested Laplace
    # approximation draws no random numbers, so there is nothing for it to
    # pin. The INLA arm's randomness is in the posterior sampling, and that
    # seed is passed where it is used.
    seed = if (identical(provider, "brms")) seed else NULL,
    backend = provider,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )
  # The scalar prior arguments describe the model as `prior` does: a base
  # fitted with `prior_fixed_sd` or `prior_vc_sd` is refitted under them, and
  # a base fitted without them keeps the package default.
  for (which in c("fixed_sd", "vc_sd")) {
    if (.fb_prior_scalar_supplied(base$fb, which)) {
      args[[paste0("prior_", which)]] <- .fb_prior_scalar_value(base$fb, which)
    }
  }
  if (!is.null(na_action)) {
    args$na_action <- na_action
  } else if (!is.null(info$na_action)) {
    args$na_action <- info$na_action
  }
  if (!is.null(known_variance)) {
    args$known_variance <- known_variance
  }
  if (!is.null(weights)) {
    args$weights <- weights
  }
  fixed_alias <- fixed_alias %||% info$fixed_alias
  if (!is.null(fixed_alias)) {
    args$fixed_alias <- fixed_alias
  }
  # The sampler settings are absent on the INLA path, where the nested
  # Laplace approximation has no chains to run; a NULL is dropped rather
  # than passed, so the call is the one the base itself made.
  if (!is.null(n_samples)) {
    args$n_samples <- n_samples
  }
  if (!is.null(warmup)) {
    args$warmup <- warmup
  }
  if (!is.null(chains)) {
    args$chains <- chains
  }

  plan <- if (identical(provider, "inla") && !is.null(readout)) {
    .fb_fold_refit_plan(readout)
  } else {
    NULL
  }
  fit <- .fb_with_fold_refit(plan, do.call(flexybayes, args))
  # A fold refit's manifest names it as one, beside the route the refit took.
  if (!is.null(fit[["extras"]][["manifest"]])) {
    route <- fit[["extras"]][["manifest"]][["route"]][["route"]]
    fit[["extras"]][["manifest"]][["route"]] <- list(
      route = "selection_fold",
      within = route
    )
  }
  fit
}


# --- posterior predictive draws on the assessment rows ------------------ #

# .fb_sel_posterior_predict_fold() --- the predictive location and scale
# draws for a fold's assessment rows, dispatched on the fitted class.
#
# It returns the per-draw predictive MEAN and the per-draw residual scale
# rather than replicated responses. Both the predictive density and the
# predictive mean follow from those two, and they follow COHERENTLY: a second
# call to a provider's own response sampler would draw its own new-level
# effects and its own noise, so a mean taken from one call and a density from
# another would not describe the same posterior draw. The provider's response
# sampler is used as an independent check on this construction in the test
# file, which is the right place for it.
#
# Dispatch is a plain `inherits()` branch rather than an S3 generic. There
# are exactly two engines and both are in this namespace, so a generic buys
# nothing an extension point would use, and an internal UseMethod() needs a
# runtime registerS3method() in .onLoad() before it resolves from a caller's
# frame -- machinery for a two-way branch. An unknown class refuses here by
# name rather than falling through to a default method.
#
# @return list(mu = <n_draws x n_assess>, sigma = <n_draws>, provider,
#   n_draws, prediction).
.fb_sel_posterior_predict_fold <- function(
  fit,
  assess_rows,
  data,
  prediction = c("new_cluster", "known_cluster"),
  n_draws = 4000L,
  seed = NULL
) {
  if (inherits(fit, "flexybayes_brms")) {
    return(.fb_sel_predict_fold_brms(
      fit,
      assess_rows,
      data,
      prediction = prediction,
      n_draws = n_draws,
      seed = seed
    ))
  }
  if (inherits(fit, "flexybayes_inla")) {
    return(.fb_sel_predict_fold_inla(
      fit,
      assess_rows,
      data,
      prediction = prediction,
      n_draws = n_draws,
      seed = seed
    ))
  }
  stop(
    "no posterior predictive route for a fit of class <",
    paste(class(fit), collapse = ", "),
    ">.",
    call. = FALSE
  )
}


# The brms method. The fit was made on the TRAINING rows only, so the
# assessment rows arrive as `newdata`. `allow_new_levels = TRUE` with
# `sample_new_levels = "gaussian"` is what draws a new cluster's effect from
# its fitted distribution (contract section 6.3, probe C7); under
# `known_cluster` the labels are the fitted ones and the argument is inert.
#
# `draw_ids` selects the draws deterministically and by position, so the
# `sigma` column read from the same draws matrix is aligned draw for draw.
# A random subsample would not be.
.fb_sel_predict_fold_brms <- function(
  fit,
  assess_rows,
  data,
  prediction = c("new_cluster", "known_cluster"),
  n_draws = 4000L,
  seed = NULL
) {
  prediction <- match.arg(prediction)
  newdata <- data[assess_rows, , drop = FALSE]
  total <- brms::ndraws(fit$brms)
  ids <- if (total <= n_draws) {
    seq_len(total)
  } else {
    unique(as.integer(round(seq(1, total, length.out = n_draws))))
  }
  draws <- as.matrix(fit$brms)
  if (!"sigma" %in% colnames(draws)) {
    stop(
      "the brms fit carries no `sigma` draw, so no predictive scale can be ",
      "read from it.",
      call. = FALSE
    )
  }
  .fb_sel_with_seed(seed %||% 1L, {
    mu <- brms::posterior_epred(
      fit$brms,
      newdata = newdata,
      allow_new_levels = TRUE,
      sample_new_levels = "gaussian",
      draw_ids = ids
    )
  })
  list(
    mu = mu,
    sigma = as.numeric(draws[ids, "sigma"]),
    provider = "brms",
    n_draws = length(ids),
    prediction = prediction
  )
}


# The INLA method. INLA has no `newdata` surface --- the package's own
# `predict.flexybayes_inla()` is fixed-effects only --- so the assessment rows
# were carried INTO the fit with the response set to NA, and their linear
# predictor is a latent quantity the fit already sampled. `fit` here is
# therefore the fit over the WHOLE fold, and `assess_rows` indexes its data.
#
# The predictive mean is the "Predictor:" block at the assessment rows; the
# scale is the Gaussian observation precision from the same draw, so the two
# describe one posterior configuration. Where `new_cluster` was asked for, the
# assessment rows already carry a relabelled grouping level that no training
# row shares, so the random effect in that Predictor entry is a draw from the
# random-effect distribution rather than a conditional mean --- the marginal
# predictive the contract's section 7 asks for.
.fb_sel_predict_fold_inla <- function(
  fit,
  assess_rows,
  data,
  prediction = c("new_cluster", "known_cluster"),
  n_draws = 4000L,
  seed = NULL
) {
  prediction <- match.arg(prediction)
  if (nrow(fit$data) != nrow(data)) {
    stop(
      "the INLA refit returned ",
      nrow(fit$data),
      " rows against the ",
      nrow(data),
      " it was given, so an assessment row index no longer names the row it ",
      "was built from.",
      call. = FALSE
    )
  }
  sample_seed <- as.integer(abs(seed %||% 1L) %% .Machine$integer.max)
  if (sample_seed == 0L) {
    sample_seed <- 1L
  }
  draws <- .fb_inla_posterior_sample(
    fit,
    as.integer(n_draws),
    seed = sample_seed,
    num.threads = "1:1"
  )

  names_latent <- rownames(draws[[1L]]$latent)
  predictor <- grep("^Predictor:", names_latent)
  if (length(predictor) != nrow(data)) {
    stop(
      "the sampled latent field carries ",
      length(predictor),
      " Predictor entries for ",
      nrow(data),
      " data rows, so the assessment rows cannot be indexed into it.",
      call. = FALSE
    )
  }
  wanted <- predictor[assess_rows]
  mu <- t(vapply(
    draws,
    function(z) as.numeric(z$latent)[wanted],
    numeric(length(wanted))
  ))

  precision_name <- grep(
    "^Precision for the Gaussian observations",
    names(draws[[1L]]$hyperpar),
    value = TRUE
  )
  if (length(precision_name) != 1L) {
    stop(
      "the sampled hyperparameters carry ",
      length(precision_name),
      " Gaussian observation precisions; this route reads exactly one.",
      call. = FALSE
    )
  }
  sigma <- vapply(
    draws,
    function(z) 1 / sqrt(z$hyperpar[[precision_name]]),
    numeric(1L)
  )

  list(
    mu = mu,
    sigma = sigma,
    provider = "inla",
    n_draws = length(draws),
    prediction = prediction
  )
}


# --- scoring the draws --------------------------------------------------- #

# .fb_sel_logmeanexp() --- log(mean(exp(x))), without leaving log space.
.fb_sel_logmeanexp <- function(x) {
  top <- max(x)
  top + log(mean(exp(x - top)))
}


# .fb_sel_bayes_score_fold() --- the per-row and per-unit records of one fold,
# from predictive draws.
#
# The counterpart of `.fb_sel_score_fold()`, and deliberately the same table
# shape so the two routes aggregate through the same function. Two log
# densities are carried because they answer different questions and the
# contract wants both:
#
#   `elpd`  the sum over the unit's rows of the POINTWISE log predictive
#           density, each row's density obtained by log-mean-exp over draws.
#           This is section 7's definition and the quantity `loo` reports.
#   `nlpd`  the negative JOINT log predictive density of the unit's rows,
#           divided by the number of rows. This is the direct counterpart of
#           the benchmark route's `nlpd`, which uses the joint block of
#           V_S_given_T rather than the marginal variance of each row, and it
#           is the one that sees a wrongly-conditioned cluster effect.
#
# They differ by exactly the within-unit dependence the draws carry, so the
# pair is also a diagnostic: on a unit whose rows are conditionally
# independent given the draw, `-elpd / n` and `nlpd` agree.
.fb_sel_bayes_score_fold <- function(
  y_assess,
  mu,
  sigma,
  cluster,
  label = "1",
  repeat_id = 1L,
  model = "full"
) {
  y_assess <- as.numeric(y_assess)
  n_assess <- length(y_assess)
  n_draws <- length(sigma)
  if (!identical(dim(mu), c(n_draws, n_assess))) {
    stop(
      "the predictive draws are ",
      paste(dim(mu), collapse = " x "),
      " against ",
      n_draws,
      " draws and ",
      n_assess,
      " assessment rows.",
      call. = FALSE
    )
  }

  # log N(y_i | mu_i^(d), sigma^(d)), one row per draw, one column per
  # assessment row.
  logdens <- stats::dnorm(
    matrix(y_assess, nrow = n_draws, ncol = n_assess, byrow = TRUE),
    mean = mu,
    sd = sigma,
    log = TRUE
  )
  lpd_row <- apply(logdens, 2L, .fb_sel_logmeanexp)
  mu_bar <- colMeans(mu)
  residual <- y_assess - mu_bar

  rows <- data.table::data.table(
    model = model,
    fold = label,
    repeat_id = repeat_id,
    cluster = as.character(cluster),
    residual = residual,
    squared_error = residual^2,
    lpd_row = lpd_row,
    nlpd_row = -lpd_row
  )

  groups <- split(
    seq_len(n_assess),
    factor(as.character(cluster), levels = unique(as.character(cluster)))
  )
  clusters <- data.table::rbindlist(lapply(names(groups), function(key) {
    index <- groups[[key]]
    joint <- .fb_sel_logmeanexp(rowSums(logdens[, index, drop = FALSE]))
    data.table::data.table(
      model = model,
      fold = label,
      repeat_id = repeat_id,
      cluster = key,
      n = length(index),
      mse = mean(residual[index]^2),
      rmse = sqrt(mean(residual[index]^2)),
      elpd = sum(lpd_row[index]),
      nlpd = -joint / length(index)
    )
  }))

  list(rows = rows, clusters = clusters)
}


# .fb_sel_aggregate_bayes_scores() --- the reported score and its standard
# error, at the declared independent unit.
#
# `rmse`, `mse` and `nlpd` go straight through `.fb_sel_aggregate_scores()`:
# its input shape is one row per (model, fold, cluster) with the metric as a
# column, which is exactly what `.fb_sel_bayes_score_fold()` returns, so the
# two routes are aggregated by one function and cannot drift apart. `elpd` is
# mirrored rather than reused because that function's `metric` argument is a
# closed vocabulary of the three loss-oriented metrics; the mirror is the same
# average-within-unit-first rule and the same refusal, with the sign flipped
# so the reported `score` is a loss on every metric.
.fb_sel_aggregate_bayes_scores <- function(
  clusters,
  metric = c("elpd", "rmse", "mse", "nlpd"),
  unit = c("cluster", "resample")
) {
  metric <- match.arg(metric)
  unit <- match.arg(unit)
  if (!identical(metric, "elpd")) {
    summary <- .fb_sel_aggregate_scores(clusters, metric = metric, unit = unit)
    summary$orientation <- "loss"
    return(summary)
  }

  key <- if (identical(unit, "cluster")) {
    as.character(clusters$cluster)
  } else {
    as.character(clusters$repeat_id)
  }
  per_unit <- vapply(
    split(clusters$elpd, factor(key, levels = unique(key))),
    mean,
    numeric(1L)
  )
  if (length(per_unit) < 2L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_scoring_unit_undeclared",
      message = paste0(
        "fb_select(): the scoring unit `",
        unit,
        "` yields ",
        length(per_unit),
        " independent unit(s). A standard error needs at least two, and a ",
        "one-standard-error rule without one is not defined."
      )
    ))
  }
  list(
    metric = metric,
    unit = unit,
    # Negated so that `score` is a loss whichever metric was asked for, which
    # is what lets `.fb_sel_select_model()`'s which.min and the increment's
    # "positive is better" sign hold without a per-metric special case. The
    # unnegated per-unit elpd is kept alongside.
    score = -mean(per_unit),
    se = stats::sd(per_unit) / sqrt(length(per_unit)),
    n_units = length(per_unit),
    per_unit = -per_unit,
    elpd_per_unit = per_unit,
    orientation = "loss"
  )
}


# --- the runtime budget -------------------------------------------------- #

# .fb_sel_cost_estimate() --- what the Bayesian route will cost, from one
# timed fit rather than from a guess.
#
# The probe refits the base's own protected model on ALL of its rows, which is
# a little more work than a training fold, so the estimate errs high. It runs
# once; the loop runs `fits_per_fold` times per fold. `probe_seconds` lets a
# caller supply a measurement it already has, in which case nothing is fitted
# at all --- that is the path the budget refusal has to be provable on, since
# a refusal that costs a model fit before it fires is not a refusal that came
# first.
.fb_sel_cost_estimate <- function(
  base,
  resamples,
  n_folds = NULL,
  fits_per_fold = 2L,
  probe_fn = NULL,
  probe_seconds = NULL,
  data = NULL
) {
  if (is.null(n_folds)) {
    n_folds <- if (identical(resamples$scheme, "explicit")) {
      length(unique(resamples$fold_id))
    } else {
      as.integer(resamples$v) * as.integer(resamples$repeats)
    }
  }
  source <- "supplied"
  if (is.null(probe_seconds)) {
    source <- "timed_probe"
    base <- .fb_sel_bayes_normalise_base(base, data)
    spec <- .fb_sel_bayes_base_spec(base)
    provider <- .fb_sel_bayes_provider(base)
    if (is.null(probe_fn)) {
      probe_fn <- function() {
        .fb_sel_bayes_refit(
          base,
          spec$data,
          .fb_sel_bayes_fixed(
            spec$response,
            spec$protected_labels,
            character(0),
            spec$offsets
          ),
          provider,
          seed = 1L
        )
      }
    }
    t0 <- proc.time()
    probe_fn()
    probe_seconds <- unname((proc.time() - t0)["elapsed"])
  }
  n_fits <- as.integer(n_folds) * as.integer(fits_per_fold)
  list(
    per_fit_seconds = probe_seconds,
    n_folds = as.integer(n_folds),
    fits_per_fold = as.integer(fits_per_fold),
    n_fits = n_fits,
    seconds = probe_seconds * n_fits,
    source = source
  )
}


# .fb_sel_bayes_check_budget() --- refuse before the loop, or return silently.
.fb_sel_bayes_check_budget <- function(cost, budget) {
  if (is.null(budget)) {
    return(invisible(cost))
  }
  if (!is.numeric(budget) || length(budget) != 1L || !is.finite(budget)) {
    stop("`budget` must be a single finite number of seconds.", call. = FALSE)
  }
  if (cost$seconds > budget) {
    stop(.fb_refusal_condition(
      reason_code = "penalty_runtime_budget_exceeded",
      message = paste0(
        "fb_select(): the Bayesian route needs ",
        cost$n_fits,
        " engine refits (",
        cost$n_folds,
        " folds x ",
        cost$fits_per_fold,
        " models) at ",
        format(round(cost$per_fit_seconds, 2), nsmall = 2),
        "s each, an estimated ",
        format(round(cost$seconds, 1), nsmall = 1),
        "s against a budget of ",
        format(round(budget, 1), nsmall = 1),
        "s. Nothing was fitted. Raise `budget`, cut the fold count, or use ",
        "the benchmark workflow, whose per-fold cost is a closed-form solve."
      )
    ))
  }
  invisible(cost)
}


# --- the driver ----------------------------------------------------------- #

# .fb_sel_bayes_cv() --- folds, per-fold engine refits, posterior predictive
# draws and scoring, in the order section 4.2 of the v3 contract fixes and
# under the estimand section 2 of v5 names.
#
# The return has the shape `.fb_sel_cv()` returns, so the two workflows are
# comparable field for field, plus the five things that distinguish this one:
# `workflow`, `provider`, `covariance_estimand`, `fits_per_fold` and `timing`.
# `aggregated$score` is a LOSS on every metric --- `elpd` is carried negated,
# and the per-unit elpd is kept unnegated on the cluster table --- so an
# increment written `protected_only - full` is positive when the candidates
# help, whichever metric was asked for.
#
# Two optional hooks serve the penalised reference, whose prior is stated on
# standardised candidate columns:
#
#   `fold_transform(analysis_data)`  learns a preprocessing on a fold's
#       analysis rows and returns `list(apply, record, ...)`. It is handed
#       those rows and nothing else, so no constant it learns can have been
#       read from an assessment row. `apply` is then put on both halves of
#       the fold, the state is passed to `refit_fn` as a fourth argument, and
#       `record` is kept on the fold's substrate.
#   `derived_models`  a named list of `function(fits, inputs, fit_seed)`
#       returning `list(mu, sigma)` for the assessment rows from the fold's
#       fitted models --- a decision product, such as a projection, that is
#       read off a fit rather than fitted itself. Each is scored and
#       aggregated like a fitted model.
.fb_sel_bayes_cv <- function(
  base,
  candidates,
  resamples = fb_resamples(),
  prediction = c("new_cluster", "known_cluster"),
  metric = c("elpd", "rmse", "nlpd"),
  engine = NULL,
  seed,
  budget = NULL,
  data = NULL,
  n_draws = 4000L,
  candidate_sets = NULL,
  cluster = NULL,
  refit_fn = NULL,
  probe_seconds = NULL,
  provider = NULL,
  fold_transform = NULL,
  derived_models = NULL
) {
  prediction <- match.arg(prediction)
  metric <- match.arg(metric)
  if (missing(seed)) {
    seed <- NULL
  }

  base <- .fb_sel_bayes_normalise_base(base, data)
  .fb_select_gate(base)
  # The engine that will do the FITTING, which is not always the engine that
  # fitted the base: a reference model carrying a joint shrinkage prior is
  # emitted on brms whatever the base was fitted by, and the fold loop has to
  # take the brms branch for it -- brms fits the training rows and predicts
  # on the assessment rows, while INLA carries the assessment rows into the
  # fit with an NA response. Left to the base's class, an INLA base scored
  # through a brms reference would read a Predictor block that is not there.
  provider <- provider %||% .fb_sel_bayes_provider(base)
  spec <- .fb_sel_bayes_base_spec(base)
  data <- spec$data
  # Before the cost probe, which fits: a term that cannot take a new level
  # refuses without a fit.
  .fb_sel_require_relabelable(base$fb, prediction)

  # The refit rebuilds the base model on every fold, so every variable the
  # base model names has to be in the data the folds are cut from -- not only
  # the candidates.
  .fb_sel_require_columns(
    data,
    unique(c(
      spec$response,
      all.vars(spec$protected),
      spec$random_vars,
      .fb_offset_vars(spec$offsets)
    )),
    "data"
  )

  folds <- .fb_sel_folds(resamples, data, seed = seed)
  # Before the cost probe, which fits: a formula variable whose value on a row
  # depends on other rows refuses by name.
  .fb_sel_require_fold_local_formula(
    list(spec$protected, candidates),
    data,
    folds
  )
  # Before the cost probe, which fits: the terms each fold can estimate, read
  # off the covariates on its training rows. A protected design they do not
  # identify, a level only in assessment rows and a term aliased in part
  # refuse here (R/selection_candidate_space.R). The penalised reference
  # compiles one prior over the declared block, so it refuses any drop.
  space <- .fb_sel_candidate_space(
    spec$protected,
    candidates,
    data,
    folds,
    groupings = unique(c(spec$random_vars, resamples$unit)),
    route = "bayes"
  )
  if (!is.null(fold_transform)) {
    .fb_sel_require_no_candidate_drop(space)
  }

  cluster_vars <- if (is.null(cluster)) resamples$unit else cluster
  .fb_sel_require_columns(data, cluster_vars, "cluster")
  cluster_key <- if (length(cluster_vars) > 0L) {
    .fb_sel_key(data, cluster_vars)
  } else {
    as.character(seq_len(nrow(data)))
  }

  # The candidate labels are fixed once, before the folds, so no fold can
  # quietly change what "the full model" means.
  candidate_labels <- attr(
    stats::terms(candidates, data = data),
    "term.labels"
  )
  sets <- if (is.null(candidate_sets)) {
    list(protected_only = integer(0), full = seq_along(candidate_labels))
  } else {
    candidate_sets
  }
  if (is.null(names(sets)) || any(!nzchar(names(sets)))) {
    stop(
      "`candidate_sets` must be a NAMED list of term indices.",
      call. = FALSE
    )
  }

  # Before the loop, and before anything is fitted when the caller already
  # knows the per-fit cost.
  cost <- .fb_sel_cost_estimate(
    base,
    resamples,
    n_folds = length(folds),
    fits_per_fold = length(sets),
    probe_seconds = probe_seconds
  )
  .fb_sel_bayes_check_budget(cost, budget)

  if (is.null(refit_fn)) {
    # A decision product reads the fold's fits itself, so they are fitted
    # in full when there is one.
    readout <- if (is.null(derived_models)) "draws" else NULL
    refit_fn <- function(fold_data, fixed, fit_seed) {
      .fb_sel_bayes_refit(
        base,
        fold_data,
        fixed,
        provider,
        engine = engine,
        seed = fit_seed,
        na_action = if (identical(provider, "inla")) "augment" else NULL,
        readout = readout
      )
    }
  }

  rows <- list()
  clusters <- list()
  substrates <- list()
  derived_terms <- list()
  fit_seconds <- numeric(0)
  fit_seeds <- integer(0)
  base_seed <- as.integer(seed %||% 0L)

  # One fold's whole share of the work, as a function of its position: the
  # engine inputs, each model's refit, its predictive draws and its scores.
  # It reads the fold's own rows and the constants fixed before the loop, and
  # it writes nothing outside itself, which is what lets a worker run it
  # (R/selection_parallel.R).
  fold_job <- function(f) {
    split <- folds[[f]]
    inputs <- .fb_sel_bayes_fold_inputs(
      data,
      split,
      response = spec$response,
      protected = spec$protected,
      candidates = candidates,
      random = spec$random_vars,
      prediction = prediction
    )

    # The fold's preprocessing learns on its analysis rows alone, and the
    # same constants are applied to both halves of the fold.
    fold_state <- NULL
    if (!is.null(fold_transform)) {
      fold_state <- fold_transform(data[split$train_rows, , drop = FALSE])
      inputs$engine_data <- fold_state$apply(inputs$engine_data)
      inputs$train_data <- fold_state$apply(inputs$train_data)
      inputs$assess_data <- fold_state$apply(inputs$assess_data)
      inputs$preprocessing <- fold_state$record
    }
    # The fold's candidate space: every model on it leaves out the terms its
    # training rows do not identify, and the record says which.
    inputs$candidate_space <- space$folds[[split$label]]
    fold_fits <- list()
    fold_scored <- list()
    fold_terms <- list()
    fold_seconds <- numeric(0)
    fold_seeds <- integer(0)

    for (s in seq_along(sets)) {
      name <- names(sets)[[s]]
      fixed <- .fb_sel_bayes_fixed(
        spec$response,
        spec$protected_labels,
        setdiff(
          candidate_labels[sets[[name]]],
          inputs$candidate_space$dropped_terms
        ),
        spec$offsets
      )
      # Deterministic per fold and per model, so one run reproduces another
      # and two models on one fold are not handed the same stream.
      fit_seed <- base_seed + 1000L * f + s

      # brms fits the training rows; INLA carries the assessment rows into
      # the fit with an NA response, because it has no newdata surface.
      fold_data <- if (identical(provider, "inla")) {
        inputs$engine_data
      } else {
        inputs$train_data
      }
      t_fit <- proc.time()
      fit <- if (is.null(fold_state)) {
        refit_fn(fold_data, fixed, fit_seed)
      } else {
        refit_fn(fold_data, fixed, fit_seed, fold_state)
      }
      fold_seconds <- c(fold_seconds, unname((proc.time() - t_fit)["elapsed"]))
      fold_seeds <- c(fold_seeds, fit_seed)
      fold_fits[[name]] <- fit

      predictive <- if (identical(provider, "inla")) {
        .fb_sel_posterior_predict_fold(
          fit,
          split$assess_rows,
          inputs$engine_data,
          prediction = prediction,
          n_draws = n_draws,
          seed = fit_seed
        )
      } else {
        .fb_sel_posterior_predict_fold(
          fit,
          seq_len(nrow(inputs$assess_data)),
          inputs$assess_data,
          prediction = prediction,
          n_draws = n_draws,
          seed = fit_seed
        )
      }

      scored <- .fb_sel_bayes_score_fold(
        inputs$assessment$y,
        predictive$mu,
        predictive$sigma,
        cluster_key[split$assess_rows],
        label = split$label,
        repeat_id = split$repeat_id,
        model = name
      )
      fold_scored[[length(fold_scored) + 1L]] <- list(
        model = name,
        rows = scored$rows,
        clusters = scored$clusters
      )
    }

    # Decision products read off this fold's fits, scored on the same
    # assessment rows under the same prediction target and the same scorer.
    for (k in seq_along(derived_models)) {
      name <- names(derived_models)[[k]]
      fit_seed <- base_seed + 1000L * f + length(sets) + k
      predictive <- derived_models[[k]](fold_fits, inputs, fit_seed)
      if (!is.null(predictive$terms)) {
        fold_terms <- c(
          fold_terms,
          stats::setNames(list(predictive$terms), name)
        )
      }
      scored <- .fb_sel_bayes_score_fold(
        inputs$assessment$y,
        predictive$mu,
        predictive$sigma,
        cluster_key[split$assess_rows],
        label = split$label,
        repeat_id = split$repeat_id,
        model = name
      )
      fold_scored[[length(fold_scored) + 1L]] <- list(
        model = name,
        rows = scored$rows,
        clusters = scored$clusters
      )
    }

    list(
      label = split$label,
      inputs = inputs,
      scored = fold_scored,
      derived_terms = fold_terms,
      fit_seconds = fold_seconds,
      fit_seeds = fold_seeds
    )
  }

  # The folds do not read each other, so they may run at the same time. The
  # plan is the machine's and not the model's: it reads one option and the
  # cores, divides the engine threads out of them, and leaves every quantity
  # the fold computes alone (R/selection_parallel.R).
  plan <- .fb_sel_worker_plan(length(folds), provider = provider)
  t_start <- proc.time()
  results <- .fb_sel_fold_lapply(
    seq_along(folds),
    fold_job,
    plan,
    # brms caches a compiled Stan program by the hash of that program, and
    # children compiling one shape at the same time write one cache entry.
    # The first fold compiles every shape the loop needs, in the parent.
    warm_first = identical(provider, "brms")
  )
  elapsed <- unname((proc.time() - t_start)["elapsed"])

  # Assembled in FOLD order, never in the order the workers finished, so
  # every table is built as a serial run builds it.
  for (res in results) {
    substrates[[res$label]] <- res$inputs
    for (item in res$scored) {
      rows[[length(rows) + 1L]] <- item$rows
      clusters[[item$model]] <- c(
        clusters[[item$model]],
        list(item$clusters)
      )
    }
    if (length(res$derived_terms) > 0L) {
      derived_terms[[res$label]] <- res$derived_terms
    }
    fit_seconds <- c(fit_seconds, res$fit_seconds)
    fit_seeds <- c(fit_seeds, res$fit_seeds)
  }

  by_model <- lapply(clusters, data.table::rbindlist)
  aggregated <- data.table::rbindlist(lapply(names(by_model), function(name) {
    summary <- .fb_sel_aggregate_bayes_scores(by_model[[name]], metric = metric)
    data.table::data.table(
      model = name,
      score = summary$score,
      se = summary$se,
      n_units = summary$n_units,
      se_unit = summary$unit,
      complexity = if (name %in% names(sets)) {
        length(sets[[name]])
      } else {
        NA_integer_
      }
    )
  }))

  list(
    folds = folds,
    substrates = substrates,
    # The design every fold was cut from, recorded once: the values the
    # score contract identifies.
    design = .fb_sel_scored_design(spec$protected, candidates, data),
    rows = data.table::rbindlist(rows),
    clusters = data.table::rbindlist(by_model),
    by_model = by_model,
    aggregated = aggregated,
    metric = metric,
    prediction = prediction,
    seed = seed,
    workflow = "bayes",
    provider = provider,
    # Read from the fold object, never from an argument: the folds integrate
    # the covariance because the engine did, and nothing upstream can declare
    # otherwise (contract section 9.2).
    covariance_estimand = unique(vapply(
      substrates,
      function(z) z$covariance_estimand,
      character(1L)
    )),
    score_orientation = "loss",
    # Everything this route learns inside a fold -- the engine refit and any
    # fold-local preprocessing -- reads the analysis rows alone, and the
    # per-fold preprocessing records say what was learned there.
    preprocessing = list(
      learned_on = "analysis_rows",
      fold_local_transform = !is.null(fold_transform),
      by_fold = lapply(substrates, function(z) z$preprocessing)
    ),
    n_draws = n_draws,
    # The declared candidate space, and per fold the terms its training rows
    # left unidentified and every model on that fold left out.
    candidate_space = space,
    # Per fold, the terms each decision product was built from where it
    # chose them on that fold.
    derived_terms = derived_terms,
    fits_per_fold = length(sets),
    # The seed each refit ran under, in fold order and one per entry of
    # `timing$fit_seconds`. Two models on one fold, and one model on two
    # folds, are never handed the same stream, and a run reproduces another.
    fit_seeds = fit_seeds,
    # How the folds were run: a property of the machine rather than of the
    # estimate, recorded here and deliberately absent from the certificate,
    # which must compare equal between a serial and a concurrent run.
    parallel = plan,
    timing = list(
      estimate_seconds = cost$seconds,
      estimate_source = cost$source,
      per_fit_estimate = cost$per_fit_seconds,
      elapsed_seconds = elapsed,
      fit_seconds = fit_seconds,
      n_fits = length(fit_seconds)
    )
  )
}
