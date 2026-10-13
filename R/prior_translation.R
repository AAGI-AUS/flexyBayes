# prior_translation: which (target, distribution) pairs of the prior
# DSL each engine can carry, written down once.
#
# Each fact about an engine lives in exactly one table, and a test
# re-reads the engine to check the table still describes it. The engine
# translates what it can represent and refuses by name, with a typed
# condition, what it cannot. A pair that is neither translated nor
# refused would be printed by `prior_summary()` as applied while the
# engine used its own default.
#
# `verdict` is one of:
#   "translate"  the engine carries this pair. `route` names the engine
#                argument it arrives through, and an arrival test proves
#                it reaches that argument.
#   "refuse"     the engine has no representation. The fit refuses typed
#                and names the two remedies (re-express in a translatable
#                distribution, or switch backend).

# ---------------------------------------------------------------- #
# The table                                                         #
# ---------------------------------------------------------------- #

# Every (engine, target, distribution) triple the DSL can produce.
# Targets are the `fb_prior()` target types (R/fb_prior.R
# `.parse_prior_target`). Distributions are the names of
# `.fb_prior_family_table()`. Both are read from those functions in
# `.fb_prior_translation_table()` below, so a distribution added to the
# DSL without a verdict here is a build-time error, not a silent drop.
#
# INLA (grounded against INLA 25.10.19):
#
#   Variance components (`sigma`, `sd()`, `smooth()`) are parameterised
#   by theta = log(precision). `pc` is INLA's own `pc.prec`.
#   `half_normal` and `uniform` are exact expression priors on the SD
#   scale (R/fb_prior.R `.inla_halfnormal_sd_expr`,
#   `.inla_uniform_sd_expr`). `exponential(rate)` is the PC prior, since
#   `pc.prec(u, a)` is Exp(-log(a) / u) on sigma, so it translates
#   exactly with `param = c(1 / rate, exp(-1))`. The remaining six
#   distributions have no representation this package emits.
#
#   `b()` arrives through `control.fixed`, which states a Gaussian prior
#   as a mean and a precision and takes per-coefficient named lists
#   (verified live: `control.fixed = list(mean = list(x = 1, default =
#   0), prec = list(x = 0.0625, default = 0.001))`). Only `normal`
#   translates, because INLA's fixed-effect prior is Gaussian.
#
#   `cor()` has no counterpart: the INLA emit produces no correlation
#   hyperparameter.
#
#   `loading()` is the scaling of an fa() copy field, which INLA gives a
#   Gaussian prior as a mean and a precision. Only `normal` translates.
#
# brms (grounded against brms 2.23.0):
#
#   Variance components take any Stan density brms will parse, with the
#   parameter's natural lower bound applied: `normal`, `student_t`,
#   `cauchy`, `exponential`, `gamma` and `uniform` were each round-
#   tripped to `lprior +=` lines in `brms::make_stancode()`. `pc` is
#   emitted as its exponential closed form and `half_normal` /
#   `half_cauchy` as the zero-centred normal / cauchy that truncation at
#   zero makes them.
#
#   `b()` keeps `normal` and `student_t`: those are the two shapes the
#   fixed-effect rows have always emitted, and a positive-support prior
#   (`half_normal`, `exponential`, `gamma`, `pc`) does not describe an
#   unbounded coefficient.
#
#   `pooled()` on `sigma(by = )` and `sd(group = , by = )` is a log-normal
#   hierarchy over the per-level SDs (R/pooled_prior.R): prior rows on the
#   log-sigma intercept and group SD, and stanvars for diag() terms.
#
#   `cor()` and `smooth()` have no reachable parameter. The brms
#   models flexyBayes emits carry classes `b`, `Intercept`, `sd` and
#   `sigma` only (correlated random slopes are refused at ingest and
#   smooths route to INLA), so a row for either would be a prior on a
#   parameter that does not exist.
.fb_prior_translation_verdicts <- function() {
  vc <- c("sigma", "sd", "smooth")
  inla_translate <- c("half_normal", "uniform", "pc", "exponential")
  brms_translate <- c(
    "normal",
    "half_normal",
    "half_cauchy",
    "cauchy",
    "student_t",
    "exponential",
    "gamma",
    "pc",
    "uniform"
  )

  rows <- list()
  add <- function(engine, target, family, verdict, route) {
    rows[[length(rows) + 1L]] <<- data.frame(
      engine = engine,
      target = target,
      family = family,
      verdict = verdict,
      route = route,
      stringsAsFactors = FALSE
    )
  }

  for (tg in vc) {
    for (fam in .fb_prior_supported_families()) {
      if (fam %in% inla_translate) {
        add("inla", tg, fam, "translate", "hyper (theta = log precision)")
      } else {
        add("inla", tg, fam, "refuse", NA_character_)
      }
      if (fam %in% brms_translate) {
        add("brms", tg, fam, "translate", "brms prior row (lb = 0)")
      } else if (identical(fam, "pooled") && tg %in% c("sigma", "sd")) {
        add("brms", tg, fam, "translate", "log-normal hierarchy (stanvars)")
      } else {
        add("brms", tg, fam, "refuse", NA_character_)
      }
    }
  }
  for (fam in .fb_prior_supported_families()) {
    add(
      "inla",
      "b",
      fam,
      if (identical(fam, "normal")) "translate" else "refuse",
      if (identical(fam, "normal")) {
        "control.fixed (mean, prec)"
      } else {
        NA_character_
      }
    )
    add(
      "brms",
      "b",
      fam,
      if (fam %in% c("normal", "student_t")) "translate" else "refuse",
      if (fam %in% c("normal", "student_t")) {
        "brms prior row"
      } else {
        NA_character_
      }
    )
    add(
      "inla",
      "loading",
      fam,
      if (identical(fam, "normal")) "translate" else "refuse",
      if (identical(fam, "normal")) {
        "copy hyper beta (mean, precision)"
      } else {
        NA_character_
      }
    )
    add("brms", "loading", fam, "refuse", NA_character_)
    add("inla", "cor", fam, "refuse", NA_character_)
    add("brms", "cor", fam, "refuse", NA_character_)
    add("brms", "smooth", fam, "refuse", NA_character_)
  }
  do.call(rbind, rows)
}

# The table as a data frame, built once per session. `smooth` on brms is
# added by the loop above and overrides the variance-component row, so
# the de-duplication keeps the last verdict written for each triple.
.fb_prior_translation_table <- function() {
  tbl <- .fb_prior_translation_verdicts()
  key <- paste(tbl$engine, tbl$target, tbl$family, sep = "|")
  tbl[!duplicated(key, fromLast = TRUE), , drop = FALSE]
}

# The verdict for one triple. An engine outside the table is
# "translate", so the table never refuses a backend it does not describe.
.fb_prior_translation_lookup <- function(engine, target_type, family) {
  if (!engine %in% c("inla", "brms")) {
    return("translate")
  }
  tbl <- .fb_prior_translation_table()
  hit <- tbl$engine == engine &
    tbl$target == target_type &
    tbl$family == family
  if (!any(hit)) {
    return("refuse")
  }
  tbl$verdict[which(hit)[[1L]]]
}


# ---------------------------------------------------------------- #
# The fit-entry check                                               #
# ---------------------------------------------------------------- #

# Refuse, by name, every spec in `prior` the engine cannot carry.
# Returns invisible TRUE when every spec translates.
#
# The message lists the untranslatable rows one per line and names both
# remedies: re-express in a distribution the engine does carry, or move
# to the other backend.
.fb_check_prior_translation <- function(prior, engine) {
  if (!inherits(prior, "fb_prior")) {
    return(invisible(TRUE))
  }
  bad <- list()
  for (s in prior$specs) {
    tt <- s$target$type %||% "<unknown>"
    fam <- s$spec$family %||% "<unknown>"
    if (identical(.fb_prior_translation_lookup(engine, tt, fam), "refuse")) {
      bad[[length(bad) + 1L]] <- paste0(
        "  ",
        .format_prior_target(s$target),
        " ~ ",
        .format_prior_distribution(s$spec)
      )
    }
  }
  if (!length(bad)) {
    return(invisible(TRUE))
  }
  stop(.fb_refusal_condition(
    reason_code = "prior_not_translatable_for_backend",
    message = paste0(
      "flexyBayes does not translate ",
      length(bad),
      " of the supplied prior specification",
      if (length(bad) == 1L) "" else "s",
      " for backend = \"",
      engine,
      "\":\n",
      paste(unlist(bad), collapse = "\n"),
      "\n\n",
      .fb_prior_translation_menu(engine),
      "\n",
      "Either re-express the prior in a distribution listed above, or ",
      "pass backend = \"",
      if (identical(engine, "inla")) "brms" else "inla",
      "\",\n  for which flexyBayes translates a different set."
    ),
    family_class = "flexybayes_prior_translation_refusal",
    engine = engine
  ))
}

# The per-engine menu, rendered from the table, not re-listed.
.fb_prior_translation_menu <- function(engine) {
  tbl <- .fb_prior_translation_table()
  tbl <- tbl[tbl$engine == engine & tbl$verdict == "translate", , drop = FALSE]
  if (!nrow(tbl)) {
    return(paste0(
      "flexyBayes translates no prior for backend = \"", engine, "\".\n"
    ))
  }
  label <- c(
    sigma = "sigma",
    sd = "sd(group = ...)",
    smooth = "smooth(\"var\")",
    b = "b(\"name\")",
    cor = "cor(group = ...)",
    loading = "loading(group = ...)"
  )
  out <- character(0)
  for (tg in unique(tbl$target)) {
    fams <- tbl$family[tbl$target == tg]
    out <- c(
      out,
      paste0(
        "  ",
        label[[tg]] %||% tg,
        " ~ ",
        paste(fams, collapse = " | ")
      )
    )
  }
  paste0(
    "flexyBayes translates these priors for backend = \"",
    engine,
    "\":\n",
    paste(out, collapse = "\n"),
    "\n"
  )
}


# ---------------------------------------------------------------- #
# Target existence                                                  #
# ---------------------------------------------------------------- #

# Refuse a prior that names a term the model does not have.
#
# `available` is the engine's own list of the parameters a prior can
# attach to, not a name list this package synthesises. On brms it is
# `brms::default_prior()`'s (class, coef, group) triples, on INLA the
# design-matrix column names and the f()-term labels.
.fb_stop_prior_target_absent <- function(
  target_label,
  kind,
  available,
  engine
) {
  available <- unique(available[nzchar(available)])
  stop(.fb_refusal_condition(
    reason_code = "prior_target_not_in_model",
    message = paste0(
      "Prior target ",
      target_label,
      " names a ",
      kind,
      " this model does not have.\n",
      "The model's ",
      kind,
      "s on backend = \"",
      engine,
      "\": ",
      if (length(available)) {
        paste(paste0("\"", available, "\""), collapse = ", ")
      } else {
        "(none)"
      },
      ".\n",
      "Check the spelling, or drop the row from fb_prior()."
    ),
    family_class = "flexybayes_prior_target_absent",
    engine = engine,
    target_label = target_label
  ))
}
