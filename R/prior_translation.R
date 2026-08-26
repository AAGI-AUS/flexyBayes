# prior_translation -- which (target, distribution) pairs of the prior
# DSL each engine can carry, written down once.
#
# The sibling of R/family_traits.R and the same discipline: a fact about
# an engine lives in exactly one table, and a test re-reads the engine to
# check the table still describes it. family_traits.R answers "what
# hyperparameter does this likelihood declare"; this file answers "can
# this engine represent this prior on this parameter".
#
# The rule exists because it was broken. `priors_to_inla()` translated
# `pc`, `half_normal` and `uniform` on `sigma` / `sd()` / `smooth()` and
# silently returned nothing for every other pair -- `b()` targets,
# `cor()` targets, and five of the ten distributions -- while
# `prior_summary()` printed each of them as applied. A prior-sensitivity
# analysis on INLA could vary a `student_t` scale and compare the
# engine's default with itself (field-sweep FS-21). The brms side had the
# mirror-image defect: its refusals were real but untyped, so no caller
# could tell "outside the translation table" from "R fell over"
# (field-sweep FS-15).
#
# A blanket refusal was the wrong fix: it would have turned six green
# cells red, because the documented full specification legitimately
# carries a `b("x") ~ normal(...)` row that INLA's `control.fixed` can
# take. So the fix is the table -- translate what the engine represents,
# refuse by name what it does not.
#
# `verdict` is one of:
#   "translate"  the engine carries this pair; `route` names the engine
#                argument it arrives through, and an arrival test proves
#                it reaches that argument.
#   "refuse"     the engine has no representation; the fit refuses typed
#                and names the two remedies (re-express in a translatable
#                distribution, or switch backend).

# ---------------------------------------------------------------- #
# The table                                                         #
# ---------------------------------------------------------------- #

# Every (engine, target, distribution) triple the DSL can produce.
# Targets are the `fb_prior()` target types (R/fb_prior.R
# `.parse_prior_target`); distributions are the names of
# `.fb_prior_family_table()`. Both are read from those functions in
# `.fb_prior_translation_table()` below, so a distribution added to the
# DSL without a verdict here is a build-time error rather than a silent
# drop.
#
# INLA (grounded against INLA 25.10.19):
#
#   Variance components (`sigma`, `sd()`, `smooth()`) are parameterised
#   by theta = log(precision). `pc` is INLA's own `pc.prec`;
#   `half_normal` and `uniform` are exact expression priors on the SD
#   scale (R/fb_prior.R `.inla_halfnormal_sd_expr`,
#   `.inla_uniform_sd_expr`); `exponential(rate)` IS the PC prior --
#   `pc.prec(u, a)` is Exp(-log(a) / u) on sigma -- so it translates
#   exactly with `param = c(1 / rate, exp(-1))`. The remaining six
#   distributions have no representation this package emits.
#
#   `b()` arrives through `control.fixed`, which states a Gaussian prior
#   as a mean and a precision and takes per-coefficient named lists
#   (verified live: `control.fixed = list(mean = list(x = 1, default =
#   0), prec = list(x = 0.0625, default = 0.001))`). Only `normal`
#   translates; INLA's fixed-effect prior is Gaussian by construction.
#
#   `cor()` has no counterpart: the INLA emit produces no correlation
#   hyperparameter.
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
#   `cor()` and `smooth()` have no reachable parameter -- the brms
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

# The verdict for one triple. An unknown engine is "translate" so that
# any backend outside the table -- past or future -- is not refused by
# a table that never described it.
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
# remedies -- re-express in a distribution the engine does carry, or
# move to the other backend -- because a refusal that says only "no" is
# the M5 class the field register files separately.
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
      "backend = \"",
      engine,
      "\" cannot carry ",
      length(bad),
      " of the supplied prior specification",
      if (length(bad) == 1L) "" else "s",
      ":\n",
      paste(unlist(bad), collapse = "\n"),
      "\n\n",
      .fb_prior_translation_menu(engine),
      "\n",
      "Two remedies: re-express the prior in a distribution this ",
      "engine carries\n  (listed above), or pass backend = \"",
      if (identical(engine, "inla")) "brms" else "inla",
      "\", which carries a different set.\n",
      "flexyBayes refuses rather than dropping the specification: ",
      "before 0.9.2 an\n  untranslatable row was silently discarded and ",
      "`prior_summary()` still printed it\n  as applied."
    ),
    family_class = "flexybayes_prior_translation_refusal",
    engine = engine
  ))
}

# The per-engine menu, rendered from the table rather than re-listed.
.fb_prior_translation_menu <- function(engine) {
  tbl <- .fb_prior_translation_table()
  tbl <- tbl[tbl$engine == engine & tbl$verdict == "translate", , drop = FALSE]
  if (!nrow(tbl)) {
    return(paste0("backend = \"", engine, "\" carries no prior DSL row.\n"))
  }
  label <- c(
    sigma = "sigma",
    sd = "sd(group = ...)",
    smooth = "smooth(\"var\")",
    b = "b(\"name\")",
    cor = "cor(group = ...)"
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
    "backend = \"",
    engine,
    "\" carries:\n",
    paste(out, collapse = "\n"),
    "\n"
  )
}


# ---------------------------------------------------------------- #
# Target existence (field-sweep FS-20)                              #
# ---------------------------------------------------------------- #

# Refuse a prior that names a term the model does not have.
#
# `available` is the engine's own list of the parameters a prior can
# attach to, so the check is grounded on what the engine will see rather
# than on a name list this package synthesises: on brms it is
# `brms::default_prior()`'s (class, coef, group) triples, on INLA the
# design-matrix column names and the f()-term labels.
#
# Before 0.9.2 the mismatch was delegated to brms's parser, which
# reported a synthesised Stan parameter name (`b_nonexistent_term`) the
# user never wrote, untyped; on INLA it was not reported at all.
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
