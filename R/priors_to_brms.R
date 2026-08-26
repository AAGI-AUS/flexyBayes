# priors_to_brms -- flexyBayes fb_prior -> brms prior translation.
#
# Stan-passthrough emit on fb_brms(backend = "brms"). The translation
# follows a ten-row table -- the original cross-engine table
# (Intercept normal, b normal, sigma uniform / half_normal / pc,
# sd half_normal / pc / uniform on the intercept-variance keyed by
# group) plus a tenth row for the uncorrelated random-slope
# variance keyed by (class = "sd", coef = "<slope_var>", group =
# "<grouping_factor>"). brms-native prior shape: a data.frame built
# by stacking brms::prior_string() rows, keyed by the brms
# class / coef / group triple. brms's parser then attaches each row
# to the matching parameter inside the generated Stan code.
#
# Two entry points:
#
#   .priors_to_brms_specs(prior, fb, prior_fixed_sd, prior_vc_sd)
#     Pure list output. Each element is a named list with fields
#     `string` (the brms density expression as a character scalar),
#     `class` (one of "b", "Intercept", "sigma", "sd"),
#     `coef` (NA_character_ when not coef-keyed),
#     `group` (NA_character_ when not sd-group-keyed). Suitable for
#     unit testing without brms installed.
#
#   .priors_to_brms(prior, fb, prior_fixed_sd, prior_vc_sd)
#     Calls the pure specs path and rbinds the rows into a single
#     brms prior data.frame via brms::prior_string(). Requires brms.
#     Returns NULL when the spec list is empty (brms then applies its
#     own default flat priors -- caller may want to inject a fall-
#     back; the Stan emit path always feeds the legacy-scalar bridge
#     so the result is never empty in normal use).
#
# Which (target, distribution) pairs the table covers is declared once,
# in `.fb_prior_translation_table()` (R/prior_translation.R), and read
# from there by both this emit and the refusal messages -- a menu
# re-listed in prose is a menu that drifts. At 0.9.2 the variance-
# component targets gained `normal`, `student_t`, `cauchy`,
# `half_cauchy`, `exponential` and `gamma`, each round-tripped through
# `brms::make_stancode()`; `cor()` and `smooth()` are refusals rather
# than deferrals, because no brms model this package emits carries
# either parameter.
#
# Anything outside the table raises a typed refusal (condition class
# `flexybayes_refusal_prior_not_translatable_for_backend`) naming both
# the unsupported spec and the brms-native fallback: pass `prior`
# directly through `...` to brms::brm() via fb_brms()'s pass-through
# `...` argument.

# ---------------------------------------------------------------- #
# Pure spec list                                                    #
# ---------------------------------------------------------------- #

.priors_to_brms_specs <- function(
  prior,
  fb,
  prior_fixed_sd = 100,
  prior_vc_sd = 1
) {
  family_has_sigma <- .brms_family_has_sigma(fb$family)

  if (
    is.null(prior) ||
      (is.list(prior) && isTRUE(prior$legacy))
  ) {
    # Legacy-scalar bridge: normal(0, prior_fixed_sd) on every b
    # coefficient (including Intercept) + lognormal(0, prior_vc_sd)
    # on sigma (only for families that carry one) and every named
    # sd group.
    return(.brms_legacy_specs(
      fb,
      prior_fixed_sd = prior_fixed_sd,
      prior_vc_sd = prior_vc_sd,
      family_has_sigma = family_has_sigma
    ))
  }

  .check_fb_prior(
    prior,
    "`.priors_to_brms()` expects an `fb_prior` object or NULL ",
    "(legacy bridge). Got: ",
    paste(class(prior), collapse = "/"),
    "."
  )

  # Build a lookup of (slope_var, grouping_factor) pairs for
  # every simple_slope_uncor random term in the IR. The flexyBayes
  # default-prior expansion synthesises a group name `<slope_var>_<g>`
  # for the slope-variance hyperparameter so the per-group uniform-on-
  # SD machinery reaches it; on the brms passthrough path we must
  # translate that synthesised name back into brms's
  # (class = "sd", coef = "<slope_var>", group = "<g>") row -- a bare
  # sd(group = "<slope_var>_<g>") row would not match any Stan-side
  # parameter and brms would refuse with "priors do not correspond
  # to any model parameter".
  slope_lookup <- list()
  for (t in fb$random_terms %||% list()) {
    if (
      identical(t$type, "simple_slope_uncor") &&
        !is.null(t$slope_var) &&
        nzchar(t$slope_var) &&
        !is.null(t$var) &&
        nzchar(t$var)
    ) {
      syn <- paste0(t$slope_var, "_", t$var)
      slope_lookup[[syn]] <- list(coef = t$slope_var, group = t$var)
    }
  }

  specs <- list()
  for (s in prior$specs) {
    # Silently drop sigma specs on families that have no residual
    # sigma parameter (Bernoulli, Poisson, etc.). brms's parser
    # would otherwise reject the prior with `The following priors do
    # not correspond to any model parameter`. Mirrors the behaviour
    # already established in priors_to_inla() (emit_inla.R
    # .build_inla_control_family).
    if (
      identical(s$target$type %||% "", "sigma") &&
        !family_has_sigma
    ) {
      next
    }

    # Slope-variance unwrap (per slope_lookup above).
    if (
      identical(s$target$type %||% "", "sd") &&
        !is.null(s$target$group) &&
        s$target$group %in% names(slope_lookup)
    ) {
      pair <- slope_lookup[[s$target$group]]
      s$target$group <- pair$group
      s$target$coef <- pair$coef
    }

    row <- .one_brms_spec(s)
    if (is.null(row)) {
      next
    }
    specs[[length(specs) + 1L]] <- row
  }

  # `prior_fixed_sd` covers the fixed effects; an fb_prior() covers the
  # parameters it names. The two are not alternatives, and before 0.9.2
  # supplying the scalar alongside the auto-default variance prior --
  # which is what passing it alone produces -- silently dropped it: the
  # legacy bridge is the only consumer of the scalar and the fb_prior
  # branch never reached it. The fixed-effect rows are added here when
  # the caller wrote the argument, and any coefficient the fb_prior names
  # itself keeps its own row (brms resolves the coef-keyed row over the
  # class-wide one).
  c(.brms_fixed_scalar_specs(fb, prior_fixed_sd, specs), specs)
}

# Class-wide fixed-effect rows for an explicitly supplied
# `prior_fixed_sd`: one `Intercept` row and one blanket `b` row, each
# omitted when the fb_prior already carries a row brms would see as a
# duplicate of it.
.brms_fixed_scalar_specs <- function(fb, prior_fixed_sd, existing) {
  if (!.fb_prior_scalar_supplied(fb, "fixed_sd")) {
    return(list())
  }
  has_row <- function(cls) {
    any(vapply(
      existing,
      function(r) {
        identical(r$class, cls) && is.na(r$coef %||% NA_character_)
      },
      logical(1)
    ))
  }
  string <- sprintf("normal(0, %s)", .fmt_num(prior_fixed_sd))
  out <- list()
  if (isTRUE(fb$intercept) && !has_row("Intercept")) {
    out[[length(out) + 1L]] <- list(
      string = string,
      class = "Intercept",
      coef = NA_character_,
      group = NA_character_,
      lb = NA_real_,
      ub = NA_real_
    )
  }
  if (length(fb$fixed_terms) > 0L && !has_row("b")) {
    out[[length(out) + 1L]] <- list(
      string = string,
      class = "b",
      coef = NA_character_,
      group = NA_character_,
      lb = NA_real_,
      ub = NA_real_
    )
  }
  out
}

# Families brms parameterises with a residual `sigma` hyperparameter.
# The roster lives in one place (`.fb_brms_families_with_sigma()`,
# R/family_traits.R) because this predicate and the
# heteroscedastic-residual gate in emit_brms.R used to carry separate
# hand-maintained copies that disagreed about gamma and beta. Dropping a
# sigma prior on a family that has no sigma keeps the cross-engine prior
# surface uniform; sending one is a fit brms refuses.
.brms_family_has_sigma <- function(fam) {
  .fb_family_has_brms_sigma(fam)
}

# Translate a single fb_prior spec to a brms-prior row. Returns
# NULL when the spec is recognised but maps to no brms row (no
# such case in v0.2 -- structured refusal fires instead).
.one_brms_spec <- function(s) {
  target <- s$target
  spec <- s$spec

  # ------------------------------ b() / Intercept ------------------------ #
  if (identical(target$type, "b")) {
    nm <- target$name
    if (identical(spec$family, "normal")) {
      mean_v <- as.numeric(.named_or_positional(
        spec$args,
        "mean",
        1L,
        default = 0
      ))
      sd_v <- as.numeric(.named_or_positional(spec$args, "sd", 2L, default = 1))
      string <- sprintf("normal(%s, %s)", .fmt_num(mean_v), .fmt_num(sd_v))
      return(.brms_b_row(nm, string))
    }
    if (identical(spec$family, "student_t")) {
      df_v <- as.numeric(.named_or_positional(spec$args, "df", 1L, default = 3))
      sc_v <- as.numeric(.named_or_positional(
        spec$args,
        "scale",
        2L,
        default = 1
      ))
      string <- sprintf("student_t(%s, 0, %s)", .fmt_num(df_v), .fmt_num(sc_v))
      return(.brms_b_row(nm, string))
    }
    .brms_unsupported(spec$family, "b()")
  }

  # ------------------------------ sigma ---------------------------------- #
  if (identical(target$type, "sigma")) {
    dens <- .brms_vc_density(spec)
    if (is.null(dens)) {
      .brms_unsupported(spec$family, "sigma")
    }
    return(list(
      string = dens$string,
      class = "sigma",
      coef = NA_character_,
      group = NA_character_,
      lb = dens$lb,
      ub = dens$ub
    ))
  }

  # ------------------------------ sd(group) ------------------------------ #
  if (identical(target$type, "sd")) {
    grp <- target$group
    # Slope-variance row: when the target carries a `coef`
    # slot, the sd row is keyed on (class = "sd", coef = "<x>",
    # group = "<g>") rather than the bare intercept-only
    # (class = "sd", group = "<g>") shape.
    coef_v <- target$coef %||% NA_character_
    dens <- .brms_vc_density(spec)
    if (is.null(dens)) {
      .brms_unsupported(
        spec$family,
        paste0(
          "sd(group = \"",
          grp,
          "\"",
          if (!is.na(coef_v)) {
            paste0(", coef = \"", coef_v, "\"")
          } else {
            ""
          },
          ")"
        )
      )
    }
    return(list(
      string = dens$string,
      class = "sd",
      coef = coef_v,
      group = grp,
      lb = dens$lb,
      ub = dens$ub
    ))
  }

  # ------------------------------ cor / smooth / name -------------------- #
  if (target$type %in% c("cor", "smooth", "name")) {
    .brms_unsupported_target(target)
  }

  .brms_unsupported_target(target)
}

# The Stan density a variance-component prior becomes, plus the bounds
# brms needs on the parameter. One builder for `sigma` and `sd()`: the
# two targets differ in how the row is keyed, never in what the density
# is, and before 0.9.2 they carried two copies of the same three
# branches.
#
# NULL means the distribution has no brms row -- the caller raises the
# typed refusal, which names the target as well as the family.
#
# Every string here was round-tripped through `brms::make_stancode()`
# against brms 2.23.0 and appears in the generated program as an
# `lprior +=` line; `tests/testthat/test-prior-translation.R` re-runs
# that round trip. The DSL lives on the standard-deviation scale, so
# every row is bounded below at zero: brms renormalises a two-sided
# density over the truncated support, which is what makes
# `half_normal(scale)` and `normal(0, scale)` the same prior on a
# variance component.
.brms_vc_density <- function(spec) {
  args <- spec$args
  switch(
    spec$family,
    "uniform" = {
      lo <- as.numeric(args$lower %||% 0)
      hi <- as.numeric(args$upper)
      # brms requires lb/ub for uniform priors so the parser bounds the
      # parameter to the prior's support.
      list(
        string = sprintf("uniform(%s, %s)", .fmt_num(lo), .fmt_num(hi)),
        lb = max(0, lo),
        ub = hi
      )
    },
    "half_normal" = list(
      string = sprintf(
        "normal(0, %s)",
        .fmt_num(as.numeric(args$scale %||% 1))
      ),
      lb = 0,
      ub = NA_real_
    ),
    "half_cauchy" = list(
      string = sprintf(
        "cauchy(0, %s)",
        .fmt_num(as.numeric(args$scale %||% 1))
      ),
      lb = 0,
      ub = NA_real_
    ),
    "pc" = {
      # The PC prior on a Gaussian variance component is exponential on
      # the SD scale with rate -log(prob) / upper (Simpson et al. 2017).
      u <- as.numeric(args$upper %||% 1)
      p <- as.numeric(args$prob %||% 0.01)
      list(
        string = sprintf("exponential(%s)", .fmt_num(-log(p) / u)),
        lb = 0,
        ub = NA_real_
      )
    },
    "normal" = list(
      string = sprintf(
        "normal(%s, %s)",
        .fmt_num(as.numeric(.named_or_positional(args, "mean", 1L, 0))),
        .fmt_num(as.numeric(.named_or_positional(args, "sd", 2L, 1)))
      ),
      lb = 0,
      ub = NA_real_
    ),
    "cauchy" = list(
      string = sprintf(
        "cauchy(%s, %s)",
        .fmt_num(as.numeric(args$location %||% 0)),
        .fmt_num(as.numeric(args$scale %||% 1))
      ),
      lb = 0,
      ub = NA_real_
    ),
    "student_t" = list(
      string = sprintf(
        "student_t(%s, %s, %s)",
        .fmt_num(as.numeric(args$df %||% 3)),
        .fmt_num(as.numeric(args$location %||% 0)),
        .fmt_num(as.numeric(args$scale %||% 1))
      ),
      lb = 0,
      ub = NA_real_
    ),
    "exponential" = list(
      string = sprintf("exponential(%s)", .fmt_num(as.numeric(args$rate))),
      lb = 0,
      ub = NA_real_
    ),
    "gamma" = list(
      string = sprintf(
        "gamma(%s, %s)",
        .fmt_num(as.numeric(args$shape)),
        .fmt_num(as.numeric(args$rate))
      ),
      lb = 0,
      ub = NA_real_
    ),
    NULL
  )
}

# brms prior row for a b() / Intercept target. Intercept is class
# "Intercept" in brms (no coef); other named coefficients are class
# "b" with the coef field naming the column.
.brms_b_row <- function(coef_name, string) {
  is_intercept <- coef_name %in% c("(Intercept)", "Intercept")
  if (is_intercept) {
    list(
      string = string,
      class = "Intercept",
      coef = NA_character_,
      group = NA_character_,
      lb = NA_real_,
      ub = NA_real_
    )
  } else {
    list(
      string = string,
      class = "b",
      coef = coef_name,
      group = NA_character_,
      lb = NA_real_,
      ub = NA_real_
    )
  }
}

# Numeric formatter -- writes a brms-parseable scalar without
# scientific notation, leading whitespace, or trailing zeros.
# Integers stay integer-shaped; fractional values use 15-digit
# significand precision (matches Stan's double precision) with the
# scientific suffix disabled.
.fmt_num <- function(x) {
  if (!is.finite(x)) {
    stop("Non-finite numeric in brms prior spec: ", x, call. = FALSE)
  }
  if (x == as.integer(x) && abs(x) < 1e9) {
    return(format(as.integer(x)))
  }
  format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE, digits = 15)
}

# Fetch a value from the args list by name first, falling back to
# the positional slot when the user passed the argument unnamed.
# `default` fires when neither route resolves. Mirrors R's standard
# match-by-name-then-position calling convention so that
# `normal(0, 50)` and `normal(mean = 0, sd = 50)` yield the same
# spec.
.named_or_positional <- function(args, nm, pos, default) {
  if (!is.null(args[[nm]])) {
    return(args[[nm]])
  }
  if (length(args) >= pos) {
    val <- args[[pos]]
    # Skip positional slots that R named with the empty string and
    # whose actual name in the call differs from `nm` (rare but
    # possible if the user wrote `normal(sd = 50)`).
    nms <- names(args) %||% rep("", length(args))
    if (!nzchar(nms[[pos]] %||% "")) return(val)
  }
  default
}

# Refusal for a distribution the brms translation table does not carry
# on a recognised target.
#
# Typed since 0.9.2. The message was already informative -- it named the
# supported table -- but the condition was a bare `simpleError`, so a
# gate or a wrapper could not tell "outside the translation table" from
# "R fell over", which is the distinction the typed-refusal contract
# exists to provide (field-sweep FS-15). The supported table is now
# rendered from `.fb_prior_translation_table()` rather than re-listed
# here, so the message cannot drift from what the emit actually does.
.brms_unsupported <- function(family, target_label) {
  stop(.fb_refusal_condition(
    reason_code = "prior_not_translatable_for_backend",
    message = paste0(
      "The Stan-passthrough emit does not translate `",
      family,
      "` priors on ",
      target_label,
      ".\n",
      .fb_prior_translation_menu("brms"),
      "Two remedies: re-express the prior in a distribution brms ",
      "carries (listed\n  above), or pass a brms `prior` object ",
      "directly via `...` (for example\n  ",
      "`fb_brms(..., backend = \"brms\", prior = brms::prior(...))`), ",
      "which bypasses\n  the flexyBayes prior DSL on the Stan path."
    ),
    family_class = "flexybayes_prior_translation_refusal",
    engine = "brms"
  ))
}

# Refusal for a target the brms models this package emits do not carry.
#
# `cor()` and `smooth()` are not deferrals in the sense the pre-0.9.2
# message implied: a flexyBayes brms model carries classes `b`,
# `Intercept`, `sd` and `sigma` and nothing else -- correlated random
# slopes are refused at ingest and smooths route to INLA -- so a row for
# either would be a prior on a parameter that does not exist.
.brms_unsupported_target <- function(target) {
  ttype <- target$type %||% "<unknown>"
  label <- switch(
    ttype,
    "cor" = paste0("cor(group = \"", target$group, "\")"),
    "smooth" = paste0("smooth(\"", target$var, "\")"),
    "name" = paste0("name = \"", target$name, "\""),
    ttype
  )
  why <- switch(
    ttype,
    "cor" = paste0(
      "No brms model this package emits carries a correlation ",
      "parameter:\n  correlated random slopes are refused at ingest, so ",
      "there is no `cor` class\n  for the prior to attach to."
    ),
    "smooth" = paste0(
      "No brms model this package emits carries a smooth: `spl()` terms ",
      "route to\n  the INLA backend, where `smooth(\"var\")` priors do ",
      "translate."
    ),
    paste0(
      "The Stan-passthrough emit carries prior targets sigma, ",
      "sd(group = ...) and\n  b(\"name\") only."
    )
  )
  stop(.fb_refusal_condition(
    reason_code = "prior_not_translatable_for_backend",
    message = paste0(
      "backend = \"brms\" cannot carry the prior target ",
      label,
      ".\n",
      why,
      "\n",
      .fb_prior_translation_menu("brms"),
      "Two remedies: drop the row, or -- for a smooth -- pass ",
      "backend = \"inla\",\n  which carries it."
    ),
    family_class = "flexybayes_prior_translation_refusal",
    engine = "brms"
  ))
}

# Legacy-scalar bridge: builds the spec list that mirrors the v0.1
# legacy prior (normal(0, prior_fixed_sd) on every fixed-effect
# coefficient incl. Intercept; lognormal(0, prior_vc_sd) on sigma and
# every named sd group). Uses lognormal because the legacy greta-side
# prior on variance components is lognormal(0, vc_sd) (codegen.R
# .sigma_decl); the brms-side row preserves that shape without a
# positivity bound (brms applies the natural sigma >= 0 anyway).
.brms_legacy_specs <- function(
  fb,
  prior_fixed_sd,
  prior_vc_sd,
  family_has_sigma = TRUE
) {
  specs <- list()

  # Intercept (class = "Intercept") gets the same normal prior.
  if (isTRUE(fb$intercept)) {
    specs[[length(specs) + 1L]] <- list(
      string = sprintf("normal(0, %s)", .fmt_num(prior_fixed_sd)),
      class = "Intercept",
      coef = NA_character_,
      group = NA_character_,
      lb = NA_real_,
      ub = NA_real_
    )
  }

  # Class = "b" applies to every non-intercept fixed-effect coef.
  # A single class-only row in brms covers all `b` coefficients.
  if (
    length(fb$fixed_terms) > 0L ||
      length(fb$residual_terms %||% list()) > 0L
  ) {
    specs[[length(specs) + 1L]] <- list(
      string = sprintf("normal(0, %s)", .fmt_num(prior_fixed_sd)),
      class = "b",
      coef = NA_character_,
      group = NA_character_,
      lb = NA_real_,
      ub = NA_real_
    )
  }

  # sigma: lognormal(0, prior_vc_sd) on the natural sigma scale.
  # brms's lognormal is parsed as lognormal(meanlog, sdlog); the
  # legacy convention is lognormal(0, vc_sd) i.e. meanlog = 0.
  # Skip on families that do not parameterise a residual sigma
  # (Bernoulli, Poisson, ...) so brms's prior parser does not raise
  # "priors do not correspond to any model parameter".
  if (isTRUE(family_has_sigma)) {
    specs[[length(specs) + 1L]] <- list(
      string = sprintf("lognormal(0, %s)", .fmt_num(prior_vc_sd)),
      class = "sigma",
      coef = NA_character_,
      group = NA_character_,
      lb = 0,
      ub = NA_real_
    )
  }

  # sd(group) -- one row per named random group.
  # simple_slope_uncor adds a 10th row class to the
  # legacy table -- the slope-variance prior, keyed on
  # (class = "sd", coef = "<slope_var>", group = "<grouping_factor>").
  # brms's set_prior() discriminates intercept- and slope-variance
  # rows on the same grouping factor by the `coef` field; intercept
  # rows leave `coef` unset (NA / empty), slope rows set it to the
  # slope variable name. The double-pipe form (x || g) and the
  # double-pipe + intercept form (1 + x || g) both ship the slope
  # row; only (1 + x || g) also ships the intercept row.
  for (term in fb$random_terms %||% list()) {
    if (is.null(term$type)) {
      next
    }
    grp <- term$var
    if (is.null(grp) || !nzchar(grp)) {
      next
    }

    if (term$type %in% c("simple", "ide", "id")) {
      specs[[length(specs) + 1L]] <- list(
        string = sprintf("lognormal(0, %s)", .fmt_num(prior_vc_sd)),
        class = "sd",
        coef = NA_character_,
        group = grp,
        lb = 0,
        ub = NA_real_
      )
      next
    }

    if (identical(term$type, "simple_slope_uncor")) {
      sv <- term$slope_var
      if (isTRUE(term$with_intercept)) {
        specs[[length(specs) + 1L]] <- list(
          string = sprintf("lognormal(0, %s)", .fmt_num(prior_vc_sd)),
          class = "sd",
          coef = NA_character_,
          group = grp,
          lb = 0,
          ub = NA_real_
        )
      }
      if (!is.null(sv) && nzchar(sv)) {
        specs[[length(specs) + 1L]] <- list(
          string = sprintf("lognormal(0, %s)", .fmt_num(prior_vc_sd)),
          class = "sd",
          coef = sv,
          group = grp,
          lb = 0,
          ub = NA_real_
        )
      }
      next
    }
  }

  specs
}


# ---------------------------------------------------------------- #
# brms object construction                                          #
# ---------------------------------------------------------------- #

# Public-style internal entry: returns a single brms-prior data.frame
# (the kind brms::brm() accepts on its `prior` argument). Requires
# brms. emit_brms.R calls this after the brms namespace check.
.priors_to_brms <- function(prior, fb, prior_fixed_sd = 100, prior_vc_sd = 1) {
  specs <- .priors_to_brms_specs(
    prior,
    fb,
    prior_fixed_sd = prior_fixed_sd,
    prior_vc_sd = prior_vc_sd
  )
  specs <- .brms_retarget_sigma_for_heterogeneous_residual(specs, fb)
  .brms_specs_to_object(specs)
}

# .brms_retarget_sigma_for_heterogeneous_residual() --- move the residual
# prior onto the parameters that actually exist when sigma is modelled.
#
# With a heterogeneous residual the emit becomes `sigma ~ 0 + f`, and the
# scalar `sigma` parameter ceases to exist -- it is replaced by one
# coefficient per level, on the LOG scale. A prior aimed at the scalar then
# matches nothing, and brms refuses outright:
#   "The following priors do not correspond to any model parameter:
#    <lower=0,upper=11.4> sigma ~ uniform(0, 11.4)"
# which is how this was found. The precedent for the fix is already in this
# file: sigma specs are dropped for families that carry no sigma, for exactly
# the same reason.
#
# C6 observation weights hit the identical wall by a different door:
# .fb_to_brms_formula() lowers `weights=` as a known OFFSET on log(sigma)
# (R/emit_brms.R documents why this, and not brms's own `weights()`
# addition term, is the precision-weighting-correct mapping), which is
# ALSO a `sigma ~ ...` distributional sub-formula -- `sigma ~ 1 +
# offset(...)` when there is no heterogeneous residual alongside it. brms's
# own parameter table for that shape (checked live via brms::get_prior(),
# not assumed) carries `class = "Intercept", dpar = "sigma"` and no
# `class = "b"` row, since offset() contributes no coefficient -- a
# different target from the heterogeneous-residual case's `class = "b"`
# (its `0 + f` has no intercept). Both reasons are handled by the same
# retarget below; when both are present at once (a per-section residual
# on data that also carries observation weights) the `0 + f` shape wins,
# so the target stays `class = "b"`.
#
# THE REPLACEMENT IS A CHOICE, AND IT IS ANNOUNCED. The default uniform(0, U)
# on the residual SD has no exact counterpart on a log-scale linear predictor,
# so it cannot be carried across unchanged; any prior here is a new decision
# rather than a translation. brms's own default for these coefficients is
# flat, which is not a decision this package is willing to make silently
# either. So the scale that produced U is reused -- the uniform default sets
# U = .FB_UNIFORM_SD_MULTIPLIER * sd(y) on the identity scale, hence
# sd(y) = U / .FB_UNIFORM_SD_MULTIPLIER -- and each log-sigma coefficient gets
# normal(log(sd(y)), 1): a prior whose median residual SD is sd(y) and whose
# 95% range spans roughly a factor of seven either side. Weakly informative,
# on the data's own scale, and stated in the prior note rather than inferred
# from behaviour.
#
# The divisor was written here as a literal 2.5 until 0.9.0, which was the
# superseded PC default's multiplier rather than the uniform default's 5. The
# location it produced was log(2 * sd(y)), so the prior's median residual SD
# was twice the data SD instead of the sd(y) this paragraph claims. Reading
# the multiplier off the constant that sets U keeps the inversion exact and
# the two files from drifting again.
#
# On a user-supplied prior the inversion is a heuristic rather than an
# identity: `upper` is whatever numeric the sigma spec ends with, and only a
# uniform(0, U) built by the default machinery is guaranteed to stand in the
# documented relation to sd(y).
.brms_retarget_sigma_for_heterogeneous_residual <- function(specs, fb) {
  het <- Filter(
    function(t) identical(t$type %||% "", "at_units"),
    fb$residual_terms %||% list()
  )
  has_weight_offset <- .fb_weights_nonconstant(.fb_ir_weights(fb))
  if ((length(het) == 0L && !has_weight_offset) || length(specs) == 0L) {
    return(specs)
  }

  # Recover the residual scale from whichever sigma spec was built, so this
  # tracks the default-prior machinery instead of re-deriving it.
  upper <- NA_real_
  for (s in specs) {
    if (identical(s$class %||% "", "sigma")) {
      m <- regmatches(s$string, regexpr("[0-9.eE+-]+\\)$", s$string))
      if (length(m) == 1L) {
        upper <- suppressWarnings(as.numeric(sub("\\)$", "", m)))
      }
    }
  }
  specs <- Filter(function(s) !identical(s$class %||% "", "sigma"), specs)

  loc <- if (is.finite(upper) && upper > 0) {
    log(upper / .FB_UNIFORM_SD_MULTIPLIER)
  } else {
    0
  }
  # `0 + f` (heterogeneous residual, with or without weights alongside)
  # has no intercept -- brms::get_prior() gives it class "b". Weights
  # alone lower to `sigma ~ 1 + offset(...)`, which keeps the intercept
  # and gives class "Intercept" instead (verified live, see the function
  # banner above) -- offset() itself contributes no coefficient to
  # retarget.
  target_class <- if (length(het) >= 1L) "b" else "Intercept"
  c(
    specs,
    list(list(
      string = sprintf("normal(%.6g, 1)", loc),
      class = target_class,
      coef = NA_character_,
      group = NA_character_,
      dpar = "sigma"
    ))
  )
}

# Stack the spec list into a single brms prior object. Each spec
# maps to one brms::prior_string() call; we then rbind. Empty specs
# return NULL so the caller can decide whether to fall back to brms's
# own defaults.
.brms_specs_to_object <- function(specs) {
  if (length(specs) == 0L) {
    return(NULL)
  }
  .check_installed(
    "brms",
    "Package 'brms' is required to compile flexyBayes priors ",
    "into a brms prior object. Install via ",
    "install.packages('brms')."
  )

  rows <- lapply(specs, function(s) {
    args <- list(prior = s$string, class = s$class)
    has_coef <- !is.na(s$coef %||% NA_character_) && nzchar(s$coef)
    if (has_coef) {
      args$coef <- s$coef
    }
    if (!is.na(s$group %||% NA_character_) && nzchar(s$group)) {
      args$group <- s$group
    }
    # `dpar` targets a distributional parameter's own linear predictor --
    # class "b" with dpar "sigma" is the coefficient block created by
    # `sigma ~ 0 + f`, which is a different parameter from the scalar
    # class "sigma". Without this the row would land on the wrong block.
    if (!is.na(s$dpar %||% NA_character_) && nzchar(s$dpar)) {
      args$dpar <- s$dpar
    }
    # brms refuses `Prior argument 'coef' may not be
    # specified when using boundaries` -- when a coef-keyed sd row
    # carries lb/ub the parser rejects the spec. For sd rows the
    # natural sd >= 0 bound is implicit in brms's Stan emit anyway;
    # drop the lb/ub when coef is set.
    if (!has_coef) {
      if (!is.na(s$lb %||% NA_real_)) {
        args$lb <- s$lb
      }
      if (!is.na(s$ub %||% NA_real_)) args$ub <- s$ub
    }
    do.call(brms::prior_string, args)
  })

  do.call(rbind, rows)
}
