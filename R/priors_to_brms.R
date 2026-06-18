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
# Anything outside the supported eight-row table (cor() / smooth() /
# half_cauchy / cauchy / lkj / exponential / gamma targets) raises a
# structured refusal naming both the unsupported spec and the brms-
# native fallback: pass `prior` directly through `...` to brms::brm()
# via fb_brms()'s pass-through `...` argument.

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

  if (!inherits(prior, "fb_prior")) {
    stop(
      "`.priors_to_brms()` expects an `fb_prior` object or NULL ",
      "(legacy bridge). Got: ",
      paste(class(prior), collapse = "/"),
      ".",
      call. = FALSE
    )
  }

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
  specs
}

# Families brms parameterises with a residual `sigma` hyperparameter.
# Mirrors emit_inla's `families_with_prec` discipline -- silently
# dropping sigma priors on families that lack the parameter keeps
# the cross-engine prior surface uniform.
.brms_family_has_sigma <- function(fam) {
  if (is.null(fam)) {
    return(TRUE)
  }
  tolower(as.character(fam)) %in%
    c(
      "gaussian",
      "stdnormal",
      "lognormal",
      "gamma",
      "beta",
      "t",
      "logistic",
      "student",
      "skew_normal"
    )
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
    if (identical(spec$family, "uniform")) {
      lo <- as.numeric(spec$args$lower %||% 0)
      hi <- as.numeric(spec$args$upper)
      string <- sprintf("uniform(%s, %s)", .fmt_num(lo), .fmt_num(hi))
      # brms requires lb/ub for uniform priors so the parser bounds
      # the parameter to the prior's support; sigma must be >= 0.
      return(list(
        string = string,
        class = "sigma",
        coef = NA_character_,
        group = NA_character_,
        lb = max(0, lo),
        ub = hi
      ))
    }
    if (identical(spec$family, "half_normal")) {
      sc_v <- as.numeric(spec$args$scale %||% 1)
      # brms half-normal: normal(0, scale) plus the natural sigma
      # positivity bound (lb = 0).
      string <- sprintf("normal(0, %s)", .fmt_num(sc_v))
      return(list(
        string = string,
        class = "sigma",
        coef = NA_character_,
        group = NA_character_,
        lb = 0,
        ub = NA_real_
      ))
    }
    if (identical(spec$family, "pc")) {
      u <- as.numeric(spec$args$upper %||% 1)
      p <- as.numeric(spec$args$prob %||% 0.01)
      rate <- -log(p) / u
      string <- sprintf("exponential(%s)", .fmt_num(rate))
      return(list(
        string = string,
        class = "sigma",
        coef = NA_character_,
        group = NA_character_,
        lb = 0,
        ub = NA_real_
      ))
    }
    .brms_unsupported(spec$family, "sigma")
  }

  # ------------------------------ sd(group) ------------------------------ #
  if (identical(target$type, "sd")) {
    grp <- target$group
    # Slope-variance row: when the target carries a `coef`
    # slot, the sd row is keyed on (class = "sd", coef = "<x>",
    # group = "<g>") rather than the bare intercept-only
    # (class = "sd", group = "<g>") shape.
    coef_v <- target$coef %||% NA_character_
    if (identical(spec$family, "uniform")) {
      lo <- as.numeric(spec$args$lower %||% 0)
      hi <- as.numeric(spec$args$upper)
      string <- sprintf("uniform(%s, %s)", .fmt_num(lo), .fmt_num(hi))
      return(list(
        string = string,
        class = "sd",
        coef = coef_v,
        group = grp,
        lb = max(0, lo),
        ub = hi
      ))
    }
    if (identical(spec$family, "half_normal")) {
      sc_v <- as.numeric(spec$args$scale %||% 1)
      string <- sprintf("normal(0, %s)", .fmt_num(sc_v))
      return(list(
        string = string,
        class = "sd",
        coef = coef_v,
        group = grp,
        lb = 0,
        ub = NA_real_
      ))
    }
    if (identical(spec$family, "pc")) {
      u <- as.numeric(spec$args$upper %||% 1)
      p <- as.numeric(spec$args$prob %||% 0.01)
      rate <- -log(p) / u
      string <- sprintf("exponential(%s)", .fmt_num(rate))
      return(list(
        string = string,
        class = "sd",
        coef = coef_v,
        group = grp,
        lb = 0,
        ub = NA_real_
      ))
    }
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

  # ------------------------------ cor / smooth / name -------------------- #
  if (target$type %in% c("cor", "smooth", "name")) {
    .brms_unsupported_target(target)
  }

  .brms_unsupported_target(target)
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

# Refusal for an unsupported distribution family on a recognised
# target.
.brms_unsupported <- function(family, target_label) {
  stop(
    "Stan-passthrough emit does not yet translate `",
    family,
    "` priors on ",
    target_label,
    ". Supported targets: ",
    "b()/Intercept ~ normal(mean, sd) | student_t(df, scale); ",
    "sigma ~ uniform(lower, upper) | half_normal(scale) | pc(upper, prob); ",
    "sd(group) ~ uniform | half_normal | pc (intercept-variance row); ",
    "sd(group, coef = <slope>) ~ uniform | half_normal | pc ",
    "(slope-variance row). ",
    "For unsupported priors, pass a brms `prior` object directly ",
    "via `...` (e.g. ",
    "`fb_brms(..., backend = \"brms\", prior = brms::prior(...))` -- ",
    "this bypasses the flexyBayes prior DSL on the Stan path).",
    call. = FALSE
  )
}

# Refusal for an unsupported target (cor() / smooth() / generic name).
.brms_unsupported_target <- function(target) {
  label <- switch(
    target$type %||% "<unknown>",
    "cor" = paste0("cor(group = \"", target$group, "\")"),
    "smooth" = paste0("smooth(\"", target$var, "\")"),
    "name" = paste0("name = \"", target$name, "\""),
    target$type %||% "<unknown>"
  )
  stop(
    "Stan-passthrough emit does not yet handle prior target `",
    label,
    "`. Supported targets: b()/Intercept/sigma/sd(group); ",
    "smoothers and structured covariances are deferred ",
    "to a future release. For full brms-side control, pass a ",
    "brms `prior` object via `...` on the Stan-passthrough call.",
    call. = FALSE
  )
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
      length(fb$rcov_terms %||% list()) > 0L
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
  .brms_specs_to_object(specs)
}

# Stack the spec list into a single brms prior object. Each spec
# maps to one brms::prior_string() call; we then rbind. Empty specs
# return NULL so the caller can decide whether to fall back to brms's
# own defaults.
.brms_specs_to_object <- function(specs) {
  if (length(specs) == 0L) {
    return(NULL)
  }
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "Package 'brms' is required to compile flexyBayes priors ",
      "into a brms prior object. Install via ",
      "install.packages('brms').",
      call. = FALSE
    )
  }

  rows <- lapply(specs, function(s) {
    args <- list(prior = s$string, class = s$class)
    has_coef <- !is.na(s$coef %||% NA_character_) && nzchar(s$coef)
    if (has_coef) {
      args$coef <- s$coef
    }
    if (!is.na(s$group %||% NA_character_) && nzchar(s$group)) {
      args$group <- s$group
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
