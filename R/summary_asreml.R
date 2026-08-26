# summary_asreml.R -- one summary object for every active engine.
#
# Constructs and prints the `summary.flexybayes` object that
# summary.flexybayes() and summary.flexybayes_inla() both return. The
# object is a list of eleven slots -- fixed, varcomp, random, missing,
# converge, n_design, n_observed, na_action, model, engine, call -- so
# that `summary(fit)$varcomp` answers the same question whichever engine
# ran the fit, which it could not before: the two engines returned two
# incomparable objects (INLA's four-slot list, brms's brmssummary) and
# neither carried a variance-component table. A fit carrying an
# autoregressive latent field adds a twelfth, `spatial_field`, which is
# the one engine-native slot and the only conditional one.
#
# Three properties of the table are load-bearing:
#
#   * `estimate` is a posterior mean, not a REML component, and the
#     printed banner says so above every table.
#   * `prior` is a projection of prior_summary(), never a second
#     prior-string builder, so the summary and the prior accessor cannot
#     disagree.
#   * the diagnostics in `converge` are whatever the engine that ran the
#     fit reports. A Laplace approximation has no R-hat, and one is not
#     invented for it.
#
# `$extras$variance_comps` keeps the five column names broom's tidier
# reads (component, estimate, sd, q2.5, q97.5); the renaming to
# estimate / std.error / conf.low / conf.high happens here, on the way
# out, so the broom contract is untouched.

# ---------------------------------------------------------------- #
# Small readers over a fitted object                                #
# ---------------------------------------------------------------- #

# .fb_fit_engine() --- which engine produced this fit.
#
# Read off the slots the object carries rather than its class vector, so
# a fit assembled by hand (the test mocks) answers as well as one that
# came from an emit.
#
# @noRd
# @keywords internal
.fb_fit_engine <- function(object) {
  if (!is.null(object$inla)) {
    return("inla")
  }
  if (!is.null(object$brms)) {
    return("brms")
  }
  if (inherits(object, "flexybayes_inla")) {
    return("inla")
  }
  if (inherits(object, "flexybayes_brms")) {
    return("brms")
  }
  "unknown"
}

# .fb_engine_line() --- the engine, named as what it did.
#
# @noRd
# @keywords internal
.fb_engine_line <- function(engine) {
  switch(
    engine,
    inla = "INLA nested Laplace approximation",
    brms = "Stan HMC via brms",
    "unrecorded"
  )
}

# .fb_fit_sampled() --- did an engine actually draw samples?
#
# The sampler lines (chains, warmup, R-hat, ESS) print only when they
# describe something that happened. INLA reports a deterministic
# approximation and has no chains at all, so printing "MCMC: 4 chains"
# over one is not a cosmetic slip: it tells a reader the fit is
# stochastic when it is not.
#
# @noRd
# @keywords internal
.fb_fit_sampled <- function(object) {
  if (identical(.fb_fit_engine(object), "inla")) {
    return(FALSE)
  }
  ch <- object$extras$call_info$chains
  is.numeric(ch) && length(ch) == 1L && !is.na(ch) && ch > 0
}

# .fb_fit_ir() --- the intermediate representation a fit carries.
#
# The INLA emit keeps it at `$fb`, the shared dispatch path writes it to
# `$extras$fb_terms`, and an older object may carry only the per-side
# term lists under `$extras$parse_info`. All three are the same terms.
#
# @noRd
# @keywords internal
.fb_fit_ir <- function(object) {
  ir <- object$extras$fb_terms %||% object$fb
  if (!is.null(ir)) {
    return(ir)
  }
  pi <- object$extras$parse_info
  if (is.null(pi)) {
    return(NULL)
  }
  list(
    random_terms = pi$random %||% list(),
    residual_terms = pi$residual %||% list()
  )
}


# ---------------------------------------------------------------- #
# The model line                                                    #
# ---------------------------------------------------------------- #
#
# Derived from the IR, never from an engine's own formula string. INLA
# writes `f(row_id, model = "ar1", group = col_id, control.group =
# list(model = "ar1"))`, which is the emitted program and not the model
# the user wrote; parsing it back would make the printed description a
# property of the emitter.

# .fb_term_human_string() --- one IR term as a reader would say it.
#
# @noRd
# @keywords internal
.fb_term_human_string <- function(term) {
  ttype <- term$type %||% "unknown"
  switch(
    ttype,
    simple = paste0(term$var, " iid"),
    ide = paste0(term$var, " iid"),
    units = "units",
    at_units = paste0("units sectioned by ", term$var),
    ar1 = paste0("ar1(", term$var, ") field + nugget"),
    ar1_spatial = paste0(
      "ar1(",
      term$row_var,
      "):ar1(",
      term$col_var,
      ") field + nugget (4 parameters)"
    ),
    spline = paste0("spline in ", term$var),
    smooth_mgcv = paste0("smooth in ", term$var),
    vm = paste0("vm(", term$var, ") known covariance"),
    ped = paste0("ped(", term$var, ") pedigree"),
    nested = paste0(term$inner, " within ", term$outer),
    diag_gxe = paste0("diag(", term$outer, "):", term$inner),
    us_gxe = paste0("us(", term$outer, "):", term$inner),
    corh_gxe = paste0("corh(", term$outer, "):", term$inner),
    fa_gxe = paste0("fa(", term$outer, ", k):", term$inner),
    term$label %||% ttype
  )
}

# .fb_model_line() --- the one-line G-side / R-side description.
#
# @noRd
# @keywords internal
.fb_model_line <- function(ir) {
  if (is.null(ir)) {
    return(NA_character_)
  }
  side <- function(terms, empty) {
    if (length(terms) == 0L) {
      return(empty)
    }
    paste(
      vapply(terms, .fb_term_human_string, character(1L)),
      collapse = "; "
    )
  }
  g <- side(ir$random_terms %||% list(), "none")
  r <- side(ir$residual_terms %||% list(), "units")
  paste0("G: ", g, "; R: ", r)
}

# .fb_agg_line() --- the N-to-K sentence, or NULL on a per-row fit.
#
# One phrasing of the compression story, shared by the summary object's
# model line and the printed header, so the two cannot drift.
#
# @noRd
# @keywords internal
.fb_agg_line <- function(object) {
  am <- object$extras$aggregation_meta
  if (is.null(am) || is.null(am$N) || is.null(am$K)) {
    return(NULL)
  }
  if (!is.finite(am$N) || !is.finite(am$K) || am$K <= 0) {
    return(NULL)
  }
  sprintf(
    "N = %s rows -> K = %s cells (ratio %.0f:1)",
    format(am$N, big.mark = " ", scientific = FALSE),
    format(am$K, big.mark = " ", scientific = FALSE),
    am$N / am$K
  )
}

# .fb_model_line_for_fit() --- the model line a fitted object earns.
#
# The G-side / R-side description is a property of the model, and the
# aggregated representation does not change the model -- it changes the
# rows the engine was handed. That belongs on the same line, because a
# reader comparing the fixed-effect table of an aggregated fit against a
# per-row one otherwise has nothing on the summary object saying the two
# were computed over different numbers of rows.
#
# @noRd
# @keywords internal
.fb_model_line_for_fit <- function(object) {
  model <- .fb_model_line(.fb_fit_ir(object))
  agg <- .fb_agg_line(object)
  if (is.null(agg) || is.na(model)) {
    return(model)
  }
  paste0(model, " [aggregated: ", agg, "]")
}

# .fb_ir_has_field() --- does the model carry an autoregressive field?
#
# @noRd
# @keywords internal
.fb_ir_has_field <- function(ir) {
  if (is.null(ir)) {
    return(FALSE)
  }
  terms_all <- c(ir$random_terms %||% list(), ir$residual_terms %||% list())
  any(vapply(
    terms_all,
    function(t) (t$type %||% "") %in% c("ar1", "ar1_spatial"),
    logical(1L)
  ))
}

# .fb_is_spatial_varcomp() --- does this table describe a latent field?
#
# Read off the canonical component names rather than the IR, so the
# honesty sentence appears exactly when the table it explains does.
#
# @noRd
# @keywords internal
.fb_is_spatial_varcomp <- function(varcomp) {
  if (is.null(varcomp) || nrow(varcomp) == 0L) {
    return(FALSE)
  }
  any(varcomp$component %in% "sd_spatial") ||
    any(startsWith(varcomp$component, "rho_"))
}


# ---------------------------------------------------------------- #
# The variance-component table                                      #
# ---------------------------------------------------------------- #

# .fb_varcomp_priors() --- the `prior` column, projected from
# prior_summary().
#
# One resolved-prior source for the package. The cell for a component the
# package priored is that prior's own canonical string; for a component
# left to the engine the cell says so in two words and prior_summary()
# carries the full sentence, which is too long for a table cell. The
# legacy scalar bridge counts as priored wherever it reaches the engine --
# see .fb_legacy_bridge_priors() below for why that needs saying.
#
# @noRd
# @keywords internal
.fb_varcomp_priors <- function(object, components) {
  na_out <- rep(NA_character_, length(components))
  ps <- tryCatch(prior_summary(object), error = function(e) NULL)
  if (is.null(ps)) {
    return(na_out)
  }

  recorded <- character(0)
  if (identical(ps$kind, "fb_prior") && !is.null(ps$fb_prior)) {
    for (s in ps$fb_prior$specs %||% list()) {
      key <- .fb_prior_spec_parameter(s)
      if (!is.na(key)) {
        recorded[[key]] <- .fb_prior_spec_string(s)
      }
    }
  }
  engine_default <- names(ps$engine_default %||% character(0))
  legacy <- identical(ps$kind, "legacy_scalar")
  bridged <- if (legacy) {
    .fb_legacy_bridge_priors(ps)
  } else {
    character(0)
  }
  # What the bridge priored where there is no engine table to read it
  # off -- the INLA route. prior_summary() carries the record; this is a
  # projection of it, in the order the two sources deserve: the engine's
  # own table first, the declaration second.
  bridge_declared <- if (legacy) {
    ps$legacy_vc_applied %||% character(0)
  } else {
    character(0)
  }
  sectioned <- .fb_sectioned_residual_prior(ps)

  vapply(
    components,
    function(cmp) {
      if (cmp %in% names(recorded)) {
        return(unname(recorded[[cmp]]))
      }
      # A per-level residual row. The declared uniform on the SD scale
      # was retargeted onto the log-sigma coefficients at emit time, so
      # `sigma` names no parameter this model has and the spec keyed to
      # it would be the wrong cell. What reached the sampler is read off
      # the engine's own table, exactly as for the legacy bridge.
      if (nzchar(sectioned) && startsWith(cmp, "sigma_")) {
        return(sectioned)
      }
      # Ahead of the engine-default branch on purpose: a component the
      # legacy bridge priored is listed as an engine default by
      # .fb_prior_record(), which classifies for triangulate()'s
      # matched-prior gate rather than for this cell.
      if (cmp %in% names(bridged)) {
        return(unname(bridged[[cmp]]))
      }
      if (cmp %in% names(bridge_declared)) {
        return(unname(bridge_declared[[cmp]]))
      }
      if (cmp %in% engine_default) {
        return("engine default")
      }
      if (legacy) {
        return("legacy scalar bridge")
      }
      # Nothing recorded, nothing bridged, and the walker did not list it
      # among the parameters it hands to the engine on purpose. That last
      # gap is reachable: a partial user prior -- `fb_prior(sigma ~ ...)`
      # on a model with a random term the prior never names -- leaves the
      # term to the engine's own hyperprior and appears in neither list,
      # so the cell rendered blank where the engine default was in force.
      # It is the same component class the branch above names, reached by
      # a different route, and it gets the same two words. A blank cell
      # tells the reader nothing, and both surfaces that check this table
      # already assert no cell is NA.
      "engine default"
    },
    character(1L),
    USE.NAMES = FALSE
  )
}

# .fb_legacy_bridge_priors() --- what the legacy scalar bridge actually
# put on each variance component.
#
# `engine default` is the wrong cell for a bridged component on the Stan
# route, and it is the cell the classification alone produces:
# .fb_prior_record() marks every variance component `engine_default` under
# the bridge because it is answering triangulate()'s question -- do two
# engines share a prior -- and under the bridge they do not. Projected
# into a single-engine table the same word says this package set no
# prior, which on Stan is false. .brms_legacy_specs() writes
# lognormal(0, prior_vc_sd) onto sigma and onto every named sd group.
#
# On the INLA route there is no engine prior table to read, so this
# returns nothing and the branch below answers from the recorded prior --
# which since 0.9.2 carries the same lognormal, because the scalar now
# reaches INLA as the expression prior that writes that density on the SD
# scale. Before then the bridge handed INLA nothing and `engine default`
# was the true cell there.
#
# The strings are read off the engine's own prior table rather than
# rebuilt from the scalar, so the cell states what reached the sampler and
# this file stays a projection of prior_summary() rather than a second
# prior-string builder. A table whose shape this does not recognise
# yields nothing and the cell falls back to what it said before.
#
# @noRd
# @keywords internal
.fb_legacy_bridge_priors <- function(ps) {
  tab <- ps$engine_prior_table
  needed <- c("prior", "class", "coef", "group")
  if (!is.data.frame(tab) || nrow(tab) == 0L) {
    return(character(0))
  }
  if (!all(needed %in% names(tab))) {
    return(character(0))
  }

  # Rows brms attributes to us. A default row is the engine's own choice
  # and belongs in the engine-default branch, not here.
  keep <- nzchar(as.character(tab$prior))
  if ("source" %in% names(tab)) {
    keep <- keep & as.character(tab$source) %in% "user"
  }
  tab <- tab[keep, , drop = FALSE]

  cls <- as.character(tab$class)
  coefs <- as.character(tab$coef)
  groups <- as.character(tab$group)
  strings <- as.character(tab$prior)

  out <- character(0)
  for (i in seq_len(nrow(tab))) {
    # A coef-keyed row prices one coefficient of a group rather than the
    # group's own scale, and the variance-component table has no cell for
    # it.
    if (nzchar(coefs[[i]])) {
      next
    }
    if (identical(cls[[i]], "sigma")) {
      out[["sigma"]] <- strings[[i]]
    } else if (identical(cls[[i]], "sd") && nzchar(groups[[i]])) {
      out[[paste0("sd_", groups[[i]])]] <- strings[[i]]
    }
  }
  out
}

# .fb_sectioned_residual_prior() --- the prior the per-level residuals
# actually ran under.
#
# `dsum(~ units | f)` makes the residual a distributional predictor, and
# the package retargets the declared uniform-on-SD onto the log-sigma
# coefficients rather than leaving a spec pointing at a `sigma` the model
# no longer has. The retargeted spec is what reached the sampler, so it
# is read back off the engine's own prior table -- the same discipline
# .fb_legacy_bridge_priors() follows, and for the same reason: this file
# projects prior_summary(), it does not rebuild prior strings.
#
# Returns "" when the fit has no such row, which leaves the cell to the
# branches below it.
#
# @noRd
# @keywords internal
.fb_sectioned_residual_prior <- function(ps) {
  tab <- ps$engine_prior_table
  needed <- c("prior", "class", "dpar")
  if (!is.data.frame(tab) || nrow(tab) == 0L || !all(needed %in% names(tab))) {
    return("")
  }
  keep <- as.character(tab$class) == "b" &
    as.character(tab$dpar) == "sigma" &
    nzchar(as.character(tab$prior))
  if ("coef" %in% names(tab)) {
    keep <- keep & !nzchar(as.character(tab$coef))
  }
  if ("source" %in% names(tab)) {
    keep <- keep & as.character(tab$source) %in% "user"
  }
  keep[is.na(keep)] <- FALSE
  if (!any(keep)) {
    return("")
  }
  as.character(tab$prior)[which(keep)[[1L]]]
}

# .fb_varcomp_notes() --- the boundary-collapse display flag.
#
# One pinned heuristic, and a display heuristic only: a component is
# flagged `collapsed` when its 97.5% quantile on the SD scale sits below
# 1% of the posterior-median residual SD, which is where the posterior
# has piled against zero and the component is doing no work. It is not a
# test, it has no null hypothesis, and nothing downstream reads it.
#
# The quantiles come from the same marginals or draws the row itself was
# summarised from, so the flag cannot disagree with the row above it.
# Correlations are skipped: they are not on a scale where "small" means
# collapsed. With no residual scale in the model -- a sectioned residual
# leaves none -- there is nothing to compare against and every cell is
# blank.
#
# The latent field carries the wider spatial threshold. A collapsed field
# does not merely sit near zero, it hands its variance to the nugget it
# is being compared against, so the same 1% cut that suits an independent
# grouping factor reads a genuinely lost field as healthy: measured on a
# holed grid, a field SD of 0.0099 beside a nugget of 0.92 -- a field the
# fit had plainly not identified -- cleared the 1% bar on its upper bound
# and the cell stayed blank while the fit-time warning fired. The two
# surfaces now agree, which is the point of having both.
#
# @noRd
# @keywords internal
.FB_COLLAPSE_FRACTION <- 0.01

.fb_varcomp_notes <- function(components, upper, medians) {
  blank <- rep("", length(components))
  ref <- if ("sigma" %in% names(medians)) {
    medians[["sigma"]]
  } else {
    NA_real_
  }
  if (!is.numeric(ref) || length(ref) != 1L || is.na(ref) || ref <= 0) {
    return(blank)
  }
  is_correlation <- startsWith(components, "rho_") |
    startsWith(components, "cor_")
  fraction <- ifelse(
    components == "sd_spatial",
    .FB_SPATIAL_COLLAPSE_FRACTION,
    .FB_COLLAPSE_FRACTION
  )
  flagged <- !is_correlation &
    is.finite(upper) &
    components != "sigma" &
    upper < fraction * ref
  ifelse(flagged, "collapsed", "")
}

# .fb_agg_varcomp_from_means() --- the last-resort aggregated table.
#
# The aggregated emits record their variance components as posterior
# MEANS only -- `$extras$summary$sigma_means` and `$tau_means` -- with no
# dispersion and no quantiles. On the INLA route that is not a loss,
# because the raw fit carries the hyperparameter marginals and
# .inla_variance_comps() reads a full row off them. On a route that kept
# nothing but the means there is no honest interval to report, so the
# three interval columns are NA rather than a spread invented from a
# point estimate, and the caller's `note` column says which rows they
# are.
#
# Component names come from the recorded names where the emit kept them
# and from the model representation otherwise, so a row answers to the
# same canonical name it would carry on the per-row route.
#
# @noRd
# @keywords internal
.fb_agg_varcomp_from_means <- function(object) {
  ps <- object$extras$summary
  if (is.null(ps)) {
    return(NULL)
  }
  sigma <- ps$sigma_means %||% numeric(0)
  tau <- ps$tau_means %||% numeric(0)

  sigma_nm <- names(sigma)
  if (is.null(sigma_nm) || any(!nzchar(sigma_nm))) {
    sigma_nm <- if (length(sigma) == 1L) {
      "sigma"
    } else {
      paste0("sigma_", seq_along(sigma))
    }
  }

  tau_nm <- names(tau)
  if (is.null(tau_nm) || any(!nzchar(tau_nm))) {
    terms_g <- .fb_fit_ir(object)$random_terms %||% list()
    tau_nm <- if (length(terms_g) == length(tau)) {
      vapply(
        terms_g,
        function(t) paste0("sd_", t$var %||% t$label %||% ""),
        character(1L)
      )
    } else {
      character(0)
    }
    if (length(tau_nm) != length(tau) || any(!nzchar(tau_nm))) {
      tau_nm <- paste0("sd_", seq_along(tau))
    }
  }

  component <- c(sigma_nm, tau_nm)
  if (length(component) == 0L) {
    return(NULL)
  }
  data.frame(
    component = component,
    estimate = c(as.numeric(sigma), as.numeric(tau)),
    sd = NA_real_,
    q2.5 = NA_real_,
    q97.5 = NA_real_,
    stringsAsFactors = FALSE
  )
}

# .fb_variance_comps() --- the stored table, whichever route wrote it.
#
# The per-row emits write `$extras$variance_comps` as the five-column
# broom-shaped data frame. The aggregated emits write a two-element list
# of posterior means under the same name, which is why every reader of
# that field -- summary(), tidy(effects = "random"), plot(type =
# "variance") -- errored on an aggregated fit rather than returning
# something. Aggregation is the DEFAULT route for an ordinary Gaussian
# call with a random term, so that was the ordinary case, not a corner.
#
# One reader normalises the two shapes. An aggregated INLA fit is routed
# to the same marginal-transform builder the per-row INLA route uses, so
# the aggregated table is the engine's own posterior on the
# standard-deviation scale rather than a reciprocal square root of a
# tabulated mean. Anything else falls back to the recorded means with
# absent intervals.
#
# @noRd
# @keywords internal
.fb_variance_comps <- function(object) {
  vc <- object$extras$variance_comps
  if (is.data.frame(vc)) {
    return(vc)
  }
  # Only the aggregated representation is rebuilt here. A per-row fit
  # whose emit-time build failed carries a NULL field on purpose, warned
  # about at the time, and re-running the build inside summary() would
  # re-raise that warning on every print.
  if (!inherits(object, "flexybayes_aggregated")) {
    return(NULL)
  }

  if (identical(.fb_fit_engine(object), "inla") && !is.null(object$inla)) {
    built <- tryCatch(.inla_variance_comps(object), error = function(e) NULL)
    if (is.data.frame(built) && nrow(built) > 0L) {
      return(built)
    }
  }

  .fb_agg_varcomp_from_means(object)
}

# .fb_varcomp_residual_by_level() --- the sectioned residual's rows.
#
# A `dsum(~ units | f)` residual leaves no scalar `sigma` behind. What the
# model carries instead is one coefficient on LOG sigma per level of f,
# which the emit-time variance-component builder does not collect -- it
# reads `sd_*` rows and `sigma`, and `b_sigma_*` is neither. So the entire
# point of fitting dsum() was printed by summary() and absent from the
# object it returned: `summary(fit)$varcomp` came back with the grouping
# factors and no residual at all, and reaching the per-level values meant
# opening `$brms` and transforming the draws by hand.
#
# The rows are taken from the same builder the printed block reads, so
# the table and the panel cannot state different numbers for one fit.
# That also fixes the point summary: the printed panel reports posterior
# MEDIANS, because exp() is convex and the median is the one summary a
# monotone transform carries through exactly. These rows therefore carry
# the median where the sampled `sd_*` rows carry a posterior mean. The
# alternative -- a mean the panel never displays -- would reintroduce the
# disagreement this closes.
#
# Naming follows the canonical vocabulary: `sigma_<level>`, on the
# standard-deviation scale like every other row of the table.
#
# C6 (observation weights) hits the same wall by a different door: with
# no sectioned residual but weights present, brms's sigma is `~ 1 +
# offset(...)` rather than a scalar -- one Intercept coefficient, not
# one per level. .brms_weights_sigma_row() recovers that single row;
# it is folded in here, named plainly "sigma" (there is no level to
# suffix), so a weighted fit's varcomp table reads exactly like an
# unweighted one's.
#
# NULL on every fit without a sectioned residual and without weights,
# which is the signal to leave the table alone.
#
# @noRd
# @keywords internal
.fb_varcomp_residual_by_level <- function(object) {
  tab <- tryCatch(
    .brms_residual_by_level_table(object),
    error = function(e) NULL
  )
  if (is.data.frame(tab) && nrow(tab) > 0L) {
    return(data.frame(
      component = paste0("sigma_", as.character(tab$level)),
      estimate = as.numeric(tab$sd),
      std.error = as.numeric(tab$sd_se %||% NA_real_),
      conf.low = as.numeric(tab$sd_lower),
      conf.high = as.numeric(tab$sd_upper),
      stringsAsFactors = FALSE
    ))
  }
  wrow <- tryCatch(.brms_weights_sigma_row(object), error = function(e) NULL)
  if (!is.data.frame(wrow) || nrow(wrow) == 0L) {
    return(NULL)
  }
  data.frame(
    component = "sigma",
    estimate = as.numeric(wrow$sd),
    std.error = as.numeric(wrow$sd_se %||% NA_real_),
    conf.low = as.numeric(wrow$sd_lower),
    conf.high = as.numeric(wrow$sd_upper),
    stringsAsFactors = FALSE
  )
}

# .fb_summary_varcomp() --- `$extras$variance_comps` as `$varcomp`.
#
# The five stored columns are renamed on the way out (sd -> std.error,
# q2.5 -> conf.low, q97.5 -> conf.high) and the two derived columns are
# appended. Nothing is renamed in place: the stored table is a broom
# contract and six test files plus three internal callers read it.
#
# @noRd
# @keywords internal
.fb_summary_varcomp <- function(object) {
  vc <- .fb_variance_comps(object)
  empty <- data.frame(
    component = character(0),
    estimate = numeric(0),
    std.error = numeric(0),
    conf.low = numeric(0),
    conf.high = numeric(0),
    prior = character(0),
    note = character(0),
    stringsAsFactors = FALSE
  )
  resid_rows <- .fb_varcomp_residual_by_level(object)
  if (is.null(vc) || nrow(vc) == 0L) {
    if (is.null(resid_rows)) {
      return(empty)
    }
    vc <- data.frame(
      component = character(0),
      estimate = numeric(0),
      sd = numeric(0),
      q2.5 = numeric(0),
      q97.5 = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  medians <- attr(vc, "posterior_median") %||% numeric(0)
  out <- data.frame(
    component = as.character(vc$component),
    estimate = as.numeric(vc$estimate),
    std.error = as.numeric(vc$sd),
    conf.low = as.numeric(vc$q2.5),
    conf.high = as.numeric(vc$q97.5),
    stringsAsFactors = FALSE
  )
  if (!is.null(resid_rows)) {
    out <- rbind(out, resid_rows)
  }
  out$prior <- .fb_varcomp_priors(object, out$component)
  out$note <- .fb_varcomp_notes(out$component, out$conf.high, medians)
  # A row the source could summarise to a point and no further says so.
  # Blank cells beside three NAs read as a rendering slip; the words say
  # the fit recorded no interval for that component.
  no_interval <- is.na(out$std.error) &
    is.na(out$conf.low) &
    is.na(out$conf.high) &
    !nzchar(out$note)
  out$note[no_interval] <- "no interval recorded"
  rownames(out) <- NULL
  out
}


# ---------------------------------------------------------------- #
# The remaining slots                                               #
# ---------------------------------------------------------------- #

# .fb_summary_fixed() --- the fixed-effect table.
#
# Delegated to tidy(), which already answers per engine from the right
# place: INLA's own marginal summary on an INLA fit, coef() / vcov() /
# confint() on the sampled engines. One route, so the summary and the
# tidier cannot report different numbers for the same coefficient.
#
# @noRd
# @keywords internal
.fb_summary_fixed <- function(object) {
  out <- tryCatch(
    tidy(object, conf.int = TRUE),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(data.frame(
      term = character(0),
      estimate = numeric(0),
      std.error = numeric(0),
      conf.low = numeric(0),
      conf.high = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------- #
# The random-effect table                                           #
# ---------------------------------------------------------------- #
#
# One data frame per grouping factor, in the same six columns on every
# engine: group, level, estimate, std.error, conf.low, conf.high. INLA
# reports a data frame per term keyed by `ID`; brms a three-dimensional
# array per grouping factor. A fit whose engine `.fb_fit_engine()`
# cannot identify falls through to `.fb_random_blups()`, which reads
# whatever `object$extras$blups` holds (a bare named vector of
# posterior means, with no uncertainty columns, on the one record shape
# that ever populated it -- see that function) -- and a reader
# comparing two fits should not have to learn the engine-specific
# shapes to use the table.
#
# A grouping factor carrying more than one effect (a random slope beside
# its intercept) contributes one row per effect and level, with the
# effect named in `level`. With a single effect, which is the ordinary
# case, `level` is the level label alone.

.FB_RANDOM_COLS <- c(
  "group",
  "level",
  "estimate",
  "std.error",
  "conf.low",
  "conf.high"
)

# .fb_random_empty() --- the typed zero-row frame.
#
# @noRd
# @keywords internal
.fb_random_empty <- function() {
  data.frame(
    group = character(0),
    level = character(0),
    estimate = numeric(0),
    std.error = numeric(0),
    conf.low = numeric(0),
    conf.high = numeric(0),
    stringsAsFactors = FALSE
  )
}

# .fb_random_rows() --- one block in the frozen column order.
#
# @noRd
# @keywords internal
.fb_random_rows <- function(
  group,
  level,
  estimate,
  std_error,
  conf_low,
  conf_high
) {
  level <- as.character(level)
  if (length(level) == 0L) {
    return(.fb_random_empty())
  }
  data.frame(
    group = rep(as.character(group), length(level)),
    level = level,
    estimate = as.numeric(estimate),
    std.error = as.numeric(std_error),
    conf.low = as.numeric(conf_low),
    conf.high = as.numeric(conf_high),
    stringsAsFactors = FALSE
  )
}

# .fb_random_column() --- one column of an engine table, or NA.
#
# INLA names its quantile columns "0.025quant" / "0.975quant", which are
# not syntactic, and omits them entirely when the fit was run without
# marginals. Reading them by name with a typed fallback keeps a missing
# column out of the arithmetic.
#
# @noRd
# @keywords internal
.fb_random_column <- function(tab, nm) {
  if (nm %in% names(tab)) {
    return(as.numeric(tab[[nm]]))
  }
  rep(NA_real_, nrow(tab))
}

# .fb_random_inla() --- INLA's per-term summaries, lowered.
#
# @noRd
# @keywords internal
.fb_random_inla <- function(object) {
  sr <- object$inla$summary.random
  if (!is.list(sr) || length(sr) == 0L) {
    return(list())
  }
  out <- lapply(names(sr), function(nm) {
    tab <- sr[[nm]]
    if (!is.data.frame(tab) || nrow(tab) == 0L) {
      return(NULL)
    }
    # C4/FS-26: `tab$ID` carries the legalised level ("X1", "low.N", ...)
    # for a group INLA's f() indexes by a relabelled factor; restore the
    # user's own label the same way coef.flexybayes_inla() does for the
    # fixed table, so ranef() -- which reaches this function via
    # coef(what = "random") -- prints what the user wrote.
    level_vals <- tab$ID %||% seq_len(nrow(tab))
    if (is.character(level_vals)) {
      level_vals <- .inla_restore_level_labels(
        object$level_labels,
        nm,
        level_vals
      )
    }
    .fb_random_rows(
      group = nm,
      level = level_vals,
      estimate = .fb_random_column(tab, "mean"),
      std_error = .fb_random_column(tab, "sd"),
      conf_low = .fb_random_column(tab, "0.025quant"),
      conf_high = .fb_random_column(tab, "0.975quant")
    )
  })
  names(out) <- names(sr)
  out[!vapply(out, is.null, logical(1L))]
}

# .fb_random_brms() --- brms::ranef()'s arrays, lowered.
#
# `brms::ranef()` returns one `levels x statistic x effect` array per
# grouping factor, with the statistics named Estimate / Est.Error /
# Q2.5 / Q97.5 (verified against brms 2.23.0 on a live fit). The array
# is stacked effect by effect so the result is one data frame per
# grouping factor whatever the effect count.
#
# @noRd
# @keywords internal
.fb_random_brms <- function(object) {
  if (is.null(object$brms) || !requireNamespace("brms", quietly = TRUE)) {
    return(list())
  }
  re <- tryCatch(brms::ranef(object$brms), error = function(e) NULL)
  if (!is.list(re) || length(re) == 0L) {
    return(list())
  }

  out <- lapply(names(re), function(nm) {
    arr <- re[[nm]]
    dn <- dimnames(arr)
    if (length(dn) != 3L || length(dn[[1L]]) == 0L) {
      return(NULL)
    }
    levs <- dn[[1L]]
    effects <- dn[[3L]]
    stat <- function(s, k) {
      if (s %in% dn[[2L]]) {
        return(as.numeric(arr[, s, k]))
      }
      rep(NA_real_, length(levs))
    }
    blocks <- lapply(effects, function(k) {
      lab <- if (length(effects) > 1L) paste0(levs, ":", k) else levs
      .fb_random_rows(
        group = nm,
        level = lab,
        estimate = stat("Estimate", k),
        std_error = stat("Est.Error", k),
        conf_low = stat("Q2.5", k),
        conf_high = stat("Q97.5", k)
      )
    })
    do.call(rbind, blocks)
  })
  names(out) <- names(re)
  out[!vapply(out, is.null, logical(1L))]
}

# .fb_random_blups() --- the fallback `.fb_summary_random()` reaches
# when `.fb_fit_engine()` cannot identify the fit as inla or brms.
#
# Reads `object$extras$blups`, a bare named vector of posterior means
# with no uncertainty attached (the shape a since-withdrawn native
# engine populated; see NEWS.md, 0.9.3). Neither active engine
# populates this slot, so on an inla/brms fit this path is never
# reached; it remains as the typed no-op for a fit
# `.fb_fit_engine()` cannot identify, returning `list()` rather than
# erroring on a missing slot. The three uncertainty columns are `NA`
# rather than filled with a number the object never held.
#
# @noRd
# @keywords internal
.fb_random_blups <- function(object) {
  bl <- object$extras$blups
  if (!is.list(bl) || length(bl) == 0L) {
    return(list())
  }
  out <- lapply(names(bl), function(nm) {
    v <- bl[[nm]]
    if (length(v) == 0L) {
      return(NULL)
    }
    .fb_random_rows(
      group = nm,
      level = names(v) %||% as.character(seq_along(v)),
      estimate = as.numeric(v),
      std_error = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_
    )
  })
  names(out) <- names(bl)
  out[!vapply(out, is.null, logical(1L))]
}

# .fb_summary_random() --- the random-effect table, one shape per engine.
#
# The single construction behind `summary(fit)$random` and
# `coef(fit, what = "random")`, so the two accessors cannot report
# different level sets for the same fit.
#
# @noRd
# @keywords internal
.fb_summary_random <- function(object) {
  engine <- .fb_fit_engine(object)
  out <- tryCatch(
    switch(
      engine,
      inla = .fb_random_inla(object),
      brms = .fb_random_brms(object),
      .fb_random_blups(object)
    ),
    error = function(e) {
      warning(
        "flexyBayes: the random-effect table could not be built from this ",
        "fit (",
        conditionMessage(e),
        "); returning none.",
        call. = FALSE
      )
      list()
    }
  )
  if (length(out) == 0L) {
    return(list())
  }
  out
}


# ---------------------------------------------------------------- #
# The coefficient accessor's shared body                            #
# ---------------------------------------------------------------- #

# .fb_coef_what() --- resolve `coef(what = )` against one fit.
#
# The fixed vector is passed in rather than recomputed, because the two
# engines reach it by different routes: a brms fit keeps it on
# `$glm$coefficients` and an INLA fit has no `$glm` slot at all. Every
# other branch is engine-independent, and lives here so the two coef()
# methods cannot answer `what = "random"` differently.
#
# @noRd
# @keywords internal
.fb_coef_what <- function(object, what, fixed) {
  switch(
    what,
    fixed = fixed,
    random = .fb_summary_random(object),
    missing = .fb_summary_missing(object),
    all = list(
      fixed = fixed,
      random = .fb_summary_random(object),
      missing = .fb_summary_missing(object)
    )
  )
}

# ---------------------------------------------------------------- #
# The unobserved design cells                                       #
# ---------------------------------------------------------------- #
#
# ASReml's `mv` factor, in this package's terms: the posterior for each
# response the design says exists and nobody observed. Both engines
# already compute it -- INLA treats an NA response as a latent
# prediction target, brms samples it through `mi()` -- so this reads what
# the fit holds rather than adding a second call to either engine.
#
# The five frozen columns come first (`row`, `estimate`, `std.error`,
# `conf.low`, `conf.high`); the design index variables the fit recorded
# are appended after them, so the column positions of the frozen five do
# not move with the model.

# .fb_missing_empty() --- the typed zero-row frame.
#
# @noRd
# @keywords internal
.fb_missing_empty <- function() {
  data.frame(
    row = integer(0),
    estimate = numeric(0),
    std.error = numeric(0),
    conf.low = numeric(0),
    conf.high = numeric(0),
    stringsAsFactors = FALSE
  )
}

# .fb_missing_rows() --- which rows of the fitted data carry no response.
#
# The engine was handed the augmented design, so these indices are rows
# of the data the fit carries, which is the same numbering the brms
# `Ymi[i]` parameter names use.
#
# @noRd
# @keywords internal
.fb_missing_rows <- function(object) {
  dat <- .fb_fit_data(object)
  resp <- .fb_fit_ir(object)$response %||%
    object$extras$parse_info$fixed$response
  if (is.null(dat) || is.null(resp) || !resp %in% names(dat)) {
    return(integer(0))
  }
  which(is.na(dat[[resp]]))
}

# .fb_missing_index_cols() --- the design index columns at those rows.
#
# @noRd
# @keywords internal
.fb_missing_index_cols <- function(object, rows) {
  vars <- object$extras$na_action$design_index_vars %||% character(0)
  dat <- .fb_fit_data(object)
  if (length(vars) == 0L || is.null(dat)) {
    return(NULL)
  }
  vars <- vars[vars %in% names(dat)]
  # An index variable sharing a name with one of the frozen five would
  # collide on the way into the frame; it is suffixed rather than
  # dropped, so the value is still reachable.
  frozen <- names(.fb_missing_empty())
  out <- lapply(vars, function(v) dat[[v]][rows])
  names(out) <- ifelse(vars %in% frozen, paste0(vars, ".index"), vars)
  if (length(out) == 0L) {
    return(NULL)
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

# .fb_missing_inla() --- INLA's fitted-value summary at the NA rows.
#
# Read from `$inla$summary.fitted.values`, which carries one row per
# observation. `$extras$summary$fitted` is only `head()` of the same
# table and would silently answer for the first six rows alone.
#
# Restricted to the identity link. INLA says so itself on any fit with a
# missing response and another link -- "otherwise the identity link will
# be used to compute the fitted values for NA data" -- so on a Poisson or
# binomial fit the rows for the unobserved cells are on the
# linear-predictor scale while the rows around them are on the response
# scale. Reporting the two under one column heading would put a number in
# the wrong dress, so the table is withheld and the reason said out loud.
# Making it available means passing `control.predictor = list(link = 1)`
# at emit time, which changes the engine call for every INLA fit.
#
# @noRd
# @keywords internal
.fb_missing_inla <- function(object, rows) {
  link <- object$extras$model_info$link %||% "identity"
  if (!identical(link, "identity")) {
    .fb_warn_missing_link_scale(link)
    return(NULL)
  }
  sfv <- object$inla$summary.fitted.values
  if (is.null(sfv) || nrow(sfv) == 0L) {
    return(NULL)
  }
  rows <- rows[rows <= nrow(sfv)]
  if (length(rows) == 0L) {
    return(NULL)
  }
  tab <- sfv[rows, , drop = FALSE]
  data.frame(
    row = as.integer(rows),
    estimate = .fb_random_column(tab, "mean"),
    std.error = .fb_random_column(tab, "sd"),
    conf.low = .fb_random_column(tab, "0.025quant"),
    conf.high = .fb_random_column(tab, "0.975quant"),
    stringsAsFactors = FALSE
  )
}

# .fb_warn_missing_link_scale() --- say why the table is empty, once.
#
# Once per session in the manner of the high-missingness note: a script
# summarising a hundred fits of the same family should say this once.
#
# @noRd
# @keywords internal
.fb_warn_missing_link_scale <- function(link) {
  if (.emit_state_get("missing_link_scale_warning")) {
    return(invisible(NULL))
  }
  .emit_state_set("missing_link_scale_warning", TRUE)
  warning(
    "flexyBayes: the unobserved-cell table is empty on this fit even ",
    "though it carries missing responses. INLA computes a fitted value ",
    "for an NA response on the identity link, and this fit uses the '",
    link,
    "' link, so those rows would be on the linear-predictor scale ",
    "while every row around them is on the response scale. The count of ",
    "unobserved cells is still nobs(fit) - nobs(fit, type = ",
    "\"observed\"), and the posterior itself is on the fit at ",
    "$inla$summary.fitted.values.",
    call. = FALSE
  )
  invisible(NULL)
}

# .fb_missing_brms() --- the `Ymi` draws a Gaussian mi() fit carries.
#
# brms names each sampled missing response `Ymi[i]`, where `i` is the row
# of the data the fit was handed (verified against brms 2.23.0 on a live
# fit: a 48-row frame with responses missing at rows 1, 4, 23, 34, 39 and
# 43 produced exactly those six parameter names). A non-Gaussian fit has
# no `mi()` term and therefore no such parameter, which is why the
# missing-response layer refuses that combination rather than deleting
# rows: there is nothing to read here in that case, and nothing was fit.
#
# @noRd
# @keywords internal
.fb_missing_brms <- function(object) {
  if (is.null(object$brms) || !requireNamespace("brms", quietly = TRUE)) {
    return(NULL)
  }
  ps <- tryCatch(brms::posterior_summary(object$brms), error = function(e) NULL)
  if (is.null(ps) || nrow(ps) == 0L) {
    return(NULL)
  }
  keep <- grep("^Ymi\\[[0-9]+\\]$", rownames(ps))
  if (length(keep) == 0L) {
    return(NULL)
  }
  ps <- ps[keep, , drop = FALSE]
  idx <- as.integer(sub("^Ymi\\[([0-9]+)\\]$", "\\1", rownames(ps)))
  ord <- order(idx)
  data.frame(
    row = idx[ord],
    estimate = as.numeric(ps[ord, "Estimate"]),
    std.error = as.numeric(ps[ord, "Est.Error"]),
    conf.low = as.numeric(ps[ord, "Q2.5"]),
    conf.high = as.numeric(ps[ord, "Q97.5"]),
    stringsAsFactors = FALSE
  )
}

# .fb_summary_missing() --- the unobserved design cells.
#
# The single construction behind `summary(fit)$missing` and
# `coef(fit, what = "missing")`. Zero rows -- typed, never `NULL` --
# whenever every response was observed, or the engine holds no posterior
# for the ones that were not.
#
# @noRd
# @keywords internal
.fb_summary_missing <- function(object) {
  out <- tryCatch(
    {
      rows <- .fb_missing_rows(object)
      body <- if (length(rows) == 0L) {
        NULL
      } else if (identical(.fb_fit_engine(object), "inla")) {
        .fb_missing_inla(object, rows)
      } else {
        .fb_missing_brms(object)
      }
      if (is.null(body) || nrow(body) == 0L) {
        NULL
      } else {
        idx <- .fb_missing_index_cols(object, body$row)
        if (is.null(idx)) body else cbind(body, idx)
      }
    },
    error = function(e) {
      warning(
        "flexyBayes: the missing-cell table could not be built from this ",
        "fit (",
        conditionMessage(e),
        "); returning none.",
        call. = FALSE
      )
      NULL
    }
  )
  if (is.null(out)) {
    return(.fb_missing_empty())
  }
  rownames(out) <- NULL
  out
}

# .fb_summary_converge() --- diagnostics in the engine's own terms.
#
# A sampled fit reports R-hat, effective sample sizes and divergent
# transitions. A Laplace approximation reports whether its mode search
# converged, whether the marginal likelihood is finite, and the largest
# symmetric Kullback-Leibler divergence between each Gaussian
# approximation and its corrected marginal -- INLA's own accuracy signal.
# Neither set is translated into the other's vocabulary.
#
# @noRd
# @keywords internal
.fb_summary_converge <- function(object) {
  engine <- .fb_fit_engine(object)

  if (identical(engine, "inla")) {
    kld <- c(
      object$inla$summary.fixed$kld,
      unlist(lapply(
        object$inla$summary.random %||% list(),
        function(r) r$kld
      ))
    )
    kld <- kld[is.finite(kld)]
    mlik <- if (is.null(object$inla$mlik)) {
      NA_real_
    } else {
      object$inla$mlik[1L, 1L]
    }
    # NA, not FALSE, where no numerical confirm ran. The aggregated emit
    # records no `num_check` at all, and isTRUE(NULL) is FALSE, so the
    # slot reported a failed check on a fit that never had one and the
    # printed line said "Numerical confirm: FAIL". Absent and failed are
    # different answers, and only one of them is true here.
    confirm <- if (is.null(object$num_check)) {
      NA
    } else {
      isTRUE(object$num_check$pass)
    }
    return(list(
      engine = "INLA nested Laplace approximation",
      mode_status = object$inla$mode$mode.status %||% NA_real_,
      mlik = mlik,
      kld_max = if (length(kld) == 0L) NA_real_ else max(kld),
      numerical_confirm = confirm,
      numerical_confirm_reasons = object$num_check$reasons %||% character(0)
    ))
  }

  conv <- object$extras$convergence
  rhat <- if (!is.null(conv$gelman)) {
    conv$gelman$psrf[, "Point est."]
  } else {
    numeric(0)
  }
  rhat <- rhat[is.finite(rhat)]
  ess_bulk <- conv$n_eff %||% numeric(0)
  ess_bulk <- ess_bulk[is.finite(ess_bulk)]
  ess_tail <- conv$n_eff_tail %||% numeric(0)
  ess_tail <- ess_tail[is.finite(ess_tail)]

  list(
    engine = .fb_engine_line(engine),
    max_rhat = if (length(rhat) == 0L) NA_real_ else max(rhat),
    min_ess_bulk = if (length(ess_bulk) == 0L) NA_real_ else min(ess_bulk),
    min_ess_tail = if (length(ess_tail) == 0L) NA_real_ else min(ess_tail),
    n_divergent = conv$n_divergent %||% NA_integer_
  )
}

# .fb_summary_counts() --- the design and observation counts.
#
# One source: the record the missing-response layer left on the fit. A
# fit assembled without that layer (a direct emit_*() call in a test)
# falls back to the emit-time row count, where design and observed are
# the same number because nothing was augmented.
#
# @noRd
# @keywords internal
.fb_summary_counts <- function(object) {
  rec <- object$extras$na_action
  n_design <- rec$n_design %||% object$extras$model_info$n_obs %||% NA_integer_
  n_observed <- rec$n_observed %||% n_design
  list(
    n_design = as.integer(n_design),
    n_observed = as.integer(n_observed)
  )
}


# ---------------------------------------------------------------- #
# The constructor                                                   #
# ---------------------------------------------------------------- #

#' Build the unified summary object for a flexyBayes fit
#'
#' Assembles the eleven-slot `summary.flexybayes` object that every
#' engine's `summary()` method returns, plus one conditional twelfth slot
#' on a fit carrying an autoregressive latent field. Exposed as an
#' internal helper so the two methods share one construction and cannot
#' drift apart.
#'
#' The fit itself is attached as the attribute `fit`. The print method
#' renders the engine-native blocks -- INLA's hyperparameter table and
#' spatial-field panel, brms's per-level residual panel -- through the
#' same helpers `print()` uses, and reading them off the summary object
#' would mean maintaining a second copy of each.
#'
#' `spatial_field` is the exception, and it is present only where the
#' model has a field. The 0.9.0 INLA summary returned that table as a
#' slot, the unified object dropped it, and the removal was silent: a
#' caller reading `summary(fit)$spatial_field` got `NULL` and whatever
#' arithmetic followed. It is restored here as an engine-native optional
#' slot rather than folded into the frozen eleven, because a fit without
#' a field has no such table and an all-`NA` row would be worse than an
#' absent one.
#'
#' @param object A fitted `flexybayes` object of any active engine.
#' @returns An object of class `c("summary.flexybayes", "list")`.
#'
#' @noRd
#' @keywords internal
.fb_summary_object <- function(object) {
  counts <- .fb_summary_counts(object)
  varcomp <- .fb_summary_varcomp(object)

  out <- structure(
    list(
      fixed = .fb_summary_fixed(object),
      varcomp = varcomp,
      random = .fb_summary_random(object),
      missing = .fb_summary_missing(object),
      converge = .fb_summary_converge(object),
      n_design = counts$n_design,
      n_observed = counts$n_observed,
      na_action = object$extras$na_action,
      model = .fb_model_line_for_fit(object),
      engine = .fb_fit_engine(object),
      call = object$extras$the_call
    ),
    class = c("summary.flexybayes", "list")
  )

  # The same table the printed field panel renders, unrounded and under
  # its own column names. `.inla_spatial_hyper_table()` returns NULL on
  # any fit without a field, and on any engine that is not INLA, so the
  # slot appears exactly where 0.9.0 put it.
  field <- .inla_spatial_hyper_table(object)
  if (!is.null(field)) {
    out$spatial_field <- field
  }

  attr(out, "fit") <- object
  out
}


# ---------------------------------------------------------------- #
# Printing                                                          #
# ---------------------------------------------------------------- #

# .fb_print_header() --- the block both print() and summary() open with.
#
# @noRd
# @keywords internal
.fb_print_header <- function(object, title, rule) {
  ci <- object$extras$call_info
  mi <- object$extras$model_info
  engine <- .fb_fit_engine(object)
  counts <- .fb_summary_counts(object)

  # An aggregated fit names its representation where a per-row fit names
  # its engine: `aggregated-gaussian` is the token the package uses to
  # signal exactness, and the aggregated print has carried it since
  # v0.3.8. The two surfaces read the label off one helper.
  label <- if (inherits(object, "flexybayes_aggregated")) {
    .agg_fit_label(mi)
  } else {
    engine
  }
  cat(sprintf("%s  [flexyBayes / %s]\n", title, label))
  cat(strrep(rule, 62L), "\n")

  if (!is.null(ci$fixed)) {
    cat("  Fixed    :", deparse(ci$fixed), "\n")
  }
  if (!is.null(ci$random)) {
    cat("  Random   :", deparse(ci$random), "\n")
  }
  if (!is.null(ci$residual) && !identical(deparse(ci$residual), "~units")) {
    cat("  Residual :", deparse(ci$residual), "\n")
  }
  ir <- .fb_fit_ir(object)
  model <- .fb_model_line(ir)
  if (!is.na(model)) {
    cat("  Model    :", model, "\n")
    # The two models share a spelling and differ by a parameter. Saying
    # so on the line that describes the model is cheaper for the reader
    # than leaving it to be discovered in the variance components.
    if (.fb_ir_has_field(ir)) {
      cat(
        "             this is not ASReml residual = ~ ar1:ar1 ",
        "(3 parameters)\n",
        sep = ""
      )
    }
  }
  # The aggregated fit's `prior_parametrization` label is deliberately
  # NOT reproduced here. This object's prior surface is the `prior`
  # column of the variance-component table, projected from
  # prior_summary(), which states the prior each component actually
  # carried rather than classifying the call. The label stays on
  # print.flexybayes_aggregated(), where it has lived since v0.3.8.
  agg <- .fb_agg_line(object)
  if (!is.null(agg)) {
    cat("  aggregation: ", agg, "\n", sep = "")
  }
  if (!is.null(mi$family)) {
    cat("  Family   :", mi$family, "/", mi$link %||% "identity", "\n")
  }

  # The representation / engine pair is the package's standing truth
  # display (v0.3.8): the representation regime on one line, the
  # inference engine on the next. A fit that predates it, or one
  # assembled without a routing decision, still gets an engine line.
  if (!is.null(object$exactness)) {
    bd <- object$extras$backend_decision
    cat("  Representation: ", .repr_label_for_fit(object, bd), "\n", sep = "")
    cat("  Engine:         ", .engine_label_for_fit(object, bd), "\n", sep = "")
  } else {
    cat("  Engine   :", .fb_engine_line(engine), "\n")
  }

  if (!is.na(counts$n_design)) {
    cat(sprintf(
      "  N        : %d design rows, %d observed responses\n",
      counts$n_design,
      counts$n_observed
    ))
  }
  rec <- object$extras$na_action
  if (!is.null(rec)) {
    cat(sprintf(
      "  na_action: %s; %d missing response(s); %d cell(s) completed\n",
      rec$na_action,
      rec$n_missing_response %||% 0L,
      rec$n_cells_completed %||% 0L
    ))
  }

  # Sampler settings describe something that happened only on an engine
  # that sampled.
  if (.fb_fit_sampled(object)) {
    cat(sprintf(
      "  Sampler  : %s chain(s) x %s samples (warmup %s)\n",
      format(ci$chains),
      format(ci$n_samples),
      format(ci$warmup)
    ))
  }
  if (!is.null(object$extras$run_time)) {
    cat("  Runtime  :", round(object$extras$run_time, 1), "sec\n")
  }
  invisible(NULL)
}

# .fb_round_prior_cell() --- the prior string at display precision.
#
# The resolved prior carries a scale to twelve significant digits, which
# is the number the projection identity is checked on and roughly nine
# digits more than a table cell can be read at. The rounding happens
# here, on the way to the console, so the stored cell stays byte-equal to
# what prior_summary() reports and nothing downstream has to know that
# the printed table was abbreviated.
#
# @noRd
# @keywords internal
.fb_round_prior_cell <- function(x, digits = 4L) {
  if (!is.character(x) || length(x) == 0L) {
    return(x)
  }
  hits <- gregexpr(
    "(?<=[=(, ])[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?",
    x,
    perl = TRUE
  )
  regmatches(x, hits) <- lapply(regmatches(x, hits), function(v) {
    if (length(v) == 0L) {
      return(v)
    }
    num <- suppressWarnings(as.numeric(v))
    ok <- !is.na(num)
    # Formatted one at a time: format() on a vector pads every element to
    # a common width, which would turn a bound of 0 into "0.000".
    v[ok] <- vapply(
      num[ok],
      function(z) format(signif(z, digits), trim = TRUE, scientific = FALSE),
      character(1L)
    )
    v
  })
  x
}

# .fb_print_converge() --- the diagnostics block, in the engine's terms.
#
# @noRd
# @keywords internal
.fb_print_converge <- function(cv) {
  cat("\n-- Convergence ", strrep("-", 47), "\n", sep = "")
  if (is.null(cv)) {
    cat("  (none recorded)\n")
    return(invisible(NULL))
  }
  if (!is.null(cv$mode_status)) {
    cat("  Engine    : ", cv$engine, "\n", sep = "")
    cat(
      "  Mode status: ",
      format(cv$mode_status),
      " (0 = converged)\n",
      sep = ""
    )
    cat(
      "  Marginal log-likelihood: ",
      format(round(cv$mlik, 3)),
      "\n",
      sep = ""
    )
    if (!is.na(cv$kld_max)) {
      cat(
        "  Largest symmetric KLD between the Gaussian approximation and ",
        "the\n  corrected marginal: ",
        format(signif(cv$kld_max, 3)),
        "\n",
        sep = ""
      )
    }
    # No line at all where no check ran: a printed PASS/FAIL is a claim
    # about a confirmation step, and the aggregated route runs none.
    if (!is.na(cv$numerical_confirm)) {
      cat(
        "  Numerical confirm: ",
        if (isTRUE(cv$numerical_confirm)) "PASS" else "FAIL",
        if (length(cv$numerical_confirm_reasons) > 0L) {
          paste0(
            " (",
            paste(cv$numerical_confirm_reasons, collapse = "; "),
            ")"
          )
        } else {
          ""
        },
        "\n",
        sep = ""
      )
    }
    return(invisible(NULL))
  }
  cat("  Engine    : ", cv$engine, "\n", sep = "")
  if (!is.na(cv$max_rhat)) {
    cat("  Max R-hat : ", format(round(cv$max_rhat, 3)), "\n", sep = "")
  }
  if (!is.na(cv$min_ess_bulk)) {
    cat("  Min ESS (bulk): ", format(round(cv$min_ess_bulk)), "\n", sep = "")
  }
  if (!is.na(cv$min_ess_tail)) {
    cat("  Min ESS (tail): ", format(round(cv$min_ess_tail)), "\n", sep = "")
  }
  if (!is.na(cv$n_divergent)) {
    cat("  Divergent transitions: ", format(cv$n_divergent), "\n", sep = "")
  }
  invisible(NULL)
}

#' Print a flexyBayes summary
#'
#' Prints the fixed-effect table, the variance components with the prior
#' each one carried, the engine's own convergence diagnostics, and
#' whatever engine-native panels the fit earns -- INLA's hyperparameter
#' and spatial-field tables, brms's per-level residual table.
#'
#' Two sentences on the output are deliberate. The banner above the
#' variance components names the estimator, because the table looks like
#' an ASReml variance-component table and is not one: every number in it
#' is a posterior summary. On a fit carrying a latent field, the sentence
#' below the table names the fourth parameter, because a separable AR1
#' field plus an independent nugget is a different model from ASReml's
#' three-parameter nugget-free residual and the tables look alike.
#'
#' The numbers inside the `prior` cell are rounded for display only. The
#' stored cell carries the resolved prior's own string at full precision,
#' which is what [prior_summary()] reports and what the two are compared
#' on.
#'
#' @param x A `summary.flexybayes` object, as returned by [summary()] on
#'   a fitted model.
#' @param digits Number of significant digits for the printed tables,
#'   including the scales inside the `prior` cell.
#' @param ... Ignored, present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the tables it prints.
#' @export
print.summary.flexybayes <- function(x, digits = 4L, ...) {
  fit <- attr(x, "fit")
  if (!is.null(fit)) {
    .fb_print_header(fit, "Bayesian mixed model summary", "=")
  }

  cat(
    "\n-- Fixed effects (posterior mean, 95% credible interval) ",
    strrep("-", 8),
    "\n",
    sep = ""
  )
  if (nrow(x$fixed) > 0L) {
    body <- x$fixed
    num <- vapply(body, is.numeric, logical(1L))
    body[num] <- lapply(body[num], round, digits = digits)
    print(body, row.names = FALSE)
  } else {
    cat("  (none)\n")
  }

  cat(
    "\n-- Variance components (95% credible interval) ",
    strrep("-", 18),
    "\n",
    sep = ""
  )
  cat(
    "Estimate = posterior mean (not REML). ",
    "Prior column is the prior that was used.\n",
    sep = ""
  )
  if (nrow(x$varcomp) > 0L) {
    body <- x$varcomp
    num <- vapply(body, is.numeric, logical(1L))
    body[num] <- lapply(body[num], round, digits = digits)
    body$prior[is.na(body$prior)] <- ""
    body$prior <- .fb_round_prior_cell(body$prior, digits = digits)
    print(body, row.names = FALSE)
    if (.fb_is_spatial_varcomp(x$varcomp)) {
      cat(
        "The field and the nugget are separate parameters: this is not ",
        "ASReml's nugget-free residual.\n",
        sep = ""
      )
    }
  } else {
    cat("  (none available)\n")
  }

  # Engine-native panels, rendered by the helpers print() uses so the two
  # surfaces cannot drift.
  if (!is.null(fit)) {
    .print_inla_hyperpar_table(fit)
    .print_inla_spatial_hypers(fit)
    .print_brms_residual_by_level(fit)
  }

  .fb_print_converge(x$converge)
  invisible(x)
}
