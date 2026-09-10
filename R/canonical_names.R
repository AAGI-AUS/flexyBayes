# canonical_names -- per-backend parameter-name resolver for flexyBayes.
#
# The `canonical_names()` S3 generic returns a list `(map, transform,
# source)` keyed by backend-native parameter names; values give the
# canonical (brms-style) parameter names; the optional transform
# slot carries value-transform functions (used by INLA's
# precision-to-SD `sqrt(1/prec)` mapping).
#
# Per-backend mapper functions register at package load via the
# internal list `.canonical_mappers`. Each mapper accepts (fit,
# fb_terms) and returns `list(map, transform)`. Users can register
# additional mappers via `register_canonical_mapper()` (internal in
# v0.2; exported in v0.3+ when third-party backend authors might use
# it).
#
# The INLA-side mapping uses summary.fixed rownames for fixed effects
# (already brms-canonical) and translates the "Precision for ..."
# hyperparameter names plus the precision-to-SD value transform. The
# brms-side mapping delegates to Stan's own canonical naming.

# ---------------------------------------------------------------- #
# Internal registry                                                  #
# ---------------------------------------------------------------- #
# Populated at the bottom of this file with the active-engine mappers
# (inla + brms). Additional mappers (e.g., for a future NIMBLE backend)
# register via `register_canonical_mapper(backend, mapper)`.

.canonical_mappers <- new.env(parent = emptyenv())


# ---------------------------------------------------------------- #
# Internal helper -- coalescing operator                            #
# ---------------------------------------------------------------- #

`%||%` <- function(x, y) if (is.null(x)) y else x


# ---------------------------------------------------------------- #
# Generic + methods                                                  #
# ---------------------------------------------------------------- #

#' Canonical parameter-name view for a flexyBayes fit
#'
#' Returns the backend-native -> canonical parameter-name map for a
#' fit, plus per-parameter value transforms where applicable (e.g.,
#' the INLA precision-to-SD `sqrt(1/prec)` transform applied to
#' hyperparameters before triangulation). The canonical convention
#' follows brms (`(Intercept)`, `<term>`, `sd_<group>`, `sigma`,
#' `r_<group>[<level>]`).
#'
#' On `flexybayes_inla` and `flexybayes_brms` fits, the per-backend
#' mapper registered at package load (`inla` or `brms`) drives the
#' resolution; the returned map is cached on
#' `fit$extras$canonical_map` for fast repeated access.
#'
#' @param fit A `flexybayes_inla` or `flexybayes_brms` object.
#' @param drop Logical: if `TRUE` (default `FALSE`), drop
#'   backend-native names that are not in the registered map
#'   (e.g., INLA's `Predictor.<i>` latent-predictor draws). When
#'   `FALSE`, un-mapped names appear in the returned `$unmapped`
#'   element.
#' @param ... Additional arguments (ignored by current methods).
#' @returns A list holding the canonical-name map and its provenance.
#'
#'   \describe{
#'     \item{`map`}{Named character vector keyed by backend-native
#'       parameter name with canonical name as the value.}
#'     \item{`transform`}{Named list of `function(x) -> x'`
#'       transforms keyed by canonical name. Empty list when no
#'       transforms apply.}
#'     \item{`source`}{Character: `"registry"`, `"user"`,
#'       `"registry_fallback_verbatim"`, or `"legacy_inferred"`.}
#'     \item{`unmapped`}{Character vector of backend-native names
#'       not in the map (when `drop = FALSE`).}
#'     \item{`prior_parametrization`}{Character, present only on
#'       aggregated fits, naming which prior the fit ran under.
#'       `"per_row_equivalent"` is the legacy scalar bridge, whose
#'       precision prior makes the aggregated posterior match the
#'       per-row posterior to numerical precision.
#'       `"package_default"` is the automatic bounded-uniform-on-SD
#'       prior, which claims no such equivalence on this route.
#'       `"custom"` is an explicit prior from the caller (see the
#'       "Matched priors" note on [triangulate()]).}
#'   }
#' @examples
#' # The map is engine-specific, so it needs a fit to read.
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(y = rnorm(60), x = rnorm(60), g = factor(rep(1:6, 10)))
#'   fit <- flexybayes(y ~ x + (1 | g), data = d, backend = "inla",
#'                     verbose = FALSE)
#'   cn <- canonical_names(fit)
#'   head(cn$map)
#' }
#' }
#' @export
canonical_names <- function(fit, drop = FALSE, ...) {
  UseMethod("canonical_names")
}

#' @rdname canonical_names
#' @export
canonical_names.flexybayes_inla <- function(fit, drop = FALSE, ...) {
  if (!is.null(fit$extras$canonical_map)) {
    return(.maybe_attach_unmapped(
      fit$extras$canonical_map,
      fit,
      drop,
      backend = "inla"
    ))
  }
  fb <- fit$extras$fb_terms %||% fit$fb
  mapper <- .canonical_mappers$inla
  if (is.null(mapper)) {
    return(list(map = character(0), transform = list(), source = "no_mapper"))
  }
  res <- mapper(fit, fb)
  res$source <- res$source %||% "registry"
  .maybe_attach_unmapped(res, fit, drop, backend = "inla")
}

#' @rdname canonical_names
#' @export
canonical_names.flexybayes_brms <- function(fit, drop = FALSE, ...) {
  if (!is.null(fit$extras$canonical_map)) {
    return(.maybe_attach_unmapped(
      fit$extras$canonical_map,
      fit,
      drop,
      backend = "brms"
    ))
  }
  fb <- fit$extras$fb_terms %||% fit$extras$parse_info
  mapper <- .canonical_mappers$brms
  if (is.null(mapper)) {
    return(list(map = character(0), transform = list(), source = "no_mapper"))
  }
  res <- mapper(fit, fb)
  res$source <- res$source %||% "registry"
  .maybe_attach_unmapped(res, fit, drop, backend = "brms")
}


# ---------------------------------------------------------------- #
# Registration helper                                               #
# ---------------------------------------------------------------- #

# Internal in v0.2 -- exported in v0.3 once the third-party backend
# author surface is documented. The helper attaches a mapper to the
# .canonical_mappers registry under the given backend key.
register_canonical_mapper <- function(backend, mapper) {
  if (!is.character(backend) || length(backend) != 1L || !nzchar(backend)) {
    stop(
      "`backend` must be a non-empty length-1 character string.",
      call. = FALSE
    )
  }
  if (!is.function(mapper)) {
    stop(
      "`mapper` must be a function (fit, fb_terms) -> ",
      "list(map, transform).",
      call. = FALSE
    )
  }
  assign(backend, mapper, envir = .canonical_mappers)
  invisible(NULL)
}


# ---------------------------------------------------------------- #
# Helpers                                                            #
# ---------------------------------------------------------------- #

# Attach an `$unmapped` element to the canonical-names result when
# requested. The `unmapped` set is the difference between (the
# names present in the fit's posterior draws) and (the keys of the
# canonical map). This surfaces backend-side draws (e.g., INLA's
# Predictor.* linear-predictor samples) that the registry doesn't
# translate.
.maybe_attach_unmapped <- function(res, fit, drop, backend) {
  # Surface the aggregated-fit prior parametrization on the result so a
  # consumer can tell whether the matched-prior (per-row-equivalent)
  # guarantee holds without reaching into the fit internals. Absent on
  # non-aggregated fits (the meta slot is NULL), so the field simply does
  # not appear there. Attached on both the drop = TRUE and drop = FALSE
  # paths since it is independent of the unmapped-draws logic.
  pp <- fit$extras$aggregation_meta$prior_parametrization
  if (!is.null(pp)) {
    res$prior_parametrization <- pp
  }
  if (isTRUE(drop)) {
    return(res)
  }
  draws_names <- tryCatch(
    .backend_draw_names(fit, backend),
    error = function(e) character(0)
  )
  if (length(draws_names) == 0L) {
    res$unmapped <- character(0)
    return(res)
  }
  res$unmapped <- setdiff(draws_names, names(res$map))
  res
}

# Backend-native draw names for the unmapped set. INLA fits expose them
# via rownames(summary.fixed) + rownames(summary.hyperpar); brms fits via
# posterior::as_draws_matrix().
.backend_draw_names <- function(fit, backend) {
  if (identical(backend, "inla")) {
    if (is.null(fit$inla)) {
      return(character(0))
    }
    c(
      rownames(fit$inla$summary.fixed) %||% character(0),
      rownames(fit$inla$summary.hyperpar) %||% character(0)
    )
  } else if (identical(backend, "brms")) {
    if (is.null(fit$brms)) {
      return(character(0))
    }
    if (!requireNamespace("posterior", quietly = TRUE)) {
      return(character(0))
    }
    m <- tryCatch(
      as.matrix(posterior::as_draws_matrix(fit$brms)),
      error = function(e) NULL
    )
    if (is.null(m)) character(0) else colnames(m) %||% character(0)
  } else {
    character(0)
  }
}

# Parse an asreml-emit factor-coefficient name like "tau_env[3,1]"
# into its components: tag = "env", level_index = 3L. Returns NULL
# when the pattern does not match.
.parse_tau_factor <- function(nm) {
  m <- regmatches(nm, regexec("^tau_([^\\[]+)\\[(\\d+),(\\d+)\\]$", nm))[[1]]
  if (length(m) < 4L) {
    return(NULL)
  }
  list(tag = m[2], level = as.integer(m[3]), col = as.integer(m[4]))
}

# Parse an indexed-slope deviation name like
# "slope_dev_f_x_raw[2,1]" against the slope_dev_lookup keyed by
# `<factor>_<continuous>`. Returns NULL when the name does not look
# like an indexed-slope draw, or when no IR term matches the parsed
# tag. The raw-vector index `i` maps to factor level `i + 1` because
# the reference level (level 1) is fixed at zero and absent from the
# raw parameter vector.
# .fb_inla_group_labels() --- generated index name -> the group a user wrote.
#
# emit_inla() builds an index column per structured random term --
# `rep_block_id` for `rep:block`, `gen_at_env_E1` for one level of a diag(),
# `gen_by_env_id` for an unstructured covariance -- and INLA names each
# summary.random entry after the column it was indexed by. The variance
# COMPONENTS were mapped back to user labels when that machinery was built;
# the random-effect table was not, so `ranef()` reported a grouping factor
# called `rep_block_id`, which is a name the user never wrote and cannot look
# up. Same defect class as the one closed for summary()$varcomp, reached by a
# different accessor.
#
# Returns a named character vector; a group with no entry keeps its own name,
# which is right for a plain `~ g` where the index IS the variable.
.fb_inla_group_labels <- function(fb) {
  map <- character(0)
  put <- function(key, value) map[[key]] <<- value
  for (t in fb$random_terms %||% list()) {
    tt <- t$type %||% ""
    if (tt %in% c("nested", "combo")) {
      vars <- if (identical(tt, "nested")) {
        c(t$outer, t$inner)
      } else {
        as.character(t$vars)
      }
      put(
        make.names(paste0(paste(vars, collapse = "_"), "_id")),
        paste(vars, collapse = ":")
      )
      next
    }
    if (identical(tt, "at_simple")) {
      for (lev in t$outer_levels %||% character(0)) {
        put(
          make.names(paste0(t$inner, "_at_", t$outer, "_", lev)),
          paste0(t$inner, " (", t$outer, " = ", lev, ")")
        )
      }
      next
    }
    if (identical(tt, "simple_slope_uncor")) {
      put(
        make.names(paste0(t$var, "_for_slope_", t$slope_var)),
        paste0(t$slope_var, " | ", t$var)
      )
      next
    }
    if (identical(tt, "simple_slope_cor")) {
      lab <- paste0(t$slope_var, " | ", t$var)
      put(make.names(paste0(t$var, "_cor_", t$slope_var, "_id")), lab)
      put(make.names(paste0(t$var, "_cor_", t$slope_var, "_slope_id")), lab)
      next
    }
    if (identical(tt, "us_gxe")) {
      put(
        make.names(paste0(t$inner, "_by_", t$outer, "_id")),
        paste0(t$inner, ":", t$outer)
      )
      next
    }
    if (identical(tt, "fa_gxe")) {
      levs <- t$outer_levels %||% character(0)
      for (j in seq_len(t$k %||% 1L)) {
        for (lev in levs) {
          put(
            make.names(paste0(t$inner, "_fa", j, "_", t$outer, "_", lev)),
            paste0(t$inner, " factor ", j, " (", t$outer, " = ", lev, ")")
          )
        }
      }
      for (lev in levs) {
        put(
          make.names(paste0(t$inner, "_psi_", t$outer, "_", lev)),
          paste0(t$inner, " specific (", t$outer, " = ", lev, ")")
        )
      }
      next
    }
  }
  map
}

.parse_slope_dev_raw <- function(nm, lookup) {
  m <- regmatches(
    nm,
    regexec("^slope_dev_(.+)_raw\\[(\\d+),(\\d+)\\]$", nm)
  )[[1]]
  if (length(m) < 4L) {
    return(NULL)
  }
  tag <- m[2]
  raw_idx <- as.integer(m[3])
  if (is.null(lookup[[tag]])) {
    return(NULL)
  }
  lvls <- lookup[[tag]]
  # raw index i -> factor level i + 1 (treatment-coded; reference
  # level pinned at zero and absent from the raw vector).
  level_idx <- raw_idx + 1L
  level_label <- if (level_idx >= 1L && level_idx <= length(lvls)) {
    lvls[level_idx]
  } else {
    as.character(level_idx)
  }
  list(
    canonical = paste0("slope_", tag, "[", level_label, "]"),
    tag = tag,
    raw_idx = raw_idx,
    level_label = level_label
  )
}


# ---------------------------------------------------------------- #
# v0.2 mapper: INLA backend                                         #
# ---------------------------------------------------------------- #
# INLA-side names exposed via fit$inla$summary.{fixed,hyperpar}:
#   summary.fixed rownames    -- already brms-canonical:
#       "(Intercept)", "<term>" (treatment-coded factor levels)
#   summary.hyperpar rownames:
#       "Precision for the Gaussian observations" -> "sigma" + sqrt(1/prec)
#       "Precision for <group>"                   -> "sd_<group>" + sqrt(1/prec)

.mapper_inla <- function(fit, fb_terms) {
  if (is.null(fit$inla)) {
    return(list(map = character(0), transform = list()))
  }

  map <- character(0)
  transform <- list()

  # Fixed effects -- identity mapping (already brms-canonical). The
  # summary.fixed table uses the bare name ("(Intercept)") while the
  # inla.posterior.sample() latent rows use the ":1"-suffixed form
  # ("(Intercept):1"). Register both so the rename fires regardless
  # of which entry point feeds the consumer.
  fixed_names <- rownames(fit$inla$summary.fixed) %||% character(0)
  for (nm in fixed_names) {
    map[nm] <- nm
    map[paste0(nm, ":1")] <- nm
  }

  # Hyperparameters -- precision-side names; rename + sqrt(1/prec). The
  # transform is keyed by the *native* hyperpar name (not the canonical
  # target), because triangulate() applies transforms to each fit's
  # draws while they still carry their backend-native names, before any
  # renaming (see triangulate.R "Apply per-parameter transforms before
  # any renaming"). Keying by the canonical name silently no-ops, which
  # would leave INLA precision draws un-converted and make a cross-engine
  # comparison pit precision against standard deviation.
  # Hyperpar names do not carry the ":1" suffix.
  prec_to_sd <- function(prec) sqrt(1 / prec)
  # WP16: AR1 / separable AR1xAR1 spatial latent-field hyperparameters. INLA
  # names them "Rho for <idx>", "GroupRho for <idx>", "Precision for <idx>"
  # where <idx> is the f() index column (<row_var>_id). Map to canonical
  # rho_row / rho_col / sd_spatial; the two Rho values are correlations
  # (identity, no transform), the Precision is the field precision -> SD.
  # Matched BEFORE the generic "Precision for <group>" rule so the field
  # precision resolves to sd_spatial, not sd_<idx>.
  spatial_idx <- character(0)
  for (t in c(
    fb_terms$random_terms %||% list(),
    fb_terms$residual_terms %||% list()
  )) {
    tt <- t$type %||% ""
    if (identical(tt, "ar1_spatial")) {
      spatial_idx <- c(spatial_idx, paste0(t$row_var, "_id"))
    } else if (identical(tt, "ar1")) {
      spatial_idx <- c(spatial_idx, paste0(t$var, "_id"))
    }
  }
  # Nested / interaction random effects. The INLA emit fits these as an iid
  # field over a generated interaction index (<a>_<b>_id), so INLA names the
  # hyperparameter after that index. The generated name must not reach the
  # user's summary: map it back to the colon-joined group, which is also what
  # brms calls it, so the two engines line up on one canonical name.
  # The section levels of a dsum() residual, for naming its per-section
  # sigma rows below.
  resid_levels <- NULL
  for (t in fb_terms$residual_terms %||% list()) {
    if (identical(t$type %||% "", "at_units") && is.null(t$level)) {
      resid_levels <- t$var_levels
    }
  }

  nested_idx <- list()
  us_prec <- list()
  us_cor <- list()
  fa_load <- list()
  for (t in fb_terms$random_terms %||% list()) {
    tt <- t$type %||% ""
    # Uncorrelated random slope: the emit writes the slope half as a second
    # iid field over a generated index (<g>_for_slope_<x>). Map it back to
    # the package's own name for that component, <slope>_<group>, which is
    # also the key its prior is registered under.
    # Heterogeneous variance by factor level: one generated index per level
    # of the outer factor. Map each back to brms's name for the same
    # component, sd_<inner>__<outer><level>, so diag() reports identically
    # on both engines.
    # us(f):g -- one iidKd field. INLA names its parts "Precision for <idx>
    # (component j)" and "RhoI:J for <idx>". Map to what brms calls the same
    # quantities for `(0 + f | g)`: sd_<g>__<f><level> and
    # cor_<g>__<f><levelI>__<f><levelJ>.
    # fa(f, k):g -- three kinds of hyperparameter. The anchor field's
    # precision is the common factor's SD; each copy's Beta is a LOADING (a
    # scaling, no transform); each psi field's precision is that level's
    # specific SD, which is the same quantity diag(f):g reports, so it takes
    # the same name.
    if (identical(tt, "fa_gxe")) {
      levs <- t$outer_levels %||% character(0)
      for (j in seq_len(t$k %||% 0L)) {
        anchor <- make.names(paste0(t$inner, "_fa", j, "_", t$outer, "_",
                                    levs[[j]]))
        us_prec[[paste0("Precision for ", anchor)]] <- paste0("sd_fa", j)
        for (e in seq_along(levs)) {
          idx <- make.names(paste0(t$inner, "_fa", j, "_", t$outer, "_",
                                   levs[[e]]))
          fa_load[[paste0("Beta for ", idx)]] <-
            paste0("fa", j, "_", t$outer, levs[[e]])
        }
      }
      for (e in seq_along(levs)) {
        psi <- make.names(paste0(t$inner, "_psi_", t$outer, "_", levs[[e]]))
        us_prec[[paste0("Precision for ", psi)]] <-
          paste0("sd_", t$inner, "__", t$outer, levs[[e]])
      }
      next
    }
    if (identical(tt, "us_gxe")) {
      idx <- make.names(paste0(t$inner, "_by_", t$outer, "_id"))
      levs <- t$outer_levels %||% character(0)
      for (j in seq_along(levs)) {
        us_prec[[paste0("Precision for ", idx, " (component ", j, ")")]] <-
          paste0("sd_", t$inner, "__", t$outer, levs[[j]])
      }
      for (a in seq_along(levs)) {
        for (b in seq_along(levs)) {
          if (b <= a) next
          us_cor[[paste0("Rho", a, ":", b, " for ", idx)]] <-
            paste0("cor_", t$inner, "__", t$outer, levs[[a]],
                   "__", t$outer, levs[[b]])
        }
      }
      next
    }
    if (identical(tt, "at_simple")) {
      for (lev in t$outer_levels %||% character(0)) {
        nested_idx[[make.names(paste0(t$inner, "_at_", t$outer, "_", lev))]] <-
          paste0(t$inner, "__", t$outer, lev)
      }
      next
    }
    if (identical(tt, "simple_slope_uncor")) {
      nested_idx[[paste0(t$var, "_for_slope_", t$slope_var)]] <-
        paste0(t$slope_var, "_", t$var)
      next
    }
    # (x | g) is one iid2d field over a two-block index, so its
    # hyperparameters arrive in the same shape us() uses: two component
    # precisions and one Rho. Component 1 is the intercept block and
    # component 2 the slope block, which fixes the reading of both. The
    # standard deviations take the names the UNCORRELATED sibling already
    # reports, so the two models are read side by side; the correlation
    # takes brms's own cor_<g>__Intercept__<x>.
    if (identical(tt, "simple_slope_cor")) {
      idx <- make.names(paste0(t$var, "_cor_", t$slope_var, "_id"))
      us_prec[[paste0("Precision for ", idx, " (component 1)")]] <-
        paste0("sd_", t$var)
      us_prec[[paste0("Precision for ", idx, " (component 2)")]] <-
        paste0("sd_", t$slope_var, "_", t$var)
      us_cor[[paste0("Rho1:2 for ", idx)]] <-
        paste0("cor_", t$var, "__Intercept__", t$slope_var)
      next
    }
    if (tt %in% c("nested", "combo")) {
      vars <- if (identical(tt, "nested")) {
        c(t$outer, t$inner)
      } else {
        as.character(t$vars)
      }
      nested_idx[[paste0(paste(vars, collapse = "_"), "_id")]] <-
        paste(vars, collapse = ":")
    }
  }

  hyper_names <- rownames(fit$inla$summary.hyperpar) %||% character(0)
  for (nm in hyper_names) {
    sp_hit <- FALSE
    if (!is.null(us_prec[[nm]])) {
      map[nm] <- us_prec[[nm]]
      transform[[nm]] <- prec_to_sd
      next
    }
    if (!is.null(fa_load[[nm]])) {
      # A loading, not a precision: no transform.
      map[nm] <- fa_load[[nm]]
      next
    }
    if (!is.null(us_cor[[nm]])) {
      # A correlation, not a precision: no transform.
      map[nm] <- us_cor[[nm]]
      next
    }
    if (startsWith(nm, "Precision for ")) {
      idx <- sub("^Precision for ", "", nm)
      if (!is.null(nested_idx[[idx]])) {
        map[nm] <- paste0("sd_", nested_idx[[idx]])
        transform[[nm]] <- prec_to_sd
        next
      }
    }
    for (idx in spatial_idx) {
      if (identical(nm, paste0("Rho for ", idx))) {
        map[nm] <- "rho_row"
        sp_hit <- TRUE
        break
      }
      if (identical(nm, paste0("GroupRho for ", idx))) {
        map[nm] <- "rho_col"
        sp_hit <- TRUE
        break
      }
      if (identical(nm, paste0("Precision for ", idx))) {
        map[nm] <- "sd_spatial"
        transform[[nm]] <- prec_to_sd
        sp_hit <- TRUE
        break
      }
    }
    if (sp_hit) {
      next
    }
    # A sectioned residual gives INLA one Gaussian likelihood per section.
    # It names the first "Precision for the Gaussian observations" and the
    # rest "...[2]", "...[3]". Map them to sigma_<level>, the same
    # vocabulary the brms path already uses for dsum(), so the two engines
    # report the sectioned residual under one set of names.
    if (length(resid_levels) > 1L &&
          startsWith(nm, "Precision for the Gaussian observations")) {
      suffix <- sub("^Precision for the Gaussian observations", "", nm)
      j <- if (nzchar(suffix)) {
        as.integer(gsub("[^0-9]", "", suffix))
      } else {
        1L
      }
      if (!is.na(j) && j >= 1L && j <= length(resid_levels)) {
        map[nm] <- paste0("sigma_", resid_levels[[j]])
        transform[[nm]] <- prec_to_sd
        next
      }
    }
    if (identical(nm, "Precision for the Gaussian observations")) {
      map[nm] <- "sigma"
      transform[[nm]] <- prec_to_sd
      next
    }
    if (startsWith(nm, "Precision for ")) {
      grp <- sub("^Precision for ", "", nm)
      map[nm] <- paste0("sd_", grp)
      transform[[nm]] <- prec_to_sd
      next
    }
    # Other hyperparameters (e.g., rho-for-ar1) left un-mapped for
    # v0.2; future ADR will register them when the structured-cov
    # work lands.
  }

  list(map = map, transform = transform)
}


# ---------------------------------------------------------------- #
# v0.2 mapper: brms-via-Stan backend (Stan passthrough)             #
# ---------------------------------------------------------------- #
# brms parameter names are the canonical names by definition (the
# brms convention IS the flexyBayes canonical convention). The
# mapper is identity over the draw-matrix column names with one
# prefix-strip rule:
#
#   b_<term>           -> <term>             (e.g. b_Days -> Days)
#   b_Intercept        -> "(Intercept)"
#   sd_<group>__Intercept -> sd_<group>      (drop the __Intercept tail)
#   sigma              -> sigma              (identity)
#
# Hyperparameters that brms exposes but flexyBayes does not yet
# canonicalise (e.g. `lp__` log-posterior, `Intercept` -- the
# scaled intercept brms emits alongside b_Intercept) are left
# un-mapped and surface via `$unmapped`.

.mapper_brms_via_stan <- function(fit, fb_terms) {
  if (is.null(fit$brms)) {
    return(list(map = character(0), transform = list()))
  }

  if (!requireNamespace("posterior", quietly = TRUE)) {
    return(list(map = character(0), transform = list()))
  }

  draw_names <- tryCatch(
    colnames(as.matrix(posterior::as_draws_matrix(fit$brms))),
    error = function(e) character(0)
  )
  if (length(draw_names) == 0L) {
    return(list(map = character(0), transform = list()))
  }

  # Build a {factor, continuous, levels} lookup for every
  # factor_numeric_interaction term in the IR so we can map brms's
  # `b_<fac><level>:<con>` (or `b_<con>:<fac><level>`) parameters
  # to the canonical `slope_<fac>_<con>[<level>]` slot.
  fixed_terms_fb <- if (!is.null(fb_terms)) {
    fb_terms$fixed_terms %||% list()
  } else {
    list()
  }
  fni_terms <- Filter(
    function(t) identical(t$type, "factor_numeric_interaction"),
    fixed_terms_fb
  )

  map <- character(0)
  for (nm in draw_names) {
    if (identical(nm, "b_Intercept")) {
      map[nm] <- "(Intercept)"
      next
    }
    # Factor:continuous indexed-slope translation.
    # brms emits these as `b_<term1>:<term2>` where `<term1>` is
    # `<fac><level>` for the non-reference level. Try both
    # (factor-first, continuous-first) orientations.
    fni_canonical <- .brms_factor_numeric_interaction(nm, fni_terms)
    if (!is.null(fni_canonical)) {
      map[nm] <- fni_canonical
      next
    }
    if (startsWith(nm, "b_")) {
      bare <- sub("^b_", "", nm)
      map[nm] <- bare
      next
    }
    if (startsWith(nm, "sd_")) {
      bare <- sub("^sd_", "", nm)
      # brms emits sd_<g>__<x> for uncorrelated random slopes
      # ((x || g) / (1 + x || g)) and sd_<g>__Intercept for the
      # corresponding intercept-deviation. Canonical name is
      # sd_<x>_<g> for slopes; sd_<g> for the intercept. The
      # __Intercept rule comes first because it is a subset of the
      # __<x> rule.
      if (grepl("__Intercept$", bare)) {
        bare <- sub("__Intercept$", "", bare)
        if (nzchar(bare)) {
          map[nm] <- paste0("sd_", bare)
        }
        next
      }
      m <- regmatches(bare, regexec("^([^_]+(?:_[^_]+)*?)__(.+)$", bare))[[1]]
      if (length(m) == 3L && nzchar(m[2L]) && nzchar(m[3L])) {
        # bare matched <g>__<x>; canonicalise to sd_<x>_<g>.
        g_part <- m[2L]
        x_part <- m[3L]
        map[nm] <- paste0("sd_", x_part, "_", g_part)
        next
      }
      if (nzchar(bare)) {
        map[nm] <- paste0("sd_", bare)
      }
      next
    }
    if (identical(nm, "sigma")) {
      map[nm] <- "sigma"
      next
    }
    # Other names (lp__, Intercept-centered, r_<group>[<lvl>,Intercept],
    # cor_<group>__*, etc.) left un-mapped for v0.2; surfaced via
    # $unmapped when drop = FALSE.
  }

  list(map = map, transform = list())
}

# Helper. Translate a brms parameter name like
# `b_fb:x` or `b_x:fb` (treatment-coded factor:continuous interaction)
# to the canonical `slope_<fac>_<con>[<level>]` slot. Returns NULL
# when the name does not match any factor_numeric_interaction term
# in the IR. Reference-level rows (level == lvls[1]) are absent from
# brms's draws by design (treatment coding); they never match.
.brms_factor_numeric_interaction <- function(nm, fni_terms) {
  if (length(fni_terms) == 0L) {
    return(NULL)
  }
  if (!startsWith(nm, "b_")) {
    return(NULL)
  }
  body <- sub("^b_", "", nm)
  if (!grepl(":", body, fixed = TRUE)) {
    return(NULL)
  }
  parts <- strsplit(body, ":", fixed = TRUE)[[1L]]
  if (length(parts) != 2L) {
    return(NULL)
  }
  for (t in fni_terms) {
    fac <- t$factor
    con <- t$continuous
    lvls <- t$levels %||% NULL
    if (is.null(fac) || is.null(con) || is.null(lvls)) {
      next
    }
    # Try factor-first, then continuous-first.
    for (orient in list(c(1L, 2L), c(2L, 1L))) {
      fac_part <- parts[[orient[[1L]]]]
      con_part <- parts[[orient[[2L]]]]
      if (!identical(con_part, con)) {
        next
      }
      if (!startsWith(fac_part, fac)) {
        next
      }
      lvl <- substr(fac_part, nchar(fac) + 1L, nchar(fac_part))
      if (!nzchar(lvl)) {
        next
      }
      if (lvl %in% lvls && !identical(lvl, lvls[[1L]])) {
        return(paste0("slope_", fac, "_", con, "[", lvl, "]"))
      }
    }
  }
  NULL
}


# ---------------------------------------------------------------- #
# Package-load registration                                         #
# ---------------------------------------------------------------- #
# Run at file-source time. R's package-build process sources the R/
# files in alphabetical order; canonical_names.R sources before
# triangulate.R (the consumer), so the registry is populated by the
# time triangulate() runs.

register_canonical_mapper("inla", .mapper_inla)
register_canonical_mapper("brms", .mapper_brms_via_stan)
