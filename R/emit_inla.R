# emit_inla -- flexyBayes INLA backend emit.
#
# v0.1 MINIMUM SUBSET: takes an fb_terms object (the IR --
# intermediate representation; see R/fb_terms.R) that has passed
# lgm_gate() and runs INLA -- Integrated Nested Laplace
# Approximations -- via INLA::inla() to produce a posterior fit.
#
# v0.1 supported model classes:
#   - Gaussian / binomial / binary / poisson / Gamma / beta /
#     lognormal / nbinomial families.
#   - Linear fixed effects (factor / continuous / interaction /
#     factor_interaction / I() expressions).
#   - Random intercepts via term type "simple" / "ide" / "id"
#     emitted as `f(group, model = "iid")`.
#   - Random walk / spline smoothers (term type "spline") emitted
#     as `f(var, model = "rw2")`.
#   - Heterogeneous residual via `at_units` -- v0.1 refused;
#     re-route via lgm_gate or backend = "brms".
#   - Spatial / GxE structured covariances -- v0.1 refused; route
#     via brms.
#
# Anything not in the supported set raises an INLA-side refusal
# pointing the user back to backend = "brms".
#
# The post-fit numerical-confirm gate lives here: assert
# mode.status == 0 and a finite marginal likelihood (mlik). Failures
# escalate to a structured warning.
#
# INLA is not on CRAN; ships via inla.r-inla-download.org. Added to
# DESCRIPTION:Suggests with Additional_repositories declaration.

# ---------------------------------------------------------------- #
# S2 -- no engine death without a typed refusal (C2, FS-25)         #
# ---------------------------------------------------------------- #

# .inla_call() --- thin, mockable wrapper around the one INLA engine
# call site. Exists so (a) there is exactly one place the death-
# pattern classification below has to route through, and (b) tests
# can substitute a controlled failure via
# testthat::local_mocked_bindings(.inla_call = ...) without needing a
# real INLA program crash (slow, and not reliably reproducible on
# demand). Carries no error handling of its own -- the caller's
# tryCatch is what turns a raw engine death into a typed refusal, so
# mocking this function alone exercises that translation exactly as a
# real crash would.
.inla_call <- function(args) {
  do.call(INLA::inla, args)
}

# .inla_is_program_death_message() --- TRUE when an error message
# names the INLA program dying, rather than a flexyBayes-side argument
# or formula error. Deliberately narrow: FS-25's observed failure text
# is "The inla-program exited with an error. Unless you interrupted it
# yourself, ..."; the alternate spelling `inla.core.safe` is INLA's
# other reporting path for the same class of event. Only these
# patterns route to `inla_program_failed` -- any other message is
# flexyBayes's own code and must propagate unwrapped (own-code errors
# must not be misnamed as an engine death).
.inla_is_program_death_message <- function(msg) {
  if (!is.character(msg) || length(msg) != 1L || is.na(msg)) {
    return(FALSE)
  }
  grepl("inla-program|inla\\.core\\.safe|\\bexited\\b", msg, ignore.case = TRUE)
}

# .inla_design_latent_effect_summary() --- best-effort tally of the
# failed design's random-effect "nodes", read directly from the random
# terms' grouping variable(s) against the data already in scope at the
# emit call. This is NOT a re-invocation of .fb_preflight()'s byte-cost
# model: that object is computed upstream in .dispatch_backend() (only
# above the 1e5-row threshold) and is not threaded down to emit_inla()
# through the backend-registry dispatch table, so it is not available
# here without new cross-file plumbing (recorded as a Handoff). This
# tally is a plain sum of observed level counts per random term,
# best-effort -- a term this function cannot size (no resolvable
# grouping variable in `data`) contributes nothing rather than aborting
# the refusal it is building.
#
# Returns list(total = <double or NA>, binding_term = <character or
# NA>, binding_n = <double or NA>) -- `binding_term` names the single
# random term contributing the most nodes, mirroring the preflight's
# own "binding term" vocabulary (R/fb_preflight.R
# design_memory_exceeds_ceiling).
.inla_design_latent_effect_summary <- function(fb, data) {
  empty <- list(
    total = NA_real_, binding_term = NA_character_, binding_n = NA_real_
  )
  terms <- fb$random_terms
  if (!length(terms) || !is.data.frame(data)) {
    return(empty)
  }
  total <- 0
  counted_any <- FALSE
  binding_term <- NA_character_
  binding_n <- -Inf
  for (t in terms) {
    vars <- if (!is.null(t$vars)) {
      as.character(t$vars)
    } else if (!is.null(t$var)) {
      as.character(t$var)
    } else {
      character(0L)
    }
    vars <- vars[nzchar(vars) & vars %in% names(data)]
    if (!length(vars)) {
      next
    }
    n_levels <- tryCatch(
      as.numeric(nrow(unique(data[vars]))),
      error = function(e) NA_real_
    )
    if (is.na(n_levels)) {
      next
    }
    counted_any <- TRUE
    total <- total + n_levels
    if (n_levels > binding_n) {
      binding_n <- n_levels
      binding_term <- if (!is.null(t$label) && nzchar(t$label)) {
        t$label
      } else {
        paste(vars, collapse = ":")
      }
    }
  }
  if (!counted_any) {
    return(empty)
  }
  list(total = total, binding_term = binding_term, binding_n = binding_n)
}

# .inla_largest_verified() --- the largest per-row INLA rung this
# package has verified to complete, read from
# inst/validation/benchmark_scaling.csv (S4/C2, FS-25 -- see the
# "inla-ceilings" study rows added there). Returns
# list(n =, random_effects =, run_date =) for the largest `n` among
# `method == "per_row"` rows carrying `outcome == "success"` AND a
# recorded `random_effects` count.
#
# The `random_effects` filter matters: this artefact also carries the
# `boundary` / `boundary-confirm` studies, which grow N on a FIXED,
# small random-effect structure (`random = ~ geno`, 6 genotypes) to
# measure the aggregated-vs-per-row row-count boundary -- a different
# question from the one an inla_program_failed refusal is answering.
# Their per_row rows reach a larger raw N (up to 1e6) than the
# genotype-growth ceilings study, but they say nothing about how large
# a LATENT FIELD this package has verified, which is what actually
# killed the FS-25 fit (a sparse-Cholesky solver cost, not a design-
# memory cost). Restricting to rows that recorded a random-effect
# count isolates the study that measured what matters here.
#
# NULL if the artefact is missing, unreadable, or carries no such row
# -- a refusal message then degrades to "not yet measured" rather than
# erroring on a missing CSV.
.inla_largest_verified <- function() {
  path <- system.file(
    "validation", "benchmark_scaling.csv",
    package = "flexyBayes"
  )
  if (!nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  df <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || !nrow(df)) {
    return(NULL)
  }
  needed <- c("method", "outcome", "n", "random_effects")
  if (!all(needed %in% names(df))) {
    return(NULL)
  }
  done <- df$method == "per_row" &
    df$outcome == "success" &
    !is.na(df$n)
  if (!any(done)) {
    return(NULL)
  }
  # Prefer a completed fit that recorded its latent-field size; fall back to
  # the largest completed per-row fit on record, flagging that its field
  # size was not recorded. Since the 2026-08-26 correction of the ceilings
  # study no completed fit on the genotype-growth ladder exists, so the
  # fallback is what a refusal reports today (1,000,000 rows, 360 cells).
  with_field <- done & !is.na(df$random_effects)
  ok <- if (any(with_field)) with_field else done
  df_ok <- df[ok, , drop = FALSE]
  best <- df_ok[which.max(df_ok$n), ]
  list(
    n = as.numeric(best$n),
    random_effects = if (
      "random_effects" %in% names(best) && !is.na(best$random_effects)
    ) {
      as.numeric(best$random_effects)
    } else {
      NA_real_
    },
    run_date = if ("run_date" %in% names(best)) {
      as.character(best$run_date)
    } else {
      NA_character_
    }
  )
}

# .inla_program_failed_condition() --- builds (does not raise) the
# `inla_program_failed` typed refusal condition object (S2/C2, FS-25),
# mirroring .fb_refusal_condition()'s own build-then-stop() idiom --
# every call site wraps this in stop(). Carries the design size, the
# largest size this package has verified, the binding quantity if
# known, and the remedies -- so the user never receives a raw engine-
# death string with no reason code and no next step
# (SCALE_STRATEGY_2026-08-22 S2).
.inla_program_failed_condition <- function(fb, data, engine_message) {
  n_rows <- if (is.data.frame(data)) as.numeric(nrow(data)) else NA_real_
  latent <- .inla_design_latent_effect_summary(fb, data)
  largest <- .inla_largest_verified()

  largest_txt <- if (is.null(largest)) {
    "not yet measured"
  } else {
    paste0(
      format(largest$n, big.mark = ",", scientific = FALSE),
      " rows",
      if (!is.na(largest$random_effects)) {
        paste0(
          " / ",
          format(largest$random_effects, big.mark = ",", scientific = FALSE),
          " random effects"
        )
      } else {
        ""
      },
      " (", largest$run_date, ")"
    )
  }

  binding_txt <- if (is.na(latent$binding_term)) {
    "not determined from the random terms"
  } else {
    paste0(
      latent$binding_term,
      " (", format(latent$binding_n, big.mark = ",", scientific = FALSE),
      " levels)"
    )
  }

  .fb_refusal_condition(
    reason_code = "inla_program_failed",
    message = paste0(
      "The INLA engine failed while fitting this design: \"",
      engine_message,
      "\" This design has ",
      format(n_rows, big.mark = ",", scientific = FALSE),
      " rows",
      if (!is.na(latent$total)) {
        paste0(
          " and an estimated ",
          format(latent$total, big.mark = ",", scientific = FALSE),
          " random-effect levels"
        )
      } else {
        ""
      },
      ". The largest design this package has verified an INLA per-row ",
      "fit to complete is ", largest_txt, ". The binding term (largest ",
      "single random-effect contributor) is ", binding_txt, ". ",
      "Remedies, in order of least to most disruptive: reduce the model ",
      "(fewer random-effect terms or coarser groupings); aggregate where ",
      "eligible (fb_plan() reports whether this design compresses); ",
      "switch engine (backend = \"brms\"); or add memory and retry."
    ),
    family_class = "flexybayes_inla_engine_refusal",
    n_rows = n_rows,
    n_latent_effects_est = latent$total,
    binding_term = latent$binding_term,
    binding_n_levels = if (is.finite(latent$binding_n)) {
      latent$binding_n
    } else {
      NA_real_
    },
    largest_verified_n = if (is.null(largest)) NA_real_ else largest$n,
    largest_verified_random_effects = if (is.null(largest)) {
      NA_real_
    } else {
      largest$random_effects
    },
    engine_message = engine_message
  )
}


# ---------------------------------------------------------------- #
# C4 -- legalise non-syntactic factor levels (FS-26)                 #
# ---------------------------------------------------------------- #

# .inla_legalise_factor_levels() --- rewrite every factor level the
# model touches to a syntactic name on a COPY of `data`, and return the
# copy alongside a per-variable map so the caller can print the user's
# own labels back through the accessors. A level was already syntactic
# in the common case (breeding-trial factor labels usually are, and
# level counts are small) -- this is a no-op there -- the map is still
# recorded uniformly so the accessors have one lookup path rather than
# a conditional one.
#
# The syntax check is NOT bare make.names(level) == level. INLA names a
# fixed-effect coefficient by concatenating the VARIABLE name onto the
# level with no separator (e.g. `env` + `2` -> `env2`), and that
# concatenated string is what has to survive INLA's own formula
# round-trip -- not the level alone. A purely-numeric level like `"2"`
# fails make.names() in isolation (a name cannot start with a digit),
# but `env2` is a perfectly valid, unambiguous identifier; legalising
# `"2"` to `"X2"` here anyway is a FALSE POSITIVE that broke the
# per-row vs. streamed-aggregated equivalence test the moment a design
# used integer-coded factor levels
# (test-aggregation-equivalence-backend.R uses `env <- factor(env)` on
# an integer vector) -- the streamed/aggregated route
# (R/emit_gaussian_aggregated.R) does not run this legalisation at all,
# so the two routes' coefficient names silently diverged. The check
# below probes with a throwaway letter-led prefix (any real `t$var` is
# itself already a legal R name, i.e. letter- or dot-led, so a letter
# prefix always reproduces the actual concatenation context) and only
# legalises a level whose PREFIXED form still changes under
# make.names() -- true for a space or other symbol-breaking character,
# false for a level that merely starts with a digit.
#
# Scans fixed_terms, random_terms and residual_terms for every
# variable name the model actually references (`$var` scalar terms,
# `$vars` interaction / at() terms), not just the response's immediate
# predictors -- a level living only inside a random or residual term
# (like FS-26's own `TRIAL`) is exactly as unrepresentable to INLA as
# one on the fixed side.
#
# Returns list(data = <copy of `data`, factor columns relabelled>,
# maps = <named list, one entry per relabelled variable, each a
# character vector keyed by the LEGAL label with the ORIGINAL label as
# its value -- .inla_restore_level_labels() reads it in that
# direction>).
.inla_legalise_factor_levels <- function(fb, data) {
  if (!is.data.frame(data)) {
    return(list(data = data, maps = list()))
  }
  vars <- unique(unlist(lapply(
    c(fb$fixed_terms, fb$random_terms, fb$residual_terms),
    function(t) {
      v <- if (!is.null(t$vars)) {
        as.character(t$vars)
      } else if (!is.null(t$var)) {
        as.character(t$var)
      } else {
        character(0L)
      }
      v[nzchar(v)]
    }
  )))
  vars <- vars[vars %in% names(data)]
  vars <- vars[vapply(vars, function(v) is.factor(data[[v]]), logical(1L))]

  if (!length(vars)) {
    return(list(data = data, maps = list()))
  }

  # A fixed, unlikely-to-collide letter-led prefix. Its only job is to
  # make every probed string start with a letter, so make.names()'s
  # verdict reflects whether the LEVEL's own characters are a problem,
  # not whether it happens to start with a digit -- concatenation onto
  # a real (letter-led) variable name already solves that part.
  probe_prefix <- "Q9zZ"

  maps <- list()
  for (v in vars) {
    original <- levels(data[[v]])
    probed <- paste0(probe_prefix, original)
    legal_probed <- make.names(probed, unique = TRUE)
    legal <- substring(legal_probed, nchar(probe_prefix) + 1L)
    if (!identical(original, legal)) {
      levels(data[[v]]) <- legal
    }
    maps[[v]] <- stats::setNames(original, legal)
  }

  list(data = data, maps = maps)
}

# .inla_restore_level_labels() --- translate a character vector of
# legal labels for variable `var` back to the user's own, using the
# map .inla_legalise_factor_levels() recorded on the fit
# (`object$level_labels`). Values not present in the map (a level the
# legalisation pass never touched, or a variable it never carried) pass
# through unchanged -- this is a display-time convenience, not a gate,
# so an unmapped value is not an error.
.inla_restore_level_labels <- function(level_labels, var, values) {
  map <- level_labels[[var]]
  if (is.null(map) || !length(values)) {
    return(values)
  }
  hit <- values %in% names(map)
  values[hit] <- unname(map[values[hit]])
  values
}

# .inla_restore_term_labels() --- the composite-name sibling of
# .inla_restore_level_labels(), for fixed-effect / interaction TERM
# strings rather than bare factor values. INLA (via R's own default
# treatment-contrast dummy coding) names a fixed-effect coefficient by
# pasting the variable name directly onto its level with no separator
# (e.g. `NAPPLIED2X.low.N`... in practice `NAPPLIED2low.N`, the legal
# level with no punctuation between var and level), so restoring the
# user's own label is a substring replace, once per relabelled
# variable x legal-level pair, over every term string. Order does not
# matter across variables (a legal level string is always prefixed by
# its own variable name in a term, so distinct variables cannot
# collide on the same substring); handles interaction terms (`a:b`)
# the same way, because the replace runs on the whole string.
#
# `object$level_labels` is `NULL` / empty for a fit this legalisation
# never touched (every non-INLA fit, and an INLA fit with no
# non-syntactic level) -- returns `terms` unchanged in that case, so
# call sites do not need their own NULL guard.
.inla_restore_term_labels <- function(level_labels, terms) {
  if (!length(level_labels) || !length(terms)) {
    return(terms)
  }
  for (var in names(level_labels)) {
    map <- level_labels[[var]]
    for (legal in names(map)) {
      original <- map[[legal]]
      if (identical(legal, original)) {
        next
      }
      needle <- paste0(var, legal)
      hit <- grepl(needle, terms, fixed = TRUE)
      if (any(hit)) {
        terms[hit] <- gsub(
          needle,
          paste0(var, original),
          terms[hit],
          fixed = TRUE
        )
      }
    }
  }
  terms
}

# .inla_restore_data_original_levels() --- the reverse of
# .inla_legalise_factor_levels(): a copy of `data` with every
# relabelled factor's LEVELS reverted from legal back to the user's
# own, via the same `level_labels` map. Level relabelling is 1:1 and
# order-preserving (make.names(..., unique = TRUE) never reorders or
# merges levels), so this round-trips exactly.
#
# The one consumer this exists for: .fb_fit_data()
# (R/ecosystem_support.R) is the single shared read-point every
# downstream-ecosystem integration (emmeans, marginaleffects, plot(),
# the missing-cell summary) uses for "the data to build a design
# matrix / reference grid from". Once coef.flexybayes_inla() restores
# the user's own coefficient names (C4/FS-26), a design matrix built
# from the still-legalised data would name its columns by the legal
# labels and no longer line up with the fit's own (now-original)
# coefficient names -- .fb_fixef_model_matrix()'s own reconciliation
# check exists precisely to catch a mismatch like that. Feeding it the
# original-labelled data keeps both sides speaking the same labels
# again, and every consumer of .fb_fit_data() gets the user's own
# labels for free rather than needing its own restoration call.
.inla_restore_data_original_levels <- function(level_labels, data) {
  if (!length(level_labels) || !is.data.frame(data)) {
    return(data)
  }
  for (v in names(level_labels)) {
    if (v %in% names(data) && is.factor(data[[v]])) {
      levels(data[[v]]) <- .inla_restore_level_labels(
        level_labels,
        v,
        levels(data[[v]])
      )
    }
  }
  data
}


# ---------------------------------------------------------------- #
# Public-style entry -- internal in v0.1                           #
# ---------------------------------------------------------------- #

emit_inla <- function(
  fb,
  data,
  known_matrices = list(),
  weights = NULL,
  n_samples = NULL,
  warmup = NULL,
  chains = NULL,
  seed = NULL,
  control = NULL,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  verbose = TRUE,
  return_code = FALSE,
  the_call = NULL,
  fixed = NULL,
  random = NULL,
  residual = NULL,
  family = NULL,
  link = NULL,
  data_name = NA_character_,
  na_action = NULL,
  requested_backend = NULL,
  requested_aggregate = NULL,
  ...
) {
  .check_fb_terms(fb, "`fb` must be an fb_terms object.")

  # Legalise non-syntactic factor levels (C4/FS-26). INLA expands a
  # factor fixed term to per-level names built as paste0(var, level)
  # and evaluates them internally; a level containing a space (or
  # anything else make.names() would touch, e.g. "low N") then fails
  # deep inside the retried inla-program call with an untyped
  # `object 'NAPPLIED2low N' not found`. Legalise every factor the
  # model uses on a COPY of `data` before anything below builds a
  # formula, an index column, or the data list INLA reads -- every
  # subsequent use of `data` in this function sees the legalised copy.
  # The map travels on the fit so summary() / coef() / ranef() /
  # predict(classify =) can translate the legal labels back to the
  # user's own before printing.
  level_map <- .inla_legalise_factor_levels(fb, data)
  data <- level_map$data

  # `seed` and `control` are Stan sampler settings threaded through from
  # the public entry point. INLA's nested Laplace approximation has no
  # random stream and no adaptation phase to tune, so both are genuine
  # no-ops here. Saying so once beats letting a user believe a fit was
  # pinned by a seed that never reached an engine.
  .note_sampler_args_ignored("INLA", seed = seed, control = control)

  .check_installed(
    "INLA",
    "Package 'INLA' is required for backend = \"inla\". ",
    "Install via:\n",
    "  install.packages('INLA',\n",
    "    repos = c(CRAN = 'https://cran.r-project.org',\n",
    "             INLA = 'https://inla.r-inla-download.org/R/stable'),\n",
    "    dep = TRUE)"
  )

  # INLA verification gate. simple_slope_uncor
  # passes lgm_gate (the allowlist includes it) so dispatch routes
  # here, but the INLA mapping for (x || g) is registered only when
  # the three-arbitrator verification test on a simple fixture has
  # passed and the artefact at
  # inst/extdata/inla-verification/simple_slope_uncor.rds records
  # pass = TRUE. Until then, refuse explicitly with a deferral
  # message naming the workaround (backend = "brms"). This is the
  # "no silent translation of an unverified mapping" policy.
  if (.has_simple_slope_uncor(fb)) {
    .check_inla_verification_simple_slope_uncor()
  }

  # Build INLA formula and family. priors_to_inla() returns a list
  # keyed by group / smoother var ("sigma" for residual). We thread
  # the per-term entries into f(..., hyper = list(prec = ...)) and
  # the "sigma" entry into control.family (residual precision).
  hyper_ctrl <- if (inherits(fb$priors, "fb_prior")) {
    priors_to_inla(fb$priors)
  } else if (.fb_prior_scalar_supplied(fb, "vc_sd")) {
    # The legacy scalar route. Until 0.9.2 this branch did not exist: the
    # INLA path read only an fb_prior, so a fit passed `prior_vc_sd` ran
    # under INLA's own log-gamma precision default while
    # `prior_summary()` printed the lognormal in its header -- the
    # accessor built to answer "what prior did this fit use?" naming a
    # prior the fit did not use.
    .priors_legacy_to_inla(fb, prior_vc_sd)
  } else {
    list()
  }
  # An `sd(group = ...)` row whose group is not a variance component of
  # this model would key a `hyper` entry no f() term reads -- the INLA
  # half of field-sweep FS-20, where the same mistake on brms surfaced as
  # brms's own parser error and here was not reported at all.
  .check_inla_hyper_keys_reachable(fb, hyper_ctrl)
  # A variable that keys both a fixed term and an f() term is a formula
  # INLA refuses -- its own message is "Key [x] is used twice and that is
  # not allowed". Caught here so the user is told which variable and what
  # to do, rather than reading the generic subprocess-exit message
  # (field-sweep FS-17).
  .check_inla_variable_not_used_twice(fb)
  # `prior_fixed_sd` reaches INLA through control.fixed, which is where
  # INLA keeps the fixed-effect prior, and so do the `b()` rows of an
  # fb_prior(). Empty unless the caller wrote one of them, so an
  # unsupplied call keeps INLA's own defaults.
  control_fixed <- .build_inla_control_fixed(
    fb,
    prior_fixed_sd,
    data = data
  )
  # v0.3.10: the formula builder consults
  # known_matrices to count blocks per blocks-format vm/ped term,
  # so it must run AFTER the data_inla setup loop (which validates
  # the blocks structure). The previous v0.3.7 ordering --- formula
  # build first, then data_inla setup --- is preserved for legacy
  # paths that do not carry blocks-format terms via the
  # known_matrices default empty list.
  inla_family <- .resolve_inla_family(fb)
  control_family <- .build_inla_control_family(hyper_ctrl, family = inla_family)

  if (return_code) {
    return(invisible(list(
      formula = .build_inla_formula(
        fb,
        hyper_ctrl = hyper_ctrl,
        known_matrices = known_matrices
      ),
      family = inla_family,
      hyper = hyper_ctrl,
      control_family = control_family,
      control_fixed = control_fixed
    )))
  }

  # When a vm/ped term routes via
  # f(<var>_id, model = "generic0", Cmatrix = <symbol>), INLA
  # resolves both the integer index column AND the Cmatrix symbol
  # from the `data` argument. Convert data.frame to a list, attach
  # an integer index column per vm/ped term, and merge
  # known_matrices so the precision matrices are in INLA's lookup
  # scope at fit time. Defensive collision check refuses when a
  # known_matrices entry shadows a data column, preserving the
  # package's no-silent-failure identity.
  data_inla <- as.list(data)
  # v0.3.10: blocks-format vm/ped terms generate
  # K per-block precision matrices that need to live in `data_inla`
  # alongside the integer index columns, since INLA looks up Cmatrix
  # symbols from the data list at fit time. The slot is also keyed by
  # the term's variable name so multiple blocks-format terms do not
  # collide on a generic name.
  blocks_precision_for_inla <- list()
  for (term in fb$random_terms) {
    if (term$type %in% c("vm", "ped")) {
      idx_name <- paste0(term$var, "_id")
      f <- factor(data[[term$var]])
      data_inla[[idx_name]] <- as.integer(f)

      # v0.3.8: enforce known-matrix dim + level
      # alignment on the INLA emit path. Pre-v0.3.8 the
      # .build_inla_formula() body deferred this to the user (it
      # documented `f(<var>_id, model =
      # "generic0", Cmatrix = Q)` and the user-responsibility
      # caveat). The validator refuses cleanly when (a) the matrix
      # dim does not match the level count, (b) the matrix carries
      # dimnames that do not match levels(<var>), or (c) the matrix
      # carries dimnames that match the level set but in a different
      # order from levels(<var>) -- INLA's generic0 requires
      # positional alignment with the integer index column built
      # above. When dimnames are absent the alignment downgrades to
      # the dim check + a future <fb_plan>-surfaced caution
      # (planned at v0.3.8).
      cov <- term$cov_representation
      if (
        !is.null(cov) &&
          cov$format %in% c("precision", "pedigree_sparse_precision")
      ) {
        Q <- known_matrices[[cov$data]]
        .validate_precision_input(
          Q,
          name = cov$data,
          group_var = term$var,
          expected_n = nlevels(f),
          fit_levels = levels(f)
        )
      }

      # v0.3.10: blocks-format vm/ped routes
      # through K independent f(<var>_id_block_<k>, model =
      # "generic0", Cmatrix = <symbol>_Q_<k>) calls. Pre-compute one
      # integer index column per block (NA outside the block) and one
      # per-block precision matrix (Q_k = solve(V_k)). Convention:
      # levels(<var>) are partitioned into K consecutive chunks of
      # sizes n_1, ..., n_K --- the validator's `block_sizes` slot is
      # the source of truth. The within-block position is the
      # 1-based offset of that level within its block.
      if (!is.null(cov) && identical(cov$format, "blocks")) {
        blocks_meta <- .validate_blocks_input(
          known_matrices[[cov$data]],
          name = cov$data,
          group_var = term$var,
          expected_n = nlevels(f)
        )
        block_sizes <- blocks_meta$block_sizes
        block_ends <- cumsum(block_sizes)
        block_starts <- c(1L, utils::head(block_ends, -1L) + 1L)
        level_ids <- as.integer(f)
        for (k in seq_along(block_sizes)) {
          idx_block_name <- paste0(term$var, "_id_block_", k)
          q_block_name <- paste0(cov$data, "_Q_", k)
          within_block <- level_ids - block_starts[[k]] + 1L
          within_block[
            level_ids < block_starts[[k]] |
              level_ids > block_ends[[k]]
          ] <- NA_integer_
          data_inla[[idx_block_name]] <- within_block

          V_k <- blocks_meta$blocks[[k]]
          Q_k <- if (inherits(V_k, "Matrix")) {
            Matrix::solve(Matrix::forceSymmetric(V_k))
          } else {
            solve(as.matrix(V_k))
          }
          blocks_precision_for_inla[[q_block_name]] <- Q_k
        }
      }
    }
  }
  if (length(blocks_precision_for_inla)) {
    conflicts <- intersect(names(data_inla), names(blocks_precision_for_inla))
    if (length(conflicts)) {
      .stop_structured_cov_refusal(
        reason_code = "known_matrices_data_name_collision",
        message = paste0(
          "blocks-format precision matrix name(s) '",
          paste(conflicts, collapse = "', '"),
          "' collide with data column names. Rename the blocks ",
          "carrier (the per-block precision matrices are named ",
          "<carrier>_Q_<k>) to disambiguate."
        ),
        conflicts = conflicts
      )
    }
    data_inla <- c(data_inla, blocks_precision_for_inla)
  }
  if (length(known_matrices)) {
    # v0.3.10: blocks-format carriers are
    # lists-of-matrices, not single matrices --- they belong in
    # blocks_precision_for_inla (in their per-block Q_k form), not
    # in data_inla under their user-facing symbol. Drop them before
    # the bulk append so INLA does not see a list column.
    blocks_symbols <- character(0L)
    for (term in fb$random_terms) {
      if (
        term$type %in%
          c("vm", "ped") &&
          !is.null(term$cov_representation) &&
          identical(term$cov_representation$format, "blocks")
      ) {
        blocks_symbols <- c(blocks_symbols, term$cov_representation$data)
      }
    }
    km_for_inla <- known_matrices[
      setdiff(names(known_matrices), blocks_symbols)
    ]
    if (length(km_for_inla)) {
      conflicts <- intersect(names(data_inla), names(km_for_inla))
      if (length(conflicts)) {
        .stop_structured_cov_refusal(
          reason_code = "known_matrices_data_name_collision",
          message = paste0(
            "known_matrices entry/entries '",
            paste(conflicts, collapse = "', '"),
            "' collide with data column names; rename to disambiguate ",
            "(INLA's data-list lookup is name-keyed)."
          ),
          conflicts = conflicts
        )
      }
      data_inla <- c(data_inla, km_for_inla)
    }
  }

  # AR1 / separable AR1xAR1 spatial latent-field index columns. INLA's ar1
  # model + grouped-AR1 need consecutive-integer indices; build <var>_id for
  # the 1D term and <row>_id / <col>_id for the separable term.
  #
  # 0.9.0: the field is read from the random terms only. It used to be read
  # from the residual terms as well, reinterpreting an ASReml residual as a
  # latent field plus the Gaussian observation nugget -- four parameters
  # under a three-parameter name. The residual spelling now refuses at
  # dispatch and the field is written where it belongs.
  #
  # The field is faithful only with one observation per grid node: with a
  # gap or a replicate the latent index no longer matches the design, so the
  # instance is validated here (the class is admitted by the gate) and an
  # incomplete or replicated grid refuses rather than being approximated.
  for (term in fb$random_terms) {
    ttype <- term$type %||% ""
    if (!ttype %in% c("ar1", "ar1_spatial")) {
      next
    }
    idx_vars <- if (identical(ttype, "ar1")) {
      term$var
    } else {
      c(term$row_var, term$col_var)
    }
    for (v in idx_vars) {
      data_inla[[paste0(v, "_id")]] <- as.integer(factor(data[[v]]))
    }
    if (identical(ttype, "ar1_spatial")) {
      # C5/FS-27: the per-trial spelling's replicate index. Built the
      # same way as the row / col index columns; INLA's replicate =
      # reads it to fit one independent realisation of the field per
      # distinct value.
      if (!is.null(term$at_var)) {
        data_inla[[paste0(term$at_var, "_id")]] <-
          as.integer(factor(data[[term$at_var]]))
      }
      n_row <- length(unique(data[[term$row_var]]))
      n_col <- length(unique(data[[term$col_var]]))
      n_nodes <- n_row * n_col
      if (is.null(term$at_var)) {
        combos <- interaction(
          data[[term$row_var]],
          data[[term$col_var]],
          drop = TRUE
        )
        one_obs_per_node <- nrow(data) == n_nodes && !any(duplicated(combos))
        bad_at_levels <- if (one_obs_per_node) character(0L) else "<the design>"
      } else {
        # C5/FS-27: "the complete-lattice check runs per level" -- one
        # field per level of at_var means each level's OWN subset must
        # independently form a complete row x col lattice; the global
        # row count check the single-field case uses does not
        # generalise (every level legitimately re-uses the same row /
        # col labels, so the whole table is n_at times n_nodes, not
        # n_nodes).
        at_levels <- unique(data[[term$at_var]])
        bad_at_levels <- character(0L)
        for (lev in at_levels) {
          sub <- data[data[[term$at_var]] == lev, , drop = FALSE]
          combos <- interaction(
            sub[[term$row_var]],
            sub[[term$col_var]],
            drop = TRUE
          )
          ok <- nrow(sub) == n_nodes && !any(duplicated(combos))
          if (!ok) {
            bad_at_levels <- c(bad_at_levels, as.character(lev))
          }
        }
      }
      if (length(bad_at_levels)) {
        stop(.fb_refusal_condition(
          reason_code = "ar1_spatial_requires_complete_grid",
          message = paste0(
            if (!is.null(term$at_var)) {
              paste0(
                "The level(s) of `", term$at_var, "` that do not ",
                "independently form a complete ", n_row, " x ", n_col,
                " (", term$row_var, ", ", term$col_var, ") lattice: ",
                paste(bad_at_levels, collapse = ", "), ". "
              )
            } else {
              paste0(
                "The data carry ", nrow(data), " rows for the ", n_row,
                " x ", n_col, " (", term$row_var, ", ", term$col_var,
                ") array, which has ", n_nodes, " nodes. "
              )
            },
            .ar1_term_spelling(term),
            " is fitted as a separable ",
            "autoregressive field indexed by that array, so each ",
            "combination of the two factors has to identify exactly one ",
            "unit",
            if (!is.null(term$at_var)) {
              paste0(", WITHIN each level of `", term$at_var, "`")
            } else {
              ""
            },
            " -- the same requirement ASReml states for ",
            "`residual = ~ ar1:ar1`, and the same reason it refuses a ",
            "trial with plots deleted. Keep the unobserved plots as ",
            "design cells with na_action = \"augment\" (the default), ",
            "which is ASReml's na.method(y = \"include\"); supply the ",
            "field-book rows for any plot absent from the data, with the ",
            "design columns filled in and the response set to NA, or pad ",
            "the array with fb_complete_grid(); or aggregate replicates ",
            "to node means."
          ),
          family_class = "flexybayes_ar1_spatial_refusal"
        ))
      }
    }
  }

  # Formula build runs here (post data_inla setup) so blocks-format
  # vm/ped terms can resolve their K-count from known_matrices.
  inla_form <- .build_inla_formula(
    fb,
    hyper_ctrl = hyper_ctrl,
    known_matrices = known_matrices
  )
  if (verbose) {
    cat(
      "\n-- flexyBayes: INLA fit ",
      paste(rep("-", 40), collapse = ""),
      "\n",
      sep = ""
    )
    cat("  formula: ", .fb_display_inla_formula(inla_form), "\n", sep = "")
    cat("  family:  ", inla_family, "\n", sep = "")
    cat(paste(rep("-", 60), collapse = ""), "\n\n")
  }

  # Fit
  #
  # control.fixed is spliced in only when it carries something. An empty
  # list means the caller did not supply `prior_fixed_sd`, and the
  # argument is then left out of the call entirely rather than passed
  # empty, so a fit that asks for nothing is byte-identical to the call
  # this emit made before the argument was wired.
  t0 <- proc.time()
  inla_args <- list(
    formula = inla_form,
    family = inla_family,
    data = data_inla,
    control.compute = list(
      config = TRUE,
      return.marginals = TRUE,
      dic = TRUE,
      waic = TRUE
    ),
    control.family = control_family,
    ...
  )
  if (length(control_fixed)) {
    inla_args$control.fixed <- control_fixed
  }
  # C6: weights lowered for the Gaussian family, identity link only --
  # INLA's scale = is a per-observation precision multiplier on the
  # Gaussian / Student-T likelihood (?INLA::inla: "Fixed (optional)
  # scale parameters of the precision ... default rep(1, n.data)"),
  # giving the same Var(y_i) = sigma^2 / w_i semantics brms's sigma
  # distributional offset lowers to on that engine (R/emit_brms.R's
  # .FB_BRMS_WEIGHTS_OFFSET_COL -- NOT brms's own weights() addition
  # term, a different quantity; that constant's banner carries the
  # grounding trail). Everything this emit cannot honour (non-gaussian,
  # non-identity link, the aggregated route) is already refused
  # upstream by .refuse_unsupported_weights() / the aggregation gate
  # (R/dispatch.R) before this call is reached -- this site only
  # carries the value through. Spliced in only when weights are
  # actually non-constant, for the same "byte-identical to the
  # unweighted call" reason control.fixed follows above -- a constant
  # vector (any value, not only 1) is the unweighted model under a
  # different spelling, and scale = <constant != 1> would silently
  # rescale the reported sigma rather than leaving it unweighted.
  ir_weights <- .fb_ir_weights(fb) %||% weights
  if (.fb_weights_nonconstant(ir_weights)) {
    inla_args$scale <- ir_weights
  }
  fit <- tryCatch(
    .inla_call(inla_args),
    error = function(e) {
      if (.inla_is_program_death_message(conditionMessage(e))) {
        stop(.inla_program_failed_condition(
          fb = fb,
          data = data,
          engine_message = conditionMessage(e)
        ))
      }
      # Not a recognised engine-death signature -- this is flexyBayes's
      # own code (a bad argument the emit built, a formula it mis-
      # constructed) or an INLA-side argument-validation error triggered
      # by one. Re-raise verbatim: wrapping it as an engine-death refusal
      # would misname the failure and hide it from anyone debugging the
      # emit itself (S2 draws this line explicitly).
      stop(e)
    }
  )
  if (is.null(fit) || !length(fit)) {
    stop(.inla_program_failed_condition(
      fb = fb,
      data = data,
      engine_message = paste0(
        "INLA::inla() returned NULL / an empty result with no error ",
        "condition."
      )
    ))
  }
  elapsed <- unname((proc.time() - t0)["elapsed"])

  # Post-fit numerical-confirm gate
  num_check <- .lgm_check_numerical(fit)
  if (!num_check$pass) {
    warning(
      "flexyBayes: INLA numerical-confirm gate flagged: ",
      paste(num_check$reasons, collapse = "; "),
      ". Treat the fit with caution.",
      call. = FALSE
    )
  }

  # Wrap as a flexybayes_inla object
  out <- structure(
    list(
      inla = fit,
      fb = fb,
      # `data` carries the level-legalised copy (C4/FS-26) -- internal
      # consumers (predict()'s design-matrix reconstruction chief among
      # them) need the same labels the INLA fit itself was built
      # against. `level_labels` is the map back to the user's own,
      # which summary() / coef() / ranef() / predict(classify =) read
      # at the display boundary.
      data = data,
      level_labels = level_map$maps,
      num_check = num_check,
      extras = list(
        summary = list(
          fixed = fit$summary.fixed,
          random = fit$summary.random,
          hyperpar = fit$summary.hyperpar,
          fitted = if (!is.null(fit$summary.fitted.values)) {
            utils::head(fit$summary.fitted.values)
          } else {
            NULL
          }
        ),
        model_info = list(
          n_obs = nrow(data),
          n_fixed = if (!is.null(fit$summary.fixed)) {
            nrow(fit$summary.fixed)
          } else {
            0L
          },
          n_random = length(fit$summary.random),
          n_hyper = if (!is.null(fit$summary.hyperpar)) {
            nrow(fit$summary.hyperpar)
          } else {
            0L
          },
          family = inla_family,
          link = fb$link
        ),
        # The argument record update() rebuilds the call from. Before
        # 0.9.1 this held six fields against brms's fifteen, and the
        # difference -- not anything about the engine -- is why update()
        # refused every INLA fit as `update_call_not_reconstructable`.
        #
        # `n_samples`, `warmup` and `chains` are recorded as what INLA
        # actually used, which is nothing: the nested Laplace
        # approximation has no chains to run and no warmup to discard, so
        # the honest record is NULL rather than the sampler settings the
        # call happened to carry. update() passes each recorded field
        # straight back to flexybayes(), where a NULL sampler setting is
        # a no-op on this engine.
        call_info = list(
          fixed = fixed,
          random = random,
          residual = residual,
          data_name = data_name,
          family = family,
          link = link,
          known_matrices = known_matrices,
          weights = weights,
          n_samples = NULL,
          warmup = NULL,
          chains = NULL,
          seed = seed,
          control = control,
          prior_fixed_sd = prior_fixed_sd,
          prior_vc_sd = prior_vc_sd,
          na_action = na_action,
          # The engine and the representation the call ASKED for, and the
          # reporting it ran under. update() re-issues all three; without
          # them an identity re-fit of a per-row fit came back aggregated,
          # whose summary speaks a different dialect and carries no
          # $varcomp. `backend` is the request ("auto"), not the engine
          # the request landed on: it is a policy, and a re-fit of a
          # changed model must be free to route again.
          backend = requested_backend,
          aggregate = requested_aggregate,
          verbose = verbose
        ),
        # Thread the IR's smooth-objects slot through to the
        # extras so predict() on the INLA path can also use
        # mgcv::Predict.matrix() when smooths are present. The INLA
        # backend itself does not currently fit s() smooths via mgcv
        # (it uses INLA's own rw2 path), so this slot will normally be
        # an empty list; threaded for shape uniformity with emit_brms.
        parse_info = list(
          smooths = .collect_smooths(fb$random_terms),
          # Threaded for shape uniformity with emit_brms so
          # cross-backend accessors (genomic_summary(), fb_structured_cov())
          # can locate the random terms the same way on every engine.
          random = fb$random_terms
        ),
        run_time = elapsed,
        the_call = the_call,
        formula = inla_form,
        # What has to match before this fit can be compared with another
        # engine's -- the same slot, built by the same function, as the
        # brms emit fills. See R/model_fingerprint.R.
        fingerprint = .fb_model_fingerprint(fb, data, prior_vc_sd = prior_vc_sd)
      )
    ),
    # The INLA fit wears the shared `flexybayes` parent, as the brms fit
    # does. Before 0.9.0 it did not, which split the class graph in two and
    # forced a parallel S3 method for every generic. The parent is safe only
    # because the five parent methods with no INLA sibling
    # (confint, model.matrix, update, anova, logLik) now resolve their inputs
    # from the slots an object actually carries and refuse by name when it
    # carries none -- see R/methods.R.
    class = c("flexybayes_inla", "flexybayes", "list")
  )

  # Variance components in the shape every other engine writes, so
  # summary(fit)$varcomp, tidy(effects = "random") and the variance plot
  # answer on INLA too. Built after the wrap because the canonical-name
  # mapper reads the fitted object, not the raw INLA return.
  #
  # A summary table is a reading of the fit, not part of it, so a failure
  # here must not destroy a posterior the engine computed successfully.
  # It is surfaced as a warning naming the cause rather than swallowed:
  # the fit comes back with an empty variance-component table and the
  # user is told why.
  out$extras$variance_comps <- tryCatch(
    .inla_variance_comps(out),
    error = function(e) {
      warning(
        "flexyBayes: the INLA fit succeeded but its variance-component ",
        "table could not be built (",
        conditionMessage(e),
        "). ",
        "summary(fit)$varcomp will be empty; the posterior itself is ",
        "unaffected and INLA's own hyperparameter table is at ",
        "fit$inla$summary.hyperpar.",
        call. = FALSE
      )
      NULL
    }
  )

  # A field that lost its identification is a fact about this fit that
  # nothing else on the printed surface states. Raised after the
  # variance-component table is built, because it is read off that table.
  .fb_warn_spatial_field_collapsed(out)
  .fb_warn_boundary_collapse(out)

  out
}

# ---------------------------------------------------------------- #
# fb_terms -> INLA formula                                         #
# ---------------------------------------------------------------- #

.build_inla_formula <- function(
  fb,
  hyper_ctrl = list(),
  known_matrices = list()
) {
  rhs_terms <- character(0)

  # Intercept
  rhs_terms <- c(rhs_terms, if (fb$intercept) "1" else "0")

  # Fixed-effect terms
  for (term in fb$fixed_terms) {
    contrib <- switch(
      term$type,
      "factor" = term$var,
      "continuous" = term$var,
      "interaction" = paste(term$vars, collapse = ":"),
      "factor_interaction" = paste(term$vars, collapse = ":"),
      # v0.2.6 -- factor:continuous indexed
      # interaction. INLA's native `f:x` notation produces the same
      # treatment-coded indexed-slope shape as base R model.matrix();
      # the verification gate
      # (.lgm_check_factor_numeric_interaction_inla_verified) confirmed
      # posterior agreement with lme4 (and a since-withdrawn engine) on a
      # gaussian-identity fixture before this branch became reachable.
      "factor_numeric_interaction" = paste(term$vars, collapse = ":"),
      "expression" = term$label,
      # This stop() is an internal contract-violation
      # assertion. The gate's .lgm_check_fixed_term_inla_support()
      # refuses any fixed term type outside the allowlist before
      # we reach here. A fired assertion means the gate and the
      # emit are out of sync -- a flexyBayes-side bug, not a
      # user-facing refusal.
      stop(
        "lgm_gate broken contract: fixed term type \"",
        term$type,
        "\" reached emit_inla() despite a passing ",
        "gate. This is a flexyBayes internal bug -- the gate's ",
        ".lgm_check_fixed_term_inla_support() allowlist is out ",
        "of sync with .build_inla_formula()'s switch(). Please ",
        "file an issue.",
        call. = FALSE
      )
    )
    rhs_terms <- c(rhs_terms, contrib)
  }

  # Random-effect terms via INLA's f(). When the user-supplied
  # priors_to_inla() output names a hyperparameter for this term,
  # splice it in as `hyper = list(prec = list(prior = ..., param = ...))`
  # so the user's prior actually reaches INLA.
  for (term in fb$random_terms) {
    # Uncorrelated random slope `(x || g)`
    # maps to two independent IID hyperparameter slots --
    #   f(<g>, model = "iid") + f(<g>_for_slope, <x>, model = "iid")
    # `<g>_for_slope` is a fresh INLA index keyed off the same
    # grouping factor so the two precisions remain independent
    # (INLA aggregates two `f()` calls under the same name).
    # Only reachable when the verification artefact records
    # pass = TRUE -- the emit_inla() entry refuses upfront otherwise.
    if (identical(term$type, "simple_slope_uncor")) {
      sv <- term$slope_var
      gtag <- term$var
      gslope <- paste0(gtag, "_for_slope_", sv)
      int_key <- gtag
      slope_key <- paste0(sv, "_", gtag)
      if (isTRUE(term$with_intercept)) {
        hyper_int <- .inla_hyper_arg(hyper_ctrl[[int_key]])
        rhs_terms <- c(
          rhs_terms,
          paste0(
            "f(",
            gtag,
            ", model = \"iid\"",
            if (nzchar(hyper_int)) paste0(", ", hyper_int) else "",
            ")"
          )
        )
      }
      hyper_slope <- .inla_hyper_arg(hyper_ctrl[[slope_key]])
      rhs_terms <- c(
        rhs_terms,
        paste0(
          "f(",
          gslope,
          ", ",
          sv,
          ", model = \"iid\"",
          if (nzchar(hyper_slope)) paste0(", ", hyper_slope) else "",
          ")"
        )
      )
      next
    }
    # vm/ped with format = "precision" or
    # "pedigree_sparse_precision" route through INLA's user-defined
    # precision interface, f(idx, model = "generic0", Cmatrix = Q).
    # INLA's generic0 model requires the first arg to be an integer
    # index into the rows/cols of Cmatrix; we pre-compute the index
    # column <var>_id in data_inla (see emit_inla() body above) and
    # reference it here so summary.random[[<var>_id]] carries the
    # natural name in the fit object. v0.3.8:
    # the level ordering of the factor against the row ordering of Q
    # is now enforced by .validate_precision_input() (called from
    # emit_inla() body above with expected_n = nlevels(f) and
    # fit_levels = levels(f)) when Q carries dimnames; the
    # pre-v0.3.8 user-responsibility caveat is no longer load-bearing.
    if (term$type %in% c("vm", "ped")) {
      cov <- term$cov_representation
      key <- term$var
      hyper_str <- .inla_hyper_arg(hyper_ctrl[[key]])

      # v0.3.10: blocks-format vm/ped emits K
      # independent f() calls, one per block, each with its own
      # generic0 precision and an INLA-default per-block precision
      # hyperparameter. The K integer index columns + K per-block
      # precision matrices are pre-computed in emit_inla()'s
      # data_inla setup loop; here we only need the block count,
      # which the resolved known_matrices entry surfaces directly.
      # The per-block hyperparameter is left at INLA's
      # loggamma(1, 5e-5) default --- no shared per-block prior at
      # v0.3.10 (the priors-to-inla hyper_ctrl entry is keyed by
      # term$var, so all K blocks share its parameters when the user
      # supplies one via fb_prior()).
      if (!is.null(cov) && identical(cov$format, "blocks")) {
        blocks_value <- known_matrices[[cov$data]]
        K <- if (is.list(blocks_value)) length(blocks_value) else 0L
        if (K == 0L) {
          stop(
            "lgm_gate broken contract: blocks-format vm/ped term ",
            "reached .build_inla_formula() but known_matrices[[\"",
            cov$data,
            "\"]] is empty or not a list. This is a ",
            "flexyBayes internal bug --- the gate or emit_inla()'s ",
            "validator should have refused first.",
            call. = FALSE
          )
        }
        contribs <- character(K)
        for (k in seq_len(K)) {
          idx_block_name <- paste0(term$var, "_id_block_", k)
          q_block_name <- paste0(cov$data, "_Q_", k)
          contribs[[k]] <- paste0(
            "f(",
            idx_block_name,
            ", model = \"generic0\", Cmatrix = ",
            q_block_name,
            if (nzchar(hyper_str)) paste0(", ", hyper_str) else "",
            ")"
          )
        }
        rhs_terms <- c(rhs_terms, contribs)
        next
      }

      idx_name <- paste0(term$var, "_id")
      contrib <- paste0(
        "f(",
        idx_name,
        ", model = \"generic0\", Cmatrix = ",
        cov$data,
        if (nzchar(hyper_str)) paste0(", ", hyper_str) else "",
        ")"
      )
      rhs_terms <- c(rhs_terms, contrib)
      next
    }
    # WP16: AR1 (1D) / separable AR1xAR1 spatial random field.
    if (term$type %in% c("ar1", "ar1_spatial")) {
      rhs_terms <- c(rhs_terms, .inla_ar1_field_term(term))
      next
    }
    key <- switch(
      term$type,
      "simple" = ,
      "ide" = ,
      "id" = term$var,
      "spline" = term$var,
      NULL
    )
    model_name <- switch(
      term$type,
      "simple" = ,
      "ide" = ,
      "id" = "iid",
      "spline" = "rw2",
      # Internal contract-violation assertion. See the
      # fixed-term site above for the rationale; the gate's
      # .lgm_check_random_term_inla_support() owns this guard.
      stop(
        "lgm_gate broken contract: random term type \"",
        term$type,
        "\" reached emit_inla() despite a passing ",
        "gate. This is a flexyBayes internal bug -- the gate's ",
        ".lgm_check_random_term_inla_support() allowlist is ",
        "out of sync with .build_inla_formula()'s switch(). ",
        "Please file an issue.",
        call. = FALSE
      )
    )
    hyper_str <- .inla_hyper_arg(hyper_ctrl[[key]])
    contrib <- paste0(
      "f(",
      term$var,
      ", model = \"",
      model_name,
      "\"",
      if (nzchar(hyper_str)) paste0(", ", hyper_str) else "",
      ")"
    )
    rhs_terms <- c(rhs_terms, contrib)
  }

  # Residual terms. INLA folds the residual variance into the likelihood, so
  # the only residual form that reaches here is the homogeneous one, which
  # contributes nothing to the formula. The separable autoregressive
  # residual used to be reinterpreted here as a latent field; it now refuses
  # at dispatch and the field is written on the random side.
  for (term in fb$residual_terms) {
    if (term$type != "units") {
      # Internal contract-violation assertion. See the
      # fixed-term site above for the rationale; the gate's
      # .lgm_check_residual_term_inla_support() owns this guard.
      stop(
        "lgm_gate broken contract: residual term type \"",
        term$type,
        "\" reached emit_inla() despite a passing ",
        "gate. This is a flexyBayes internal bug -- the gate's ",
        ".lgm_check_residual_term_inla_support() allowlist is out ",
        "of sync with .build_inla_formula()'s residual guard. ",
        "Please file an issue.",
        call. = FALSE
      )
    }
  }

  rhs <- paste(rhs_terms, collapse = " + ")
  stats::as.formula(paste0(fb$response, " ~ ", rhs))
}

# Build the `hyper = list(prec = list(prior = ..., param = ...))`
# substring spliced into a random-effect f() call. Returns "" when
# no entry is present so default INLA priors apply.
.inla_hyper_arg <- function(entry) {
  if (is.null(entry)) {
    return("")
  }
  prior <- entry$prior
  param <- entry$param
  if (is.null(prior)) {
    return("")
  }
  body <- paste0("prior = ", deparse(prior))
  if (length(param)) {
    body <- paste0(body, ", param = c(", paste(param, collapse = ", "), ")")
  }
  paste0("hyper = list(prec = list(", body, "))")
}

# ---------------------------------------------------------------- #
# Displaying the emitted formula                                    #
# ---------------------------------------------------------------- #
#
# A uniform-on-SD prior reaches INLA as a C expression evaluated per
# hyperparameter -- roughly a hundred characters of `theta`, `-1.0e10`
# and `return(...)` spliced into the formula. That is the right thing to
# send the engine and the wrong thing to show a reader: it buries the
# model in the encoding of one prior, and the two displays that print the
# formula (the emit-time banner and print()) were unreadable on any fit
# carrying a non-default variance prior.
#
# The compression is for those two displays only. `$extras$formula` and
# the formula INLA itself holds keep every character, because they are
# what the fit was actually run from.

# .fb_signif_string() --- a number at display precision.
#
# @noRd
# @keywords internal
.fb_signif_string <- function(x, digits = 4L) {
  format(signif(as.numeric(x), digits), trim = TRUE, scientific = FALSE)
}

# .fb_hyper_prior_tag() --- name the prior a hyper block encodes.
#
# Read off the block itself, which is the resolved prior rendered for the
# engine: the uniform-on-SD and half-normal-on-SD expressions carry their
# own bounds and scale, and a named INLA prior carries its name and
# parameters. An encoding this does not recognise is tagged as an
# expression rather than guessed at.
#
# @noRd
# @keywords internal
.fb_hyper_prior_tag <- function(body) {
  num <- "([-+0-9.eE]+)"
  m <- regmatches(
    body,
    regexec(paste0("expression: *L=", num, "; *U=", num), body)
  )[[1L]]
  if (length(m) == 3L) {
    return(sprintf(
      "<prior: uniform-SD(%s, %s)>",
      .fb_signif_string(m[[2L]]),
      .fb_signif_string(m[[3L]])
    ))
  }
  m <- regmatches(
    body,
    regexec(paste0("expression: *U=", num), body)
  )[[1L]]
  if (length(m) == 2L) {
    return(sprintf("<prior: uniform-SD(0, %s)>", .fb_signif_string(m[[2L]])))
  }
  m <- regmatches(
    body,
    regexec(paste0("expression: *s=", num), body)
  )[[1L]]
  if (length(m) == 2L) {
    # Both SD-scale expression priors open with the same `s=` scale, so
    # they are told apart on the term only one of them carries: the
    # lognormal is quadratic in theta, the half-normal is not.
    if (grepl("theta*theta", body, fixed = TRUE)) {
      return(sprintf(
        "<prior: lognormal-SD(0, %s)>",
        .fb_signif_string(m[[2L]])
      ))
    }
    return(sprintf(
      "<prior: half-normal-SD(%s)>",
      .fb_signif_string(m[[2L]])
    ))
  }

  nm <- regmatches(body, regexec("prior *= *\"([^\"]+)\"", body))[[1L]]
  if (length(nm) == 2L && !startsWith(nm[[2L]], "expression:")) {
    par <- regmatches(body, regexec("param *= *c\\(([^)]*)\\)", body))[[1L]]
    if (length(par) == 2L) {
      vals <- trimws(strsplit(par[[2L]], ",", fixed = TRUE)[[1L]])
      vals <- vapply(vals, .fb_signif_string, character(1L), USE.NAMES = FALSE)
      return(sprintf(
        "<prior: %s(%s)>",
        nm[[2L]],
        paste(vals, collapse = ", ")
      ))
    }
    return(sprintf("<prior: %s>", nm[[2L]]))
  }
  "<prior: expression>"
}

# .fb_compress_hyper_blobs() --- swap each hyper block for its tag.
#
# The blocks are found by scanning rather than by regular expression: a
# prior expression contains its own parentheses and quotes, so a balanced
# match has to know which of them are inside a string literal.
#
# @noRd
# @keywords internal
.fb_compress_hyper_blobs <- function(txt) {
  marker <- "hyper = list("
  repeat {
    start <- regexpr(marker, txt, fixed = TRUE)
    if (start < 0L) {
      return(txt)
    }
    chars <- strsplit(txt, "", fixed = TRUE)[[1L]]
    open_at <- start + nchar(marker) - 1L
    depth <- 0L
    in_string <- FALSE
    end_at <- NA_integer_
    for (i in seq(open_at, length(chars))) {
      ch <- chars[[i]]
      if (in_string) {
        if (ch == "\"" && !identical(chars[[max(i - 1L, 1L)]], "\\")) {
          in_string <- FALSE
        }
        next
      }
      if (ch == "\"") {
        in_string <- TRUE
      } else if (ch == "(") {
        depth <- depth + 1L
      } else if (ch == ")") {
        depth <- depth - 1L
        if (depth == 0L) {
          end_at <- i
          break
        }
      }
    }
    if (is.na(end_at)) {
      return(txt)
    }
    body <- paste(chars[seq(open_at, end_at)], collapse = "")
    txt <- paste0(
      paste(chars[seq_len(start - 1L)], collapse = ""),
      .fb_hyper_prior_tag(body),
      paste(chars[seq(end_at + 1L, length(chars))], collapse = "")
    )
  }
}

# .fb_display_inla_formula() --- the formula as a reader should see it.
#
# @noRd
# @keywords internal
.fb_display_inla_formula <- function(form) {
  if (is.null(form)) {
    return("(none recorded)")
  }
  txt <- paste(deparse(form), collapse = " ")
  txt <- gsub("[[:space:]]+", " ", txt)
  .fb_compress_hyper_blobs(txt)
}

# The spelling a user wrote for an AR1 field term, rebuilt from the IR so a
# refusal quotes the formula back rather than a canonical form the user never
# typed. ar1(env):geno and ar1(row):ar1(col) parse to the same IR type and
# differ only in `col_ar1`, so the second dimension is quoted plain when it
# carries no correlation of its own.
.ar1_term_spelling <- function(term) {
  if (identical(term$type, "ar1")) {
    return(paste0("ar1(", term$var, ")"))
  }
  inner <- if (isTRUE(term$col_ar1)) {
    paste0("ar1(", term$col_var, ")")
  } else {
    term$col_var
  }
  base <- paste0("ar1(", term$row_var, "):", inner)
  # C5/FS-27: the per-trial spelling wraps the single-field form in
  # at(<at_var>): -- quoting it back is what makes a refusal on the
  # per-trial shape name the formula the user actually wrote rather
  # than the single-field one.
  if (!is.null(term$at_var)) {
    return(paste0("at(", term$at_var, "):", base))
  }
  base
}

# WP16: build the INLA f() latent-field term for an AR1 / separable
# AR1xAR1 spatial term. A 1D ar1(t) maps to f(<t>_id, model = "ar1"); a
# separable ar1(row):ar1(col) maps to INLA's grouped-AR1 idiom --
# f(<row>_id, model = "ar1", group = <col>_id,
#   control.group = list(model = "ar1" | "iid")) -- whose precision is the
# Kronecker AR1(row) (x) AR1(col) (col_ar1 = FALSE gives AR1 x iid). The
# idiom's faithfulness was validated against a hand-built AR1(x)AR1 GLS/REML
# oracle (rebuild/wp16_inla_spatial_oracle.R). The <..>_id integer index
# columns are pre-built in emit_inla()'s data_inla setup.
#
# C5/FS-27: at(trial):ar1(row):ar1(col) -- one separable field per
# level of `trial` -- adds `replicate = <at_var>_id` to the same f().
# INLA's replicate mechanism fits K independent realisations of the
# identical model (same node count, same hyperparameters shared across
# every replicate) -- exactly "one field per level, shared
# hyperparameters", with no second f() call and no per-replicate
# hyperparameter block to keep in sync. `<at_var>_id` is pre-built in
# emit_inla()'s data_inla setup alongside the row / col index columns.
.inla_ar1_field_term <- function(term) {
  if (identical(term$type, "ar1")) {
    return(paste0("f(", term$var, "_id, model = \"ar1\")"))
  }
  col_model <- if (isTRUE(term$col_ar1)) "ar1" else "iid"
  replicate_clause <- if (!is.null(term$at_var)) {
    paste0(", replicate = ", term$at_var, "_id")
  } else {
    ""
  }
  paste0(
    "f(",
    term$row_var,
    "_id, model = \"ar1\", group = ",
    term$col_var,
    "_id, control.group = list(model = \"",
    col_model,
    "\")",
    replicate_clause,
    ")"
  )
}

# Build the control.family list passed to INLA::inla(). When the user
# supplies an fb_prior() spec keyed under "sigma", attach it as the
# residual-scale hyperprior under the keyword the likelihood actually
# declares -- `prec` for gaussian, lognormal, logistic, t and gamma, and
# `phi` for beta. The keyword roster lives in
# `.fb_inla_residual_hyper()` (R/family_traits.R) and is read from
# INLA's own `inla.models()` declaration; writing `prec` for beta is what
# made a beta fit on INLA fail with a raw engine error rather than fit.
#
# A likelihood with no residual-scale hyperparameter -- poisson,
# binomial, and the overdispersion-parameterised nbinomial and
# betabinomial -- takes no hyper input at all, so the sigma prior is
# dropped, as it is on the brms side for a family with no sigma. Returns
# an empty list when no residual prior is specified, so INLA's own
# default applies.
.build_inla_control_family <- function(hyper_ctrl, family) {
  entry <- hyper_ctrl[["sigma"]]
  if (is.null(entry)) {
    return(list())
  }
  keyword <- .fb_inla_hyper_keyword(family)
  if (is.null(keyword)) {
    return(list())
  }
  body <- entry[c("prior", "param")]
  body <- body[!vapply(body, is.null, logical(1))]
  stats::setNames(list(stats::setNames(list(body), keyword)), "hyper")
}

# Build the control.fixed list passed to INLA::inla(). INLA states the
# fixed-effect prior as a precision, so the documented SD maps to
# prec = 1 / sd^2, applied to the slopes and to the intercept alike --
# the argument documents itself as covering the intercept, the factor
# contrasts and the continuous slopes uniformly, and INLA's own defaults
# treat the intercept separately (prec.intercept = 0, i.e. flat).
#
# Empty unless the caller wrote `prior_fixed_sd` or supplied a `b()`
# prior row. Applying the documented default of 100 unconditionally would
# replace INLA's flat intercept prior on every fit, which is a different
# model for a response far from zero; supplied means honoured, and
# prior_summary() names the engine default otherwise.
#
# The two routes compose the way INLA composes them: the scalar sets the
# blanket `mean` / `prec`, a `b()` row sets a per-coefficient entry, and
# INLA resolves the named entry over the `default` one. A `b()` row alone
# leaves the unnamed coefficients on INLA's own defaults, which is the
# same supplied-ness rule the scalar follows.
.build_inla_control_fixed <- function(
  fb,
  prior_fixed_sd,
  data = NULL,
  coef_names = NULL
) {
  scalar_supplied <- .fb_prior_scalar_supplied(fb, "fixed_sd")
  per_coef <- .priors_to_inla_control_fixed(
    fb$priors,
    available = coef_names %||% .fb_inla_fixed_coef_names(fb, data)
  )
  if (!scalar_supplied) {
    return(per_coef)
  }
  if (
    !is.numeric(prior_fixed_sd) ||
      length(prior_fixed_sd) != 1L ||
      !is.finite(prior_fixed_sd) ||
      prior_fixed_sd <= 0
  ) {
    stop(.fb_refusal_condition(
      reason_code = "prior_hyperparameter_out_of_domain",
      message = paste0(
        "`prior_fixed_sd` must be a single positive number -- it is the ",
        "standard deviation of the fixed-effect normal prior. Got ",
        paste(format(prior_fixed_sd), collapse = ", "),
        "."
      )
    ))
  }
  prec <- 1 / (prior_fixed_sd^2)
  out <- list(
    mean = 0,
    prec = prec,
    mean.intercept = 0,
    prec.intercept = prec
  )
  if (!length(per_coef)) {
    return(out)
  }
  # Per-coefficient rows win over the blanket scalar, expressed the way
  # INLA reads it: a named list whose `default` entry is the scalar.
  if (!is.null(per_coef$mean)) {
    out$mean <- c(per_coef$mean, list(default = 0))
    out$prec <- c(per_coef$prec, list(default = prec))
  }
  for (nm in c("mean.intercept", "prec.intercept")) {
    if (!is.null(per_coef[[nm]])) {
      out[[nm]] <- per_coef[[nm]]
    }
  }
  out
}

# The fixed-effect coefficient names INLA will build from this model:
# the design-matrix columns of the fixed part, minus the intercept.
# Used to check a `b()` prior row against the model the fit will run
# (field-sweep FS-20). NULL when the names cannot be derived, which the
# caller reads as "no check possible" rather than "no coefficients".
.fb_inla_fixed_coef_names <- function(fb, data) {
  if (is.null(data) || !is.data.frame(data)) {
    return(NULL)
  }
  labels <- vapply(
    fb$fixed_terms %||% list(),
    function(term) {
      term$label %||%
        term$var %||%
        paste(term$vars %||% character(0), collapse = ":")
    },
    character(1)
  )
  labels <- labels[nzchar(labels)]
  if (!length(labels)) {
    return(character(0))
  }
  form <- tryCatch(
    stats::as.formula(paste0("~ ", paste(labels, collapse = " + "))),
    error = function(e) NULL
  )
  if (is.null(form)) {
    return(NULL)
  }
  cols <- tryCatch(
    colnames(stats::model.matrix(form, data = data)),
    error = function(e) NULL
  )
  if (is.null(cols)) {
    return(NULL)
  }
  setdiff(cols, "(Intercept)")
}


# Refuse an `sd(group = ...)` / `smooth("var")` prior whose key names no
# variance component of this model.
#
# The reachable set is the one the default-prior walker uses
# (`.fb_default_prior_targets()`), widened by every random term's own
# variable so a term outside the walker is not refused for being outside
# the walker. `sigma` is always reachable: the residual scale is dropped
# by family, not by term, and that drop is deliberate and documented in
# R/family_traits.R.
.check_inla_hyper_keys_reachable <- function(fb, hyper_ctrl) {
  keys <- setdiff(names(hyper_ctrl), "sigma")
  if (!length(keys)) {
    return(invisible(TRUE))
  }
  targets <- .fb_default_prior_targets(fb)
  term_vars <- unlist(lapply(
    fb$random_terms %||% list(),
    function(term) {
      c(term$var, term$inner, term$outer, term$slope_var)
    }
  ))
  available <- unique(c(targets$shared, targets$vm_ped, term_vars))
  available <- available[nzchar(available %||% character(0))]
  missing_keys <- setdiff(keys, available)
  if (!length(missing_keys)) {
    return(invisible(TRUE))
  }
  .fb_stop_prior_target_absent(
    target_label = paste0("sd(group = \"", missing_keys[[1L]], "\")"),
    kind = "variance component",
    available = available,
    engine = "inla"
  )
}


# ---------------------------------------------------------------- #
# Duplicate-key guard (field-sweep FS-17)                           #
# ---------------------------------------------------------------- #

# Refuse a model whose emitted INLA formula would key the same variable
# twice -- once as a fixed term and once as the index of an f() term.
#
# INLA's own message for this is "Key [x] is used twice and that is not
# allowed. A typical example where this happens is: y ~ x + f(x). Change
# this formula into: y ~ x + f(x2) where you define x2 = x", but it
# reaches the caller only as the generic "the inla-program exited with an
# error" wrapper, which names neither the variable nor the remedy. Both
# of INLA's remedies are repeated here, plus the one that is specific to
# a spline: an rw2 smooth already carries a linear trend in its null
# space, so the fixed copy of the same variable is usually redundant
# rather than needed.
.check_inla_variable_not_used_twice <- function(fb) {
  fixed_keys <- unlist(lapply(
    fb$fixed_terms %||% list(),
    function(term) {
      switch(
        term$type %||% "",
        "factor" = ,
        "continuous" = term$var,
        character(0)
      )
    }
  ))
  latent_keys <- unlist(lapply(
    fb$random_terms %||% list(),
    function(term) {
      switch(
        term$type %||% "",
        "simple" = ,
        "ide" = ,
        "id" = ,
        "spline" = term$var,
        character(0)
      )
    }
  ))
  clash <- intersect(
    fixed_keys %||% character(0),
    latent_keys %||% character(0)
  )
  if (!length(clash)) {
    return(invisible(TRUE))
  }
  is_spline <- any(vapply(
    fb$random_terms %||% list(),
    function(term) {
      identical(term$type, "spline") && term$var %in% clash
    },
    logical(1)
  ))
  stop(.fb_refusal_condition(
    reason_code = "inla_variable_used_twice",
    message = paste0(
      "backend = \"inla\": ",
      paste(paste0("`", clash, "`"), collapse = ", "),
      " would index both a fixed-effect term and a\n",
      "latent f() term in the emitted formula, and INLA refuses a key ",
      "used twice.\n\n",
      if (is_spline) {
        paste0(
          "Two remedies. Drop the fixed copy -- `y ~ 1, random = ~ spl(",
          clash[[1L]],
          ")` -- because\n  an rw2 smooth already carries a linear ",
          "trend in its null space, so the\n  fixed term is redundant ",
          "rather than additional. Or duplicate the column\n  (`",
          clash[[1L]],
          "2 <- ",
          clash[[1L]],
          "`) and smooth the copy, which is what INLA's own ",
          "message\n  suggests and keeps the two effects separately ",
          "identified.\n"
        )
      } else {
        paste0(
          "Two remedies. Drop the term from one side, or duplicate the ",
          "column (`",
          clash[[1L]],
          "2 <- ",
          clash[[1L]],
          "`)\n  and use the copy on the random side, which is what ",
          "INLA's own message suggests.\n"
        )
      },
      "backend = \"brms\" carries the same model without the ",
      "duplicate-key restriction."
    ),
    family_class = "flexybayes_inla_emit_refusal",
    variables = clash
  ))
}

# ---------------------------------------------------------------- #
# fb$family -> INLA family name                                    #
# ---------------------------------------------------------------- #

.resolve_inla_family <- function(fb) {
  fam <- if (inherits(fb$family, "family")) {
    fb$family$family
  } else {
    as.character(fb$family)
  }
  switch(
    tolower(fam),
    "gaussian" = "gaussian",
    "stdnormal" = "stdnormal",
    "binomial" = "binomial",
    "binary" = "binomial",
    "poisson" = "poisson",
    "negative_binomial" = "nbinomial",
    "negbinom" = "nbinomial",
    "nbinomial" = "nbinomial",
    "gamma" = "gamma",
    "beta" = "beta",
    "lognormal" = "lognormal",
    "exponential" = "exponential",
    fam
  )
}

# ---------------------------------------------------------------- #
# Post-fit numerical-confirm gate                                  #
# ---------------------------------------------------------------- #

.lgm_check_numerical <- function(fit) {
  pass <- TRUE
  reasons <- character(0)

  if (
    !is.null(fit$mode$mode.status) &&
      !isTRUE(fit$mode$mode.status == 0)
  ) {
    pass <- FALSE
    reasons <- c(
      reasons,
      paste0(
        "INLA mode.status = ",
        fit$mode$mode.status,
        " (non-zero indicates non-convergence)"
      )
    )
  }

  if (!is.null(fit$mlik) && nrow(fit$mlik) > 0L) {
    if (!is.finite(fit$mlik[1, 1])) {
      pass <- FALSE
      reasons <- c(reasons, "marginal log-likelihood is non-finite")
    }
  }

  list(pass = pass, reasons = reasons)
}

# ---------------------------------------------------------------- #
# INLA verification gate for (x || g)                             #
# ---------------------------------------------------------------- #

# Predicate: is any random term in the IR an uncorrelated random
# slope? Drives the per-fit verification check at the top of
# emit_inla().
.has_simple_slope_uncor <- function(fb) {
  if (is.null(fb$random_terms) || length(fb$random_terms) == 0L) {
    return(FALSE)
  }
  for (term in fb$random_terms) {
    if (identical(term$type, "simple_slope_uncor")) return(TRUE)
  }
  FALSE
}

# Refuse the (x || g) INLA mapping, which is deferred on every host.
#
# The three-arbitrator verification named a since-withdrawn engine as
# one arbitrator, so the criterion cannot be re-run as designed.
# The artefact at inst/extdata/inla-verification/simple_slope_uncor.rds
# survives as a developer rehearsal hook: it is excluded from the build
# and is consulted only when a developer sets
# options(flexyBayes.dev_inla_verification_artefacts = TRUE) by hand, so
# shipped behaviour is the same whether or not a host carries one.
#
# The refusal names backend = "brms" as the workaround (auto already
# falls back to it) and raises a structured condition of class
# flexybayes_inla_simple_slope_uncor_deferred so downstream tooling can
# pattern-match.
.check_inla_verification_simple_slope_uncor <- function() {
  pass <- FALSE
  if (.fb_dev_verification_artefacts_enabled()) {
    artefact_path <- system.file(
      "extdata",
      "inla-verification",
      "simple_slope_uncor.rds",
      package = "flexyBayes"
    )
    if (nzchar(artefact_path) && file.exists(artefact_path)) {
      art <- tryCatch(readRDS(artefact_path), error = function(e) NULL)
      if (is.list(art) && isTRUE(art$pass)) pass <- TRUE
    }
  }
  if (isTRUE(pass)) {
    return(invisible(TRUE))
  }

  msg <- paste0(
    "INLA mapping for uncorrelated random slopes (x || g) is ",
    "deferred to a future release.\n",
    "The INLA mapper is registered ",
    "only when the\nthree-arbitrator verification test (INLA vs an ",
    "engine since withdrawn vs lme4 on a simple fixture\nat J = 20 ",
    "groups) passes within the Wasserstein-1 \u2264 0.20 * tau_true\n",
    "tolerance on both sd_<g> and sd_<x>_<g>. One of the three ",
    "arbitrators has since been withdrawn as a fitting engine\n",
    "(see NEWS.md), so this INLA-native mapping stays deferred until ",
    "the\n",
    "verification is re-designed around active backends. The refusal ",
    "is\n",
    "host-independent -- no on-disk artefact lifts it.\n\n",
    "Workaround: fit via brms instead, which already represents ",
    "(x || g)\n",
    "natively. For example:\n",
    "  fit <- flexybayes(y ~ x + (x || g), data = d, ",
    "backend = \"brms\")\n",
    "  fit <- fb_brms   (y ~ x + (x || g), data = d)\n\n",
    "backend = \"auto\" already does this for you: an INLA-eligible ",
    "model\n",
    "that hits this deferral falls back to brms automatically.\n\n",
    "This refusal is mandatory: no silent translation of an ",
    "unverified mapping ships."
  )
  cond <- structure(
    class = c(
      "flexybayes_inla_simple_slope_uncor_deferred",
      "error",
      "condition"
    ),
    list(
      message = msg,
      call = NULL,
      deferral_target = "a future release",
      workaround = "backend = \"brms\""
    )
  )
  stop(cond)
}

# ---------------------------------------------------------------- #
# Separable AR1 field: the four parameters, named                  #
# ---------------------------------------------------------------- #
#
# INLA reports an autoregressive field on the precision scale and names its
# correlations after the integer index column the emit built, so the raw
# hyperparameter table reads `Precision for row_id`, `Rho for row_id`,
# `GroupRho for row_id` -- three lines that do not say which correlation
# runs along rows, which along columns, or which variance belongs to the
# field rather than to the observation nugget.
#
# All four are the model. The whole difference between what flexyBayes fits
# (a latent field plus an independent nugget) and ASReml's separable
# residual (one correlated residual, no nugget) is the fourth parameter, so
# a spatial fit prints all four under names that carry their meaning. A user
# who acts on the correlations can then see what they are acting on.

# ---------------------------------------------------------------- #
# Variance components, on the scale a reader thinks in               #
# ---------------------------------------------------------------- #
#
# INLA reports its hyperparameters as PRECISIONS. A standard deviation is
# 1 / sqrt(precision), a nonlinear and decreasing transform, so the
# posterior mean of the SD is not 1 / sqrt() of the posterior mean
# precision -- and the gap grows with the posterior spread, which is
# exactly where a variance component matters. Every SD-scale row below is
# therefore computed from the precision MARGINAL:
#
#   estimate   inla.emarginal(1 / sqrt(x), m)   -- posterior mean of the
#                                                  SD, integrated over
#                                                  the marginal
#   std.error  sqrt(E[1/x] - E[1/sqrt(x)]^2)    -- posterior SD of the SD
#   bounds     1 / sqrt(inla.qmarginal(1 - p))  -- the quantile identity
#
# The bounds use the identity rather than a numerically transformed
# density because a strictly monotone transform carries quantiles
# EXACTLY: the 2.5% quantile of the SD is 1 / sqrt() of the 97.5%
# quantile of the precision, with the order reversed because the
# transform is decreasing. This is not the forbidden operation, which is
# transforming a posterior MEAN -- a mean is not equivariant under a
# nonlinear transform and a quantile is. It is also the route the
# spatial-field table has taken since 0.9.0.
#
# Building the transformed density with INLA::inla.tmarginal() and
# reading quantiles off that was tried first and rejected on evidence:
# it evaluates the transform at inla.qmarginal((1:2048)/2049, m), whose
# extreme-left values fall slightly below the marginal's own support on a
# wide-support precision marginal -- 17 non-positive precisions on an
# ordinary random-intercept fit -- so 1 / sqrt() returns NaN and the call
# errors out. Where it does succeed it agrees with the identity above to
# about 1e-3 relative, being the interpolated answer to the exact one.
#
# Correlations (`Rho`, `GroupRho`) are not transformed: they are already
# on their own scale.
#
# The table is written in the five-column shape the brms path writes
# (component, estimate, sd, q2.5, q97.5) because broom's tidier reads
# those names (R/tidiers.R). The canonical component names come from the
# same registry mapper triangulate() and the draws accessors use, so one
# component answers to one name everywhere.
.inla_variance_comps <- function(object) {
  empty <- data.frame(
    component = character(0),
    estimate = numeric(0),
    sd = numeric(0),
    q2.5 = numeric(0),
    q97.5 = numeric(0),
    stringsAsFactors = FALSE
  )
  fit <- object$inla
  hp <- fit$summary.hyperpar
  if (is.null(hp) || nrow(hp) == 0L) {
    return(empty)
  }

  mapper <- tryCatch(
    .mapper_inla(object, object$fb),
    error = function(e) list(map = character(0))
  )
  marginals <- fit$marginals.hyperpar %||% list()
  need <- c("mean", "sd", "0.025quant", "0.975quant")
  have_cols <- all(need %in% colnames(hp))

  rows <- list()
  medians <- numeric(0)
  unavailable <- character(0)

  for (nm in rownames(hp)) {
    canonical <- unname(mapper$map[nm] %||% NA_character_)
    if (is.na(canonical)) {
      canonical <- nm
    }
    is_precision <- startsWith(nm, "Precision for ")

    if (is_precision) {
      m <- marginals[[nm]]
      if (is.null(m)) {
        # No marginal, no honest SD-scale summary. Transforming the
        # tabulated point estimate would produce a number that looks
        # computed and is not, so the row is dropped and named below.
        unavailable <- c(unavailable, nm)
        next
      }
      est <- INLA::inla.emarginal(function(x) 1 / sqrt(x), m)
      second <- INLA::inla.emarginal(function(x) 1 / x, m)
      sdv <- sqrt(max(second - est^2, 0))
      # Reversed probabilities: 1 / sqrt() is decreasing, so the SD's
      # lower bound is carried by the precision's upper quantile.
      #
      # Clamped to the marginal's own support. inla.qmarginal()
      # interpolates the quantile function, and on a precision whose
      # density is piled near zero the interpolation can return a value
      # slightly outside the grid it was built from -- occasionally a
      # negative precision, which has no square root. The clamp reads
      # the bound off the edge of the support instead, which is the
      # widest value the marginal actually carries; it never invents a
      # bound the density does not reach.
      support <- range(m[, 1L])
      prec_q <- INLA::inla.qmarginal(c(0.975, 0.5, 0.025), m)
      prec_q <- pmin(pmax(prec_q, support[[1L]]), support[[2L]])
      if (
        !all(is.finite(c(est, sdv, prec_q, support))) ||
          any(prec_q <= 0)
      ) {
        unavailable <- c(unavailable, nm)
        next
      }
      sd_q <- 1 / sqrt(prec_q)
      rows[[length(rows) + 1L]] <- data.frame(
        component = canonical,
        estimate = est,
        sd = sdv,
        q2.5 = sd_q[[1L]],
        q97.5 = sd_q[[3L]],
        stringsAsFactors = FALSE
      )
      medians[[canonical]] <- sd_q[[2L]]
      next
    }

    if (!have_cols) {
      unavailable <- c(unavailable, nm)
      next
    }
    rows[[length(rows) + 1L]] <- data.frame(
      component = canonical,
      estimate = hp[nm, "mean"],
      sd = hp[nm, "sd"],
      q2.5 = hp[nm, "0.025quant"],
      q97.5 = hp[nm, "0.975quant"],
      stringsAsFactors = FALSE
    )
    if ("0.5quant" %in% colnames(hp)) {
      medians[[canonical]] <- hp[nm, "0.5quant"]
    }
  }

  if (length(unavailable) > 0L) {
    warning(
      "flexyBayes: the hyperparameter(s) ",
      paste(unavailable, collapse = ", "),
      " carry no usable posterior marginal, so no variance-component row ",
      "was written for them and summary() will not report them. The fit ",
      "itself is unaffected. Check that control.compute = ",
      "list(return.marginals = TRUE), which is the flexyBayes default, ",
      "reached the engine.",
      call. = FALSE
    )
  }
  if (length(rows) == 0L) {
    return(empty)
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  # Posterior medians travel as an attribute rather than a sixth column:
  # the five column names are a broom contract. The boundary-collapse
  # display flag in summary() reads them.
  attr(out, "posterior_median") <- medians
  out
}


# ---------------------------------------------------------------- #
# The spatial-field collapse warning                                #
# ---------------------------------------------------------------- #
#
# An autoregressive field fitted on an incomplete grid can lose its
# signal entirely: the field standard deviation runs to a floor, both
# correlations sit at approximately zero with credible intervals spanning
# almost the whole of [-1, 1], and the variance the field should have
# carried is absorbed into the nugget. The fit's own convergence block
# reports a converged mode and a passing numerical confirm throughout,
# and the augmentation record says it completed the design, so nothing on
# the printed surface says the field was not identified. Measured on ten
# hole patterns over one 12 x 10 grid at 10 per cent missingness, five
# lost the field.
#
# The intervals are the only honest tell and they are easy to read past,
# so the fit says it out loud instead.
#
# Scope is the FIELD, on purpose. A grouping factor whose variance
# component collapses is a different fact about a different structure --
# often the thing a hierarchical model is there to show -- and it keeps
# the quieter `collapsed` cell in the variance-component table.

# A field standard deviation whose upper credible bound sits below this
# fraction of the nugget's posterior median carries no spatial signal at
# the scale of the noise it sits beside. Looser than the display flag's
# .FB_COLLAPSE_FRACTION because a collapsed field is compared against the
# nugget that absorbed its variance, not against an independent residual.
.FB_SPATIAL_COLLAPSE_FRACTION <- 0.05

# A correlation whose credible interval reaches beyond +/- this value at
# both ends spans nearly the whole parameter space: the data placed no
# constraint on it at all.
.FB_SPATIAL_RHO_UNIDENTIFIED <- 0.9

# .fb_spatial_collapse_reasons() --- which field parameters lost their
# identification, as a character vector of phrases. Empty when the field
# is healthy or when the fit carries no field.
#
# Read off the canonical variance-component table rather than INLA's own
# hyperparameter matrix, so the warning and the `collapsed` cell in
# summary(fit)$varcomp are readings of one table.
#
# @noRd
# @keywords internal
.fb_spatial_collapse_reasons <- function(object) {
  if (length(.inla_spatial_field_terms(object)) == 0L) {
    return(character(0L))
  }
  vc <- object$extras$variance_comps
  if (!is.data.frame(vc) || nrow(vc) == 0L) {
    return(character(0L))
  }
  medians <- attr(vc, "posterior_median") %||% numeric(0)
  cmp <- as.character(vc$component)
  reasons <- character(0L)

  i_field <- match("sd_spatial", cmp)
  nugget <- if ("sigma" %in% names(medians)) medians[["sigma"]] else NA_real_
  if (
    !is.na(i_field) &&
      is.numeric(nugget) &&
      length(nugget) == 1L &&
      is.finite(nugget) &&
      nugget > 0 &&
      is.finite(vc$q97.5[[i_field]]) &&
      vc$q97.5[[i_field]] < .FB_SPATIAL_COLLAPSE_FRACTION * nugget
  ) {
    reasons <- c(
      reasons,
      sprintf(
        paste0(
          "the field standard deviation is at the boundary (upper credible ",
          "bound %.4g against a nugget of %.4g)"
        ),
        vc$q97.5[[i_field]],
        nugget
      )
    )
  }

  rho_rows <- which(startsWith(cmp, "rho"))
  for (i in rho_rows) {
    lo <- vc$q2.5[[i]]
    hi <- vc$q97.5[[i]]
    if (
      is.finite(lo) &&
        is.finite(hi) &&
        lo < -.FB_SPATIAL_RHO_UNIDENTIFIED &&
        hi > .FB_SPATIAL_RHO_UNIDENTIFIED
    ) {
      reasons <- c(
        reasons,
        sprintf(
          "%s is unidentified (interval [%.3f, %.3f])",
          cmp[[i]],
          lo,
          hi
        )
      )
    }
  }
  reasons
}

# .fb_warn_spatial_field_collapsed() --- say it at fit time.
#
# @noRd
# @keywords internal
.fb_warn_spatial_field_collapsed <- function(object) {
  if (isTRUE(getOption("flexyBayes.silence_spatial_collapse_warning", FALSE))) {
    return(invisible(NULL))
  }
  reasons <- .fb_spatial_collapse_reasons(object)
  if (length(reasons) == 0L) {
    return(invisible(NULL))
  }
  warning(
    "flexyBayes: the autoregressive field in this fit did not identify -- ",
    paste(reasons, collapse = ", "),
    ". The field carries no spatial signal and the fit is effectively an ",
    "independent-errors model with the field's variance absorbed into the ",
    "residual; the convergence block reports a converged mode regardless, ",
    "because the optimiser did converge -- to a solution with no field in ",
    "it. An incomplete grid raises the risk: in repeated runs over one ",
    "12 x 10 field at 10 per cent missing responses, half the hole ",
    "patterns lost the field and half did not, on the same data. Three ",
    "routes: complete the grid, either by supplying the field-book rows ",
    "for every sown plot with the response set to NA or by padding the ",
    "array with fb_complete_grid(); fit the same model on the brms engine ",
    "for a second reading; or give the field an informative prior with ",
    "fb_prior() so it is not left to run to a boundary. Silence via ",
    "options(flexyBayes.silence_spatial_collapse_warning = TRUE).",
    call. = FALSE
  )
  invisible(NULL)
}


# The autoregressive field terms carried by a fitted IR.
.inla_spatial_field_terms <- function(object) {
  Filter(
    function(t) (t$type %||% "") %in% c("ar1", "ar1_spatial"),
    object$fb$random_terms %||% list()
  )
}

# One labelled row per field parameter, on the correlation and standard-
# deviation scales.
#
# The point estimate is the posterior MEDIAN, not the mean, because the
# standard deviations are read off INLA's precision marginals and only a
# monotone summary survives the 1 / sqrt() transform exactly. The interval
# is INLA's own 2.5% / 97.5% quantiles, transformed (and reversed, since
# 1 / sqrt() is decreasing). Nothing here is a re-estimate.
.inla_spatial_hyper_table <- function(object) {
  field_terms <- .inla_spatial_field_terms(object)
  if (length(field_terms) == 0L) {
    return(NULL)
  }
  hp <- object$inla$summary.hyperpar
  if (is.null(hp) || nrow(hp) == 0L) {
    return(NULL)
  }
  need <- c("0.5quant", "0.025quant", "0.975quant")
  if (!all(need %in% colnames(hp))) {
    return(NULL)
  }

  as_corr <- function(row_name, label) {
    i <- match(row_name, rownames(hp))
    if (is.na(i)) {
      return(NULL)
    }
    data.frame(
      parameter = label,
      median = hp[i, "0.5quant"],
      lower = hp[i, "0.025quant"],
      upper = hp[i, "0.975quant"],
      stringsAsFactors = FALSE
    )
  }
  as_sd <- function(row_name, label) {
    i <- match(row_name, rownames(hp))
    if (is.na(i)) {
      return(NULL)
    }
    data.frame(
      parameter = label,
      median = 1 / sqrt(hp[i, "0.5quant"]),
      lower = 1 / sqrt(hp[i, "0.975quant"]),
      upper = 1 / sqrt(hp[i, "0.025quant"]),
      stringsAsFactors = FALSE
    )
  }

  out <- list()
  for (term in field_terms) {
    if (identical(term$type, "ar1")) {
      idx <- paste0(term$var, "_id")
      out <- c(
        out,
        list(
          as_corr(
            paste0("Rho for ", idx),
            paste0("correlation along ", term$var, " (rho)")
          ),
          as_sd(
            paste0("Precision for ", idx),
            "field SD (sigma_field)"
          )
        )
      )
    } else {
      idx <- paste0(term$row_var, "_id")
      out <- c(
        out,
        list(
          as_corr(
            paste0("Rho for ", idx),
            paste0(
              "correlation along ",
              term$row_var,
              " (rho_",
              term$row_var,
              ")"
            )
          ),
          as_corr(
            paste0("GroupRho for ", idx),
            paste0(
              "correlation along ",
              term$col_var,
              " (rho_",
              term$col_var,
              ")"
            )
          ),
          as_sd(
            paste0("Precision for ", idx),
            "field SD (sigma_field)"
          )
        )
      )
    }
  }
  out <- c(
    out,
    list(as_sd(
      "Precision for the Gaussian observations",
      "nugget SD (sigma_e)"
    ))
  )
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0L) {
    return(NULL)
  }
  tab <- do.call(rbind, out)
  rownames(tab) <- NULL
  tab
}

# INLA's own hyperparameter table, on the scale INLA reports it.
#
# The variance-component table above it is the same information on the
# standard-deviation scale, which is the scale a reader thinks in. This
# block stays because it is the engine's own answer, unmodified, and a
# user checking flexyBayes against a plain INLA::inla() run needs
# something to check against. A no-op on every other engine.
.print_inla_hyperpar_table <- function(object) {
  hp <- object$extras$summary$hyperpar
  if (is.null(object$inla) || is.null(hp) || nrow(hp) == 0L) {
    return(invisible(NULL))
  }
  cat(
    "\n-- Hyperparameters (INLA, precision scale) ",
    strrep("-", 21),
    "\n",
    sep = ""
  )
  print(round(hp, 4))
  invisible(hp)
}

# Render the field block. Shared by print() and summary() so the two cannot
# drift apart.
.print_inla_spatial_hypers <- function(object) {
  tab <- .inla_spatial_hyper_table(object)
  if (is.null(tab)) {
    return(invisible(NULL))
  }
  separable <- any(vapply(
    .inla_spatial_field_terms(object),
    function(t) identical(t$type, "ar1_spatial"),
    logical(1L)
  ))
  cat(
    "\n",
    if (separable) "Separable AR1 field" else "AR1 field",
    " (posterior median, 95% interval):\n",
    sep = ""
  )
  body <- data.frame(
    parameter = tab$parameter,
    median = round(tab$median, 4),
    `2.5%` = round(tab$lower, 4),
    `97.5%` = round(tab$upper, 4),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  print(body, row.names = FALSE)
  cat(
    "  The field and the nugget are separate parameters: this is a ",
    "latent\n  autoregressive field plus independent observation noise, ",
    "not ASReml's\n  nugget-free separable residual. Variances are the ",
    "squares of the SDs.\n",
    sep = ""
  )
  invisible(tab)
}


# ---------------------------------------------------------------- #
# Print + summary methods for the flexybayes_inla wrapper          #
# ---------------------------------------------------------------- #

#' Print method for flexybayes_inla
#'
#' Internal S3 method. Brief one-screen description of an INLA fit
#' produced via `fb(... backend = "inla")` or `emit_inla()`.
#'
#' Opens with the header every engine's print shares, so the three prints
#' cannot disagree about what the fit is or how many rows it saw, then
#' adds what belongs to this engine alone: the formula as it reached
#' `INLA::inla()`, the post-fit numerical-confirm verdict, and -- on a fit
#' carrying a latent autoregressive field -- the field's own parameters
#' on the correlation and standard-deviation scales.
#'
#' @param x   A `flexybayes_inla` object, the fit an INLA run returns.
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the one-screen summary
#'   it prints.
#' @keywords internal
#' @export
print.flexybayes_inla <- function(x, ...) {
  mi <- x$extras$model_info
  .fb_print_header(x, "Bayesian mixed model", "-")
  # The full formula, prior expressions and all, stays on the object at
  # `$extras$formula`; what prints is the same formula with each prior
  # block named rather than spelled out.
  cat(
    "  INLA formula: ",
    .fb_display_inla_formula(x$extras$formula),
    "\n",
    sep = ""
  )
  cat(
    "  Params   : ",
    mi$n_fixed,
    " fixed, ",
    mi$n_random,
    " random, ",
    mi$n_hyper,
    " hyperparameter(s)\n",
    sep = ""
  )
  cat(
    "  numerical confirm: ",
    if (isTRUE(x$num_check$pass)) "PASS" else "FAIL",
    if (!isTRUE(x$num_check$pass)) {
      paste0(" (", paste(x$num_check$reasons, collapse = "; "), ")")
    } else {
      ""
    },
    "\n",
    sep = ""
  )
  .print_inla_spatial_hypers(x)
  cat(strrep("-", 62), "\n")
  cat("  Object of class <flexybayes_inla>\n")
  cat("  $inla -- raw INLA fit (use INLA's summary, plot, etc.)\n")
  cat("  $fb   -- the fb_terms IR used for dispatch\n")
  invisible(x)
}

#' Summary method for flexybayes_inla
#'
#' Internal S3 method. Builds and prints the same `summary.flexybayes`
#' object every active engine returns, so `summary(fit)$varcomp` answers
#' on an INLA fit as it does on a brms one. Before 0.9.1 this method had
#' its own dialect -- a bare four-slot list of INLA's own tables, not
#' comparable with what the other engine returned and carrying no
#' variance-component table at all.
#'
#' The variance components reach the standard-deviation scale through
#' INLA's precision marginals rather than by transforming a tabulated
#' point estimate; see `.inla_variance_comps()`. INLA's own
#' precision-scale hyperparameter table still prints, unmodified,
#' beneath them.
#'
#' @param object A `flexybayes_inla` object, the fit an INLA run
#'   returns.
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, an object of class
#'   `c("summary.flexybayes", "list")`. See [summary.flexybayes()] for
#'   the slots.
#' @keywords internal
#' @export
summary.flexybayes_inla <- function(object, ...) {
  out <- .fb_summary_object(object)
  print(out)
  invisible(out)
}
