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
#     re-route via lgm_gate or backend = "greta".
#   - Spatial / GxE structured covariances -- v0.1 refused; route
#     via greta.
#
# Anything not in the supported set raises an INLA-side refusal
# pointing the user back to backend = "greta".
#
# The post-fit numerical-confirm gate lives here: assert
# mode.status == 0 and a finite marginal likelihood (mlik). Failures
# escalate to a structured warning.
#
# INLA is not on CRAN; ships via inla.r-inla-download.org. Added to
# DESCRIPTION:Suggests with Additional_repositories declaration.

# ---------------------------------------------------------------- #
# Public-style entry -- internal in v0.1                           #
# ---------------------------------------------------------------- #

emit_inla <- function(
  fb,
  data,
  known_matrices = list(),
  seed = NULL,
  control = NULL,
  verbose = TRUE,
  return_code = FALSE,
  the_call = NULL,
  fixed = NULL,
  random = NULL,
  residual = NULL,
  family = NULL,
  link = NULL,
  data_name = NA_character_,
  ...
) {
  .check_fb_terms(fb, "`fb` must be an fb_terms object.")

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
  } else {
    list()
  }
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
      control_family = control_family
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
      n_row <- length(unique(data[[term$row_var]]))
      n_col <- length(unique(data[[term$col_var]]))
      combos <- interaction(
        data[[term$row_var]], data[[term$col_var]],
        drop = TRUE
      )
      one_obs_per_node <- nrow(data) == n_row * n_col &&
        !any(duplicated(combos))
      if (!one_obs_per_node) {
        stop(.fb_refusal_condition(
          reason_code = "ar1_spatial_requires_complete_grid",
          message = paste0(
            .ar1_term_spelling(term), " on INLA is ",
            "emitted as a separable autoregressive latent field over the ",
            n_row, " x ", n_col, " grid, faithful only with one observation ",
            "per (", term$row_var, ", ", term$col_var, ") node. The data ",
            "have ", nrow(data), " rows for ", n_row * n_col, " nodes ",
            "(incomplete or replicated grid), so the field's latent index ",
            "would not correspond to the design. Keep the unobserved nodes ",
            "as design cells with na_action = \"augment\" (the default), ",
            "supply one observation per node, or aggregate to node means."
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
    cat("  formula: ", deparse(inla_form), "\n", sep = "")
    cat("  family:  ", inla_family, "\n", sep = "")
    cat(paste(rep("-", 60), collapse = ""), "\n\n")
  }

  # Fit
  t0 <- proc.time()
  fit <- tryCatch(
    INLA::inla(
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
    ),
    error = function(e) {
      stop("INLA fit failed: ", conditionMessage(e), call. = FALSE)
    }
  )
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
  structure(
    list(
      inla = fit,
      fb = fb,
      data = data,
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
        call_info = list(
          fixed = fixed,
          random = random,
          residual = residual,
          data_name = data_name,
          family = family,
          link = link
        ),
        # Thread the IR's smooth-objects slot through to the
        # extras so predict() on the INLA path can also use
        # mgcv::Predict.matrix() when smooths are present. The INLA
        # backend itself does not currently fit s() smooths via mgcv
        # (it uses INLA's own rw2 path), so this slot will normally be
        # an empty list; threaded for shape uniformity with emit_greta.
        parse_info = list(
          smooths = .collect_smooths(fb$random_terms),
          # Threaded for shape uniformity with emit_greta / emit_brms so
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
        fingerprint = .fb_model_fingerprint(fb, data)
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
      # posterior agreement with both greta and lme4 on a gaussian-
      # identity fixture before this branch became reachable.
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
  paste0("ar1(", term$row_var, "):", inner)
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
.inla_ar1_field_term <- function(term) {
  if (identical(term$type, "ar1")) {
    return(paste0("f(", term$var, "_id, model = \"ar1\")"))
  }
  col_model <- if (isTRUE(term$col_ar1)) "ar1" else "iid"
  paste0(
    "f(", term$row_var, "_id, model = \"ar1\", group = ",
    term$col_var, "_id, control.group = list(model = \"", col_model, "\"))"
  )
}

# Build the control.family list passed to INLA::inla(). When the
# user supplies an fb_prior() spec keyed under "sigma", attach it as
# the residual-precision hyperprior -- but only for families whose
# likelihood actually carries a precision hyperparameter (Gaussian,
# lognormal, gamma, beta, T). For Poisson, binomial, exponential
# etc. the likelihood has no `prec` hyperparameter and INLA refuses
# any control.family$hyper input, so we silently drop the sigma prior
# in that case (it is not a no-op for the model -- there is no such
# parameter to constrain). Returns an empty list when no residual
# prior is specified, so INLA's default loggamma applies for the
# Gaussian-family case where prec is the lone hyperparameter.
.build_inla_control_family <- function(hyper_ctrl, family) {
  entry <- hyper_ctrl[["sigma"]]
  if (is.null(entry)) {
    return(list())
  }
  families_with_prec <- c(
    "gaussian",
    "stdnormal",
    "lognormal",
    "gamma",
    "beta",
    "T",
    "logistic"
  )
  if (!tolower(family) %in% families_with_prec) {
    return(list())
  }
  body <- entry[c("prior", "param")]
  body <- body[!vapply(body, is.null, logical(1))]
  list(hyper = list(prec = body))
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
# The three-arbitrator verification named greta as one arbitrator and
# greta is quarantined, so the criterion cannot be re-run as designed.
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
    "only when the\nthree-arbitrator verification test (INLA vs ",
    "greta vs lme4 on a simple fixture\nat J = 20 groups) passes ",
    "within the Wasserstein-1 \u2264 0.20 * tau_true\n",
    "tolerance on both sd_<g> and sd_<x>_<g>. greta served as one of\n",
    "the three arbitrators and has since been withdrawn as a fitting ",
    "engine\n",
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
      out <- c(out, list(
        as_corr(
          paste0("Rho for ", idx),
          paste0("correlation along ", term$var, " (rho)")
        ),
        as_sd(
          paste0("Precision for ", idx),
          "field SD (sigma_field)"
        )
      ))
    } else {
      idx <- paste0(term$row_var, "_id")
      out <- c(out, list(
        as_corr(
          paste0("Rho for ", idx),
          paste0("correlation along ", term$row_var, " (rho_",
                 term$row_var, ")")
        ),
        as_corr(
          paste0("GroupRho for ", idx),
          paste0("correlation along ", term$col_var, " (rho_",
                 term$col_var, ")")
        ),
        as_sd(
          paste0("Precision for ", idx),
          "field SD (sigma_field)"
        )
      ))
    }
  }
  out <- c(out, list(as_sd(
    "Precision for the Gaussian observations",
    "nugget SD (sigma_e)"
  )))
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0L) {
    return(NULL)
  }
  tab <- do.call(rbind, out)
  rownames(tab) <- NULL
  tab
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
    "\n", if (separable) "Separable AR1 field" else "AR1 field",
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
#' Internal S3 method. Brief one-screen summary of an INLA fit
#' produced via `fb(... backend = "inla")` or `emit_inla()`.
#'
#' @param x   A `flexybayes_inla` object, the fit an INLA run returns.
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the one-screen summary
#'   it prints.
#' @keywords internal
#' @export
print.flexybayes_inla <- function(x, ...) {
  mi <- x$extras$model_info
  cat("Bayesian fit  [flexybayes_inla / INLA backend]\n")
  cat(strrep("-", 55), "\n")
  cat("  formula: ", deparse(x$extras$formula), "\n", sep = "")
  cat("  family:  ", mi$family, "\n", sep = "")
  cat("  n_obs:   ", mi$n_obs, "\n", sep = "")
  cat("  fixed:   ", mi$n_fixed, "\n", sep = "")
  cat("  random:  ", mi$n_random, "\n", sep = "")
  cat("  hyper:   ", mi$n_hyper, "\n", sep = "")
  cat("  runtime: ", round(x$extras$run_time, 2), " sec\n", sep = "")
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
  cat(strrep("-", 55), "\n")
  cat("  $inla -- raw INLA fit (use INLA's summary, plot, etc.)\n")
  cat("  $fb   -- the fb_terms IR used for dispatch\n")
  invisible(x)
}

#' Summary method for flexybayes_inla
#'
#' Returns the fixed / random / hyperpar posterior summary tables
#' produced by INLA. A fit carrying an autoregressive latent field also
#' gets a `spatial_field` element: the field's own parameters on the
#' correlation and standard-deviation scales, which is the form a reader
#' needs and INLA's precision-scale hyperparameter table is not.
#' Internal S3 method.
#'
#' @param object A `flexybayes_inla` object, the fit an INLA run
#'   returns.
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, the summary list, carrying a `spatial_field` data
#'   frame when the model holds an AR1 field.
#' @keywords internal
#' @export
summary.flexybayes_inla <- function(object, ...) {
  cat("Bayesian fit summary  [flexybayes_inla / INLA backend]\n")
  cat(strrep("-", 60), "\n")
  cat("Fixed effects:\n")
  if (
    !is.null(object$extras$summary$fixed) &&
      nrow(object$extras$summary$fixed) > 0L
  ) {
    print(round(object$extras$summary$fixed, 4))
  } else {
    cat("  (none)\n")
  }
  cat("\nHyperparameters:\n")
  if (
    !is.null(object$extras$summary$hyperpar) &&
      nrow(object$extras$summary$hyperpar) > 0L
  ) {
    print(round(object$extras$summary$hyperpar, 4))
  } else {
    cat("  (none)\n")
  }
  .print_inla_spatial_hypers(object)
  cat("\nRandom effects:\n")
  if (length(object$extras$summary$random) > 0L) {
    cat(
      "  groups: ",
      paste(names(object$extras$summary$random), collapse = ", "),
      "\n",
      sep = ""
    )
  } else {
    cat("  (none)\n")
  }
  cat(strrep("-", 60), "\n")
  out <- object$extras$summary
  field <- .inla_spatial_hyper_table(object)
  if (!is.null(field)) {
    out$spatial_field <- field
  }
  invisible(out)
}
