# emit_greta -- flexyBayes greta backend emit shim.
#
# A no-semantic-change shim that takes an fb_terms IR
# (intermediate representation) object and runs the full greta
# backend pipeline that flexybayes() previously ran inline. This
# lifts greta dispatch out of flexybayes() and behind a single
# emit_*() function, mirroring what emit_brms and emit_inla do.
#
# Behaviour invariant: for every input previously valid for
# flexybayes(), the migrated path flexybayes() -> fb_from_asreml() ->
# emit_greta() produces byte-identical generated greta code to the
# earlier flexybayes() inline pipeline. Snapshot-tested in
# tests/testthat/test-emit-greta.R; transparency-tested in
# tests/testthat/test-fb-from-asreml.R.
#
# Internal -- not exported. Called from flexybayes() and from
# triangulate() when greta is one of the requested backends.

emit_greta <- function(
  fb,
  data,
  known_matrices = list(),
  weights = NULL,
  n_samples = 1000,
  warmup = 500,
  chains = 4,
  prior_fixed_sd = 100,
  prior_vc_sd = 1,
  verbose = TRUE,
  mcmc_verbose = TRUE,
  return_code = FALSE,
  the_call = NULL,
  fixed = NULL,
  random = NULL,
  rcov = NULL,
  family = NULL,
  link = NULL,
  data_name = NA_character_
) {
  if (!is_fb_terms(fb)) {
    stop("`fb` must be an fb_terms object (see fb_from_asreml).", call. = FALSE)
  }

  # ---------------------------------------------------------------- #
  # Unwrap the IR back into the shape the existing helpers expect.   #
  # fixed_info, random_terms, rcov_terms,                            #
  # fam_link are byte-identical to what flexybayes() previously      #
  # produced inline via .parse_fixed / .parse_formula /              #
  # .resolve_family.                                                 #
  # ---------------------------------------------------------------- #

  fixed_info <- list(
    response = fb$response,
    intercept = fb$intercept,
    terms = fb$fixed_terms
  )
  random_terms <- fb$random_terms
  rcov_terms <- fb$rcov_terms
  fam_link <- list(family = fb$family, link = fb$link)

  # ---------------------------------------------------------------- #
  # Build evaluation environment (unchanged from earlier flexybayes).#
  # ---------------------------------------------------------------- #

  greta_ns <- asNamespace("greta")
  ev <- new.env(parent = greta_ns)
  .setup_env(
    ev,
    fixed_info,
    random_terms,
    rcov_terms,
    data,
    known_matrices,
    weights
  )

  # ---------------------------------------------------------------- #
  # Generate code (unchanged).                                       #
  # ---------------------------------------------------------------- #

  # If the user supplied an fb_prior(), derive per-variance-component
  # spec maps so codegen can emit greta::uniform(lower, upper) (the
  # v0.1 default) or greta::exponential(rate) (PC; explicit choice)
  # instead of the legacy lognormal default for the
  # matched groups.
  legacy_specs <- if (inherits(fb$priors, "fb_prior")) {
    priors_to_legacy(fb$priors)
  } else {
    list()
  }
  pc_per_vc <- legacy_specs$pc_per_vc %||% list()
  uniform_per_vc <- legacy_specs$uniform_per_vc %||% list()

  ctx <- list(
    params = character(0),
    predictor = character(0),
    code = character(0),
    fam_link = fam_link,
    prior_vc = prior_vc_sd,
    prior_fx = prior_fixed_sd,
    pc_per_vc = pc_per_vc,
    uniform_per_vc = uniform_per_vc,
    # Codegen for s() smooths binds the basis matrix into
    # this environment as B_s_<vname>; emit_greta evaluates the
    # generated code against `ev` (which is `ctx$env` here) so the
    # bound name is in scope. Other code-paths that build a ctx
    # without an env (notably return_code = TRUE branches that never
    # eval) must still work -- the assign() below is harmless on a
    # throw-away env.
    env = ev
  )

  ctx <- .code_fixed(ctx, fixed_info)
  ctx <- .code_random(ctx, random_terms, data, known_matrices)
  ctx <- .code_rcov(ctx, rcov_terms, data)
  ctx <- .code_predictor(ctx, fixed_info)
  ctx <- .code_likelihood(ctx, fixed_info, rcov_terms, data, weights)
  ctx <- .code_model(ctx, n_samples, warmup, chains, mcmc_verbose)

  code_str <- paste(ctx$code, collapse = "\n")

  if (verbose) {
    cat(
      "\n-- flexyBayes: generated greta code ",
      paste(rep("-", 40), collapse = ""),
      "\n",
      sep = ""
    )
    cat(code_str, "\n")
    cat(paste(rep("-", 60), collapse = ""), "\n\n")
  }

  if (return_code) {
    return(invisible(code_str))
  }

  # ---------------------------------------------------------------- #
  # Evaluate (unchanged).                                            #
  # ---------------------------------------------------------------- #

  t0 <- proc.time()
  tryCatch(
    eval(parse(text = code_str), envir = ev),
    error = function(e) {
      stop("flexyBayes evaluation error: ", conditionMessage(e))
    }
  )
  elapsed <- unname((proc.time() - t0)["elapsed"])

  # ---------------------------------------------------------------- #
  # Post-processing (unchanged).                                     #
  # ---------------------------------------------------------------- #

  all_obj_names <- ls(ev)
  is_greta <- vapply(
    all_obj_names,
    function(nm) {
      tryCatch(
        inherits(get(nm, envir = ev, inherits = FALSE), "greta_array"),
        error = function(e) FALSE
      )
    },
    logical(1)
  )
  greta_arrays <- if (any(is_greta)) {
    nms <- all_obj_names[is_greta]
    setNames(
      lapply(nms, function(nm) get(nm, envir = ev, inherits = FALSE)),
      nms
    )
  } else {
    list()
  }

  post_summary <- tryCatch(summary(ev$atg_draws), error = function(e) NULL)

  n_eff <- tryCatch(coda::effectiveSize(ev$atg_draws), error = function(e) NULL)
  gelman <- if (chains >= 2 && length(ev$atg_draws) >= 2) {
    tryCatch(
      coda::gelman.diag(ev$atg_draws, multivariate = FALSE),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  conv_diag <- list(n_eff = n_eff, gelman = gelman)

  # Three-part output assembly
  glm_obj <- .build_glm(
    fixed_info = fixed_info,
    data = data,
    draws = ev$atg_draws,
    ev = ev,
    fam_link = fam_link,
    the_call = the_call,
    fixed = fixed
  )

  greta_out <- structure(
    list(
      model = ev$atg_model,
      draws = ev$atg_draws,
      greta_arrays = greta_arrays,
      env = ev
    ),
    class = "flexybayes_greta"
  )

  vc_table <- .build_variance_comps(post_summary, ctx$params)
  blups <- .build_blups(post_summary, random_terms)
  preds <- .build_predictions(data, fixed_info, post_summary, ev, fam_link)

  extras <- structure(
    list(
      summary = post_summary,
      convergence = conv_diag,
      variance_comps = vc_table,
      blups = blups,
      predictions = preds,
      code = code_str,
      param_names = unique(ctx$params),
      parse_info = list(
        fixed = fixed_info,
        random = random_terms,
        rcov = rcov_terms,
        family = fam_link,
        # Collect non-NULL smooth_obj slots into a named list
        # keyed by smooth-term variable name. Empty list when the model
        # has no smooths. predict.flexybayes() reads this slot and
        # re-evaluates the basis on newdata via mgcv::Predict.matrix().
        smooths = .collect_smooths(random_terms),
        # Collect per-smooth
        # approximation metadata (scheme, rank, projection V_K, full
        # singular spectrum, realised Frobenius capture) keyed by
        # smooth variable. Empty list when no smooth was routed through
        # an approximation scheme. validate_approximation() reads this
        # slot for its verdict; predict.flexybayes() reads V_K to
        # project the newdata basis into the truncated coefficient
        # space (R/predict_kernel.R).
        approx = .collect_approx(random_terms),
        # Record which factor:continuous emit branch fired
        # (Option C primary or Option D fallback). NULL if no
        # factor_numeric_interaction term was present.
        factor_continuous_emit = ctx$factor_continuous_emit
      ),
      call_info = list(
        fixed = fixed,
        random = random,
        rcov = rcov,
        data_name = data_name,
        family = family,
        link = link,
        known_matrices = known_matrices,
        weights = weights,
        n_samples = n_samples,
        warmup = warmup,
        chains = chains,
        prior_fixed_sd = prior_fixed_sd,
        prior_vc_sd = prior_vc_sd
      ),
      run_time = elapsed,
      model_info = list(
        n_obs = nrow(data),
        n_fixed = length(fixed_info$terms) + fixed_info$intercept,
        n_random = length(random_terms),
        n_params = length(unique(ctx$params)),
        family = fam_link$family,
        link = fam_link$link
      )
    ),
    class = "flexybayes_extras"
  )

  result <- structure(
    list(
      glm = glm_obj,
      greta = greta_out,
      extras = extras
    ),
    class = "flexybayes"
  )

  invisible(result)
}

# Collect non-NULL `smooth_obj` slots from the IR's
# random_terms list into a named list keyed by the smooth-term
# variable name (e.g., `s(x)` -> `smooths[["x"]] <- <Smooth>`).
# Returns an empty named list when no smooths are present, so the
# downstream `parse_info$smooths` slot has a uniform shape.
# Internal helper; also reused by emit_inla() for the IR slot.
.collect_smooths <- function(random_terms) {
  out <- list()
  if (length(random_terms) == 0L) {
    return(out)
  }
  for (t in random_terms) {
    if (!is.null(t$smooth_obj)) {
      key <- t$var %||% t$smooth_label %||% paste0("s_", length(out) + 1L)
      out[[key]] <- t$smooth_obj
    }
  }
  out
}

# Collect the per-smooth
# approximation metadata recorded by .enrich() (R/parse_formula.R) into
# a named list keyed by smooth-term variable name. Each entry carries
# the scheme, the truncation rank, the projection V_K (the predict path
# projects the newdata basis through it), the full singular spectrum,
# and the realised Frobenius capture. Returns an empty named list when
# no smooth was routed through an approximation scheme, so the
# parse_info$approx slot has a uniform shape.
.collect_approx <- function(random_terms) {
  out <- list()
  if (length(random_terms) == 0L) {
    return(out)
  }
  for (t in random_terms) {
    if (!is.null(t$approx_spec)) {
      key <- t$var %||% t$smooth_label %||% paste0("s_", length(out) + 1L)
      out[[key]] <- t$approx_spec
    }
  }
  out
}
