# fb_greta() -- the greta engine pin, plus the native-greta fitting core.
#
# v0.5.0 backend-axis recovery. `fb_greta()` is one
# of the three engine pins (alongside `fb_inla()` and `fb_brms()`): sugar
# over the universal entry with the engine fixed --- `fb_greta(...)` ==
# `fb(..., backend = "greta")`. It accepts every grammar the universal
# entry does, including a native `greta_model` graph: a formula is lowered
# through the shared emit path, while a native graph is fit directly by
# `greta::mcmc()` via `.fit_native_greta()`.
#
# Before v0.5.0 `fb_greta(model = <greta_model>)` was a distinct
# native-graph verb. That signature is removed here; the
# native-graph fitting body it owned is extracted into the internal
# `.fit_native_greta()`, which the universal-entry dispatch
# (`.dispatch_native_greta()` in R/dispatch.R) now drives whenever the IR
# carries `source == "greta"`. The `model = ` argument is remapped to the
# universal entry's model-spec slot for call-compatibility, so
# `fb_greta(model = m)` and `fb_greta(m)` both still fit a native graph.

# .fit_native_greta() -- fit a native greta model graph by greta::mcmc()
# and assemble the classed flexybayes_direct_greta result. The IR (`fb`,
# a greta-source fb_terms built by fb_from_greta()) carries the model
# graph on `fb$greta_meta$model` plus the canonical-name map, so this
# helper is self-contained. Internal -- reached via the universal entry's
# native-greta dispatch, never called directly by users.
#
# This is the verbatim body of the pre-v0.5.0 fb_greta() native path:
# no emit_greta() lowering, a minimal GLM-compatible
# shim built from the posterior summary, and the
# c("flexybayes_direct_greta", "flexybayes", "list") class.
.fit_native_greta <- function(
  fb,
  n_samples = 1000L,
  warmup = 500L,
  chains = 4L,
  verbose = TRUE,
  mcmc_verbose = TRUE,
  the_call = NULL
) {
  if (!requireNamespace("greta", quietly = TRUE)) {
    stop(
      "Package 'greta' is required to fit a native greta model. ",
      "Install with:\n  install.packages('greta')",
      call. = FALSE
    )
  }

  model <- fb$greta_meta$model
  if (is.null(model) || !inherits(model, "greta_model")) {
    stop(
      ".fit_native_greta(): the IR does not carry a greta model graph ",
      "on `greta_meta$model`. This is a programming error -- the IR ",
      "must be built by fb_from_greta().",
      call. = FALSE
    )
  }

  # ------- 2. Direct greta::mcmc() call (no emit_greta lowering). - #
  run_time_start <- proc.time()[["elapsed"]]
  draws <- greta::mcmc(
    model = model,
    n_samples = n_samples,
    warmup = warmup,
    chains = chains,
    verbose = mcmc_verbose
  )
  run_time <- proc.time()[["elapsed"]] - run_time_start

  if (isTRUE(verbose)) {
    message(
      "fb_greta(): fit complete in ",
      round(run_time, 1),
      " sec; ",
      chains,
      " chain(s) x ",
      n_samples,
      " samples; ",
      length(fb$greta_meta$arrays),
      " target parameter(s)."
    )
  }

  # ------- 3. Posterior summary on canonical-named axes. ---------- #

  all_draws <- do.call(rbind, lapply(draws, as.matrix))
  if (is.null(colnames(all_draws))) {
    stop(
      "greta MCMC returned draws without column names; ",
      "fb_greta() requires named target parameters.",
      call. = FALSE
    )
  }

  # Rename draw columns to canonical where mapped (greta uses
  # canonical-array-name plus an element index for vector-valued
  # targets, e.g., `b[1,1]`). We only rename the scalar columns; the
  # element-suffixed columns retain the original greta name with the
  # canonical-prefix substitution. This keeps the canonical-name
  # registry compatible with both scalar and vector targets.
  greta_names <- names(fb$greta_meta$canonical_map)
  canon <- unname(fb$greta_meta$canonical_map)
  for (i in seq_along(greta_names)) {
    if (identical(greta_names[i], canon[i])) {
      next
    }
    pat <- paste0("^", greta_names[i], "(\\[|$)")
    repl <- paste0(canon[i], "\\1")
    colnames(all_draws) <- sub(pat, repl, colnames(all_draws))
  }

  post_mean <- colMeans(all_draws)
  post_vcov <- if (ncol(all_draws) > 0L) cov(all_draws) else matrix(0, 0, 0)
  post_sd <- sqrt(diag(post_vcov))
  post_q025 <- apply(all_draws, 2L, quantile, probs = 0.025, na.rm = TRUE)
  post_q975 <- apply(all_draws, 2L, quantile, probs = 0.975, na.rm = TRUE)

  # ------- 4. Convergence diagnostics. ---------------------------- #
  conv <- tryCatch(
    {
      list(
        gelman = if (length(draws) >= 2L) {
          coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)
        } else {
          NULL
        },
        n_eff = tryCatch(coda::effectiveSize(draws), error = function(e) NULL)
      )
    },
    error = function(e) list(gelman = NULL, n_eff = NULL)
  )

  # ------- 5. Minimal $glm shim. ---------------------------------- #
  # A minimal GLM-shaped shim built from the posterior
  # summary so that emmeans / marginaleffects see *something*. The
  # shim carries posterior means under canonical names as
  # `coefficients`; the family is the inferred likelihood family.
  fam_obj <- tryCatch(
    eval(call(fb$family))(),
    error = function(e) gaussian()
  )
  glm_shim <- list(
    coefficients = post_mean,
    residuals = numeric(0),
    fitted.values = numeric(0),
    rank = length(post_mean),
    family = fam_obj,
    linear.predictors = numeric(0),
    deviance = NA_real_,
    aic = NA_real_,
    null.deviance = NA_real_,
    iter = NA_integer_,
    df.residual = NA_integer_,
    df.null = NA_integer_,
    converged = TRUE,
    boundary = FALSE,
    call = the_call
  )
  attr(glm_shim, "posterior_vcov") <- post_vcov
  attr(glm_shim, "posterior_sd") <- post_sd
  attr(glm_shim, "posterior_q025") <- post_q025
  attr(glm_shim, "posterior_q975") <- post_q975
  class(glm_shim) <- c("flexybayes_glm", "glm", "lm")

  # ------- 6. Assemble the classed result. ------------------------ #
  param_names <- colnames(all_draws)
  result <- list(
    glm = glm_shim,
    greta = list(
      model = model,
      draws = draws,
      greta_arrays = fb$greta_meta$arrays,
      env = new.env(parent = emptyenv())
    ),
    extras = list(
      summary = data.frame(
        param = param_names,
        mean = post_mean,
        sd = post_sd,
        q025 = post_q025,
        q975 = post_q975,
        stringsAsFactors = FALSE,
        row.names = NULL
      ),
      convergence = conv,
      variance_comps = NULL,
      blups = NULL,
      predictions = NULL,
      code = NA_character_, # no flexyBayes-generated code
      param_names = param_names,
      parse_info = list(
        source = "greta",
        canonical_map = fb$greta_meta$canonical_map
      ),
      call_info = list(
        fixed = quote(`<greta-direct entry>`),
        random = NULL,
        rcov = NULL,
        chains = chains,
        n_samples = n_samples,
        warmup = warmup,
        call = the_call
      ),
      run_time = run_time,
      model_info = list(
        family = fb$family,
        link = "identity",
        n_params = length(post_mean),
        n_fixed = length(post_mean),
        n_random = 0L,
        n_obs = fb$greta_meta$n_data %||% NA_integer_,
        likelihood = fb$greta_meta$likelihood,
        canonical_map = fb$greta_meta$canonical_map
      ),
      fb_terms = fb,
      backend_decision = list(
        backend = "greta",
        path = "direct_entry",
        gate_checks = NULL,
        reason = paste0(
          "fb_greta() bypasses lgm_gate(); ",
          "user-built model accepted as-is."
        )
      )
    )
  )
  # S3 dispatch order: subclass first (most-specific), then parent.
  # The reverse order breaks R's S3 dispatch (the leftmost match wins):
  # `print.flexybayes` would always shadow
  # `print.flexybayes_direct_greta`. The implementation uses the
  # R-correct order.
  class(result) <- c("flexybayes_direct_greta", "flexybayes", "list")
  result
}


#' Fit a flexyBayes model via the greta engine
#'
#' Engine pin: fits the model with greta (full Hamiltonian Monte Carlo)
#' only. This is sugar for [flexybayes()]`(..., backend = "greta")` and
#' accepts the same grammars --- an ASReml `fixed` / `random` / `rcov`
#' specification, a brms-style bar-grouped formula, or a native
#' `greta_model` graph built with `greta::model()`. A formula is lowered
#' through the shared emit path; a native graph is fit directly by
#' [greta::mcmc()] and returned as a `flexybayes_direct_greta` object.
#'
#' To attach a canonical-name map to a native graph (so [triangulate()]
#' can align it with an INLA or Stan fit), build the intermediate
#' representation first and pass it on: `fb_greta(fb_from_greta(model,`
#' `canonical_names = c(...)))`.
#'
#' @param ... Arguments passed to [flexybayes()] (e.g. `fixed`, `random`,
#'   `rcov`, `data`, `family`, `prior`, `syntax`), or a native
#'   `greta_model` / greta-source IR as the model-spec slot. The
#'   `backend` argument is pinned to `"greta"`; a conflicting `backend`
#'   value raises a structured refusal. The pre-v0.5.0 `model = `
#'   native-graph argument is remapped to the model-spec slot for
#'   call-compatibility, so `fb_greta(model = m)` still fits a native
#'   graph.
#'
#' @return For a formula: an object of class `"flexybayes"` (a greta
#'   fit). For a native `greta_model`: a `flexybayes_direct_greta`
#'   object (subclass of `"flexybayes"`) carrying the MCMC draws, a
#'   GLM-compatible shim, and the canonical-name map; see
#'   [flexybayes()] for the shared structure.
#'
#' @family flexyBayes engine pins
#' @seealso [flexybayes()] and [fb()] for the universal entry that picks
#'   a backend; [fb_from_greta()] for building a greta-source IR with a
#'   canonical-name map.
#' @examples
#' \dontrun{
#' # live greta fit -- needs a working Python / TensorFlow stack
#' data(sleepstudy, package = "lme4")
#' fit <- fb_greta(Reaction ~ Days + (1 | Subject), data = sleepstudy,
#'                 n_samples = 200, warmup = 200, chains = 1,
#'                 verbose = FALSE, mcmc_verbose = FALSE)
#' coef(fit)
#' }
#' @export
fb_greta <- function(...) {
  cl <- match.call()
  # Call-compatibility: the removed `model = ` native-graph argument maps
  # to the universal entry's model-spec slot (`fixed`), so a native graph
  # passed either positionally or as `model = ` reaches the native-greta
  # fitting path.
  nm <- names(cl)
  if (!is.null(nm) && "model" %in% nm) {
    if ("fixed" %in% nm) {
      stop(
        "fb_greta(): pass the model once -- `model` (no longer ",
        "accepted) and `fixed` name the same model-spec slot.",
        call. = FALSE
      )
    }
    nm[nm == "model"] <- "fixed"
    names(cl) <- nm
  }
  .fb_engine_pin("greta", cl, parent.frame())
}


# ---------------------------------------------------------------- #
# backend_decision() -- trivial-trace helper                          #
# ---------------------------------------------------------------- #

#' Backend dispatch trace for a flexyBayes fit
#'
#' Returns the dispatch trace recorded at fit time: which backend
#' was selected, which gate checks ran, and why. On fb_greta() fits
#' the trace is trivial (the user bypassed the gate by entering on
#' the greta-direct path).
#'
#' @param fit A `flexybayes` object.
#' @return A list with the following components. The first four
#'   are present on every fit; the four routing-trace
#'   fields are present on v0.3.6+ fits and NULL
#'   on earlier fits for backward compatibility.
#'   \describe{
#'     \item{`backend`}{Character; one of `"greta"`, `"inla"`,
#'       `"brms"`, `"gretaR"`.}
#'     \item{`path`}{Character; the dispatch path token.}
#'     \item{`gate_checks`}{List or NULL; the `lgm_gate()` check
#'       trail (failures on refusal; capabilities on accept).}
#'     \item{`reason`}{Character; the dispatch-decision rationale.}
#'     \item{`preflight_summary`}{An `<fb_preflight>` object or
#'       NULL. Populated when `.fb_preflight()` ran (>1e5-row
#'       path); NULL on the small-data fast path.}
#'     \item{`representation_plan`}{Named list of slim per-term
#'       entries `(term_id, representation_class, justification)`
#'       derived from `preflight_summary`; NULL when no preflight.}
#'     \item{`rejected_routes`}{List of `(backend, reason)` pairs
#'       for the backends considered but not chosen. Empty for
#'       explicit user requests (the routing policy is bypassed
#'       when the backend is named directly).}
#'     \item{`routing_policy_version`}{Character; e.g.
#'       `"stage5a_v1"`. The audit-anchor for reproducibility -- a
#'       policy change bumps this string.}
#'   }
#' @export
backend_decision <- function(fit) {
  if (!inherits(fit, "flexybayes") && !inherits(fit, "flexybayes_inla")) {
    stop(
      "`fit` must be a `flexybayes` (or `flexybayes_inla`) object.",
      call. = FALSE
    )
  }
  decision <- fit$extras$backend_decision
  if (is.null(decision)) {
    # Older fit (no recorded decision): synthesise a minimal
    # legacy trace so downstream tooling never sees NULL.
    return(list(
      backend = if (inherits(fit, "flexybayes_inla")) "inla" else "greta",
      path = "legacy_no_recorded_decision",
      gate_checks = NULL,
      reason = paste0(
        "fit predates dispatch-trace recording; ",
        "backend inferred from class."
      ),
      preflight_summary = NULL,
      representation_plan = NULL,
      rejected_routes = list(),
      routing_policy_version = NA_character_
    ))
  }
  # v0.3.6+: surface the four new fields with
  # NULL defaults so earlier fits (those that recorded the
  # 4-field trace before v0.3.6 landed) read out with the
  # same shape downstream tooling sees on v0.3.6+ fits. Construct
  # a fresh canonical 8-field list with explicit NULL values --
  # `decision[["x"]] <- NULL` removes the slot rather than setting
  # it to NULL, which would silently drop the field on every read.
  rejected <- decision$rejected_routes
  if (is.null(rejected)) {
    rejected <- list()
  }
  version <- decision$routing_policy_version
  if (is.null(version)) {
    version <- NA_character_
  }
  list(
    backend = decision$backend,
    path = decision$path,
    gate_checks = decision$gate_checks,
    reason = decision$reason,
    preflight_summary = decision$preflight_summary,
    representation_plan = decision$representation_plan,
    rejected_routes = rejected,
    routing_policy_version = version
  )
}
