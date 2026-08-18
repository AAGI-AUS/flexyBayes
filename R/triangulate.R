# triangulate -- cross-engine posterior comparison on a gated overlap.
#
# Takes two model fits produced on different backends (INLA and the brms
# / Stan passthrough) and computes per-parameter posterior comparison
# metrics: the 1D Wasserstein-1 distance, the difference in posterior
# means, the ratio of posterior SDs, and the tail drift (the 0.025 and
# 0.975 quantile differences). Each fit's own posterior mean and SD are
# reported alongside, so the size of a discrepancy can be read against
# the scale of the parameter.
#
# Three gates run before any number is reported, because a distance
# between two posteriors means nothing until the two engines were asked
# the same question:
#
#   (1) Comparability. Each fit carries a model fingerprint built at fit
#       time (R/model_fingerprint.R). Two fingerprints that disagree are
#       refused by name -- there is no override argument, because the
#       function's contract is same model, same data.
#   (2) Diagnostics. A fit whose own sampler or approximation did not
#       converge cannot support a verdict about another engine, so a
#       failing fit makes the overall status inconclusive.
#   (3) Matched priors. A variance component is compared only when both
#       fits record the same prior for it. The rest are reported as
#       not_compared with the reason, rather than compared silently.
#
# Backends use different internal parameter naming conventions:
#   - brms: e.g., "b_Intercept", "b_x", "sd_g__Intercept", "sigma"
#   - INLA: e.g., "(Intercept):1", "x:1", "g:1", "Precision for g"
# canonical_names() reconciles the common cases automatically. Users
# supply a `name_map` (named character vector or list) to align fit_b's
# parameter names to fit_a's canonical names where it does not. Without
# alignment, only literal-match parameters are compared.
#
# Internal posterior extraction is via the fb_as_draws_simple() generic
# with methods for `flexybayes` (the brms path) and `flexybayes_inla`
# (the INLA path). Both return a named list whose values are numeric
# vectors of posterior draws.

# ---------------------------------------------------------------- #
# Public-facing entry                                              #
# ---------------------------------------------------------------- #

#' Cross-engine posterior triangulation
#'
#' Compute per-parameter posterior comparison metrics across two
#' Bayesian fits produced on different backends (INLA and the brms /
#' Stan passthrough).
#'
#' For each parameter present in both fits (post `name_map`), the
#' returned table reports each fit's posterior mean and SD, the
#' difference in means, the ratio of SDs, the Q2.5 and Q97.5 differences
#' (tail drift), and the 1D empirical Wasserstein-1 distance. Parameters
#' present in only one fit are reported in `only_a` / `only_b`.
#'
#' `triangulate()` is a diagnostic on a gated overlap, not the package's
#' validation contract. Agreement is not validation: the parser, the
#' intermediate representation, the prior interlingua and the data
#' preparation are shared, so two engines can agree on a mistranslated
#' model. What agreement does support is the narrower claim that neither
#' engine's approximation is distorting the posterior of the compared
#' parameters, and it supports even that only inside the three gates
#' described below.
#'
#' @section The three gates:
#' **Comparability.** Fits made by [flexybayes()] carry a model
#' fingerprint: the canonical formula triple, the family and link, the
#' data dimensions, column names and a content digest, and the prior
#' recorded for each variance component. Two fingerprints that disagree
#' raise a typed refusal (`flexybayes_refusal_triangulate_incomparable_fits`)
#' naming the first element that differs. There is no override argument:
#' comparing two models is a different operation, and it is not this one.
#' A fit built outside the package carries no fingerprint, so
#' comparability cannot be verified and the overall status is
#' `inconclusive`.
#'
#' **Diagnostics.** A brms fit passes when every model parameter has
#' R-hat at or below 1.01 and bulk effective sample size at or above 400,
#' with no divergent transitions; an INLA fit passes its own
#' numerical-confirm gate (mode status and a finite marginal likelihood).
#' A failing fit yields status `inconclusive` and no parameter verdicts:
#' disagreement between an unconverged fit and a converged one is not a
#' finding.
#'
#' **Matched priors.** A variance component (`sigma`, `sd_<group>`,
#' `cor_<group>`, `rho_*`) is compared only when both fits record the
#' same prior for it -- the shared default this package injects before
#' dispatch, or an identical explicit [fb_prior()]. Components outside
#' that record are reported as `not_compared` with the reason, most often
#' that the term is outside the default-prior walker and each engine
#' chose its own hyperprior. Fixed-effect coefficients are compared
#' without a prior match: both engines put a vague normal on them, and
#' the fixed structure itself is already part of the fingerprint.
#'
#' @section Thresholds and what they are not:
#' A parameter is called concordant when the absolute mean difference is
#' at most 0.1 posterior SD, the Wasserstein-1 distance is at most 0.1
#' posterior SD, and the SD ratio lies in \[0.9, 1.1\]; the reference
#' scale is the larger of the two posterior SDs. The overall status is
#' `discordant` when any compared parameter is discordant, `concordant`
#' when all are, and `inconclusive` when a gate blocked the comparison or
#' nothing was left to compare.
#'
#' These cut-offs are **heuristics**, in the sense Gelman et al. (2020)
#' use for cross-implementation checks in the Bayesian workflow: they are
#' chosen to flag differences a reader should look at, not derived from a
#' sampling distribution, and they carry no error rate. A `concordant`
#' verdict is not a calibration statement. flexyBayes ships no
#' simulation-based calibration of its own fits (Talts et al., 2018), so
#' whether either engine is calibrated for a given model and dataset
#' remains an open question and future work.
#'
#' Read `discordant` as "look at this", not as "one engine is wrong".
#' Wasserstein-1 compares the whole marginal, so two posteriors whose
#' means and SDs agree can still be flagged on shape alone -- which is
#' the expected signature of a Laplace marginal for a variance component
#' next to an HMC one. All three metrics are estimated from finite draws,
#' and the Wasserstein-1 estimate in particular is biased upward at small
#' `n_samples`; raise `n_samples` and re-run before reading a borderline
#' flag as a finding.
#'
#' Backends use different parameter naming conventions; supply
#' `name_map = c(<fit_b name> = <canonical name>, ...)` to align
#' them. Without alignment, only literal-match parameter names are
#' compared.
#'
#' Backends also use different *parameter scales* for variance
#' components -- INLA reports precision (`Precision for g`), brms
#' reports standard deviation (`sd_g__Intercept`). Supply
#' `transform_a` / `transform_b` -- named lists of one-argument
#' functions keyed by parameter name -- to put the two posteriors on
#' a common scale before comparison. Transforms are applied first;
#' then `name_map` aligns the (already-transformed) fit_b names to
#' fit_a's canonical names. Names in `transform_b` therefore refer to
#' fit_b's *original* parameter names, not the post-`name_map`
#' canonical names.
#'
#' @section Aggregated fits and matched priors:
#' When one of the inputs is an aggregated-gaussian fit (cell-level
#' sufficient statistics rather than the per-row likelihood), the
#' posteriors being compared are only directly comparable when the two
#' fits share priors. The aggregated path combines the cell-mean
#' likelihood with a precision prior carrying a closed-form correction
#' that absorbs the within-cell sum-of-squares. Under the *legacy scalar*
#' prior this recovers the per-row posterior to numerical precision, so
#' the aggregated fit is tagged `prior_parametrization =
#' "per_row_equivalent"` (visible in the aggregated `print()` / `summary()`
#' and in [canonical_names()]). Under the package's own auto-default it
#' is tagged `"package_default"`, which names the prior without claiming
#' that equivalence. When an explicit prior is supplied the
#' fit is tagged `"custom"`: the equivalence against a *default-prior*
#' per-row fit no longer holds, and on the aggregated INLA path the
#' observation-precision prior is not plumbed through, so a custom
#' residual prior is silently not applied there. Before reading the
#' agreement metrics on a custom-prior aggregated fit, confirm both
#' inputs carry the same prior with [prior_summary()].
#'
#' @param fit_a A fit object carrying a `fb_as_draws_simple` method,
#'   which means `flexybayes_brms` from the brms path or
#'   `flexybayes_inla` from the INLA path.
#' @param fit_b A second fit object of the same kind, typically produced
#'   by the other backend on the same model and data.
#' @param name_map An optional named character vector or list mapping
#'   `fit_b`'s parameter names (left) to canonical names matching
#'   `fit_a` (right).
#' @param transform_a,transform_b An optional named list of functions,
#'   each taking a numeric vector of posterior draws and returning a
#'   numeric vector of the same length. Names key parameters in
#'   `fit_a` / `fit_b` (using each fit's *original* parameter
#'   names). Common use: pass
#'   `transform_b = list("Precision for g" = function(x) 1 / sqrt(x))`
#'   to convert INLA's precision draws to standard-deviation scale
#'   so they line up with brms's `sd_g__Intercept`. Optional, and
#'   normally unnecessary -- [canonical_names()] supplies the common
#'   transforms.
#' @param n_samples A single integer giving the number of posterior
#'   samples to draw for fits whose extractor needs sampling, such as
#'   INLA via `INLA::inla.posterior.sample`.
#' @param data_independence A single logical declaring whether the two fits were
#'   built on independently-sourced data. `triangulate()` measures inter-fit
#'   *agreement*, and the backend-independence registry certifies code (not
#'   data) independence -- so if both fits share the same upstream data, a
#'   fabricated data fact is common-mode and their agreement cannot detect it.
#'   `TRUE` declares the data independent (no caveat); `FALSE` (same data) or
#'   `NA` (the default, undeclared) attach a `shared_upstream_caveat` field to
#'   the result, surfaced prominently by the print method, so agreement is never
#'   silently mistaken for corroboration of a shared data fact (Independent
#'   Oracle Principle).
#' @returns A `triangulate_result` S3 object (a list). Key fields:
#'   `metrics` (a data.frame with one row per common parameter and the
#'   columns `param`, `mean_a`, `mean_b`, `mean_diff`, `sd_a`, `sd_b`,
#'   `sd_ratio`, `q025_diff`, `q975_diff`, `wasserstein_1`,
#'   `mean_shift_sd`, `w1_sd`, `verdict`, `reason`), `status` (one of
#'   `"concordant"`, `"discordant"`, `"inconclusive"`),
#'   `status_reasons` (character, why a status is inconclusive or
#'   discordant), `comparability` (the fingerprint verdict),
#'   `diagnostics` (the per-fit diagnostic summaries), `common`
#'   (character), `only_a`, `only_b`, `n_common`, `n_compared`,
#'   `source_a`, `source_b`.
#' @references
#' Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter,
#' B., Yao, Y., et al. (2020). Bayesian workflow. *arXiv preprint*
#' arXiv:2011.01808.
#'
#' Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A.
#' (2018). Validating Bayesian inference algorithms with
#' simulation-based calibration. *arXiv preprint* arXiv:1804.06788.
#'
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., et al. (2021).
#' Rank-normalization, folding, and localization: an improved R-hat for
#' assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667--718.
#' @examples
#' # Live INLA posterior sampling can fail in restricted-process
#' # check environments (the `inla.posterior.sample` parallelism
#' # check trips); the example uses `\dontrun{}` deliberately. On
#' # an interactive install with INLA + brms + sn available it
#' # runs in a few seconds.
#' \dontrun{
#' if (requireNamespace("brms", quietly = TRUE) &&
#'     requireNamespace("INLA", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(y = rnorm(40), x = rnorm(40),
#'                   g = factor(rep(1:5, 8)))
#'   fit_b <- fb(y ~ x + (1 | g), data = d, backend = "brms",
#'               n_samples = 500, warmup = 500, chains = 2,
#'               verbose = FALSE, mcmc_verbose = FALSE)
#'   fit_i <- fb(y ~ x + (1 | g), data = d, backend = "inla",
#'               verbose = FALSE)
#'   # canonical_names() aligns the names and the precision-to-SD
#'   # transform, so no name_map is needed in the common cases.
#'   print(triangulate(fit_b, fit_i))
#' }
#' }
#' @export
triangulate <- function(
  fit_a,
  fit_b,
  name_map = NULL,
  transform_a = NULL,
  transform_b = NULL,
  n_samples = 1000L,
  data_independence = NA
) {
  if (!is.logical(data_independence) || length(data_independence) != 1L) {
    stop("`data_independence` must be a single logical (TRUE / FALSE / NA).",
         call. = FALSE)
  }

  # Gate 1 -- comparability. Ahead of the draws because extracting an
  # INLA posterior means sampling it, and there is no point sampling two
  # posteriors that are not answers to the same question.
  comparability <- .triangulate_comparability(fit_a, fit_b)

  # Gate 2 -- each fit's own diagnostics.
  diag_a <- .fb_fit_diagnostics(fit_a)
  diag_b <- .fb_fit_diagnostics(fit_b)

  draws_a <- fb_as_draws_simple(fit_a, n_samples = n_samples)
  draws_b <- fb_as_draws_simple(fit_b, n_samples = n_samples)

  if (!is.list(draws_a) || is.null(names(draws_a))) {
    stop(
      "`fb_as_draws_simple(fit_a)` must return a named list of ",
      "numeric vectors.",
      call. = FALSE
    )
  }
  if (!is.list(draws_b) || is.null(names(draws_b))) {
    stop(
      "`fb_as_draws_simple(fit_b)` must return a named list of ",
      "numeric vectors.",
      call. = FALSE
    )
  }

  # When transforms / name_map are not user-supplied, fall
  # back to the per-backend canonical-name registry. User-supplied
  # values always win over the registry per the explicit-over-implicit
  # principle. Resolve canonical_names() defensively (returns empty
  # list when no mapper is registered for the backend / class).
  cn_a <- tryCatch(canonical_names(fit_a, drop = TRUE), error = function(e) {
    list(map = character(0), transform = list())
  })
  cn_b <- tryCatch(canonical_names(fit_b, drop = TRUE), error = function(e) {
    list(map = character(0), transform = list())
  })

  if (is.null(transform_a) && length(cn_a$transform) > 0L) {
    transform_a <- cn_a$transform
  }
  if (is.null(transform_b) && length(cn_b$transform) > 0L) {
    transform_b <- cn_b$transform
  }

  # Apply per-parameter transforms before any renaming so callers can
  # key transforms by each fit's native parameter names.
  draws_a <- .triangulate_apply_transform(draws_a, transform_a, "transform_a")
  draws_b <- .triangulate_apply_transform(draws_b, transform_b, "transform_b")

  if (!is.null(name_map)) {
    if (
      !(is.list(name_map) || is.character(name_map)) ||
        is.null(names(name_map))
    ) {
      stop(
        "`name_map` must be a named character vector or list ",
        "mapping fit_b's parameter names to canonical names.",
        call. = FALSE
      )
    }
  }

  # Combined-precedence rename: user-supplied name_map > registry
  # > identity. The user map applies only to fit_b (existing
  # contract); the registry applies to both fits. For each native
  # name we check name_map first (if applicable) then fall back to
  # the canonical-name registry; un-mapped names retain their
  # backend-native form.
  .resolve_canonical <- function(nm, reg_map, user_map) {
    if (!is.null(user_map) && nm %in% names(user_map)) {
      return(as.character(user_map[[nm]]))
    }
    if (!is.null(reg_map) && nm %in% names(reg_map)) {
      return(as.character(reg_map[[nm]]))
    }
    nm
  }

  names(draws_a) <- vapply(
    names(draws_a),
    .resolve_canonical,
    character(1),
    reg_map = cn_a$map,
    user_map = NULL
  )
  names(draws_b) <- vapply(
    names(draws_b),
    .resolve_canonical,
    character(1),
    reg_map = cn_b$map,
    user_map = name_map
  )

  common <- intersect(names(draws_a), names(draws_b))
  only_a <- setdiff(names(draws_a), common)
  only_b <- setdiff(names(draws_b), common)

  metrics <- if (length(common) > 0L) {
    rows <- lapply(common, function(p) {
      .triangulate_one(p, draws_a[[p]], draws_b[[p]])
    })
    do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  } else {
    .triangulate_empty_table()
  }

  # Gate 3 -- matched priors, then the per-parameter verdicts and the
  # overall status the two gates above may already have decided.
  metrics <- .triangulate_verdicts(
    metrics,
    comparability = comparability,
    diag_a = diag_a,
    diag_b = diag_b
  )
  overall <- .triangulate_status(
    metrics,
    comparability = comparability,
    diag_a = diag_a,
    diag_b = diag_b
  )

  # Label the pair with its independence axis vocabulary.
  # The lookup is symmetric (sorted pair key); an unregistered pair
  # (same backend, or a backend with no registered claims) yields an
  # empty axis set and the report renders without the label rather than
  # refusing.
  src_a <- .triangulate_source(fit_a)
  src_b <- .triangulate_source(fit_b)
  indep <- .lookup_pair_independence(c(src_a, src_b))

  # FX-10 (Independent Oracle Principle). triangulate measures inter-fit
  # *agreement* only; the backend-independence registry certifies CODE
  # independence (algorithmic / implementation / specification), never DATA-fact
  # independence. If both fits consumed the same upstream `data`, a fabricated
  # data fact is common-mode: every backend ingests it, every backend agrees,
  # and tight agreement is reported as strong consensus. Unless the caller
  # declares the fits used independently-sourced data, we attach a caveat (and
  # warn) so agreement is not mistaken for corroboration of a shared data fact.
  shared_upstream_caveat <- if (isTRUE(data_independence)) {
    NA_character_
  } else {
    paste0(
      "triangulate measures inter-fit agreement, not correspondence: the ",
      "backend-independence registry certifies code independence, not ",
      "data-fact independence. ",
      if (identical(data_independence, FALSE)) {
        "Both fits consumed the SAME data, so agreement is common-mode and "
      } else {
        "Data independence was not declared, so if both fits share the same "
      },
      "a fabricated upstream data fact would not be detected by their ",
      "agreement. Declare data_independence = TRUE only when the fits used ",
      "independently-sourced data."
    )
  }

  structure(
    list(
      metrics = metrics,
      status = overall$status,
      status_reasons = overall$reasons,
      comparability = comparability,
      diagnostics = list(a = diag_a, b = diag_b),
      common = common,
      only_a = only_a,
      only_b = only_b,
      n_common = length(common),
      n_compared = sum(metrics$verdict %in% c("concordant", "discordant")),
      source_a = src_a,
      source_b = src_b,
      independence = if (is.null(indep)) character(0) else indep$axes,
      axis_justification = if (is.null(indep)) {
        NA_character_
      } else {
        indep$justification
      },
      data_independence = data_independence,
      shared_upstream_caveat = shared_upstream_caveat
    ),
    class = c("triangulate_result", "list")
  )
}

is_triangulate_result <- function(x) inherits(x, "triangulate_result")

# Apply a named list of one-argument transforms to a draws list.
# Names not present in `draws` are silently ignored (so a single
# transform spec can be reused across models that share only some
# parameters). Each function must return a vector of the same length
# as the input draws.
.triangulate_apply_transform <- function(draws, transform, label) {
  if (is.null(transform)) {
    return(draws)
  }
  if (
    !is.list(transform) ||
      is.null(names(transform)) ||
      any(!nzchar(names(transform)))
  ) {
    stop(
      "`",
      label,
      "` must be a named list of one-argument ",
      "functions keyed by parameter name.",
      call. = FALSE
    )
  }
  for (nm in names(transform)) {
    fn <- transform[[nm]]
    if (!is.function(fn)) {
      stop("`", label, "[[\"", nm, "\"]]` must be a function.", call. = FALSE)
    }
    if (nm %in% names(draws)) {
      x_in <- draws[[nm]]
      x_out <- fn(x_in)
      if (!is.numeric(x_out) || length(x_out) != length(x_in)) {
        stop(
          "`",
          label,
          "[[\"",
          nm,
          "\"]]` must return a numeric ",
          "vector of length ",
          length(x_in),
          ".",
          call. = FALSE
        )
      }
      draws[[nm]] <- as.numeric(x_out)
    }
  }
  draws
}

# ---------------------------------------------------------------- #
# Per-parameter metrics                                            #
# ---------------------------------------------------------------- #

.triangulate_one <- function(name, a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  m_a <- mean(a, na.rm = TRUE)
  m_b <- mean(b, na.rm = TRUE)
  s_a <- stats::sd(a, na.rm = TRUE)
  s_b <- stats::sd(b, na.rm = TRUE)
  q025_a <- stats::quantile(a, 0.025, names = FALSE, na.rm = TRUE)
  q025_b <- stats::quantile(b, 0.025, names = FALSE, na.rm = TRUE)
  q975_a <- stats::quantile(a, 0.975, names = FALSE, na.rm = TRUE)
  q975_b <- stats::quantile(b, 0.975, names = FALSE, na.rm = TRUE)

  list(
    param = name,
    mean_a = m_a,
    mean_b = m_b,
    mean_diff = m_a - m_b,
    sd_a = s_a,
    sd_b = s_b,
    sd_ratio = if (is.finite(s_b) && s_b > 0) s_a / s_b else NA_real_,
    q025_diff = q025_a - q025_b,
    q975_diff = q975_a - q975_b,
    wasserstein_1 = .wasserstein1_1d(a, b)
  )
}

.triangulate_empty_table <- function() {
  data.frame(
    param = character(0),
    mean_a = numeric(0),
    mean_b = numeric(0),
    mean_diff = numeric(0),
    sd_a = numeric(0),
    sd_b = numeric(0),
    sd_ratio = numeric(0),
    q025_diff = numeric(0),
    q975_diff = numeric(0),
    wasserstein_1 = numeric(0),
    mean_shift_sd = numeric(0),
    w1_sd = numeric(0),
    verdict = character(0),
    reason = character(0),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------- #
# Gate 1 -- comparability                                          #
# ---------------------------------------------------------------- #

# Compare the two fits' model fingerprints. Refuses on a difference;
# otherwise returns the verdict, which is either verified (both fits
# carried a fingerprint and they agree) or unverified (at least one fit
# carries none, so nothing was checked and the caller must be told so
# rather than reassured by silence).
.triangulate_comparability <- function(fit_a, fit_b) {
  fp_a <- .fb_fit_fingerprint(fit_a)
  fp_b <- .fb_fit_fingerprint(fit_b)

  if (is.null(fp_a) || is.null(fp_b)) {
    missing_on <- c(
      if (is.null(fp_a)) "fit_a" else NULL,
      if (is.null(fp_b)) "fit_b" else NULL
    )
    return(list(
      verified = FALSE,
      fingerprint_a = fp_a,
      fingerprint_b = fp_b,
      note = paste0(
        "no model fingerprint on ",
        paste(missing_on, collapse = " and "),
        ", so it could not be checked that the two fits share a model, a ",
        "dataset and a prior. Fits made by flexybayes() carry one; fits ",
        "imported or built by hand do not."
      )
    ))
  }

  diff <- .fb_fingerprint_first_difference(fp_a, fp_b)
  if (!is.null(diff)) {
    stop(.fb_refusal_condition(
      reason_code = "triangulate_incomparable_fits",
      message = paste0(
        "triangulate() compares two fits of the same model on the same ",
        "data. These two differ on ", diff$label, ":\n",
        "  fit_a: ", diff$a, "\n",
        "  fit_b: ", diff$b, "\n",
        "A distance between posteriors of two different questions is not a ",
        "disagreement between engines. Re-fit both on the same model, data ",
        "and prior, or compare them by reading the two summaries -- there ",
        "is no argument that lifts this check."
      ),
      differing_element = diff$element,
      value_a = diff$a,
      value_b = diff$b
    ))
  }

  list(
    verified = TRUE,
    fingerprint_a = fp_a,
    fingerprint_b = fp_b,
    note = NA_character_
  )
}


# ---------------------------------------------------------------- #
# Gate 2 -- per-fit diagnostics                                    #
# ---------------------------------------------------------------- #

# Thresholds. R-hat 1.01 and bulk ESS 400 are the rank-normalised
# workflow defaults (Vehtari et al. 2021); they are heuristics with no
# error rate, and the roxygen says so.
.FB_TRIANGULATE_RHAT_MAX <- 1.01
.FB_TRIANGULATE_ESS_MIN <- 400

# Parameters excluded from the sampler diagnostics: Stan's own bookkeeping
# quantities are not model parameters, and their mixing is not what the
# gate is about.
.FB_TRIANGULATE_DIAG_EXCLUDE <- c("lp__", "lprior", "lp")

# Drop the bookkeeping quantities from a named diagnostic vector. An
# unnamed vector is returned whole: subsetting on `names(v) %in% ...`
# when `names(v)` is NULL empties the vector, which would silently turn
# a failing fit into one with no diagnostics at all.
.fb_drop_bookkeeping <- function(v) {
  if (is.null(v) || is.null(names(v))) {
    return(v)
  }
  v[!(names(v) %in% .FB_TRIANGULATE_DIAG_EXCLUDE)]
}

# Name of the element at `idx`, or a placeholder when the vector has no
# names (a fit whose diagnostics arrived without variable labels).
.fb_diag_label <- function(v, idx) {
  nms <- names(v)
  if (is.null(nms) || is.na(nms[[idx]])) "a model parameter" else nms[[idx]]
}

# Normalise whatever diagnostics an engine records into one shape:
# ok (TRUE / FALSE / NA when the object carries none), the reasons a
# FALSE was reached, and the numbers behind them.
.fb_fit_diagnostics <- function(fit) {
  if (inherits(fit, "flexybayes_inla")) {
    return(.fb_inla_diagnostics(fit))
  }
  if (inherits(fit, "flexybayes_brms")) {
    return(.fb_brms_diagnostics(fit))
  }
  list(
    engine = .triangulate_source(fit),
    ok = NA,
    reasons = "this fit records no convergence diagnostics",
    max_rhat = NA_real_,
    min_ess_bulk = NA_real_,
    n_divergent = NA_integer_
  )
}

.fb_inla_diagnostics <- function(fit) {
  nc <- fit$num_check
  # The aggregated INLA path builds its own fit object and does not run
  # the post-fit numerical confirm, so read it off the raw INLA object
  # rather than reporting the diagnostics as unavailable.
  if (is.null(nc) && !is.null(fit$inla)) {
    nc <- tryCatch(.lgm_check_numerical(fit$inla), error = function(e) NULL)
  }
  if (is.null(nc)) {
    return(list(
      engine = "inla",
      ok = NA,
      reasons = "this fit records no numerical-confirm result",
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_,
      n_divergent = NA_integer_
    ))
  }
  list(
    engine = "inla",
    ok = isTRUE(nc$pass),
    reasons = if (isTRUE(nc$pass)) character(0) else nc$reasons,
    max_rhat = NA_real_,
    min_ess_bulk = NA_real_,
    n_divergent = NA_integer_
  )
}

.fb_brms_diagnostics <- function(fit) {
  conv <- fit$extras$convergence
  out <- list(
    engine = "brms",
    ok = NA,
    reasons = character(0),
    max_rhat = NA_real_,
    min_ess_bulk = NA_real_,
    n_divergent = NA_integer_
  )
  if (is.null(conv)) {
    out$reasons <- "this fit records no sampler diagnostics"
    return(out)
  }

  rhat <- tryCatch(
    conv$gelman$psrf[, "Point est.", drop = FALSE],
    error = function(e) numeric(0)
  )
  if (is.matrix(rhat)) {
    rhat <- stats::setNames(as.numeric(rhat[, 1L]), rownames(rhat))
  }
  rhat <- .fb_drop_bookkeeping(rhat)
  rhat <- rhat[is.finite(rhat)]
  ess <- .fb_drop_bookkeeping(conv$n_eff)
  ess <- ess[is.finite(ess)]
  ndiv <- conv$n_divergent

  reasons <- character(0)
  ok <- TRUE
  if (length(rhat) == 0L && length(ess) == 0L) {
    out$reasons <- "this fit records no sampler diagnostics"
    return(out)
  }
  if (length(rhat) > 0L) {
    out$max_rhat <- max(rhat)
    if (out$max_rhat > .FB_TRIANGULATE_RHAT_MAX) {
      ok <- FALSE
      reasons <- c(reasons, sprintf(
        "max R-hat %.4f on `%s` exceeds %.2f",
        out$max_rhat,
        .fb_diag_label(rhat, which.max(rhat)),
        .FB_TRIANGULATE_RHAT_MAX
      ))
    }
  }
  if (length(ess) > 0L) {
    out$min_ess_bulk <- min(ess)
    if (out$min_ess_bulk < .FB_TRIANGULATE_ESS_MIN) {
      ok <- FALSE
      reasons <- c(reasons, sprintf(
        "min bulk ESS %.0f on `%s` is below %d",
        out$min_ess_bulk,
        .fb_diag_label(ess, which.min(ess)),
        as.integer(.FB_TRIANGULATE_ESS_MIN)
      ))
    }
  }
  if (!is.null(ndiv) && is.finite(ndiv)) {
    out$n_divergent <- as.integer(ndiv)
    if (out$n_divergent > 0L) {
      ok <- FALSE
      reasons <- c(reasons, sprintf(
        "%d divergent transition%s",
        out$n_divergent,
        if (out$n_divergent == 1L) "" else "s"
      ))
    }
  }
  out$ok <- ok
  out$reasons <- reasons
  out
}


# ---------------------------------------------------------------- #
# Gate 3 -- matched priors, and the verdicts                       #
# ---------------------------------------------------------------- #

# Verdict thresholds, on the scale of the wider of the two posteriors.
.FB_TRIANGULATE_MEAN_SHIFT_MAX <- 0.1
.FB_TRIANGULATE_W1_MAX <- 0.1
.FB_TRIANGULATE_SD_RATIO_BAND <- c(0.9, 1.1)

# Canonical names of the parameters the matched-prior gate governs. The
# fixed-effect coefficients are deliberately not in it: both engines put
# a vague normal on them, the difference between two vague normals does
# not move an identified coefficient, and the fixed structure is already
# compared through the fingerprint.
.fb_is_variance_component_param <- function(param) {
  grepl("^(sigma$|sd_|cor_|rho_|Precision )", param)
}

# Attach the per-parameter prior verdict and the concordance verdict.
.triangulate_verdicts <- function(metrics, comparability, diag_a, diag_b) {
  if (nrow(metrics) == 0L) {
    return(metrics)
  }
  fp_a <- comparability$fingerprint_a
  fp_b <- comparability$fingerprint_b

  metrics$mean_shift_sd <- NA_real_
  metrics$w1_sd <- NA_real_
  verdict <- character(nrow(metrics))
  reason <- rep(NA_character_, nrow(metrics))

  # A fit that failed its own diagnostics supports no verdict about
  # anything. The distances are still reported -- they are what a reader
  # needs to see to judge how bad the failure is -- but no parameter is
  # called concordant or discordant on the strength of a posterior the
  # sampler did not explore.
  failed <- c(
    if (isFALSE(diag_a$ok)) "fit_a" else NULL,
    if (isFALSE(diag_b$ok)) "fit_b" else NULL
  )
  if (length(failed) > 0L) {
    metrics$verdict <- "not_compared"
    metrics$reason <- paste0(
      paste(failed, collapse = " and "),
      " failed the convergence gate, so no parameter verdict is supportable"
    )
    return(metrics)
  }

  for (i in seq_len(nrow(metrics))) {
    p <- metrics$param[[i]]
    excl <- .triangulate_prior_exclusion(p, fp_a, fp_b)
    if (!is.na(excl)) {
      verdict[[i]] <- "not_compared"
      reason[[i]] <- excl
      next
    }
    s_ref <- max(metrics$sd_a[[i]], metrics$sd_b[[i]], na.rm = TRUE)
    if (!is.finite(s_ref) || s_ref <= 0) {
      verdict[[i]] <- "not_compared"
      reason[[i]] <- paste0(
        "no posterior spread to scale the comparison against (both SDs are ",
        "zero or not finite)"
      )
      next
    }
    metrics$mean_shift_sd[[i]] <- abs(metrics$mean_diff[[i]]) / s_ref
    metrics$w1_sd[[i]] <- metrics$wasserstein_1[[i]] / s_ref
    ratio <- metrics$sd_ratio[[i]]
    breach <- c(
      if (metrics$mean_shift_sd[[i]] > .FB_TRIANGULATE_MEAN_SHIFT_MAX) {
        sprintf("mean shift %.2f SD", metrics$mean_shift_sd[[i]])
      },
      if (metrics$w1_sd[[i]] > .FB_TRIANGULATE_W1_MAX) {
        sprintf("Wasserstein-1 %.2f SD", metrics$w1_sd[[i]])
      },
      if (
        !is.finite(ratio) ||
          ratio < .FB_TRIANGULATE_SD_RATIO_BAND[[1L]] ||
          ratio > .FB_TRIANGULATE_SD_RATIO_BAND[[2L]]
      ) {
        sprintf("SD ratio %.2f", ratio)
      }
    )
    if (length(breach) > 0L) {
      verdict[[i]] <- "discordant"
      reason[[i]] <- paste(breach, collapse = "; ")
    } else {
      verdict[[i]] <- "concordant"
    }
  }

  metrics$verdict <- verdict
  metrics$reason <- reason
  metrics
}

# Why a parameter is outside the matched-prior overlap, or NA when it is
# inside it. Non-variance-component parameters are always inside.
.triangulate_prior_exclusion <- function(param, fp_a, fp_b) {
  if (!.fb_is_variance_component_param(param)) {
    return(NA_character_)
  }
  if (is.null(fp_a) || is.null(fp_b)) {
    return(NA_character_)
  }

  eng <- c(fp_a$engine_default_params, fp_b$engine_default_params)
  if (param %in% names(eng)) {
    return(paste0("no shared prior: ", unname(eng[[param]])))
  }

  pa <- fp_a$priors %||% character(0)
  pb <- fp_b$priors %||% character(0)
  in_a <- param %in% names(pa)
  in_b <- param %in% names(pb)
  if (in_a && in_b) {
    # An unequal recorded prior is caught by the comparability gate and
    # never reaches here, so equality is the only remaining case.
    return(NA_character_)
  }
  paste0(
    "no shared prior: recorded on ",
    if (in_a) "fit_a only" else if (in_b) "fit_b only" else "neither fit"
  )
}


# ---------------------------------------------------------------- #
# Overall status                                                   #
# ---------------------------------------------------------------- #

.triangulate_status <- function(metrics, comparability, diag_a, diag_b) {
  reasons <- character(0)

  if (!isTRUE(comparability$verified)) {
    reasons <- c(reasons, paste0("comparability unverified: ",
                                 comparability$note))
  }
  for (nm in c("a", "b")) {
    d <- if (identical(nm, "a")) diag_a else diag_b
    if (isFALSE(d$ok)) {
      reasons <- c(reasons, sprintf(
        "fit_%s (%s) failed its own diagnostics: %s",
        nm,
        d$engine,
        paste(d$reasons, collapse = "; ")
      ))
    } else if (is.na(d$ok)) {
      reasons <- c(reasons, sprintf(
        "fit_%s (%s) diagnostics unavailable: %s",
        nm,
        d$engine,
        paste(d$reasons, collapse = "; ")
      ))
    }
  }

  compared <- metrics$verdict %in% c("concordant", "discordant")
  if (sum(compared) == 0L) {
    reasons <- c(
      reasons,
      paste0(
        "no parameter survived the comparability, diagnostic and ",
        "matched-prior gates"
      )
    )
  }

  status <- if (length(reasons) > 0L) {
    "inconclusive"
  } else if (any(metrics$verdict == "discordant")) {
    "discordant"
  } else {
    "concordant"
  }

  if (identical(status, "discordant")) {
    disc <- metrics[metrics$verdict == "discordant", , drop = FALSE]
    reasons <- sprintf("%s: %s", disc$param, disc$reason)
  }

  list(status = status, reasons = reasons)
}

# Empirical 1D Wasserstein-1 distance via quantile interpolation:
# W1 = integral over u in [0, 1] of |F_a^-1(u) - F_b^-1(u)| du.
# Approximated by a 99-point quantile grid (u in 0.01..0.99).
.wasserstein1_1d <- function(a, b) {
  if (length(a) == 0L || length(b) == 0L) {
    return(NA_real_)
  }
  qs <- seq(0.01, 0.99, by = 0.01)
  qa <- stats::quantile(a, qs, names = FALSE, na.rm = TRUE)
  qb <- stats::quantile(b, qs, names = FALSE, na.rm = TRUE)
  mean(abs(qa - qb))
}

# (An "R-hat-on-means" between-engine statistic was removed in the development
# version: pooling two different engines' posteriors as "chains" conflates
# genuine between-engine approximation bias with within-sampler non-convergence,
# so it is not a valid convergence diagnostic. Cross-engine discrepancy is
# reported by the distributional metrics -- wasserstein_1, sd_ratio, mean_diff,
# and the quantile differences.)

.triangulate_source <- function(fit) {
  if (inherits(fit, "flexybayes_brms")) {
    return("brms")
  }
  if (inherits(fit, "flexybayes_inla")) {
    return("inla")
  }
  if (inherits(fit, "flexybayes")) {
    return("greta")
  }
  paste(class(fit), collapse = "/")
}

# ---------------------------------------------------------------- #
# Per-fit posterior extraction (S3 generic)                        #
# ---------------------------------------------------------------- #

#' Extract per-parameter posterior draws from a model fit
#'
#' S3 generic used by `triangulate()` to extract a named list of
#' numeric posterior-draw vectors from each fit. The methods that reach
#' an active engine are `flexybayes_brms` and `flexybayes_inla`; the
#' `flexybayes` and `flexybayes_gretaR` methods read the greta-shaped
#' draws slot and are retained for objects lifted in by
#' [fb_from_greta()], since greta is quarantined as a fitting engine.
#' User-defined methods can extend the generic.
#'
#' @param fit A model fit object carrying a posterior, dispatched on by
#'   the methods listed above.
#' @param ... Method-specific arguments, such as `n_samples` for the
#'   INLA method.
#' @returns A named list of numeric vectors, one element per parameter,
#'   each element holding that parameter's posterior draws.
#' @export
fb_as_draws_simple <- function(fit, ...) UseMethod("fb_as_draws_simple")

#' @rdname fb_as_draws_simple
#' @keywords internal
#' @export
fb_as_draws_simple.flexybayes <- function(fit, ...) {
  if (is.null(fit$greta) || is.null(fit$greta$draws)) {
    stop(
      "Cannot extract draws: fit$greta$draws is missing. This method ",
      "reads the greta slot, which a fit on an active engine does not ",
      "carry -- a brms fit dispatches on `flexybayes_brms` and an INLA ",
      "fit on `flexybayes_inla`. Check the fit completed ",
      "(return_code = FALSE) and carries the slot its engine writes.",
      call. = FALSE
    )
  }
  m <- as.matrix(fit$greta$draws)
  cols <- colnames(m)
  if (is.null(cols)) {
    cols <- paste0("V", seq_len(ncol(m)))
  }
  setNames(lapply(seq_len(ncol(m)), function(j) as.numeric(m[, j])), cols)
}

#' @rdname fb_as_draws_simple
#' @keywords internal
#' @export
fb_as_draws_simple.flexybayes_inla <- function(fit, n_samples = 1000L, ...) {
  .check_installed(
    "INLA",
    "Package 'INLA' is required to extract draws from a ",
    "flexybayes_inla object."
  )
  if (is.null(fit$inla)) {
    stop("Cannot extract draws: fit$inla is missing.", call. = FALSE)
  }

  n_samples <- as.integer(n_samples)
  # Force single-threaded sampling, matching the summary path in
  # methods_inla.R. inla.posterior.sample() otherwise spawns one process
  # per core, which fails under the two-core limit `R CMD check --as-cran`
  # enforces (the worker count is reported as "N simultaneous processes
  # spawned").
  samples <- tryCatch(
    INLA::inla.posterior.sample(n_samples, fit$inla, num.threads = "1:1"),
    error = function(e) {
      stop(
        "INLA::inla.posterior.sample() failed: ",
        conditionMessage(e),
        ". Re-fit with control.compute = list(config = TRUE).",
        call. = FALSE
      )
    }
  )

  if (length(samples) == 0L) {
    stop("INLA::inla.posterior.sample() returned an empty list.", call. = FALSE)
  }

  latent_names <- rownames(samples[[1]]$latent)
  hyperpar_names <- names(samples[[1]]$hyperpar)

  out <- list()
  for (pn in latent_names) {
    out[[pn]] <- vapply(
      samples,
      function(s) as.numeric(s$latent[pn, 1]),
      numeric(1)
    )
  }
  for (pn in hyperpar_names) {
    out[[pn]] <- vapply(
      samples,
      function(s) as.numeric(s$hyperpar[pn]),
      numeric(1)
    )
  }
  out
}

#' @rdname fb_as_draws_simple
#' @keywords internal
#' @export
fb_as_draws_simple.flexybayes_brms <- function(fit, ...) {
  .check_installed(
    "posterior",
    "Package 'posterior' is required to extract draws from a ",
    "flexybayes_brms object."
  )
  if (is.null(fit$brms)) {
    stop("Cannot extract draws: fit$brms is missing.", call. = FALSE)
  }
  m <- tryCatch(
    as.matrix(posterior::as_draws_matrix(fit$brms)),
    error = function(e) {
      stop(
        "posterior::as_draws_matrix() failed on the ",
        "brms fit: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  cols <- colnames(m)
  if (is.null(cols)) {
    cols <- paste0("V", seq_len(ncol(m)))
  }
  setNames(lapply(seq_len(ncol(m)), function(j) as.numeric(m[, j])), cols)
}

#' @rdname fb_as_draws_simple
#' @keywords internal
#' @export
fb_as_draws_simple.default <- function(fit, ...) {
  stop(
    "triangulate() / fb_as_draws_simple() does not know how to ",
    "extract draws from an object of class ",
    paste(class(fit), collapse = "/"),
    ". Define an `fb_as_draws_simple.<class>` method.",
    call. = FALSE
  )
}

# ---------------------------------------------------------------- #
# Print method                                                     #
# ---------------------------------------------------------------- #

#' Print method for triangulate_result
#'
#' Internal S3 method. Brief summary header followed by the per-
#' parameter metrics table.
#'
#' @param x   A `triangulate_result` object as returned by
#'   [triangulate()].
#' @param ... Ignored. Present for compatibility with the generic.
#' @returns Invisibly, `x` unchanged. Called for the header and metrics
#'   table it prints.
#' @keywords internal
#' @export
print.triangulate_result <- function(x, ...) {
  cat("<triangulate_result>\n")
  cat("  status:   ", x$status %||% "unknown", "\n", sep = "")
  for (r in x$status_reasons %||% character(0)) {
    cat("    - ", r, "\n", sep = "")
  }
  cat("  source_a: ", x$source_a, "\n", sep = "")
  cat("  source_b: ", x$source_b, "\n", sep = "")
  # The independence axis labels the kind of convergence
  # claim this pair underwrites. Absent (empty) for same-backend or
  # unregistered pairs.
  if (length(x$independence) > 0L) {
    cat(
      "  independence: ",
      .format_independence_axes(x$independence),
      "\n",
      sep = ""
    )
    if (!is.null(x$axis_justification) && !is.na(x$axis_justification)) {
      cat("    (", x$axis_justification, ")\n", sep = "")
    }
  }
  cat("  n_common: ", x$n_common, "\n", sep = "")
  if (!is.null(x$n_compared)) {
    cat("  compared: ", x$n_compared, " of ", x$n_common,
        " (matched-prior overlap)\n", sep = "")
  }
  .print_only(x$only_a, "only_a")
  .print_only(x$only_b, "only_b")
  if (!is.null(x$shared_upstream_caveat) && !is.na(x$shared_upstream_caveat)) {
    cat("  [!] common-mode caveat: agreement does NOT test a shared upstream",
        "data fact\n")
    cat("      (data_independence = ",
        if (is.na(x$data_independence)) "undeclared" else
          as.character(x$data_independence), ")\n", sep = "")
  }

  if (x$n_common > 0L) {
    cat("\nMetrics (per common parameter):\n")
    print(.round_metrics(.triangulate_print_columns(x$metrics), 4))
    .print_not_compared(x$metrics)
  } else {
    cat(
      "\n  No common parameters. Supply `name_map` to align ",
      "fit_b's parameter names to fit_a's canonical names.\n",
      sep = ""
    )
  }
  invisible(x)
}

# Render the independence axes as a " + "-joined string.
# Colour grades match the methodological strength of the
# convergence claim -- algorithmic (strongest) blue, implementation
# plain, specification dim. Colour is additive: a monochrome / dumb
# terminal sees the same labels in plain text (cli degrades to identity
# off a dynamic TTY), so the label is never colour-load-bearing.
.format_independence_axes <- function(axes) {
  graded <- vapply(
    axes,
    function(a) {
      switch(
        a,
        algorithmic = cli::col_blue(a),
        specification = cli::col_grey(a),
        a
      )
    },
    character(1L)
  )
  paste(graded, collapse = " + ")
}

.print_only <- function(names_vec, label) {
  n <- length(names_vec)
  cat(
    "  ",
    label,
    ":   ",
    n,
    " parameter",
    if (n != 1L) "s" else "",
    if (n > 0L) {
      paste0(
        " (",
        paste(utils::head(names_vec, 5), collapse = ", "),
        if (n > 5L) ", ..." else "",
        ")"
      )
    } else {
      ""
    },
    "\n",
    sep = ""
  )
}

# The printed view. The full table stays on the object -- the tail drift
# and the raw Wasserstein-1 are there for anyone who wants them -- but a
# fourteen-column data.frame wraps into unreadability on a console, so the
# print method shows the columns a verdict is read from.
.triangulate_print_columns <- function(df) {
  keep <- c(
    "param",
    "mean_a",
    "mean_b",
    "mean_shift_sd",
    "sd_ratio",
    "w1_sd",
    "verdict"
  )
  keep <- keep[keep %in% names(df)]
  df[, keep, drop = FALSE]
}

# Parameters the gates excluded, with the reason. Printed under the table
# because an exclusion is a finding: it says the two engines were never
# asked the same question about that parameter.
.print_not_compared <- function(df) {
  if (is.null(df$verdict)) {
    return(invisible(NULL))
  }
  nc <- df[df$verdict == "not_compared", , drop = FALSE]
  if (nrow(nc) == 0L) {
    return(invisible(NULL))
  }
  cat("\nNot compared:\n")
  for (i in seq_len(nrow(nc))) {
    cat("  ", nc$param[[i]], " -- ", nc$reason[[i]], "\n", sep = "")
  }
  invisible(NULL)
}

.round_metrics <- function(df, digits = 4L) {
  num_cols <- vapply(df, is.numeric, logical(1))
  df[num_cols] <- lapply(df[num_cols], round, digits = digits)
  df
}
