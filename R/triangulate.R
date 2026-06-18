# triangulate -- cross-engine posterior comparison.
#
# The signature feature of flexyBayes. Takes two model fits produced
# via different backends (greta / INLA / brms passthrough) and
# computes per-parameter posterior comparison metrics: Wasserstein-1
# distance (1D), R-hat-on-means, tail drift (0.025 / 0.975 quantile
# difference), and SD ratio.
#
# Backends use different internal parameter naming conventions:
#   - greta: e.g., "mu_atg", "beta_x", "tau_env[1]", "sigma_geno"
#   - INLA: e.g., "(Intercept):1", "x:1", "g:1", "Precision for g"
# Users supply a `name_map` (named character vector or list) to
# align fit_b's parameter names to fit_a's canonical names. Without
# alignment, only literal-match parameters are compared.
#
# Internal posterior extraction is via the fb_as_draws_simple()
# generic with methods for `flexybayes` (greta backend) and
# `flexybayes_inla` (INLA backend). Both return a named list whose
# values are numeric vectors of posterior draws.

# ---------------------------------------------------------------- #
# Public-facing entry                                              #
# ---------------------------------------------------------------- #

#' Cross-engine posterior triangulation
#'
#' Compute per-parameter posterior comparison metrics across two
#' Bayesian fits produced via different backends (greta / INLA /
#' brms passthrough). The signature feature of flexyBayes.
#'
#' For each parameter present in both fits (post `name_map`), the
#' returned table reports: posterior means, posterior SDs, Q2.5 /
#' Q97.5 differences (tail drift), Wasserstein-1 distance (1D
#' empirical), the SD ratio, and an R-hat-on-means scalar that
#' compares between-engine vs within-engine variance. Parameters
#' present in only one fit are reported in `only_a` / `only_b`.
#'
#' Backends use different parameter naming conventions; supply
#' `name_map = c(<fit_b name> = <canonical name>, ...)` to align
#' them. Without alignment, only literal-match parameter names are
#' compared.
#'
#' Backends also use different *parameter scales* for variance
#' components -- INLA reports precision (`Precision for g`), greta
#' reports standard deviation (`sigma_g`). Supply
#' `transform_a` / `transform_b` -- named lists of one-argument
#' functions keyed by parameter name -- to put the two posteriors on
#' a common scale before comparison. Transforms are applied first;
#' then `name_map` aligns the (already-transformed) fit_b names to
#' fit_a's canonical names. Names in `transform_b` therefore refer to
#' fit_b's *original* parameter names, not the post-`name_map`
#' canonical names.
#'
#' @section Matched priors:
#' When one of the inputs is an aggregated-gaussian fit (cell-level
#' sufficient statistics rather than the per-row likelihood), the
#' posteriors being compared are only directly comparable when the two
#' fits share priors. The aggregated path combines the cell-mean
#' likelihood with a precision prior carrying a closed-form correction
#' that absorbs the within-cell sum-of-squares; under the *default*
#' prior this recovers the per-row posterior to numerical precision, so
#' the aggregated fit is tagged `prior_parametrization =
#' "per_row_equivalent"` (visible in the aggregated `print()` / `summary()`
#' and in [canonical_names()]). When an explicit prior is supplied the
#' fit is tagged `"custom"`: the equivalence against a *default-prior*
#' per-row fit no longer holds, and on the aggregated INLA path the
#' observation-precision prior is not plumbed through, so a custom
#' residual prior is silently not applied there. Before reading the
#' agreement metrics on a custom-prior aggregated fit, confirm both
#' inputs carry the same prior with [prior_summary()].
#'
#' @param fit_a a fit object with a `fb_as_draws_simple` method
#'   (e.g., `flexybayes` from greta, `flexybayes_inla` from INLA).
#' @param fit_b a second fit object, typically from a different
#'   backend.
#' @param name_map named character vector or list mapping fit_b's
#'   parameter names (left) to canonical names matching fit_a
#'   (right). Optional.
#' @param transform_a,transform_b named list of functions, each
#'   taking a numeric vector of posterior draws and returning a
#'   numeric vector of the same length. Names key parameters in
#'   `fit_a` / `fit_b` (using each fit's *original* parameter
#'   names). Common use: pass
#'   `transform_b = list("Precision for g" = function(x) 1 / sqrt(x))`
#'   to convert INLA's precision draws to standard-deviation scale
#'   so they line up with greta's `sigma_g`. Optional.
#' @param n_samples integer: number of posterior samples to draw
#'   for fits whose extractor needs sampling (e.g., INLA via
#'   `INLA::inla.posterior.sample`).
#' @param data_independence single logical declaring whether the two fits were
#'   built on independently-sourced data. `triangulate()` measures inter-fit
#'   *agreement*, and the backend-independence registry certifies code (not
#'   data) independence -- so if both fits share the same upstream data, a
#'   fabricated data fact is common-mode and their agreement cannot detect it.
#'   `TRUE` declares the data independent (no caveat); `FALSE` (same data) or
#'   `NA` (the default, undeclared) attach a `shared_upstream_caveat` field to
#'   the result, surfaced prominently by the print method, so agreement is never
#'   silently mistaken for corroboration of a shared data fact (Independent
#'   Oracle Principle).
#' @return a `triangulate_result` S3 object (list). Key fields:
#'   `metrics` (data.frame, one row per common parameter),
#'   `common` (character), `only_a`, `only_b`, `n_common`,
#'   `source_a`, `source_b`.
#' @examples
#' # Live INLA posterior sampling can fail in restricted-process
#' # check environments (the `inla.posterior.sample` parallelism
#' # check trips); the example uses `\dontrun{}` deliberately. On
#' # an interactive install with greta + INLA + sn available it
#' # runs in a few seconds.
#' \dontrun{
#' if (requireNamespace("greta", quietly = TRUE) &&
#'     requireNamespace("INLA",  quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(y = rnorm(40), x = rnorm(40),
#'                   g = factor(rep(1:5, 8)))
#'   fit_g <- fb(y ~ x + (1 | g), data = d, backend = "greta",
#'               n_samples = 100, warmup = 100, chains = 1,
#'               verbose = FALSE)
#'   fit_i <- fb(y ~ x + (1 | g), data = d, backend = "inla",
#'               verbose = FALSE)
#'   prec_to_sd <- function(x) 1 / sqrt(x)
#'   tri <- triangulate(
#'     fit_g, fit_i,
#'     transform_b = list(
#'       "Precision for g"                         = prec_to_sd,
#'       "Precision for the Gaussian observations" = prec_to_sd
#'     ),
#'     name_map = c(
#'       "(Intercept):1"                           = "mu_atg",
#'       "x:1"                                     = "beta_x",
#'       "Precision for g"                         = "sigma_g",
#'       "Precision for the Gaussian observations" = "sigma_e_atg"
#'     )
#'   )
#'   print(tri)
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
      common = common,
      only_a = only_a,
      only_b = only_b,
      n_common = length(common),
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
    stringsAsFactors = FALSE
  )
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
#' numeric posterior-draw vectors from each fit. Methods exist for
#' the `flexybayes` (greta backend) and `flexybayes_inla` (INLA
#' backend) classes; user-defined methods can extend the generic.
#'
#' @param fit a model fit object.
#' @param ... method-specific arguments (e.g., `n_samples` for INLA).
#' @return a named list of numeric vectors.
#' @export
fb_as_draws_simple <- function(fit, ...) UseMethod("fb_as_draws_simple")

#' @rdname fb_as_draws_simple
#' @keywords internal
#' @export
fb_as_draws_simple.flexybayes <- function(fit, ...) {
  if (is.null(fit$greta) || is.null(fit$greta$draws)) {
    stop(
      "Cannot extract draws: fit$greta$draws is missing. ",
      "Did the fit complete (return_code = FALSE) and use ",
      "backend = \"greta\"?",
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
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop(
      "Package 'INLA' is required to extract draws from a ",
      "flexybayes_inla object.",
      call. = FALSE
    )
  }
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
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop(
      "Package 'posterior' is required to extract draws from a ",
      "flexybayes_brms object.",
      call. = FALSE
    )
  }
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
#' @param x   a `triangulate_result` object.
#' @param ... unused.
#' @return invisibly returns `x`.
#' @keywords internal
#' @export
print.triangulate_result <- function(x, ...) {
  cat("<triangulate_result>\n")
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
    print(.round_metrics(x$metrics, 4))
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

.round_metrics <- function(df, digits = 4L) {
  num_cols <- vapply(df, is.numeric, logical(1))
  df[num_cols] <- lapply(df[num_cols], round, digits = digits)
  df
}
