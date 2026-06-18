# fb_log_posterior() -- exposes a fitted model's log-posterior as a callable.
#
# This is the one inference-result outflow from flexyBayes: a fitted
# flexyBayes object is turned into a vectorised, domain-safe, unnormalised
# log-posterior callable that a downstream tool (for example the proxymix
# package, via `proxymix::from_fb_posterior()`) can compress into a
# closed-form Gaussian-mixture proxy. flexyBayes owns the producer; the
# consumer never depends on flexyBayes at run time (an acyclic dependency
# invariant), so this file emits the *callable*, read by the consumer as a
# bare callable carrying its metadata as attributes.
#
# The greta backend is the canonical real producer: greta retains the
# model graph (`fit$greta$model`), whose directed acyclic graph exposes an
# unadjusted joint-density function. Evaluated at a model's free
# (unconstrained) state, that function returns
# `log p(data | theta) + log p(theta)` on the *natural* (constrained) scale
# -- exactly the unnormalised log-posterior the contract requires, with
# no Jacobian adjustment (the adjustment belongs to the free-state
# parameterisation, not the natural-scale density). This producer therefore
# maps natural-scale inputs to the free state, calls the dag's unadjusted
# log-prob, and returns the result.
#
# The brms and INLA backends abstain. brms's log-density (`rstan::log_prob`)
# lives on the Stan unconstrained scale, and the parameter-name to
# Stan-unconstrained-index mapping is version-fragile (transformed
# parameters, declaration order); a wrong mapping would yield a wrong
# log-density silently, so an honest abstain beats a plausible-but-wrong
# producer (Independent Oracle Principle). INLA's posterior is a
# deterministic Laplace/grid approximation, not a sampling log-density at
# arbitrary natural-scale points, so it abstains too. Both raise an
# informative, classed condition rather than guessing.

# ---------------------------------------------------------------- #
# The generic                                                      #
# ---------------------------------------------------------------- #

#' Emit a flexyBayes posterior as a log-density producer
#'
#' Turns a fitted flexyBayes object into a log-posterior producer:
#' a vectorised, domain-safe, unnormalised log-posterior callable that
#' `proxymix::from_fb_posterior()` compresses into a closed-form
#' Gaussian-mixture proxy. It is the single inference-result outflow from
#' flexyBayes; the contract is the *log-density*, not the draws, so the
#' returned object is addressed purely through its callable.
#'
#' The returned value is a **bare callable**
#' `function(theta_matrix) -> numeric`. Its input is a numeric matrix whose
#' rows index independent parameter draws and whose columns index
#' parameters, in `attr(., "parameter_names")` order, on the natural
#' (constrained) scale. Its output is a length-`nrow(theta_matrix)` numeric
#' vector of `log p(theta | data) + const` (unnormalised). The callable is
#' vectorised, side-effect free, and domain-safe: a row outside the
#' parameters' support returns `-Inf` rather than raising an error (the
#' consumer probes it at construction). It carries, as attributes:
#'
#' \describe{
#'   \item{`parameter_names`}{Character vector naming the parameters; its
#'     length fixes the proxy's ambient dimension. Vector-valued targets
#'     are flattened in column-major order with index suffixes (e.g.
#'     `beta[1,1]`, `beta[2,1]`).}
#'   \item{`log_normalizer`}{The additive correction that would normalise
#'     the density, i.e. `-log Z`. For a posterior the marginal likelihood
#'     is generally unknown, so this is `NA_real_` -- honest, and the
#'     consumer reports a shifted (not absolute) divergence.}
#'   \item{`support_lower`, `support_upper`}{Length-`n_dim` numeric support
#'     bounds taken from the model's parameter constraints (`NA` for an
#'     unbounded coordinate). A variance / scale parameter is bounded below
#'     by zero, for instance. Used only to centre and scale the consumer's
#'     default importance proposal.}
#'   \item{`draws`}{An `n` by `n_dim` numeric matrix of the fit's posterior
#'     draws on the natural scale, column-aligned to `parameter_names`.
#'     Used only to seed the consumer's default proposal; never required.}
#' }
#'
#' Backend support. The **greta** backend is the canonical real producer:
#' it evaluates the model graph's unadjusted joint density at the
#' free-state image of the supplied natural-scale parameters, which is the
#' unnormalised natural-scale log-posterior exactly. The **brms** and
#' **INLA** backends abstain with an informative condition -- brms's
#' log-density lives on the Stan unconstrained scale with a version-fragile
#' name mapping, and INLA's posterior is a deterministic approximation, not
#' a sampling log-density; an honest abstain is preferred to a
#' plausible-but-wrong log-density.
#'
#' Acyclic note. A consumer such as proxymix uses this callable without
#' depending on flexyBayes; flexyBayes does not list proxymix in `Imports`
#' or `Suggests`. The cross-package demonstration lives in a separate
#' integration harness, not in this package, preserving the acyclic
#' dependency graph.
#'
#' @param fit A fitted flexyBayes object. The greta classes
#'   (`flexybayes` / `flexybayes_direct_greta`) produce a real callable;
#'   the brms and INLA classes abstain.
#' @param ... Reserved for future producer options; currently unused.
#'
#' @return For a greta fit, a bare callable `function(theta_matrix)` with
#'   the attributes described above, ready to pass to
#'   `proxymix::from_fb_posterior()`. For a brms or INLA fit, the function
#'   does not return: it raises a classed `fb_c4_unavailable` condition.
#'
#' @family flexyBayes interop
#' @seealso `proxymix::from_fb_posterior()` for the consumer (compresses
#'   the returned callable into a Gaussian-mixture proxy).
#' @examples
#' \dontrun{
#' library(greta)
#' n <- 30
#' y <- rnorm(n, 1.5, 2)
#' mu <- normal(0, 5)
#' sigma <- normal(0, 5, truncation = c(0, Inf))
#' yd <- as_data(y)
#' distribution(yd) <- normal(mu, sigma)
#' m <- model(mu, sigma)
#' fit <- fb_greta(fb_from_greta(m), n_samples = 500, warmup = 500,
#'                 chains = 2, verbose = FALSE, mcmc_verbose = FALSE)
#' producer <- fb_log_posterior(fit)
#' attr(producer, "parameter_names")
#' producer(matrix(c(1.5, 2.0), nrow = 1)) # natural-scale log-posterior
#' ## Compress with proxymix (in a separate integration harness):
#' ## proxymix::from_fb_posterior(producer, N = 2)
#' }
#' @export
fb_log_posterior <- function(fit, ...) {
  UseMethod("fb_log_posterior")
}

# ---------------------------------------------------------------- #
# Abstain: default + non-greta backends                            #
# ---------------------------------------------------------------- #

#' @rdname fb_log_posterior
#' @export
fb_log_posterior.default <- function(fit, ...) {
  .fb_c4_abstain(
    backend = "this object class",
    detail = paste0(
      "fb_log_posterior() produces a C4 log-density only for a fitted ",
      "flexyBayes object. Got class: ",
      paste(class(fit), collapse = "/"),
      "."
    )
  )
}

#' @rdname fb_log_posterior
#' @export
fb_log_posterior.flexybayes_brms <- function(fit, ...) {
  .fb_c4_abstain(
    backend = "the brms backend",
    detail = paste0(
      "C4 log-density producer not available for the brms backend in ",
      "this version. brms's log-density (rstan::log_prob) is defined on ",
      "the Stan unconstrained scale, and the mapping from flexyBayes ",
      "parameter names to Stan unconstrained indices is version-fragile ",
      "(transformed parameters, declaration order). An honest abstain is ",
      "preferred to a plausible-but-wrong log-density. Refit the model ",
      "via fb_greta() to obtain a C4 producer."
    )
  )
}

#' @rdname fb_log_posterior
#' @export
fb_log_posterior.flexybayes_inla <- function(fit, ...) {
  .fb_c4_abstain(
    backend = "the INLA backend",
    detail = paste0(
      "C4 log-density producer not available for the INLA backend. INLA's ",
      "posterior is a deterministic Laplace / grid approximation, not a ",
      "sampling log-density evaluable at arbitrary natural-scale points, ",
      "so there is no faithful unnormalised log-posterior to emit. Refit ",
      "the model via fb_greta() to obtain a C4 producer."
    )
  )
}

# Raise the classed abstention condition. Never returns.
.fb_c4_abstain <- function(backend, detail) {
  cond <- structure(
    class = c("fb_c4_unavailable", "error", "condition"),
    list(
      message = paste0(
        "No C4 log-density producer is available for ", backend, ".\n",
        detail
      ),
      call = sys.call(-1L)
    )
  )
  stop(cond)
}

# ---------------------------------------------------------------- #
# Real producer: greta backend                                     #
# ---------------------------------------------------------------- #

#' @rdname fb_log_posterior
#' @export
fb_log_posterior.flexybayes <- function(fit, ...) {
  if (!requireNamespace("greta", quietly = TRUE)) {
    stop(
      "Package 'greta' is required to build a C4 log-density producer ",
      "from a greta-backed flexyBayes fit. Install with:\n",
      "  install.packages('greta')",
      call. = FALSE
    )
  }

  greta_slot <- fit$greta
  if (is.null(greta_slot) || is.null(greta_slot$model)) {
    .fb_c4_abstain(
      backend = "this flexyBayes fit",
      detail = paste0(
        "The fit does not carry a greta model graph on `$greta$model`. ",
        "fb_log_posterior() produces a C4 log-density only for greta-backed ",
        "fits (fb_greta(...) or a native greta_model). For brms / INLA ",
        "fits the producer abstains by design."
      )
    )
  }
  model <- greta_slot$model
  if (!inherits(model, "greta_model")) {
    .fb_c4_abstain(
      backend = "this flexyBayes fit",
      detail = paste0(
        "`$greta$model` is not a `greta_model` (class: ",
        paste(class(model), collapse = "/"),
        "). A C4 producer needs the retained greta graph."
      )
    )
  }

  # The model's directed acyclic graph and its ordered target nodes. The
  # free-state column order is the target-node order, with each node
  # flattened in column-major order; this is also the column order of the
  # posterior draws. Both are derived once, here, from the same dag.
  dag <- model$dag
  if (is.null(dag) || !is.environment(dag)) {
    .fb_c4_abstain(
      backend = "this flexyBayes fit",
      detail = paste0(
        "`$greta$model$dag` is NULL or not an environment; the greta graph ",
        "is unavailable. Confirm the fit was produced by a recent greta ",
        "(validated against greta 0.5.x)."
      )
    )
  }
  target_nodes <- dag$target_nodes
  if (is.null(target_nodes) || length(target_nodes) == 0L) {
    .fb_c4_abstain(
      backend = "this flexyBayes fit",
      detail = "The greta graph has no target nodes to infer over."
    )
  }

  # Per-target element metadata: flatten each node to its scalar
  # coordinates, recording the natural-scale name, lower / upper support
  # bound, and the natural -> free transform per coordinate.
  meta <- .fb_greta_param_meta(target_nodes)
  parameter_names <- meta$names
  # Consumer convention: NA marks an unbounded coordinate. The transform
  # maths below uses meta$lower / meta$upper (kept as +/-Inf).
  support_lower <- meta$support_lower
  support_upper <- meta$support_upper
  n_dim <- length(parameter_names)

  # The dag's unadjusted joint-density function, built once. Evaluated at a
  # free-state matrix (rows = independent draws), it returns one
  # natural-scale unnormalised log-posterior per row. greta needs its tf
  # environment built before this function is usable; a tiny calculate()
  # on the targets primes it without perturbing the fit.
  lp_fun <- .fb_greta_build_logprob(model, target_nodes)

  # Posterior draws on the natural scale, column-aligned to
  # parameter_names, to seed the consumer's default proposal. Best-effort:
  # the contract is the log-density, not the draws, so a failure here is
  # silent (draws = NULL).
  draws <- .fb_greta_draws_matrix(greta_slot$draws, parameter_names)

  # The producer callable. Domain-safe: rows outside support map to -Inf
  # without erroring; a non-finite free-state image (transform of an
  # out-of-support natural value) is likewise -Inf. Side-effect free.
  log_density <- function(theta_matrix) {
    if (is.null(dim(theta_matrix))) {
      theta_matrix <- matrix(theta_matrix, nrow = 1L)
    }
    theta_matrix <- as.matrix(theta_matrix)
    storage.mode(theta_matrix) <- "double"
    if (ncol(theta_matrix) != n_dim) {
      stop(
        "C4 log-density: expected ", n_dim, " columns (",
        paste(parameter_names, collapse = ", "), "), got ",
        ncol(theta_matrix), ".",
        call. = FALSE
      )
    }
    n_row <- nrow(theta_matrix)
    out <- rep(-Inf, n_row)
    if (n_row == 0L) {
      return(out)
    }

    # Natural -> free per coordinate; rows that fall outside any
    # coordinate's support (or whose free image is non-finite) are flagged
    # and left at -Inf, never sent to the graph.
    free <- matrix(NA_real_, nrow = n_row, ncol = n_dim)
    in_support <- rep(TRUE, n_row)
    for (j in seq_len(n_dim)) {
      # The transform deliberately maps out-of-support naturals to NaN
      # (e.g. log of a non-positive value); that is the domain-safety
      # signal, not an anomaly, so its warning is suppressed here.
      tj <- suppressWarnings(.fb_natural_to_free(
        theta_matrix[, j], meta$lower[j], meta$upper[j]
      ))
      bad <- !is.finite(tj)
      in_support[bad] <- FALSE
      free[, j] <- tj
    }

    if (!any(in_support)) {
      return(out)
    }
    free_ok <- free[in_support, , drop = FALSE]
    vals <- tryCatch(
      as.numeric(lp_fun(free_ok)),
      error = function(e) rep(NA_real_, nrow(free_ok))
    )
    # A graph that returns NaN / NA at an in-support point is treated as
    # zero density there (-Inf), never propagated as a numeric.
    vals[!is.finite(vals)] <- -Inf
    out[in_support] <- vals
    out
  }

  attr(log_density, "parameter_names") <- parameter_names
  attr(log_density, "log_normalizer") <- NA_real_
  attr(log_density, "support_lower") <- support_lower
  attr(log_density, "support_upper") <- support_upper
  attr(log_density, "draws") <- draws
  class(log_density) <- c("fb_log_posterior_producer", "function")
  log_density
}

#' @export
print.fb_log_posterior_producer <- function(x, ...) {
  pn <- attr(x, "parameter_names")
  dr <- attr(x, "draws")
  cat(sprintf(
    "<fb_log_posterior_producer> (C4 greta producer) in p = %d dimensions\n",
    length(pn)
  ))
  cat(sprintf("  parameters    : %s\n", paste(pn, collapse = ", ")))
  lo <- attr(x, "support_lower")
  hi <- attr(x, "support_upper")
  bounded <- (!is.null(lo) && any(is.finite(lo))) ||
    (!is.null(hi) && any(is.finite(hi)))
  cat(sprintf("  support bounds: %s\n",
              if (bounded) "supplied (from parameter constraints)"
              else "unbounded"))
  cat(sprintf("  log Z         : unknown (NA -- posterior marginal likelihood)\n"))
  cat(sprintf("  proposal seed : %s\n",
              if (!is.null(dr)) sprintf("%d posterior draws", nrow(dr))
              else "<none>"))
  invisible(x)
}

# ---------------------------------------------------------------- #
# Internal helpers -- greta graph introspection + transforms        #
# ---------------------------------------------------------------- #

# Flatten the ordered target nodes into per-coordinate metadata: the
# natural-scale parameter name (scalar nodes keep their name; vector /
# matrix nodes gain a column-major index suffix matching greta's draw
# column names, e.g. `beta[1,1]`), and the lower / upper support bound per
# coordinate. The order is the free-state / draw column order.
.fb_greta_param_meta <- function(target_nodes) {
  names_out <- character(0)
  lower_out <- numeric(0)
  upper_out <- numeric(0)
  node_names <- names(target_nodes)
  for (k in seq_along(target_nodes)) {
    nd <- target_nodes[[k]]
    nm <- node_names[k]
    d <- nd$dim
    n_elem <- if (is.null(d)) 1L else as.integer(prod(d))
    # greta stores per-element lower / upper; recycle a scalar bound.
    lo <- nd$lower
    hi <- nd$upper
    lo <- if (is.null(lo)) rep(-Inf, n_elem) else rep_len(as.numeric(lo), n_elem)
    hi <- if (is.null(hi)) rep(Inf, n_elem) else rep_len(as.numeric(hi), n_elem)

    if (n_elem == 1L) {
      coord_names <- nm
    } else {
      # Column-major index labels matching greta's draw column names. For a
      # 2-D target greta labels as `nm[i,j]`; for a 1-D vector greta still
      # carries the matrix form `nm[i,1]` (a greta_array is at least 2-D).
      dim_use <- if (is.null(d) || length(d) == 1L) c(n_elem, 1L) else d
      idx <- expand.grid(lapply(dim_use, seq_len), KEEP.OUT.ATTRS = FALSE)
      coord_names <- paste0(
        nm, "[", apply(idx, 1L, paste, collapse = ","), "]"
      )
    }
    names_out <- c(names_out, coord_names)
    lower_out <- c(lower_out, lo)
    upper_out <- c(upper_out, hi)
  }
  # Expose NA (not +/-Inf) for unbounded coordinates in the support
  # attributes -- that is the consumer's "unbounded" convention.
  support_lower <- lower_out
  support_upper <- upper_out
  support_lower[!is.finite(support_lower)] <- NA_real_
  support_upper[!is.finite(support_upper)] <- NA_real_
  list(
    names = names_out,
    lower = lower_out, # kept as +/-Inf for the transform maths
    upper = upper_out,
    support_lower = support_lower,
    support_upper = support_upper
  )
}

# Build the dag's unadjusted joint-density function, priming greta's tf
# environment first. The unadjusted density is the natural-scale
# log-posterior (no free-state Jacobian), which is what C4 wants.
.fb_greta_build_logprob <- function(model, target_nodes) {
  # Prime the tf environment via a no-op calculate() on the targets. This
  # is the public entry that builds the graph; without it the dag's
  # log-prob generator references an undefined tf environment.
  primed <- tryCatch({
    do.call(
      greta::calculate,
      c(unname(target_nodes), list(nsim = 1L))
    )
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(primed)) {
    # calculate() can refuse on some graphs; fall back to building the dag
    # directly. Either path leaves a usable log-prob function.
    tryCatch(model$dag$define_tf(), error = function(e) NULL)
  }
  model$dag$generate_log_prob_function("unadjusted")
}

# Reorder the posterior draws to the parameter_names (free-state) column
# order and return as a plain natural-scale numeric matrix. Best-effort:
# returns NULL if the draws are absent or cannot be aligned.
.fb_greta_draws_matrix <- function(draws, parameter_names) {
  if (is.null(draws)) {
    return(NULL)
  }
  mat <- tryCatch(as.matrix(draws), error = function(e) NULL)
  if (is.null(mat) || is.null(colnames(mat))) {
    return(NULL)
  }
  if (!all(parameter_names %in% colnames(mat))) {
    return(NULL)
  }
  out <- mat[, parameter_names, drop = FALSE]
  storage.mode(out) <- "double"
  if (anyNA(out)) {
    return(NULL)
  }
  unname(out)
}

# Natural -> free (unconstrained) transform for a single coordinate,
# vectorised over a numeric vector of natural values. Mirrors greta's
# default bijectors:
#   * both bounds infinite      -> identity
#   * lower finite, upper Inf   -> log(x - lower)
#   * lower -Inf, upper finite  -> log(upper - x)
#   * both finite               -> logit((x - lower) / (upper - lower))
# A natural value outside the open support maps to a non-finite free value
# (NaN / Inf), which the caller reads as "outside support" -> -Inf density.
.fb_natural_to_free <- function(x, lower, upper) {
  lo_fin <- is.finite(lower)
  hi_fin <- is.finite(upper)
  if (!lo_fin && !hi_fin) {
    return(x)
  }
  if (lo_fin && !hi_fin) {
    z <- x - lower
    out <- ifelse(z > 0, log(z), NaN)
    return(out)
  }
  if (!lo_fin && hi_fin) {
    z <- upper - x
    out <- ifelse(z > 0, log(z), NaN)
    return(out)
  }
  # Both finite: logit on the open interval.
  u <- (x - lower) / (upper - lower)
  ifelse(u > 0 & u < 1, stats::qlogis(u), NaN)
}
