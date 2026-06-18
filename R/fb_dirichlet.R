# fb_dirichlet.R -- Dirichlet family + compositional concentration fitter.
#
# The Dirichlet distribution is the natural model for compositional data:
# rows that are proportions summing to one (soil texture fractions, dietary
# composition, the relative abundance of species in a community, allele
# frequencies at a locus). Like the GEV it does not fit the single-column
# GLM-link mixed-model emit path -- the response is a K-column simplex, not a
# scalar with a mean-link -- so it is provided as a dedicated, self-contained
# fitter (`fb_dirichlet()`) following the genomic-fitter pattern, with a
# matching family-object constructor (`fb_family_dirichlet()`) registered
# against the family system so `family = "dirichlet"` is recognised and the
# user routed here rather than refused generically.
#
# The default method is dependency-free maximum likelihood (base `optim` on
# the log concentrations); an optional Bayesian route uses greta's native
# `dirichlet` distribution (already in the package's Suggests). Both routes
# return the same `fb_dirichlet_fit` shape.
# Reference: Minka (2000), *Estimating a Dirichlet distribution*.

# ---- Family-object constructor -----------------------------------------

#' Dirichlet family object
#'
#' Constructs the family descriptor for the Dirichlet distribution, the
#' natural model for compositional (simplex) data. The object mirrors the
#' shape of the other flexyBayes family descriptors -- a named list carrying
#' the family name, the parameter description, and the support -- while
#' signalling that the Dirichlet is fitted through the dedicated
#' `fb_dirichlet()` entry point rather than the GLM-link emit path.
#'
#' The Dirichlet on the \eqn{K}-simplex is parameterised by a vector of
#' \eqn{K} positive concentration parameters
#' \eqn{\alpha_1, \ldots, \alpha_K}. Their sum controls the concentration
#' (large sum gives compositions tightly clustered around the mean), and the
#' normalised vector \eqn{\alpha / \sum_k \alpha_k} is the mean composition.
#'
#' @returns An object of class `c("fb_family_dirichlet", "fb_family")`: a
#'   list with `family` (the canonical string `"dirichlet"`), `parameters`
#'   (`"concentration (alpha), one per simplex component"`), `support`
#'   (`"K-simplex"`), and `fitter` (`"fb_dirichlet"`).
#'
#' @seealso [fb_dirichlet()]
#' @examples
#' fam <- fb_family_dirichlet()
#' fam$family
#' @export
fb_family_dirichlet <- function() {
  structure(
    list(
      family = "dirichlet",
      parameters = "concentration (alpha), one per simplex component",
      support = "K-simplex",
      fitter = "fb_dirichlet"
    ),
    class = c("fb_family_dirichlet", "fb_family")
  )
}

# ---- Density / simulation primitives ------------------------------------

# Dirichlet log-density for a matrix `x` (rows = observations, columns =
# simplex components) under concentration vector `alpha`. Returns one
# log-density per row.
.ddirichlet_log <- function(x, alpha) {
  if (any(alpha <= 0)) {
    return(rep(-Inf, nrow(x)))
  }
  const <- lgamma(sum(alpha)) - sum(lgamma(alpha))
  const + colSums((alpha - 1) * t(log(x)))
}

#' Simulate from a Dirichlet distribution
#'
#' Draws compositional rows from a Dirichlet with the given concentration
#' vector, via the gamma-normalisation construction.
#'
#' @param n Integer. The number of compositions (rows) to draw.
#' @param alpha Numeric vector of positive concentration parameters; its
#'   length sets the number of simplex components \eqn{K}.
#'
#' @returns A numeric `n` by `K` matrix whose rows sum to one.
#'
#' @seealso [fb_dirichlet()]
#' @examples
#' set.seed(1)
#' X <- rdirichlet(5L, alpha = c(2, 5, 3))
#' rowSums(X)
#' @export
rdirichlet <- function(n, alpha) {
  if (!is.numeric(alpha) || length(alpha) < 2L || any(alpha <= 0)) {
    stop(
      "`alpha` must be a numeric vector of at least two positive values.",
      call. = FALSE
    )
  }
  k <- length(alpha)
  g <- matrix(
    stats::rgamma(n * k, shape = rep(alpha, each = n)),
    nrow = n
  )
  g / rowSums(g)
}

# ---- Fitter -------------------------------------------------------------

#' Fit a Dirichlet distribution to compositional data
#'
#' Estimates the concentration vector \eqn{\alpha} of a Dirichlet
#' distribution
#' from a matrix of compositional rows (each row a composition on the
#' simplex). The default method is maximum likelihood via base `optim()`
#' (dependency-free, deterministic); a Bayesian alternative uses greta's
#' native `dirichlet` distribution. Both routes return the same
#' `fb_dirichlet_fit` shape, so a `tidy()` summary is identical across
#' methods.
#'
#' Rows are renormalised to sum to one, and exact zeros or ones are nudged
#' inside the open simplex by a small `eps` (a Dirichlet density is `-Inf` on
#' the simplex boundary). The maximum-likelihood route optimises the log
#' concentrations to keep them positive and recovers standard errors from the
#' Hessian by the delta method. The Bayesian route places independent
#' half-normal priors on the concentrations and reports posterior means with
#' quantile-based credible intervals.
#'
#' The fitted mean composition (`alpha / sum(alpha)`) is reported on
#' `fit$mean_composition`.
#'
#' @param x A numeric matrix or `data.frame` of compositional rows: at least
#'   two columns, all entries non-negative, with at least four rows.
#' @param method Character. `"ml"` for maximum likelihood (the default) or
#'   `"greta"` for the Bayesian route.
#' @param conf_level Numeric in `(0, 1)`. The interval level for the parameter
#'   summary. Defaults to `0.95`.
#' @param eps Numeric. The boundary nudge applied to exact zeros / ones.
#'   Defaults to `1e-06`.
#' @param n_samples,warmup,chains Integer MCMC controls for the `"greta"`
#'   method (ignored by `"ml"`). Default to `2000`, `1000`, and `4`.
#' @param verbose Logical. Whether the greta route prints sampler progress.
#'   Defaults to `FALSE`.
#'
#' @returns An object of class `c("fb_dirichlet_fit", "fb_family_fit")`: a
#'   list with `estimates` (a `data.frame` of `term` / `estimate` /
#'   `std.error` / `conf.low` / `conf.high`, one row per component),
#'   `mean_composition` (numeric, summing to one), `method`, `n_obs`,
#'   `n_components`, `logLik` (ML only), and `draws` (greta only).
#'
#' @seealso [fb_family_dirichlet()], [rdirichlet()], [tidy.fb_dirichlet_fit()]
#' @examples
#' set.seed(1)
#' X <- rdirichlet(300L, alpha = c(2, 5, 3))
#' fit <- fb_dirichlet(X)
#' fit$estimates
#' @export
fb_dirichlet <- function(
  x,
  method = c("ml", "greta"),
  conf_level = 0.95,
  eps = 1e-06,
  n_samples = 2000L,
  warmup = 1000L,
  chains = 4L,
  verbose = FALSE
) {
  method <- match.arg(method)
  x <- .check_composition(x, eps)
  labels <- colnames(x)

  fit <- if (method == "ml") {
    .fb_dirichlet_ml(x, conf_level, labels)
  } else {
    .fb_dirichlet_greta(
      x, conf_level, labels, n_samples, warmup, chains, verbose
    )
  }

  alpha <- fit$estimates$estimate
  fit$mean_composition <- stats::setNames(alpha / sum(alpha), labels)
  fit$method <- method
  fit$n_obs <- nrow(x)
  fit$n_components <- ncol(x)
  structure(fit, class = c("fb_dirichlet_fit", "fb_family_fit"))
}

# Maximum-likelihood Dirichlet fit. Optimises log(alpha) so the
# concentrations stay positive, and delta-corrects the Hessian SEs back to
# the natural scale.
.fb_dirichlet_ml <- function(x, conf_level, labels) {
  nll <- function(log_alpha) {
    alpha <- exp(log_alpha)
    ll <- .ddirichlet_log(x, alpha)
    if (any(!is.finite(ll))) {
      return(1e10)
    }
    -sum(ll)
  }

  start <- log(colMeans(x) * 10 + 0.5)
  opt <- suppressWarnings(stats::optim(
    start,
    nll,
    method = "BFGS",
    hessian = TRUE
  ))
  if (!isTRUE(opt$convergence == 0L)) {
    warning(
      "fb_dirichlet(): the maximum-likelihood optimiser did not converge ",
      "(convergence code ", opt$convergence, "). Treat the estimates as ",
      "provisional and inspect the data.",
      call. = FALSE
    )
  }

  alpha <- exp(opt$par)
  se_log <- .safe_hessian_se(opt$hessian)
  # Delta method: SE(alpha) = SE(log_alpha) * d exp(log_alpha) = SE * alpha.
  se_natural <- se_log * alpha

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  estimates <- data.frame(
    term = labels,
    estimate = alpha,
    std.error = se_natural,
    conf.low = pmax(0, alpha - z * se_natural),
    conf.high = alpha + z * se_natural,
    stringsAsFactors = FALSE
  )

  list(estimates = estimates, logLik = -opt$value)
}

# Bayesian Dirichlet fit via greta's native `dirichlet` distribution. The
# concentration vector is a positive variable; the simplex rows are observed
# as `dirichlet(alpha, n_realisations = nrow(x))`.
.fb_dirichlet_greta <- function(
  x,
  conf_level,
  labels,
  n_samples,
  warmup,
  chains,
  verbose
) {
  if (!requireNamespace("greta", quietly = TRUE)) {
    stop(
      "fb_dirichlet(method = \"greta\") needs the greta package. Install it ",
      "with install.packages(\"greta\", repos = c(\"https://greta-dev.r-univ",
      "erse.dev\", getOption(\"repos\"))) and run greta::install_greta_deps",
      "(), or use method = \"ml\".",
      call. = FALSE
    )
  }

  draws <- .dirichlet_greta_sample(x, n_samples, warmup, chains, verbose)
  probs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
  est <- vapply(
    seq_len(ncol(draws)),
    function(j) {
      col <- draws[, j]
      c(
        mean(col),
        stats::sd(col),
        stats::quantile(col, probs[1L]),
        stats::quantile(col, probs[2L])
      )
    },
    numeric(4L)
  )
  estimates <- data.frame(
    term = labels,
    estimate = est[1L, ],
    std.error = est[2L, ],
    conf.low = est[3L, ],
    conf.high = est[4L, ],
    stringsAsFactors = FALSE
  )

  list(estimates = estimates, draws = draws)
}

# Build and sample the greta Dirichlet model. Kept separate so the
# parameterisation is testable in isolation.
.dirichlet_greta_sample <- function(x, n_samples, warmup, chains, verbose) {
  ge <- asNamespace("greta")
  variable <- get("variable", envir = ge)
  dirichlet_dist <- get("dirichlet", envir = ge)
  as_data <- get("as_data", envir = ge)
  mcmc <- get("mcmc", envir = ge)
  model <- get("model", envir = ge)
  distribution_setter <- get("distribution<-", envir = ge)

  k <- ncol(x)
  alpha <- variable(lower = 0, dim = c(1L, k))
  x_data <- as_data(x)
  distribution_setter(x_data, value = dirichlet_dist(
    alpha,
    n_realisations = nrow(x)
  ))
  m <- model(alpha)
  sims <- mcmc(
    m,
    n_samples = n_samples,
    warmup = warmup,
    chains = chains,
    verbose = isTRUE(verbose)
  )
  mat <- as.matrix(sims)
  # greta names the alpha columns alpha[1,1], alpha[1,2], ...; order is
  # column-major and stable, so the K trailing columns are the K
  # concentrations in component order.
  mat[, seq_len(k), drop = FALSE]
}

# ---- tidy / print -------------------------------------------------------

#' Tidy a Dirichlet fit
#'
#' Returns the concentration-parameter summary as a `broom`-style
#' `data.frame`, one row per simplex component, with the canonical columns.
#'
#' @param x An `fb_dirichlet_fit` object from [fb_dirichlet()].
#' @param ... Currently unused; present for generic compatibility.
#'
#' @returns A `data.frame` with `term`, `estimate`, `std.error`, `conf.low`,
#'   and `conf.high`.
#'
#' @seealso [fb_dirichlet()]
#' @examples
#' set.seed(1)
#' fit <- fb_dirichlet(rdirichlet(300L, c(2, 5, 3)))
#' tidy(fit)
#' @export
tidy.fb_dirichlet_fit <- function(x, ...) {
  x$estimates
}

#' Print a Dirichlet fit
#'
#' @param x An `fb_dirichlet_fit` object.
#' @param ... Currently unused; present for generic compatibility.
#' @returns `x`, invisibly.
#' @export
print.fb_dirichlet_fit <- function(x, ...) {
  cat("Dirichlet (compositional) fit  [flexyBayes]\n")
  cat(strrep("-", 50L), "\n")
  cat(
    "  Method :", x$method, "  N =", x$n_obs,
    "  components =", x$n_components, "\n"
  )
  est <- x$estimates
  for (i in seq_len(nrow(est))) {
    cat(sprintf(
      "  alpha[%-10s] %8.3f  (%.3f, %.3f)\n",
      est$term[i],
      est$estimate[i],
      est$conf.low[i],
      est$conf.high[i]
    ))
  }
  cat(
    "  Mean composition:",
    paste(sprintf("%.3f", x$mean_composition), collapse = " "),
    "\n"
  )
  cat(strrep("-", 50L), "\n")
  invisible(x)
}

# ---- Validation helper --------------------------------------------------

# Validate and clean a compositional matrix: numeric, >= 2 columns, >= 4
# rows, non-negative, renormalised to sum to one, with boundary nudging.
.check_composition <- function(x, eps) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(
      "`x` must be a numeric matrix or data.frame of compositional rows.",
      call. = FALSE
    )
  }
  if (ncol(x) < 2L) {
    stop(
      "`x` must have at least two columns (simplex components); got ",
      ncol(x), ".",
      call. = FALSE
    )
  }
  if (nrow(x) < 4L) {
    stop(
      "`x` must have at least four rows to identify the concentrations; got ",
      nrow(x), ".",
      call. = FALSE
    )
  }
  if (any(!is.finite(x)) || any(x < 0)) {
    stop("`x` must contain only finite, non-negative values.", call. = FALSE)
  }

  labels <- colnames(x)
  if (is.null(labels)) {
    labels <- paste0("c", seq_len(ncol(x)))
  }
  # Renormalise rows, then nudge off the boundary and renormalise again so
  # the Dirichlet density is finite.
  rs <- rowSums(x)
  if (any(rs <= 0)) {
    stop("`x` has a row that sums to zero; cannot normalise.", call. = FALSE)
  }
  x <- x / rs
  x <- pmin(pmax(x, eps), 1 - eps)
  x <- x / rowSums(x)
  colnames(x) <- labels
  x
}
