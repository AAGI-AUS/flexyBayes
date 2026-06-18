# Genome-wide association scan (G3) --- the specialised fast path.
#
# A whole-genome scan is not "fit a Bayesian model per marker" (infeasible
# at 10^4-10^6 markers). It is the EMMAX / P3D mixed-model scan (Kang et
# al. 2010): fit the polygenic null mixed model once to estimate the
# variance components, then test each marker by generalised least squares
# under those fixed components. The shared spectral primitive (G0a) makes
# this an O(n) test per marker after a single O(n^3) eigendecomposition of
# the relationship matrix --- the rotated model has a diagonal residual
# covariance, so the per-marker GLS is an exact weighted least squares.
#
# This file is engine-agnostic base R: the scan needs no greta / INLA /
# brms backend (the backends enter only at the optional top-hit Bayesian
# refinement, where the genome has been reduced to a handful of loci and a
# full multi-backend fit is affordable). It composes `.fb_spectral()`,
# `.fb_reml_vc()` (REML variance components via the rotated restricted
# likelihood), and the vectorised EMMAX score test.

# --- REML variance components (EMMA-style) ------------------------ #

# `.fb_reml_vc()` --- restricted-maximum-likelihood estimates of the
# genetic and residual variances for the polygenic null mixed model
# y = X beta + g + e, g ~ N(0, sigma_g^2 K), e ~ N(0, sigma_e^2 I), via a
# one-dimensional optimisation of the rotated REML profile likelihood
# over the variance ratio delta = sigma_e^2 / sigma_g^2 (Kang et al. 2008).
#
# In the eigenbasis of K the covariance is diagonal,
# Cov(U'y) = sigma_g^2 diag(lambda_i + delta), so the REML log-likelihood
# is a fast scalar function of delta. Returns sigma_g^2, sigma_e^2, delta,
# the implied heritability 1 / (1 + delta), and the rotated quantities the
# scan reuses.
.fb_reml_vc <- function(spec, y, X, delta_log_range = c(-10, 10)) {
  .fb_spectral_check(spec)
  if (spec$rank != spec$n) {
    stop(
      "`.fb_reml_vc()` needs the full-rank spectral decomposition; got ",
      "rank ", spec$rank, " of ", spec$n, ".",
      call. = FALSE
    )
  }
  n <- spec$n
  X <- as.matrix(X)
  p <- ncol(X)
  if (length(y) != n || nrow(X) != n) {
    stop("`y` and `X` must have ", n, " rows (the dimension of K).",
      call. = FALSE
    )
  }
  if (n - p < 1L) {
    stop(
      "the relationship dimension (", n, ") must exceed the number of ",
      "fixed-effect columns (", p, ") for REML.",
      call. = FALSE
    )
  }

  y_star <- as.vector(.fb_spectral_rotate(spec, y))
  X_star <- .fb_spectral_rotate(spec, X)
  lambda <- spec$values

  # REML profile negative log-likelihood (constants dropped) as a
  # function of log(delta). The weighted design Xtil = X* / sqrt(d) and
  # ytil = y* / sqrt(d) put the problem in homoscedastic form.
  neg_ll <- function(log_delta) {
    delta <- exp(log_delta)
    d <- lambda + delta
    sq <- sqrt(d)
    ytil <- y_star / sq
    Xtil <- X_star / sq
    xtx <- crossprod(Xtil)
    beta <- tryCatch(
      solve(xtx, crossprod(Xtil, ytil)),
      error = function(e) NULL
    )
    if (is.null(beta)) {
      return(Inf)
    }
    r <- ytil - Xtil %*% beta
    rss <- sum(r^2)
    ld_xtx <- determinant(xtx, logarithm = TRUE)$modulus
    0.5 * ((n - p) * log(rss) + sum(log(d)) + as.numeric(ld_xtx))
  }

  opt <- stats::optimize(
    neg_ll, interval = delta_log_range, tol = 1e-8
  )
  delta <- exp(opt$minimum)
  d <- lambda + delta
  sq <- sqrt(d)
  ytil <- y_star / sq
  Xtil <- X_star / sq
  xtx <- crossprod(Xtil)
  beta <- solve(xtx, crossprod(Xtil, ytil))
  r <- ytil - Xtil %*% beta
  rss <- sum(r^2)
  var_g <- rss / (n - p)
  var_e <- delta * var_g

  list(
    var_g = var_g,
    var_e = var_e,
    delta = delta,
    h2 = 1 / (1 + delta),
    beta = as.vector(beta),
    y_star = y_star,
    X_star = X_star,
    weights = 1 / d,
    converged = is.finite(opt$objective)
  )
}

# --- the scan ----------------------------------------------------- #

#' Genome-wide association scan (EMMAX / P3D)
#'
#' Test each marker for association with a phenotype while correcting for
#' polygenic background and population / family structure via a genomic
#' relationship matrix. `fb_gwas()` fits the polygenic null mixed model
#' once to estimate the variance components by REML, then --- holding
#' those components fixed (the P3D approximation of Kang et al. 2010) ---
#' tests every marker by exact generalised least squares in the eigenbasis
#' of the relationship matrix. The single eigendecomposition is the shared
#' spectral primitive; each marker is then an `O(n)` weighted-least-squares
#' score test, so a whole-genome scan is feasible without a per-marker
#' model fit.
#'
#' The backends (greta / INLA / brms) are not used by the scan itself ---
#' it is a deterministic frequentist fast path. They enter only when a
#' handful of significant loci are re-fit as full Bayesian models for
#' credible effect sizes, which is affordable at that reduced scale.
#'
#' @param formula A two-sided formula `y ~ covariates` for the response
#'   and fixed-effect covariates (an intercept is included unless removed
#'   with `- 1`). Use `y ~ 1` for an intercept-only background.
#' @param data A data frame with the response and covariates; one row per
#'   individual (the association unit).
#' @param markers A numeric `n x m` matrix of marker genotypes (allele
#'   dosages), one row per individual aligned to `data`, one column per
#'   marker. Column names are used as marker identifiers.
#' @param K An `n x n` genomic (or pedigree) relationship matrix. When
#'   `NULL` (default) it is built from the centred markers as
#'   \eqn{G = Z_c Z_c^\top / m}, the average-allele-frequency genomic
#'   relationship (Astle & Balding 2009). This normalises by the marker
#'   count `m`, unlike VanRaden's (2008) first method, which divides by
#'   \eqn{2\sum_j p_j(1 - p_j)}; the two coincide only when every marker has
#'   minor-allele frequency 0.5 and otherwise differ by a constant scale that
#'   the REML variance ratio absorbs, so the scan statistics are unaffected.
#'   Note that building `K` from the same markers being tested causes *proximal
#'   contamination*: a marker's own signal enters the polygenic term and
#'   deflates its test statistic (conservative, but a power loss). For a
#'   well-calibrated scan supply `K` from an independent background panel
#'   or a leave-one-chromosome-out construction.
#' @param marker_map Optional data frame with columns `marker`, `chr`,
#'   `pos` to annotate the results and order the Manhattan plot.
#' @param tol PSD / symmetry tolerance passed to the spectral primitive.
#'
#' @return An `fb_gwas` object: a list with `results` (a data frame of
#'   `marker`, `effect`, `se`, `statistic` (the chi-square), `p_value`,
#'   `p_bonferroni`, `q_value` (Benjamini-Hochberg FDR), plus any map
#'   columns), `lambda_gc` (the genomic-control inflation factor),
#'   `var_g` / `var_e` / `h2` (the null REML components), and metadata.
#'
#' @references Kang, H. M., et al. (2010). Variance component model to
#'   account for sample structure in genome-wide association studies.
#'   \emph{Nature Genetics}, 42(4), 348-354. Astle, W., & Balding, D. J.
#'   (2009). Population structure and cryptic relatedness in genetic
#'   association studies. \emph{Statistical Science}, 24(4), 451-471.
#'   VanRaden, P. M. (2008). Efficient methods to compute genomic
#'   predictions. \emph{Journal of Dairy Science}, 91(11), 4414-4423.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 100L
#' m <- 200L
#' M <- matrix(rbinom(n * m, 2L, 0.3), n, m)
#' colnames(M) <- paste0("snp", seq_len(m))
#' y <- 2 * scale(M[, 50L]) + rnorm(n)
#' dat <- data.frame(y = y)
#' scan <- fb_gwas(y ~ 1, data = dat, markers = M)
#' head(scan$results[order(scan$results$p_value), ])
#' scan$lambda_gc
#' }
#' @export
fb_gwas <- function(
  formula,
  data,
  markers,
  K = NULL,
  marker_map = NULL,
  tol = 1e-8
) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a two-sided formula, e.g. y ~ 1.", call. = FALSE)
  }
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, mf)
  n <- length(y)

  markers <- as.matrix(markers)
  if (nrow(markers) != n) {
    stop(
      "`markers` must have one row per individual (", n, "); got ",
      nrow(markers), ".",
      call. = FALSE
    )
  }
  m <- ncol(markers)
  marker_ids <- colnames(markers)
  if (is.null(marker_ids)) {
    marker_ids <- paste0("marker", seq_len(m))
  }

  spec <- if (is.null(K)) {
    # Build the full n x n genomic relationship (avg-allele-frequency, /m;
    # Astle & Balding 2009) and take its *complete*
    # eigendecomposition. REML and the EMMAX rotation need the full
    # n-dimensional eigenbasis -- including the zero-eigenvalue directions
    # that exist when there are fewer markers than individuals (a
    # rank-deficient G) -- which the marker-SVD shortcut, truncated to rank
    # min(n, m), does not provide. The covariance V = sigma_g^2 G +
    # sigma_e^2 I is full rank regardless, so .fb_reml_vc() needs every
    # eigenvector.
    Zc <- scale(markers, center = TRUE, scale = FALSE)
    Zc[!is.finite(Zc)] <- 0
    G <- tcrossprod(Zc) / ncol(markers)
    .fb_spectral(G, tol = tol, name = "markers (genomic relationship, /m)")
  } else {
    K <- as.matrix(K)
    if (nrow(K) != n || ncol(K) != n) {
      stop(
        "`K` must be an ", n, " x ", n, " relationship matrix; got ",
        nrow(K), " x ", ncol(K), ".",
        call. = FALSE
      )
    }
    .fb_spectral(K, tol = tol, name = "K")
  }

  reml <- .fb_reml_vc(spec, y, X)

  scan <- .fb_emmax_scan(spec, reml, markers)

  chisq <- scan$statistic
  lambda_gc <- stats::median(chisq, na.rm = TRUE) / stats::qchisq(0.5, 1L)
  p_value <- scan$p_value

  results <- data.frame(
    marker = marker_ids,
    effect = scan$effect,
    se = scan$se,
    statistic = chisq,
    p_value = p_value,
    p_bonferroni = pmin(1, p_value * m),
    q_value = stats::p.adjust(p_value, method = "BH"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  if (!is.null(marker_map)) {
    results <- .fb_gwas_join_map(results, marker_map)
  }

  structure(
    list(
      results = results,
      lambda_gc = lambda_gc,
      var_g = reml$var_g,
      var_e = reml$var_e,
      h2 = reml$h2,
      n = n,
      n_markers = m,
      n_fixed = ncol(X),
      call = match.call()
    ),
    class = c("fb_gwas", "list")
  )
}

# Vectorised EMMAX score test. Rotate the markers into the eigenbasis,
# scale every quantity by 1/sqrt(lambda + delta) to homoscedastic form,
# residualise the markers and response against the covariates, then test
# each marker by ordinary least squares with the residual variance fixed
# at the null REML estimate sigma_g^2 (P3D). All markers are processed by
# matrix algebra, no per-marker loop over model fits.
.fb_emmax_scan <- function(spec, reml, markers) {
  sq <- sqrt(1 / reml$weights) # = sqrt(lambda + delta)
  y_til <- reml$y_star / sq
  X_til <- reml$X_star / sq

  # Projection onto the orthogonal complement of the (scaled) covariates.
  xtx_inv <- solve(crossprod(X_til))
  hat <- X_til %*% xtx_inv %*% t(X_til)
  resid_op <- function(M) M - hat %*% M

  y_perp <- resid_op(y_til)

  # Rotate + scale + residualise every marker at once.
  M_star <- crossprod(spec$vectors, markers) # rank x m
  M_til <- M_star / sq
  M_perp <- resid_op(M_til)

  ss <- colSums(M_perp^2)
  xy <- as.vector(crossprod(M_perp, y_perp))
  # Markers with no residual variance (e.g. collinear with covariates)
  # cannot be tested; flag as NA rather than divide by zero.
  ok <- ss > .Machine$double.eps^0.5
  effect <- rep(NA_real_, length(ss))
  se <- rep(NA_real_, length(ss))
  effect[ok] <- xy[ok] / ss[ok]
  se[ok] <- sqrt(reml$var_g / ss[ok])
  stat <- (effect / se)^2
  p_value <- stats::pchisq(stat, df = 1L, lower.tail = FALSE)

  list(effect = effect, se = se, statistic = stat, p_value = p_value)
}

# --- display ------------------------------------------------------ #

#' @exportS3Method print fb_gwas
print.fb_gwas <- function(x, ...) {
  cat("<fb_gwas>  EMMAX / P3D genome scan\n")
  cat(
    "  ", x$n, " individuals x ", x$n_markers, " markers; ",
    x$n_fixed, " fixed-effect column(s)\n",
    sep = ""
  )
  cat(
    "  null REML: h^2 = ", format(round(x$h2, 3L)),
    "  (var_g = ", format(round(x$var_g, 4L)),
    ", var_e = ", format(round(x$var_e, 4L)), ")\n",
    sep = ""
  )
  cat("  genomic-control lambda_GC: ", format(round(x$lambda_gc, 3L)), "\n",
    sep = ""
  )
  n_sig <- sum(x$results$p_bonferroni < 0.05, na.rm = TRUE)
  cat("  Bonferroni-significant markers (5%): ", n_sig, "\n", sep = "")
  top <- x$results[order(x$results$p_value), , drop = FALSE]
  top <- utils::head(top, 5L)
  cat("  top markers:\n")
  for (i in seq_len(nrow(top))) {
    cat(
      "    ", top$marker[i], "  p = ", format(signif(top$p_value[i], 3L)),
      "  effect = ", format(round(top$effect[i], 3L)), "\n",
      sep = ""
    )
  }
  invisible(x)
}

#' @exportS3Method plot fb_gwas
plot.fb_gwas <- function(x, type = c("manhattan", "qq"), ...) {
  type <- match.arg(type)
  res <- x$results
  if (identical(type, "qq")) {
    obs <- -log10(sort(res$p_value))
    expd <- -log10(stats::ppoints(length(obs)))
    graphics::plot(
      expd, obs,
      xlab = "expected -log10(p)", ylab = "observed -log10(p)",
      main = paste0("QQ (lambda_GC = ", round(x$lambda_gc, 3L), ")"),
      pch = 16L, col = "#3366aa", ...
    )
    graphics::abline(0, 1, col = "grey50", lty = 2L)
    return(invisible(x))
  }
  pos <- if (!is.null(res$pos)) res$pos else seq_len(nrow(res))
  chr <- if (!is.null(res$chr)) as.integer(factor(res$chr)) else 1L
  graphics::plot(
    pos, -log10(res$p_value),
    xlab = "position", ylab = "-log10(p)", main = "Manhattan",
    pch = 16L, col = c("#3366aa", "#aa6633")[1L + (chr %% 2L)], ...
  )
  graphics::abline(
    h = -log10(0.05 / nrow(res)), col = "red", lty = 2L
  )
  invisible(x)
}

# --- internal helpers --------------------------------------------- #

# Join a marker map (marker / chr / pos) onto the results, preserving the
# scan order.
.fb_gwas_join_map <- function(results, marker_map) {
  marker_map <- as.data.frame(marker_map)
  if (!"marker" %in% names(marker_map)) {
    stop("`marker_map` must have a `marker` column.", call. = FALSE)
  }
  idx <- match(results$marker, marker_map$marker)
  for (col in setdiff(names(marker_map), "marker")) {
    results[[col]] <- marker_map[[col]][idx]
  }
  results
}
