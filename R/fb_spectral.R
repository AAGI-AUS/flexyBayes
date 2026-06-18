# Spectral primitive --- the shared efficiency machinery (G-I4).
#
# `.fb_spectral()` eigendecomposes a relationship / kinship matrix K
# once (K = U Lambda U', symmetric PSD) and packages the rotation that
# turns a structured genetic random effect into an independent one. It
# is the single piece of new numerical machinery the genomics /
# MET expansion adds; everything downstream composes it:
#
#   * GWAS null scan (R/fb_gwas.R): one O(n^3) decomposition, then the
#     rotated model y* = U'y, X* = U'X has a diagonal residual
#     covariance D = sigma_g^2 Lambda + sigma_e^2 I, so each marker is
#     an O(n) generalised-least-squares score test (EMMAX / P3D) rather
#     than a fresh O(n^3) mixed-model fit.
#   * GBLUP scaling: the genetic random effect u ~ N(0, sigma_g^2 K)
#     rotates to alpha = U'u ~ N(0, sigma_g^2 Lambda), n independent
#     coordinates a backend fits cheaply.
#   * REML variance-component estimation: the rotated profile
#     log-likelihood is a fast low-dimensional optimisation in
#     (sigma_g^2, sigma_e^2).
#
# The object IS the cache: it carries U and Lambda, so a caller
# computes it once per K and reuses it across every marker, fold, or
# variance-component evaluation. There is no hidden global state ---
# correctness over magic; consumers hold the returned object (e.g. on
# `fit$extras$spectral`).
#
# This file is engine-agnostic base R + `Matrix`: no greta / INLA /
# brms dependency, so it is unit-testable to exact numerical tolerance
# without a backend. The user-facing fit routes that consume it
# (vm() scaling, fb_gwas()) raise the typed structured-covariance
# refusals; this internal primitive validates its own arguments with
# plain base-R stop() (the humane caller-facing form, r_style E9).

# --- constructor -------------------------------------------------- #

# `.fb_spectral()` --- eigendecompose a symmetric PSD matrix.
#
# Arguments
#   K      A symmetric positive-semidefinite matrix (base-R matrix or a
#          `Matrix` object). The genomic relationship matrix G, the
#          pedigree relationship matrix A, or any kinship matrix. Sparse
#          `Matrix` inputs are densified for `eigen()` (eigendecom-
#          position is inherently dense); the O(n^2) materialisation is
#          accepted deliberately --- relationship matrices are dense.
#   rank   Optional integer: keep only the leading `rank` eigenpairs
#          (the rank-K truncation used for very large n, reusing the
#          existing approximation vocabulary). NULL (default) keeps the
#          full decomposition.
#   tol    Relative tolerance for the PSD contract. Eigenvalues in
#          [-tol * lambda_max, 0) are numerical noise and are clamped to
#          0 (recorded honestly on `negative_clamped`); an eigenvalue
#          below -tol * lambda_max is genuine indefiniteness and is
#          refused (K is not a valid covariance).
#   name   Symbol name for messages (e.g. the known-matrix label).
#
# Returns an `fb_spectral` object: a classed list with
#   $values         kept eigenvalues, descending, clamped non-negative.
#   $vectors        U (n x rank), orthonormal columns, K %*% U = U Lambda.
#   $rank           number of kept eigenpairs.
#   $rank_full      n (the matrix dimension).
#   $n              n.
#   $trace          sum of all non-negative eigenvalues (full spectrum).
#   $capture_trace  sum(kept) / sum(all) --- proportion of genetic
#                   variance retained (the breeder's "variance captured").
#   $capture_frobenius
#                   sum(kept^2) / sum(all^2) --- the Frobenius capture
#                   that matches the package's low_rank_smooth bias
#                   bound (1 - this is the relative squared truncation
#                   error).
#   $negative_clamped
#                   list(count, min) --- how many eigenvalues were
#                   clamped from small-negative to 0 and the most
#                   negative value seen (NA min when none). Surfaced so
#                   PSD repair is never silent.
#   $tol, $name     echoed inputs.
.fb_spectral <- function(K, rank = NULL, tol = 1e-8, name = "K") {
  K <- .fb_spectral_as_dense(K, name)
  d <- dim(K)
  if (is.null(d) || length(d) != 2L || d[[1L]] != d[[2L]]) {
    stop(
      "`", name, "` must be a square matrix for spectral decomposition; ",
      "got dim ", paste(d, collapse = " x "), ".",
      call. = FALSE
    )
  }
  n <- d[[1L]]
  if (n == 0L) {
    stop("`", name, "` is a 0 x 0 matrix; nothing to decompose.", call. = FALSE)
  }

  if (!.fb_spectral_is_symmetric(K, tol)) {
    stop(
      "`", name, "` must be symmetric for spectral decomposition. A ",
      "relationship / kinship matrix is symmetric by construction; ",
      "isSymmetric() returned FALSE beyond tolerance ", format(tol), ".",
      call. = FALSE
    )
  }

  if (!is.null(rank)) {
    if (
      length(rank) != 1L || is.na(rank) || rank < 1L || rank > n ||
        rank != as.integer(rank)
    ) {
      stop(
        "`rank` must be a single integer in 1:", n, " (the dimension of `",
        name, "`); got ", paste(rank, collapse = ", "), ".",
        call. = FALSE
      )
    }
    rank <- as.integer(rank)
  }

  e <- eigen(K, symmetric = TRUE)
  vals_raw <- e$values
  lambda_max <- max(vals_raw)
  if (!is.finite(lambda_max)) {
    stop(
      "spectral decomposition of `", name, "` produced non-finite ",
      "eigenvalues; check the matrix for NA / Inf entries.",
      call. = FALSE
    )
  }

  # PSD contract: clamp numerical-noise negatives, refuse genuine
  # indefiniteness. lambda_max <= 0 means K is not PSD at all.
  neg_floor <- -abs(tol) * max(lambda_max, .Machine$double.eps)
  worst_neg <- min(vals_raw)
  if (lambda_max <= 0 || worst_neg < neg_floor) {
    stop(
      "`", name, "` is not positive-semidefinite: its smallest eigenvalue ",
      format(worst_neg, digits = 4L), " is below the tolerance floor ",
      format(neg_floor, digits = 4L), " (lambda_max = ",
      format(lambda_max, digits = 4L), "). A covariance / relationship ",
      "matrix must be PSD; check the matrix construction.",
      call. = FALSE
    )
  }
  clamped <- vals_raw < 0
  n_clamped <- sum(clamped)
  vals <- vals_raw
  vals[clamped] <- 0

  trace_all <- sum(vals)
  frob_all <- sum(vals^2)
  keep <- if (is.null(rank)) n else rank
  idx <- seq_len(keep)
  vals_keep <- vals[idx]
  vec_keep <- e$vectors[, idx, drop = FALSE]

  capture_trace <- if (trace_all > 0) sum(vals_keep) / trace_all else 1
  capture_frob <- if (frob_all > 0) sum(vals_keep^2) / frob_all else 1

  structure(
    list(
      values = vals_keep,
      vectors = vec_keep,
      rank = keep,
      rank_full = n,
      n = n,
      trace = trace_all,
      capture_trace = capture_trace,
      capture_frobenius = capture_frob,
      negative_clamped = list(
        count = n_clamped,
        min = if (n_clamped > 0L) worst_neg else NA_real_
      ),
      tol = tol,
      name = name
    ),
    class = c("fb_spectral", "list")
  )
}

# `.fb_spectral_from_markers()` --- build the spectral object directly
# from a marker matrix via the SVD of the centred genotypes, avoiding
# the explicit n x n product K = Z Z' / m when that is less stable or
# unnecessary.
#
# For centred markers Z_c (n individuals x m markers) the VanRaden
# (2008) genomic relationship is G = Z_c Z_c' / m. The SVD
# Z_c = P S Q' gives G = P (S^2 / m) P', so the eigenvectors of G are
# the left singular vectors P and the eigenvalues are diag(S^2) / m ---
# computed without ever forming G. Downstream code receives the same
# `fb_spectral` object shape as `.fb_spectral(G)`.
#
# Arguments
#   Z       n x m marker matrix (allele dosages 0/1/2 or any numeric
#           coding).
#   m       Normalising constant for the relationship scale. Default
#           ncol(Z); the VanRaden 2 sum p (1 - p) normaliser is the
#           caller's declared choice (kinship scaling shifts the
#           heritability posterior by a constant --- see the MET /
#           genomics vignette, "Kinship matrix scaling"). Stated here,
#           not assumed.
#   center  Centre each marker column to mean 0 (the additive-relation-
#           ship convention). Default TRUE.
#   scale   Scale each marker column to unit variance (Astle & Balding
#           2009 style). Default FALSE.
#   rank    Optional leading-rank truncation, as `.fb_spectral()`.
#   tol     PSD / symmetry tolerance, as `.fb_spectral()`.
#   name    Symbol name for messages.
.fb_spectral_from_markers <- function(
  Z,
  m = NULL,
  center = TRUE,
  scale = FALSE,
  rank = NULL,
  tol = 1e-8,
  name = "Z"
) {
  Z <- .fb_spectral_as_dense(Z, name)
  d <- dim(Z)
  if (is.null(d) || length(d) != 2L) {
    stop(
      "`", name, "` must be an n x m marker matrix; got ",
      paste(class(Z), collapse = "/"), ".",
      call. = FALSE
    )
  }
  n <- d[[1L]]
  n_markers <- d[[2L]]
  if (is.null(m)) {
    m <- n_markers
  }
  if (length(m) != 1L || is.na(m) || m <= 0) {
    stop(
      "`m` (the relationship normaliser) must be a single positive ",
      "number; got ", paste(m, collapse = ", "), ".",
      call. = FALSE
    )
  }
  Zc <- scale(Z, center = center, scale = scale)
  # scale() carries centring / scaling attributes; drop them so the
  # downstream matrix algebra sees a plain numeric matrix.
  attr(Zc, "scaled:center") <- NULL
  attr(Zc, "scaled:scale") <- NULL
  bad <- !is.finite(Zc)
  if (any(bad)) {
    stop(
      "centring / scaling `", name, "` produced non-finite values (",
      sum(bad), " entries) --- a marker column with zero variance cannot ",
      "be unit-scaled. Drop monomorphic markers before forming the ",
      "relationship matrix.",
      call. = FALSE
    )
  }

  if (!is.null(rank)) {
    max_rank <- min(n, n_markers)
    if (
      length(rank) != 1L || is.na(rank) || rank < 1L || rank > max_rank ||
        rank != as.integer(rank)
    ) {
      stop(
        "`rank` must be a single integer in 1:", max_rank,
        " (min(n, m) for the marker SVD); got ", paste(rank, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    rank <- as.integer(rank)
  }

  sv <- svd(Zc)
  eig <- sv$d^2 / m
  keep_full <- length(eig)
  keep <- if (is.null(rank)) keep_full else rank
  idx <- seq_len(keep)

  trace_all <- sum(eig)
  frob_all <- sum(eig^2)
  vals_keep <- eig[idx]
  vec_keep <- sv$u[, idx, drop = FALSE]

  structure(
    list(
      values = vals_keep,
      vectors = vec_keep,
      rank = keep,
      rank_full = keep_full,
      n = n,
      trace = trace_all,
      capture_trace = if (trace_all > 0) sum(vals_keep) / trace_all else 1,
      capture_frobenius = if (frob_all > 0) sum(vals_keep^2) / frob_all else 1,
      negative_clamped = list(count = 0L, min = NA_real_),
      tol = tol,
      name = name
    ),
    class = c("fb_spectral", "list")
  )
}

# --- rotation operators ------------------------------------------- #

# `.fb_spectral_rotate()` --- left-multiply by U', mapping data into the
# eigenbasis where the genetic random effect is independent. For a
# response vector y this returns y* = U'y; for a design matrix X it
# returns X* = U'X. The result has `rank` rows.
.fb_spectral_rotate <- function(spec, M) {
  .fb_spectral_check(spec)
  M <- .fb_spectral_conform(spec, M, "rotate")
  crossprod(spec$vectors, M)
}

# `.fb_spectral_backrotate()` --- left-multiply by U, mapping a quantity
# in the eigenbasis back to the original observation basis (U A). The
# input must have `rank` rows.
.fb_spectral_backrotate <- function(spec, A) {
  .fb_spectral_check(spec)
  A <- as.matrix(A)
  if (nrow(A) != spec$rank) {
    stop(
      "back-rotation expects ", spec$rank, " rows (the spectral rank); ",
      "got ", nrow(A), ".",
      call. = FALSE
    )
  }
  spec$vectors %*% A
}

# `.fb_spectral_sqrt()` --- the matrix square root B = U Lambda^{1/2}
# with B B' = K (rank-truncated when rank < n). The decorrelating
# factor: if alpha ~ N(0, sigma^2 I_rank) then u = B alpha ~
# N(0, sigma^2 K). Used to simulate genetic effects with known truth
# (parameter-recovery cells) and as the spectral alternative to the
# Cholesky K = L L' rotation.
.fb_spectral_sqrt <- function(spec) {
  .fb_spectral_check(spec)
  sweep(spec$vectors, 2L, sqrt(spec$values), `*`)
}

# `.fb_spectral_dvar()` --- the diagonal residual variances of the
# rotated mixed model: d_i = sigma_g^2 lambda_i + sigma_e^2, the
# per-coordinate variances of y* = U'y under
# Cov(y) = sigma_g^2 K + sigma_e^2 I. Length `rank`. The GLS / EMMAX
# weights are 1 / d_i. For a truncated spectrum the n - rank directions
# in the orthogonal complement carry residual variance sigma_e^2 only;
# that complement is the consumer's responsibility (documented), this
# helper returns the captured-subspace variances.
.fb_spectral_dvar <- function(spec, var_g, var_e) {
  .fb_spectral_check(spec)
  .fb_spectral_check_vc(var_g, var_e)
  var_g * spec$values + var_e
}

# `.fb_spectral_logdet()` --- log|sigma_g^2 K + sigma_e^2 I| via the
# spectrum: sum_i log(sigma_g^2 lambda_i + sigma_e^2). Requires the
# full-rank decomposition (rank == n); the truncated determinant would
# need the n - rank residual-only directions added back, so it is
# refused here rather than returned wrong. This is the determinant term
# of the mixed-model REML / GLS profile log-likelihood the GWAS null
# fit and variance-component estimation evaluate.
.fb_spectral_logdet <- function(spec, var_g, var_e) {
  .fb_spectral_check(spec)
  .fb_spectral_check_vc(var_g, var_e)
  if (spec$rank != spec$n) {
    stop(
      "`.fb_spectral_logdet()` requires the full-rank decomposition ",
      "(rank == n); this object is truncated to rank ", spec$rank, " of ",
      spec$n, ". The truncated determinant needs the residual-only ",
      "complement and is not computed here.",
      call. = FALSE
    )
  }
  sum(log(var_g * spec$values + var_e))
}

# --- predicates + display ----------------------------------------- #

# `is_fb_spectral()` --- TRUE for an `.fb_spectral()` object.
is_fb_spectral <- function(x) inherits(x, "fb_spectral")

#' @exportS3Method print fb_spectral
print.fb_spectral <- function(x, ...) {
  cat("<fb_spectral> ", x$name, ": ", x$n, " x ", x$n, "\n", sep = "")
  cat(
    "  rank: ", x$rank, "/", x$rank_full,
    if (x$rank < x$rank_full) " (truncated)" else " (full)", "\n",
    sep = ""
  )
  cat(
    "  variance captured: ", format(round(x$capture_trace, 4L)),
    " (trace), ", format(round(x$capture_frobenius, 4L)), " (Frobenius)\n",
    sep = ""
  )
  nc <- x$negative_clamped
  if (!is.null(nc) && nc$count > 0L) {
    cat(
      "  PSD repair: ", nc$count, " small-negative eigenvalue",
      if (nc$count != 1L) "s" else "", " clamped to 0 (min ",
      format(nc$min, digits = 3L), ")\n",
      sep = ""
    )
  }
  invisible(x)
}

# --- internal helpers --------------------------------------------- #

# Coerce a `Matrix` object to a dense base-R matrix for the dense
# LAPACK routines (eigen / svd), leaving base matrices untouched.
.fb_spectral_as_dense <- function(M, name) {
  if (is.null(M)) {
    stop("`", name, "` is NULL; expected a matrix.", call. = FALSE)
  }
  if (inherits(M, "Matrix")) {
    return(as.matrix(M))
  }
  if (!is.matrix(M)) {
    M <- tryCatch(
      as.matrix(M),
      error = function(e) {
        stop(
          "`", name, "` could not be coerced to a matrix: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }
  if (!is.numeric(M)) {
    stop("`", name, "` must be numeric; got ", typeof(M), ".", call. = FALSE)
  }
  M
}

# Symmetry probe tolerant of floating-point asymmetry, on the values
# only (dimnames are irrelevant to the numerical contract).
.fb_spectral_is_symmetric <- function(K, tol) {
  K_un <- K
  dimnames(K_un) <- NULL
  isSymmetric(K_un, tol = max(100 * .Machine$double.eps, abs(tol)))
}

# Guard: the object is a well-formed spectral decomposition.
.fb_spectral_check <- function(spec) {
  if (!is_fb_spectral(spec)) {
    stop(
      "expected an `fb_spectral` object from `.fb_spectral()`; got ",
      paste(class(spec), collapse = "/"), ".",
      call. = FALSE
    )
  }
  invisible(spec)
}

# Guard: variance components are single non-negative finite numbers.
.fb_spectral_check_vc <- function(var_g, var_e) {
  ok <- function(v) {
    length(v) == 1L && is.numeric(v) && is.finite(v) && v >= 0
  }
  if (!ok(var_g) || !ok(var_e)) {
    stop(
      "`var_g` and `var_e` must each be a single non-negative finite ",
      "number (the genetic and residual variances).",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Conform an input to an n-row matrix for rotation, accepting a vector
# (treated as a single column) or a matrix with n rows.
.fb_spectral_conform <- function(spec, M, op) {
  if (is.null(dim(M))) {
    M <- matrix(M, ncol = 1L)
  } else {
    M <- as.matrix(M)
  }
  if (nrow(M) != spec$n) {
    stop(
      op, " expects ", spec$n, " rows (the matrix dimension of `",
      spec$name, "`); got ", nrow(M), ".",
      call. = FALSE
    )
  }
  M
}
