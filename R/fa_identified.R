# Factor-analytic genotype-by-environment covariance, fa(f, k):g: the choice
# of anchor levels, the identified covariance of a fit, and the agreement of
# several fits on it.
#
# The INLA emit fixes the loading of factor j to one at its anchor level and
# enters factor j only at the anchor and the levels placed after it. That is
# a reparameterisation of G = Lambda Lambda' + Psi only when the anchors'
# loadings are not zero, so the anchors are chosen from the data rather than
# taken from factor-level order.


# .fb_fa_levels() -- the levels of f in the order the emit enters factors:
# the k anchors first, then the other levels in factor order.
.fb_fa_levels <- function(term, levs) {
  if (is.null(term$fa_order)) {
    return(levs)
  }
  levs[term$fa_order]
}


# .fb_fa_choose_anchors() -- set `fa_order` on every fa(f, k):g term.
#
# `response` is the IR's response part, read by .fb_working_response().
.fb_fa_choose_anchors <- function(random_terms, data, response) {
  y <- .fb_working_response(response, data)
  lapply(random_terms, function(term) {
    if (!identical(term$type, "fa_gxe")) {
      return(term)
    }
    term$fa_order <- .fb_fa_anchor_order(term, data, y)
    term
  })
}


# .fb_fa_anchor_order() -- a permutation of the levels of f that puts the
# anchors first.
#
# The genotype-by-environment table of cell means gives a moment estimate S
# of G. Its common part is estimated by replacing each diagonal entry with
# the variance the other environments predict, S_ee - 1 / (S^-1)_ee (the
# squared multiple correlation, Guttman 1956, doi:10.1007/BF02289137), so a
# level whose variance is all specific variance or error is not mistaken for
# one that carries the factors. A pivoted Cholesky step then takes the level
# with the largest common variance not yet explained by the earlier anchors.
# NULL, which keeps factor order, when the data give no estimate.
.fb_fa_anchor_order <- function(term, data, y) {
  if (!all(c(term$outer, term$inner) %in% names(data)) || !is.numeric(y)) {
    return(NULL)
  }
  ok <- is.finite(y)
  outer_f <- factor(data[[term$outer]])[ok]
  inner_f <- factor(data[[term$inner]])[ok]
  mat_m <- tapply(y[ok], list(inner_f, outer_f), mean)
  mat_s <- suppressWarnings(stats::cov(mat_m, use = "pairwise.complete.obs"))
  if (anyNA(diag(mat_s))) {
    return(NULL)
  }
  # A pair of environments with fewer than two genotypes in common.
  mat_s[is.na(mat_s)] <- 0

  mat_c <- mat_s
  if (rcond(mat_s) > sqrt(.Machine$double.eps)) {
    diag(mat_c) <- diag(mat_s) - 1 / diag(solve(mat_s))
  }

  s1 <- ncol(mat_c)
  vec_ord <- integer(0)
  for (j1 in seq_len(min(term$k, s1))) {
    vec_d <- diag(mat_c)
    vec_d[vec_ord] <- -Inf
    p1 <- unname(which.max(vec_d))
    vec_ord <- c(vec_ord, p1)
    if (mat_c[p1, p1] > 0) {
      vec_l <- mat_c[, p1] / sqrt(mat_c[p1, p1])
      mat_c <- mat_c - tcrossprod(vec_l)
    }
  }
  c(vec_ord, setdiff(seq_len(s1), vec_ord))
}


# .fb_fa_g() -- G = Lambda diag(sd_factor^2) Lambda' + diag(psi^2).
.fb_fa_g <- function(mat_l, vec_sf, vec_psi) {
  mat_l %*% diag(vec_sf^2, nrow = length(vec_sf)) %*% t(mat_l) +
    diag(vec_psi^2, nrow = length(vec_psi))
}


# .fb_fa_identified() -- the genetic covariance G and the genetic
# correlations of a fa(f, k):g fit on INLA.
#
# Both are computed on each of `n` joint draws of the hyperparameters
# (fb_variance_draws()) and summarised by their posterior medians, entry by
# entry. G and the correlations do not change when the loadings are rotated
# or change sign; the loadings themselves do. Rows and columns follow the
# level order of f.
.fb_fa_identified <- function(fit, n = 2000L, seed = 1L) {
  is_fa <- vapply(
    fit$fb$random_terms,
    function(t) identical(t$type, "fa_gxe"),
    logical(1L)
  )
  if (!any(is_fa)) {
    stop("`fit` carries no fa(f, k):g term.", call. = FALSE)
  }
  term <- fit$fb$random_terms[[which(is_fa)[[1L]]]]
  levs <- term$outer_levels
  levs_fa <- .fb_fa_levels(term, levs)
  k <- term$k
  s1 <- length(levs)
  dr <- fb_variance_draws(fit, n, seed = seed)

  nm_psi <- paste0("sd_", term$inner, "__", term$outer, levs)
  nm_sf <- paste0("sd_fa", seq_len(k))
  # The loading of factor j at each level: one at its anchor, a draw at the
  # levels after it, zero before it.
  lst_l <- lapply(seq_len(k), function(j1) {
    mat <- matrix(0, nrow(dr), s1, dimnames = list(NULL, levs))
    mat[, levs_fa[[j1]]] <- 1
    for (lev in levs_fa[-seq_len(j1)]) {
      mat[, lev] <- dr[[paste0("fa", j1, "_", term$outer, lev)]]
    }
    mat
  })

  arr_g <- array(NA_real_, c(nrow(dr), s1, s1))
  arr_r <- arr_g
  for (i1 in seq_len(nrow(dr))) {
    mat_l <- vapply(lst_l, function(m) m[i1, ], numeric(s1))
    mat_l <- matrix(mat_l, nrow = s1)
    mat_g <- .fb_fa_g(
      mat_l,
      unlist(dr[i1, nm_sf]),
      unlist(dr[i1, nm_psi])
    )
    arr_g[i1, , ] <- mat_g
    arr_r[i1, , ] <- stats::cov2cor(mat_g)
  }

  mat_g <- apply(arr_g, c(2L, 3L), stats::median)
  mat_r <- apply(arr_r, c(2L, 3L), stats::median)
  dimnames(mat_g) <- list(levs, levs)
  dimnames(mat_r) <- list(levs, levs)
  list(G = mat_g, cor = mat_r, anchors = levs_fa[seq_len(k)], k = k,
       n_draws = nrow(dr))
}


# .fb_fa_agreement() -- how far several fits of one fa(f, k):g model are
# apart on the identified covariance.
#
# `fits` is a list of at least two results of .fb_fa_identified() (or any
# lists carrying `G` and `cor` over the same levels). Returned:
#   g_rel    per level, (max - min) / min of the genetic variance G_ee
#   cor_abs  per pair of levels, max - min of the genetic correlation
# and the largest of each. Covariances off the diagonal are compared through
# the correlations, because a covariance near zero has no stable relative
# difference.
.fb_fa_agreement <- function(fits) {
  if (length(fits) < 2L) {
    stop("Agreement needs at least two fits.", call. = FALSE)
  }
  levs <- rownames(fits[[1L]]$G)
  same <- vapply(
    fits,
    function(x) identical(rownames(x$G), levs) && identical(rownames(x$cor), levs),
    logical(1L)
  )
  if (!all(same)) {
    stop("The fits do not cover the same levels in the same order.",
         call. = FALSE)
  }
  mat_g <- vapply(fits, function(x) diag(x$G), numeric(length(levs)))
  mat_g <- matrix(mat_g, nrow = length(levs), dimnames = list(levs, NULL))
  g_min <- apply(mat_g, 1L, min)
  g_rel <- (apply(mat_g, 1L, max) - g_min) / g_min

  idx <- which(upper.tri(fits[[1L]]$cor), arr.ind = TRUE)
  nm_pair <- paste0(levs[idx[, 1L]], ":", levs[idx[, 2L]])
  mat_r <- vapply(fits, function(x) x$cor[idx], numeric(nrow(idx)))
  mat_r <- matrix(mat_r, nrow = nrow(idx), dimnames = list(nm_pair, NULL))
  cor_abs <- apply(mat_r, 1L, max) - apply(mat_r, 1L, min)

  list(
    g_rel_max = max(g_rel),
    cor_abs_max = if (length(cor_abs)) max(cor_abs) else 0,
    g_rel = g_rel,
    cor_abs = cor_abs,
    n_fits = length(fits)
  )
}


# Thresholds of the k-fit agreement gate, from U-FA-1: the largest relative
# spread of a genetic variance and the largest absolute spread of a genetic
# correlation across the k fits.
.FB_FA_GATE_G_REL <- 0.10
.FB_FA_GATE_COR_ABS <- 0.05


# .fb_fa_term() -- the first fa(f, k):g term of an IR, or NULL.
.fb_fa_term <- function(fb) {
  for (t in fb$random_terms %||% list()) {
    if (identical(t$type, "fa_gxe")) {
      return(t)
    }
  }
  NULL
}


# .fb_fa_mode_identified() -- G and the genetic correlations at an INLA
# fit's joint hyperparameter mode, `inla_fit$mode$theta`.
#
# Read on the coordinates INLA ran on. A standardised response scales every
# entry of G by the same factor, which leaves the relative spreads and the
# correlations the gate compares unchanged.
.fb_fa_mode_identified <- function(inla_fit, term) {
  th <- inla_fit$mode$theta
  levs <- term$outer_levels
  levs_fa <- .fb_fa_levels(term, levs)
  k <- term$k
  sd_of <- function(idx) exp(-th[[paste0("Log precision for ", idx)]] / 2)
  mat_l <- matrix(0, length(levs), k, dimnames = list(levs, NULL))
  vec_sf <- numeric(k)
  for (j1 in seq_len(k)) {
    vec_sf[[j1]] <- sd_of(.fb_inla_fa_index(term$outer, term$inner, j1,
                                            levs_fa[[j1]]))
    mat_l[levs_fa[[j1]], j1] <- 1
    for (lev in levs_fa[-seq_len(j1)]) {
      idx <- .fb_inla_fa_index(term$outer, term$inner, j1, lev)
      mat_l[lev, j1] <- th[[paste0("Beta_intern for ", idx)]]
    }
  }
  vec_psi <- vapply(
    levs,
    function(lev) sd_of(.fb_inla_fa_specific(term$outer, term$inner, lev)),
    numeric(1L)
  )
  mat_g <- .fb_fa_g(mat_l, vec_sf, vec_psi)
  dimnames(mat_g) <- list(levs, levs)
  list(G = mat_g, cor = stats::cov2cor(mat_g))
}


# .fb_fa_gate() -- the k-fit agreement record of a fa(f, k):g fit.
#
# `fits` holds the k INLA fits. `agree` is TRUE when both spreads are within
# the thresholds.
.fb_fa_gate <- function(fits, term) {
  agr <- .fb_fa_agreement(lapply(fits, .fb_fa_mode_identified, term = term))
  c(
    agr,
    list(
      g_tol = .FB_FA_GATE_G_REL,
      cor_tol = .FB_FA_GATE_COR_ABS,
      agree = agr$g_rel_max <= .FB_FA_GATE_G_REL &&
        agr$cor_abs_max <= .FB_FA_GATE_COR_ABS
    )
  )
}
