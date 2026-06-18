# Identified-quantity reporting for factor-analytic structured covariance.
#
# A factor-analytic term fa(outer, k) decomposes the outer-factor
# covariance as G = Lambda Lambda' + diag(psi), with Lambda an
# (n_outer x k) loadings matrix. The loadings themselves are only
# identified up to an orthogonal rotation and a per-column sign flip, so
# the marginal posterior of an individual Lambda entry is multimodal and
# its Rhat is meaningless -- a large Rhat on a raw loading is expected and
# is not evidence of non-convergence. The implied covariance G (and the
# correlation derived from it) is invariant to that rotation and sign, so
# G is the identified quantity: its posterior is interpretable and its
# Rhat is a genuine convergence diagnostic.
#
# fb_structured_cov() reconstructs G per posterior draw from the Lambda
# and psi draws, and returns the posterior summary plus an entrywise Rhat.
# Reporting G sidesteps loadings sign-alignment entirely.

#' Identified covariance for factor-analytic structured-covariance terms
#'
#' For each `fa(outer, k)` term in a greta fit, reconstruct the implied
#' outer-factor covariance \eqn{G = \Lambda\Lambda^\top + \mathrm{diag}(\psi)}
#' from the posterior draws and summarise it. Unlike the raw loadings
#' \eqn{\Lambda} -- which are identified only up to rotation and sign, so
#' their per-entry Rhat is meaningless -- the covariance \eqn{G} and the
#' correlation derived from it are rotation- and sign-invariant. They are
#' therefore the identified quantities whose posterior is interpretable
#' and whose Rhat is a genuine convergence diagnostic. Consult this in
#' preference to the raw-loading Rhat when judging whether a
#' factor-analytic fit has converged.
#'
#' @param fit A `flexybayes` fit produced on the greta backend with at
#'   least one `fa()` structured-covariance term.
#'
#' @return A named list with one entry per factor-analytic term (named by
#'   the term's outer factor). Each entry is a list with: `levels` (the
#'   outer-factor levels labelling the rows/columns), `cov_mean`,
#'   `cov_lower`, `cov_upper` (posterior mean and 95% interval of
#'   \eqn{G}), `cor_mean` (posterior-mean correlation), `rhat` (entrywise
#'   Rhat of \eqn{G}, `NA` if fewer than two chains), `max_rhat`, and
#'   `k` (the number of factors). Returns an empty list (with a message)
#'   when the fit carries no factor-analytic term. Non-factor-analytic
#'   structured terms (`us`, `ar1`) are reported as not-yet-reconstructed.
#'
#' @examples
#' \dontrun{
#' fit <- flexybayes(
#'   y ~ 1, random = ~ fa(env, 2):id(geno), data = met_data,
#'   family = "gaussian", backend = "greta"
#' )
#' sc <- fb_structured_cov(fit)
#' sc$env$cov_mean   # identified genetic covariance across environments
#' sc$env$max_rhat   # convergence of the identified quantity
#' }
#'
#' @export
fb_structured_cov <- function(fit) {
  if (!inherits(fit, "flexybayes")) {
    stop("`fit` must be a flexybayes object.", call. = FALSE)
  }
  rt <- fit$extras$parse_info$random %||% list()
  fa_terms <- Filter(function(t) identical(t$type %||% "", "fa_gxe"), rt)
  other_struct <- Filter(
    function(t) (t$type %||% "") %in% c("us_gxe", "ar1_spatial"),
    rt
  )

  if (!length(fa_terms)) {
    if (length(other_struct)) {
      message(
        "flexyBayes: fb_structured_cov() reconstructs the identified ",
        "covariance for factor-analytic fa() terms; reconstruction for ",
        "us()/ar1() terms is not yet implemented."
      )
    } else {
      message(
        "flexyBayes: this fit carries no factor-analytic structured-",
        "covariance term; nothing to report."
      )
    }
    return(invisible(list()))
  }

  draws <- tryCatch(fit$greta$draws, error = function(e) NULL)
  if (is.null(draws)) {
    stop(
      "fb_structured_cov() requires the greta posterior draws; this fit ",
      "does not carry them (only the greta backend is supported).",
      call. = FALSE
    )
  }

  out <- list()
  for (term in fa_terms) {
    res <- .fb_fa_identified_cov(draws, term)
    out[[term$outer]] <- res
  }
  if (length(other_struct)) {
    message(
      "flexyBayes: fb_structured_cov() reported the factor-analytic ",
      "term(s); reconstruction for us()/ar1() terms is not yet ",
      "implemented and was skipped."
    )
  }
  out
}

# Reconstruct G = Lambda Lambda' + diag(psi) per draw for one fa term and
# summarise it. `draws` is a coda mcmc.list (one mcmc per chain); columns
# are named Lambda_<tag>[i,j] and psi_<tag>[i,1].
.fb_fa_identified_cov <- function(draws, term) {
  tag <- paste0(term$inner, "_", term$outer, "_fa", term$k)
  no <- term$n_outer
  k <- term$k

  lambda_cols <- as.vector(vapply(
    seq_len(k),
    function(j) sprintf("Lambda_%s[%d,%d]", tag, seq_len(no), j),
    character(no)
  ))
  psi_cols <- sprintf("psi_%s[%d,1]", tag, seq_len(no))

  chain_mats <- lapply(draws, as.matrix)
  avail <- colnames(chain_mats[[1]])
  if (!all(c(lambda_cols, psi_cols) %in% avail)) {
    stop(
      "fb_structured_cov(): expected Lambda/psi draw columns for term '",
      tag, "' are absent from the posterior; the fit may predate the ",
      "current fa() codegen.",
      call. = FALSE
    )
  }

  # Per-chain entrywise G draws for Rhat; pooled draws for the summary.
  # Upper-triangle (incl. diagonal) entries indexed once.
  idx <- which(upper.tri(matrix(0, no, no), diag = TRUE), arr.ind = TRUE)
  entry_label <- sprintf("G[%d,%d]", idx[, 1], idx[, 2])

  g_one_draw <- function(row) {
    Lambda <- matrix(row[lambda_cols], nrow = no, ncol = k)
    psi <- row[psi_cols]
    Lambda %*% t(Lambda) + diag(psi, nrow = no)
  }

  per_chain_entries <- lapply(chain_mats, function(m) {
    t(apply(m, 1L, function(r) {
      G <- g_one_draw(r)
      G[idx]
    }))
  })

  pooled <- do.call(rbind, per_chain_entries)
  mean_vec <- colMeans(pooled)
  lo_vec <- apply(pooled, 2L, stats::quantile, probs = 0.025, names = FALSE)
  hi_vec <- apply(pooled, 2L, stats::quantile, probs = 0.975, names = FALSE)

  fill_sym <- function(vec) {
    M <- matrix(NA_real_, no, no)
    M[idx] <- vec
    M[lower.tri(M)] <- t(M)[lower.tri(M)]
    M
  }
  cov_mean <- fill_sym(mean_vec)
  cov_lower <- fill_sym(lo_vec)
  cov_upper <- fill_sym(hi_vec)

  # Correlation from the posterior-mean covariance.
  d <- sqrt(diag(cov_mean))
  cor_mean <- cov_mean / outer(d, d)

  rhat_vec <- rep(NA_real_, length(entry_label))
  if (length(per_chain_entries) >= 2L) {
    ml <- coda::as.mcmc.list(lapply(per_chain_entries, coda::as.mcmc))
    psrf <- tryCatch(
      coda::gelman.diag(ml, multivariate = FALSE, autoburnin = FALSE)$psrf[, 1],
      error = function(e) rep(NA_real_, length(entry_label))
    )
    rhat_vec <- psrf
  }
  rhat <- matrix(NA_real_, no, no)
  rhat[idx] <- rhat_vec
  rhat[lower.tri(rhat)] <- t(rhat)[lower.tri(rhat)]

  levels <- term$outer_levels %||% as.character(seq_len(no))
  dimnames(cov_mean) <- list(levels, levels)
  dimnames(cov_lower) <- list(levels, levels)
  dimnames(cov_upper) <- list(levels, levels)
  dimnames(cor_mean) <- list(levels, levels)
  dimnames(rhat) <- list(levels, levels)

  max_rhat <- if (all(is.na(rhat_vec))) {
    NA_real_
  } else {
    max(rhat_vec, na.rm = TRUE)
  }
  list(
    levels = levels,
    k = k,
    cov_mean = cov_mean,
    cov_lower = cov_lower,
    cov_upper = cov_upper,
    cor_mean = cor_mean,
    rhat = rhat,
    max_rhat = max_rhat
  )
}
