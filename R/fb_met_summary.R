# Breeder-facing MET summary (G2). A factor-analytic G x E fit
# fa(env, k):gen answers the two questions a plant breeder asks of every
# multi-environment trial:
#
#   * overall performance (OP): which genotypes are best on average across
#     environments? -- the across-environment mean of the realised
#     genotype-by-environment effects.
#   * stability: which genotypes hold their performance across
#     environments, and which trade mean for sensitivity? -- the
#     across-environment spread of those effects.
#
# plus the environment genetic-correlation matrix (which environments rank
# genotypes alike, and which reverse them -- the crossover structure a
# reverse-U / G x E paradox investigation reads directly), and the
# genotype-by-environment BLUPs themselves.
#
# These are computed from the *realised* effects g_mat = F Lambda' + delta
# (monitored by the fa codegen), which are identified -- rotation- and
# sign-invariant -- unlike the raw loadings. The environment covariance
# reuses the same identified reconstruction fb_structured_cov() reports.
# Greta-only, mirroring fb_structured_cov() (the factor-analytic term fits
# on the greta backend).

# --- the summary -------------------------------------------------- #

#' Breeder summary of a factor-analytic multi-environment-trial fit
#'
#' For a `fa(env, k):gen` factor-analytic G x E fit, summarise the
#' quantities a plant breeder acts on: each genotype's overall performance
#' (the across-environment mean of its realised effects) and stability (the
#' across-environment spread), the genotype-by-environment BLUPs, and the
#' environment genetic-correlation matrix (the crossover structure). The
#' realised effects are identified -- invariant to the rotation and sign
#' ambiguity of the raw loadings -- so their posterior summaries are
#' interpretable; judge convergence on these and on the identified
#' covariance ([fb_structured_cov()]) rather than on the raw loadings.
#'
#' @param fit A `flexybayes` greta fit with a `fa()` factor-analytic G x E
#'   term.
#' @param genotype_levels,environment_levels Optional character labels for
#'   the inner (genotype) and outer (environment) factors; default to
#'   positional labels.
#'
#' @return An `fb_met_summary` object (one entry is built per `fa()` term,
#'   the function returns the first / named): `op` (data frame of overall
#'   performance per genotype with credible interval), `stability` (data
#'   frame of across-environment spread per genotype), `gxe_blup` (the
#'   posterior-mean genotype-by-environment effect matrix), `env_cor` (the
#'   environment genetic-correlation matrix), `loadings` (posterior-mean
#'   factor loadings), and metadata.
#'
#' @seealso [fb_structured_cov()] for the identified environment covariance
#'   and its convergence diagnostic.
#' @examples
#' \dontrun{
#' fit <- flexybayes(
#'   yield ~ env, random = ~ fa(env, 2):gen, data = met, backend = "greta"
#' )
#' ms <- fb_met_summary(fit)
#' head(ms$op[order(-ms$op$mean), ]) # best genotypes on average
#' ms$env_cor # environment crossover structure
#' }
#' @export
fb_met_summary <- function(
  fit,
  genotype_levels = NULL,
  environment_levels = NULL
) {
  if (inherits(fit, c("flexybayes_inla", "flexybayes_brms"))) {
    stop(
      "fb_met_summary() needs a greta factor-analytic fit; the supplied fit ",
      "uses the ", if (inherits(fit, "flexybayes_inla")) "INLA" else "brms",
      " backend. ",
      "Breeder summaries are computed from the realised factor-analytic ",
      "effects, which fit on the greta backend; the INLA / brms MET path ",
      "reports variance components via summary() / fb_structured_cov(). ",
      "Refit fa(env, k):gen with backend = \"greta\" for fb_met_summary().",
      call. = FALSE
    )
  }
  if (!inherits(fit, "flexybayes")) {
    stop("`fit` must be a flexybayes (greta) object.", call. = FALSE)
  }
  rt <- fit$extras$parse_info$random %||% list()
  fa_terms <- Filter(function(t) identical(t$type %||% "", "fa_gxe"), rt)
  if (!length(fa_terms)) {
    stop(
      "fb_met_summary(): this fit carries no factor-analytic fa() G x E ",
      "term; nothing to summarise. Fit, e.g., random = ~ fa(env, 2):gen.",
      call. = FALSE
    )
  }
  draws <- tryCatch(fit$greta$draws, error = function(e) NULL)
  if (is.null(draws)) {
    stop(
      "fb_met_summary() requires the greta posterior draws; this fit does ",
      "not carry them (the factor-analytic term is greta-only).",
      call. = FALSE
    )
  }

  term <- fa_terms[[1L]]
  pooled <- do.call(rbind, lapply(draws, as.matrix))
  res <- .fb_met_one(pooled, draws, term, genotype_levels, environment_levels)
  structure(res, class = c("fb_met_summary", "list"))
}

# Summarise one factor-analytic term.
.fb_met_one <- function(pooled, draws, term, geno_levels, env_levels) {
  tag <- paste0(term$inner, "_", term$outer, "_fa", term$k)
  ni <- term$n_inner
  no <- term$n_outer

  geno_levels <- if (!is.null(geno_levels)) {
    as.character(geno_levels)
  } else {
    paste0(term$inner, seq_len(ni))
  }
  env_levels <- if (!is.null(env_levels)) {
    as.character(env_levels)
  } else {
    paste0(term$outer, seq_len(no))
  }

  # Map each g_mat column to its (genotype i, environment j) index.
  prefix <- paste0("g_mat_", tag)
  cols <- grep(paste0("^", prefix, "\\["), colnames(pooled), value = TRUE)
  if (length(cols) != ni * no) {
    stop(
      "fb_met_summary(): expected ", ni * no, " g_mat draw columns for ",
      "term '", tag, "'; found ", length(cols), ". The fit may predate the ",
      "g_mat monitoring (re-fit with the current fa() codegen).",
      call. = FALSE
    )
  }
  ij <- .fb_met_parse_ij(cols, prefix)

  # Per-draw genotype overall performance (mean across environments) and
  # stability (sd across environments), vectorised over draws.
  op_draws <- matrix(NA_real_, nrow(pooled), ni)
  stab_draws <- matrix(NA_real_, nrow(pooled), ni)
  for (i in seq_len(ni)) {
    ci <- cols[ij[, 1L] == i]
    block <- pooled[, ci, drop = FALSE]
    op_draws[, i] <- rowMeans(block)
    stab_draws[, i] <- apply(block, 1L, stats::sd)
  }

  op <- .fb_met_summ_df(op_draws, geno_levels, "genotype")
  stability <- .fb_met_summ_df(stab_draws, geno_levels, "genotype")
  names(stability)[names(stability) == "mean"] <- "mean"

  # G x E BLUP: posterior-mean realised effects, n_inner x n_outer.
  blup_vec <- colMeans(pooled[, cols, drop = FALSE])
  gxe_blup <- matrix(NA_real_, ni, no, dimnames = list(geno_levels, env_levels))
  gxe_blup[cbind(ij[, 1L], ij[, 2L])] <- blup_vec

  # Environment genetic correlation (the identified covariance fb_structured_cov
  # reports) -- the crossover / reverse-U structure.
  ident <- tryCatch(
    .fb_fa_identified_cov(draws, term),
    error = function(e) NULL
  )
  env_cor <- if (!is.null(ident)) {
    cm <- ident$cor_mean
    dimnames(cm) <- list(env_levels, env_levels)
    cm
  } else {
    NULL
  }

  # Posterior-mean loadings (n_outer x k), for the biplot.
  lambda_cols <- as.vector(vapply(
    seq_len(term$k),
    function(jj) sprintf("Lambda_%s[%d,%d]", tag, seq_len(no), jj),
    character(no)
  ))
  loadings <- NULL
  if (all(lambda_cols %in% colnames(pooled))) {
    loadings <- matrix(
      colMeans(pooled[, lambda_cols, drop = FALSE]), no, term$k,
      dimnames = list(env_levels, paste0("factor", seq_len(term$k)))
    )
  }

  list(
    op = op,
    stability = stability,
    gxe_blup = gxe_blup,
    env_cor = env_cor,
    loadings = loadings,
    k = term$k,
    n_genotypes = ni,
    n_environments = no,
    inner = term$inner,
    outer = term$outer
  )
}

# --- display ------------------------------------------------------ #

#' @exportS3Method print fb_met_summary
print.fb_met_summary <- function(x, ...) {
  cat("<fb_met_summary>  factor-analytic MET (", x$k, " factor",
    if (x$k != 1L) "s" else "", ")\n", sep = "")
  cat(
    "  ", x$n_genotypes, " genotypes x ", x$n_environments, " environments\n",
    sep = ""
  )
  best <- x$op[order(-x$op$mean), , drop = FALSE]
  cat("  top overall performers:\n")
  for (i in seq_len(min(3L, nrow(best)))) {
    cat(
      "    ", best$genotype[i], "  OP = ", format(round(best$mean[i], 3L)),
      " [", format(round(best$q2.5[i], 2L)), ", ",
      format(round(best$q97.5[i], 2L)), "]\n",
      sep = ""
    )
  }
  if (!is.null(x$env_cor)) {
    off <- x$env_cor[upper.tri(x$env_cor)]
    cat(
      "  environment genetic correlation: range [",
      format(round(min(off), 2L)), ", ", format(round(max(off), 2L)),
      "]", if (min(off) < 0) "  (negative = crossover G x E)" else "", "\n",
      sep = ""
    )
  }
  invisible(x)
}

# --- internal helpers --------------------------------------------- #

# Parse the (i, j) indices out of g_mat_<tag>[i,j] column names.
.fb_met_parse_ij <- function(cols, prefix) {
  inside <- sub(paste0("^", prefix, "\\[(.*)\\]$"), "\\1", cols)
  parts <- strsplit(inside, ",", fixed = TRUE)
  i <- as.integer(vapply(parts, `[[`, character(1), 1L))
  j <- as.integer(vapply(parts, `[[`, character(1), 2L))
  cbind(i = i, j = j)
}

# Per-genotype posterior summary data frame from a draws-by-genotype
# matrix.
.fb_met_summ_df <- function(draws_mat, labels, id_name) {
  out <- data.frame(
    label = labels,
    mean = colMeans(draws_mat),
    sd = apply(draws_mat, 2L, stats::sd),
    q2.5 = apply(draws_mat, 2L, stats::quantile, 0.025, names = FALSE),
    q97.5 = apply(draws_mat, 2L, stats::quantile, 0.975, names = FALSE),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "label"] <- id_name
  out
}
