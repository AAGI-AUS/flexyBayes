# predict_posterior.R -- per-draw posterior aggregation + sampled-RE
# helpers used by the v0.3.5 file-output / sample-mode predict paths.
#
# v0.3.8 audit close: the legacy .predict_per_draw()
# entry point was superseded by .predict_linear_draws() in
# R/predict_kernel.R, which is the single per-draw kernel routed
# through by every prediction call site (in-memory smooth-aware,
# file-backed, sampled-RE). The helpers retained here are the
# downstream aggregation primitives and the sampled-RE realisation
# generator.
#
#   .predict_intervals(per_draw_mat, level = 0.95)
#     Returns list(point, lower, upper) with point = posterior mean
#     (rowMeans) and lower / upper = posterior quantiles at
#     (1 - level) / 2 and 1 - (1 - level) / 2.
#
#   .predict_sample_re_for_unknowns(object, unknown_summary,
#                                    n_draws_eff)
#     For each unknown level in each column with allow_new_levels =
#     "sample", samples (n_unknown_levels x n_draws_eff) random
#     effect realisations rnorm(0, tau_per_draw) where tau_per_draw
#     reads off the codegen-named tau_<group> draw column on
#     $greta$draws. Returns a list keyed by row index (global; not
#     chunk-relative), each entry a length-n_draws_eff vector. Uses
#     the caller's RNG state (set.seed() upstream is the
#     reproducibility hook).
#
#   .predict_in_memory_sample(object, newdata, type, se.fit)
#     In-memory wrapper for allow_new_levels = "sample". Builds the
#     per-row sampled-RE realisations via
#     .predict_sample_re_for_unknowns(), then routes through the
#     shared kernel (.predict_linear_draws) with
#     include = c("fixed", "smooth", "random_sampled"). Posterior
#     aggregation (rowMeans point + MC se.fit) at this layer.
#
# Not exported.

.predict_intervals <- function(per_draw_mat, level = 0.95) {
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop(
      ".predict_intervals(): `level` must be a numeric scalar in ",
      "(0, 1); got: ",
      deparse(level),
      call. = FALSE
    )
  }
  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)
  point <- rowMeans(per_draw_mat)
  q_mat <- apply(
    per_draw_mat,
    1L,
    stats::quantile,
    probs = probs,
    names = FALSE
  )
  # apply() on a single-row input returns a vector, not a matrix.
  if (is.null(dim(q_mat))) {
    q_mat <- matrix(q_mat, nrow = 2L)
  }
  list(
    point = unname(point),
    lower = unname(q_mat[1L, ]),
    upper = unname(q_mat[2L, ])
  )
}


.predict_sample_re_for_unknowns <- function(
  object,
  unknown_summary,
  n_draws_eff,
  n_draws_cap = 4000L
) {
  if (length(unknown_summary) == 0L) {
    return(list())
  }
  if (is.null(object$greta) || is.null(object$greta$draws)) {
    stop(
      ".predict_sample_re_for_unknowns(): allow_new_levels = ",
      "\"sample\" requires posterior draws on $greta$draws. ",
      "Fit does not carry them. INLA fits do not currently ",
      "expose per-draw access; an INLA draws adapter is queued ",
      "for a future release.",
      call. = FALSE
    )
  }

  draws_mat <- do.call(rbind, lapply(object$greta$draws, as.matrix))
  n_total <- nrow(draws_mat)
  if (n_total > n_draws_cap) {
    idx <- as.integer(seq.int(1L, n_total, length.out = n_draws_cap))
    draws_mat <- draws_mat[idx, , drop = FALSE]
  }

  out <- list()
  for (col in names(unknown_summary)) {
    info <- unknown_summary[[col]]
    sigma_col <- paste0("sigma_", col)
    tau_pat <- paste0("^tau_", col, "(\\[|$)")
    tau_cols <- grep(tau_pat, colnames(draws_mat), value = TRUE)

    # Resolve a length-n_draws_eff vector of group-level RE SDs. The
    # codegen names the group's RE precision-scale as `sigma_<col>`
    # for the per-row path and `tau_<col>[k]` for the indexed path.
    # The standard deviation for the unknown level is best
    # approximated by the prior SD `sigma_<col>` when available
    # (this is the population-level posterior on the RE
    # distribution's scale); fall back to the empirical SD of the
    # known-level tau draws if not.
    if (sigma_col %in% colnames(draws_mat)) {
      tau_per_draw <- draws_mat[, sigma_col]
    } else if (length(tau_cols) > 0L) {
      # Per-draw empirical SD across the known levels.
      tau_mat <- draws_mat[, tau_cols, drop = FALSE]
      tau_per_draw <- apply(tau_mat, 1L, stats::sd)
    } else {
      stop(
        ".predict_sample_re_for_unknowns(): no RE-scale draws ",
        "found for column `",
        col,
        "`. Looked for `",
        sigma_col,
        "` and `tau_",
        col,
        "[k]`. allow_new_levels = \"sample\" ",
        "requires the fit's codegen to expose the RE-scale ",
        "posterior; this fit does not.",
        call. = FALSE
      )
    }

    if (length(tau_per_draw) != n_draws_eff) {
      # Stride-subsample tau_per_draw to match n_draws_eff.
      n_t <- length(tau_per_draw)
      i <- as.integer(seq.int(1L, n_t, length.out = n_draws_eff))
      tau_per_draw <- tau_per_draw[i]
    }

    # For each unknown row index in this column, sample one RE
    # realisation per draw: rnorm(n_draws_eff, 0, tau_per_draw[k]).
    # Caller's set.seed() determines the realisations.
    row_indices <- info$row_indices
    for (i in row_indices) {
      key <- as.character(i)
      out[[key]] <- stats::rnorm(n_draws_eff, mean = 0, sd = tau_per_draw)
    }
  }

  out
}


# In-memory predict() helper for allow_new_levels = "sample" mode.
# Computes per-draw predictions including sampled-RE contributions
# on the unknown rows, then returns the posterior-mean point
# estimate (and optional posterior-SD se.fit) as a numeric vector
# (or list when se.fit = TRUE). Single-pass; the in-memory sample
# path does not chunk -- chunk-aware sample is in the file-output
# path (.predict_to_file).
.predict_in_memory_sample <- function(object, newdata, type, se.fit) {
  fam_link <- object$extras$parse_info$family

  sample_re_summary <- attr(newdata, "sample_re_summary")
  draws_count <- nrow(do.call(rbind, lapply(object$greta$draws, as.matrix)))
  n_draws_eff <- min(draws_count, 4000L)

  sampled_re <- if (!is.null(sample_re_summary)) {
    .predict_sample_re_for_unknowns(
      object = object,
      unknown_summary = sample_re_summary,
      n_draws_eff = n_draws_eff,
      n_draws_cap = 4000L
    )
  } else {
    NULL
  }

  # v0.3.8 audit Critical Fix #1: route through the
  # shared kernel so sample-mode fits with s() smooths layer the
  # smooth-basis contribution onto the linear predictor (pre-v0.3.8
  # .predict_per_draw() built mm from fixed_formula alone, silently
  # dropping the smooth basis).
  per_draw_lin <- tryCatch(
    .predict_linear_draws(
      object = object,
      newdata = newdata,
      n_draws_cap = 4000L,
      sampled_re = sampled_re,
      include = c("fixed", "smooth", "random_sampled")
    ),
    error = function(e) {
      if (
        grepl(
          "no fixed-effect draw columns",
          conditionMessage(e),
          fixed = TRUE
        ) ||
          grepl("no overlap between newdata", conditionMessage(e), fixed = TRUE)
      ) {
        return(NULL)
      }
      stop(e)
    }
  )
  if (is.null(per_draw_lin)) {
    return(
      if (isTRUE(se.fit)) {
        list(
          fit = rep(NA_real_, nrow(newdata)),
          se.fit = rep(NA_real_, nrow(newdata))
        )
      } else {
        rep(NA_real_, nrow(newdata))
      }
    )
  }

  per_draw <- if (identical(type, "response")) {
    .apply_link_inverse(per_draw_lin, fam_link$link)
  } else {
    per_draw_lin
  }

  point <- rowMeans(per_draw)
  if (isTRUE(se.fit)) {
    se <- apply(per_draw, 1L, stats::sd)
    return(list(fit = point, se.fit = se))
  }
  point
}
