# predict_file_output.R -- file-output for predict().
#
# Three internal helpers driving the v0.3.5 file-output extension on
# predict.flexybayes():
#
#   .fst_available()
#     requireNamespace("fst", quietly = TRUE) wrapper. Mockable in
#     tests via testthat::local_mocked_bindings().
#
#   .resolve_format(format, N, interop, fst_available)
#     Pure path function per the resolution table:
#       format = "auto" + interop          -> "csv"
#       format = "auto" + N >= 1e6 + fst   -> "fst"
#       format = "auto" + N >= 1e6 + !fst  -> "rds"
#       format = "auto" + N <  1e6         -> "rds"
#       format = "csv" / "rds"             -> identity
#       format = "fst"  + !fst             -> structured stop()
#
#   .predict_write_file(predictions, newdata, output_file, format)
#     Assembles the file-output data.table (point/lower/upper +
#     optional se.fit + newdata columns) and dispatches to
#     data.table::fwrite() / base::saveRDS() / fst::write_fst() per
#     the resolved format. Refuses if output_file exists -- no
#     silent overwrite, consistent with the package's broader
#     no-silent-refusal posture.
#
# Not exported.

.fst_available <- function() {
  requireNamespace("fst", quietly = TRUE)
}


.resolve_format <- function(format, N, interop, fst_available) {
  if (!is.character(format) || length(format) != 1L) {
    stop(
      ".resolve_format(): `format` must be a length-1 character; ",
      "got: ",
      deparse(format),
      call. = FALSE
    )
  }
  if (!format %in% c("auto", "csv", "rds", "fst")) {
    stop(
      ".resolve_format(): `format` must be one of ",
      "'auto', 'csv', 'rds', 'fst'; got: '",
      format,
      "'.",
      call. = FALSE
    )
  }
  if (!is.numeric(N) || length(N) != 1L || is.na(N) || N < 0) {
    stop(
      ".resolve_format(): `N` must be a non-negative numeric ",
      "scalar; got: ",
      deparse(N),
      call. = FALSE
    )
  }

  if (identical(format, "csv")) {
    return("csv")
  }
  if (identical(format, "rds")) {
    return("rds")
  }
  if (identical(format, "fst")) {
    if (isTRUE(fst_available)) {
      return("fst")
    }
    stop(
      "format = \"fst\" was requested but the 'fst' package is ",
      "not installed.\n\n",
      "Install with:\n",
      "  install.packages(\"fst\")\n\n",
      "Or override:\n",
      "  predict(fit, newdata, output_file = \"preds.rds\", ",
      "format = \"rds\")",
      call. = FALSE
    )
  }

  # format == "auto"
  if (isTRUE(interop)) {
    return("csv")
  }
  if (N >= 1e6 && isTRUE(fst_available)) {
    return("fst")
  }
  "rds"
}


.predict_write_file <- function(predictions, newdata, output_file, format) {
  if (file.exists(output_file)) {
    stop(
      ".predict_write_file(): `output_file` already exists at ",
      "`",
      output_file,
      "`. Refusing to overwrite silently. ",
      "Remove the file first (unlink(\"",
      output_file,
      "\")) ",
      "or pick a different path.",
      call. = FALSE
    )
  }

  out <- data.table::as.data.table(predictions)
  nd <- data.table::as.data.table(newdata)
  # Newdata columns are bound alongside predictions so the file
  # carries the prediction grid plus the per-row predictors. rds /
  # fst preserve factor type; csv coerces factor -> character per
  # csv's design.
  out <- cbind(out, nd)

  switch(
    format,
    "csv" = data.table::fwrite(out, file = output_file),
    "rds" = base::saveRDS(out, file = output_file),
    "fst" = {
      if (!requireNamespace("fst", quietly = TRUE)) {
        stop(
          ".predict_write_file(): format = \"fst\" but the 'fst' ",
          "package is not installed. .resolve_format() should ",
          "have caught this; check upstream.",
          call. = FALSE
        )
      }
      fst::write_fst(out, path = output_file)
    },
    stop(
      ".predict_write_file(): unsupported format `",
      format,
      "`.",
      call. = FALSE
    )
  )

  invisible(output_file)
}


# Top-level file-output dispatcher invoked by predict.flexybayes()
# when output_file is supplied. Splits newdata into chunks of
# `chunk_size` rows (or one pass when chunk_size is NULL),
# per-chunk computes the model.matrix + per-draw linear predictor
# (via .predict_per_draw) + posterior intervals (via
# .predict_intervals), accumulates point/lower/upper vectors
# across chunks, optionally accumulates se.fit, applies link
# transform per `type`, writes to disk via .predict_write_file
# under the resolved format. Returns invisible(output_file).
#
# For allow_new_levels = "sample": pre-samples sampled-RE
# realisations once at the outer level via
# .predict_sample_re_for_unknowns(), then dispatches sampled_re
# per-chunk by translating global row indices to chunk-relative
# indices. Same set of pre-sampled values used across chunks --
# preserves the bitwise-identical chunked-vs-
# unchunked file-output equivalence.
.predict_to_file <- function(
  object,
  newdata,
  output_file,
  format,
  interop,
  type,
  se.fit,
  chunk_size,
  allow_new_levels
) {
  N <- nrow(newdata)
  resolved_format <- .resolve_format(
    format = format,
    N = N,
    interop = interop,
    fst_available = .fst_available()
  )

  # Pre-sample RE realisations for unknown rows under
  # allow_new_levels = "sample". The pre-sample happens once at
  # outer level so chunked vs unchunked produces identical
  # numerical output for the same caller seed.
  sample_re_summary <- attr(newdata, "sample_re_summary")
  sampled_re_global <- NULL
  if (
    identical(allow_new_levels, "sample") &&
      !is.null(sample_re_summary)
  ) {
    # Use min(n_total_draws, 4000) for the draws cap; this
    # determines n_draws_eff for both .predict_sample_re_for_unknowns
    # and .predict_per_draw below.
    draws_count <- nrow(do.call(rbind, lapply(object$greta$draws, as.matrix)))
    n_draws_eff <- min(draws_count, 4000L)
    sampled_re_global <- .predict_sample_re_for_unknowns(
      object = object,
      unknown_summary = sample_re_summary,
      n_draws_eff = n_draws_eff,
      n_draws_cap = 4000L
    )
  }

  fam_link <- object$extras$parse_info$family

  # Chunk index sequence. chunk_size = NULL or chunk_size >= N
  # means a single pass over all rows.
  chunks <- if (is.null(chunk_size) || chunk_size >= N) {
    list(seq_len(N))
  } else {
    split(seq_len(N), ceiling(seq_len(N) / as.integer(chunk_size)))
  }

  point_acc <- numeric(N)
  lower_acc <- numeric(N)
  upper_acc <- numeric(N)
  se_acc <- if (isTRUE(se.fit)) numeric(N) else NULL

  for (rng in chunks) {
    chunk <- newdata[rng, , drop = FALSE]

    # Translate global row indices for sampled RE to chunk-relative
    # row indices. NULL when no sampled-RE is in play.
    sampled_re_chunk <- NULL
    if (!is.null(sampled_re_global)) {
      keys_global <- as.integer(names(sampled_re_global))
      hits <- which(keys_global %in% rng)
      if (length(hits) > 0L) {
        sampled_re_chunk <- setNames(
          sampled_re_global[hits],
          as.character(match(keys_global[hits], rng))
        )
      }
    }

    # v0.3.8 audit Critical Fix #1: route through the
    # shared kernel so fits with s() smooth random terms include the
    # smooth-basis contribution per draw. Pre-v0.3.8 the file-backed
    # path built model.matrix(fixed_formula, chunk) without consulting
    # extras$parse_info$smooths -- silently producing predictions
    # missing the smooth contribution on disk. The kernel routes both
    # paths through the same .linear_formula() / mgcv::PredictMat()
    # logic the in-memory smooth branch already used.
    per_draw_lin <- tryCatch(
      .predict_linear_draws(
        object = object,
        newdata = chunk,
        n_draws_cap = 4000L,
        sampled_re = sampled_re_chunk,
        include = c("fixed", "smooth", "random_sampled")
      ),
      error = function(e) {
        # No-overlap case: kernel found no fixed-effect overlap; mirror
        # the legacy NA-filled chunk semantics rather than aborting the
        # whole file-output.
        if (
          grepl(
            "no fixed-effect draw columns",
            conditionMessage(e),
            fixed = TRUE
          ) ||
            grepl(
              "no overlap between newdata",
              conditionMessage(e),
              fixed = TRUE
            )
        ) {
          return(NULL)
        }
        stop(e)
      }
    )
    if (is.null(per_draw_lin)) {
      point_acc[rng] <- NA_real_
      lower_acc[rng] <- NA_real_
      upper_acc[rng] <- NA_real_
      if (isTRUE(se.fit)) {
        se_acc[rng] <- NA_real_
      }
      next
    }

    # Apply link transform per draw before computing posterior
    # intervals. This gives an honest posterior expectation on the
    # response scale rather than inv_link(mean(linear_predictor)),
    # which differs for non-identity links.
    per_draw <- if (identical(type, "response")) {
      .apply_link_inverse(per_draw_lin, fam_link$link)
    } else {
      per_draw_lin
    }

    intervals <- .predict_intervals(per_draw, level = 0.95)
    point_acc[rng] <- intervals$point
    lower_acc[rng] <- intervals$lower
    upper_acc[rng] <- intervals$upper
    if (isTRUE(se.fit)) {
      se_acc[rng] <- apply(per_draw, 1L, stats::sd)
    }
  }

  predictions <- if (isTRUE(se.fit)) {
    list(
      point = point_acc,
      lower = lower_acc,
      upper = upper_acc,
      se.fit = se_acc
    )
  } else {
    list(point = point_acc, lower = lower_acc, upper = upper_acc)
  }

  .predict_write_file(
    predictions = predictions,
    newdata = newdata,
    output_file = output_file,
    format = resolved_format
  )
}


# Per-draw link-inverse on the linear-predictor matrix produced by
# .predict_per_draw(). Mirrors the family-link table in the legacy
# posterior-mean predict body.
.apply_link_inverse <- function(per_draw_lin, link) {
  switch(
    link,
    "identity" = per_draw_lin,
    "log" = exp(per_draw_lin),
    "logit" = 1 / (1 + exp(-per_draw_lin)),
    "probit" = stats::pnorm(per_draw_lin),
    per_draw_lin
  )
}
