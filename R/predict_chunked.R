# predict_chunked.R -- chunked-prediction core +
# file-output "sample" activation.
#
# Three internal helpers driving the v0.3.4 core + v0.3.5
# sample-activation extensions on predict.flexybayes():
#
#   .predict_resolve_factors(newdata, fb_dataset, allow_new_levels)
#     Dictionary-backed factor handling. Coerces each
#     factor / character column in newdata that has a dictionary on
#     <fb_dataset>$dictionaries to factor(levels = dictionary).
#     Detects unknown levels (values present in newdata that are not
#     in the dictionary). Routes per the allow_new_levels mode:
#       - "population" -- coerce to NA; rows surface in the predict()
#                         return with a warning naming the count.
#                         Population-level prediction (no RE
#                         contribution from the unknown group) is the
#                         downstream behaviour the existing predict
#                         body already implements when a factor level
#                         maps to NA in the model.matrix() expansion.
#       - "refuse"     -- structured stop naming column + first 5
#                         unknown levels + total count.
#       - "sample"     -- coerce to NA (same as "population") and
#                         attach attr(newdata, "sample_re_summary")
#                         carrying the unknown-row index information
#                         so downstream .predict_per_draw() can
#                         layer a sampled random-effect realisation
#                         onto the linear predictor for unknown
#                         rows. No warning fires (sampling is the
#                         requested semantics, not a fallback).
#                         Activated at v0.3.5 (was reserved at
#                         v0.3.4 with a forward-pointer refusal).
#
#   .predict_chunked_iterate(object, newdata, chunk_size, ...)
#     Chunk_size iteration. Splits newdata into chunks
#     of `chunk_size` rows, calls predict.flexybayes() per chunk
#     (with chunk_size = NULL so the recursion bottoms out), and
#     concatenates the per-chunk numeric vectors (or se.fit lists)
#     into a single result of the same shape as the all-in-one call.
#
# Used by:
#   - predict.flexybayes() (R/methods.R) -- the chunk_size +
#     allow_new_levels + output_file arguments route through these
#     helpers.
#
# Not exported.

.predict_resolve_factors <- function(newdata, fb_dataset, allow_new_levels) {
  if (
    is.null(fb_dataset) ||
      is.null(fb_dataset$dictionaries) ||
      length(fb_dataset$dictionaries) == 0L
  ) {
    return(newdata)
  }

  if (!is.data.frame(newdata)) {
    stop(
      ".predict_resolve_factors() expects a data.frame newdata; ",
      "got: ",
      paste(class(newdata), collapse = "/"),
      call. = FALSE
    )
  }

  unknown_summary <- list()
  for (col in intersect(names(fb_dataset$dictionaries), names(newdata))) {
    fit_levels <- fb_dataset$dictionaries[[col]]
    col_vals <- as.character(newdata[[col]])
    unknown_idx <- which(
      !(col_vals %in% fit_levels) &
        !is.na(col_vals)
    )

    if (length(unknown_idx) > 0L) {
      unknown_levels <- unique(col_vals[unknown_idx])
      unknown_summary[[col]] <- list(
        n_rows = length(unknown_idx),
        unknown_levels = unknown_levels,
        row_indices = unknown_idx
      )

      if (identical(allow_new_levels, "refuse")) {
        show <- if (length(unknown_levels) > 5L) {
          c(unknown_levels[1:5], "...")
        } else {
          unknown_levels
        }
        stop(
          "predict.flexybayes(allow_new_levels = \"refuse\"): ",
          "column `",
          col,
          "` carries level(s) not present in the ",
          "fit-time dictionary. ",
          length(unknown_idx),
          " row(s) affected across ",
          length(unknown_levels),
          " unknown level(s): ",
          paste(show, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
      # "population" and "sample" branches both set unknown rows to
      # NA in this column. Downstream model.matrix() drops NA rows
      # from the design matrix (na.action default); the predict
      # body fills the corresponding prediction slots with NA. For
      # "sample", the sample_re_summary attribute below carries the
      # row indices so downstream per-draw computation can add a
      # sampled RE realisation per unknown row.
      col_vals[unknown_idx] <- NA_character_
    }

    newdata[[col]] <- factor(col_vals, levels = fit_levels)
  }

  if (length(unknown_summary) > 0L) {
    total_rows <- sum(vapply(
      unknown_summary,
      function(x) x$n_rows,
      integer(1L)
    ))
    n_cols <- length(unknown_summary)

    if (identical(allow_new_levels, "sample")) {
      # No warning -- sampling is the requested semantics. The
      # downstream per-draw computation (R/predict_posterior.R)
      # reads the attribute and adds a fresh
      # Normal(0, tau_<col>_per_draw) realisation per unknown row
      # per draw.
      attr(newdata, "sample_re_summary") <- unknown_summary
    } else {
      # "population" branch: surface the count via warning + a
      # caller-introspectable attribute. Population-level
      # prediction returns NA for these rows.
      warning(
        "predict.flexybayes(allow_new_levels = \"",
        allow_new_levels,
        "\"): ",
        total_rows,
        " row(s) carry unknown factor level(s) across ",
        n_cols,
        " column(s) (",
        paste(names(unknown_summary), collapse = ", "),
        "); set to NA. ",
        "Population-level prediction returns NA for these rows. ",
        "Pass `allow_new_levels = \"refuse\"` to stop ",
        "explicitly on unknown levels, or ",
        "`allow_new_levels = \"sample\"` to layer a sampled ",
        "random-effect realisation onto each unknown row.",
        call. = FALSE
      )
      attr(newdata, "n_unknown_levels") <- unknown_summary
    }
  }

  newdata
}


.predict_chunked_iterate <- function(
  object,
  newdata,
  chunk_size,
  type,
  se.fit,
  allow_new_levels,
  ...
) {
  n <- nrow(newdata)
  if (n <= chunk_size) {
    return(predict.flexybayes(
      object,
      newdata = newdata,
      type = type,
      se.fit = se.fit,
      chunk_size = NULL,
      allow_new_levels = allow_new_levels,
      ...
    ))
  }

  n_chunks <- ceiling(n / chunk_size)
  fit_acc <- vector("list", n_chunks)
  se_acc <- if (isTRUE(se.fit)) vector("list", n_chunks) else NULL

  for (k in seq_len(n_chunks)) {
    lo <- (k - 1L) * chunk_size + 1L
    hi <- min(k * chunk_size, n)
    chunk <- newdata[lo:hi, , drop = FALSE]
    out <- predict.flexybayes(
      object,
      newdata = chunk,
      type = type,
      se.fit = se.fit,
      chunk_size = NULL,
      allow_new_levels = allow_new_levels,
      ...
    )
    if (isTRUE(se.fit)) {
      fit_acc[[k]] <- out$fit
      se_acc[[k]] <- out$se.fit
    } else {
      fit_acc[[k]] <- out
    }
  }

  if (isTRUE(se.fit)) {
    list(
      fit = unlist(fit_acc, use.names = FALSE),
      se.fit = unlist(se_acc, use.names = FALSE)
    )
  } else {
    unlist(fit_acc, use.names = FALSE)
  }
}
