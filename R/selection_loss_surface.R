# The held-out loss over the relative penalty and the elastic-net mixing
# parameter.
#
# One optimiser run walks one path on s = lambda / lambda_max at the mixing
# value its penalty declares, so a surface over (s, alpha) is several runs.
# fb_select() is not given a grid of alpha: choosing alpha would make it a
# second tuning coordinate of the decision, and every read-out built on one
# path (entrance, stability, the refit, the certificate) would need a rule
# for it. The display instead takes the runs the caller made, one per mixing
# value, and draws them only when they were scored on the same folds and
# the same score contract, walked the same grid on s and chose their point
# by the same rule. The chosen point marked is the one each run selected on
# its own path; the display chooses nothing.


# .fb_plot_plotly_installed() --- is plotly available.
#
# One line, kept as its own binding so the absent-plotly refusal has a
# triggering test on a machine where plotly is installed.
.fb_plot_plotly_installed <- function() {
  requireNamespace("plotly", quietly = TRUE)
}


#' The path one optimiser result walked on s, or the reason it has none
#'
#' @param x An `fb_select` object.
#' @return `list(ok, reason)`, and with `ok = TRUE` the fields the surface
#'   reads: `alpha`, `scores` (ordered by `lambda_index`), `selection`,
#'   `contract_hash`, `metric`, `rule`, `scoring_unit`, `columns` and
#'   `covariance`.
#' @noRd
#' @keywords internal
.fb_sel_surface_path <- function(x) {
  none <- function(reason) list(ok = FALSE, reason = reason)
  if (!inherits(x, "fb_select")) {
    return(none("it is not an fb_select() result"))
  }
  arm <- x$benchmark$glmnet
  scores <- arm$selection$scores
  if (!isTRUE(arm$selection$ran) || !is.data.frame(scores) ||
        nrow(scores) == 0L) {
    return(none("no penalised optimiser path was scored out of fold"))
  }
  if (!any(is.finite(scores$s))) {
    return(none(paste0(
      "this penalty has no finite null anchor, so its path is indexed by `",
      arm$selection$tuning_scale %||% "the penalty",
      "` and has no s = lambda / lambda_max coordinate"
    )))
  }
  alpha <- arm$support$alpha
  if (length(alpha) != 1L || !is.finite(alpha)) {
    return(none("its penalty has no elastic-net mixing parameter"))
  }
  hash <- unique(x$assessment_predictions$contract_hash)
  if (length(hash) != 1L || is.na(hash)) {
    return(none("it carries no verified score contract for its decision"))
  }
  scores <- as.data.frame(scores)
  list(
    ok = TRUE,
    alpha = as.numeric(alpha),
    scores = scores[order(scores$lambda_index), , drop = FALSE],
    selection = arm$selection,
    contract_hash = hash,
    metric = arm$selection$metric,
    rule = arm$selection$rule,
    scoring_unit = arm$scoring_unit,
    columns = arm$selection$candidate_availability$column,
    covariance = arm$covariance_source
  )
}


#' Refuse a set of results that does not span a loss surface
#'
#' @param why The reason, one sentence without its full stop.
#' @noRd
#' @keywords internal
.fb_sel_surface_refuse <- function(why) {
  stop(.fb_refusal_condition(
    reason_code = "selection_loss_surface_undefined",
    message = paste0(
      "plot(type = \"loss_surface\") draws the out-of-fold loss of optimiser ",
      "results at two or more elastic-net mixing values, scored on one score ",
      "contract: ", why, ". Run fb_select() once per mixing value with the ",
      "same base, candidates, resamples, seed and decision, and pass the ",
      "other results as `alongside = list(...)`."
    )
  ))
}


#' The loss surface over s and the elastic-net mixing parameter
#'
#' @param x An `fb_select` object; its selected point is the highlighted one.
#' @param alongside A list of `fb_select` objects at other mixing values.
#' @param ... Ignored.
#' @return The plotly widget, returned visibly, with attribute
#'   `values`: `surface` (one row per mixing value and grid point: `alpha`,
#'   `lambda_index`, `s`, `loss`, `se`), `chosen` (one row per result:
#'   `alpha`, `lambda_index`, `s`, `loss`, `primary`), `metric`,
#'   `scoring_unit` and `contract_hash`. `NULL` with a message when `x` has
#'   no path on s.
#' @noRd
#' @keywords internal
.fb_select_plot_loss_surface <- function(x, alongside = list(), ...) {
  if (!.fb_plot_plotly_installed()) {
    stop(.fb_refusal_condition(
      reason_code = "loss_surface_requires_plotly",
      message = paste0(
        "plot(type = \"loss_surface\") draws a three-dimensional surface ",
        "with plotly, which is not installed. Install it with ",
        "install.packages(\"plotly\"), or read the same values from ",
        "tidy(x) and plot(x, type = \"entrance_path\"), which draw one ",
        "mixing value on the base graphics device."
      )
    ))
  }
  primary <- .fb_sel_surface_path(x)
  if (!isTRUE(primary$ok)) {
    return(.fb_select_plot_unavailable("the loss surface", primary$reason))
  }
  if (!is.list(alongside) || inherits(alongside, "fb_select")) {
    .fb_sel_surface_refuse(
      "`alongside` must be a list of fb_select() results"
    )
  }

  # Every result on one contract, grid and rule ---------------------------
  paths <- list(primary)
  for (k in seq_along(alongside)) {
    path <- .fb_sel_surface_path(alongside[[k]])
    if (!isTRUE(path$ok)) {
      .fb_sel_surface_refuse(paste0(
        "entry ", k, " of `alongside` has no path on s: ", path$reason
      ))
    }
    fields <- c(
      "contract_hash", "metric", "rule", "scoring_unit", "columns",
      "covariance"
    )
    differing <- fields[!vapply(
      fields, function(f) identical(path[[f]], primary[[f]]), logical(1L)
    )]
    if (!identical(path$scores$s, primary$scores$s)) {
      differing <- c(differing, "grid on s")
    }
    if (length(differing) > 0L) {
      .fb_sel_surface_refuse(paste0(
        "entry ", k, " of `alongside` differs from `x` in ",
        paste(differing, collapse = ", ")
      ))
    }
    paths[[length(paths) + 1L]] <- path
  }
  alphas <- vapply(paths, function(p) p$alpha, numeric(1L))
  if (anyDuplicated(alphas) > 0L) {
    .fb_sel_surface_refuse(paste0(
      "the mixing value ", format(alphas[anyDuplicated(alphas)]),
      " appears more than once"
    ))
  }
  if (length(alphas) < 2L) {
    .fb_sel_surface_refuse("only one mixing value was given")
  }

  # The values, one row per mixing value ----------------------------------
  order_alpha <- order(alphas)
  paths <- paths[order_alpha]
  alphas <- alphas[order_alpha]
  is_primary <- order_alpha == 1L
  s_grid <- primary$scores$s
  z <- do.call(rbind, lapply(paths, function(p) p$scores$score))
  surface <- do.call(rbind, lapply(paths, function(p) {
    data.frame(
      alpha = p$alpha,
      lambda_index = as.integer(p$scores$lambda_index),
      s = p$scores$s,
      loss = p$scores$score,
      se = p$scores$se,
      stringsAsFactors = FALSE
    )
  }))
  rownames(surface) <- NULL
  chosen <- do.call(rbind, lapply(seq_along(paths), function(k) {
    p <- paths[[k]]
    at <- match(p$selection$selected_lambda_index, p$scores$lambda_index)
    data.frame(
      alpha = p$alpha,
      lambda_index = as.integer(p$scores$lambda_index[[at]]),
      s = p$scores$s[[at]],
      loss = p$scores$score[[at]],
      primary = is_primary[[k]],
      stringsAsFactors = FALSE
    )
  }))
  rownames(chosen) <- NULL

  # The drawing -------------------------------------------------------------
  unit <- paste(primary$scoring_unit %||% "unit", collapse = "+")
  subtitle <- paste0(
    "out-of-fold ", primary$metric, ", the mean over ", unit,
    "s of each ", unit, "'s value; diamond: the point this result chose"
  )
  hover <- paste0(
    "s = %{x:.3g}<br>alpha = %{y}<br>", primary$metric, " = %{z:.4g}",
    "<extra></extra>"
  )
  widget <- plotly::plot_ly()
  widget <- plotly::add_trace(
    widget,
    type = "surface",
    x = s_grid,
    y = alphas,
    z = z,
    colorscale = "Viridis",
    # Slightly translucent, so a chosen point on the surface is not hidden
    # by the surface in front of it.
    opacity = 0.9,
    colorbar = list(title = list(text = primary$metric)),
    hovertemplate = hover,
    name = "loss"
  )
  others <- chosen[!chosen$primary, , drop = FALSE]
  widget <- plotly::add_trace(
    widget,
    type = "scatter3d",
    mode = "markers",
    x = others$s,
    y = others$alpha,
    z = others$loss,
    marker = list(size = 5, color = "white",
                  line = list(color = "black", width = 2)),
    hovertemplate = hover,
    name = "chosen at the other alphas"
  )
  mine <- chosen[chosen$primary, , drop = FALSE]
  widget <- plotly::add_trace(
    widget,
    type = "scatter3d",
    mode = "markers",
    x = mine$s,
    y = mine$alpha,
    z = mine$loss,
    marker = list(size = 8, color = "black", symbol = "diamond"),
    hovertemplate = hover,
    name = "chosen point"
  )
  widget <- plotly::layout(
    widget,
    title = list(text = paste0(
      "Out-of-fold loss over s and the elastic-net mixing parameter",
      "<br><sub>", subtitle, "</sub>"
    )),
    scene = list(
      # The grid is geometric in s, so a logarithmic axis spaces it evenly;
      # strong penalty to weak, as on the entrance path.
      xaxis = list(
        title = list(text = "s = lambda / lambda_max"),
        type = "log",
        autorange = "reversed",
        dtick = 1
      ),
      yaxis = list(title = list(text = "alpha")),
      zaxis = list(title = list(text = primary$metric)),
      # Strong penalty on the left, alpha rising away from the viewer.
      camera = list(eye = list(x = 0.6, y = -2.0, z = 0.7))
    ),
    # One row under the subtitle, clear of it and of the axes.
    legend = list(
      orientation = "h", x = 0.5, xanchor = "center", y = 0.94,
      yanchor = "top"
    )
  )
  attr(widget, "values") <- list(
    surface = surface,
    chosen = chosen,
    metric = primary$metric,
    scoring_unit = primary$scoring_unit,
    contract_hash = primary$contract_hash
  )
  widget
}
