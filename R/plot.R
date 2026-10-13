# Plot methods for flexyBayes

#' Plot diagnostics for a flexyBayes model
#'
#' Draws one of eight standard displays for a fitted model, chosen by
#' `type`. Every backend reaches the same method, and a display that needs
#' a slot a given engine does not populate (posterior draws on an INLA
#' fit, for instance) prints a message naming what
#' is missing and returns invisibly instead of erroring or drawing an
#' empty panel. The one exception is `"pp_check"`, which raises the
#' refusal [pp_check.flexybayes()] raises, so a caller can catch it by
#' class.
#'
#' The default display depends on what the fit carries. A sampled
#' posterior gets `"diagnostics"`, because the first question about a
#' sampler is whether it converged. A deterministic approximation has no
#' chains to trace, so an INLA fit gets `"residuals"` instead of a
#' message declining to draw the display it was asked for.
#'
#' @param x A fitted `flexybayes` object of any backend. Aggregated and
#'   generalised-linear fits reach the same method
#'   through their own registrations.
#' @param type A single string naming the display to draw. One of:
#'   `"diagnostics"`, trace plots and marginal densities per parameter
#'   (needs a sampled posterior), `"residuals"`, residuals against
#'   fitted values beside a normal quantile-quantile plot, `"effects"`, a
#'   forest plot of the fixed effects with their credible intervals,
#'   available on every backend that supplies [coef()] and [confint()].
#'   `"variance"`, a bar chart of the variance components with credible
#'   intervals, `"blups"`, a caterpillar plot of the random-effect
#'   predictions ordered by magnitude, `"pp_check"`, a posterior
#'   predictive check overlaying replicated datasets on the observed
#'   response, which needs predictive draws and so refuses by name on a
#'   fit that carries none (see [pp_check.flexybayes()]), and
#'   `"variogram"`, the empirical semivariance of the residuals over the
#'   design index, and `"spatial_field"`, the posterior mean of a separable
#'   autoregressive field, `random = ~ ar1(row):ar1(col)`, drawn on the
#'   row-by-column array of the trial on the response scale. `NULL` (the
#'   default) resolves to `"diagnostics"` on a sampled fit and
#'   `"residuals"` otherwise.
#' @param ... Further arguments passed to the underlying plotting call,
#'   for example `variable` to restrict which parameters a diagnostic
#'   display covers, or `type` and `ndraws` for `"pp_check"`. The spatial
#'   field takes `trial`, the identifier printed in its title and the only
#'   label it adds to the array's own row and column labels. On a fit with
#'   one field per trial, `random = ~ id(trial):ar1(row):ar1(col)`, it names
#'   the level whose field is drawn.
#' @returns Invisibly, the object the underlying plotting call returns:
#'   a \pkg{ggplot2} object for the displays built with it, the
#'   semivariance table for `"variogram"`, for `"spatial_field"` a list
#'   holding `values` (the field, one row per row level and one column per
#'   column level), `row_var`, `col_var`, `trial`, `title` and `spelling`
#'   (the field's term as written), and `NULL`
#'   for those drawn on the base graphics device. A fit with no separable
#'   field refuses `"spatial_field"` with
#'   `spatial_field_requires_ar1_field`. Called for the plot it draws.
#' @seealso [pp_check.flexybayes()] for the posterior predictive check
#'   `type = "pp_check"` draws, and the refusal it raises where a fit has
#'   no predictive draws.
#' @export
plot.flexybayes <- function(x, type = NULL, ...) {
  type <- if (is.null(type)) {
    .fb_default_plot_type(x)
  } else {
    match.arg(type, .FB_PLOT_TYPES)
  }

  switch(
    type,
    "diagnostics" = .plot_diagnostics(x, ...),
    "residuals" = .plot_residuals(x, ...),
    "effects" = .plot_effects(x, ...),
    "variance" = .plot_variance(x, ...),
    "blups" = .plot_blups(x, ...),
    "pp_check" = .plot_pp_check(x, ...),
    "variogram" = .plot_variogram(x, ...),
    "spatial_field" = .plot_spatial_field(x, ...)
  )
}

# The closed vocabulary of displays.
.FB_PLOT_TYPES <- c(
  "diagnostics",
  "residuals",
  "effects",
  "variance",
  "blups",
  "pp_check",
  "variogram",
  "spatial_field"
)

# .fb_has_sampler_draws(): did an engine draw samples for this fit?
#
# A brms fit keeps its draws inside the brmsfit at `$brms`. A nested
# Laplace approximation has none at all.
#
# @noRd
# @keywords internal
.fb_has_sampler_draws <- function(x) {
  inherits(x$brms, "brmsfit")
}

# .fb_default_plot_type(): the display a fit earns by what it holds.
#
# @noRd
# @keywords internal
.fb_default_plot_type <- function(x) {
  if (.fb_has_sampler_draws(x)) {
    return("diagnostics")
  }
  "residuals"
}

# Registered for each fit class instead of relying on inheritance. The
# brms and INLA fits do share the "flexybayes" parent as of 0.9.0, so the
# parent registration would reach them, but the aggregated and glm
# classes carry their own dispatch order and an explicit
# registration keeps the set visible in one place. The shared body above
# is backend-aware: a display that reads a slot a given backend does not
# populate (MCMC draws on an INLA fit, for instance) degrades to an
# informative message via .plot_unavailable() instead of erroring.
# "effects" works on every backend that exposes coef() and confint().

#' @rdname plot.flexybayes
#' @export
plot.flexybayes_inla <- function(x, ...) plot.flexybayes(x, ...)

#' @rdname plot.flexybayes
#' @export
plot.flexybayes_brms <- function(x, ...) plot.flexybayes(x, ...)

#' @rdname plot.flexybayes
#' @export
plot.flexybayes_aggregated <- function(x, ...) plot.flexybayes(x, ...)

#' @rdname plot.flexybayes
#' @export
plot.flexybayes_glm <- function(x, ...) plot.flexybayes(x, ...)

# Emit a non-silent, non-erroring notice that a plot type is not
# available for this fit's backend, then return invisibly. Backends
# differ in which slots they expose (brms carries MCMC draws inside its
# brmsfit, INLA does not), so a plot type
# that reads a slot the backend never populates degrades to a message
# instead of crashing through to graphics::plot.default().
.plot_unavailable <- function(type, reason) {
  message("plot(type = \"", type, "\") is not available for this fit: ", reason)
  invisible(NULL)
}

# MCMC diagnostics: trace + density
.plot_diagnostics <- function(x, ...) {
  if (!.fb_has_sampler_draws(x)) {
    return(.plot_unavailable(
      "diagnostics",
      paste0(
        "trace and density plots describe a sampler, and this fit carries ",
        "no sampler draws. A nested Laplace approximation reports mode ",
        "status and a Kullback-Leibler divergence instead (see ",
        "summary(fit)$converge) and plot(fit) defaults to the residual ",
        "panels on such a fit."
      )
    ))
  }

  # A brms fit keeps its draws inside the brmsfit, and brms draws its own
  # panels from them. Forwarding is what makes plot(brms_fit) show the
  # trace and density it always claimed to. This is the only sampler
  # this branch can reach (.fb_has_sampler_draws() already confirmed
  # `x$brms` is a `brmsfit`), so brms is required here, not merely
  # preferred.
  if (!requireNamespace("brms", quietly = TRUE)) {
    return(.plot_unavailable(
      "diagnostics",
      "the fit's draws live in a brmsfit and brms is not installed."
    ))
  }
  dots <- list(...)
  if (is.null(dots$combo)) {
    dots$combo <- c("dens", "trace")
  }
  if (is.null(dots$ask)) {
    dots$ask <- FALSE
  }
  invisible(do.call(
    graphics::plot,
    c(list(x = x$brms), dots)
  ))
}

# Residual diagnostics
.plot_residuals <- function(x, ...) {
  fitted_vals <- tryCatch(stats::fitted(x), error = function(e) NULL)
  resid_vals <- tryCatch(stats::residuals(x), error = function(e) NULL)

  if (
    is.null(fitted_vals) ||
      is.null(resid_vals) ||
      length(fitted_vals) == 0L ||
      length(resid_vals) == 0L
  ) {
    return(.plot_unavailable(
      "residuals",
      "fitted values and residuals are not available for this backend."
    ))
  }

  if (any(is.na(fitted_vals)) || any(is.na(resid_vals))) {
    message("Cannot plot residuals: fitted values contain NA.")
    return(invisible(NULL))
  }

  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  on.exit(par(old_par))

  # 1. Residuals vs Fitted
  plot(
    fitted_vals,
    resid_vals,
    xlab = "Fitted values",
    ylab = "Residuals",
    main = "Residuals vs Fitted",
    pch = 16,
    col = "#00000060"
  )
  abline(h = 0, lty = 2, col = "red")
  lines(lowess(fitted_vals, resid_vals), col = "blue")

  # 2. QQ plot
  qqnorm(resid_vals, main = "Normal Q-Q", pch = 16, col = "#00000060")
  qqline(resid_vals, col = "red")

  # 3. Scale-Location
  plot(
    fitted_vals,
    sqrt(abs(resid_vals)),
    xlab = "Fitted values",
    ylab = "sqrt(|Residuals|)",
    main = "Scale-Location",
    pch = 16,
    col = "#00000060"
  )
  lines(lowess(fitted_vals, sqrt(abs(resid_vals))), col = "blue")

  # 4. Histogram of residuals
  hist(
    resid_vals,
    breaks = 30,
    main = "Residual Distribution",
    xlab = "Residuals",
    col = "lightblue",
    border = "white"
  )

  invisible(NULL)
}

# Forest plot of fixed effects
.plot_effects <- function(x, level = 0.95, ...) {
  beta <- tryCatch(stats::coef(x), error = function(e) NULL)
  if (is.null(beta) || length(beta) == 0L) {
    message("No fixed effects to plot.")
    return(invisible(NULL))
  }

  ci <- tryCatch(stats::confint(x, level = level), error = function(e) NULL)
  if (is.null(ci) || nrow(ci) != length(beta)) {
    return(.plot_unavailable(
      "effects",
      "fixed-effect credible intervals are not available for this backend."
    ))
  }

  n <- length(beta)
  y_pos <- seq_len(n)

  old_par <- par(mar = c(4, max(nchar(names(beta))) * 0.6 + 2, 2, 1))
  on.exit(par(old_par))

  xlim <- range(c(ci, 0)) * 1.1

  plot(
    beta,
    y_pos,
    xlim = xlim,
    yaxt = "n",
    xlab = "Estimate",
    ylab = "",
    main = paste0("Fixed effects (", round(level * 100), "% CrI)"),
    pch = 16,
    cex = 1.2
  )
  segments(ci[, 1], y_pos, ci[, 2], y_pos, lwd = 2)
  abline(v = 0, lty = 2, col = "grey50")
  axis(2, at = y_pos, labels = names(beta), las = 1, cex.axis = 0.8)

  invisible(NULL)
}

# Variance components plot
.plot_variance <- function(x, ...) {
  # Through the shared reader, because the aggregated emits record the
  # field as a list of posterior means instead of the five-column table
  # this panel needs, and `nrow()` of that list is NULL.
  vc <- .fb_variance_comps(x)
  if (is.null(vc) || nrow(vc) == 0) {
    message("No variance components to plot.")
    return(invisible(NULL))
  }

  n <- nrow(vc)
  y_pos <- seq_len(n)

  old_par <- par(mar = c(4, max(nchar(vc$component)) * 0.5 + 2, 2, 1))
  on.exit(par(old_par))

  # The estimates join the upper bounds in the range so a table that
  # carries no intervals still gets an axis, instead of -Inf.
  x_upper <- suppressWarnings(
    max(c(vc$q97.5, vc$estimate), na.rm = TRUE)
  )
  if (!is.finite(x_upper) || x_upper <= 0) {
    x_upper <- 1
  }
  # A correlation row can be negative, and an axis pinned at zero would
  # draw it off the panel.
  x_lower <- suppressWarnings(min(c(0, vc$q2.5, vc$estimate), na.rm = TRUE))
  xlim <- c(x_lower * 1.1, x_upper * 1.1)

  plot(
    vc$estimate,
    y_pos,
    xlim = xlim,
    yaxt = "n",
    xlab = "Estimate",
    ylab = "",
    main = "Variance components (95% CrI)",
    pch = 16,
    cex = 1.2
  )
  segments(vc$q2.5, y_pos, vc$q97.5, y_pos, lwd = 2)
  axis(2, at = y_pos, labels = vc$component, las = 1, cex.axis = 0.7)

  invisible(NULL)
}

# Caterpillar plot of BLUPs
.plot_blups <- function(x, ...) {
  blups <- x$extras$blups
  if (length(blups) == 0) {
    message("No BLUPs available.")
    return(invisible(NULL))
  }

  n_terms <- length(blups)
  old_par <- par(mfrow = c(1, n_terms), mar = c(4, 4, 2, 1))
  on.exit(par(old_par))

  for (nm in names(blups)) {
    vals <- blups[[nm]]
    n <- length(vals)
    ord <- order(vals)

    plot(
      vals[ord],
      seq_len(n),
      pch = 16,
      cex = 0.6,
      xlab = "BLUP",
      ylab = "Rank",
      main = nm,
      col = "#00000080"
    )
    abline(v = 0, lty = 2, col = "red")
  }

  invisible(NULL)
}

# Posterior predictive check.
#
# The display this type names is a check of the model against its own
# replicated data, and it is drawn wherever the fit carries predictive
# draws to replicate from: which is the brms engine. The seam is
# pp_check.flexybayes() (R/pp_check_support.R), so the plot() entry
# point and the bayesplot generic cannot disagree about what a check
# is, and a fit with no predictive draws refuses by name here too.
#
# Before 0.9.1 this drew observed values against fitted values, plus a
# density of the fitted values labelled Posterior Predictive, on every
# backend that carried a response and a fitted vector. No replicated
# dataset was ever involved. The panel is gone instead of retitled:
# what it showed, residuals against fitted values, is
# plot(fit, type = "residuals").
.plot_pp_check <- function(x, ...) {
  p <- .fb_pp_check_impl(x, ...)
  print(p)
  invisible(p)
}


# ---------------------------------------------------------------- #
# The empirical residual variogram                                  #
# ---------------------------------------------------------------- #
#
# The display a field trialist reaches for after fitting a spatial
# model: how the residual semivariance behaves with separation along the
# design index. It is an EMPIRICAL variogram of the residuals, no
# fitted variogram is overlaid and none is claimed, so a flat surface
# here says the fitted structure absorbed what was there and nothing
# stronger.
#
# Observed rows only. A row carried as a latent design cell has no
# observed response and therefore no residual, and pairing an NA into a
# squared difference would silently drop whole lags instead of one
# pair. The subtitle prints both counts so the reader can see how much
# of the array the picture rests on.

# .fb_variogram_index(): the design index a fit is built over.
#
# @noRd
# @keywords internal
.fb_variogram_index <- function(x) {
  vars <- x$extras$na_action$design_index_vars %||% character(0)
  if (length(vars) == 0L) {
    vars <- .fb_design_index_vars(.fb_fit_ir(x) %||% list())
  }
  unique(vars[!is.na(vars) & nzchar(vars)])
}

# .fb_variogram_position(): one index column as an ordinal position.
#
# A design index is a lattice: what a lag counts is steps along the
# array, not the distance between two label strings. A factor is coded
# by its sorted level order, a numeric column is taken at face value.
#
# @noRd
# @keywords internal
.fb_variogram_position <- function(v) {
  if (is.numeric(v)) {
    return(as.numeric(v))
  }
  as.numeric(factor(as.character(v), levels = sort(unique(as.character(v)))))
}

# .fb_variogram_table(): semivariance by lag.
#
# gamma(h) = mean over pairs at separation h of half the squared
# residual difference. Returned as a data frame so a caller can read the
# numbers the picture is drawn from.
#
# @noRd
# @keywords internal
.fb_variogram_table <- function(resid, positions, index_vars) {
  n <- length(resid)
  diffs <- outer(resid, resid, "-")
  half_sq <- (diffs^2) / 2
  lags <- lapply(positions, function(p) abs(outer(p, p, "-")))

  upper <- upper.tri(half_sq)
  key <- do.call(
    paste,
    c(lapply(lags, function(m) round(m[upper], 6L)), list(sep = "\r"))
  )
  gamma <- tapply(half_sq[upper], key, mean)
  count <- tapply(half_sq[upper], key, length)

  parts <- strsplit(names(gamma), "\r", fixed = TRUE)
  out <- as.data.frame(
    lapply(seq_along(index_vars), function(k) {
      as.numeric(vapply(parts, function(p) p[[k]], character(1L)))
    }),
    stringsAsFactors = FALSE
  )
  names(out) <- paste0("lag_", index_vars)
  out$semivariance <- as.numeric(gamma)
  out$n_pairs <- as.integer(count)
  out <- out[do.call(order, out[paste0("lag_", index_vars)]), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "n_observed") <- n
  out
}

# .fb_variogram_draw(): the picture, on the base device this file uses.
#
# The bottom margin is five lines, not the four this file's other panels
# use, because `sub =` is drawn at line `par("mgp")[1] + 1` of the bottom
# margin, line 4 with the default mgp, and a four-line margin does not
# reach it. The counts were computed, carried on the returned table, and
# drawn outside the device, so the caption never appeared for a reader.
#
# @noRd
# @keywords internal
.fb_variogram_draw <- function(tab, index_vars, n_observed, n_design) {
  title <- "Empirical residual variogram"
  subtitle <- sprintf(
    "%d observed / %d design rows", n_observed, n_design
  )
  # Carried on the returned table as well as printed, so the two cannot
  # be checked separately and drift.
  attr(tab, "title") <- title
  attr(tab, "subtitle") <- subtitle
  attr(tab, "n_design") <- as.integer(n_design)
  lag_cols <- paste0("lag_", index_vars)

  if (length(index_vars) == 1L) {
    old_par <- graphics::par(mar = c(5, 4, 3, 1))
    on.exit(graphics::par(old_par))
    plot(
      tab[[lag_cols[[1L]]]],
      tab$semivariance,
      type = "b",
      pch = 16,
      xlab = paste0("lag (", index_vars[[1L]], ")"),
      ylab = "Semivariance",
      main = title,
      sub = subtitle
    )
    return(invisible(tab))
  }

  # Two or more index variables: the first two carry the surface, which
  # is the row x column array a field trial is laid out on.
  ux <- sort(unique(tab[[lag_cols[[1L]]]]))
  uy <- sort(unique(tab[[lag_cols[[2L]]]]))
  z <- matrix(NA_real_, nrow = length(ux), ncol = length(uy))
  z[cbind(
    match(tab[[lag_cols[[1L]]]], ux),
    match(tab[[lag_cols[[2L]]]], uy)
  )] <- tab$semivariance

  old_par <- graphics::par(mar = c(5, 4, 3, 1))
  on.exit(graphics::par(old_par))
  # The palette is graphics::image()'s own sequential default, which
  # keeps grDevices off the Imports list for one call.
  graphics::image(
    x = ux,
    y = uy,
    z = z,
    xlab = paste0("lag (", index_vars[[1L]], ")"),
    ylab = paste0("lag (", index_vars[[2L]], ")"),
    main = title,
    sub = subtitle
  )
  ok <- is.finite(z)
  if (sum(ok) > 3L && diff(range(z[ok])) > 0) {
    graphics::contour(x = ux, y = uy, z = z, add = TRUE, col = "#00000060")
  }
  invisible(tab)
}

# .plot_variogram(): the type = "variogram" entry point.
#
# @noRd
# @keywords internal
.FB_VARIOGRAM_MAX_OBS <- 4000L

.plot_variogram <- function(x, ...) {
  index_vars <- .fb_variogram_index(x)
  dat <- .fb_fit_data(x)

  if (length(index_vars) == 0L || is.null(dat) ||
    !all(index_vars %in% names(dat))) {
    stop(.fb_refusal_condition(
      reason_code = "variogram_requires_design_index",
      message = paste0(
        "plot(type = \"variogram\") needs a design index to measure ",
        "separation along, and this fit carries none",
        if (length(index_vars) > 0L) {
          paste0(
            " (the recorded index ", paste(index_vars, collapse = ", "),
            " is not a column of the fitted data)"
          )
        } else {
          ": its terms name no row or column array"
        },
        ". A residual variogram is a picture of how far apart two plots ",
        "are, so a model with no such array has no lag to plot against. ",
        "Fit a spatial term (random = ~ ar1(row):ar1(col)) or use ",
        "plot(fit, type = \"residuals\") for the fitted-value diagnostics."
      ),
      family_class = "flexybayes_variogram_requires_design_index"
    ))
  }

  resid_all <- tryCatch(stats::residuals(x), error = function(e) NULL)
  if (is.null(resid_all) || length(resid_all) != nrow(dat)) {
    return(.plot_unavailable(
      "variogram",
      "residuals are not available for this backend."
    ))
  }

  # Observed rows only: an augmented design cell has no observed response
  # and so no residual.
  keep <- which(!is.na(resid_all))
  if (length(keep) < 3L) {
    return(.plot_unavailable(
      "variogram",
      paste0(
        "only ", length(keep), " row(s) carry an observed residual, which ",
        "is too few to form a semivariance."
      )
    ))
  }
  if (length(keep) > .FB_VARIOGRAM_MAX_OBS) {
    return(.plot_unavailable(
      "variogram",
      paste0(
        length(keep), " observed rows would need ",
        format(length(keep)^2, big.mark = " ", scientific = FALSE),
        " pairwise separations to be held at once. Plot a subset of the ",
        "array, or read the fitted correlations from summary(fit)$varcomp."
      )
    ))
  }

  positions <- lapply(index_vars, function(v) {
    .fb_variogram_position(dat[[v]][keep])
  })
  tab <- .fb_variogram_table(resid_all[keep], positions, index_vars)
  .fb_variogram_draw(tab, index_vars, length(keep), nrow(dat))
}


# ---------------------------------------------------------------- #
# The posterior mean of a separable spatial field                   #
# ---------------------------------------------------------------- #
#
# A Stage 1 fit of one trial carries its spatial trend as the latent field
# of random = ~ ar1(row):ar1(col), which INLA writes as
# f(<row>_id, model = "ar1", group = <col>_id) (R/emit_inla.R). The field's
# posterior mean is stored node by node in summary.random[["<row>_id"]],
# mapped back to the response scale when the fit was built. Node k holds
# row index ((k - 1) mod n_row) + 1 and column index
# floor((k - 1) / n_row) + 1, and with one field per trial,
# id(trial):ar1(row):ar1(col), each trial's nodes follow the previous
# trial's. The indices are the positions of the row and column labels in
# factor(), as the emitter builds them, so the display is labelled with the
# array's own labels and with nothing read from any other column.

# .fb_spatial_field(): the field of one trial as a row x column matrix.
#
# @noRd
# @keywords internal
.fb_spatial_field <- function(x, trial = NULL) {
  terms <- Filter(
    function(t) identical(t$type %||% "", "ar1_spatial"),
    .fb_fit_ir(x)$random_terms %||% list()
  )
  if (length(terms) == 0L) {
    stop(.fb_refusal_condition(
      reason_code = "spatial_field_requires_ar1_field",
      message = paste0(
        "plot(type = \"spatial_field\") draws the posterior mean of a ",
        "separable spatial field, and this fit carries none: its random ",
        "terms name no row-by-column array. Fit one with random = ~ ",
        "ar1(row):ar1(col), or use plot(fit, type = \"variogram\") on a ",
        "fit with a design index."
      )
    ))
  }
  if (length(terms) > 1L) {
    stop(
      "plot(type = \"spatial_field\") draws one field, and this fit carries ",
      length(terms), " separable fields (",
      paste(vapply(terms, .ar1_term_spelling, character(1L)), collapse = ", "),
      ").",
      call. = FALSE
    )
  }
  term <- terms[[1L]]
  if (!is.null(trial) &&
        (!is.character(trial) || length(trial) != 1L || is.na(trial))) {
    stop("`trial` must be a single string.", call. = FALSE)
  }

  # The index each label was given, as the emitter gave it -------------------
  data <- x$data
  level_of <- function(v) {
    .inla_restore_level_labels(
      x$level_labels, v, levels(factor(data[[v]]))
    )
  }
  row_levels <- level_of(term$row_var)
  col_levels <- level_of(term$col_var)
  n_row <- length(row_levels)
  n_col <- length(col_levels)
  offset <- 0L
  if (!is.null(term$at_var)) {
    at_levels <- level_of(term$at_var)
    if (is.null(trial) || !trial %in% at_levels) {
      stop(
        "`trial` must name the level of `", term$at_var, "` whose field ",
        "is drawn: this fit carries one field for each of ",
        paste(at_levels, collapse = ", "), ".",
        call. = FALSE
      )
    }
    offset <- n_row * n_col * (match(trial, at_levels) - 1L)
  }

  # The stored posterior mean, node by node ---------------------------------
  stored <- x$inla$summary.random[[paste0(term$row_var, "_id")]]
  values <- matrix(
    stored$mean[offset + seq_len(n_row * n_col)],
    nrow = n_row,
    ncol = n_col,
    dimnames = list(row_levels, col_levels)
  )
  list(
    values = values,
    row_var = term$row_var,
    col_var = term$col_var,
    trial = trial,
    title = if (is.null(trial)) {
      "Spatial field"
    } else {
      paste0("Spatial field, ", trial)
    },
    spelling = .ar1_term_spelling(term)
  )
}

# The viridis sequence at eleven steps, as grDevices::hcl.colors(11,
# "viridis") gives it, written out so the display needs no graphics-device
# package.
.FB_HEAT_COLOURS <- c(
  "#4B0055", "#44286E", "#274983", "#006791", "#008298", "#009B95",
  "#00B28A", "#00C475", "#7ED357", "#C2DE34", "#FDE333"
)

# .plot_spatial_field(): the type = "spatial_field" entry point.
#
# The palette is .FB_HEAT_COLOURS, ordered in lightness, so it reads in grey
# and to the common colour-vision deficiencies.
#
# @noRd
# @keywords internal
.plot_spatial_field <- function(x, trial = NULL, ...) {
  field <- .fb_spatial_field(x, trial)
  values <- field$values
  n_row <- nrow(values)
  n_col <- ncol(values)
  colours <- .FB_HEAT_COLOURS
  span <- range(values)
  breaks <- seq(span[[1L]], span[[2L]], length.out = length(colours) + 1L)

  old_par <- graphics::par(mar = c(5, 5, 4, 6))
  on.exit(graphics::par(old_par))
  graphics::image(
    x = seq_len(n_col),
    y = seq_len(n_row),
    z = t(values),
    col = colours,
    breaks = breaks,
    axes = FALSE,
    xlab = field$col_var,
    ylab = field$row_var,
    main = field$title,
    ...
  )
  graphics::axis(1L, at = seq_len(n_col), labels = colnames(values),
                 cex.axis = 0.8)
  graphics::axis(2L, at = seq_len(n_row), labels = rownames(values),
                 las = 1L, cex.axis = 0.8)
  graphics::box()
  graphics::mtext(
    paste0("posterior mean of ", field$spelling, ", response scale"),
    side = 3L, line = 0.4, cex = 0.8
  )

  # The colour key, in the right margin -------------------------------------
  usr <- graphics::par("usr")
  key_x <- usr[[2L]] + 0.04 * (usr[[2L]] - usr[[1L]])
  key_w <- 0.05 * (usr[[2L]] - usr[[1L]])
  key_y <- usr[[3L]] + (breaks - span[[1L]]) / diff(span) *
    (usr[[4L]] - usr[[3L]])
  graphics::rect(key_x, key_y[-length(key_y)], key_x + key_w, key_y[-1L],
                 col = colours, border = NA, xpd = NA)
  graphics::rect(key_x, key_y[[1L]], key_x + key_w, key_y[[length(key_y)]],
                 border = "black", xpd = NA)
  key_at <- c(span[[1L]], 0, span[[2L]])
  key_at <- key_at[key_at >= span[[1L]] & key_at <= span[[2L]]]
  graphics::text(
    key_x + key_w,
    usr[[3L]] + (key_at - span[[1L]]) / diff(span) * (usr[[4L]] - usr[[3L]]),
    labels = format(signif(key_at, 3L)), pos = 4L, cex = 0.7, xpd = NA
  )
  invisible(field)
}
