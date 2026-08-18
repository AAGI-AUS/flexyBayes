# Plot methods for flexyBayes

#' Plot diagnostics for a flexyBayes model
#'
#' Draws one of seven standard displays for a fitted model, chosen by
#' `type`. Every backend reaches the same method, and a display that needs
#' a slot a given engine does not populate -- posterior draws on the
#' deterministic INLA path, for instance -- prints a message naming what
#' is missing and returns invisibly rather than erroring or drawing an
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
#' @param x A fitted `flexybayes` object of any backend. Aggregated,
#'   generalised-linear and direct-greta fits reach the same method
#'   through their own registrations.
#' @param type A single string naming the display to draw. One of:
#'   `"diagnostics"`, trace plots and marginal densities per parameter
#'   (needs a sampled posterior); `"residuals"`, residuals against
#'   fitted values beside a normal quantile-quantile plot; `"effects"`, a
#'   forest plot of the fixed effects with their credible intervals,
#'   available on every backend that supplies [coef()] and [confint()];
#'   `"variance"`, a bar chart of the variance components with credible
#'   intervals; `"blups"`, a caterpillar plot of the random-effect
#'   predictions ordered by magnitude; `"pp_check"`, a posterior
#'   predictive check overlaying replicated datasets on the observed
#'   response, which needs predictive draws and so refuses by name on a
#'   fit that carries none (see [pp_check.flexybayes()]); and
#'   `"variogram"`, the empirical semivariance of the residuals over the
#'   design index. `NULL` (the default) resolves to `"diagnostics"` on a
#'   sampled fit and `"residuals"` otherwise.
#' @param ... Further arguments passed to the underlying plotting call,
#'   for example `variable` to restrict which parameters a diagnostic
#'   display covers, or `type` and `ndraws` for `"pp_check"`.
#' @returns Invisibly, the object the underlying plotting call returns --
#'   a \pkg{ggplot2} object for the displays built with it, the
#'   semivariance table for `"variogram"`, and `NULL` for those drawn on
#'   the base graphics device. Called for the plot it draws.
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
    "variogram" = .plot_variogram(x, ...)
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
  "variogram"
)

# .fb_has_sampler_draws() --- did an engine draw samples for this fit?
#
# Written fresh rather than lifted from .plot_diagnostics(), which tested
# `x$greta$draws` alone: that slot is NULL on a brms fit, so plot() on
# one declined to draw its diagnostics and named brms as a supported
# backend in the same sentence. A brms fit keeps its draws inside the
# brmsfit at `$brms`, and a nested Laplace approximation has none at all.
#
# @noRd
# @keywords internal
.fb_has_sampler_draws <- function(x) {
  if (!is.null(x$greta$draws)) {
    return(TRUE)
  }
  inherits(x$brms, "brmsfit")
}

# .fb_default_plot_type() --- the display a fit earns by what it holds.
#
# @noRd
# @keywords internal
.fb_default_plot_type <- function(x) {
  if (.fb_has_sampler_draws(x)) {
    return("diagnostics")
  }
  "residuals"
}

# Registered for each fit class rather than relying on inheritance. The
# brms and INLA fits do share the "flexybayes" parent as of 0.9.0, so the
# parent registration would reach them, but the aggregated, glm and
# direct-greta classes carry their own dispatch order and an explicit
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
plot.flexybayes_direct_greta <- function(x, ...) plot.flexybayes(x, ...)

#' @rdname plot.flexybayes
#' @export
plot.flexybayes_glm <- function(x, ...) plot.flexybayes(x, ...)

# Emit a non-silent, non-erroring notice that a plot type is not
# available for this fit's backend, then return invisibly. Backends
# differ in which slots they expose (greta carries MCMC draws, fitted
# values and a variance-component table; INLA does not), so a plot type
# that reads a slot the backend never populates degrades to a message
# rather than crashing through to graphics::plot.default().
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
        "status and a Kullback-Leibler divergence instead -- see ",
        "summary(fit)$converge -- and plot(fit) defaults to the residual ",
        "panels on such a fit."
      )
    ))
  }

  # A brms fit keeps its draws inside the brmsfit, and brms draws its own
  # panels from them. Forwarding is what makes plot(brms_fit) show the
  # trace and density it always claimed to.
  if (is.null(x$greta$draws)) {
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
    return(invisible(do.call(
      graphics::plot,
      c(list(x = x$brms), dots)
    )))
  }

  if (requireNamespace("bayesplot", quietly = TRUE)) {
    draws <- x$greta$draws
    p1 <- bayesplot::mcmc_trace(draws, ...)
    print(p1)
  } else {
    # Base R fallback
    draws <- x$greta$draws
    all_draws <- do.call(rbind, lapply(draws, as.matrix))
    n_params <- min(ncol(all_draws), 6)
    par_names <- colnames(all_draws)[seq_len(n_params)]

    old_par <- par(mfrow = c(n_params, 2), mar = c(3, 3, 2, 1))
    on.exit(par(old_par))

    for (i in seq_len(n_params)) {
      nm <- par_names[i]
      # Trace
      for (ch in seq_along(draws)) {
        vals <- as.matrix(draws[[ch]])[, nm]
        if (ch == 1) {
          plot(
            vals,
            type = "l",
            main = paste("Trace:", nm),
            xlab = "",
            ylab = "",
            col = ch
          )
        } else {
          lines(vals, col = ch)
        }
      }
      # Density
      vals <- all_draws[, nm]
      plot(density(vals), main = paste("Density:", nm), xlab = "", ylab = "")
    }
  }
  invisible(NULL)
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
  # field as a list of posterior means rather than the five-column table
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
  # carries no intervals still gets an axis, rather than -Inf.
  x_upper <- suppressWarnings(
    max(c(vc$q97.5, vc$estimate), na.rm = TRUE)
  )
  if (!is.finite(x_upper) || x_upper <= 0) {
    x_upper <- 1
  }
  xlim <- c(0, x_upper * 1.1)

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
# draws to replicate from -- which is the brms engine. The seam is
# pp_check.flexybayes() (R/pp_check_support.R), so the plot() entry
# point and the bayesplot generic cannot disagree about what a check
# is, and a fit with no predictive draws refuses by name here too.
#
# Before 0.9.1 this drew observed values against fitted values, plus a
# density of the fitted values labelled Posterior Predictive, on every
# backend that carried a response and a fitted vector. No replicated
# dataset was ever involved. The panel is gone rather than retitled:
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
# design index. It is an EMPIRICAL variogram of the residuals -- no
# fitted variogram is overlaid and none is claimed, so a flat surface
# here says the fitted structure absorbed what was there and nothing
# stronger.
#
# Observed rows only. A row carried as a latent design cell has no
# observed response and therefore no residual, and pairing an NA into a
# squared difference would silently drop whole lags rather than one
# pair. The subtitle prints both counts so the reader can see how much
# of the array the picture rests on.

# .fb_variogram_index() --- the design index a fit is built over.
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

# .fb_variogram_position() --- one index column as an ordinal position.
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

# .fb_variogram_table() --- semivariance by lag.
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

# .fb_variogram_draw() --- the picture, on the base device this file uses.
#
# The bottom margin is five lines, not the four this file's other panels
# use, because `sub =` is drawn at line `par("mgp")[1] + 1` of the bottom
# margin -- line 4 with the default mgp -- and a four-line margin does not
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

# .plot_variogram() --- the type = "variogram" entry point.
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
        "Fit a spatial term -- random = ~ ar1(row):ar1(col) -- or use ",
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
