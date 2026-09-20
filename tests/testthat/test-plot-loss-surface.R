# =============================================================================
# plot(<fb_select>, type = "loss_surface"): the out-of-fold loss over the
# relative penalty s = lambda / lambda_max and the elastic-net mixing
# parameter alpha, drawn with plotly from optimiser results that were scored
# on one score contract, one per mixing value.
#
# The oracle is each result's stored per-unit out-of-fold table,
# `benchmark$glmnet$clusters_by_model$benchmark_<k>`: this file averages each
# scoring unit's rows and then the units with its own code, and compares every
# cell of the drawn surface with that value. The grid coordinates are read
# from the grid's definition and the mixing value off the support record; the
# chosen point is read off each result's own selection record.
# None of it calls the package's summarisers.
#
# What each block holds:
#
#   Values. Every cell of the drawn surface equals the unit-averaged stored
#   out-of-fold loss at its (alpha, s), in whichever order the results are
#   given, and the axes are the stored grid and the mixing values.
#
#   Chosen point. The highlighted marker is the result the plot was called
#   on, at its selected s and its alpha, at the stored loss there; every
#   other result's selected point is marked too.
#
#   Refusals. plotly absent, one mixing value, a repeated mixing value,
#   another partition of the units, another grid, another decision rule, and
#   an entry that is not an optimiser result on the s coordinate each refuse
#   by name before anything is drawn.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

# Sites of rows, one protected covariate and candidates z1..z<p>, three of
# which carry signal, with a site random effect.
fbls_data <- function(n_site = 24L, per = 5L, p = 8L, seed = 7L) {
  .fb_sel_with_seed(seed, {
    n <- n_site * per
    site <- factor(rep(sprintf("s%02d", seq_len(n_site)), each = per))
    d <- data.frame(site = site, x = stats::rnorm(n))
    for (j in seq_len(p)) {
      d[[paste0("z", j)]] <- stats::rnorm(n)
    }
    d$y <- 2 + 0.5 * d$x + 0.9 * d$z1 - 0.8 * d$z2 + 0.6 * d$z3 +
      stats::rnorm(n_site)[as.integer(site)] + stats::rnorm(n)
    d
  })
}

# A base carrying the parsed model, the data and the call record the fold
# routes read, with no engine behind it.
fbls_stub_base <- function(d) {
  fb <- suppressWarnings(fb_from_asreml(y ~ x, random = ~site, data = d))
  fb$family <- "gaussian"
  fb$link <- "identity"
  structure(
    list(
      fb = fb,
      data = d,
      extras = list(
        call_info = list(fixed = y ~ x, random = ~site, family = "gaussian")
      )
    ),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
}

fbls_quiet <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_high_missingness_warning = TRUE,
    flexyBayes.silence_default_prior_note = TRUE,
    .local_envir = envir
  )
}

fbls_select <- function(d, alpha = 1, seed = 1L, v = 6L,
                        decision = fb_decision(method = "full",
                                               metric = "rmse"),
                        penalty = NULL) {
  if (is.null(penalty)) {
    penalty <- fb_penalty(
      term(everything()) ~ elastic_net(lambda = 0.2, alpha = alpha)
    )
  }
  fb_select(
    fbls_stub_base(d),
    ~ z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8,
    penalty = penalty,
    workflow = "penalized_optimizer",
    resamples = fb_resamples(unit = ~site, v = v),
    decision = decision,
    seed = seed
  )
}

fbls_cache <- new.env(parent = emptyenv())

# Three results on one contract, at mixing values 0.25, 0.5 and 1.
fbls_trio <- function() {
  if (is.null(fbls_cache$trio)) {
    d <- fbls_data()
    fbls_cache$trio <- lapply(
      c(a025 = 0.25, a05 = 0.5, a1 = 1),
      function(alpha) fbls_select(d, alpha = alpha)
    )
  }
  fbls_cache$trio
}

# The drawn traces, as plotly writes them for the browser.
fbls_traces <- function(widget) {
  plotly::plotly_build(widget)$x$data
}

# The loss at every grid point of one result, from its stored per-unit
# out-of-fold tables: each unit's rows averaged, then the units.
fbls_hand_losses <- function(sel, metric = "rmse") {
  arm <- sel$benchmark$glmnet
  models <- grep("^benchmark_", names(arm$clusters_by_model), value = TRUE)
  index <- sort(as.integer(sub("^benchmark_", "", models)))
  vapply(index, function(k) {
    tab <- as.data.frame(arm$clusters_by_model[[paste0("benchmark_", k)]])
    per_unit <- vapply(
      split(tab[[metric]], as.character(tab$cluster)), mean, numeric(1L)
    )
    mean(per_unit)
  }, numeric(1L))
}

# The grid on s, as the optimiser defines it: data-free, geometric, twenty
# points from 1 to 1e-3.
fbls_hand_grid <- function() {
  exp(seq(0, log(1e-3), length.out = 20L))
}


# --- values -------------------------------------------------------------- #

test_that("every cell of the loss surface is the stored out-of-fold loss at its s and alpha", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("plotly")
  fbls_quiet()
  trio <- fbls_trio()

  widget <- plot(trio$a05, type = "loss_surface",
                 alongside = list(trio$a1, trio$a025))
  expect_s3_class(widget, "plotly")
  traces <- fbls_traces(widget)
  surface <- traces[[which(vapply(traces, function(t) identical(t$type, "surface"), logical(1L)))]]

  # The axes: the stored grid along x, the mixing values ascending along y.
  grid <- fbls_hand_grid()
  expect_identical(as.numeric(surface$x), grid)
  expect_identical(as.numeric(surface$y), c(0.25, 0.5, 1))

  # Each row of z is one mixing value, each column one grid point.
  z <- matrix(unlist(surface$z), nrow = 3L, byrow = !is.matrix(surface$z))
  if (is.matrix(surface$z)) {
    z <- surface$z
  }
  expect_identical(dim(z), c(3L, length(grid)))
  for (name in names(trio)) {
    row <- match(trio[[name]]$benchmark$glmnet$support$alpha, c(0.25, 0.5, 1))
    expect_equal(
      as.numeric(z[row, ]), fbls_hand_losses(trio[[name]]),
      tolerance = 1e-12, label = paste0("drawn losses at alpha of ", name)
    )
  }

  # The values the widget carries for a reader equal the drawn ones.
  values <- attr(widget, "values")
  expect_identical(sort(unique(values$surface$alpha)), c(0.25, 0.5, 1))
  for (name in names(trio)) {
    alpha <- trio[[name]]$benchmark$glmnet$support$alpha
    part <- values$surface[values$surface$alpha == alpha, , drop = FALSE]
    part <- part[order(part$lambda_index), , drop = FALSE]
    expect_identical(part$s, grid)
    expect_equal(part$loss, fbls_hand_losses(trio[[name]]), tolerance = 1e-12)
  }
  expect_identical(values$metric, "rmse")
  expect_identical(
    values$contract_hash,
    unique(trio$a05$assessment_predictions$contract_hash)
  )

  # The order the results are given in does not change the surface.
  again <- plot(trio$a05, type = "loss_surface",
                alongside = list(trio$a025, trio$a1))
  expect_identical(attr(again, "values")$surface, values$surface)
})


# --- chosen point -------------------------------------------------------- #

test_that("the highlighted marker is the selected point of the result plotted", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("plotly")
  fbls_quiet()
  trio <- fbls_trio()

  for (primary in names(trio)) {
    others <- unname(trio[setdiff(names(trio), primary)])
    widget <- plot(trio[[primary]], type = "loss_surface", alongside = others)
    chosen <- attr(widget, "values")$chosen
    selection <- trio[[primary]]$benchmark$glmnet$selection
    losses <- fbls_hand_losses(trio[[primary]])

    # One row per result; exactly one is the plotted result's.
    expect_identical(nrow(chosen), 3L)
    expect_identical(sum(chosen$primary), 1L)
    mine <- chosen[chosen$primary, , drop = FALSE]
    expect_identical(mine$alpha, trio[[primary]]$benchmark$glmnet$support$alpha)
    expect_identical(mine$s, selection$selected_s)
    expect_identical(mine$lambda_index, as.integer(selection$selected_lambda_index))
    expect_equal(mine$loss, losses[[selection$selected_lambda_index]],
                 tolerance = 1e-12)

    # The drawn highlight sits there.
    traces <- fbls_traces(widget)
    marker_traces <- traces[vapply(traces, function(t) identical(t$type, "scatter3d"), logical(1L))]
    highlight <- Filter(function(t) identical(t$name, "chosen point"), marker_traces)
    expect_length(highlight, 1L)
    expect_identical(as.numeric(highlight[[1L]]$x), selection$selected_s)
    expect_identical(as.numeric(highlight[[1L]]$y), mine$alpha)
    expect_equal(as.numeric(highlight[[1L]]$z), mine$loss, tolerance = 1e-12)

    # Every other result's own selected point is marked at its own alpha.
    for (other in others) {
      row <- chosen[chosen$alpha == other$benchmark$glmnet$support$alpha, , drop = FALSE]
      expect_identical(row$s, other$benchmark$glmnet$selection$selected_s)
      expect_false(row$primary)
    }
  }
})


# --- refusals ------------------------------------------------------------ #

test_that("the loss surface refuses by name without plotly", {
  skip_if_not_installed("glmnet")
  fbls_quiet()
  trio <- fbls_trio()
  local_mocked_bindings(.fb_plot_plotly_installed = function() FALSE)
  expect_error(
    plot(trio$a05, type = "loss_surface", alongside = list(trio$a1)),
    class = "flexybayes_refusal_loss_surface_requires_plotly"
  )
  expect_error(
    plot(trio$a05, type = "loss_surface", alongside = list(trio$a1)),
    regexp = "plotly", fixed = TRUE
  )
})

test_that("results that do not span a surface on one contract refuse by name", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("plotly")
  fbls_quiet()
  trio <- fbls_trio()
  d <- fbls_data()
  code <- "flexybayes_refusal_selection_loss_surface_undefined"

  # One mixing value, alone or repeated.
  expect_error(plot(trio$a05, type = "loss_surface"), class = code)
  expect_error(
    plot(trio$a05, type = "loss_surface", alongside = list(fbls_select(d, alpha = 0.5))),
    class = code
  )
  # Another assignment of the sites to folds is another score contract.
  other_folds <- fbls_select(d, alpha = 1, seed = 2L)
  expect_false(identical(
    unique(other_folds$assessment_predictions$contract_hash),
    unique(trio$a05$assessment_predictions$contract_hash)
  ))
  expect_error(
    plot(trio$a05, type = "loss_surface", alongside = list(other_folds)),
    class = code
  )
  # Another grid on s. `fb_decision()` fixes the grid, so a result from
  # another grid is made by dropping the weakest point of a stored one.
  other_grid <- trio$a1
  scores <- other_grid$benchmark$glmnet$selection$scores
  other_grid$benchmark$glmnet$selection$scores <-
    scores[scores$lambda_index < max(scores$lambda_index), ]
  expect_error(
    plot(trio$a05, type = "loss_surface", alongside = list(other_grid)),
    class = code
  )
  # Another decision rule chooses its point by another rule.
  other_rule <- fbls_select(d, alpha = 1, decision = fb_decision(
    method = "full", metric = "rmse", rule = "min"
  ))
  expect_error(
    plot(trio$a05, type = "loss_surface", alongside = list(other_rule)),
    class = code
  )
  # An entry that is not a selection, and a ridge path, which has no s.
  expect_error(
    plot(trio$a05, type = "loss_surface", alongside = list(data.frame(x = 1))),
    class = code
  )
  ridge <- fbls_select(d, penalty = fb_penalty(term(everything()) ~ ridge(scale = 1)))
  expect_error(
    plot(trio$a05, type = "loss_surface", alongside = list(ridge)),
    class = code
  )
  # `alongside` is a list of results.
  expect_error(
    plot(trio$a05, type = "loss_surface", alongside = trio$a1),
    class = code
  )
})

test_that("a result with no path on s says so instead of drawing", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("plotly")
  fbls_quiet()
  d <- fbls_data()
  ridge <- fbls_select(d, penalty = fb_penalty(term(everything()) ~ ridge(scale = 1)))
  expect_message(
    out <- plot(ridge, type = "loss_surface", alongside = list(fbls_trio()$a1)),
    "no s = lambda / lambda_max coordinate", fixed = TRUE
  )
  expect_null(out)
})
