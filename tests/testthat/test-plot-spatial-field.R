# =============================================================================
# plot(<flexybayes fit>, type = "spatial_field"): the posterior mean of a
# separable autoregressive field, ar1(row):ar1(col), on the row-by-column
# array of one trial, labelled with the trial identifier the caller supplies
# and the array's own row and column labels.
#
# Oracles, all read from the engine's own output and none from the package's
# index code:
#
#   The plot-to-node map. INLA's posterior samples hold, for every plot i,
#   the linear predictor, the fixed effects and the field node by node. In
#   each sample the predictor equals the fixed part plus the field at the
#   node the plot was indexed to, exactly, so the node of every plot is read
#   off the engine's own draws, using the index columns INLA received
#   (`$inla$.args$data`) and INLA's documented grouped ordering
#   (node = row index + n_row (column index - 1) + n_row n_col (trial - 1)).
#
#   The values. The drawn cell of plot i equals the stored posterior mean of
#   the field at that node, `$inla$summary.random`.
#
#   The scale. The stored mean is mapped back from the standardised response
#   INLA was run on. On the engine's own coordinates the posterior mean of
#   the field is the weight-averaged conditional mean over INLA's
#   hyperparameter configurations; times the response standard deviation it
#   must equal the drawn field.
# =============================================================================

# --- fixtures ------------------------------------------------------------ #

.fbsf_ar1 <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))

# One observation per node of an n_row x n_col array, a separable AR1 field
# plus a nugget; rows and columns labelled with the caller's own strings.
# The plots are in shuffled order and the column factor's levels run from
# the last column to the first, so neither the order of the rows nor a sort
# of the labels gives the index the emitter wrote.
.fbsf_grid <- function(n_row = 8L, n_col = 6L, seed = 2026L, shift = 0) {
  set.seed(seed)
  g <- expand.grid(col = seq_len(n_col), row = seq_len(n_row))
  chol_l <- t(chol(1.44 * kronecker(.fbsf_ar1(n_row, 0.6), .fbsf_ar1(n_col, 0.3))))
  g$y <- 20 + shift + as.numeric(chol_l %*% stats::rnorm(n_row * n_col)) +
    stats::rnorm(n_row * n_col, 0, 0.4)
  g$row <- factor(sprintf("R%02d", g$row))
  g$range <- factor(sprintf("C%02d", g$col), levels = sprintf("C%02d", n_col:1))
  g$col <- NULL
  g[sample.int(nrow(g)), , drop = FALSE]
}

.fbsf_silence <- function(envir = parent.frame()) {
  withr::local_options(
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_uniform_inla_approx = TRUE,
    .local_envir = envir
  )
}

.fbsf_null_device <- function(envir = parent.frame()) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  withr::defer(
    {
      grDevices::dev.off()
      unlink(path)
    },
    envir = envir
  )
  invisible(path)
}

.fbsf_cache <- new.env(parent = emptyenv())

.fbsf_field_fit <- function() {
  if (is.null(.fbsf_cache$field)) {
    .fbsf_cache$field <- suppressMessages(flexybayes(
      y ~ 1,
      random = ~ ar1(row):ar1(range),
      data = .fbsf_grid(),
      backend = "inla",
      na_action = "augment",
      verbose = FALSE
    ))
  }
  .fbsf_cache$field
}

# Two trials of 6 x 5, one field per trial.
.fbsf_at_fit <- function() {
  if (is.null(.fbsf_cache$at)) {
    g1 <- .fbsf_grid(6L, 5L, seed = 11L)
    g1$trial <- "T1"
    g2 <- .fbsf_grid(6L, 5L, seed = 12L, shift = 3)
    g2$trial <- "T2"
    g <- rbind(g1, g2)
    g$trial <- factor(g$trial)
    .fbsf_cache$at <- suppressMessages(flexybayes(
      y ~ trial,
      random = ~ at(trial):ar1(row):ar1(range),
      data = g,
      backend = "inla",
      na_action = "augment",
      verbose = FALSE
    ))
  }
  .fbsf_cache$at
}

# The node of every plot, from the index columns INLA received, confirmed on
# the engine's own posterior draws: in each draw, the predictor of plot i is
# the fixed part plus the field at node i. Returns the nodes.
.fbsf_engine_nodes <- function(fit, row_var, col_var, fixed = ~1,
                               at_var = NULL, n_draws = 100L) {
  engine_data <- fit$inla$.args$data
  row_id <- engine_data[[paste0(row_var, "_id")]]
  col_id <- engine_data[[paste0(col_var, "_id")]]
  n_row <- max(row_id)
  n_col <- max(col_id)
  node <- row_id + n_row * (col_id - 1L)
  if (!is.null(at_var)) {
    node <- node + n_row * n_col * (engine_data[[paste0(at_var, "_id")]] - 1L)
  }
  draws <- INLA::inla.posterior.sample(
    n_draws, fit$inla, num.threads = "1:1", seed = 1L
  )
  latent_names <- rownames(draws[[1L]]$latent)
  predictor <- grep("^Predictor:", latent_names)
  field <- grep(paste0("^", row_var, "_id:"), latent_names)
  # INLA keeps its data as a list; model.matrix() of an intercept-only
  # formula over a list would have no rows.
  X <- stats::model.matrix(fixed, as.data.frame(engine_data))
  gap <- vapply(draws, function(s) {
    lat <- s$latent[, 1L]
    beta <- lat[paste0(colnames(X), ":1")]
    difference <- lat[predictor] - drop(X %*% beta) - lat[field][node]
    # Every plot is compared, or the gap is not a check.
    if (length(difference) != length(node) || length(node) == 0L) {
      return(Inf)
    }
    max(abs(difference))
  }, numeric(1L))
  list(node = node, max_gap = max(gap), n_nodes = length(field))
}


# --- the values ------------------------------------------------------------ #

test_that("each drawn cell is the stored field mean at the node the engine indexed that plot to", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .fbsf_silence()
  .fbsf_null_device()
  fit <- .fbsf_field_fit()
  d <- fit$data

  out <- plot(fit, type = "spatial_field", trial = "trial 07")
  expect_type(out, "list")
  field <- out$values
  expect_identical(dim(field), c(8L, 6L))
  # Rows and columns in the order of the caller's factor levels.
  expect_identical(rownames(field), sprintf("R%02d", 1:8))
  expect_identical(colnames(field), sprintf("C%02d", 6:1))

  # The map from plot to node, read off the engine's draws: exact up to
  # floating-point rounding in every draw.
  engine <- .fbsf_engine_nodes(fit, "row", "range")
  expect_lt(engine$max_gap, 1e-10)
  expect_identical(engine$n_nodes, 48L)

  stored <- fit$inla$summary.random$row_id$mean
  drawn <- field[cbind(as.character(d$row), as.character(d$range))]
  expect_identical(drawn, stored[engine$node])
  # Every node is drawn once.
  expect_identical(sort(engine$node), seq_len(48L))
})

test_that("the drawn field is the posterior mean on the response scale, read against the engine's own coordinates", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .fbsf_silence()
  .fbsf_null_device()
  fit <- .fbsf_field_fit()
  d <- fit$data
  out <- plot(fit, type = "spatial_field")
  engine <- .fbsf_engine_nodes(fit, "row", "range", n_draws = 10L)

  # INLA was run on (y - mean(y)) / sd(y); read that scale off the data.
  scale <- stats::sd(d$y)
  configs <- fit$inla$misc$configs
  tag <- which(configs$contents$tag == "row_id")
  # A configuration holds the latent field without the linear predictor that
  # `contents` lists first, when INLA pre-optimises the model.
  offset <- sum(configs$contents$length) - length(configs$config[[1L]]$mean)
  index <- configs$contents$start[[tag]] - offset - 1L +
    seq_len(configs$contents$length[[tag]])
  log_post <- vapply(configs$config, function(k) k$log.posterior, numeric(1L))
  weight <- exp(log_post - max(log_post))
  weight <- weight / sum(weight)
  engine_mean <- Reduce(`+`, Map(
    function(k, w) w * k$improved.mean[index], configs$config, weight
  ))
  drawn <- out$values[cbind(as.character(d$row), as.character(d$range))]
  # The stored means come from INLA's integrated marginals, the weighted
  # conditional means from its configurations. Measured on this fixture, in
  # units of the field's smallest posterior sd: the two differ by 4.4e-4; a
  # field left on the engine's coordinates would sit 0.54 away, and a
  # transposed map 5.1 away.
  sd_min <- min(fit$inla$summary.random$row_id$sd)
  expect_lt(max(abs(drawn - scale * engine_mean[engine$node])), 1e-2 * sd_min)
})

test_that("one field per trial draws the trial asked for", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .fbsf_silence()
  .fbsf_null_device()
  fit <- .fbsf_at_fit()
  d <- fit$data
  engine <- .fbsf_engine_nodes(fit, "row", "range", fixed = ~trial,
                               at_var = "trial")
  expect_lt(engine$max_gap, 1e-10)
  stored <- fit$inla$summary.random$row_id$mean

  for (trial in c("T1", "T2")) {
    out <- plot(fit, type = "spatial_field", trial = trial)
    rows <- which(as.character(d$trial) == trial)
    drawn <- out$values[cbind(as.character(d$row[rows]), as.character(d$range[rows]))]
    expect_identical(drawn, stored[engine$node[rows]])
    expect_identical(dim(out$values), c(6L, 5L))
  }
  # The field of a fit with one field per trial needs the trial named.
  expect_error(plot(fit, type = "spatial_field"), "`trial`", fixed = TRUE)
  expect_error(plot(fit, type = "spatial_field", trial = "T3"), "`trial`", fixed = TRUE)
})

test_that("the display carries only the identifier the caller supplied and the array's labels", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .fbsf_silence()
  .fbsf_null_device()
  fit <- .fbsf_field_fit()

  out <- plot(fit, type = "spatial_field", trial = "trial 07")
  expect_identical(out$title, "Spatial field, trial 07")
  expect_identical(out$row_var, "row")
  expect_identical(out$col_var, "range")
  bare <- plot(fit, type = "spatial_field")
  expect_identical(bare$title, "Spatial field")
  # The drawn values do not depend on the label.
  expect_identical(bare$values, out$values)
  expect_error(plot(fit, type = "spatial_field", trial = c("a", "b")), "`trial`", fixed = TRUE)
})


# --- refusals ------------------------------------------------------------- #

test_that("a fit with no separable field refuses by name", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  .fbsf_silence()
  set.seed(88L)
  d <- data.frame(g = factor(rep(letters[1:6], each = 8L)), x = stats::rnorm(48L))
  d$y <- 2 + 0.5 * d$x + stats::rnorm(48L)
  fit <- suppressMessages(flexybayes(
    y ~ x, random = ~g, data = d, backend = "inla", aggregate = FALSE,
    verbose = FALSE
  ))
  expect_error(
    plot(fit, type = "spatial_field"),
    class = "flexybayes_refusal_spatial_field_requires_ar1_field"
  )
  expect_error(plot(fit, type = "spatial_field"), "ar1(row):ar1(col)", fixed = TRUE)
})

test_that("a fit carrying two separable fields asks which one", {
  terms <- list(
    list(type = "ar1_spatial", row_var = "row", col_var = "col", col_ar1 = TRUE),
    list(type = "ar1_spatial", row_var = "row2", col_var = "col2", col_ar1 = TRUE)
  )
  stub <- structure(
    list(fb = list(random_terms = terms), data = data.frame(row = 1)),
    class = c("flexybayes_inla", "flexybayes", "list")
  )
  expect_error(plot(stub, type = "spatial_field"), "2 separable fields", fixed = TRUE)
})
