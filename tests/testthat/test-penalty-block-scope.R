# regularized_horseshoe(scope = "block"): the candidate block as its own brms
# non-linear part, y ~ a + b, with the protected effects in part a.
#
# The oracle for the fitted posterior is a brms fit written out by hand here,
# with the same formula, prior, seed and chains. Its Monte Carlo error is
# measured from three seeds of that direct fit: for each coefficient, the
# larger of the spread of its posterior mean over the seeds and the median
# spread over all coefficients as a share of the posterior standard deviation.
# Agreement is a gap from the mean of the three seeds within four standard
# errors of that difference.


# --- fixtures --------------------------------------------------------------- #

# 48 rows, 4 regions, 21 clusters nested in regions, 30 standardised
# candidates of which x01 and x02 carry real effects, and a known variance per
# row.
hsb_data <- function(seed = 20261001L) {
  set.seed(seed)
  n <- 48L
  n_cluster <- 21L
  cluster_region <- rep(c("R1", "R2", "R3", "R4"), length.out = n_cluster)
  cluster <- c(
    seq_len(n_cluster),
    sample(n_cluster, n - n_cluster, replace = TRUE)
  )
  d <- data.frame(
    cluster = factor(sprintf("C%02d", cluster)),
    region = factor(cluster_region[cluster])
  )
  for (j in seq_len(30L)) {
    d[[sprintf("x%02d", j)]] <- as.numeric(scale(stats::rnorm(n)))
  }
  u <- stats::rnorm(n_cluster, sd = 0.3)
  shift <- c(R1 = 0, R2 = -1, R3 = 0.5, R4 = 1)
  d$y <- 4 + shift[as.character(d$region)] + 0.8 * d$x01 - 0.6 * d$x02 +
    u[cluster] + stats::rnorm(n, sd = 0.3)
  d$v <- stats::runif(n, 0.01, 0.06)
  d
}

hsb_candidates <- sprintf("x%02d", seq_len(30L))

hsb_candidate_formula <- function() {
  stats::as.formula(paste("~", paste(hsb_candidates, collapse = " + ")))
}

hsb_penalty <- function() {
  fb_penalty(
    term(everything()) ~ regularized_horseshoe(p0 = 2, scope = "block")
  )
}

# par_ratio = p0 / (p - p0), calibrated from 48 rows to 21 clusters.
hsb_prior_string <- function() {
  sprintf(
    paste0(
      "horseshoe(df = 1, df_global = 1, scale_slab = 2, df_slab = 4, ",
      "par_ratio = %s)"
    ),
    format(2 / 28 * sqrt(48 / 21), digits = 15L)
  )
}


# --- grammar, targeting and compilation ------------------------------------ #

test_that("scope = 'block' is accepted and other values still refuse", {
  pen <- hsb_penalty()
  expect_identical(pen$specs[[1L]]$args$scope, "block")
  expect_error(
    fb_penalty(
      term(everything()) ~ regularized_horseshoe(p0 = 2, scope = "x")
    ),
    class = "flexybayes_refusal_penalty_hyperparameter_invalid"
  )
})

test_that("a block-scoped horseshoe resolves beside protected effects", {
  d <- hsb_data()
  map <- .fb_sel_term_map(~ region + x01 + x02 + x03, d)
  blocks <- .fb_pen_resolve(hsb_penalty(), map, protected = "region")
  expect_identical(blocks[[1L]]$coef_names, c("x01", "x02", "x03"))

  # The same penalty with the default scope is the refusal it replaces.
  expect_error(
    .fb_pen_resolve(
      fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)),
      map,
      protected = "region"
    ),
    class = "flexybayes_refusal_penalty_block_targeting_unsupported"
  )
})

test_that("a second block beside a block-scoped prior refuses", {
  d <- hsb_data()
  map <- .fb_sel_term_map(~ region + x01 + x02 + x03, d)
  expect_error(
    .fb_pen_resolve(
      fb_penalty(
        term(everything()) ~ regularized_horseshoe(p0 = 1, scope = "block"),
        term(x01) ~ ridge(scale = 1)
      ),
      map,
      protected = "region"
    ),
    class = "flexybayes_refusal_penalty_block_targeting_unsupported"
  )
})

test_that("the compiled prior row sits on part b and names the block", {
  skip_if_not_installed("brms")
  d <- hsb_data()
  map <- .fb_sel_term_map(~ region + x01 + x02 + x03, d)
  blocks <- .fb_pen_resolve(hsb_penalty(), map, protected = "region")
  blocks <- .fb_pen_scaling(blocks, map$X[, 5:7], train_rows = seq_len(48L))
  compiled <- .fb_pen_compile_brms(
    blocks,
    c("regionR2", "regionR3", "regionR4", "x01", "x02", "x03"),
    n_train = 48
  )
  expect_identical(compiled$prior$nlpar, "b")
  expect_identical(compiled$prior$class, "b")
  expect_match(compiled$prior$prior, "^horseshoe\\(")
  expect_identical(compiled$nl_block$protected_part, "a")
  expect_identical(compiled$nl_block$block_part, "b")
  expect_identical(compiled$nl_block$coef_names, c("x01", "x02", "x03"))
})

test_that("the linear formula splits into the protected part and the block", {
  skip_if_not_installed("brms")
  nl_block <- list(
    protected_part = "a",
    block_part = "b",
    candidate_terms = c("x01", "x02")
  )
  form <- .fb_nl_block_formula(
    y | se(s, sigma = TRUE) ~ region + x01 + x02 + (1 | cluster),
    nl_block
  )
  expect_true(isTRUE(attr(form$formula, "nl")))
  expect_identical(
    deparse(form$formula),
    "y | se(s, sigma = TRUE) ~ a + b"
  )
  expect_identical(deparse(form$pforms$a), "a ~ 1 + region + (1 | cluster)")
  expect_identical(deparse(form$pforms$b), "b ~ 0 + x01 + x02")
})

test_that("a variable named as a part refuses by name", {
  skip_if_not_installed("brms")
  nl_block <- list(
    protected_part = "a",
    block_part = "b",
    candidate_terms = "x"
  )
  expect_error(
    .fb_nl_block_formula(y ~ a + x, nl_block),
    class = "flexybayes_refusal_penalty_block_part_name_taken"
  )
})

test_that("protected prior rows move to part a and sigma stays", {
  skip_if_not_installed("brms")
  prior <- c(
    brms::set_prior("normal(0, 5)", class = "b", coef = "regionR2"),
    brms::set_prior("normal(4, 2)", class = "Intercept"),
    brms::set_prior("normal(0, 1)", class = "sd", group = "cluster"),
    brms::set_prior("normal(0, 2)", class = "sigma")
  )
  moved <- .fb_nl_block_prior(prior, list(protected_part = "a"))
  expect_identical(moved$class, c("b", "b", "sd", "sigma"))
  expect_identical(moved$coef, c("regionR2", "Intercept", "", ""))
  expect_identical(moved$nlpar, c("a", "a", "a", ""))
})

test_that("draw names lose the part prefix brms writes", {
  nm <- c(
    "b_a_Intercept", "b_a_regionR2", "b_b_x01", "b_b_a_soil", "sdb_b_x01",
    "hs_global_b", "sd_cluster__a_Intercept", "r_cluster__a[C01,Intercept]",
    "sigma"
  )
  expect_identical(
    .fb_nl_canonical_names(nm, c("a", "b")),
    c(
      "b_Intercept", "b_regionR2", "b_x01", "b_a_soil", "sdb_x01",
      "hs_global", "sd_cluster__Intercept", "r_cluster[C01,Intercept]",
      "sigma"
    )
  )
})


# --- fb_select() before any fit -------------------------------------------- #

hsb_inla_base <- function(d) {
  flexybayes(
    y ~ region + (1 | cluster),
    data = d,
    backend = "inla",
    verbose = FALSE
  )
}

test_that("the plan beside a protected region is a fit, decided in full", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.inla_threads = "1:1")
  d <- hsb_data()
  base <- hsb_inla_base(d)
  plan <- fb_select(
    base,
    candidates = hsb_candidate_formula(),
    penalty = hsb_penalty(),
    resamples = fb_resamples(unit = ~cluster, v = 3L),
    plan = TRUE,
    seed = 1L
  )
  expect_identical(plan$selection$decision, "full")
  expect_identical(plan$selection$decision_resolved_from, "auto")

  # Asked for by name, the projection refuses before anything is fitted.
  expect_error(
    fb_select(
      base,
      candidates = hsb_candidate_formula(),
      penalty = hsb_penalty(),
      resamples = fb_resamples(unit = ~cluster, v = 3L),
      decision = fb_decision(method = "projection"),
      plan = TRUE,
      seed = 1L
    ),
    class = "flexybayes_refusal_projection_structure_unsupported"
  )
})

test_that("a candidate outside the block refuses", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.inla_threads = "1:1")
  d <- hsb_data()
  base <- hsb_inla_base(d)
  expect_error(
    fb_select(
      base,
      candidates = ~ x01 + x02 + x03,
      penalty = fb_penalty(
        term(x01, x02) ~ regularized_horseshoe(p0 = 1, scope = "block")
      ),
      resamples = fb_resamples(unit = ~cluster, v = 3L),
      plan = TRUE,
      seed = 1L
    ),
    class = "flexybayes_refusal_penalty_block_targeting_unsupported"
  )
})


# --- the fitted reference against brms directly ---------------------------- #

hsb_cache <- new.env(parent = emptyenv())

hsb_engine <- list(
  chains = 2L,
  n_samples = 2000L,
  warmup = 1000L,
  adapt_delta = 0.995
)

# The reference fit and three seeds of the direct brms fit, on one dataset.
hsb_fits <- function(known_variance = FALSE) {
  key <- if (known_variance) "kv" else "plain"
  if (!is.null(hsb_cache[[key]])) {
    return(hsb_cache[[key]])
  }
  d <- hsb_data()
  base <- flexybayes(
    y ~ region + (1 | cluster),
    data = d,
    backend = "brms",
    known_variance = if (known_variance) d$v else NULL,
    chains = 2L,
    n_samples = 1000L,
    warmup = 1000L,
    seed = 11L,
    verbose = FALSE,
    mcmc_verbose = FALSE
  )
  sel <- fb_select(
    base,
    candidates = hsb_candidate_formula(),
    penalty = hsb_penalty(),
    resamples = fb_resamples(unit = ~cluster, v = 3L),
    engine = hsb_engine,
    seed = 20261001L
  )

  # The direct fit, written out. The base carries brms's defaults on the
  # coefficients and the package's half-normal on both standard deviations.
  scale <- 5 * stats::sd(d$y)
  lhs <- if (known_variance) "y | se(s, sigma = TRUE)" else "y"
  d$s <- sqrt(d$v)
  form <- brms::bf(
    stats::as.formula(paste(lhs, "~ a + b")),
    a ~ 1 + region + (1 | cluster),
    stats::as.formula(
      paste("b ~ 0 +", paste(hsb_candidates, collapse = " + "))
    ),
    nl = TRUE
  )
  prior <- c(
    brms::set_prior(hsb_prior_string(), class = "b", nlpar = "b"),
    brms::set_prior(
      sprintf("normal(0, %s)", format(scale, digits = 15L)),
      class = "sd", group = "cluster", nlpar = "a"
    ),
    brms::set_prior(
      sprintf("normal(0, %s)", format(scale, digits = 15L)),
      class = "sigma"
    )
  )
  direct <- list(brms::brm(
    form,
    data = d,
    prior = prior,
    chains = 2L,
    iter = 3000L,
    warmup = 1000L,
    seed = 20261001L,
    control = list(adapt_delta = 0.995),
    refresh = 0L,
    silent = 2L
  ))
  for (s in c(20261002L, 20261003L)) {
    direct[[length(direct) + 1L]] <- stats::update(
      direct[[1L]],
      seed = s,
      refresh = 0L,
      silent = 2L
    )
  }
  out <- list(data = d, base = base, sel = sel, direct = direct)
  assign(key, out, envir = hsb_cache)
  out
}

# Posterior means of the named coefficients, one column per direct seed.
hsb_direct_means <- function(direct, names) {
  vapply(
    direct,
    function(fit) {
      m <- colMeans(.fb_brms_draws_matrix(fit))
      m[paste0("b_", names)]
    },
    numeric(length(names))
  )
}

# The gap between the package's posterior means and the mean over the direct
# seeds, in Monte Carlo standard errors of that difference.
hsb_agreement <- function(fits, names) {
  means <- hsb_direct_means(fits$direct, names)
  post_sd <- apply(
    .fb_brms_draws_matrix(fits$direct[[1L]])[, paste0("b_", names)],
    2L,
    stats::sd
  )
  spread <- apply(means, 1L, stats::sd)
  share <- stats::median(spread / post_sd)
  mcse <- pmax(spread, share * post_sd)
  gap <- abs(stats::coef(fits$sel$reference)[names] - rowMeans(means))
  se_gap <- mcse * sqrt(1 + 1 / ncol(means))
  list(z = gap / se_gap, share = share, gap = gap, mcse = mcse)
}

test_that("block-scoped reference agrees with the direct brms fit", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_long_tier()
  fits <- hsb_fits()
  ref <- fits$sel$reference

  # Same formula, prior, data, seed and chains.
  expect_identical(
    deparse(ref$brms$formula$pforms$a),
    deparse(fits$direct[[1L]]$formula$pforms$a)
  )

  names_all <- c("regionR2", "regionR3", "regionR4", hsb_candidates)
  agreement <- hsb_agreement(fits, names_all)
  package_means <- stats::coef(ref)[names_all]

  # The candidates under their own names, within Monte Carlo error.
  expect_true(all(hsb_candidates %in% names(stats::coef(ref))))
  expect_false(any(grepl("^b_", names(stats::coef(ref)))))
  expect_true(all(agreement$z[hsb_candidates] < 4))

  # The real effects are found and the region effects are left unshrunk.
  expect_gt(package_means[["x01"]], 0.5)
  expect_lt(package_means[["x02"]], -0.4)
  region <- c("regionR2", "regionR3", "regionR4")
  expect_true(all(agreement$z[region] < 4))
  expect_lt(package_means[["regionR2"]], -0.6)
  expect_gt(package_means[["regionR4"]], 0.5)

  # tidy() and the posterior-threshold reader use the same names.
  tab <- tidy(fits$sel)
  expect_equal(
    tab$estimate[match(hsb_candidates, tab$term)],
    unname(package_means[hsb_candidates]),
    tolerance = 1e-10
  )
  newdata <- fits$data[1:3, ]
  expect_length(stats::predict(ref, newdata = newdata), 3L)
})

test_that("prior read-back shows the horseshoe on the block only", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_long_tier()
  ref <- hsb_fits()$sel$reference
  summary <- prior_summary(ref)
  table <- as.data.frame(summary$engine_prior_table)

  # One horseshoe row, on class b of part b, over exactly the candidates.
  hs <- grepl("^horseshoe\\(", table$prior)
  expect_identical(sum(hs), 1L)
  expect_identical(table$nlpar[hs], "b")
  expect_identical(table$prior[hs], hsb_prior_string())
  block_rows <- table$nlpar == "b" & nzchar(table$coef)
  expect_setequal(table$coef[block_rows], hsb_candidates)

  # Part a keeps the base's priors: no global-local prior on any row.
  protected <- table$nlpar == "a"
  expect_setequal(
    table$coef[protected & table$class == "b" & nzchar(table$coef)],
    c("Intercept", "regionR2", "regionR3", "regionR4")
  )
  expect_false(any(grepl("horseshoe", table$prior[protected])))
  expect_true(any(
    protected & table$class == "sd" & grepl("^normal\\(0, ", table$prior)
  ))
  expect_match(summary$fixed_sd_engine_default, "uncentred")

  # The manifest reads the block's row back as sent.
  slot <- summary$slots[summary$slots$engine_slot == "b|||||b", ]
  expect_identical(nrow(slot), 1L)
  expect_true(slot$readback_agrees)
})

test_that("block-scoped reference carries a known variance per row", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  skip_if_not_long_tier()
  fits <- hsb_fits(known_variance = TRUE)
  ref <- fits$sel$reference
  expect_match(deparse(ref$brms$formula$formula), "se\\(", fixed = FALSE)

  names_all <- c("regionR2", "regionR3", "regionR4", hsb_candidates)
  expect_true(all(hsb_agreement(fits, names_all)$z < 4))

  # Scored out of fold, the protected baseline refit cannot carry it.
  expect_error(
    fb_select(
      fits$base,
      candidates = hsb_candidate_formula(),
      penalty = hsb_penalty(),
      resamples = fb_resamples(unit = ~cluster, v = 3L),
      benchmark = TRUE,
      seed = 1L,
      plan = TRUE
    ),
    class = "flexybayes_refusal_known_variance_not_supported"
  )
})
