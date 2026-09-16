# The default variance-component prior is a half-normal on both engines.
#
# Every variance component the package default reaches carries
# half_normal(scale = s) on its standard deviation, on INLA and on brms, with
# s = 5 sd(y) / d on an identity link, 3 / d on a log link and 5 / d on a logit
# link, and s / sd(x) on an uncorrelated slope on a numeric covariate x. d is
# the divisor adopted for "matched to the former bound". Before 0.13.0 the
# default was uniform(0, U) with U = 5 sd(y), 3 or 5; INLA's reported interval
# for a weakly identified component could run past U, which the support audit
# refuses (agridat::john.alpha, `rep + rep:block`, 8 of 8 repeats).
#
# The checks read what each engine received. On INLA: INLA's own record of its
# priors (`all.hyper`, through the manifest's read-back), with each expression
# evaluated against the closed-form half-normal density on the log precision.
# On brms: the Stan program and brms::prior_summary(). Every expected scale is
# computed here from the data; none is read from the package.

.hn1_divisor <- 1

# --- fixtures ------------------------------------------------------------ #

# The capability grid's fixtures (tools/execution_grid.R), in the same
# generation order, so each class runs on the design its grid cell runs on.
.hn1_env <- new.env(parent = emptyenv())

.hn1_data <- function() {
  if (!is.null(.hn1_env$data)) {
    return(.hn1_env$data)
  }
  set.seed(11L)
  d_cap <- expand.grid(
    rep = factor(seq_len(2L)),
    gen = factor(seq_len(6L)),
    env = factor(seq_len(3L))
  )
  d_cap$x <- stats::rnorm(nrow(d_cap))
  d_cap$y <- stats::rnorm(nrow(d_cap), 10, 2)
  d_cap$bin <- stats::rbinom(nrow(d_cap), 1L, 0.5)
  set.seed(21L)
  d_cap$hg <- ifelse(
    stats::rbinom(nrow(d_cap), 1L, 0.7) == 1L,
    stats::rgamma(nrow(d_cap), shape = 2, rate = 0.5),
    0
  )
  w_cap <- as.numeric(seq_len(nrow(d_cap)))

  set.seed(12L)
  d_capgrid <- expand.grid(row = factor(seq_len(6L)), col = factor(seq_len(5L)))
  d_capgrid$y <- stats::rnorm(nrow(d_capgrid))

  k_cap <- diag(6)
  dimnames(k_cap) <- list(as.character(seq_len(6L)), as.character(seq_len(6L)))
  q_cap <- Matrix::Matrix(k_cap, sparse = TRUE)
  dimnames(q_cap) <- dimnames(k_cap)

  set.seed(13L)
  d_agg <- expand.grid(
    rep = factor(seq_len(20L)),
    gen = factor(seq_len(6L)),
    env = factor(seq_len(3L))
  )
  d_agg$y <- stats::rnorm(nrow(d_agg), 10, 2)
  d_agg$cnt <- stats::rpois(nrow(d_agg), 3)

  ar1_mat <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))
  d_ar1_at_trial <- do.call(rbind, lapply(seq_len(3L), function(i) {
    set.seed(1000L + i)
    grid <- expand.grid(col = seq_len(8L), row = seq_len(6L))
    lc <- t(chol(kronecker(ar1_mat(6L, 0.5), ar1_mat(8L, 0.3))))
    u <- as.numeric(lc %*% stats::rnorm(48L))
    data.frame(
      y = 5 + u + stats::rnorm(48L, 0, 0.5),
      row = factor(grid$row), col = factor(grid$col),
      trial = factor(paste0("T", i))
    )
  }))

  set.seed(2001L)
  grid_tg <- expand.grid(t = seq_len(12L), g = seq_len(10L))
  lt <- t(chol(ar1_mat(12L, 0.6)))
  u_tg <- as.numeric(vapply(
    seq_len(10L),
    function(i) as.numeric(lt %*% stats::rnorm(12L)),
    numeric(12L)
  ))
  d_ar1_tg <- data.frame(
    y = 5 + u_tg + stats::rnorm(120L, 0, 0.5),
    t = factor(grid_tg$t), g = factor(grid_tg$g)
  )
  d_ar1_tg$y[c(7L, 50L)] <- NA

  .hn1_env$data <- list(
    d_cap = d_cap, w_cap = w_cap, d_capgrid = d_capgrid, k_cap = k_cap,
    q_cap = q_cap, d_agg = d_agg, d_ar1_at_trial = d_ar1_at_trial,
    d_ar1_tg = d_ar1_tg
  )
  .hn1_env$data
}

# The default's scale in original units, by link.
.hn1_scale <- function(link, y = NULL) {
  switch(
    link,
    identity = 5 * stats::sd(y, na.rm = TRUE) / .hn1_divisor,
    log = 3 / .hn1_divisor,
    logit = 5 / .hn1_divisor
  )
}

# One INLA expression prior evaluated at `theta`, as R arithmetic.
.hn1_eval_expression <- function(text, theta) {
  body <- sub("return\\(", "(", sub("^expression:", "", text))
  env <- new.env()
  env$theta <- theta
  eval(parse(text = body), envir = env)
}

# The half-normal density on the SD, written on INLA's log precision:
# sigma = exp(-theta / 2), |d sigma / d theta| = sigma / 2.
.hn1_closed_form <- function(theta, s) {
  sigma <- exp(-theta / 2)
  log(2) + stats::dnorm(sigma, 0, s, log = TRUE) + log(sigma / 2)
}

# --- the model classes ----------------------------------------------------- #

# Each capability-matrix row that fits on INLA: the call, the link, the
# response the scale is read from, the default slots (manifest input
# identifiers) and the slope slot with its covariate.
.hn1_inla_classes <- function(f) {
  fb_inla <- function(...) {
    suppressMessages(flexybayes(..., backend = "inla", verbose = FALSE))
  }
  list(
    "Gaussian LMM, simple random intercept" = list(
      fit = function() fb_inla(y ~ env, random = ~ gen, data = f$d_cap),
      link = "identity", y = f$d_cap$y, slots = c("sigma:sigma", "sd:gen")
    ),
    "GLMM (binomial, Poisson, negative binomial, gamma, beta), simple random effect" = list(
      fit = function() {
        fb_inla(bin ~ env, random = ~ gen, data = f$d_cap, family = "binomial")
      },
      link = "logit", y = NULL, slots = "sd:gen"
    ),
    "Uncorrelated random slope" = list(
      fit = function() fb_inla(y ~ x + (x || gen), data = f$d_cap),
      link = "identity", y = f$d_cap$y,
      slots = c("sigma:sigma", "sd:gen", "sd:x_gen"),
      slope = list(id = "sd:x_gen", x = f$d_cap$x)
    ),
    "Factor-by-numeric fixed interaction" = list(
      fit = function() fb_inla(y ~ env * x, random = ~ gen, data = f$d_cap),
      link = "identity", y = f$d_cap$y, slots = c("sigma:sigma", "sd:gen")
    ),
    "Nested / interaction random effects, multi-stratum" = list(
      fit = function() fb_inla(y ~ env, random = ~ gen:env, data = f$d_cap),
      link = "identity", y = f$d_cap$y, slots = c("sigma:sigma", "sd:gen:env")
    ),
    "Heterogeneous variance by factor level" = list(
      fit = function() fb_inla(y ~ env, random = ~ diag(env):gen, data = f$d_cap),
      link = "identity", y = f$d_cap$y, slots = c("sigma:sigma", "sd:gen")
    ),
    "Heterogeneous residual by factor level" = list(
      fit = function() {
        fb_inla(y ~ env, random = ~ gen, residual = ~ dsum(~ units | env),
                data = f$d_cap)
      },
      link = "identity", y = f$d_cap$y, slots = c("sigma:sigma", "sd:gen")
    ),
    "Baseline heteroscedastic MET: iid genotype, iid genotype-by-environment and a sectioned residual" = list(
      fit = function() {
        fb_inla(y ~ env, random = ~ gen + gen:env,
                residual = ~ dsum(~ units | env), data = f$d_cap)
      },
      link = "identity", y = f$d_cap$y,
      slots = c("sigma:sigma", "sd:gen", "sd:gen:env")
    ),
    "Factor-analytic genotype-by-environment covariance" = list(
      fit = function() fb_inla(y ~ env, random = ~ fa(env, 1):gen, data = f$d_cap),
      link = "identity", y = f$d_cap$y, slots = c("sigma:sigma", "sd:gen")
    ),
    "Known-covariance genomic / pedigree random effect" = list(
      fit = function() {
        fb_inla(y ~ env, random = ~ vm(gen, cov = fb_cov(Q, type = "precision")),
                data = f$d_cap, known_matrices = list(Q = f$q_cap))
      },
      link = "identity", y = f$d_cap$y, slots = c("sigma:sigma", "sd:gen")
    ),
    "Separable AR1 spatial field" = list(
      fit = function() fb_inla(y ~ 1, random = ~ ar1(row):ar1(col), data = f$d_capgrid),
      link = "identity", y = f$d_capgrid$y, slots = "sigma:sigma"
    ),
    "Univariate P-spline" = list(
      fit = function() fb_inla(y ~ 1, random = ~ spl(x), data = f$d_cap),
      link = "identity", y = f$d_cap$y, slots = "sigma:sigma"
    ),
    "Observation weights (Gaussian, identity link)" = list(
      fit = function() {
        fb_inla(y ~ env, random = ~ gen, data = f$d_cap, weights = f$w_cap)
      },
      link = "identity", y = f$d_cap$y, slots = c("sigma:sigma", "sd:gen")
    ),
    "Exact sufficient-statistic aggregation (binomial, Poisson)" = list(
      fit = function() {
        fb_inla(cnt ~ env, random = ~ gen, data = f$d_agg, aggregate = TRUE,
                family = "poisson")
      },
      link = "log", y = NULL, slots = "sd:gen"
    ),
    "Per-trial separable AR1 field" = list(
      fit = function() {
        fb_inla(y ~ 1, random = ~ at(trial):ar1(row):ar1(col),
                data = f$d_ar1_at_trial)
      },
      link = "identity", y = f$d_ar1_at_trial$y, slots = "sigma:sigma"
    ),
    "Grouped temporal AR1 field" = list(
      fit = function() fb_inla(y ~ 1, random = ~ ar1(t):g, data = f$d_ar1_tg),
      link = "identity", y = f$d_ar1_tg$y, slots = "sigma:sigma"
    )
  )
}

# Each capability-matrix row that fits on brms through flexybayes(): the Stan
# program's call, the link, the response, whether the slope slot is divided,
# and whether the residual is a log-sigma predictor.
.hn1_brms_classes <- function(f) {
  fb_code <- function(...) {
    suppressMessages(flexybayes(..., backend = "brms", return_code = TRUE))
  }
  list(
    "Gaussian LMM, simple random intercept" = list(
      code = function() fb_code(y ~ env, random = ~ gen, data = f$d_cap),
      link = "identity", y = f$d_cap$y
    ),
    "GLMM (binomial, Poisson, negative binomial, gamma, beta), simple random effect" = list(
      code = function() {
        fb_code(bin ~ env, random = ~ gen, data = f$d_cap, family = "binomial")
      },
      link = "logit", y = NULL
    ),
    "Hurdle gamma (zero mass plus a positive gamma part)" = list(
      code = function() {
        fb_code(hg ~ env, random = ~ gen, data = f$d_cap, family = "hurdle_gamma")
      },
      link = "log", y = NULL
    ),
    "Uncorrelated random slope" = list(
      code = function() fb_code(y ~ x + (x || gen), data = f$d_cap),
      link = "identity", y = f$d_cap$y, slope_x = f$d_cap$x
    ),
    "Factor-by-numeric fixed interaction" = list(
      code = function() fb_code(y ~ env * x, random = ~ gen, data = f$d_cap),
      link = "identity", y = f$d_cap$y
    ),
    # (x | g) has no INLA counterpart, so its slope SD takes the scale
    # undivided, as the intercept SD does.
    "Correlated random slope" = list(
      code = function() fb_code(y ~ x + (x | gen), data = f$d_cap),
      link = "identity", y = f$d_cap$y
    ),
    "Nested / interaction random effects, multi-stratum" = list(
      code = function() fb_code(y ~ env, random = ~ gen:env, data = f$d_cap),
      link = "identity", y = f$d_cap$y
    ),
    "Heterogeneous variance by factor level" = list(
      code = function() fb_code(y ~ env, random = ~ diag(env):gen, data = f$d_cap),
      link = "identity", y = f$d_cap$y
    ),
    "Unstructured genotype-by-environment covariance" = list(
      code = function() fb_code(y ~ env, random = ~ us(env):gen, data = f$d_cap),
      link = "identity", y = f$d_cap$y
    ),
    "Heterogeneous residual by factor level" = list(
      code = function() {
        fb_code(y ~ env, random = ~ gen, residual = ~ dsum(~ units | env),
                data = f$d_cap)
      },
      link = "identity", y = f$d_cap$y, log_sigma = TRUE
    ),
    "Baseline heteroscedastic MET: iid genotype, iid genotype-by-environment and a sectioned residual" = list(
      code = function() {
        fb_code(y ~ env, random = ~ gen + gen:env,
                residual = ~ dsum(~ units | env), data = f$d_cap)
      },
      link = "identity", y = f$d_cap$y, log_sigma = TRUE
    ),
    "Known-covariance genomic / pedigree random effect" = list(
      code = function() {
        fb_code(y ~ env, random = ~ vm(gen, K), data = f$d_cap,
                known_matrices = list(K = f$k_cap))
      },
      link = "identity", y = f$d_cap$y
    ),
    "Observation weights (Gaussian, identity link)" = list(
      code = function() {
        fb_code(y ~ env, random = ~ gen, data = f$d_cap, weights = f$w_cap)
      },
      link = "identity", y = f$d_cap$y, log_sigma = TRUE
    )
  )
}


# --- the default object ---------------------------------------------------- #

test_that("the default is a half-normal whose scale is the former bound over the adopted divisor", {
  expect_identical(flexyBayes:::.FB_DEFAULT_SCALE_DIVISOR, .hn1_divisor)
  set.seed(5L)
  d <- data.frame(y = stats::rnorm(30L, 4, 3), g = factor(rep(1:5, 6L)))

  p <- flexyBayes:::.default_half_normal_prior(
    data = d, response = "y", family = "gaussian",
    random_groups = "g", vm_ped_groups = "k"
  )
  families <- vapply(p$specs, function(s) s$spec$family, character(1L))
  scales <- vapply(p$specs, function(s) s$spec$args$scale, numeric(1L))
  expect_identical(unique(families), "half_normal")
  expect_equal(scales, rep(5 * stats::sd(d$y) / .hn1_divisor, 3L), tolerance = 1e-12)
  expect_identical(attr(p, "fb_prior_default_basis"), "identity_sd_half_normal")
  expect_equal(attr(p, "fb_prior_default_scale"), 5 * stats::sd(d$y) / .hn1_divisor,
               tolerance = 1e-12)

  links <- list(
    poisson = c(3, "log_link_half_normal"),
    negative_binomial = c(3, "log_link_half_normal"),
    gamma = c(3, "log_link_half_normal"),
    binomial = c(5, "logit_half_normal"),
    beta = c(5, "logit_half_normal")
  )
  for (fam in names(links)) {
    s <- flexyBayes:::.default_half_normal_scale(NULL, "y", fam)
    expect_equal(s$scale, as.numeric(links[[fam]][[1L]]) / .hn1_divisor,
                 label = paste(fam, "scale"))
    expect_identical(s$basis, links[[fam]][[2L]], label = paste(fam, "basis"))
  }
})


# --- INLA ------------------------------------------------------------------- #

test_that("every default slot on INLA carries the closed-form half-normal, by model class", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  f <- .hn1_data()
  classes <- .hn1_inla_classes(f)

  # The classes are the capability matrix's INLA "fits" rows, every one.
  cap <- flexyBayes:::.fb_capability_matrix()
  expect_setequal(names(classes), cap$model_class[cap$inla == "fits"])

  theta <- c(-4, -1, 0, 1, 4)
  for (cls in names(classes)) {
    spec <- classes[[cls]]
    fit <- spec$fit()
    slots <- fit$extras$manifest$slots
    std <- fit$extras$response_standardisation$scale %||% 1
    if (!is.null(fit$extras$response_standardisation)) {
      expect_equal(std, stats::sd(spec$y, na.rm = TRUE), tolerance = 1e-12,
                   label = paste(cls, "standardisation scale"))
    }
    package <- slots[
      slots$source == "package" & !is.na(slots$input_id) &
        grepl("^(sigma|sd):", slots$input_id),
      ,
      drop = FALSE
    ]
    expect_setequal(unique(package$input_id), spec$slots)
    expect_false(any(grepl("U=", slots$readback, fixed = TRUE)), label = cls)
    expect_identical(NROW(fit$extras$support_audit), 0L, label = cls)

    for (i1 in seq_len(nrow(package))) {
      id <- package$input_id[[i1]]
      text <- package$readback[[i1]]
      lab <- paste(cls, id, package$engine_slot[[i1]])
      expect_true(package$readback_agrees[[i1]], label = lab)
      expect_match(text, "^expression: s=", label = lab)
      expected <- .hn1_scale(spec$link, spec$y)
      if (!is.null(spec$slope) && identical(id, spec$slope$id)) {
        observed_x <- spec$slope$x[!is.na(spec$y)]
        expected <- expected / stats::sd(observed_x)
      }
      s_received <- as.numeric(sub("^expression: s=([^;]+);.*$", "\\1", text))
      expect_lt(abs(s_received * std / expected - 1), 1e-12, label = paste(lab, "scale"))
      expect_lt(
        max(abs(.hn1_eval_expression(text, theta) -
          .hn1_closed_form(theta, s_received))),
        1e-10,
        label = paste(lab, "density")
      )
    }
  }
})


# --- brms ------------------------------------------------------------------- #

test_that("every default slot on brms carries normal(0, s) on a standard deviation bounded only below", {
  skip_on_cran()
  skip_if_not_installed("brms")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  f <- .hn1_data()
  classes <- .hn1_brms_classes(f)

  # The classes are the capability matrix's brms "fits" rows that a
  # flexybayes() call reaches; the joint shrinkage rows are fb_select()'s.
  cap <- flexyBayes:::.fb_capability_matrix()
  reachable <- cap$model_class[
    cap$brms == "fits" & !grepl("shrinkage prior", cap$model_class)
  ]
  expect_setequal(names(classes), reachable)

  num <- "([-+0-9.eE]+)"
  for (cls in names(classes)) {
    spec <- classes[[cls]]
    lines <- strsplit(as.character(spec$code()), "\n", fixed = TRUE)[[1L]]
    s <- .hn1_scale(spec$link, spec$y)
    expected <- s
    if (!is.null(spec$slope_x)) {
      expected <- c(s, s / stats::sd(spec$slope_x))
    }

    declared <- grep(
      "^ *(real|vector)<[^>]*>(\\[[^]]*\\])? *(sigma|sd_[0-9]+);",
      lines,
      value = TRUE
    )
    expect_gt(length(declared), 0L)
    expect_false(any(grepl("upper=", declared, fixed = TRUE)), label = cls)
    expect_false(any(grepl("uniform_lpdf\\((sigma|sd_)", lines)), label = cls)

    priors <- grep(
      paste0("lprior \\+= normal_lpdf\\((sigma|sd_[0-9]+)(\\[[0-9]+\\])? \\| 0, ", num, "\\)"),
      lines,
      value = TRUE
    )
    got <- as.numeric(sub(
      paste0(".*normal_lpdf\\([^|]+\\| 0, ", num, "\\).*"), "\\1", priors
    ))
    expect_gt(length(got), 0L)
    for (v in got) {
      expect_lt(min(abs(v / expected - 1)), 1e-6, label = paste(cls, "scale", v))
    }
    for (v in expected) {
      expect_lt(min(abs(got / v - 1)), 1e-6, label = paste(cls, "expected scale", v))
    }
    if (isTRUE(spec$log_sigma)) {
      loc_line <- grep(
        paste0("normal_lpdf\\((b_sigma|Intercept_sigma) \\| ", num, ", 1\\)"),
        lines,
        value = TRUE
      )
      expect_length(loc_line, 1L)
      loc <- as.numeric(sub(paste0(".*\\| ", num, ", 1\\).*"), "\\1", loc_line))
      expect_lt(abs(loc - log(stats::sd(spec$y))), 1e-5, label = paste(cls, "log-sigma location"))
    }
  }
})


# --- both engines ----------------------------------------------------------- #

test_that("the default reaches INLA and brms with the same scale on every slot", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  f <- .hn1_data()

  fit_i <- suppressMessages(flexybayes(
    y ~ x + (x || gen), data = f$d_cap, backend = "inla", verbose = FALSE
  ))
  # The block reads the priors each engine received, not the posterior, so
  # the short chain's sampler warnings carry nothing it checks.
  fit_b <- suppressWarnings(suppressMessages(flexybayes(
    y ~ x + (x || gen), data = f$d_cap, backend = "brms",
    chains = 1L, warmup = 200L, n_samples = 200L, seed = 1L, verbose = FALSE
  )))

  s <- 5 * stats::sd(f$d_cap$y) / .hn1_divisor
  oracle <- c(sigma = s, sd_gen = s, sd_x_gen = s / stats::sd(f$d_cap$x))

  slots <- fit_i$extras$manifest$slots
  std <- fit_i$extras$response_standardisation$scale
  inla_scale <- function(id) {
    text <- slots$readback[which(slots$input_id == id)[[1L]]]
    as.numeric(sub("^expression: s=([^;]+);.*$", "\\1", text)) * std
  }
  inla <- c(
    sigma = inla_scale("sigma:sigma"),
    sd_gen = inla_scale("sd:gen"),
    sd_x_gen = inla_scale("sd:x_gen")
  )

  ps <- as.data.frame(brms::prior_summary(fit_b$brms))
  ps <- ps[ps$source == "user", , drop = FALSE]
  brms_scale <- function(cls, coef = "", group = "") {
    row <- ps$prior[ps$class == cls & ps$coef == coef & ps$group == group]
    expect_length(row, 1L)
    as.numeric(sub("^normal\\(0, ([^)]+)\\)$", "\\1", row))
  }
  brms <- c(
    sigma = brms_scale("sigma"),
    sd_gen = brms_scale("sd", group = "gen"),
    sd_x_gen = brms_scale("sd", coef = "x", group = "gen")
  )

  for (nm in names(oracle)) {
    expect_lt(abs(inla[[nm]] / oracle[[nm]] - 1), 1e-12, label = paste("INLA", nm))
    expect_lt(abs(brms[[nm]] / oracle[[nm]] - 1), 1e-6, label = paste("brms", nm))
    expect_lt(abs(brms[[nm]] / inla[[nm]] - 1), 1e-6, label = paste("engines", nm))
  }

  keys <- names(oracle)
  expect_identical(
    fit_b$extras$fingerprint$priors[keys],
    fit_i$extras$fingerprint$priors[keys]
  )
  expect_identical(
    fit_b$extras$fingerprint$prior_sources[keys],
    fit_i$extras$fingerprint$prior_sources[keys]
  )
  brms_slots <- fit_b$extras$manifest$slots
  expect_true(all(brms_slots$readback_agrees[brms_slots$source == "package"]))
})


# --- the design that tripped the uniform ----------------------------------- #

test_that("john.alpha rep + rep:block fits on INLA under the default, REML inside every interval", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
  skip_if_not_installed("lme4")
  withr::local_options(flexyBayes.silence_default_prior_note = TRUE)
  d <- agridat::john.alpha

  reml <- lme4::lmer(yield ~ gen + (1 | rep) + (1 | rep:block), data = d, REML = TRUE)
  vc <- as.data.frame(lme4::VarCorr(reml))
  reml_sd <- c(
    sigma = vc$sdcor[vc$grp == "Residual"],
    sd_rep = vc$sdcor[vc$grp == "rep"],
    "sd_rep:block" = vc$sdcor[vc$grp == "rep:block"]
  )

  # Under the uniform default every one of eight repeats refused with
  # prior_support_violated (`rep` interval 0.17 to 3.65 against 3.05).
  for (i1 in seq_len(8L)) {
    fit <- suppressMessages(flexybayes(
      yield ~ gen, random = ~ rep + rep:block, data = d,
      backend = "inla", verbose = FALSE
    ))
    expect_identical(NROW(fit$extras$support_audit), 0L)
    v <- summary(fit)$varcomp
    for (nm in names(reml_sd)) {
      row <- v[v$component == nm, , drop = FALSE]
      expect_length(row$estimate, 1L)
      expect_true(
        row$conf.low <= reml_sd[[nm]] && reml_sd[[nm]] <= row$conf.high,
        label = sprintf("repeat %d: REML %s = %.4f inside [%.4f, %.4f]",
                        i1, nm, reml_sd[[nm]], row$conf.low, row$conf.high)
      )
    }
  } # ends i1, over repeats
})
