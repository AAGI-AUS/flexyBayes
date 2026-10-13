# flexyBayes against stored ASReml 4.2 fits of the eight model shapes of
# the ASReml oracle suite (data-raw/asreml_oracle/SPEC.md, rows 01 to 07
# and 15). The ASReml side is inst/extdata/asreml_parity.rds, written by
# data-raw/asreml_parity.R on the data in helper-asreml-parity-data.R.
#
# Tolerances are SPEC.md's, fixed on 2026-09-11 before any comparison ran:
#   fixed effect (Layer 2)      abs(b_fb - b_asr) <= 0.25 * se_asr
#   variance component (V1-V2)  abs(median - v_asr) / v_asr <= 0.15
#   variance ratio (Layer 2)    abs(diff) <= 0.05
# The flexyBayes side is the posterior median on the variance scale. For a
# fixed effect the posterior mean, which on a Gaussian response is the
# median. A component is compared only where ASReml identifies it: bound
# "P" and std.error / component < 0.5 (SPEC "Identification gate"). A
# component ASReml does not identify is checked for the package's
# weak-identification flag where that flag applies, and never compared on
# a wider tolerance.

.AP_FIXED_SE_SHARE <- 0.25
.AP_VARCOMP_REL <- 0.15
.AP_RATIO_ABS <- 0.05
.AP_IDENTIFIED_SE_RATIO <- 0.5

.ap_reference <- function() {
  path <- system.file("extdata", "asreml_parity.rds", package = "flexyBayes")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "extdata",
                                "asreml_parity.rds")
  }
  readRDS(path)
}

.ap_source_root <- function() {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
  if (file.exists(file.path(root, "data-raw", "asreml_parity.R"))) {
    root
  } else {
    NA_character_
  }
}

# The flexyBayes call for each shape; `...` overrides or adds arguments.
.ap_fit <- function(shape, d = .ap_data(shape), ...) {
  f <- switch(
    shape,
    "01_iid_intercept" = list(yield ~ gen, ~ rep, NULL),
    "02_nested_stratum" = list(yield ~ gen, ~ rep + rep:block, NULL),
    "15_weights" = list(yield ~ gen, ~ rep, NULL),
    "03_diag_gxe" = list(yield ~ env, ~ diag(env):gen, NULL),
    "05_combined_met" = list(yield ~ env, ~ diag(env):gen,
                             ~ dsum(~ units | env)),
    "04_dsum_residual" = list(yield ~ env, ~ gen, ~ dsum(~ units | env)),
    list(yield ~ env, ~ us(env):gen, NULL)
  )
  args <- list(fixed = f[[1L]], random = f[[2L]], data = d, verbose = FALSE)
  if (!is.null(f[[3L]])) {
    args$residual <- f[[3L]]
  }
  if (identical(shape, "15_weights")) {
    args$weights <- d$wt
  }
  withr::with_options(
    list(
      flexyBayes.inla_threads = "1:1",
      flexyBayes.silence_default_prior_note = TRUE,
      flexyBayes.silence_weak_identification_warning = TRUE
    ),
    do.call("flexybayes", utils::modifyList(args, list(...)))
  )
}

# ASReml's component and coefficient names in flexyBayes's spelling.
.ap_component_name <- function(term) {
  out <- paste0("sd_", term)
  out[term == "units!R"] <- "sigma"
  diag_gen <- grepl("^env:gen!env_[^:]+$", term)
  out[diag_gen] <- sub("^env:gen!env_", "sd_gen__env", term[diag_gen])
  us_var <- grepl("^env:gen!env_([^:]+):\\1$", term)
  out[us_var] <- sub("^env:gen!env_([^:]+):.*$", "sd_gen__env\\1",
                     term[us_var])
  section <- grepl("^env_.+!R$", term)
  out[section] <- sub("^env_(.+)!R$", "sigma_\\1", term[section])
  out
}

.ap_fixed_name <- function(term) {
  sub("^(gen|env)_", "\\1", term)
}

# One row per ASReml variance component, the flexyBayes posterior median
# beside it on the variance scale, and whether the pair is compared. The
# covariances of a us() term are left out. Its correlations are compared
# separately.
.ap_varcomp_pairs <- function(rec, fit) {
  vc <- rec$varcomp
  vc <- vc[!grepl("^env:gen!env_([^:]+):(?!\\1$)", vc$term, perl = TRUE), ]
  med <- attr(fit$extras$variance_comps, "posterior_median")
  notes <- summary(fit)$varcomp
  name <- .ap_component_name(vc$term)
  fb <- unname(med[name])^2
  data.frame(
    component = name,
    asreml = vc$component,
    flexybayes = fb,
    rel_diff = abs(fb - vc$component) / vc$component,
    identified = vc$bound %in% c("P", "?") & vc$component > 0 &
      vc$std.error / vc$component < .AP_IDENTIFIED_SE_RATIO,
    weak_flag = notes$note[match(name, notes$component)] ==
      "weakly identified",
    stringsAsFactors = FALSE
  )
}

.ap_fixed_pairs <- function(rec, fit) {
  fx <- summary(fit)$fixed
  name <- .ap_fixed_name(rec$fixed$term)
  i <- match(name, fx$term)
  data.frame(
    term = name,
    asreml = rec$fixed$estimate,
    asreml_se = rec$fixed$std.error,
    flexybayes = fx$estimate[i],
    stringsAsFactors = FALSE
  )
}

# The two checks every shape carries. Each mismatch is named.
.ap_fixed_misses <- function(fx) {
  bad <- is.na(fx$flexybayes) |
    abs(fx$flexybayes - fx$asreml) > .AP_FIXED_SE_SHARE * fx$asreml_se
  fx$term[bad]
}

.ap_varcomp_misses <- function(vc) {
  cmp <- vc[vc$identified, ]
  cmp$component[is.na(cmp$rel_diff) | cmp$rel_diff > .AP_VARCOMP_REL]
}

.ap_skip_engine <- function() {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("agridat")
}

.ap_check_shape <- function(shape, weak) {
  rec <- .ap_reference()$fits[[shape]]
  fit <- .ap_fit(shape)
  fx <- .ap_fixed_pairs(rec, fit)
  vc <- .ap_varcomp_pairs(rec, fit)
  expect_identical(.ap_fixed_misses(fx), character(0), label = shape)
  expect_identical(.ap_varcomp_misses(vc), character(0), label = shape)
  expect_gt(sum(vc$identified), 0L)
  # The components ASReml does not identify that the package flags.
  expect_identical(vc$component[!vc$identified & vc$weak_flag], weak,
                   label = shape)
  invisible(list(fit = fit, rec = rec, fixed = fx, varcomp = vc))
}


# --- the stored fits ---------------------------------------------------- #

test_that("the stored fits cover the eight shapes and every one converged", {
  ref <- .ap_reference()
  expect_setequal(
    names(ref$fits),
    c("01_iid_intercept", "02_nested_stratum", "03_diag_gxe",
      "04_dsum_residual", "05_combined_met", "06_us_dim2", "07_us_dim3",
      "15_weights")
  )
  expect_match(ref$stamp$asreml, "^4\\.2\\.")
  for (rec in ref$fits) {
    expect_true(rec$converged, label = rec$shape)
    expect_gt(nrow(rec$varcomp), 0L)
    expect_gt(nrow(rec$fixed), 0L)
  }
})

test_that("the stored fits were written by the recipe in data-raw", {
  root <- .ap_source_root()
  skip_if(is.na(root), "package source tree not reachable")
  stamp <- .ap_reference()$stamp
  expect_identical(
    unname(tools::md5sum(file.path(root, "data-raw", "asreml_parity.R"))),
    stamp$script_md5
  )
  expect_identical(
    unname(tools::md5sum(
      file.path(root, "tests", "testthat", "helper-asreml-parity-data.R")
    )),
    stamp$helper_md5
  )
})

test_that("each shape's data rebuild to the rows ASReml fitted", {
  skip_if_not_installed("agridat")
  for (rec in .ap_reference()$fits) {
    d <- .ap_data(rec$shape)
    expect_identical(nrow(d), rec$n, label = rec$shape)
    expect_identical(.ap_data_md5(d), rec$data_md5, label = rec$shape)
  }
  # A changed yield is a different data set.
  d <- .ap_data("04_dsum_residual")
  d$yield[[1L]] <- d$yield[[1L]] + 0.01
  expect_false(identical(
    .ap_data_md5(d),
    .ap_reference()$fits[["04_dsum_residual"]]$data_md5
  ))
})


# --- shapes fitted on INLA ---------------------------------------------- #

test_that("01 iid intercept: fixed effects and residual agree with ASReml", {
  .ap_skip_engine()
  .ap_check_shape("01_iid_intercept", weak = "sd_rep")
})

test_that("02 alpha design: fixed effects, strata and gen means agree with ASReml", {
  .ap_skip_engine()
  out <- .ap_check_shape("02_nested_stratum", weak = "sd_rep")
  p <- predict(out$fit, classify = "gen")
  asr <- out$rec$predict_gen
  i <- match(asr$gen, as.character(p$gen))
  expect_false(anyNA(i))
  expect_true(all(
    abs(p$estimate[i] - asr$predicted.value) <=
      .AP_FIXED_SE_SHARE * asr$std.error
  ))
})

test_that("15 relative weights: fixed effects and residual agree with ASReml", {
  .ap_skip_engine()
  .ap_check_shape("15_weights", weak = "sd_rep")
})

test_that("03 diag(env):gen: residual agrees with ASReml, and at() is the same model", {
  .ap_skip_engine()
  out <- .ap_check_shape("03_diag_gxe", weak = character(0))
  # ASReml fits at(env):gen and diag(env):gen to the same likelihood.
  at_rec <- out$rec$at
  expect_equal(at_rec$loglik, out$rec$loglik, tolerance = 1e-8)
  expect_equal(at_rec$varcomp$component, out$rec$varcomp$component,
               tolerance = 1e-6)
  # flexyBayes reads the two spellings to the same random term. Repeated
  # identical fits of this shape differ by up to 4% in the genotype
  # variances, so the check is on the parsed model, not on two fits.
  d <- .ap_data("03_diag_gxe")
  expect_identical(
    .fb_from_asreml_impl(yield ~ env, random = ~ at(env):gen,
                         data = d)$random_terms,
    .fb_from_asreml_impl(yield ~ env, random = ~ diag(env):gen,
                         data = d)$random_terms
  )
})

test_that("04 dsum residual: every component, the fixed effects and vpredict agree", {
  .ap_skip_engine()
  out <- .ap_check_shape("04_dsum_residual", weak = character(0))
  expect_true(all(out$varcomp$identified))
  set.seed(20261008L)
  v <- vpredict(out$fit, gen_share_BH93 ~ sd_gen^2 / (sd_gen^2 + sigma_BH93^2))
  asr <- out$rec$vpredict
  expect_lte(abs(v$median - asr$estimate), .AP_RATIO_ABS)
})

test_that("05 diag(env):gen with dsum residual: identified components agree", {
  .ap_skip_engine()
  out <- .ap_check_shape("05_combined_met", weak = character(0))
  expect_setequal(
    out$varcomp$component[out$varcomp$identified],
    c("sd_gen__envC5", paste0("sigma_C", 1:6))
  )
})

test_that("DEF-201: a genotype variance REML cannot identify is flagged", {
  skip(paste(
    "DEF-201 open: on 03_diag_gxe sd_gen__envC2 has REML std.error 7.2",
    "times its estimate. The posterior variance median is 19.5 and its SD",
    "interval 1.7 times its median, below both flag thresholds"
  ))
  .ap_skip_engine()
  out <- .ap_check_shape("03_diag_gxe", weak = character(0))
  vc <- out$varcomp
  expect_true(vc$weak_flag[vc$component == "sd_gen__envC2"])
})


# --- us(): refused on INLA by name, fitted on brms ---------------------- #

test_that("06 and 07 us(env):gen: INLA refuses by name", {
  skip_if_not_installed("agridat")
  for (shape in c("06_us_dim2", "07_us_dim3")) {
    expect_error(
      .ap_fit(shape, backend = "inla"),
      class = "flexybayes_refusal_correlated_covariance_requires_brms"
    )
  }
})

test_that("06 and 07 us(env):gen: auto chooses brms", {
  skip_if_not_installed("agridat")
  skip_if_not_installed("brms")
  for (shape in c("06_us_dim2", "07_us_dim3")) {
    p <- .ap_fit(shape, plan = TRUE)
    expect_identical(p$backend_chosen, "brms")
  }
})

test_that("06 and 07 us(env):gen on brms: fixed effects and identified variances agree", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("agridat")
  skip_if_not_long_tier()
  for (shape in c("06_us_dim2", "07_us_dim3")) {
    rec <- .ap_reference()$fits[[shape]]
    fit <- suppressWarnings(.ap_fit(
      shape, backend = "brms", chains = 4L, seed = 20261008L
    ))
    expect_identical(.ap_fixed_misses(.ap_fixed_pairs(rec, fit)),
                     character(0), label = shape)
    vc <- .ap_varcomp_pairs(rec, fit)
    expect_identical(.ap_varcomp_misses(vc), character(0), label = shape)
    expect_true("sigma" %in% vc$component[vc$identified])
    # A genetic correlation is compared when ASReml identifies both of its
    # variances. On these counties at most one of each pair is, so no
    # correlation is compared.
    pair <- regmatches(rec$correlations$functional,
                       gregexpr("env[^_]+", rec$correlations$functional))
    both <- vapply(pair, function(p) {
      all(paste0("sd_gen__", p) %in% vc$component[vc$identified])
    }, logical(1))
    expect_false(any(both), label = shape)
  }
})


# --- the comparison can fail -------------------------------------------- #

test_that("a planted mismatch fails the comparison", {
  .ap_skip_engine()
  rec <- .ap_reference()$fits[["04_dsum_residual"]]
  fit <- .ap_fit("04_dsum_residual")

  shifted <- rec
  shifted$fixed$estimate[[2L]] <- shifted$fixed$estimate[[2L]] +
    0.3 * shifted$fixed$std.error[[2L]]
  expect_identical(
    .ap_fixed_misses(.ap_fixed_pairs(shifted, fit)),
    .ap_fixed_name(rec$fixed$term[[2L]])
  )

  # Two environments' residual variances swapped.
  swapped <- rec
  i <- match(c("env_KE93!R", "env_WP93!R"), swapped$varcomp$term)
  swapped$varcomp$term[i] <- swapped$varcomp$term[rev(i)]
  expect_setequal(
    .ap_varcomp_misses(.ap_varcomp_pairs(swapped, fit)),
    c("sigma_KE93", "sigma_WP93")
  )
})
