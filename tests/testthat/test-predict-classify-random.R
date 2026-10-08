# =============================================================================
# predict(classify = ) with random terms, against ASReml 4.2.
#
# The expected values are ASReml's predicted values and standard errors on
# the same rows, stored by data-raw/asreml_classify_random.R in
# fixtures/asreml_classify_random.rds. The data are stored with them.
#
# Tolerances, per level (written before the comparison was run):
#   * |posterior mean - predicted value| <= 0.25 ASReml SE. The BLUP part
#     is shrunk by a variance ratio that ASReml fixes at REML and the
#     posterior averages over; a ratio a few percent off moves a BLUP of
#     up to 2 SE by 0.1-0.2 SE.
#   * Spearman correlation with the predicted values >= 0.99.
#   * posterior SD / ASReml SE within [0.95, 1.30]. ASReml's SE treats the
#     variances as known; the posterior SD also carries their uncertainty.
#
# The gen:county table broke the point bound: up to 0.36 SE on the MET fit
# and 0.29 SE with a cell removed. With six counties, INLA's posterior mean
# of the gen:county variance is 21.3 against REML's 15.0, so the cell
# effects are shrunk less. The test with the variances held at their REML
# values shows the construction itself agrees to 0.006 SE, so the
# gen:county bound is set at 0.40 SE and the reason stated here.
# =============================================================================

.cr_ref <- function() {
  readRDS(test_path("fixtures", "asreml_classify_random.rds"))
}

.cr_cache <- new.env(parent = emptyenv())

.cr_fit <- function(key, formula, random, data, backend = "inla", ...) {
  id <- paste(key, backend, sep = "_")
  if (is.null(.cr_cache[[id]])) {
    withr::local_options(
      flexyBayes.inla_threads = "1:1",
      flexyBayes.silence_default_prior_note = TRUE
    )
    .cr_cache[[id]] <- suppressMessages(suppressWarnings(flexybayes(
      formula,
      random = random,
      data = data,
      backend = backend,
      verbose = FALSE,
      mcmc_verbose = FALSE,
      ...
    )))
  }
  .cr_cache[[id]]
}

# One message per violated bound; character(0) when the table agrees.
.cr_check <- function(tab, ref, by, point = 0.25, width = c(0.95, 1.30),
                      rank = 0.99) {
  m <- merge(as.data.frame(tab), ref, by = by, suffixes = c("", ".asreml"))
  if (nrow(m) != nrow(ref)) {
    return(sprintf("%d of %d levels matched", nrow(m), nrow(ref)))
  }
  dz <- abs(m$estimate - m$predicted.value) / m$std.error.asreml
  ratio <- m$std.error / m$std.error.asreml
  rho <- stats::cor(m$estimate, m$predicted.value, method = "spearman")
  msgs <- character(0)
  if (max(dz) > point) {
    msgs <- c(msgs, sprintf("point: max %.3f SE > %.2f", max(dz), point))
  }
  if (min(ratio) < width[[1L]] || max(ratio) > width[[2L]]) {
    msgs <- c(msgs, sprintf("width: SD/SE in [%.3f, %.3f]",
                            min(ratio), max(ratio)))
  }
  if (rho < rank) {
    msgs <- c(msgs, sprintf("rank: Spearman %.4f < %.2f", rho, rank))
  }
  msgs
}

.cr_rcb <- function() {
  .cr_fit("rcb", yield ~ rep, ~gen, .cr_ref()$data$rcb)
}

.cr_met <- function() {
  .cr_fit("met", yield ~ county + county:rep, ~ gen + gen:county,
          .cr_ref()$data$met)
}


# ---------------------------------------------------------------- #
# 1. The check can fail                                              #
# ---------------------------------------------------------------- #

test_that("the oracle check flags a shifted, narrowed or reordered table", {
  ref <- .cr_ref()$rcb$gen
  good <- data.frame(
    gen = ref$gen,
    estimate = ref$predicted.value,
    std.error = 1.05 * ref$std.error
  )
  expect_identical(.cr_check(good, ref, "gen"), character(0))

  shifted <- good
  shifted$estimate[[3L]] <- shifted$estimate[[3L]] + 0.3 * ref$std.error[[3L]]
  expect_match(.cr_check(shifted, ref, "gen"), "^point")

  narrow <- good
  narrow$std.error <- 0.8 * ref$std.error
  expect_match(.cr_check(narrow, ref, "gen"), "^width")

  reordered <- good
  reordered$estimate <- rev(reordered$estimate)
  expect_true(any(grepl("^rank", .cr_check(reordered, ref, "gen"))))
})


# ---------------------------------------------------------------- #
# 2. INLA against ASReml                                             #
# ---------------------------------------------------------------- #

test_that("RCB, gen random: classify = 'gen' agrees with ASReml on INLA", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  ref <- .cr_ref()
  tab <- predict(.cr_rcb(), classify = "gen")

  expect_s3_class(tab, "fb_predict_classify")
  expect_identical(
    names(tab),
    c("gen", "estimate", "std.error", "conf.low", "conf.high")
  )
  expect_identical(nrow(tab), 24L)
  expect_identical(.cr_check(tab, ref$rcb$gen, "gen"), character(0))
  expect_identical(attr(tab, "random"), "gen")
  expect_identical(attr(tab, "average"), "rep")
  expect_true(all(tab$conf.low < tab$estimate & tab$estimate < tab$conf.high))
})

test_that("a random term naming a factor outside classify is ignored", {
  # rep:block names block, which is in neither classify nor the fixed part.
  skip_if_not_installed("INLA")
  skip_on_cran()
  ref <- .cr_ref()
  fit <- .cr_fit("rcb_blk", yield ~ rep, ~ gen + rep:block, ref$data$rcb)
  tab <- predict(fit, classify = "gen")

  expect_identical(.cr_check(tab, ref$rcb_blk$gen, "gen"), character(0))
  expect_identical(attr(tab, "ignored"), "block")
})

test_that("MET: classify = 'gen' averages gen:county over county, as ASReml", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  ref <- .cr_ref()
  tab <- predict(.cr_met(), classify = "gen")

  expect_identical(.cr_check(tab, ref$met$gen, "gen"), character(0))
  expect_identical(attr(tab, "random"), c("gen", "gen:county"))
  expect_identical(attr(tab, "average"), c("county", "rep"))
})

test_that("MET: classify = 'gen:county' agrees with ASReml on INLA", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  ref <- .cr_ref()
  tab <- predict(.cr_met(), classify = "gen:county")

  expect_identical(nrow(tab), 96L)
  expect_identical(
    .cr_check(tab, ref$met$gen_county, c("gen", "county"), point = 0.40),
    character(0)
  )
})

test_that("a gen:county cell with no data is predicted as ASReml does", {
  # Every plot of G01 in C1 is removed. The cell's effect has no node in
  # INLA's field and is drawn from the term's distribution.
  skip_if_not_installed("INLA")
  skip_on_cran()
  ref <- .cr_ref()
  fit <- .cr_fit("met_gap", yield ~ county + county:rep, ~ gen + gen:county,
                 ref$data$met_gap)

  expect_identical(
    .cr_check(predict(fit, classify = "gen"), ref$met_gap$gen, "gen"),
    character(0)
  )
  tab <- predict(fit, classify = "gen:county")
  expect_identical(
    .cr_check(tab, ref$met_gap$gen_county, c("gen", "county"), point = 0.40),
    character(0)
  )
  # The empty cell is wider than G01's other cells, as on ASReml's side.
  g01 <- tab[tab$gen == "G01", ]
  expect_gt(g01$std.error[g01$county == "C1"],
            max(g01$std.error[g01$county != "C1"]))
})

test_that("with the variances held at REML, the table reproduces ASReml", {
  # Each standard deviation is held within 0.5% of ASReml's REML value, so
  # SEs move by at most 0.5% and shrinkage by about 1%: 0.02 SE and a
  # width ratio within 1%.
  skip_if_not_installed("INLA")
  skip_on_cran()
  ref <- .cr_ref()
  vc <- sqrt(ref$met$varcomp)
  pin <- function(s) call("uniform", lower = 0.995 * s, upper = 1.005 * s)
  prior <- eval(bquote(fb_prior(
    sigma ~ .(pin(vc[["units!R"]])),
    sd(group = "gen") ~ .(pin(vc[["gen"]])),
    sd(group = "gen:county") ~ .(pin(vc[["gen:county"]]))
  )))
  fit <- .cr_fit("met_pinned", yield ~ county + county:rep,
                 ~ gen + gen:county, ref$data$met, prior = prior)

  expect_identical(
    .cr_check(predict(fit, classify = "gen"), ref$met$gen, "gen",
              point = 0.02, width = c(0.99, 1.01)),
    character(0)
  )
  expect_identical(
    .cr_check(predict(fit, classify = "gen:county"), ref$met$gen_county,
              c("gen", "county"), point = 0.02, width = c(0.99, 1.01)),
    character(0)
  )
})

test_that("DEF-201: gen fixed, gen:county random includes gen:county", {
  # ASReml's rules use every term built from the classify and averaged
  # variables, so gen:county enters predict(classify = "gen") averaged over
  # county. The table used to be the fixed part alone, 11% wider than
  # ASReml's SE on this fit. The variances are held within 0.5% of REML, as
  # in the test above, so the width ratio must be within 1%.
  skip_if_not_installed("INLA")
  skip_on_cran()
  ref <- .cr_ref()
  vc <- sqrt(ref$met_fixed$varcomp)
  pin <- function(s) call("uniform", lower = 0.995 * s, upper = 1.005 * s)
  prior <- eval(bquote(fb_prior(
    sigma ~ .(pin(vc[["units!R"]])),
    sd(group = "gen:county") ~ .(pin(vc[["gen:county"]]))
  )))
  fit <- .cr_fit("met_fixed_pinned", yield ~ gen + county + county:rep,
                 ~gen:county, ref$data$met, prior = prior)
  tab <- predict(fit, classify = "gen")

  expect_identical(
    .cr_check(tab, ref$met_fixed$gen, "gen", point = 0.02,
              width = c(0.99, 1.01)),
    character(0)
  )
  expect_identical(attr(tab, "random"), "gen:county")
})

test_that("gen fixed, gen:county random agrees with ASReml by default", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  ref <- .cr_ref()
  fit <- .cr_fit("met_fixed", yield ~ gen + county + county:rep,
                 ~gen:county, ref$data$met)
  expect_identical(
    .cr_check(predict(fit, classify = "gen"), ref$met_fixed$gen, "gen"),
    character(0)
  )
})

test_that("a fixed classify with no random term keeps the emmeans table", {
  skip_if_not_installed("INLA")
  skip_if_not_installed("emmeans")
  skip_on_cran()
  ref <- .cr_ref()
  tab <- suppressMessages(predict(.cr_met(), classify = "county"))

  expect_null(attr(tab, "random"))
  expect_identical(
    .cr_check(tab, ref$met$county, "county", point = 0.01,
              width = c(0.95, 1.30)),
    character(0)
  )
})

test_that("the banner names the random terms, the averaged and the ignored", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  ref <- .cr_ref()
  fit <- .cr_fit("rcb_blk", yield ~ rep, ~ gen + rep:block, ref$data$rcb)
  out <- utils::capture.output(print(predict(fit, classify = "gen")))

  expect_match(out[[1L]], "Predicted values: gen", fixed = TRUE)
  expect_match(out[[2L]], "fixed part plus gen. Not an SED table.",
               fixed = TRUE)
  expect_identical(out[[3L]], "Averaged over: rep")
  expect_identical(out[[4L]], "Ignored: block")
})


# ---------------------------------------------------------------- #
# 3. brms against ASReml                                             #
# ---------------------------------------------------------------- #
#
# The INLA bounds plus three Monte Carlo errors: 2,000 draws put about 2%
# on an SD and 0.03 SE on a mean, so 0.10 SE is added to the point bound.

test_that("RCB and an empty MET cell agree with ASReml on brms", {
  skip_if_not_installed("brms")
  skip_on_cran()
  skip_on_ci()
  skip_if_not_long_tier()
  ref <- .cr_ref()
  args <- list(n_samples = 1000L, warmup = 1000L, chains = 2L, seed = 1L)

  rcb <- do.call(.cr_fit, c(
    list("rcb", yield ~ rep, ~gen, ref$data$rcb, backend = "brms"), args
  ))
  expect_identical(
    .cr_check(predict(rcb, classify = "gen"), ref$rcb$gen, "gen",
              point = 0.35, width = c(0.93, 1.32)),
    character(0)
  )

  gap <- do.call(.cr_fit, c(
    list("met_gap", yield ~ county + county:rep, ~ gen + gen:county,
         ref$data$met_gap, backend = "brms"),
    args
  ))
  expect_identical(
    .cr_check(predict(gap, classify = "gen:county"), ref$met_gap$gen_county,
              c("gen", "county"), point = 0.50, width = c(0.93, 1.32)),
    character(0)
  )
})


# ---------------------------------------------------------------- #
# 4. Refusals                                                        #
# ---------------------------------------------------------------- #

test_that("an argument of ASReml's predict() is refused by name", {
  obj <- structure(list(), class = "flexybayes")
  for (arg in c("present", "average", "ignore", "sed", "vcov", "levels")) {
    args <- c(list(obj, "gen"), stats::setNames(list(TRUE), arg))
    cnd <- tryCatch(
      do.call(flexyBayes:::.fb_predict_classify, args),
      error = function(e) e
    )
    expect_s3_class(cnd, "flexybayes_refusal_classify_argument_not_supported")
    expect_identical(cnd$argument, arg)
    expect_match(conditionMessage(cnd), paste0("`", arg, "`"), fixed = TRUE)
  }
})

test_that("the ASReml-argument refusal is reached through predict()", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  expect_error(
    predict(.cr_rcb(), classify = "gen", sed = TRUE),
    class = "flexybayes_refusal_classify_argument_not_supported"
  )
})

test_that("a non-identity link with a random classify term is refused", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  d <- .cr_ref()$data$rcb
  d$count <- round(10 * d$yield)
  fit <- .cr_fit("rcb_pois", count ~ rep, ~gen, d, family = "poisson")
  cnd <- tryCatch(predict(fit, classify = "gen"), error = function(e) e)

  expect_s3_class(
    cnd,
    "flexybayes_refusal_classify_random_factor_not_supported"
  )
  expect_identical(cnd$factor, "gen")
  expect_match(conditionMessage(cnd), "log link", fixed = TRUE)
})

test_that("a random factor with no usable term is refused", {
  # rep is not fixed here, so rep:block names a variable the rules ignore
  # and there is no term to predict block from.
  skip_if_not_installed("INLA")
  skip_on_cran()
  fit <- .cr_fit("blk_only", yield ~ 1, ~ gen + rep:block, .cr_ref()$data$rcb)
  cnd <- tryCatch(predict(fit, classify = "block"), error = function(e) e)

  expect_s3_class(
    cnd,
    "flexybayes_refusal_classify_random_factor_not_supported"
  )
  expect_identical(cnd$factor, "block")
  msg <- conditionMessage(cnd)
  expect_match(msg, "classify = \"rep:block\"", fixed = TRUE)
  expect_match(msg, "ranef(fit)", fixed = TRUE)
})

test_that("the classify refusal codes are registered", {
  entries <- ls(flexyBayes:::.refusal_registry)
  expect_true("classify_argument_not_supported" %in% entries)
  expect_true("classify_random_factor_not_supported" %in% entries)
})
