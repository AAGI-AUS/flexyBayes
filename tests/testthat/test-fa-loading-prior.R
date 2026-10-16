# The prior on the loadings of fa(f, k):g on INLA.
#
# Each loading is the scaling of a copy of the anchor field, so it is the
# ratio of a level's loading to the anchor's. INLA's own prior on a copy
# scaling is normal(mean 1, precision 10), an SD of 0.32 on that ratio. On
# agridat::dasilva.maize it moved the common genetic variance at the posterior
# mode by up to 0.72 relative to ASReml 4.2 (DEF-229). The prior is written as
# `loading(group = <inner>) ~ normal(mean, sd)`, sent on every copy, recorded
# in the model manifest and read back from INLA's own record (`all.hyper`).

.flp_fixture <- function(n_gen = 40L, n_rep = 3L, scale = 1) {
  set.seed(51L)
  vec_l <- c(1.0, 1.6, 0.5, -0.9)
  vec_f <- rnorm(n_gen, 0, 0.9)
  mat_u <- outer(vec_f, vec_l) + matrix(rnorm(n_gen * 4L, 0, 0.25), ncol = 4L)
  tab <- expand.grid(
    gen = factor(seq_len(n_gen)),
    env = factor(seq_len(4L)),
    rep = seq_len(n_rep)
  )
  tab$y <- scale *
    (5 +
      mat_u[cbind(as.integer(tab$gen), as.integer(tab$env))] +
      rnorm(nrow(tab), 0, 0.3))
  tab
}

.flp_copy_text <- function(tab, prior = NULL) {
  fb <- fb_from_asreml(
    fixed = y ~ env,
    random = ~ fa(env, 1):gen,
    data = tab,
    family = "gaussian",
    prior = prior
  )
  code <- suppressMessages(emit_inla(
    fb,
    tab,
    return_code = TRUE,
    verbose = FALSE
  ))
  txt <- paste(deparse(code$formula, width.cutoff = 500L), collapse = " ")
  gsub("\\s+", " ", txt)
}

test_that("loading(group = ) ~ normal() is a prior target", {
  p <- fb_prior(loading(group = "gen") ~ normal(mean = 1, sd = 1))
  expect_identical(p$specs[[1L]]$target, list(type = "loading", group = "gen"))
  entry <- priors_to_inla(p)[["loading:gen"]]
  expect_identical(entry$prior, "normal")
  expect_equal(entry$param, c(1, 1))
})

test_that("a loading prior other than normal() is refused on INLA by name", {
  skip_if_not_installed("INLA")
  tab <- .flp_fixture(n_gen = 10L, n_rep = 2L)
  p <- fb_prior(loading(group = "gen") ~ half_normal(scale = 1))
  expect_error(
    .flp_copy_text(tab, p),
    class = "flexybayes_refusal_prior_not_translatable_for_backend"
  )
})

test_that("every copy carries the loading prior, the default when none is stated", {
  skip_if_not_installed("INLA")
  tab <- .flp_fixture(n_gen = 10L, n_rep = 2L)
  txt <- .flp_copy_text(tab)
  n_copy <- lengths(regmatches(txt, gregexpr("copy = ", txt, fixed = TRUE)))
  expect_identical(n_copy, 3L)
  n_prior <- lengths(regmatches(
    txt,
    gregexpr(
      "hyper = list(beta = list(prior = \"normal\", param = c(1, 10)))",
      txt,
      fixed = TRUE
    )
  ))
  expect_identical(n_prior, 3L)
  # A stated prior replaces it, as mean and precision.
  txt <- .flp_copy_text(tab, fb_prior(loading(group = "gen") ~ normal(1, 1)))
  expect_match(txt, "param = c(1, 1)))", fixed = TRUE)
  expect_false(grepl("param = c(1, 10)", txt, fixed = TRUE))
})

test_that("the loading prior does not move with the response scale", {
  skip_if_not_installed("INLA")
  # A loading is a ratio of two standard deviations, so the standardised
  # fit carries the same prior on y and on 1000 y.
  a <- .flp_copy_text(.flp_fixture(n_gen = 10L, n_rep = 2L))
  b <- .flp_copy_text(.flp_fixture(n_gen = 10L, n_rep = 2L, scale = 1000))
  pat <- "param = c(1, 10)"
  expect_identical(
    lengths(regmatches(b, gregexpr(pat, b, fixed = TRUE))),
    lengths(regmatches(a, gregexpr(pat, a, fixed = TRUE)))
  )
})

test_that("DEF-229: the loading prior is in the manifest and read back from INLA", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.inla_threads = "1:1",
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE
  )
  fit <- .fa_fit_past_refusal(y ~ env, ~ fa(env, 1):gen, .flp_fixture())
  man <- fit$extras$manifest
  inp <- man$inputs[man$inputs$input_id == "loading:gen", , drop = FALSE]
  expect_identical(nrow(inp), 1L)
  expect_identical(inp$source, "filled_default")
  expect_identical(inp$distribution, "normal")
  sl <- man$slots[man$slots$input_id %in% "loading:gen", , drop = FALSE]
  # One slot per copy: three levels after the anchor.
  expect_identical(nrow(sl), 3L)
  expect_true(all(sl$readback_agrees))
  expect_true(all(sl$source == "package"))
  expect_true(all(grepl("normal(1, 10)", sl$readback, fixed = TRUE)))
  txt <- paste(
    utils::capture.output(print(prior_summary(fit))),
    collapse = "\n"
  )
  expect_match(txt, "loading_gen ~ normal(mean=1, sd=0.316", fixed = TRUE)
  expect_match(txt, "-> random:gen_fa1_env_1:theta", fixed = TRUE)
  expect_false(grepl("flexyBayes sets no prior on them", txt, fixed = TRUE))
})

test_that("a stated loading prior is read back as stated", {
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.inla_threads = "1:1",
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE
  )
  testthat::local_mocked_bindings(
    .refuse_unrepresentable_structures = function(fb) invisible(NULL)
  )
  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~ fa(env, 1):gen,
    data = .flp_fixture(),
    prior = fb_prior(loading(group = "gen") ~ normal(mean = 1, sd = 1)),
    backend = "inla",
    verbose = FALSE
  ))
  man <- fit$extras$manifest
  expect_identical(
    man$inputs$source[man$inputs$input_id == "loading:gen"],
    "explicit"
  )
  sl <- man$slots[man$slots$input_id %in% "loading:gen", , drop = FALSE]
  expect_identical(nrow(sl), 3L)
  expect_true(all(sl$readback_agrees))
  expect_true(all(grepl("normal(1, 1)", sl$readback, fixed = TRUE)))
})

test_that("the read-back catches a loading prior that differs from the one recorded", {
  # Planted defect: the formula INLA receives carries precision 2 while the
  # manifest records the precision 10 the package meant to send.
  skip_if_not_installed("INLA")
  skip_on_cran()
  withr::local_options(
    flexyBayes.inla_threads = "1:1",
    flexyBayes.silence_default_prior_note = TRUE,
    flexyBayes.silence_weak_identification_warning = TRUE
  )
  real <- .build_inla_formula
  plant <- function(f) {
    txt <- paste(deparse(f, width.cutoff = 500L), collapse = " ")
    txt <- gsub("param = c(1, 10)", "param = c(1, 2)", txt, fixed = TRUE)
    stats::as.formula(txt, env = environment(f))
  }
  testthat::local_mocked_bindings(
    .refuse_unrepresentable_structures = function(fb) invisible(NULL),
    .build_inla_formula = function(...) {
      out <- real(...)
      if (inherits(out, "formula")) {
        return(plant(out))
      }
      out$formula <- plant(out$formula)
      out
    }
  )
  fit <- suppressMessages(flexybayes(
    y ~ env,
    random = ~ fa(env, 1):gen,
    data = .flp_fixture(),
    backend = "inla",
    verbose = FALSE
  ))
  sl <- fit$extras$manifest$slots
  sl <- sl[sl$input_id %in% "loading:gen", , drop = FALSE]
  expect_identical(nrow(sl), 3L)
  expect_false(any(sl$readback_agrees))
  expect_true(all(grepl("normal(1, 2)", sl$readback, fixed = TRUE)))
})
