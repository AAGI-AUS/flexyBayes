# =============================================================================
# The per-(target, distribution) translation capability table.
#
# `priors_to_inla()` used to translate three distributions on three targets
# and return nothing at all for every other pair, while `prior_summary()`
# printed each of them as applied: five of the ten DSL distributions and two
# of the five targets were inert on INLA, silently (field-sweep FS-21). The
# brms side had the mirror-image defect -- real refusals, untyped, so no
# caller could tell "outside the translation table" from "R fell over"
# (FS-15).
#
# The tests below hold three things:
#
#   1. the table is total over the DSL -- every (engine, target,
#      distribution) triple has a verdict, so a distribution added without
#      one fails here rather than dropping silently;
#   2. every "translate" cell ARRIVES -- the hyper string is read back out
#      of the fitted INLA object, or the density out of the Stan program;
#   3. every "refuse" cell refuses TYPED at the fit, naming the row.
#
# The arrival half is the point. A translation test that asserts what
# `priors_to_inla()` returns proves the function agrees with itself; these
# read the engine.
# =============================================================================

.pt_data <- function(seed = 20260819L) {
  set.seed(seed)
  n <- 60L
  g <- factor(rep(seq_len(10L), each = 6L))
  x <- stats::rnorm(n)
  u <- stats::rnorm(10L, 0, 0.7)
  eta <- 0.5 + 0.4 * x + u[as.integer(g)]
  data.frame(y = eta + stats::rnorm(n), x = x, g = g)
}

.pt_silence <- function() {
  options(flexyBayes.silence_default_prior_note = TRUE)
}

.pt_brms_args <- list(
  chains = 1L,
  n_samples = 200L,
  warmup = 100L,
  seed = 1L,
  verbose = FALSE
)

# One fit at reachability budget, or the condition it raised.
.pt_fit <- function(prior, backend) {
  .pt_silence()
  d <- .pt_data()
  args <- list(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    backend = backend,
    prior = prior,
    verbose = FALSE
  )
  if (identical(backend, "brms")) {
    args <- c(args, .pt_brms_args[c("chains", "n_samples", "warmup", "seed")])
  }
  tryCatch(
    suppressWarnings(suppressMessages(do.call(flexybayes, args))),
    error = function(e) e
  )
}


# --- (1) the table is total over the DSL -------------------------------- #

test_that("every (engine, target, distribution) triple carries a verdict", {
  tbl <- flexyBayes:::.fb_prior_translation_table()
  families <- flexyBayes:::.fb_prior_supported_families()
  targets <- c("sigma", "sd", "smooth", "b", "cor")
  expect_setequal(unique(tbl$engine), c("inla", "brms"))
  expect_setequal(unique(tbl$target), targets)
  expect_setequal(unique(tbl$family), families)
  expect_true(all(tbl$verdict %in% c("translate", "refuse")))
  expect_equal(nrow(tbl), 2L * length(targets) * length(families))
  # No duplicated triple: the lookup returns the first hit, so a duplicate
  # would make the verdict depend on row order.
  key <- paste(tbl$engine, tbl$target, tbl$family, sep = "|")
  expect_false(any(duplicated(key)))
  # Every translate row names the engine argument it arrives through.
  expect_true(all(nzchar(tbl$route[tbl$verdict == "translate"])))
  expect_true(all(is.na(tbl$route[tbl$verdict == "refuse"])))
})

test_that("the rendered menu is derived from the table, not re-listed", {
  for (engine in c("inla", "brms")) {
    tbl <- flexyBayes:::.fb_prior_translation_table()
    tbl <- tbl[tbl$engine == engine & tbl$verdict == "translate", ]
    menu <- flexyBayes:::.fb_prior_translation_menu(engine)
    for (i in seq_len(nrow(tbl))) {
      expect_match(menu, tbl$family[[i]], fixed = TRUE)
    }
    # A refused family must not appear as if it were on offer.
    refused <- setdiff(
      flexyBayes:::.fb_prior_supported_families(),
      unique(tbl$family)
    )
    for (fam in refused) {
      expect_false(grepl(fam, menu, fixed = TRUE), label = paste(engine, fam))
    }
  }
})


# --- (2) arrival: INLA -------------------------------------------------- #

test_that("every translatable INLA variance-component cell reaches the hyper string", {
  skip_if_not_installed("INLA")
  .pt_silence()
  cases <- list(
    list(
      "half_normal",
      fb_prior(sd(group = "g") ~ half_normal(scale = 1.5)),
      "s=1.5"
    ),
    list(
      "uniform",
      fb_prior(sd(group = "g") ~ uniform(lower = 0, upper = 5)),
      "5"
    ),
    list(
      "pc",
      fb_prior(sd(group = "g") ~ pc(upper = 2, prob = 0.05)),
      "pc.prec"
    ),
    list(
      "exponential",
      fb_prior(sd(group = "g") ~ exponential(rate = 2)),
      "pc.prec"
    )
  )
  for (case in cases) {
    fit <- .pt_fit(case[[2L]], "inla")
    expect_false(inherits(fit, "condition"), label = case[[1L]])
    form <- paste(deparse(fit$inla$.args$formula), collapse = "")
    expect_match(form, "hyper", label = case[[1L]])
    expect_match(form, case[[3L]], fixed = TRUE, label = case[[1L]])
  }
})

test_that("the sigma cell reaches control.family rather than the formula", {
  skip_if_not_installed("INLA")
  fit <- .pt_fit(fb_prior(sigma ~ exponential(rate = 3)), "inla")
  expect_false(inherits(fit, "condition"))
  # INLA rewrites the supplied `hyper` list into its own theta-indexed
  # slots, so the entry is found by the short name the likelihood declares
  # (R/family_traits.R) rather than by the key we wrote.
  # INLA stamps an `inla.read.only` attribute on each field, so the
  # comparisons strip attributes rather than using identical().
  hyper <- fit$inla$.args$control.family[[1L]]$hyper
  slot <- Filter(function(h) isTRUE(h$short.name == "prec"), hyper)[[1L]]
  expect_true(grepl("pc.prec", slot$prior, fixed = TRUE))
  # -log(alpha) / U must recover the requested rate exactly.
  expect_equal(-log(slot$param[[2L]]) / slot$param[[1L]], 3)
})

test_that("exponential on an SD target IS the PC prior, checked numerically", {
  # The identity that makes the exponential row translatable: pc.prec(U,
  # alpha) states P(sigma > U) = alpha, i.e. Exp(-log(alpha) / U) on the
  # standard-deviation scale. Grounded against the definition rather than
  # against the emit.
  tr <- flexyBayes:::priors_to_inla(
    fb_prior(sd(group = "g") ~ exponential(rate = 2.5))
  )
  param <- tr[["g"]]$param
  rate <- -log(param[[2L]]) / param[[1L]]
  expect_equal(rate, 2.5)
  # And the tail statement the parameters encode is the exponential's own.
  expect_equal(stats::pexp(param[[1L]], rate, lower.tail = FALSE), param[[2L]])
})

test_that("a b() prior reaches INLA's control.fixed by coefficient name", {
  skip_if_not_installed("INLA")
  fit <- .pt_fit(fb_prior(b("x") ~ normal(mean = 1, sd = 4)), "inla")
  expect_false(inherits(fit, "condition"))
  cf <- fit$inla$.args$control.fixed
  expect_equal(cf$mean$x, 1)
  expect_equal(cf$prec$x, 1 / 16)
})

test_that("the b() route composes with the prior_fixed_sd scalar", {
  skip_if_not_installed("INLA")
  .pt_silence()
  d <- .pt_data()
  fit <- suppressMessages(flexybayes(
    y ~ x,
    random = ~g,
    data = d,
    family = "gaussian",
    backend = "inla",
    prior_fixed_sd = 7,
    prior = fb_prior(b("x") ~ normal(mean = 0, sd = 2)),
    verbose = FALSE
  ))
  cf <- fit$inla$.args$control.fixed
  # The named entry wins for `x`; the scalar is the default for the rest.
  expect_equal(cf$prec$x, 1 / 4)
  expect_equal(cf$prec$default, 1 / 49)
})

test_that("the six full-specification INLA cells stay green", {
  # The regression pin. A blanket refusal of untranslatable rows would have
  # turned these red, because the documented full specification carries a
  # b() row -- which is why the fix is a table and not a refusal.
  skip_if_not_installed("INLA")
  .pt_silence()
  full <- fb_prior(
    sigma ~ half_normal(scale = 2),
    sd(group = "g") ~ half_normal(scale = 1.5),
    b("x") ~ normal(mean = 0, sd = 4)
  )
  for (fam in c(
    "gaussian",
    "poisson",
    "negative_binomial",
    "binomial",
    "gamma",
    "beta"
  )) {
    d <- .pt_data()
    d$y <- switch(
      fam,
      "poisson" = ,
      "negative_binomial" = stats::rpois(nrow(d), 2),
      "binomial" = stats::rbinom(nrow(d), 1L, 0.5),
      "gamma" = stats::rgamma(nrow(d), shape = 2, rate = 1),
      "beta" = stats::rbeta(nrow(d), 2, 2),
      d$y
    )
    fit <- tryCatch(
      suppressWarnings(suppressMessages(flexybayes(
        y ~ x,
        random = ~g,
        data = d,
        family = fam,
        backend = "inla",
        prior = full,
        verbose = FALSE
      ))),
      error = function(e) e
    )
    expect_false(inherits(fit, "condition"), label = fam)
    expect_equal(fit$inla$.args$control.fixed$prec$x, 1 / 16, label = fam)
  }
})


# --- (2) arrival: brms -------------------------------------------------- #

test_that("every translatable brms variance-component cell reaches the Stan program", {
  skip_if_not_installed("brms")
  .pt_silence()
  d <- .pt_data()
  cases <- list(
    list(
      quote(sd(group = "g") ~ normal(mean = 0, sd = 2)),
      "normal_lpdf\\(sd_1 \\| 0, 2\\)"
    ),
    list(
      quote(sd(group = "g") ~ half_normal(scale = 2)),
      "normal_lpdf\\(sd_1 \\| 0, 2\\)"
    ),
    list(
      quote(sd(group = "g") ~ half_cauchy(scale = 2)),
      "cauchy_lpdf\\(sd_1 \\| 0, 2\\)"
    ),
    list(
      quote(sd(group = "g") ~ cauchy(location = 1, scale = 2)),
      "cauchy_lpdf\\(sd_1 \\| 1, 2\\)"
    ),
    list(
      quote(sd(group = "g") ~ student_t(df = 3, scale = 2)),
      "student_t_lpdf\\(sd_1 \\| 3, 0, 2\\)"
    ),
    list(
      quote(sd(group = "g") ~ exponential(rate = 2)),
      "exponential_lpdf\\(sd_1 \\| 2\\)"
    ),
    list(
      quote(sd(group = "g") ~ gamma(shape = 1, rate = 2)),
      "gamma_lpdf\\(sd_1 \\| 1, 2\\)"
    ),
    list(
      quote(sd(group = "g") ~ uniform(lower = 0, upper = 5)),
      "uniform_lpdf\\(sd_1 \\| 0, 5\\)"
    ),
    list(
      quote(sigma ~ pc(upper = 1, prob = 0.05)),
      "exponential_lpdf\\(sigma \\| 2.99"
    )
  )
  for (case in cases) {
    p <- eval(bquote(fb_prior(.(case[[1L]]))))
    code <- suppressMessages(flexybayes(
      y ~ x,
      random = ~g,
      data = d,
      family = "gaussian",
      backend = "brms",
      prior = p,
      return_code = TRUE,
      verbose = FALSE
    ))
    expect_match(
      paste(code, collapse = "\n"),
      case[[2L]],
      label = deparse(case[[1L]])
    )
  }
})

test_that("a variance-component row is bounded below at zero", {
  # Which is what makes `normal(0, s)` and `half_normal(s)` the same prior
  # on an SD: brms renormalises over the truncated support.
  sp <- flexyBayes:::.brms_vc_density(
    list(family = "normal", args = list(mean = 0, sd = 2))
  )
  expect_identical(sp$string, "normal(0, 2)")
  expect_identical(sp$lb, 0)
  hn <- flexyBayes:::.brms_vc_density(
    list(family = "half_normal", args = list(scale = 2))
  )
  expect_identical(hn$string, sp$string)
  expect_identical(hn$lb, sp$lb)
})


# --- (3) refusal: every "refuse" cell refuses typed ---------------------- #

test_that("an untranslatable INLA row refuses typed at the fit, naming the row", {
  skip_if_not_installed("INLA")
  cases <- list(
    fb_prior(sd(group = "g") ~ student_t(df = 3, scale = 1)),
    fb_prior(sd(group = "g") ~ half_cauchy(scale = 1)),
    fb_prior(sd(group = "g") ~ gamma(shape = 1, rate = 2)),
    fb_prior(sigma ~ normal(mean = 0, sd = 2)),
    fb_prior(cor(group = "g") ~ lkj(eta = 2)),
    fb_prior(b("x") ~ student_t(df = 3, scale = 1))
  )
  for (p in cases) {
    fit <- .pt_fit(p, "inla")
    expect_s3_class(
      fit,
      "flexybayes_refusal_prior_not_translatable_for_backend"
    )
    msg <- conditionMessage(fit)
    expect_match(msg, p$specs[[1L]]$spec$family, fixed = TRUE)
    # Both remedies, named.
    expect_match(msg, "re-express")
    expect_match(msg, "backend = \"brms\"", fixed = TRUE)
  }
})

test_that("an untranslatable brms row refuses typed at the fit, naming the row", {
  skip_if_not_installed("brms")
  cases <- list(
    fb_prior(sigma ~ lkj(eta = 2)),
    fb_prior(b("x") ~ half_normal(scale = 1)),
    fb_prior(b("x") ~ pc(upper = 1, prob = 0.05)),
    fb_prior(cor(group = "g") ~ lkj(eta = 2)),
    fb_prior(smooth("x", basis = "rw2") ~ pc(upper = 1, prob = 0.05))
  )
  for (p in cases) {
    fit <- .pt_fit(p, "brms")
    expect_s3_class(
      fit,
      "flexybayes_refusal_prior_not_translatable_for_backend"
    )
    expect_match(conditionMessage(fit), "brms", fixed = TRUE)
  }
})

test_that("the refusal is registered and discoverable", {
  reg <- fb_refusals(reason_code = "prior_not_translatable_for_backend")
  expect_equal(nrow(reg), 1L)
  expect_equal(reg$since_version, "0.9.2")
})

test_that("no prior specification is accepted, printed, and then dropped", {
  # The FS-21 contract in one line: whatever prior_summary() shows on an
  # INLA fit reached the engine, because anything that could not have is
  # refused before the fit exists.
  skip_if_not_installed("INLA")
  fit <- .pt_fit(
    fb_prior(
      sd(group = "g") ~ half_normal(scale = 1.5),
      b("x") ~ normal(mean = 0, sd = 4)
    ),
    "inla"
  )
  expect_false(inherits(fit, "condition"))
  ps <- prior_summary(fit)
  form <- paste(deparse(fit$inla$.args$formula), collapse = "")
  cf <- fit$inla$.args$control.fixed
  for (s in ps$fb_prior$specs) {
    arrived <- switch(
      s$target$type,
      "sd" = grepl("hyper", form, fixed = TRUE),
      "b" = !is.null(cf$mean[[s$target$name]]),
      TRUE
    )
    expect_true(arrived, label = flexyBayes:::.format_prior_target(s$target))
  }
})


# --- (4) target existence (field-sweep FS-20) ---------------------------- #

test_that("a prior naming a term the model does not have refuses typed on brms", {
  skip_if_not_installed("brms")
  fit <- .pt_fit(
    fb_prior(b("nonexistent_term") ~ normal(mean = 0, sd = 2)),
    "brms"
  )
  expect_s3_class(fit, "flexybayes_refusal_prior_target_not_in_model")
  expect_match(conditionMessage(fit), "nonexistent_term", fixed = TRUE)
  # The model's actual terms, in the user's vocabulary -- not brms's
  # synthesised `b_nonexistent_term`.
  expect_match(conditionMessage(fit), "\"x\"", fixed = TRUE)
  expect_false(grepl("b_nonexistent_term", conditionMessage(fit), fixed = TRUE))
})

test_that("a prior naming a group the model does not have refuses typed on brms", {
  skip_if_not_installed("brms")
  fit <- .pt_fit(
    fb_prior(sd(group = "nonexistent_group") ~ half_normal(scale = 1)),
    "brms"
  )
  expect_s3_class(fit, "flexybayes_refusal_prior_target_not_in_model")
  expect_match(conditionMessage(fit), "nonexistent_group", fixed = TRUE)
  expect_match(conditionMessage(fit), "\"g\"", fixed = TRUE)
})

test_that("the same two mistakes refuse typed on INLA, where they were ignored", {
  skip_if_not_installed("INLA")
  for (p in list(
    fb_prior(b("nonexistent_term") ~ normal(mean = 0, sd = 2)),
    fb_prior(sd(group = "nonexistent_group") ~ half_normal(scale = 1))
  )) {
    fit <- .pt_fit(p, "inla")
    expect_s3_class(fit, "flexybayes_refusal_prior_target_not_in_model")
    expect_match(conditionMessage(fit), "nonexistent")
  }
})

test_that("a legitimate target on either engine is not caught by the check", {
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  p <- fb_prior(
    b("x") ~ normal(mean = 0, sd = 4),
    sd(group = "g") ~ half_normal(scale = 1.5)
  )
  expect_false(inherits(.pt_fit(p, "brms"), "condition"))
  expect_false(inherits(.pt_fit(p, "inla"), "condition"))
})

test_that("prior_summary() does not claim the engine default for a b() row", {
  # The one-story contract in the other direction. A `b()` prior now
  # reaches the engine, so the engine-default sentence has to except the
  # coefficients it names.
  skip_if_not_installed("INLA")
  fit <- .pt_fit(fb_prior(b("x") ~ normal(mean = 0, sd = 4)), "inla")
  expect_false(inherits(fit, "condition"))
  txt <- paste(
    utils::capture.output(print(prior_summary(fit))),
    collapse = " "
  )
  expect_match(txt, "except `x`", fixed = TRUE)
  fit2 <- .pt_fit(fb_prior(sd(group = "g") ~ half_normal(scale = 1.5)), "inla")
  txt2 <- paste(
    utils::capture.output(print(prior_summary(fit2))),
    collapse = " "
  )
  expect_match(txt2, "the fixed effects carry", fixed = TRUE)
  expect_false(grepl("except", txt2, fixed = TRUE))
})
