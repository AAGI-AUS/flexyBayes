# =============================================================================
# The prior mini-language refuses malformed specifications, by name.
#
# The 0.9.1 field-coverage sweep submitted 37 constructions to `fb_prior()`:
# 31 malformed and 6 deliberate controls. Twenty-three of the 31 were
# accepted without error or warning, and three of those reached a completed
# fit under a silently substituted hyperparameter -- a prior-sensitivity
# analysis run through that path compares the default prior with itself.
#
# The corpus below is that sweep section, kept verbatim as a regression
# fixture. A mini-language is a public API: an argument name it does not
# have is a specification it cannot honour, and the refusal has to carry a
# class a caller can match on rather than a message a caller must parse.
#
# Each malformed row names the reason code it must raise. Each control row
# must still construct -- the corpus is a two-sided contract, and a
# constructor that refuses everything would pass a one-sided one.
# =============================================================================

# --- malformed corpus: every row refuses, with the named reason ---------- #

.corpus_malformed <- list(
  list(
    "wrong_arg_scale",
    sd(group = "g") ~ normal(scale = 2),
    "prior_argument_unknown"
  ),
  list(
    "wrong_arg_mu_sigma",
    sd(group = "g") ~ normal(mu = 0, sigma = 2),
    "prior_argument_unknown"
  ),
  list(
    "wrong_arg_hn_sd",
    sd(group = "g") ~ half_normal(sd = 1.5),
    "prior_argument_unknown"
  ),
  list(
    "extra_arg_hn",
    sd(group = "g") ~ half_normal(scale = 2, extra = 9),
    "prior_argument_unknown"
  ),
  list(
    "hc_extra_location",
    sigma ~ half_cauchy(location = 0, scale = 2),
    "prior_argument_unknown"
  ),
  list(
    "dup_named_arg",
    sigma ~ normal(mean = 0, sd = 2, sd = 3),
    "prior_argument_duplicated"
  ),
  list("pc_missing_prob", sigma ~ pc(upper = 1), "prior_argument_missing"),
  list("pc_missing_upper", sigma ~ pc(prob = 0.05), "prior_argument_missing"),
  list("t_single_positional", sigma ~ student_t(2.5), "prior_argument_missing"),
  list(
    "gamma_missing_rate",
    sigma ~ gamma(shape = 1),
    "prior_argument_missing"
  ),
  list(
    "negative_scale",
    sd(group = "g") ~ half_normal(scale = -1),
    "prior_hyperparameter_out_of_domain"
  ),
  list(
    "zero_scale",
    sd(group = "g") ~ half_normal(scale = 0),
    "prior_hyperparameter_out_of_domain"
  ),
  list(
    "pc_negative_upper",
    sigma ~ pc(upper = -1, prob = 0.05),
    "prior_hyperparameter_out_of_domain"
  ),
  list(
    "pc_prob_gt1",
    sigma ~ pc(upper = 1, prob = 1.5),
    "prior_hyperparameter_out_of_domain"
  ),
  list(
    "t_negative_df",
    sigma ~ student_t(df = -3, location = 0, scale = 2),
    "prior_hyperparameter_out_of_domain"
  ),
  list(
    "exp_negative_rate",
    sigma ~ exponential(rate = -1),
    "prior_hyperparameter_out_of_domain"
  ),
  list(
    "lkj_negative_eta",
    cor(group = "g") ~ lkj(eta = -1),
    "prior_hyperparameter_out_of_domain"
  ),
  list(
    "uniform_reversed",
    sigma ~ uniform(lower = 2, upper = 1),
    "prior_hyperparameter_out_of_domain"
  ),
  list(
    "uniform_negative_lower",
    sigma ~ uniform(lower = -1, upper = 5),
    "prior_hyperparameter_out_of_domain"
  ),
  list(
    "zerolen_scale",
    sd(group = "g") ~ half_normal(scale = numeric(0)),
    "prior_hyperparameter_not_scalar"
  ),
  list(
    "len2_scale",
    sd(group = "g") ~ half_normal(scale = c(1, 2)),
    "prior_hyperparameter_not_scalar"
  ),
  list(
    "na_scale",
    sd(group = "g") ~ half_normal(scale = NA),
    "prior_hyperparameter_not_scalar"
  ),
  list(
    "inf_scale",
    sd(group = "g") ~ half_normal(scale = Inf),
    "prior_hyperparameter_not_scalar"
  ),
  list(
    "char_scale",
    sd(group = "g") ~ half_normal(scale = "1"),
    "prior_hyperparameter_not_scalar"
  ),
  list(
    "bare_unknown_target",
    phi ~ half_normal(scale = 1),
    "prior_target_unsupported"
  ),
  list(
    "unknown_dist_normla",
    sd(group = "g") ~ normla(0, 2),
    "prior_distribution_unknown"
  ),
  list(
    "unknown_dist_cauchy2",
    sd(group = "g") ~ cauchy2(0, 2),
    "prior_distribution_unknown"
  ),
  list("target_not_call", sigma ~ half_normal, "prior_distribution_not_a_call"),
  list(
    "one_sided_formula",
    ~ half_normal(scale = 1),
    "prior_spec_not_two_sided"
  ),
  list(
    "sd_no_group",
    sd() ~ pc(upper = 1, prob = 0.01),
    "prior_target_argument_missing"
  )
)

test_that("every malformed prior specification refuses with its reason code", {
  for (case in .corpus_malformed) {
    label <- case[[1L]]
    cond <- tryCatch(
      {
        fb_prior(case[[2L]])
        NULL
      },
      condition = function(c) c
    )
    expect_false(
      is.null(cond),
      label = paste0(label, " was accepted (silent-wrong)")
    )
    expect_s3_class(cond, "flexybayes_refusal")
    expect_s3_class(cond, paste0("flexybayes_refusal_", case[[3L]]))
  }
})

test_that("the non-formula entry shapes refuse by name", {
  expect_s3_class(
    tryCatch(fb_prior(), condition = function(c) c),
    "flexybayes_refusal_prior_spec_empty"
  )
  expect_s3_class(
    tryCatch(
      fb_prior("sigma ~ half_normal(scale = 1)"),
      condition = function(c) c
    ),
    "flexybayes_refusal_prior_spec_not_formula"
  )
  expect_s3_class(
    tryCatch(
      fb_prior(list(sigma ~ half_normal(scale = 1))),
      condition = function(c) c
    ),
    "flexybayes_refusal_prior_spec_not_formula"
  )
})

# --- controls: the corpus is two-sided ----------------------------------- #

.corpus_controls <- list(
  list("positional_normal", sd(group = "g") ~ normal(0, 2)),
  list("normal_reversed_named", sd(group = "g") ~ normal(sd = 2, mean = 0)),
  list("positional_half_normal", sd(group = "g") ~ half_normal(3)),
  list("positional_pc", sd(group = "g") ~ pc(2, 0.05)),
  list("bad_b_target", b("nonexistent_term") ~ normal(mean = 0, sd = 2)),
  list(
    "bad_sd_group",
    sd(group = "nonexistent_group") ~ half_normal(scale = 1)
  ),
  list("sd_positional_group", sd("g") ~ half_normal(scale = 1)),
  list("b_unquoted_target", b(x) ~ normal(mean = 0, sd = 2)),
  list("uniform_upper_only", sd(group = "g") ~ uniform(upper = 5)),
  list("student_t_full", b("x") ~ student_t(df = 4, scale = 2.5)),
  list("cauchy_full", sigma ~ cauchy(location = 0, scale = 2)),
  list("gamma_full", sigma ~ gamma(shape = 1, rate = 2)),
  list("lkj_full", cor(group = "g") ~ lkj(eta = 2)),
  list("smooth_target", smooth("time") ~ pc(upper = 1, prob = 0.01))
)

test_that("the deliberate controls still construct", {
  for (case in .corpus_controls) {
    p <- fb_prior(case[[2L]])
    expect_s3_class(p, "fb_prior")
    expect_length(p$specs, 1L)
  }
})

test_that("a data-driven hyperparameter still evaluates in the caller frame", {
  y <- stats::rnorm(50L)
  p <- fb_prior(sigma ~ pc(upper = 2.5 * stats::sd(y), prob = 0.05))
  expect_equal(p$specs[[1]]$spec$args$upper, 2.5 * stats::sd(y))
})

# --- the parameter contract has one source ------------------------------- #

test_that("the family parameter table is the only list of prior families", {
  tbl <- flexyBayes:::.fb_prior_family_table()
  expect_setequal(
    names(tbl),
    flexyBayes:::.fb_prior_supported_families()
  )
  for (fam in names(tbl)) {
    entry <- tbl[[fam]]
    expect_true(all(entry$required %in% entry$params), info = fam)
    expect_setequal(names(entry$domain), entry$params)
    expect_true(
      all(entry$domain %in% c("positive", "nonnegative", "real", "unit_open")),
      info = fam
    )
  }
})


# =============================================================================
# The fit-entry half of the corpus (appended at 0.9.2).
#
# The rows above are constructions the DSL can judge on its own. These are
# the ones it cannot: whether a specification is honourable depends on the
# model (does that term exist?) and on the engine (can it carry that
# distribution on that parameter?), so they are legitimate at construction
# and refuse at the fit. Kept in the same file because a caller reading
# "which malformed priors does flexyBayes catch" should find one list.
#
# Sources: field-sweep FS-15 (the brms prior-emit layer was untyped),
# FS-20 (a prior on a term the model does not have reached brms's parser),
# FS-21 (an untranslatable INLA row was accepted, printed as applied, and
# dropped). Append-only: a row here is a regression pin.
# =============================================================================

.corpus_fit_entry_data <- function() {
  set.seed(20260819L)
  n <- 60L
  g <- factor(rep(seq_len(10L), each = 6L))
  x <- stats::rnorm(n)
  data.frame(y = 0.5 + 0.4 * x + stats::rnorm(n), x = x, g = g)
}

.corpus_fit_entry <- function(prior, backend) {
  options(flexyBayes.silence_default_prior_note = TRUE)
  args <- list(
    y ~ x,
    random = ~g,
    data = .corpus_fit_entry_data(),
    family = "gaussian",
    backend = backend,
    prior = prior,
    verbose = FALSE
  )
  if (identical(backend, "brms")) {
    args <- c(
      args,
      list(chains = 1L, n_samples = 200L, warmup = 100L, seed = 1L)
    )
  }
  tryCatch(
    suppressWarnings(suppressMessages(do.call(flexybayes, args))),
    error = function(e) e
  )
}

# name, spec, backend, reason code
.corpus_fit_entry_rows <- list(
  list(
    "bad_b_target_brms",
    b("nonexistent_term") ~ normal(mean = 0, sd = 2),
    "brms",
    "prior_target_not_in_model"
  ),
  list(
    "bad_sd_group_brms",
    sd(group = "nonexistent_group") ~ half_normal(scale = 1),
    "brms",
    "prior_target_not_in_model"
  ),
  list(
    "bad_b_target_inla",
    b("nonexistent_term") ~ normal(mean = 0, sd = 2),
    "inla",
    "prior_target_not_in_model"
  ),
  list(
    "bad_sd_group_inla",
    sd(group = "nonexistent_group") ~ half_normal(scale = 1),
    "inla",
    "prior_target_not_in_model"
  ),
  list(
    "lkj_on_sigma_brms",
    sigma ~ lkj(eta = 2),
    "brms",
    "prior_not_translatable_for_backend"
  ),
  list(
    "cor_target_brms",
    cor(group = "g") ~ lkj(eta = 2),
    "brms",
    "prior_not_translatable_for_backend"
  ),
  list(
    "smooth_target_brms",
    smooth("x", basis = "rw2") ~ pc(upper = 1, prob = 0.05),
    "brms",
    "prior_not_translatable_for_backend"
  ),
  list(
    "b_half_normal_brms",
    b("x") ~ half_normal(scale = 1),
    "brms",
    "prior_not_translatable_for_backend"
  ),
  list(
    "student_t_sd_inla",
    sd(group = "g") ~ student_t(df = 3, scale = 1),
    "inla",
    "prior_not_translatable_for_backend"
  ),
  list(
    "half_cauchy_sd_inla",
    sd(group = "g") ~ half_cauchy(scale = 1),
    "inla",
    "prior_not_translatable_for_backend"
  ),
  list(
    "gamma_sd_inla",
    sd(group = "g") ~ gamma(shape = 1, rate = 2),
    "inla",
    "prior_not_translatable_for_backend"
  ),
  list(
    "normal_sigma_inla",
    sigma ~ normal(mean = 0, sd = 2),
    "inla",
    "prior_not_translatable_for_backend"
  ),
  list(
    "cor_target_inla",
    cor(group = "g") ~ lkj(eta = 2),
    "inla",
    "prior_not_translatable_for_backend"
  ),
  list(
    "b_student_t_inla",
    b("x") ~ student_t(df = 3, scale = 1),
    "inla",
    "prior_not_translatable_for_backend"
  )
)

test_that("every fit-entry corpus row refuses typed, with the named reason", {
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  for (case in .corpus_fit_entry_rows) {
    p <- fb_prior(case[[2L]])
    expect_s3_class(p, "fb_prior") # legitimate at construction
    err <- .corpus_fit_entry(p, case[[3L]])
    expect_s3_class(err, paste0("flexybayes_refusal_", case[[4L]]))
    expect_s3_class(err, "flexybayes_refusal")
  }
})

# The controls for this half: specifications that used to refuse untyped on
# the brms emit and now translate, and the INLA rows that always did.

.corpus_fit_entry_controls <- list(
  list("positional_normal_brms", sd(group = "g") ~ normal(0, 2), "brms"),
  list(
    "reversed_named_brms",
    sd(group = "g") ~ normal(sd = 2, mean = 0),
    "brms"
  ),
  list("half_cauchy_sigma_brms", sigma ~ half_cauchy(scale = 2), "brms"),
  list("exponential_sd_brms", sd(group = "g") ~ exponential(rate = 2), "brms"),
  list("gamma_sd_brms", sd(group = "g") ~ gamma(shape = 1, rate = 2), "brms"),
  list("student_t_sigma_brms", sigma ~ student_t(df = 3, scale = 2), "brms"),
  list("exponential_sd_inla", sd(group = "g") ~ exponential(rate = 2), "inla"),
  list("b_normal_inla", b("x") ~ normal(mean = 0, sd = 4), "inla"),
  list(
    "half_normal_sd_inla",
    sd(group = "g") ~ half_normal(scale = 1.5),
    "inla"
  )
)

test_that("the fit-entry controls reach a completed fit", {
  skip_if_not_installed("brms")
  skip_if_not_installed("INLA")
  skip_on_cran()
  skip_on_ci()
  for (case in .corpus_fit_entry_controls) {
    fit <- .corpus_fit_entry(fb_prior(case[[2L]]), case[[3L]])
    expect_false(inherits(fit, "condition"), label = case[[1L]])
  }
})
