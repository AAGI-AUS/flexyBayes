# Tests for fb_prior() -- the PC-canonical hybrid prior DSL
# (deliverable 3).

# ---------------------------------------------------------------- #
# Constructor                                                      #
# ---------------------------------------------------------------- #

test_that("fb_prior() builds an fb_prior object from a single spec", {
  p <- fb_prior(sigma ~ pc(upper = 2, prob = 0.05))
  expect_s3_class(p, "fb_prior")
  expect_s3_class(p, "list")
  expect_length(p$specs, 1L)
  expect_identical(p$specs[[1]]$target$type, "sigma")
  expect_identical(p$specs[[1]]$spec$family, "pc")
  expect_identical(p$specs[[1]]$spec$args$upper, 2)
  expect_identical(p$specs[[1]]$spec$args$prob, 0.05)
})

test_that("fb_prior() accepts multiple specs", {
  p <- fb_prior(
    sigma ~ pc(upper = 2, prob = 0.05),
    sd(group = "subject") ~ half_normal(scale = 1),
    b("treatment") ~ student_t(df = 4, scale = 2.5)
  )
  expect_length(p$specs, 3L)
  types <- vapply(p$specs, function(s) s$target$type, character(1))
  expect_setequal(types, c("sigma", "sd", "b"))
})

test_that("fb_prior() requires at least one spec", {
  expect_error(fb_prior(), "at least one specification")
})

test_that("fb_prior() rejects non-formula arguments", {
  expect_error(fb_prior("not_a_formula"), "must be a two-sided formula")
  expect_error(fb_prior(~ pc(1, 0.01)), "must be two-sided")
})

# ---------------------------------------------------------------- #
# Targets                                                          #
# ---------------------------------------------------------------- #

test_that("fb_prior() parses sigma target", {
  p <- fb_prior(sigma ~ pc(upper = 1, prob = 0.01))
  expect_identical(p$specs[[1]]$target, list(type = "sigma"))
})

test_that("fb_prior() parses sd(group = ...) target", {
  p <- fb_prior(sd(group = "subject") ~ pc(upper = 1, prob = 0.01))
  expect_identical(p$specs[[1]]$target, list(type = "sd", group = "subject"))
})

test_that("fb_prior() parses b(...) target", {
  p <- fb_prior(b("treatment") ~ normal(mean = 0, sd = 5))
  expect_identical(p$specs[[1]]$target, list(type = "b", name = "treatment"))
})

test_that("fb_prior() parses cor(group = ...) target", {
  p <- fb_prior(cor(group = "subject") ~ lkj(eta = 2))
  expect_identical(p$specs[[1]]$target, list(type = "cor", group = "subject"))
})

test_that("fb_prior() parses smooth(var) target", {
  p <- fb_prior(smooth("time") ~ pc(upper = 1, prob = 0.01))
  expect_identical(p$specs[[1]]$target$type, "smooth")
  expect_identical(p$specs[[1]]$target$var, "time")
  expect_identical(p$specs[[1]]$target$basis, "rw2")
})

test_that("fb_prior() rejects unsupported targets", {
  expect_error(
    fb_prior(unsupported_target() ~ pc(1, 0.01)),
    "Unsupported prior target"
  )
  expect_error(fb_prior(sd() ~ pc(1, 0.01)), "sd\\(\\) prior target requires")
})

# ---------------------------------------------------------------- #
# Distributions                                                    #
# ---------------------------------------------------------------- #

test_that("fb_prior() accepts the v0.1 distribution families", {
  for (fam in c(
    "pc",
    "half_normal",
    "half_cauchy",
    "student_t",
    "normal",
    "exponential",
    "lkj",
    "cauchy",
    "gamma"
  )) {
    expr <- as.call(c(as.name(fam), 1, 0.01))
    f <- stats::as.formula(paste0("sigma ~ ", deparse(expr)))
    p <- fb_prior(f)
    expect_identical(
      p$specs[[1]]$spec$family,
      fam,
      info = paste("family =", fam)
    )
  }
})

test_that("fb_prior() rejects unsupported distributions", {
  expect_error(
    fb_prior(sigma ~ wishart(df = 3)),
    "Unsupported prior distribution: wishart"
  )
})

# ---------------------------------------------------------------- #
# Predicate + print                                                #
# ---------------------------------------------------------------- #

test_that("is_fb_prior() recognises the class", {
  p <- fb_prior(sigma ~ pc(upper = 1, prob = 0.01))
  expect_true(flexyBayes:::is_fb_prior(p))
  expect_false(flexyBayes:::is_fb_prior(list()))
  expect_false(flexyBayes:::is_fb_prior(NULL))
})

test_that("print.fb_prior() emits a multi-line summary", {
  p <- fb_prior(
    sigma ~ pc(upper = 2, prob = 0.05),
    sd(group = "subject") ~ half_normal(scale = 1)
  )
  out <- capture.output(print(p))
  expect_true(any(grepl("fb_prior", out)))
  expect_true(any(grepl("sigma", out)))
  expect_true(any(grepl("pc\\(upper = 2", out)))
  expect_true(any(grepl("subject", out)))
  expect_true(any(grepl("half_normal", out)))
})

# ---------------------------------------------------------------- #
# Translation helpers                                              #
# ---------------------------------------------------------------- #

test_that("priors_to_legacy() extracts vc_sd from sigma pc spec", {
  p <- fb_prior(sigma ~ pc(upper = 3, prob = 0.05))
  out <- flexyBayes:::priors_to_legacy(p)
  expect_identical(out$vc_sd, 3)
  expect_false(out$legacy)
})

test_that("priors_to_legacy() extracts fixed_sd from b normal spec", {
  p <- fb_prior(b("treatment") ~ normal(mean = 0, sd = 7))
  out <- flexyBayes:::priors_to_legacy(p)
  expect_identical(out$fixed_sd, 7)
})

test_that("priors_to_legacy() falls back to defaults for non-fb_prior input", {
  out <- flexyBayes:::priors_to_legacy(NULL)
  expect_identical(out$fixed_sd, 10)
  expect_identical(out$vc_sd, 1)
  expect_true(out$legacy)
})

test_that("priors_to_inla() emits pc.prec for sigma pc spec", {
  p <- fb_prior(sigma ~ pc(upper = 2, prob = 0.05))
  out <- flexyBayes:::priors_to_inla(p)
  expect_named(out, "sigma")
  expect_identical(out$sigma$prior, "pc.prec")
  expect_identical(out$sigma$param, c(2, 0.05))
})

test_that("priors_to_inla() emits pc.prec keyed by sd group", {
  p <- fb_prior(sd(group = "subject") ~ pc(upper = 1, prob = 0.01))
  out <- flexyBayes:::priors_to_inla(p)
  expect_named(out, "subject")
})

test_that("priors_to_inla() emits an exact half_normal expression prior", {
  p <- fb_prior(sigma ~ half_normal(scale = 1))
  out <- flexyBayes:::priors_to_inla(p)
  expect_named(out, "sigma")
  expect_match(out$sigma$prior, "^expression:")
  expect_null(out$sigma$param)
  expect_null(out$sigma$approx_from)
  expect_identical(out$sigma$meta$family, "half_normal")
  expect_identical(out$sigma$meta$scale, 1)
})

# ---------------------------------------------------------------- #
# uniform() DSL                                                    #
# ---------------------------------------------------------------- #

test_that("fb_prior() accepts uniform(lower, upper)", {
  p <- fb_prior(sd(group = "g") ~ uniform(lower = 0, upper = 5))
  expect_s3_class(p, "fb_prior")
  expect_identical(p$specs[[1]]$spec$family, "uniform")
  expect_identical(p$specs[[1]]$spec$args$lower, 0)
  expect_identical(p$specs[[1]]$spec$args$upper, 5)
})

test_that("uniform() rejects missing upper, negative lower, upper <= lower", {
  expect_error(
    fb_prior(sigma ~ uniform(lower = 0)),
    regexp = "uniform\\(\\) requires `upper`"
  )
  expect_error(
    fb_prior(sigma ~ uniform(lower = -1, upper = 5)),
    regexp = "lower must be >= 0"
  )
  expect_error(
    fb_prior(sigma ~ uniform(lower = 5, upper = 5)),
    regexp = "upper must be > lower"
  )
})

test_that("priors_to_inla() emits an exact uniform-on-SD expression prior", {
  local_clean_emit_state()
  p <- fb_prior(sigma ~ uniform(lower = 0, upper = 4))
  # No approximation message: the uniform is now represented exactly.
  expect_no_message(out <- flexyBayes:::priors_to_inla(p))
  expect_match(out$sigma$prior, "^expression:")
  expect_match(out$sigma$prior, "U=4")
  expect_null(out$sigma$param)
  expect_null(out$sigma$approx_from)
  expect_identical(out$sigma$meta$family, "uniform")
  expect_identical(out$sigma$meta$lower, 0)
  expect_identical(out$sigma$meta$upper, 4)
})

test_that(".default_uniform_prior() maps to faithful INLA expression priors", {
  # The synthesised default reaches all three INLA paths via
  # priors_to_inla(); it must carry the exact uniform-on-SD expression
  # prior (matching the greta backend's flat uniform), not the former
  # PC approximation that disagreed with greta on small group counts.
  dat <- data.frame(y = stats::rnorm(60), g = factor(rep(1:5, each = 12)))
  dp <- flexyBayes:::.default_uniform_prior(
    data = dat, response = "y", family = "gaussian", random_groups = "g"
  )
  out <- flexyBayes:::priors_to_inla(dp)
  expect_setequal(names(out), c("sigma", "g"))
  expect_match(out$sigma$prior, "^expression:")
  expect_match(out$g$prior, "^expression:")
  expect_false(any(vapply(
    out, function(e) identical(e$prior, "pc.prec"), logical(1)
  )))
})

test_that("priors_to_legacy() exposes uniform-on-VC via uniform_per_vc (ADR 0004)", {
  p <- fb_prior(sd(group = "g") ~ uniform(lower = 0, upper = 5))
  out <- flexyBayes:::priors_to_legacy(p)
  expect_true("g" %in% names(out$uniform_per_vc))
  expect_identical(out$uniform_per_vc$g$lower, 0)
  expect_identical(out$uniform_per_vc$g$upper, 5)
})

test_that("priors_to_legacy() exposes uniform-on-sigma via uniform_per_vc", {
  p <- fb_prior(sigma ~ uniform(lower = 0, upper = 7))
  out <- flexyBayes:::priors_to_legacy(p)
  expect_true("__sigma__" %in% names(out$uniform_per_vc))
  expect_identical(out$uniform_per_vc[["__sigma__"]]$upper, 7)
})

# ---------------------------------------------------------------- #
# v0.1.x default-prior change (PC default + soft deprecation)      #
# ---------------------------------------------------------------- #

test_that(".default_pc_prior() builds PC specs scaled by sd(y) (Gaussian)", {
  set.seed(11)
  d <- data.frame(y = rnorm(40, sd = 2), g = factor(rep(1:5, 8)))
  p <- flexyBayes:::.default_pc_prior(
    data = d,
    response = "y",
    family = "gaussian",
    random_groups = c("g"),
    alpha = 0.05
  )
  expect_s3_class(p, "fb_prior")
  expect_identical(p$specs[[1]]$target$type, "sigma")
  expect_identical(p$specs[[1]]$spec$family, "pc")
  # Scale is 2.5 * sd(y); allow tolerance for sample sd noise.
  expect_true(abs(p$specs[[1]]$spec$args$upper - 2.5 * sd(d$y)) < 1e-9)
  expect_identical(p$specs[[2]]$target$type, "sd")
  expect_identical(p$specs[[2]]$target$group, "g")
})

test_that(".default_prior_scale() picks family-aware scale", {
  set.seed(12)
  d_g <- data.frame(y = rnorm(40, sd = 3))
  s_g <- flexyBayes:::.default_prior_scale(d_g, "y", "gaussian")
  expect_identical(s_g$basis, "identity_sd")

  d_p <- data.frame(y = rpois(40, lambda = 5))
  s_p <- flexyBayes:::.default_prior_scale(d_p, "y", "poisson")
  expect_identical(s_p$basis, "log_link_sd")

  s_b <- flexyBayes:::.default_prior_scale(NULL, "y", "binomial")
  expect_identical(s_b$basis, "logit_default")
  expect_identical(s_b$scale, 2.5)
})

# ---------------------------------------------------------------- #
# v0.1 uniform default (ADR 0004 supersedes ADR 0003 PC default)   #
# ---------------------------------------------------------------- #

test_that(".default_uniform_prior() builds uniform specs scaled by 5*sd(y) (Gaussian)", {
  set.seed(11)
  d <- data.frame(y = rnorm(40, sd = 2), g = factor(rep(1:5, 8)))
  p <- flexyBayes:::.default_uniform_prior(
    data = d,
    response = "y",
    family = "gaussian",
    random_groups = c("g")
  )
  expect_s3_class(p, "fb_prior")
  expect_identical(p$specs[[1]]$target$type, "sigma")
  expect_identical(p$specs[[1]]$spec$family, "uniform")
  expect_identical(p$specs[[1]]$spec$args$lower, 0)
  expect_true(abs(p$specs[[1]]$spec$args$upper - 5 * sd(d$y)) < 1e-9)
  expect_identical(p$specs[[2]]$target$type, "sd")
  expect_identical(p$specs[[2]]$target$group, "g")
  expect_identical(p$specs[[2]]$spec$family, "uniform")
})

# AMBITION_STAGE.md §1.4: vm() + ped() structured-cov groups gain
# uniform-on-SD defaults (extending ADR 0004). Subtests cover the
# class of the returned object, the scale calculation, the
# attribute presence on the result, and that the other structured
# forms (at, us, fa, ar1) still fall through to the legacy
# lognormal default (unchanged in this work).

test_that(".default_uniform_prior() extends to vm() + ped() per §1.4", {
  set.seed(20260523L)
  d <- data.frame(y = stats::rnorm(50L, sd = 2))
  p <- flexyBayes:::.default_uniform_prior(
    data = d,
    response = "y",
    family = "gaussian",
    random_groups = character(0),
    vm_ped_groups = c("geno_vm", "animal_ped")
  )
  # Class.
  expect_s3_class(p, "fb_prior")
  # First spec is still sigma (residual SD).
  expect_identical(p$specs[[1]]$target$type, "sigma")
  # Two additional specs for vm + ped, each as sd(group = ...).
  expect_length(p$specs, 3L)
  expect_identical(p$specs[[2]]$target$type, "sd")
  expect_identical(p$specs[[2]]$target$group, "geno_vm")
  expect_identical(p$specs[[2]]$spec$family, "uniform")
  expect_identical(p$specs[[3]]$target$group, "animal_ped")
  expect_identical(p$specs[[3]]$spec$family, "uniform")
  # Scale calculation: identity-link Gaussian -> 5 * sd(y).
  expect_true(abs(p$specs[[2]]$spec$args$upper - 5 * stats::sd(d$y)) < 1e-9)
  # Attribute trail: per-form marker + per-prior basis.
  expect_identical(
    attr(p, "fb_prior_default_vm_ped_groups"),
    c("geno_vm", "animal_ped")
  )
  expect_identical(attr(p, "fb_prior_default_basis"), "identity_sd_uniform")
  expect_identical(attr(p$specs[[2]], "_default_uniform_form"), "vm_or_ped")
  expect_identical(attr(p$specs[[3]], "_default_uniform_form"), "vm_or_ped")
})

test_that(".default_uniform_prior() handles vm/ped alone (no simple RI)", {
  set.seed(20260524L)
  d <- data.frame(y = stats::rnorm(40L, sd = 3))
  p <- flexyBayes:::.default_uniform_prior(
    data = d,
    response = "y",
    family = "gaussian",
    vm_ped_groups = "geno_vm"
  )
  expect_length(p$specs, 2L) # sigma + one vm group
  expect_identical(p$specs[[2]]$target$group, "geno_vm")
})

test_that(".default_uniform_prior() leaves at/us/fa/ar1 outside its default set", {
  set.seed(20260525L)
  d <- data.frame(y = stats::rnorm(40L))
  # The function takes only `random_groups` (simple/ide/id) and
  # `vm_ped_groups` -- other structured forms are not part of the
  # caller's enumeration in flexybayes()/fb_brms() (see code in
  # the same files), so they continue to fall through to the
  # legacy lognormal default at codegen time. We assert here that
  # passing no vm_ped_groups returns only sigma + simple-RI specs.
  p <- flexyBayes:::.default_uniform_prior(
    data = d,
    response = "y",
    family = "gaussian",
    random_groups = "g",
    vm_ped_groups = character(0)
  )
  expect_length(p$specs, 2L) # sigma + g
  expect_identical(p$specs[[2]]$target$type, "sd")
  expect_identical(p$specs[[2]]$target$group, "g")
  # No vm_ped attribute when none supplied.
  expect_null(attr(p, "fb_prior_default_vm_ped_groups"))
})

test_that(".default_uniform_scale() picks family-aware bounds", {
  set.seed(12)
  d_g <- data.frame(y = rnorm(40, sd = 3))
  s_g <- flexyBayes:::.default_uniform_scale(d_g, "y", "gaussian")
  expect_identical(s_g$basis, "identity_sd_uniform")
  expect_true(abs(s_g$scale - 5 * sd(d_g$y)) < 1e-9)

  s_p <- flexyBayes:::.default_uniform_scale(NULL, "y", "poisson")
  expect_identical(s_p$basis, "log_link_uniform")
  expect_identical(s_p$scale, 3)

  s_b <- flexyBayes:::.default_uniform_scale(NULL, "y", "binomial")
  expect_identical(s_b$basis, "logit_uniform")
  expect_identical(s_b$scale, 5)

  s_be <- flexyBayes:::.default_uniform_scale(NULL, "y", "beta")
  expect_identical(s_be$basis, "logit_uniform")
  expect_identical(s_be$scale, 5)

  s_nb <- flexyBayes:::.default_uniform_scale(NULL, "y", "negative_binomial")
  expect_identical(s_nb$basis, "log_link_uniform")
  expect_identical(s_nb$scale, 3)
})
