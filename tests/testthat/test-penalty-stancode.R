# =============================================================================
# Generated Stan code: what each compiled family actually does to which
# coefficient.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 4.4 and 12
# (the F0 gate), Unit A.
#
# These are semantic assertions, not text snapshots. A snapshot of a brms
# program goes red when brms reformats a comment and stays green when the
# index moves, which is precisely backwards: an off-by-one in `cand_idx`
# shrinks the wrong coefficient and the generated program is perfectly valid
# Stan either way. So every assertion here is about a quantity with a
# meaning -- which `b[i]` receives which density, how many contributions each
# candidate gets, which latent parameters exist and at what dimension, and
# whether the protected coefficient is touched at all.
#
# The index mapping is checked against `brms::make_standata()` rather than
# against the convention written in the compiler's header. The `b` vector
# excludes the intercept under the centred parameterisation, and the only
# authority on that is brms.
# =============================================================================


# --- fixtures --------------------------------------------------------------- #

fbps_data <- function(n = 60L, seed = 5150L) {
  set.seed(seed)
  data.frame(
    y = stats::rnorm(n),
    p1 = stats::rnorm(n),
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n),
    x3 = stats::rnorm(n),
    g = factor(rep(seq_len(6L), length.out = n)),
    stringsAsFactors = FALSE
  )
}

fbps_formula <- function() y ~ p1 + x1 + x2 + x3 + (1 | g)

fbps_map <- function(d = fbps_data()) {
  .fb_sel_term_map(~ p1 + x1 + x2 + x3, d)
}

# The b vector of the extended model, intercept excluded, in order. `p1` is
# the protected coefficient and sits at b[1].
fbps_coef_all <- c("p1", "x1", "x2", "x3")

fbps_lines <- function(code) strsplit(code, "\n", fixed = TRUE)[[1L]]

# Every `b[<i>]` that appears inside a log-density contribution, with
# repeats, so "exactly one contribution per coefficient" is a countable
# claim rather than a set membership one.
fbps_prior_b_indices <- function(code, density = NULL) {
  ln <- fbps_lines(code)
  ln <- grep("lprior \\+=|target \\+=", ln, value = TRUE)
  if (!is.null(density)) {
    ln <- grep(paste0(density, "_lpdf\\("), ln, value = TRUE)
  }
  hits <- regmatches(ln, gregexpr("b\\[[0-9]+\\]", ln))
  idx <- as.integer(gsub("[^0-9]", "", unlist(hits, use.names = FALSE)))
  sort(idx)
}

fbps_declares <- function(code, pattern) {
  any(grepl(pattern, fbps_lines(code)))
}


# --- coefficient-row families ------------------------------------------------ #

test_that("ridge puts one normal on each candidate index and none on the protected one", {
  skip_if_not_installed("brms")
  d <- fbps_data()
  blocks <- .fb_pen_resolve(
    fb_penalty(term(everything()) ~ ridge(scale = 0.4)),
    fbps_map(d),
    protected = "p1"
  )
  out <- .fb_pen_compile_brms(blocks, fbps_coef_all, n_train = nrow(d))
  code <- brms::make_stancode(fbps_formula(), d, prior = out$prior)

  expect_identical(fbps_prior_b_indices(code, "normal"), c(2L, 3L, 4L))
  # Exactly one contribution per candidate, and none of any kind on b[1].
  expect_identical(fbps_prior_b_indices(code), c(2L, 3L, 4L))
  expect_true(all(grepl(
    "normal_lpdf\\(b\\[[234]\\] \\| 0, 0.4\\)",
    grep("normal_lpdf\\(b\\[", fbps_lines(code), value = TRUE)
  )))

  # The b index is brms's, and brms is the authority on it.
  sd <- brms::make_standata(fbps_formula(), d, prior = out$prior)
  expect_equal(sd$Kc, 4)
  expect_identical(dimnames(sd$X)[[2L]], c("Intercept", fbps_coef_all))
  expect_identical(out$param_map$b_index, 2:4)
  for (k in seq_len(nrow(out$param_map))) {
    expect_equal(
      unname(sd$X[, out$param_map$b_index[k] + 1L]),
      unname(d[[out$param_map$coefficient[k]]])
    )
  }
})


test_that("laplace puts one double_exponential on each candidate index", {
  skip_if_not_installed("brms")
  d <- fbps_data()
  blocks <- .fb_pen_resolve(
    fb_penalty(term(x1, x2, x3) ~ laplace(scale = 0.25)),
    fbps_map(d),
    protected = "p1"
  )
  out <- .fb_pen_compile_brms(blocks, fbps_coef_all, n_train = nrow(d))
  code <- brms::make_stancode(fbps_formula(), d, prior = out$prior)

  expect_identical(
    fbps_prior_b_indices(code, "double_exponential"), c(2L, 3L, 4L)
  )
  expect_identical(fbps_prior_b_indices(code), c(2L, 3L, 4L))
  expect_true(fbps_declares(code, "double_exponential_lpdf\\(b\\[2\\] \\| 0, 0.25\\)"))
})


test_that("two coefficient-row families partition the candidate indices", {
  skip_if_not_installed("brms")
  d <- fbps_data()
  blocks <- .fb_pen_resolve(
    fb_penalty(
      term(everything()) ~ ridge(scale = 1),
      term(x3) ~ laplace(scale = 0.1)
    ),
    fbps_map(d),
    protected = "p1"
  )
  out <- .fb_pen_compile_brms(blocks, fbps_coef_all, n_train = nrow(d))
  code <- brms::make_stancode(fbps_formula(), d, prior = out$prior)

  expect_identical(fbps_prior_b_indices(code, "normal"), c(2L, 3L))
  expect_identical(fbps_prior_b_indices(code, "double_exponential"), 4L)
  # Still exactly one contribution each, and still nothing on the protected
  # coefficient: the union is the partition, not an overlap.
  expect_identical(fbps_prior_b_indices(code), c(2L, 3L, 4L))
})


test_that("a protected coefficient keeps its own fb_prior row alongside a block", {
  skip_if_not_installed("brms")
  d <- fbps_data()
  blocks <- .fb_pen_resolve(
    fb_penalty(term(everything()) ~ laplace(scale = 0.25)),
    fbps_map(d),
    protected = "p1"
  )
  out <- .fb_pen_compile_brms(blocks, fbps_coef_all, n_train = nrow(d))
  protected_row <- .brms_specs_to_object(
    list(.brms_b_row("p1", "normal(0, 10)"))
  )
  code <- brms::make_stancode(
    fbps_formula(), d,
    prior = rbind(protected_row, out$prior)
  )
  # b[1] now carries exactly one contribution, and it is the protected row's
  # normal(0, 10) -- not a shrinkage density.
  expect_identical(fbps_prior_b_indices(code), c(1L, 2L, 3L, 4L))
  expect_true(fbps_declares(code, "normal_lpdf\\(b\\[1\\] \\| 0, 10\\)"))
  expect_identical(
    fbps_prior_b_indices(code, "double_exponential"), c(2L, 3L, 4L)
  )
})


# --- class-wide families ------------------------------------------------------ #

test_that("the regularized horseshoe compiles to the global-local parameters", {
  skip_if_not_installed("brms")
  d <- fbps_data()
  map <- .fb_sel_term_map(~ x1 + x2 + x3, d)
  blocks <- .fb_pen_resolve(
    fb_penalty(
      term(everything()) ~ regularized_horseshoe(
        p0 = 1, slab_scale = 2, slab_df = 4, df = 1, df_global = 1
      )
    ),
    map
  )
  out <- .fb_pen_compile_brms(blocks, c("x1", "x2", "x3"), n_train = nrow(d))
  expect_identical(out$route, "class_wide")
  expect_null(out$stanvars)

  fml <- y ~ x1 + x2 + x3 + (1 | g)
  code <- brms::make_stancode(fml, d, prior = out$prior)
  sd <- brms::make_standata(fml, d, prior = out$prior)

  # Presence and dimension of every latent parameter the prior introduces.
  expect_true(fbps_declares(code, "real<lower=0> hs_global;"))
  expect_true(fbps_declares(code, "real<lower=0> hs_slab;"))
  expect_true(fbps_declares(code, "vector<lower=0>\\[Kscales\\] hs_local;"))
  expect_true(fbps_declares(code, "vector\\[Kc\\] zb;"))
  expect_true(fbps_declares(code, "b = zb \\.\\* sdb;"))
  expect_equal(sd$Kscales, 3)
  expect_equal(sd$Kc, 3)

  # There is no per-coefficient prior row at all: the coupling is what makes
  # it a shrinkage prior, and expanding it row-wise would destroy it.
  expect_length(fbps_prior_b_indices(code), 0L)

  # p0 reaches the program through par_ratio, which brms turns into the
  # global scale. The relation is brms's, measured here rather than assumed.
  expect_equal(blocks[[1L]]$hyperparameters$par_ratio, 1 / 2)
  expect_equal(sd$hs_scale_global, 0.5 / sqrt(nrow(d)))
  expect_equal(sd$hs_df, 1)
  expect_equal(sd$hs_df_global, 1)
  expect_equal(sd$hs_scale_slab, 2)
  expect_equal(sd$hs_df_slab, 4)
})


test_that("the horseshoe hyperparameters travel from the grammar to the data block", {
  skip_if_not_installed("brms")
  d <- fbps_data()
  map <- .fb_sel_term_map(~ x1 + x2 + x3, d)
  blocks <- .fb_pen_resolve(
    fb_penalty(
      term(everything()) ~ horseshoe(
        p0 = 2, slab_scale = 3, slab_df = 6, df = 2, df_global = 5
      )
    ),
    map
  )
  out <- .fb_pen_compile_brms(blocks, c("x1", "x2", "x3"), n_train = nrow(d))
  sd <- brms::make_standata(y ~ x1 + x2 + x3 + (1 | g), d, prior = out$prior)
  expect_equal(sd$hs_df, 2)
  expect_equal(sd$hs_df_global, 5)
  expect_equal(sd$hs_scale_slab, 3)
  expect_equal(sd$hs_df_slab, 6)
  expect_equal(sd$hs_scale_global, (2 / 1) / sqrt(nrow(d)))
})


test_that("R2D2 compiles to its simplex and beta parameters", {
  skip_if_not_installed("brms")
  d <- fbps_data()
  map <- .fb_sel_term_map(~ x1 + x2 + x3, d)
  blocks <- .fb_pen_resolve(
    fb_penalty(
      term(everything()) ~ r2d2(mean_R2 = 0.3, prec_R2 = 3, cons_D2 = 0.4)
    ),
    map
  )
  out <- .fb_pen_compile_brms(blocks, c("x1", "x2", "x3"), n_train = nrow(d))
  fml <- y ~ x1 + x2 + x3 + (1 | g)
  code <- brms::make_stancode(fml, d, prior = out$prior)
  sd <- brms::make_standata(fml, d, prior = out$prior)

  expect_true(fbps_declares(code, "simplex\\[Kscales\\] R2D2_phi;"))
  expect_true(fbps_declares(code, "real<lower=0,upper=1> R2D2_R2;"))
  expect_true(fbps_declares(code, "dirichlet_lpdf\\(R2D2_phi \\| R2D2_cons_D2\\)"))
  expect_equal(sd$Kscales, 3)
  expect_length(fbps_prior_b_indices(code), 0L)
  expect_equal(sd$R2D2_mean_R2, 0.3)
  expect_equal(sd$R2D2_prec_R2, 3)
  expect_equal(as.numeric(sd$R2D2_cons_D2), rep(0.4, 3))
})


test_that("the class-wide refusal fires before brms is reached", {
  skip_if_not_installed("brms")
  d <- fbps_data()
  map <- fbps_map(d)
  pen <- fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 1))

  # Resolution happens upstream of every brms call, and this is the proof
  # rather than the claim: with brms's own entry points replaced by a fault,
  # the refusal that comes back is still the targeting one.
  testthat::local_mocked_bindings(
    make_stancode = function(...) stop("brms was reached"),
    make_standata = function(...) stop("brms was reached"),
    .package = "brms"
  )
  expect_error(
    .fb_pen_resolve(pen, map, protected = "p1"),
    class = "flexybayes_refusal_penalty_block_targeting_unsupported"
  )
  # And the same gate holds one altitude down, at the compiler, where the b
  # vector rather than the term map is what knows about the protected effect.
  blocks <- .fb_pen_resolve(pen, .fb_sel_term_map(~ x1 + x2 + x3, d))
  expect_error(
    .fb_pen_compile_brms(blocks, fbps_coef_all, n_train = nrow(d)),
    class = "flexybayes_refusal_penalty_block_targeting_unsupported"
  )
})


# --- the stanvar route --------------------------------------------------------- #

test_that("the elastic net indexes brms's b vector, and the index is right", {
  skip_if_not_installed("brms")
  d <- fbps_data()
  blocks <- .fb_pen_resolve(
    fb_penalty(term(x1, x2, x3) ~ elastic_net(lambda = 0.5, alpha = 0.25)),
    fbps_map(d),
    protected = "p1"
  )
  out <- .fb_pen_compile_brms(blocks, fbps_coef_all, n_train = nrow(d))
  expect_identical(out$route, "stanvar_density")
  expect_s3_class(out$stanvars, "stanvars")

  fml <- fbps_formula()
  code <- brms::make_stancode(fml, d, prior = out$prior, stanvars = out$stanvars)
  sd <- brms::make_standata(fml, d, prior = out$prior, stanvars = out$stanvars)

  # Declared at the block's own length, in the data block.
  expect_true(fbps_declares(code, "array\\[3\\] int cand_idx_1;"))
  expect_identical(as.integer(sd$cand_idx_1), c(2L, 3L, 4L))

  # The index mapping oracle: b[i] is the coefficient of X[, i + 1] under
  # brms's centred parameterisation, so cand_idx must select exactly the
  # candidate columns of the data and never the intercept or the protected
  # column.
  expect_identical(
    dimnames(sd$X)[[2L]][sd$cand_idx_1 + 1L], c("x1", "x2", "x3")
  )
  for (k in seq_along(sd$cand_idx_1)) {
    expect_equal(
      unname(sd$X[, sd$cand_idx_1[k] + 1L]),
      unname(d[[c("x1", "x2", "x3")[k]]])
    )
  }
  expect_false(1L %in% sd$cand_idx_1)

  # Exactly one shrinkage contribution, carrying the training row count and
  # both hyperparameters.
  target_lines <- grep(
    "cand_idx_1", grep("target \\+=", fbps_lines(code), value = TRUE),
    value = TRUE
  )
  expect_length(target_lines, 1L)
  expect_match(target_lines, "-60 \\* 0.5 \\*")
  expect_match(target_lines, "0.25 \\* sum\\(abs\\(b\\[cand_idx_1\\]\\)\\)")
  expect_match(target_lines, "\\(1 - 0.25\\) / 2 \\* dot_self\\(b\\[cand_idx_1\\]\\)")
})


test_that("the elastic net leaves the protected coefficient its own row", {
  skip_if_not_installed("brms")
  d <- fbps_data()
  blocks <- .fb_pen_resolve(
    fb_penalty(term(x1, x2, x3) ~ elastic_net(lambda = 0.5, alpha = 0.25)),
    fbps_map(d),
    protected = "p1"
  )
  out <- .fb_pen_compile_brms(blocks, fbps_coef_all, n_train = nrow(d))
  protected_row <- .brms_specs_to_object(
    list(.brms_b_row("p1", "normal(0, 10)"))
  )
  code <- brms::make_stancode(
    fbps_formula(), d,
    prior = rbind(protected_row, out$prior), stanvars = out$stanvars
  )
  expect_true(fbps_declares(code, "normal_lpdf\\(b\\[1\\] \\| 0, 10\\)"))
  # The only indexed contribution is the protected row's; the shrinkage term
  # reaches b through cand_idx_1 and cannot touch b[1].
  expect_identical(fbps_prior_b_indices(code), 1L)
})


test_that("n_train enters the elastic-net density and nothing else", {
  skip_if_not_installed("brms")
  d <- fbps_data()
  blocks <- .fb_pen_resolve(
    fb_penalty(term(x1, x2, x3) ~ elastic_net(lambda = 0.5, alpha = 0.25)),
    fbps_map(d),
    protected = "p1"
  )
  small <- .fb_pen_compile_brms(blocks, fbps_coef_all, n_train = 40L)
  large <- .fb_pen_compile_brms(blocks, fbps_coef_all, n_train = 400L)
  expect_match(small$prior$prior, "-40 \\*")
  expect_match(large$prior$prior, "-400 \\*")
  expect_identical(small$stanvars$cand_idx_1$sdata, large$stanvars$cand_idx_1$sdata)
})


# --- every family reaches a valid Stan program --------------------------------- #

test_that("each family in the table compiles to parseable Stan", {
  skip_if_not_installed("brms")
  d <- fbps_data()
  candidate_map <- .fb_sel_term_map(~ x1 + x2 + x3, d)
  fml <- y ~ x1 + x2 + x3 + (1 | g)
  penalties <- list(
    ridge = fb_penalty(term(everything()) ~ ridge(scale = 0.5)),
    laplace = fb_penalty(term(everything()) ~ laplace(scale = 0.5)),
    elastic_net = fb_penalty(
      term(everything()) ~ elastic_net(lambda = 0.4, alpha = 0.6)
    ),
    regularized_horseshoe = fb_penalty(
      term(everything()) ~ regularized_horseshoe(p0 = 1)
    ),
    r2d2 = fb_penalty(term(everything()) ~ r2d2())
  )
  for (nm in names(penalties)) {
    blocks <- .fb_pen_resolve(penalties[[nm]], candidate_map)
    out <- .fb_pen_compile_brms(
      blocks, c("x1", "x2", "x3"),
      n_train = nrow(d)
    )
    code <- brms::make_stancode(
      fml, d,
      prior = out$prior, stanvars = out$stanvars
    )
    expect_true(nzchar(code), info = nm)
    expect_true(fbps_declares(code, "^parameters \\{"), info = nm)
    expect_true(fbps_declares(code, "vector\\[Kc\\] (b|zb);"), info = nm)
    expect_identical(
      out$route, .fb_pen_family_table()[[nm]]$route,
      info = nm
    )
  }
})
