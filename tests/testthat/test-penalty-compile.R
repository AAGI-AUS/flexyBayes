# =============================================================================
# The penalty IR: resolution against a model matrix, the role of every column,
# the fold-local scaling and its inverse.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 4.2, 4.3 and
# 4.4, Unit A.
#
# Two things here are external facts and are treated as such. The name brms
# gives a design column is measured against `brms::get_prior()` rather than
# recalled, because a prior keyed on a name brms does not use is a prior that
# is silently never applied -- the 0.9.2 defect, in a new place. And the
# back-transform is checked against an ordinary least-squares fit on the
# unstandardised design: the transform is exact linear algebra, so the oracle
# is exact too, and the tolerance is a contract floor rather than a fitted
# one.
#
# Every refusal in this file is shown firing on the input that should produce
# it, and the two invariants that cannot be reached through the public path
# (block disjointness, and the b-vector agreement inside the compiler) are
# fired directly on hand-built blocks. A guard nobody has watched fail is not
# a guard.
# =============================================================================


# --- fixtures --------------------------------------------------------------- #

# A design with all three column kinds: numerics, a three-level factor whose
# contrast columns are one term, and an interaction that generates two more.
fbpc_data <- function(n = 60L, seed = 4021L) {
  set.seed(seed)
  data.frame(
    y = stats::rnorm(n),
    p1 = stats::rnorm(n),
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n),
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    g = factor(rep(seq_len(6L), length.out = n)),
    stringsAsFactors = FALSE
  )
}

fbpc_map <- function(d = fbpc_data()) {
  .fb_sel_term_map(~ p1 + x1 + x2 * f, d)
}

# The candidate coefficient names of that map with `p1` protected, in
# model-matrix order. Written out rather than derived, so a change in the
# resolver shows up here as a disagreement rather than as two functions
# agreeing about the wrong thing.
fbpc_candidates <- c("x1", "x2", "fb", "fc", "x2:fb", "x2:fc")

# The whole b vector of the extended model, intercept excluded.
fbpc_coef_all <- c("p1", "x1", "x2", "fb", "fc", "x2:fb", "x2:fc")


# --- resolution against the model matrix ------------------------------------ #

test_that("everything() resolves to every candidate column, in order", {
  map <- fbpc_map()
  pen <- fb_penalty(term(everything()) ~ ridge(scale = 1))
  blocks <- .fb_pen_resolve(pen, map, protected = "p1")

  expect_s3_class(blocks, "fb_prior_blocks")
  expect_length(blocks, 1L)
  expect_s3_class(blocks[[1L]], "fb_prior_block")
  expect_identical(blocks[[1L]]$columns, seq_len(6L))
  expect_identical(blocks[[1L]]$coef_names, fbpc_candidates)
  expect_identical(blocks[[1L]]$route, "coef_rows")
  expect_identical(blocks[[1L]]$hyperparameter_mode, "fixed")
  expect_false(blocks[[1L]]$protected)
})


test_that("a factor is one term and all of its contrast columns", {
  map <- fbpc_map()
  pen <- fb_penalty(term(f) ~ ridge(scale = 1))
  blocks <- .fb_pen_resolve(pen, map, protected = "p1")
  expect_identical(blocks[[1L]]$coef_names, c("fb", "fc"))

  # And an interaction is one term and every column it generates.
  pen2 <- fb_penalty(term(c("x2:f")) ~ ridge(scale = 1))
  blocks2 <- .fb_pen_resolve(pen2, map, protected = "p1")
  expect_identical(blocks2[[1L]]$coef_names, c("x2:fb", "x2:fc"))
})


test_that("matches() runs on term labels at resolution = 'term'", {
  map <- fbpc_map()
  pen <- fb_penalty(term(matches("^x")) ~ ridge(scale = 1))
  blocks <- .fb_pen_resolve(pen, map, protected = "p1")
  # x1, x2 and x2:f -- the interaction's label starts with x as well.
  expect_identical(blocks[[1L]]$coef_names, c("x1", "x2", "x2:fb", "x2:fc"))
})


test_that("matches() runs on coefficient names at resolution = 'coefficient'", {
  map <- fbpc_map()
  pen <- fb_penalty(
    term(matches("^x1$")) ~ ridge(scale = 1),
    resolution = "coefficient"
  )
  blocks <- .fb_pen_resolve(pen, map, protected = "p1")
  expect_identical(blocks[[1L]]$coef_names, "x1")
})


test_that("later specifications override earlier ones for the coordinates they name", {
  map <- fbpc_map()
  pen <- fb_penalty(
    term(everything()) ~ ridge(scale = 1),
    term(matches("^x")) ~ laplace(scale = 0.5),
    term(x1) ~ none()
  )
  blocks <- .fb_pen_resolve(pen, map, protected = "p1")
  role <- unlist(lapply(blocks, function(b) {
    stats::setNames(rep(b$family, length(b$coef_names)), b$coef_names)
  }))
  expect_identical(role[["x1"]], "none")
  expect_identical(role[["x2"]], "laplace")
  expect_identical(role[["x2:fb"]], "laplace")
  expect_identical(role[["fb"]], "ridge")
  expect_identical(role[["fc"]], "ridge")
  expect_length(role, 6L)
})


test_that("two equally broad selectors on one coordinate refuse", {
  map <- fbpc_map()
  expect_error(
    .fb_pen_resolve(
      fb_penalty(
        term(x1, x2) ~ ridge(scale = 1),
        term(x1) ~ laplace(scale = 0.5)
      ),
      map,
      protected = "p1"
    ),
    class = "flexybayes_refusal_selection_penalty_target_ambiguous"
  )
  expect_error(
    .fb_pen_resolve(
      fb_penalty(
        term(matches("^x")) ~ ridge(scale = 1),
        term(matches("^x1")) ~ laplace(scale = 0.5)
      ),
      map,
      protected = "p1"
    ),
    class = "flexybayes_refusal_selection_penalty_target_ambiguous"
  )
  # Different breadth is the ordinary idiom and must not trip the gate.
  expect_silent(
    .fb_pen_resolve(
      fb_penalty(
        term(everything()) ~ ridge(scale = 1),
        term(x1) ~ laplace(scale = 0.5)
      ),
      map,
      protected = "p1"
    )
  )
})


test_that("a target the candidate set does not contain refuses", {
  map <- fbpc_map()
  expect_error(
    .fb_pen_resolve(
      fb_penalty(term(z9) ~ ridge(scale = 1)), map, protected = "p1"
    ),
    class = "flexybayes_refusal_selection_penalty_target_not_in_model"
  )
  # A protected term is not a candidate, and the refusal says so.
  expect_error(
    .fb_pen_resolve(
      fb_penalty(term(p1) ~ ridge(scale = 1)), map, protected = "p1"
    ),
    class = "flexybayes_refusal_selection_penalty_target_not_in_model"
  )
  # A regular expression that matches nothing is a target, not an empty set.
  expect_error(
    .fb_pen_resolve(
      fb_penalty(term(matches("^nothing")) ~ ridge(scale = 1)),
      map,
      protected = "p1"
    ),
    class = "flexybayes_refusal_selection_penalty_target_not_in_model"
  )
  # Protecting the whole model leaves nothing to shrink.
  expect_error(
    .fb_pen_resolve(
      fb_penalty(term(everything()) ~ ridge(scale = 1)),
      map,
      protected = c("p1", "x1", "x2", "f", "x2:f")
    ),
    class = "flexybayes_refusal_selection_penalty_target_not_in_model"
  )
})


test_that("one column in two blocks refuses, on hand-built blocks", {
  # Unreachable through the public path, because later-overrides-earlier is
  # what makes it unreachable. The assertion exists so that a change to the
  # resolution rule cannot quietly put two densities on one coefficient, and
  # the only way to watch it fire is to build the state it forbids.
  map <- fbpc_map()
  blocks <- .fb_pen_resolve(
    fb_penalty(
      term(x1) ~ ridge(scale = 1),
      term(x2) ~ laplace(scale = 0.5)
    ),
    map,
    protected = "p1"
  )
  blocks[[2L]]$columns <- blocks[[1L]]$columns
  blocks[[2L]]$coef_names <- blocks[[1L]]$coef_names
  expect_error(
    .fb_pen_assert_disjoint(blocks),
    class = "flexybayes_refusal_penalty_block_overlap"
  )
})


# --- hyperparameters that need the resolved design -------------------------- #

test_that("p0 outside [1, p - 1] refuses, with p counted in design columns", {
  # No protected term, so the class-wide family is admissible and the only
  # thing under test is the domain of p0. Five candidate terms generate seven
  # design columns, and p is the column count: p0 = 6 is admissible and would
  # be out of range if terms were being counted.
  map <- fbpc_map()
  expect_silent(
    .fb_pen_resolve(
      fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 6)), map
    )
  )
  expect_error(
    .fb_pen_resolve(
      fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 7)), map
    ),
    class = "flexybayes_refusal_penalty_hyperparameter_invalid"
  )
  expect_error(
    .fb_pen_resolve(
      fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 0.5)), map
    ),
    class = "flexybayes_refusal_penalty_hyperparameter_invalid"
  )
})


test_that("par_ratio is p0 / (p - p0) and is stored on the block", {
  map <- fbpc_map()
  blocks <- .fb_pen_resolve(
    fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)), map
  )
  expect_equal(blocks[[1L]]$hyperparameters$par_ratio, 2 / 5)
  expect_equal(blocks[[1L]]$hyperparameters$p, 7L)
})


# --- coordinate semantics ---------------------------------------------------- #

test_that("a coordinate-wise family on one contrast column refuses", {
  map <- fbpc_map()
  expect_error(
    .fb_pen_resolve(
      fb_penalty(
        term(fb) ~ laplace(scale = 0.5),
        resolution = "coefficient"
      ),
      map,
      protected = "p1"
    ),
    class = "flexybayes_refusal_penalty_factor_contrast_ambiguous"
  )
  expect_error(
    .fb_pen_resolve(
      fb_penalty(
        term(matches("^x2:")) ~ elastic_net(lambda = 0.5, alpha = 0.5),
        resolution = "coefficient"
      ),
      map,
      protected = "p1"
    ),
    class = "flexybayes_refusal_penalty_factor_contrast_ambiguous"
  )
  # Ridge is an L2 ball: rotate the block and the density is unchanged, so
  # the contrast matrix cannot move its answer and the gate must not fire.
  expect_silent(
    .fb_pen_resolve(
      fb_penalty(term(fb) ~ ridge(scale = 1), resolution = "coefficient"),
      map,
      protected = "p1"
    )
  )
  # Nor does it fire on a single-column term at coefficient resolution.
  expect_silent(
    .fb_pen_resolve(
      fb_penalty(term(x1) ~ laplace(scale = 0.5), resolution = "coefficient"),
      map,
      protected = "p1"
    )
  )
  # Nor at term resolution, where the whole block is the target.
  expect_silent(
    .fb_pen_resolve(
      fb_penalty(term(f) ~ laplace(scale = 0.5)), map,
      protected = "p1"
    )
  )
})


# --- class-wide targeting ----------------------------------------------------- #

test_that("a class-wide family with a protected fixed effect refuses", {
  map <- fbpc_map()
  expect_error(
    .fb_pen_resolve(
      fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)),
      map,
      protected = "p1"
    ),
    class = "flexybayes_refusal_penalty_block_targeting_unsupported"
  )
  expect_error(
    .fb_pen_resolve(
      fb_penalty(term(everything()) ~ r2d2()), map,
      protected = "p1"
    ),
    class = "flexybayes_refusal_penalty_block_targeting_unsupported"
  )
  # An intercept-only protected part is exactly the admissible case.
  expect_silent(
    .fb_pen_resolve(
      fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 2)), map
    )
  )
})


test_that("scope = 'all_fixed' is the documented opt-in", {
  map <- fbpc_map()
  blocks <- .fb_pen_resolve(
    fb_penalty(
      term(everything()) ~ regularized_horseshoe(p0 = 2, scope = "all_fixed")
    ),
    map,
    protected = "p1"
  )
  expect_identical(blocks[[1L]]$hyperparameters$scope, "all_fixed")
  expect_identical(blocks[[1L]]$coef_names, fbpc_candidates)

  # It cannot coexist with a second block: the global-local prior is on every
  # population-level coefficient, so there is nothing left to give one to.
  expect_error(
    .fb_pen_resolve(
      fb_penalty(
        term(everything()) ~ regularized_horseshoe(
          p0 = 2, scope = "all_fixed"
        ),
        term(x1) ~ ridge(scale = 1)
      ),
      map,
      protected = "p1"
    ),
    class = "flexybayes_refusal_penalty_block_targeting_unsupported"
  )
})


test_that("a class-wide family covering only some candidates refuses", {
  map <- fbpc_map()
  expect_error(
    .fb_pen_resolve(
      fb_penalty(term(x1) ~ regularized_horseshoe(p0 = 2)), map
    ),
    class = "flexybayes_refusal_penalty_block_targeting_unsupported"
  )
  # Two class-wide families cannot be aimed at two blocks either: brms
  # carries one special prior on class "b".
  expect_error(
    .fb_pen_resolve(
      fb_penalty(
        term(everything()) ~ regularized_horseshoe(p0 = 2),
        term(x1) ~ r2d2()
      ),
      map
    ),
    class = "flexybayes_refusal_penalty_block_targeting_unsupported"
  )
})


# --- brms coefficient names: the independent oracle --------------------------- #

test_that("the coefficient names are the ones brms will use", {
  skip_if_not_installed("brms")
  # An awkward design on purpose: a factor whose levels carry a space and a
  # hyphen, and an interaction with it. brms deletes the space and
  # transliterates the hyphen, and a prior keyed on the model-matrix spelling
  # would never be applied.
  set.seed(77L)
  n <- 40L
  d <- data.frame(
    y = stats::rnorm(n),
    x1 = stats::rnorm(n),
    f = factor(rep(c("a", "b c", "d-e"), length.out = n)),
    stringsAsFactors = FALSE
  )
  fml <- y ~ x1 * f
  map <- .fb_sel_term_map(~ x1 * f, d)
  ours <- .fb_pen_brms_coef_names(
    map$columns$column[map$columns$term_index > 0L]
  )
  theirs <- brms::get_prior(fml, d)
  theirs <- theirs$coef[theirs$class == "b" & nzchar(theirs$coef)]
  expect_setequal(ours, theirs)
  expect_true(all(c("fbc", "fdMe", "x1:fbc", "x1:fdMe") %in% ours))
})


# --- the compiler's b-vector arithmetic ---------------------------------------- #

test_that("param_map places every coefficient at its b index", {
  skip_if_not_installed("brms")
  map <- fbpc_map()
  blocks <- .fb_pen_resolve(
    fb_penalty(
      term(everything()) ~ ridge(scale = 1),
      term(x1) ~ laplace(scale = 0.5)
    ),
    map,
    protected = "p1"
  )
  out <- .fb_pen_compile_brms(blocks, fbpc_coef_all, n_train = 60L)
  expect_identical(
    out$param_map$coefficient, fbpc_coef_all[-1L]
  )
  expect_identical(out$param_map$b_index, 2:7)
  expect_identical(
    out$param_map$family[out$param_map$coefficient == "x1"], "laplace"
  )
  expect_setequal(out$route, "coef_rows")
  # The protected coefficient is in the b vector and in no block.
  expect_false("p1" %in% out$param_map$coefficient)
})


test_that("the compiler refuses blocks the b vector does not contain", {
  skip_if_not_installed("brms")
  map <- fbpc_map()
  blocks <- .fb_pen_resolve(
    fb_penalty(term(everything()) ~ ridge(scale = 1)), map,
    protected = "p1"
  )
  expect_error(
    .fb_pen_compile_brms(blocks, c("p1", "x1"), n_train = 60L),
    "not in the b vector"
  )
  expect_error(
    .fb_pen_compile_brms(blocks, c(fbpc_coef_all, "x1"), n_train = 60L),
    "repeats"
  )
  expect_error(
    .fb_pen_compile_brms(blocks, fbpc_coef_all, n_train = 0L),
    "positive number"
  )
})


test_that("an explicit none() block emits no prior row but keeps its map row", {
  skip_if_not_installed("brms")
  map <- fbpc_map()
  blocks <- .fb_pen_resolve(
    fb_penalty(
      term(everything()) ~ ridge(scale = 1),
      term(x1) ~ none()
    ),
    map,
    protected = "p1"
  )
  out <- .fb_pen_compile_brms(blocks, fbpc_coef_all, n_train = 60L)
  expect_equal(nrow(out$prior), 5L)
  expect_false("x1" %in% out$prior$coef)
  expect_identical(
    out$param_map$route[out$param_map$coefficient == "x1"], "none"
  )
  expect_setequal(out$route, "coef_rows")

  # A penalty that is nothing but none() rows compiles to no prior at all.
  blocks2 <- .fb_pen_resolve(
    fb_penalty(term(everything()) ~ none()), map,
    protected = "p1"
  )
  out2 <- .fb_pen_compile_brms(blocks2, fbpc_coef_all, n_train = 60L)
  expect_null(out2$prior)
  expect_null(out2$stanvars)
  expect_identical(out2$route, "none")
})


# --- fold-local scaling --------------------------------------------------------- #

fbpc_numeric_design <- function(n = 80L, p = 4L, seed = 909L) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p, mean = 3, sd = 2), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  X
}

fbpc_numeric_blocks <- function(X, scaling = "train_fold_standardize") {
  d <- as.data.frame(X)
  d$y <- stats::rnorm(nrow(X))
  map <- .fb_sel_term_map(
    stats::as.formula(paste("~", paste(colnames(X), collapse = " + "))), d
  )
  .fb_pen_resolve(
    fb_penalty(term(everything()) ~ ridge(scale = 1), scaling = scaling),
    map
  )
}


test_that("the centre and scale are learned on the training rows alone", {
  X <- fbpc_numeric_design()
  blocks <- fbpc_numeric_blocks(X)
  train <- seq_len(50L)
  blocks <- .fb_pen_scaling(blocks, X, train_rows = train)

  expect_equal(unname(blocks[[1L]]$centre), unname(colMeans(X[train, ])))
  expect_equal(
    unname(blocks[[1L]]$scale),
    unname(apply(X[train, ], 2L, stats::sd))
  )
  expect_identical(names(blocks[[1L]]$centre), colnames(X))

  Xs <- .fb_pen_apply_scaling(X, blocks)
  # Exactly standardised on the training rows.
  expect_equal(unname(colMeans(Xs[train, ])), rep(0, ncol(X)))
  expect_equal(unname(apply(Xs[train, ], 2L, stats::sd)), rep(1, ncol(X)))
  # And NOT standardised on the assessment rows, which is the point: an
  # assessment row transformed by its own moments has entered the training
  # design, and the out-of-fold score is then optimistic by an amount nobody
  # measures.
  expect_false(isTRUE(all.equal(
    unname(colMeans(Xs[-train, , drop = FALSE])), rep(0, ncol(X))
  )))
})


test_that("scaling = 'none' stores the identity transform", {
  X <- fbpc_numeric_design()
  blocks <- fbpc_numeric_blocks(X, scaling = "none")
  blocks <- .fb_pen_scaling(blocks, X, train_rows = seq_len(50L))
  expect_equal(unname(blocks[[1L]]$centre), rep(0, ncol(X)))
  expect_equal(unname(blocks[[1L]]$scale), rep(1, ncol(X)))
  expect_equal(.fb_pen_apply_scaling(X, blocks), X)
})


test_that("a column with no spread in the training fold refuses", {
  X <- fbpc_numeric_design()
  X[seq_len(50L), 2L] <- 1
  blocks <- fbpc_numeric_blocks(X)
  expect_error(
    .fb_pen_scaling(blocks, X, train_rows = seq_len(50L)),
    class = "flexybayes_refusal_penalty_scaling_column_constant"
  )
  # The same column is fine once the fold sees its variation.
  expect_silent(.fb_pen_scaling(blocks, X, train_rows = 40:80))
})


test_that("applying an unscaled block is an error, not a silent pass-through", {
  X <- fbpc_numeric_design()
  blocks <- fbpc_numeric_blocks(X)
  expect_error(
    .fb_pen_apply_scaling(X, blocks), "carries no scaling constants"
  )
  expect_error(
    .fb_pen_back_transform(matrix(0, 1L, 1L, dimnames = list(NULL, "b_x1")),
      blocks
    ),
    "carries no scaling constants"
  )
})


test_that("the argument guards on the scaling entry points fire", {
  X <- fbpc_numeric_design()
  blocks <- fbpc_numeric_blocks(X)
  expect_error(.fb_pen_scaling(blocks, as.data.frame(X), 1:50), "must be a matrix")
  expect_error(.fb_pen_scaling(blocks, X, 1L), "at least two rows")
  expect_error(.fb_pen_scaling(blocks, X, c(1L, 500L)), "outside the rows")
  expect_error(.fb_pen_apply_scaling(as.data.frame(X), blocks), "must be a matrix")
  expect_error(
    .fb_pen_back_transform(matrix(0, 1L, 1L), blocks), "column names"
  )
})


# --- the back transform, against an independent least-squares oracle ---------- #

test_that("back-transformed coefficients are the unstandardised fit", {
  # The oracle is deliberately not another flexyBayes function. Standardising
  # a design and unstandardising its coefficients is exact linear algebra:
  # the ordinary least-squares fit on the raw design is what the transform
  # has to reproduce, to machine precision and not to a tolerance anyone
  # chose.
  X <- fbpc_numeric_design()
  n <- nrow(X)
  beta <- c(1.5, -0.8, 0, 2.2)
  y <- as.numeric(4 + X %*% beta + stats::rnorm(n, sd = 0.3))

  blocks <- fbpc_numeric_blocks(X)
  blocks <- .fb_pen_scaling(blocks, X, train_rows = seq_len(n))
  Xs <- .fb_pen_apply_scaling(X, blocks)

  fit_std <- stats::lm.fit(cbind(Intercept = 1, Xs), y)$coefficients
  fit_raw <- stats::lm.fit(cbind(Intercept = 1, X), y)$coefficients

  draws <- matrix(
    fit_std,
    nrow = 1L,
    dimnames = list(NULL, paste0("b_", names(fit_std)))
  )
  back <- .fb_pen_back_transform(draws, blocks)

  expect_equal(
    unname(back[1L, paste0("b_", colnames(X))]),
    unname(fit_raw[colnames(X)]),
    tolerance = 1e-10
  )
  expect_equal(
    unname(back[1L, "b_Intercept"]),
    unname(fit_raw[["Intercept"]]),
    tolerance = 1e-10
  )
})


test_that("the back transform leaves untargeted columns and draws alone", {
  X <- fbpc_numeric_design()
  d <- as.data.frame(X)
  d$y <- stats::rnorm(nrow(X))
  map <- .fb_sel_term_map(~ x1 + x2 + x3 + x4, d)
  blocks <- .fb_pen_resolve(
    fb_penalty(
      term(x1, x2) ~ ridge(scale = 1),
      term(x3, x4) ~ none()
    ),
    map
  )
  blocks <- .fb_pen_scaling(blocks, X, train_rows = seq_len(nrow(X)))
  draws <- matrix(
    seq_len(6L) / 10,
    nrow = 1L,
    dimnames = list(NULL, c("b_x1", "b_x2", "b_x3", "b_x4", "sigma", "lp__"))
  )
  back <- .fb_pen_back_transform(draws, blocks)
  # A none() block has scale 1, so its coefficients come back unchanged, and
  # so does anything that is not a coefficient at all.
  expect_equal(back[1L, "b_x3"], draws[1L, "b_x3"])
  expect_equal(back[1L, "sigma"], draws[1L, "sigma"])
  expect_equal(back[1L, "lp__"], draws[1L, "lp__"])
  expect_false(isTRUE(all.equal(back[1L, "b_x1"], draws[1L, "b_x1"])))
})


test_that("a block object prints its family, route and resolved columns", {
  map <- fbpc_map()
  blocks <- .fb_pen_resolve(
    fb_penalty(term(everything()) ~ ridge(scale = 1)), map,
    protected = "p1"
  )
  out <- paste(utils::capture.output(print(blocks)), collapse = "\n")
  expect_match(out, "1 block\\(s\\)")
  expect_match(out, "ridge")
  expect_match(out, "route: coef_rows")
  expect_match(out, "x2:fc")
})
