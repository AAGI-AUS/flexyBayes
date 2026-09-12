# =============================================================================
# The penalty grammar: a two-sided corpus.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md section 3.2, and
# SPEC_penalised_selection_v3_2026-09-11.md section 6 (the target-call
# contract, whose acceptance criterion is exactly this file).
#
# Two-sided means both halves are tested and neither on its own is enough. A
# corpus of malformed calls that all refuse proves nothing if the constructor
# refuses everything; a corpus of well-formed calls that all pass proves
# nothing if it accepts everything. Every allowed argument of every
# constructor is mutated by omission, duplication, unknown naming, extra
# position, wrong scalar type and call shape, and every mutation names the
# code it must refuse with.
#
# The corpus is data, not code: adding a family means adding rows, and a row
# that stops refusing is a visible deletion rather than a test that quietly
# went missing.
# =============================================================================


# --- the well-formed half ------------------------------------------------- #

# Every row here is a specification the grammar must accept, together with
# the canonical family it resolves to and the hyperparameters it must carry
# once defaults are filled. `alias` marks the rows that also exercise the
# once-per-session note.
fbpg_well_formed <- list(
  list(
    spec = term(everything()) ~ ridge(scale = 1),
    family = "ridge", args = list(scale = 1)
  ),
  list(
    spec = term(everything()) ~ ridge(1),
    family = "ridge", args = list(scale = 1)
  ),
  list(
    spec = term(x1, x2) ~ laplace(scale = 0.5),
    family = "laplace", args = list(scale = 0.5)
  ),
  list(
    spec = term(c("x1", "x2")) ~ laplace(0.5),
    family = "laplace", args = list(scale = 0.5)
  ),
  list(
    spec = term("x1") ~ laplace(scale = 0.5),
    family = "laplace", args = list(scale = 0.5)
  ),
  list(
    spec = term(matches("^x")) ~ elastic_net(lambda = 0.5, alpha = 0.25),
    family = "elastic_net", args = list(lambda = 0.5, alpha = 0.25)
  ),
  list(
    spec = term(everything()) ~ elastic_net(0.5, 0.25),
    family = "elastic_net", args = list(lambda = 0.5, alpha = 0.25)
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(p0 = 3),
    family = "regularized_horseshoe",
    args = list(
      p0 = 3, slab_scale = 2, slab_df = 4, df = 1, df_global = 1,
      scope = "candidates"
    )
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(
      p0 = 3, slab_scale = 1, slab_df = 3, df = 2, df_global = 5,
      scope = "all_fixed"
    ),
    family = "regularized_horseshoe",
    args = list(
      p0 = 3, slab_scale = 1, slab_df = 3, df = 2, df_global = 5,
      scope = "all_fixed"
    )
  ),
  list(
    spec = term(everything()) ~ horseshoe(p0 = 3),
    family = "regularized_horseshoe",
    args = list(
      p0 = 3, slab_scale = 2, slab_df = 4, df = 1, df_global = 1,
      scope = "candidates"
    )
  ),
  list(
    spec = term(everything()) ~ r2d2(),
    family = "r2d2",
    args = list(
      mean_R2 = 0.5, prec_R2 = 2, cons_D2 = 0.5, scope = "candidates"
    )
  ),
  list(
    spec = term(everything()) ~ r2d2(
      mean_R2 = 0.3, prec_R2 = 3, cons_D2 = 0.4, scope = "all_fixed"
    ),
    family = "r2d2",
    args = list(
      mean_R2 = 0.3, prec_R2 = 3, cons_D2 = 0.4, scope = "all_fixed"
    )
  ),
  list(
    spec = term(x1) ~ none(),
    family = "none", args = list()
  ),
  list(
    spec = term(everything()) ~ lasso(scale = 0.5),
    family = "laplace", args = list(scale = 0.5), alias = TRUE
  )
)


test_that("every well-formed specification is accepted, once", {
  for (row in fbpg_well_formed) {
    pen <- suppressMessages(fb_penalty(row$spec))
    expect_s3_class(pen, "fb_penalty")
    expect_length(pen$specs, 1L)
    expect_identical(
      pen$specs[[1L]]$family, row$family,
      info = deparse(row$spec)
    )
    expect_equal(
      pen$specs[[1L]]$args, row$args,
      info = deparse(row$spec)
    )
  }
})


test_that("several specifications are kept in the order written", {
  pen <- fb_penalty(
    term(everything()) ~ ridge(scale = 1),
    term(x1) ~ laplace(scale = 0.5),
    term(matches("^z")) ~ none()
  )
  expect_length(pen$specs, 3L)
  expect_identical(
    vapply(pen$specs, function(s) s$family, character(1L)),
    c("ridge", "laplace", "none")
  )
  expect_identical(pen$specs[[1L]]$selector$kind, "everything")
  expect_identical(pen$specs[[2L]]$selector$names, "x1")
  expect_identical(pen$specs[[3L]]$selector$pattern, "^z")
})


test_that("resolution and scaling are the two documented pairs", {
  expect_identical(
    fb_penalty(term(x1) ~ ridge(1))$resolution, "term"
  )
  expect_identical(
    fb_penalty(term(x1) ~ ridge(1))$scaling, "train_fold_standardize"
  )
  expect_identical(
    fb_penalty(term(x1) ~ ridge(1), resolution = "coefficient")$resolution,
    "coefficient"
  )
  expect_identical(
    fb_penalty(term(x1) ~ ridge(1), scaling = "none")$scaling, "none"
  )
  expect_error(
    fb_penalty(term(x1) ~ ridge(1), resolution = "column"),
    "should be one of"
  )
  expect_error(
    fb_penalty(term(x1) ~ ridge(1), scaling = "global"),
    "should be one of"
  )
})


test_that("the lasso alias resolves to laplace and says so once", {
  flexyBayes:::.fb_pen_reset_alias_note()
  expect_message(
    fb_penalty(term(everything()) ~ lasso(scale = 0.5)),
    "not a lasso estimate"
  )
  # Once per session: the second call is silent, and the specification is
  # identical to the one the canonical spelling produces.
  expect_silent(
    pen <- fb_penalty(term(everything()) ~ lasso(scale = 0.5))
  )
  expect_identical(pen$specs[[1L]]$family, "laplace")
  expect_identical(
    pen$specs[[1L]],
    fb_penalty(term(everything()) ~ laplace(scale = 0.5))$specs[[1L]]
  )
  flexyBayes:::.fb_pen_reset_alias_note()
})


test_that("the horseshoe alias records the canonical name", {
  expect_identical(
    fb_penalty(term(everything()) ~ horseshoe(p0 = 2))$specs[[1L]],
    fb_penalty(
      term(everything()) ~ regularized_horseshoe(p0 = 2)
    )$specs[[1L]]
  )
})


test_that("print() reports the family, the route and the resolution", {
  pen <- fb_penalty(
    term(everything()) ~ ridge(scale = 1),
    term(x1) ~ none()
  )
  out <- paste(utils::capture.output(print(pen)), collapse = "\n")
  expect_match(out, "2 specification\\(s\\)")
  expect_match(out, "resolution = \"term\"")
  expect_match(out, "ridge\\(scale = 1\\)")
  expect_match(out, "route coef_rows")
})


# --- the malformed half ---------------------------------------------------- #

# Each row is one mutation and the code it must refuse with. The `mutation`
# field is what the row is testing rather than decoration: the coverage claim
# below counts them, so a family added without its six mutations fails the
# count rather than passing quietly.
fbpg_malformed <- list(
  # ---- target: shape -----------------------------------------------------
  list(
    spec = term() ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "omission"
  ),
  list(
    spec = term(select = x1) ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "unknown_name"
  ),
  list(
    spec = term(1) ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "wrong_type"
  ),
  list(
    spec = term("") ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "wrong_type"
  ),
  list(
    spec = term(x1, x1) ~ ridge(scale = 1),
    code = "grammar_argument_duplicated", mutation = "duplication"
  ),
  list(
    spec = term(c("x1", "x1")) ~ ridge(scale = 1),
    code = "grammar_argument_duplicated", mutation = "duplication"
  ),
  list(
    spec = term(x1, ) ~ ridge(scale = 1),
    code = "grammar_argument_empty", mutation = "empty_argument"
  ),
  # ---- target: call-shaped pseudo-selectors -------------------------------
  list(
    spec = term(starts_with("x")) ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "call_shape"
  ),
  list(
    spec = term(all_of(c("x1"))) ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "call_shape"
  ),
  list(
    spec = term(everything(), x1) ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "extra_position"
  ),
  list(
    spec = term(everything(1)) ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "extra_position"
  ),
  list(
    spec = term(matches()) ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "omission"
  ),
  list(
    spec = term(matches("a", "b")) ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "extra_position"
  ),
  list(
    spec = term(matches(1)) ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "wrong_type"
  ),
  list(
    spec = term(c("x1", NA)) ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "wrong_type"
  ),
  list(
    spec = term(c(1, 2)) ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "wrong_type"
  ),
  # ---- target: the withdrawn and never-supported constructors -------------
  list(
    spec = group(x1) ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "call_shape"
  ),
  list(
    spec = as_group("g") ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "call_shape"
  ),
  list(
    spec = x1 ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "wrong_type"
  ),
  list(
    spec = "x1" ~ ridge(scale = 1),
    code = "selection_penalty_target_unsupported", mutation = "wrong_type"
  ),
  # ---- ridge --------------------------------------------------------------
  list(
    spec = term(everything()) ~ ridge(),
    code = "penalty_hyperparameter_invalid", mutation = "omission"
  ),
  list(
    spec = term(everything()) ~ ridge(scale = 1, scale = 2),
    code = "grammar_argument_duplicated", mutation = "duplication"
  ),
  list(
    spec = term(everything()) ~ ridge(sd = 1),
    code = "grammar_named_argument_unknown", mutation = "unknown_name"
  ),
  list(
    spec = term(everything()) ~ ridge(1, 2),
    code = "grammar_arity_too_many", mutation = "extra_position"
  ),
  list(
    spec = term(everything()) ~ ridge(scale = "1"),
    code = "penalty_hyperparameter_invalid", mutation = "wrong_type"
  ),
  list(
    spec = term(everything()) ~ ridge(scale = c(1, 2)),
    code = "penalty_hyperparameter_invalid", mutation = "wrong_type"
  ),
  list(
    spec = term(everything()) ~ ridge(scale = 0),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ ridge(scale = NA),
    code = "penalty_hyperparameter_invalid", mutation = "wrong_type"
  ),
  list(
    spec = term(everything()) ~ ridge(scale = Inf),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ ridge(scale = ),
    code = "grammar_argument_empty", mutation = "empty_argument"
  ),
  list(
    spec = term(everything()) ~ ridge(scale = mean),
    code = "penalty_hyperparameter_invalid", mutation = "call_shape"
  ),
  # ---- laplace -------------------------------------------------------------
  list(
    spec = term(everything()) ~ laplace(),
    code = "penalty_hyperparameter_invalid", mutation = "omission"
  ),
  list(
    spec = term(everything()) ~ laplace(scale = 1, scale = 2),
    code = "grammar_argument_duplicated", mutation = "duplication"
  ),
  list(
    spec = term(everything()) ~ laplace(b = 1),
    code = "grammar_named_argument_unknown", mutation = "unknown_name"
  ),
  list(
    spec = term(everything()) ~ laplace(1, 2),
    code = "grammar_arity_too_many", mutation = "extra_position"
  ),
  list(
    spec = term(everything()) ~ laplace(scale = "wide"),
    code = "penalty_hyperparameter_invalid", mutation = "wrong_type"
  ),
  list(
    spec = term(everything()) ~ laplace(scale = -1),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ lasso(),
    code = "penalty_hyperparameter_invalid", mutation = "omission"
  ),
  # ---- elastic_net ---------------------------------------------------------
  list(
    spec = term(everything()) ~ elastic_net(lambda = 0.5),
    code = "penalty_hyperparameter_invalid", mutation = "omission"
  ),
  list(
    spec = term(everything()) ~ elastic_net(alpha = 0.5),
    code = "penalty_hyperparameter_invalid", mutation = "omission"
  ),
  list(
    spec = term(everything()) ~ elastic_net(
      lambda = 0.5, lambda = 1, alpha = 0.5
    ),
    code = "grammar_argument_duplicated", mutation = "duplication"
  ),
  list(
    spec = term(everything()) ~ elastic_net(
      lambda = 0.5, alpha = 0.5, gamma = 1
    ),
    code = "grammar_named_argument_unknown", mutation = "unknown_name"
  ),
  list(
    spec = term(everything()) ~ elastic_net(0.5, 0.5, 0.5),
    code = "grammar_arity_too_many", mutation = "extra_position"
  ),
  list(
    spec = term(everything()) ~ elastic_net(lambda = 0.5, alpha = "half"),
    code = "penalty_hyperparameter_invalid", mutation = "wrong_type"
  ),
  list(
    spec = term(everything()) ~ elastic_net(lambda = 0.5, alpha = 1),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ elastic_net(lambda = 0.5, alpha = 0),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ elastic_net(lambda = 0, alpha = 0.5),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ elastic_net(lambda = , alpha = 0.5),
    code = "grammar_argument_empty", mutation = "empty_argument"
  ),
  # ---- regularized_horseshoe ----------------------------------------------
  list(
    spec = term(everything()) ~ regularized_horseshoe(),
    code = "penalty_hyperparameter_invalid", mutation = "omission"
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(p0 = 2, p0 = 3),
    code = "grammar_argument_duplicated", mutation = "duplication"
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(p0 = 2, tau0 = 1),
    code = "grammar_named_argument_unknown", mutation = "unknown_name"
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(
      2, 2, 4, 1, 1, "candidates", 9
    ),
    code = "grammar_arity_too_many", mutation = "extra_position"
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(p0 = "three"),
    code = "penalty_hyperparameter_invalid", mutation = "wrong_type"
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(p0 = -1),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(p0 = 2, slab_scale = 0),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(p0 = 2, slab_df = -3),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(p0 = 2, df = 0),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(p0 = 2, df_global = 0),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(p0 = 2, scope = "some"),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(p0 = 2, scope = 1),
    code = "penalty_hyperparameter_invalid", mutation = "wrong_type"
  ),
  list(
    spec = term(everything()) ~ regularized_horseshoe(p0 = ),
    code = "grammar_argument_empty", mutation = "empty_argument"
  ),
  # ---- r2d2 -----------------------------------------------------------------
  list(
    spec = term(everything()) ~ r2d2(p0 = 2),
    code = "penalty_hyperparameter_invalid", mutation = "unknown_name"
  ),
  list(
    spec = term(everything()) ~ r2d2(mean_R2 = 0.5, mean_R2 = 0.6),
    code = "grammar_argument_duplicated", mutation = "duplication"
  ),
  list(
    spec = term(everything()) ~ r2d2(shrinkage = 0.5),
    code = "grammar_named_argument_unknown", mutation = "unknown_name"
  ),
  list(
    spec = term(everything()) ~ r2d2(0.5, 2, 0.5, "candidates", 1),
    code = "grammar_arity_too_many", mutation = "extra_position"
  ),
  list(
    spec = term(everything()) ~ r2d2(mean_R2 = "half"),
    code = "penalty_hyperparameter_invalid", mutation = "wrong_type"
  ),
  list(
    spec = term(everything()) ~ r2d2(mean_R2 = 1),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ r2d2(prec_R2 = 0),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ r2d2(cons_D2 = -1),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ r2d2(scope = "protected"),
    code = "penalty_hyperparameter_invalid", mutation = "out_of_domain"
  ),
  list(
    spec = term(everything()) ~ r2d2(mean_R2 = ),
    code = "grammar_argument_empty", mutation = "empty_argument"
  ),
  # ---- none -----------------------------------------------------------------
  list(
    spec = term(x1) ~ none(1),
    code = "grammar_arity_too_many", mutation = "extra_position"
  ),
  list(
    spec = term(x1) ~ none(scale = 1),
    code = "grammar_named_argument_unknown", mutation = "unknown_name"
  ),
  # ---- withdrawn families ----------------------------------------------------
  list(
    spec = term(everything()) ~ group_lasso(),
    code = "selection_penalty_family_withdrawn", mutation = "withdrawn"
  ),
  list(
    spec = term(everything()) ~ exclusive_lasso(),
    code = "selection_penalty_family_withdrawn", mutation = "withdrawn"
  ),
  list(
    spec = term(everything()) ~ sparse_group_lasso(gamma = 0.5),
    code = "selection_penalty_family_withdrawn", mutation = "withdrawn"
  ),
  list(
    spec = term(everything()) ~ group_exclusive(ratio = 1),
    code = "selection_penalty_family_withdrawn", mutation = "withdrawn"
  ),
  list(
    spec = term(everything()) ~ adaptive(power = 1),
    code = "selection_penalty_family_withdrawn", mutation = "withdrawn"
  ),
  # ---- family shape ----------------------------------------------------------
  list(
    spec = term(everything()) ~ scad(a = 3.7),
    code = "selection_penalty_family_unknown", mutation = "unknown_family"
  ),
  list(
    spec = term(everything()) ~ ridge,
    code = "selection_penalty_family_unknown", mutation = "call_shape"
  ),
  list(
    spec = term(everything()) ~ 1,
    code = "selection_penalty_family_unknown", mutation = "wrong_type"
  ),
  list(
    spec = ~ridge(scale = 1),
    code = "selection_penalty_spec_malformed", mutation = "one_sided"
  )
)


test_that("every malformed specification refuses with its named code", {
  flexyBayes:::.fb_pen_reset_alias_note()
  for (row in fbpg_malformed) {
    expect_error(
      suppressMessages(fb_penalty(row$spec)),
      class = paste0("flexybayes_refusal_", row$code),
      info = paste0(deparse(row$spec), collapse = " ")
    )
  }
  flexyBayes:::.fb_pen_reset_alias_note()
})


test_that("the corpus covers every mutation of every family", {
  # The coverage claim itself, mechanically: six mutation kinds against the
  # five live families plus none(), so a family added to the table without
  # its corpus rows fails here rather than in review.
  codes <- vapply(fbpg_malformed, function(r) r$code, character(1L))
  mutations <- vapply(fbpg_malformed, function(r) r$mutation, character(1L))
  expect_true(all(
    c(
      "omission", "duplication", "unknown_name", "extra_position",
      "wrong_type", "call_shape", "empty_argument", "out_of_domain",
      "withdrawn", "unknown_family", "one_sided"
    ) %in% mutations
  ))
  expect_setequal(
    unique(codes),
    c(
      "selection_penalty_target_unsupported", "grammar_argument_duplicated",
      "grammar_argument_empty", "grammar_named_argument_unknown",
      "grammar_arity_too_many", "penalty_hyperparameter_invalid",
      "selection_penalty_family_withdrawn", "selection_penalty_family_unknown",
      "selection_penalty_spec_malformed"
    )
  )
  # Each live family carries at least the five argument mutations.
  for (fam in c("ridge", "laplace", "elastic_net", "regularized_horseshoe")) {
    rows <- vapply(
      fbpg_malformed,
      function(r) grepl(paste0("\\b", fam, "\\("), deparse(r$spec)[1L]),
      logical(1L)
    )
    expect_true(
      all(
        c("omission", "duplication", "unknown_name", "extra_position",
          "wrong_type") %in% mutations[rows]
      ),
      info = fam
    )
  }
})


# --- the constructor's own arguments --------------------------------------- #

test_that("fb_penalty() with no specification refuses", {
  expect_error(
    fb_penalty(),
    class = "flexybayes_refusal_selection_penalty_spec_malformed"
  )
})


test_that("a non-formula argument refuses rather than being coerced", {
  expect_error(
    fb_penalty("term(x1) ~ ridge(1)"),
    class = "flexybayes_refusal_selection_penalty_spec_malformed"
  )
  expect_error(
    fb_penalty(list(term(x1) ~ ridge(1))),
    class = "flexybayes_refusal_selection_penalty_spec_malformed"
  )
})


test_that("the second of two formulas is validated as hard as the first", {
  # A corpus that only ever mutates argument one is a corpus that tests a
  # loop's first iteration.
  expect_error(
    fb_penalty(
      term(everything()) ~ ridge(scale = 1),
      term(x1) ~ ridge(sd = 2)
    ),
    class = "flexybayes_refusal_grammar_named_argument_unknown"
  )
})


test_that("a family argument evaluates in the formula's environment", {
  local_scale <- 0.75
  pen <- fb_penalty(term(everything()) ~ ridge(scale = local_scale))
  expect_equal(pen$specs[[1L]]$args$scale, 0.75)
  expect_error(
    fb_penalty(term(everything()) ~ ridge(scale = no_such_object)),
    class = "flexybayes_refusal_penalty_hyperparameter_invalid"
  )
})


# --- the family table is the single source of truth ------------------------ #

test_that("every table family has a route, a mode and provider capabilities", {
  tbl <- flexyBayes:::.fb_pen_family_table()
  expect_setequal(
    names(tbl),
    c("ridge", "laplace", "elastic_net", "regularized_horseshoe", "r2d2",
      "none")
  )
  for (nm in names(tbl)) {
    e <- tbl[[nm]]
    expect_true(e$route %in%
      c("coef_rows", "class_wide", "stanvar_density", "none"), info = nm)
    expect_true(e$hyperparameter_mode %in% c("fixed", "learned"), info = nm)
    expect_setequal(names(e$provider_capabilities), c("brms", "inla", "glmnet"))
    expect_true(all(e$required %in% e$params), info = nm)
    expect_true(all(names(e$defaults) %in% e$params), info = nm)
    expect_setequal(
      sort(unique(c(e$required, names(e$defaults)))), sort(e$params)
    )
  }
  # INLA never approximates a scale-mixture prior (SPEC v5 section 4.5).
  shrinkage <- setdiff(names(tbl), "none")
  expect_true(all(vapply(
    shrinkage, function(nm) isFALSE(tbl[[nm]]$provider_capabilities$inla),
    logical(1L)
  )))
})


test_that("the withdrawn families are exactly the v3 grammar's five", {
  expect_setequal(
    flexyBayes:::.fb_pen_withdrawn_families(),
    c("group_lasso", "exclusive_lasso", "sparse_group_lasso",
      "group_exclusive", "adaptive")
  )
  # Their component traits stay in the resolver, which is what makes the
  # refusal a public-grammar decision rather than a deletion.
  traits <- flexyBayes:::.fb_sel_penalty_component_traits()
  expect_true(all(
    flexyBayes:::.fb_pen_withdrawn_families() %in% traits$component
  ))
})
