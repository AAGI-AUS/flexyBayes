# Penalty compiler --- joint prior blocks to a brms prior object
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 4.2, 4.3 and
# 4.4, Unit A.
#
# One function turns the resolved IR into the two things brms takes: a prior
# data frame and, where a family needs one, a `stanvars` object. Protected
# coefficients keep their own `fb_prior()` rows and flow through
# R/priors_to_brms.R unchanged; the two prior data frames meet by `rbind()`
# at the fit.
#
# Everything this file asserts about brms was measured against the installed
# brms rather than recalled. The argument names of `brms::horseshoe()` and
# `brms::R2D2()` are read off `formals()` in test-penalty-stancode.R; the
# index convention of the `b` vector is read off `brms::make_standata()`,
# which is the only place the correspondence between a design column and a
# Stan coefficient is written down. That matters more here than anywhere else
# in the module: an off-by-one in `cand_idx` shrinks the wrong coefficient,
# and the generated program is perfectly valid Stan either way.
#
# The index convention, measured (brms 2.23.0, `make_standata()` on a
# Gaussian mixed model): with the default centred parameterisation the data
# block carries `X` with `K` columns of which the first is the intercept,
# `Xc` with `Kc = K - 1` columns, and `Xc[, i - 1] = X[, i] - means_X[i - 1]`.
# The parameter `b` has length `Kc`. So `b[j]` is population-level
# coefficient `j` of the model with the intercept excluded, which is exactly
# the order of `coef_names_all`.


# The two non-linear parts of a `scope = "block"` reference model,
# `y ~ a + b`: the protected part and the candidate block.
.FB_NL_PROTECTED_PART <- "a"
.FB_NL_BLOCK_PART <- "b"


# --- the compiler ---------------------------------------------------------- #

# .fb_pen_compile_brms() --- blocks to (prior, stanvars, param_map, route,
# nl_block). `nl_block` is NULL unless a class-wide block has
# `scope = "block"`, and then names the two non-linear parts and the block's
# coefficients (R/penalty_block_nl.R).
#
# `coef_names_all` is the brms `b` vector in order, intercept excluded. It is
# an argument rather than something derived here because only the caller
# knows the extended model's fixed part: the candidates are one half of it
# and the protected effects are the other, and a compiler that guessed would
# be guessing at exactly the fact `cand_idx` depends on.
#
# `n_train` is the number of training rows the reference fit will see. It
# enters only the elastic-net density, whose canonical objective is stated
# per observation.
.fb_pen_compile_brms <- function(blocks, coef_names_all, n_train) {
  .check_installed(
    "brms",
    "Package 'brms' is required to compile a flexyBayes penalty into a ",
    "brms prior object. Install via install.packages('brms')."
  )
  coef_names_all <- as.character(coef_names_all)
  if (anyDuplicated(coef_names_all) > 0L) {
    stop(
      ".fb_pen_compile_brms(): `coef_names_all` repeats ",
      paste(
        unique(coef_names_all[duplicated(coef_names_all)]),
        collapse = ", "
      ),
      ". The b vector names one coefficient per position.",
      call. = FALSE
    )
  }
  if (!is.numeric(n_train) || length(n_train) != 1L || n_train < 1) {
    stop(
      ".fb_pen_compile_brms(): `n_train` must be a single positive number.",
      call. = FALSE
    )
  }

  .fb_pen_assert_targeting(blocks, coef_names_all)

  missing_coefs <- setdiff(
    unlist(lapply(blocks, function(b) b$coef_names), use.names = FALSE),
    coef_names_all
  )
  if (length(missing_coefs) > 0L) {
    stop(
      ".fb_pen_compile_brms(): these block coefficients are not in the b ",
      "vector: ", paste(missing_coefs, collapse = ", "),
      ". The resolved blocks and the fitted design disagree.",
      call. = FALSE
    )
  }

  rows <- list()
  stanvars <- NULL
  map_rows <- list()
  nl_block <- NULL

  for (b in blocks) {
    idx <- match(b$coef_names, coef_names_all)
    map_rows[[length(map_rows) + 1L]] <- data.frame(
      coefficient = b$coef_names,
      b_index = as.integer(idx),
      block = rep(b$id, length(idx)),
      family = rep(b$family, length(idx)),
      route = rep(b$route, length(idx)),
      hyperparameter_mode = rep(b$hyperparameter_mode, length(idx)),
      stringsAsFactors = FALSE
    )

    if (identical(b$route, "none")) {
      next
    }

    if (identical(b$route, "coef_rows")) {
      density <- .fb_pen_coef_density(b)
      for (nm in b$coef_names) {
        rows[[length(rows) + 1L]] <- brms::set_prior(
          density,
          class = "b", coef = nm
        )
      }
      next
    }

    if (identical(b$route, "class_wide")) {
      block_scope <- identical(b$hyperparameters$scope %||% "", "block")
      rows[[length(rows) + 1L]] <- brms::set_prior(
        .fb_pen_class_wide_string(b),
        class = "b",
        nlpar = if (block_scope) .FB_NL_BLOCK_PART else ""
      )
      if (block_scope) {
        nl_block <- list(
          protected_part = .FB_NL_PROTECTED_PART,
          block_part = .FB_NL_BLOCK_PART,
          coef_names = b$coef_names
        )
      }
      next
    }

    if (identical(b$route, "stanvar_density")) {
      compiled <- .fb_pen_stanvar_elastic_net(b, idx, n_train)
      rows[[length(rows) + 1L]] <- compiled$prior
      stanvars <- if (is.null(stanvars)) {
        compiled$stanvars
      } else {
        stanvars + compiled$stanvars
      }
      next
    }

    # A grouped family reaches here only when the Bayesian route was asked
    # for. It has a provider -- the penalised optimiser -- so the refusal is
    # not the withdrawal code, and it names the workflow that does fit it
    # rather than leaving the caller with a route they cannot take.
    if (identical(b$route, "optimizer_only")) {
      stop(.fb_refusal_condition(
        reason_code = "selection_penalty_family_optimizer_only",
        message = paste0(
          "fb_penalty(): `", b$family, "()` is fitted by the penalised ",
          "optimiser and has no brms form. It is a penalty with a proximal ",
          "operator, not a density brms can be given, and no nearby prior ",
          "is substituted under its name. Ask for it with ",
          "workflow = \"penalized_optimizer\", or choose ridge(), ",
          "laplace() or elastic_net(), which the Bayesian route compiles ",
          "exactly."
        ),
        penalty_family = b$family
      ))
    }

    stop(
      ".fb_pen_compile_brms(): block ", b$id, " carries the unknown route `",
      b$route, "`.",
      call. = FALSE
    )
  }

  param_map <- if (length(map_rows) == 0L) {
    data.frame(
      coefficient = character(0), b_index = integer(0), block = integer(0),
      family = character(0), route = character(0),
      hyperparameter_mode = character(0), stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, map_rows)
  }
  param_map <- param_map[order(param_map$b_index), , drop = FALSE]
  rownames(param_map) <- NULL

  routes <- sort(unique(param_map$route[param_map$route != "none"]))
  if (length(routes) == 0L) {
    routes <- "none"
  }

  list(
    prior = if (length(rows) == 0L) NULL else do.call(rbind, rows),
    stanvars = stanvars,
    param_map = param_map,
    route = routes,
    nl_block = nl_block
  )
}


# --- per-route builders ---------------------------------------------------- #

# The two coefficient-row densities. Both are ordinary brms priors on a
# single coefficient, which is what makes `coef_rows` the only route that
# leaves protected coefficients entirely alone.
.fb_pen_coef_density <- function(b) {
  if (identical(b$family, "ridge")) {
    return(sprintf(
      "normal(0, %s)", .fb_pen_num(b$hyperparameters$scale)
    ))
  }
  if (identical(b$family, "laplace")) {
    return(sprintf(
      "double_exponential(0, %s)", .fb_pen_num(b$hyperparameters$scale)
    ))
  }
  stop(
    ".fb_pen_coef_density(): `", b$family, "` has no coefficient-row form.",
    call. = FALSE
  )
}


# The class-wide specials, written in brms's own grammar.
#
# `par_ratio` is brms's calibration of the horseshoe's global scale: it is
# the prior guess at the ratio of non-zero to zero coefficients, so the
# contract's `p0` out of `p` becomes `p0 / (p - p0)`. Passing `par_ratio`
# rather than `scale_global` is deliberate -- brms derives the global scale
# from `par_ratio`, `sigma` and the sample size, which is the Piironen and
# Vehtari calibration and not something to reproduce here.
.fb_pen_class_wide_string <- function(b) {
  h <- b$hyperparameters
  if (identical(b$family, "regularized_horseshoe")) {
    return(sprintf(
      paste0(
        "horseshoe(df = %s, df_global = %s, scale_slab = %s, df_slab = %s, ",
        "par_ratio = %s)"
      ),
      .fb_pen_num(h$df), .fb_pen_num(h$df_global),
      .fb_pen_num(h$slab_scale), .fb_pen_num(h$slab_df),
      .fb_pen_num(h$par_ratio)
    ))
  }
  if (identical(b$family, "r2d2")) {
    return(sprintf(
      "R2D2(mean_R2 = %s, prec_R2 = %s, cons_D2 = %s)",
      .fb_pen_num(h$mean_R2), .fb_pen_num(h$prec_R2), .fb_pen_num(h$cons_D2)
    ))
  }
  stop(
    ".fb_pen_class_wide_string(): `", b$family, "` is not class-wide.",
    call. = FALSE
  )
}


# The elastic net, as an unnormalised `target +=` contribution over an index
# vector carried in `stanvars`.
#
# Two things make this admissible. The index is into brms's `b`, which
# excludes the intercept, so the contribution reaches exactly the candidate
# coefficients and a protected coefficient keeps whatever row `fb_prior()`
# gave it. And the hyperparameters are fixed in this release, so the
# normalising constant -- which depends on `lambda` and `alpha` and on
# nothing else -- is a constant, and an unnormalised density is the same
# posterior. Learned elastic-net hyperparameters would need that constant and
# are a non-goal (SPEC v5 section 15).
.fb_pen_stanvar_elastic_net <- function(b, idx, n_train) {
  h <- b$hyperparameters
  name <- paste0("cand_idx_", b$id)
  m <- length(idx)
  sv <- brms::stanvar(
    x = as.integer(idx),
    name = name,
    block = "data",
    scode = sprintf("array[%d] int %s;", m, name)
  )
  density <- sprintf(
    paste0(
      "target += -%s * %s * (%s * sum(abs(b[%s])) + ",
      "(1 - %s) / 2 * dot_self(b[%s]))"
    ),
    .fb_pen_num(n_train), .fb_pen_num(h$lambda),
    .fb_pen_num(h$alpha), name,
    .fb_pen_num(h$alpha), name
  )
  list(
    prior = brms::set_prior(density, class = "b", check = FALSE),
    stanvars = sv
  )
}
