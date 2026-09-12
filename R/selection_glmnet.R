# The penalised-optimiser provider: `glmnet` on the whitened substrate.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 2, 8.3 and 9,
# implementing the provider contract of
# SPEC_penalised_selection_v3_2026-09-11.md section 8 with the path and
# support semantics of its section 7.
#
# This is `workflow = "penalized_optimizer"`. It is never a Bayesian answer
# and never carries a Bayesian label: `fully_bayesian` is FALSE on every
# object this file returns, and the equivalence level of the comparison is
# stored beside it rather than left to a reader's assumption.
#
# **Nothing about `glmnet` is taken from its documentation.** Its lambda
# conventions are measured here, on a fixed internal design, against a
# Cholesky closed form and a subgradient residual written in
# `R/selection_bridges.R` from the canonical objective. The measurement
# corrected the working hypothesis it started from -- see
# `.fb_sel_glmnet_convention()` -- which is the whole reason the module
# measures rather than reads.

# --- the provider gate -------------------------------------------------- #

# .fb_sel_glmnet_state --- the per-session convention cache.
#
# Keyed by the installed `glmnet` version, so an upgrade inside a session
# re-measures rather than reusing a convention that may no longer hold.
.fb_sel_glmnet_state <- new.env(parent = emptyenv())


# .fb_sel_provider_installed() --- is a provider package available.
#
# One line, and it exists as its own binding so the absent-provider refusal
# has a triggering test that does not require uninstalling anything: a
# registered code nobody has seen fire is a string.
.fb_sel_provider_installed <- function(package) {
  requireNamespace(package, quietly = TRUE)
}


# .fb_sel_require_glmnet() --- the provider and family gate.
#
# `fallback = "none"` is the default of the public verb, so an absent
# provider is a named refusal and never a substitution. The family check
# lives here too: the optimiser workflow serves three families, and a request
# for any of the v3 grammar's group, exclusive, sparse-group or adaptive
# families is refused by name rather than approximated by the nearest thing
# `glmnet` can fit.
.fb_sel_require_glmnet <- function(family = NULL) {
  if (!.fb_sel_provider_installed("glmnet")) {
    stop(.fb_refusal_condition(
      reason_code = "selection_provider_absent",
      message = paste0(
        "fb_select(): the `penalized_optimizer` workflow needs the glmnet ",
        "package, which is not installed. Install it with ",
        "install.packages(\"glmnet\"), or choose a workflow whose provider ",
        "is present. No substitute is fitted: a penalised-likelihood answer ",
        "never appears in place of the one that was asked for."
      ),
      provider = "glmnet",
      workflow = "penalized_optimizer"
    ))
  }
  if (!is.null(family)) {
    .fb_sel_glmnet_assert_family(family)
  }
  invisible(TRUE)
}


# .fb_sel_glmnet_assert_family() --- the three families with a provider.
.fb_sel_glmnet_assert_family <- function(family) {
  if (!is.character(family) || length(family) != 1L || !nzchar(family)) {
    stop("`family` must be a single non-empty string.", call. = FALSE)
  }
  supported <- .fb_sel_penalty_families()
  if (!family %in% supported) {
    stop(.fb_refusal_condition(
      reason_code = "selection_penalty_family_withdrawn",
      message = paste0(
        "fb_select(): the penalty family `",
        family,
        "` has no provider in this release. The optimiser workflow serves ",
        paste(supported, collapse = ", "),
        " only. The group, exclusive, sparse-group and adaptive families of ",
        "the v3 grammar are withdrawn until each ships with its own provider ",
        "and its own gates; their component traits remain in the resolver, ",
        "and a family without a provider is a promise rather than a ",
        "capability."
      ),
      penalty_family = family,
      supported = supported
    ))
  }
  invisible(family)
}


# .fb_sel_glmnet_fit() --- one `glmnet` call at one lambda, with the tight
# convergence settings this module needs.
#
# `thresh` and `maxit` moved into a `control` list in `glmnet` 5.0 and the
# direct arguments now warn. Which spelling the installed version wants is
# read off the callable contract rather than assumed from a version number:
# `formals()` is the interface, a version string is a label for it.
#
# One lambda per call, never a path. `glmnet`'s path fitting stops early on a
# deviance criterion and its `s =` interpolation is not the solution at `s`,
# so a grid walked one point at a time is both the exact route and the one
# whose every point carries its own KKT residual.
.fb_sel_glmnet_fit <- function(x, y, alpha, lambda, thresh, maxit) {
  arguments <- list(
    x = x,
    y = as.numeric(y),
    family = "gaussian",
    alpha = alpha,
    lambda = lambda,
    standardize = FALSE,
    intercept = FALSE
  )
  if ("control" %in% names(formals(glmnet::glmnet))) {
    arguments$control <- list(thresh = thresh, maxit = maxit)
  } else {
    arguments$thresh <- thresh
    arguments$maxit <- maxit
  }
  fit <- do.call(glmnet::glmnet, arguments)
  beta <- as.numeric(stats::coef(fit))[-1L]
  names(beta) <- colnames(x)
  beta
}


# --- the measured convention ------------------------------------------- #

# .fb_sel_glmnet_convention() --- the installed `glmnet`'s lambda maps,
# measured rather than read.
#
# The measurement runs on a fixed, seeded, well-conditioned internal design
# whose response carries a non-zero mean, because the candidate response
# scales differ from one another only when it does. For each candidate rule
# the canonical lambda is mapped to a `glmnet` lambda, `glmnet` is run, and
# the answer is graded against `.fb_sel_bridge_ridge()` for the L2 part and
# against `.fb_sel_kkt_residual()` for the L1 part. The rule with the
# smallest residual is the convention; the residuals of every rule are kept
# so a reader can see the margin between them.
#
# **What the measurement found, and what it corrected.** The L1 map is
# `lambda_glmnet = lambda_canonical / n` exactly. The L2 map is not: under
# `standardize = FALSE, intercept = FALSE` the ridge penalty is applied to a
# response `glmnet` has divided by its own scale, and the recovered rule is
# the UNCENTRED root mean square of the response, not its standard deviation.
# The standard-deviation hypothesis that this work started from leaves a
# residual near 1e-4 relative; the root-mean-square rule leaves the solver's
# own convergence floor. The two differ only when the response has a non-zero
# mean, which is why the internal design gives it one.
#
# A change in `glmnet`'s convention turns this measurement, and the test that
# pins it, red -- rather than moving every answer the module gives by a
# factor nobody looked at.
.fb_sel_glmnet_convention <- function(refresh = FALSE) {
  .fb_sel_require_glmnet()
  version <- as.character(utils::packageVersion("glmnet"))
  key <- paste0("convention_", version)
  if (!isTRUE(refresh) && !is.null(.fb_sel_glmnet_state[[key]])) {
    return(.fb_sel_glmnet_state[[key]])
  }

  thresh <- 1e-16
  maxit <- 1e6
  rules <- c("none", "sd_population", "sd_sample", "rms")

  design <- .fb_sel_with_seed(20260912L, {
    n <- 60L
    p <- 5L
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    colnames(x) <- paste0("m", seq_len(p))
    # A non-zero mean is what separates the root-mean-square rule from the
    # standard-deviation rules; without it the measurement cannot tell them
    # apart and would report whichever came first.
    y <- 2 + drop(x %*% c(1.5, -0.8, 0, 0.4, 0)) + stats::rnorm(n, sd = 0.6)
    list(x = x, y = y, n = n)
  })

  ridge_scale <- 0.8
  lambda_l2 <- 1 / ridge_scale^2
  closed <- .fb_sel_bridge_ridge(design$x, design$y, ridge_scale)
  laplace_scale <- 0.4
  lambda_l1 <- 1 / laplace_scale

  l2_deviation <- numeric(length(rules))
  l2_kkt <- numeric(length(rules))
  l1_kkt <- numeric(length(rules))
  names(l2_deviation) <- rules
  names(l2_kkt) <- rules
  names(l1_kkt) <- rules

  for (rule in rules) {
    scale <- .fb_sel_response_scale(design$y, rule)
    ridge_beta <- .fb_sel_glmnet_fit(
      design$x,
      design$y,
      alpha = 0,
      lambda = lambda_l2 * scale / design$n,
      thresh = thresh,
      maxit = maxit
    )
    l2_deviation[[rule]] <- max(abs(ridge_beta - closed))
    l2_kkt[[rule]] <- .fb_sel_kkt_residual(
      design$x,
      design$y,
      ridge_beta,
      "ridge",
      list(scale = ridge_scale)
    )
    lasso_beta <- .fb_sel_glmnet_fit(
      design$x,
      design$y,
      alpha = 1,
      lambda = lambda_l1 * scale / design$n,
      thresh = thresh,
      maxit = maxit
    )
    l1_kkt[[rule]] <- .fb_sel_kkt_residual(
      design$x,
      design$y,
      lasso_beta,
      "laplace",
      list(scale = laplace_scale)
    )
  }

  l2_rule <- rules[[which.min(l2_deviation)]]
  l1_rule <- rules[[which.min(l1_kkt)]]

  convention <- list(
    l1 = "lambda_canonical / n",
    l2 = paste0("lambda_canonical * ", l2_rule, "(response) / n"),
    l1_scale_rule = l1_rule,
    response_scale_rule = l2_rule,
    response_scale = NA_real_,
    residuals = list(
      l2_max_abs_by_rule = l2_deviation,
      l2_kkt_by_rule = l2_kkt,
      l1_kkt_by_rule = l1_kkt,
      l2_max_abs = unname(l2_deviation[[l2_rule]]),
      l2_kkt = unname(l2_kkt[[l2_rule]]),
      l1_kkt = unname(l1_kkt[[l1_rule]])
    ),
    glmnet_version = version,
    settings = list(
      standardize = FALSE,
      intercept = FALSE,
      thresh = thresh,
      maxit = maxit
    ),
    design = paste0(
      "seeded internal design, ",
      design$n,
      " rows by ",
      ncol(design$x),
      " columns, response mean offset 2, ridge scale ",
      ridge_scale,
      ", laplace scale ",
      laplace_scale
    ),
    ridge_oracle = "closed form via .fb_sel_bridge_ridge()",
    lasso_oracle = "canonical KKT residual via .fb_sel_kkt_residual()"
  )
  .fb_sel_glmnet_state[[key]] <- convention
  convention
}


# .fb_sel_convention_at() --- the convention with the response scale this
# problem's response implies.
#
# The L1 map carries no response scale and the L2 map carries exactly one,
# which is the asymmetry the measurement found. If a future `glmnet` were to
# put a scale on the L1 part as well, the composite map below would no longer
# be the right algebra, so the assumption is asserted here rather than
# assumed silently.
.fb_sel_convention_at <- function(convention, z) {
  if (!is.list(convention) || is.null(convention$response_scale_rule)) {
    stop(
      "`convention` must come from .fb_sel_glmnet_convention().",
      call. = FALSE
    )
  }
  if (!identical(convention$l1_scale_rule, "none")) {
    stop(
      "internal contract violation: the measured L1 map now carries the ",
      "response scale `",
      convention$l1_scale_rule,
      "`, so the composite canonical-to-glmnet map in ",
      ".fb_sel_canonical_to_glmnet() is no longer the right algebra and must ",
      "be re-derived before any path is fitted.",
      call. = FALSE
    )
  }
  convention$response_scale <- .fb_sel_response_scale(
    z,
    convention$response_scale_rule
  )
  convention
}


# --- one canonical point ------------------------------------------------ #

# .fb_sel_glmnet_at() --- the canonical minimiser at one (lambda, alpha),
# through `glmnet` where `glmnet` will take it and through the module's own
# proximal solver where it will not.
#
# Installed `glmnet` refuses a design with fewer than two columns, and a
# zero-penalty point is a least-squares solve rather than a penalised one.
# Both route to `.fb_sel_prox_mode_lambda()`, whose KKT certificate is the
# same one the `glmnet` answer is graded by, so the two routes are
# interchangeable in everything except which code did the arithmetic. v3
# section 8.1's one-column closed form is deliberately not used: it has not
# been independently re-derived and falsified, and that document says it does
# not ship until it has been.
.fb_sel_glmnet_at <- function(
  A,
  z,
  lambda_canonical,
  alpha,
  convention,
  thresh = 1e-16,
  maxit = 1e6
) {
  A <- as.matrix(A)
  z <- as.numeric(z)
  n <- nrow(A)
  p <- ncol(A)
  lambda_l1 <- lambda_canonical * alpha
  lambda_l2 <- lambda_canonical * (1 - alpha)

  if (p == 0L) {
    return(list(
      beta = numeric(0),
      route = "empty_candidate_block",
      lambda_glmnet = NA_real_,
      alpha_glmnet = NA_real_,
      kkt = 0,
      canonical_objective = 0.5 * sum(z^2),
      provider_objective = NA_real_
    ))
  }

  map <- .fb_sel_canonical_to_glmnet(
    lambda_canonical,
    alpha,
    n,
    convention
  )

  if (p < 2L || lambda_canonical <= 0) {
    solved <- .fb_sel_prox_mode_lambda(
      A,
      z,
      lambda_l1,
      lambda_l2,
      tol = 1e-9,
      max_iter = 2e5
    )
    beta <- solved$beta
    route <- "internal_prox"
  } else {
    beta <- .fb_sel_glmnet_fit(
      A,
      z,
      map$alpha,
      map$lambda,
      thresh = thresh,
      maxit = maxit
    )
    route <- "glmnet"
  }

  list(
    beta = beta,
    route = route,
    lambda_glmnet = map$lambda,
    alpha_glmnet = map$alpha,
    kkt = .fb_sel_kkt_residual_lambda(A, z, beta, lambda_l1, lambda_l2),
    canonical_objective = .fb_sel_canonical_objective(
      A,
      z,
      beta,
      lambda_l1,
      lambda_l2
    ),
    # `glmnet`'s own documented objective at the same point. It is recorded
    # beside the canonical one because the measured L2 convention says the
    # two are not the same function, and a certificate that showed only one
    # of them would hide exactly that.
    provider_objective = (0.5 / n) *
      sum((z - drop(A %*% beta))^2) +
      map$lambda *
        (map$alpha * sum(abs(beta)) + 0.5 * (1 - map$alpha) * sum(beta^2))
  )
}


# --- path and support semantics ----------------------------------------- #

# .fb_sel_support_semantics() --- v3 section 7's row for a family.
#
# Ridge does not select. It reaches zero only as lambda goes to infinity, so
# it has no finite null anchor and no support to report, and a ridge path
# that came back with a "selected set" would be reporting an artefact of
# where the grid happened to stop.
.fb_sel_support_semantics <- function(family) {
  .fb_sel_glmnet_assert_family(family)
  switch(
    family,
    ridge = list(
      family = "ridge",
      selects = FALSE,
      support_semantics = "none",
      path_parameter = "log_lambda",
      null_anchor = "infinite",
      selection_resolution = NA_character_
    ),
    laplace = list(
      family = "laplace",
      selects = TRUE,
      support_semantics = "coefficient",
      path_parameter = "lambda",
      null_anchor = "finite",
      selection_resolution = "coefficient"
    ),
    elastic_net = list(
      family = "elastic_net",
      selects = TRUE,
      support_semantics = "coefficient",
      path_parameter = "lambda",
      null_anchor = "finite",
      selection_resolution = "coefficient"
    )
  )
}


# .fb_sel_glmnet_decision() --- the decision settings, with the contract's
# defaults filled in.
.fb_sel_glmnet_decision <- function(decision = NULL) {
  defaults <- list(
    rule = "one_se",
    metric = "rmse",
    tolerance = NULL,
    benchmark = TRUE,
    n_lambda = 20L,
    lambda_min_ratio = 1e-3
  )
  if (is.null(decision)) {
    return(defaults)
  }
  if (!is.list(decision)) {
    stop("`decision` must be a list or NULL.", call. = FALSE)
  }
  unknown <- setdiff(names(decision), names(defaults))
  if (length(unknown) > 0L) {
    stop(
      "unknown decision setting(s): ",
      paste(unknown, collapse = ", "),
      ". The settings are: ",
      paste(names(defaults), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  out <- utils::modifyList(defaults, decision)
  out$rule <- match.arg(out$rule, c("one_se", "min", "tolerance"))
  out$metric <- match.arg(out$metric, c("rmse", "mse", "nlpd"))
  out$n_lambda <- as.integer(out$n_lambda)
  if (is.na(out$n_lambda) || out$n_lambda < 2L) {
    stop("`n_lambda` must be at least two grid points.", call. = FALSE)
  }
  out
}


# --- block normalisation ------------------------------------------------ #

# .fb_sel_glmnet_blocks() --- the resolved penalty blocks as one canonical
# penalty over one column set.
#
# `glmnet`'s `alpha` is a scalar and its `penalty.factor` is renormalised to
# sum to the column count, so two blocks with different families -- or the
# same family at different hyperparameters -- are not the problem its printed
# lambda names. The Bayesian workflow carries block-specific priors exactly;
# the optimiser benchmark refuses rather than fitting a nearby problem under
# the requested problem's label.
.fb_sel_glmnet_blocks <- function(blocks, candidate_names, n) {
  if (!is.list(blocks) || length(blocks) == 0L) {
    stop(
      "`blocks` must be a non-empty list of resolved penalty blocks.",
      call. = FALSE
    )
  }
  if (!is.null(names(blocks)) && any(names(blocks) == "family")) {
    blocks <- list(blocks)
  }

  resolved <- lapply(seq_along(blocks), function(index) {
    block <- blocks[[index]]
    if (!is.list(block) || is.null(block$family)) {
      stop(
        "block ",
        index,
        " must be a list carrying `family`, `columns` and ",
        "`hyperparameters`.",
        call. = FALSE
      )
    }
    .fb_sel_glmnet_assert_family(block$family)
    columns <- block$columns
    if (is.character(columns)) {
      missing_columns <- setdiff(columns, candidate_names)
      if (length(missing_columns) > 0L) {
        stop(
          "block ",
          index,
          " names candidate column(s) the design does not ",
          "carry: ",
          paste(missing_columns, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
      columns <- match(columns, candidate_names)
    }
    columns <- as.integer(columns)
    if (
      length(columns) == 0L ||
        anyNA(columns) ||
        any(columns < 1L | columns > length(candidate_names))
    ) {
      stop(
        "block ",
        index,
        " must index at least one candidate column inside ",
        "the design.",
        call. = FALSE
      )
    }
    list(
      family = block$family,
      columns = columns,
      column_names = candidate_names[columns],
      hyperparameters = block$hyperparameters
    )
  })

  all_columns <- unlist(lapply(resolved, function(b) b$columns))
  if (anyDuplicated(all_columns) > 0L) {
    repeated <- candidate_names[unique(all_columns[duplicated(all_columns)])]
    stop(
      "internal contract violation: candidate column(s) in more than one ",
      "penalty block: ",
      paste(repeated, collapse = ", "),
      ". Overlap is resolved before a provider is called.",
      call. = FALSE
    )
  }

  families <- vapply(resolved, function(b) b$family, character(1L))
  penalties <- lapply(resolved, function(b) {
    .fb_sel_canonical_penalty(b$family, b$hyperparameters, n = n)
  })
  distinct <- unique(vapply(
    penalties,
    function(p) paste0(p$family, ":", p$lambda_l1, ":", p$lambda_l2),
    character(1L)
  ))
  if (length(distinct) > 1L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_penalty_blocks_not_expressible",
      message = paste0(
        "fb_select(): the resolved penalty blocks name ",
        length(distinct),
        " distinct penalties (",
        paste(distinct, collapse = "; "),
        "), and `glmnet` fits one. Its `alpha` is a scalar and its ",
        "`penalty.factor` is renormalised to sum to the column count, so a ",
        "second penalty is not the problem its printed lambda names. Run the ",
        "blocks as one family at one set of hyperparameters, or take the ",
        "Bayesian workflow, which carries block-specific priors exactly."
      ),
      penalties = distinct,
      families = unique(families)
    ))
  }

  list(
    blocks = resolved,
    family = families[[1L]],
    columns = sort(all_columns),
    column_names = candidate_names[sort(all_columns)],
    penalty = penalties[[1L]]
  )
}


# --- the fold reduction ------------------------------------------------- #

# .fb_sel_glmnet_reduction() --- one fold's whitened, protected-residualised
# provider payload.
#
# The training rows are whitened by the fold's own Cholesky factor of V_TT --
# the same factor `.fb_sel_gls_fit()` uses, so the penalised path and the
# unpenalised benchmark sit in the same metric -- and the protected block is
# residualised out by Frisch-Waugh-Lovell before the provider sees anything.
# Candidate columns no block names are unpenalised by construction and join
# the protected block, which is the `term(...) ~ none()` role of section 3.2
# expressed where it belongs: in the reduction, not in a penalty factor.
.fb_sel_glmnet_reduction <- function(substrate, penalised_columns) {
  if (!inherits(substrate, "fb_sel_substrate")) {
    stop("`substrate` must be an <fb_sel_substrate>.", call. = FALSE)
  }
  train_rows <- substrate$train_rows
  factor_lower <- t(substrate$blocks$factor_tt)

  candidate_names <- colnames(substrate$x_candidate)
  free_columns <- setdiff(
    seq_len(ncol(substrate$x_candidate)),
    penalised_columns
  )
  x_protected <- cbind(
    substrate$x_protected,
    substrate$x_candidate[, free_columns, drop = FALSE]
  )
  x_penalised <- substrate$x_candidate[, penalised_columns, drop = FALSE]

  whiten <- function(m) forwardsolve(factor_lower, m)
  reduction <- .fb_sel_protected_reduction(
    y = drop(whiten(matrix(
      substrate$response[train_rows],
      ncol = 1L
    ))),
    X0 = .fb_sel_glmnet_named(
      whiten(x_protected[train_rows, , drop = FALSE]),
      colnames(x_protected)
    ),
    Xp = .fb_sel_glmnet_named(
      whiten(x_penalised[train_rows, , drop = FALSE]),
      colnames(x_penalised)
    ),
    whitener = NULL,
    on_rank_deficient = "report"
  )

  list(
    reduction = reduction,
    payload = .fb_sel_provider_payload(reduction),
    x_full = cbind(substrate$x_protected, substrate$x_candidate),
    protected_names = colnames(x_protected),
    penalised_names = candidate_names[penalised_columns],
    free_names = candidate_names[free_columns]
  )
}


# .fb_sel_glmnet_named() --- restore the column names a triangular solve
# drops, so the reduction can hold the provider to them.
.fb_sel_glmnet_named <- function(m, nms) {
  m <- as.matrix(m)
  if (ncol(m) > 0L) {
    colnames(m) <- nms
  }
  m
}


# .fb_sel_glmnet_fold_solution() --- one fold at one canonical point, from
# the provider's coefficients to the assessment-row predictions.
.fb_sel_glmnet_fold_solution <- function(
  substrate,
  parts,
  lambda_canonical,
  alpha,
  convention
) {
  payload <- parts$payload
  solved <- .fb_sel_glmnet_at(
    payload$x,
    payload$y,
    lambda_canonical,
    alpha,
    convention
  )
  beta_penalised <- solved$beta

  recovered <- .fb_sel_recover_protected(parts$reduction, beta_penalised)
  full_names <- colnames(parts$x_full)
  beta_full <- stats::setNames(rep(0, length(full_names)), full_names)
  if (length(recovered) > 0L) {
    recovered[is.na(recovered)] <- 0
    beta_full[names(recovered)] <- recovered
  }
  if (length(beta_penalised) > 0L) {
    beta_full[parts$penalised_names] <- beta_penalised
  }

  train_rows <- substrate$train_rows
  assess_rows <- substrate$assess_rows
  mu <- .fb_sel_fold_predict(
    substrate$blocks,
    beta_full,
    parts$x_full[train_rows, , drop = FALSE],
    substrate$response[train_rows],
    parts$x_full[assess_rows, , drop = FALSE]
  )

  list(
    beta_penalised = beta_penalised,
    beta_protected = recovered,
    beta_full = beta_full,
    mu_assess = mu,
    route = solved$route,
    kkt = solved$kkt,
    lambda_glmnet = solved$lambda_glmnet,
    alpha_glmnet = solved$alpha_glmnet,
    canonical_objective = solved$canonical_objective,
    provider_objective = solved$provider_objective
  )
}


# --- the lambda grid ---------------------------------------------------- #

# .fb_sel_glmnet_grid() --- the family-canonical path coordinate, fixed once
# before the folds.
#
# A grid rebuilt inside each fold would make "the same lambda" mean a
# different problem on each of them, and the scores would be averaged across
# points that are not the same point. For a family with an L1 part the
# anchor is the smallest canonical lambda at which every coefficient is zero,
# taken as the largest such value across folds; for ridge, which has no
# finite null anchor (v3 section 7), the grid is decades either side of the
# prior's own canonical lambda.
#
# The anchor reads training rows only. It is computed from each fold's
# residualised payload, which is built from that fold's training rows alone,
# so a common grid is not a common leak.
.fb_sel_glmnet_grid <- function(parts_by_fold, penalty, decision) {
  alpha <- penalty$alpha_composite
  anchor <- penalty$lambda_composite
  if (alpha > 0) {
    lambda_max <- max(vapply(
      parts_by_fold,
      function(parts) {
        payload <- parts$payload
        if (ncol(payload$x) == 0L) {
          return(0)
        }
        max(abs(drop(crossprod(payload$x, payload$y)))) / alpha
      },
      numeric(1L)
    ))
    if (!is.finite(lambda_max) || lambda_max <= 0) {
      lambda_max <- anchor
    }
    grid <- exp(seq(
      log(lambda_max),
      log(lambda_max * decision$lambda_min_ratio),
      length.out = decision$n_lambda
    ))
  } else {
    grid <- anchor * 10^seq(2, -2, length.out = decision$n_lambda)
  }
  sort(unique(c(grid, anchor)), decreasing = TRUE)
}


# --- the provider path -------------------------------------------------- #

# .fb_sel_glmnet_path() --- the v3 section 8 provider contract, on the
# substrate this module already builds.
#
# Two arms come back. The **fixed** arm is the single canonical point the
# declared hyperparameters name: it is the point at which a mode identity
# with a Bayesian fit can hold, and the one the equivalence level is about.
# The **benchmark** arm walks the family-canonical lambda grid and applies
# the one-standard-error decision at the fold unit, which is what a
# `glmnet` user would do and therefore what a head-to-head comparison has to
# put on the table.
#
# Every returned quantity that names an estimand, a provider, a convention or
# a route is copied from the object that did the work. `fully_bayesian` is
# FALSE, and it is a constant rather than an argument.
.fb_sel_glmnet_path <- function(
  substrate,
  blocks,
  folds = NULL,
  decision = NULL,
  seed = NULL,
  cluster = NULL
) {
  .fb_sel_require_glmnet()
  decision <- .fb_sel_glmnet_decision(decision)

  substrates <- if (inherits(substrate, "fb_sel_substrate")) {
    list(substrate)
  } else {
    substrate
  }
  if (
    length(substrates) == 0L ||
      !all(vapply(
        substrates,
        inherits,
        logical(1L),
        what = "fb_sel_substrate"
      ))
  ) {
    stop(
      "`substrate` must be an <fb_sel_substrate> or a list of them.",
      call. = FALSE
    )
  }
  labels <- vapply(substrates, function(s) s$label, character(1L))
  names(substrates) <- labels
  if (!is.null(folds)) {
    fold_labels <- vapply(folds, function(f) f$label, character(1L))
    if (!setequal(fold_labels, labels)) {
      stop(
        "`folds` and `substrate` must describe the same folds: folds carry ",
        paste(fold_labels, collapse = ", "),
        " and the substrates carry ",
        paste(labels, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }

  candidate_names <- colnames(substrates[[1L]]$x_candidate)
  n_train <- length(substrates[[1L]]$train_rows)
  resolved <- .fb_sel_glmnet_blocks(blocks, candidate_names, n = n_train)
  penalty <- resolved$penalty
  semantics <- .fb_sel_support_semantics(resolved$family)
  convention_base <- .fb_sel_glmnet_convention()

  run <- function() {
    parts_by_fold <- lapply(substrates, function(s) {
      .fb_sel_glmnet_reduction(s, resolved$columns)
    })
    conventions <- lapply(parts_by_fold, function(parts) {
      .fb_sel_convention_at(convention_base, parts$payload$y)
    })

    grid <- .fb_sel_glmnet_grid(parts_by_fold, penalty, decision)
    arms <- list(fixed = penalty$lambda_composite)
    if (isTRUE(decision$benchmark)) {
      arms$benchmark <- grid
    }

    coefficients <- list()
    fitted <- list()
    canonical <- list()
    provider <- list()
    kkt <- list()
    coordinates <- list()
    supports <- list()
    # Kept per model rather than as one table, for the reason
    # `.fb_sel_cv()` gives: flexyBayes is not a data.table-aware namespace,
    # so the module collects and binds rather than subsetting a data.table
    # it would then be reading under data.frame semantics.
    clusters_by_model <- list()
    supports_by_model <- list()

    for (arm in names(arms)) {
      lambdas <- arms[[arm]]
      for (index in seq_along(lambdas)) {
        lambda_canonical <- lambdas[[index]]
        model <- paste0(arm, "_", index)
        for (label in labels) {
          s <- substrates[[label]]
          parts <- parts_by_fold[[label]]
          solution <- .fb_sel_glmnet_fold_solution(
            s,
            parts,
            lambda_canonical,
            penalty$alpha_composite,
            conventions[[label]]
          )
          scale <- s$learned$preprocessing_scale
          estimates <- solution$beta_penalised
          if (length(estimates) > 0L) {
            estimates <- estimates / scale[parts$penalised_names]
          }
          coefficients[[length(coefficients) + 1L]] <-
            data.table::data.table(
              arm = arm,
              lambda_index = index,
              lambda_canonical = lambda_canonical,
              fold = label,
              role = c(
                rep("candidate", length(estimates)),
                rep("protected", length(solution$beta_protected))
              ),
              term = c(
                parts$penalised_names,
                names(solution$beta_protected)
              ),
              estimate = c(
                as.numeric(estimates),
                as.numeric(solution$beta_protected)
              )
            )
          fitted[[length(fitted) + 1L]] <- data.table::data.table(
            arm = arm,
            lambda_index = index,
            lambda_canonical = lambda_canonical,
            fold = label,
            row = s$assess_rows,
            fitted = as.numeric(solution$mu_assess)
          )
          canonical[[length(canonical) + 1L]] <- data.table::data.table(
            arm = arm,
            lambda_index = index,
            lambda_canonical = lambda_canonical,
            fold = label,
            canonical_objective = solution$canonical_objective,
            route = solution$route
          )
          provider[[length(provider) + 1L]] <- data.table::data.table(
            arm = arm,
            lambda_index = index,
            lambda_canonical = lambda_canonical,
            fold = label,
            provider_objective = solution$provider_objective,
            lambda_glmnet = solution$lambda_glmnet,
            alpha_glmnet = solution$alpha_glmnet,
            route = solution$route
          )
          kkt[[length(kkt) + 1L]] <- data.table::data.table(
            arm = arm,
            lambda_index = index,
            lambda_canonical = lambda_canonical,
            fold = label,
            kkt = solution$kkt
          )
          coordinates[[length(coordinates) + 1L]] <-
            data.table::data.table(
              arm = arm,
              lambda_index = index,
              lambda_canonical = lambda_canonical,
              alpha_canonical = penalty$alpha_composite,
              lambda_glmnet = solution$lambda_glmnet,
              alpha_glmnet = solution$alpha_glmnet,
              path_parameter = semantics$path_parameter,
              null_anchor = semantics$null_anchor,
              fold = label
            )
          nonzero <- names(
            solution$beta_penalised
          )[abs(solution$beta_penalised) > 0]
          support_row <- data.table::data.table(
            arm = arm,
            lambda_index = index,
            lambda_canonical = lambda_canonical,
            fold = label,
            n_nonzero = length(nonzero),
            effective_df = .fb_sel_glmnet_effective_df(
              parts$payload$x,
              lambda_canonical,
              penalty$alpha_composite
            ),
            terms = paste(nonzero, collapse = ", ")
          )
          supports[[length(supports) + 1L]] <- support_row
          supports_by_model[[model]] <- c(
            supports_by_model[[model]],
            list(support_row)
          )
          key <- if (is.null(cluster)) {
            rep(label, length(s$assess_rows))
          } else {
            as.character(cluster)[s$assess_rows]
          }
          scored <- .fb_sel_score_fold(
            s$assessment$y,
            solution$mu_assess,
            s$assessment$v_s_given_t,
            key,
            label = label,
            repeat_id = s$repeat_id,
            model = model
          )
          clusters_by_model[[model]] <- c(
            clusters_by_model[[model]],
            list(scored$clusters)
          )
        }
      }
    }

    list(
      coefficients = data.table::rbindlist(coefficients),
      fitted = data.table::rbindlist(fitted),
      canonical_objective = data.table::rbindlist(canonical),
      provider_objective = data.table::rbindlist(provider),
      kkt = data.table::rbindlist(kkt),
      coordinates = data.table::rbindlist(coordinates),
      supports = data.table::rbindlist(supports),
      clusters_by_model = lapply(clusters_by_model, data.table::rbindlist),
      supports_by_model = lapply(supports_by_model, data.table::rbindlist),
      grid = grid,
      arms = arms,
      conventions = conventions
    )
  }

  computed <- if (is.null(seed)) run() else .fb_sel_with_seed(seed, run())

  selection <- .fb_sel_glmnet_decide(
    computed,
    decision,
    semantics,
    labels
  )

  equivalence <- .fb_sel_equivalence_level(
    prior_family = resolved$family,
    provider = "glmnet",
    # The fold plugs its covariance in and the whitened objective carries a
    # unit residual scale by construction, so both are held fixed here. This
    # is read off what the substrate did, not off an argument.
    covariance_fixed = TRUE,
    sigma_fixed = TRUE
  )

  structure(
    list(
      workflow = "penalized_optimizer",
      provider = "glmnet",
      fully_bayesian = FALSE,
      family = resolved$family,
      blocks = resolved$blocks,
      penalty = penalty,
      coefficients_original_scale = computed$coefficients,
      fitted_values = computed$fitted,
      support_by_resolution = list(
        semantics = semantics$support_semantics,
        selection_resolution = semantics$selection_resolution,
        selects = semantics$selects,
        by_lambda = computed$supports,
        selected = selection$selected_terms
      ),
      canonical_objective = computed$canonical_objective,
      provider_objective = computed$provider_objective,
      kkt_by_lambda = computed$kkt,
      path_coordinates = computed$coordinates,
      provider_version = as.character(utils::packageVersion("glmnet")),
      provider_lambda_convention = convention_base,
      equivalence_level = equivalence,
      decision = decision,
      selection = selection,
      covariance_estimand = substrates[[1L]]$covariance_estimand,
      prediction = substrates[[1L]]$prediction,
      folds = labels,
      seed = seed,
      contract = "SPEC_bayesian_selection_v5_2026-09-12"
    ),
    class = "fb_sel_canonical_path"
  )
}


# .fb_sel_glmnet_effective_df() --- the ridge degrees of freedom at a
# canonical lambda, which is the complexity coordinate for a family that does
# not select.
#
# sum_j d_j^2 / (d_j^2 + l2) over the singular values of the residualised
# candidate block. With an L1 part present the support size is the complexity
# the decision rule uses instead, and this is reported alongside rather than
# used.
.fb_sel_glmnet_effective_df <- function(A, lambda_canonical, alpha) {
  if (ncol(A) == 0L) {
    return(0)
  }
  lambda_l2 <- lambda_canonical * (1 - alpha)
  d2 <- svd(A, nu = 0L, nv = 0L)$d^2
  sum(d2 / (d2 + lambda_l2))
}


# .fb_sel_glmnet_decide() --- the out-of-fold decision over the benchmark
# arm.
#
# The scoring unit is the fold unless the caller names a cluster key, and the
# standard error is taken at that unit, per v3 section 9: a row-level
# standard error on correlated data understates it and collapses the
# one-standard-error rule to the minimum.
#
# The complexity a ridge path is ordered by is its effective degrees of
# freedom, not a support size, because a ridge path has no support (v3
# section 7). The selected set it reports is empty for the same reason.
.fb_sel_glmnet_decide <- function(computed, decision, semantics, labels) {
  if (!"benchmark" %in% names(computed$arms)) {
    return(list(
      ran = FALSE,
      reason = "the benchmark arm was not requested",
      selected_terms = character(0)
    ))
  }
  models <- grep("^benchmark_", names(computed$clusters_by_model), value = TRUE)
  indices <- sort(as.integer(sub("^benchmark_", "", models)))

  rows <- lapply(indices, function(index) {
    model <- paste0("benchmark_", index)
    summary <- .fb_sel_aggregate_scores(
      computed$clusters_by_model[[model]],
      metric = decision$metric,
      unit = "cluster"
    )
    support_rows <- computed$supports_by_model[[model]]
    complexity <- if (isTRUE(semantics$selects)) {
      mean(support_rows$n_nonzero)
    } else {
      mean(support_rows$effective_df)
    }
    data.table::data.table(
      model = model,
      lambda_index = index,
      lambda_canonical = support_rows$lambda_canonical[[1L]],
      score = summary$score,
      se = summary$se,
      n_units = summary$n_units,
      complexity = complexity
    )
  })
  scores <- data.table::rbindlist(rows)

  chosen <- .fb_sel_select_model(
    scores,
    rule = decision$rule,
    tolerance = decision$tolerance
  )
  chosen_index <- scores$lambda_index[[chosen$index]]
  selected_terms <- if (isTRUE(semantics$selects)) {
    terms <- computed$supports_by_model[[
      paste0("benchmark_", chosen_index)
    ]]$terms
    parts <- unlist(strsplit(terms, ", ", fixed = TRUE))
    sort(unique(parts[nzchar(parts)]))
  } else {
    character(0)
  }

  list(
    ran = TRUE,
    scores = scores,
    rule = chosen$rule,
    metric = decision$metric,
    unit = "cluster",
    selected_model = chosen$selected,
    selected_lambda_index = chosen_index,
    selected_lambda_canonical = scores$lambda_canonical[[chosen$index]],
    minimum_model = chosen$minimum,
    threshold = chosen$threshold,
    selected_terms = selected_terms,
    support_semantics = semantics$support_semantics
  )
}
