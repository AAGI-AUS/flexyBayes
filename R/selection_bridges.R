# Exact-limit bridges between a shrinkage prior's posterior mode and a
# penalised optimiser.
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 4.3, 8.1 and
# 8.2, on the whitened substrate of SPEC_penalised_selection_v3_2026-09-11.md
# section 5.1.
#
# The canonical objective, and the only one this file knows about, is
#
#   0.5 * || z - A beta ||^2  +  pen(beta)
#
# on the whitened, protected-residualised scale. Every family this release
# ships writes its penalty as
#
#   pen(beta) = lambda_l1 * ||beta||_1  +  (lambda_l2 / 2) * ||beta||_2^2
#
# so one pair of numbers carries ridge, laplace and elastic net, and one
# Karush-Kuhn-Tucker condition certifies all three.
#
# Two of the oracles here are deliberately independent of any solver. The
# ridge bridge is a Cholesky solve of the normal equations; the KKT residual
# is a subgradient measurement at whatever point it is handed. Neither reads
# a provider's convergence flag, and neither shares a source with the code
# whose answer it grades: probe B4 found Stan's L-BFGS reporting convergence
# 0.16 away from a lasso optimum, which is exactly the failure a residual
# measured at the returned point catches and a status flag does not.

# --- the canonical penalty -------------------------------------------- #

# .fb_sel_penalty_families() --- the three families with a mode identity that
# this release ships to an optimiser.
.fb_sel_penalty_families <- function() {
  c(
    "ridge", "laplace", "elastic_net",
    "group_lasso", "sparse_group_lasso", "exclusive_lasso", "group_exclusive"
  )
}


# .fb_sel_canonical_penalty() --- a family and its hyperparameters as the
# (l1, l2) pair of the canonical objective.
#
# Section 4.3 fixes each map:
#
#   ridge(scale = s)              Normal(0, s^2)          -> l2 = 1 / s^2
#   laplace(scale = b)            DoubleExponential(0, b) -> l1 = 1 / b
#   elastic_net(lambda, alpha)    exp{-n lambda [alpha |beta| +
#                                   (1 - alpha)/2 beta^2]}
#                                 -> l1 = n lambda alpha,
#                                    l2 = n lambda (1 - alpha)
#
# The `n` in the elastic-net density is the number of rows entering the
# objective, so `lambda` for that family is already on the per-observation
# scale a penalised optimiser speaks; ridge's and laplace's are prior scales
# and are not.
.fb_sel_canonical_penalty <- function(family, hyper, n) {
  family <- match.arg(family, .fb_sel_penalty_families())
  if (!is.list(hyper)) {
    stop("`hyper` must be a named list of hyperparameters.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || n <= 0) {
    stop("`n` must be a single positive number of rows.", call. = FALSE)
  }

  positive <- function(value, what) {
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        !is.finite(value) ||
        value <= 0
    ) {
      stop(
        "`",
        what,
        "` must be a single finite positive number.",
        call. = FALSE
      )
    }
    as.numeric(value)
  }

  unit_open <- function(value, what) {
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        !is.finite(value) ||
        value <= 0 ||
        value >= 1
    ) {
      stop(
        "`",
        what,
        "` must be a single number strictly inside (0, 1).",
        call. = FALSE
      )
    }
    as.numeric(value)
  }

  parts <- switch(
    family,
    ridge = {
      scale <- positive(hyper$scale, "scale")
      c(l1 = 0, l2 = 1 / scale^2, lg = 0, le = 0)
    },
    laplace = {
      scale <- positive(hyper$scale, "scale")
      c(l1 = 1 / scale, l2 = 0, lg = 0, le = 0)
    },
    elastic_net = {
      lambda <- positive(hyper$lambda, "lambda")
      alpha <- unit_open(hyper$alpha, "alpha")
      c(
        l1 = n * lambda * alpha, l2 = n * lambda * (1 - alpha), lg = 0, le = 0
      )
    },
    # The grouped families. Each states its scale the way the ungrouped
    # family with the same homogeneity does, so the exact limits hold and no
    # new convention is invented: a group norm is homogeneous of degree one
    # like a Laplace scale (lg = 1 / scale), and the exclusive penalty is
    # homogeneous of degree two like a ridge scale (le = 1 / scale^2). That
    # second choice is what makes a singleton exclusive group EXACTLY a
    # ridge, which is the identity the solver tests pin.
    group_lasso = {
      scale <- positive(hyper$scale, "scale")
      c(l1 = 0, l2 = 0, lg = 1 / scale, le = 0)
    },
    sparse_group_lasso = {
      scale <- positive(hyper$scale, "scale")
      alpha <- unit_open(hyper$alpha, "alpha")
      c(l1 = alpha / scale, l2 = 0, lg = (1 - alpha) / scale, le = 0)
    },
    exclusive_lasso = {
      scale <- positive(hyper$scale, "scale")
      c(l1 = 0, l2 = 0, lg = 0, le = 1 / scale^2)
    },
    group_exclusive = {
      group_scale <- positive(hyper$group_scale, "group_scale")
      exclusive_scale <- positive(hyper$exclusive_scale, "exclusive_scale")
      c(l1 = 0, l2 = 0, lg = 1 / group_scale, le = 1 / exclusive_scale^2)
    }
  )

  grouped <- family %in% .fb_sel_grouped_families()
  composite <- parts[["l1"]] + parts[["l2"]]
  # A grouped penalty is not a point on the (lambda, alpha) line the provider
  # map speaks, so it carries a multiplier instead: the path scales all four
  # numbers together and the declared hyperparameters are the point at one.
  # `finite_anchor` says whether the family has a null at a finite multiplier
  # at all -- the exclusive penalty is smooth at the origin and has none, for
  # the same reason ridge has none.
  list(
    family = family,
    lambda_l1 = parts[["l1"]],
    lambda_l2 = parts[["l2"]],
    lambda_group = parts[["lg"]],
    lambda_exclusive = parts[["le"]],
    grouped = grouped,
    group_rule = if (grouped) hyper$group %||% "term" else NA_character_,
    finite_anchor = if (grouped) {
      parts[["l1"]] > 0 || parts[["lg"]] > 0
    } else {
      parts[["l1"]] > 0
    },
    # The composite coordinate the provider map speaks: pen(beta) =
    # lambda_composite [ alpha ||beta||_1 + (1 - alpha)/2 ||beta||_2^2 ].
    # A grouped penalty has no such coordinate; its path multiplier is one.
    lambda_composite = if (grouped) 1 else composite,
    alpha_composite = if (grouped) {
      NA_real_
    } else if (composite > 0) {
      parts[["l1"]] / composite
    } else {
      0
    },
    n = as.numeric(n),
    hyperparameters = hyper
  )
}


# .fb_sel_penalty_value() --- pen(beta) at a point.
.fb_sel_penalty_value <- function(beta, lambda_l1, lambda_l2) {
  lambda_l1 * sum(abs(beta)) + 0.5 * lambda_l2 * sum(beta^2)
}


# .fb_sel_canonical_objective() --- the objective this module minimises.
.fb_sel_canonical_objective <- function(A, z, beta, lambda_l1, lambda_l2) {
  residual <- as.numeric(z) - drop(A %*% beta)
  0.5 * sum(residual^2) + .fb_sel_penalty_value(beta, lambda_l1, lambda_l2)
}


# --- the ridge bridge -------------------------------------------------- #

# .fb_sel_bridge_ridge() --- the exact ridge mode, by Cholesky.
#
#   beta-hat = (A'A + I / scale^2)^-1 A'z
#
# written as one Cholesky factorisation and two triangular solves. The
# inverse is never formed: `solve()` on the cross-product would be the
# shorter line and the less accurate one, and at any width worth the name it
# is also the slower.
#
# This is the ridge exact-limit oracle of section 4.3, and it is deliberately
# independent of `glmnet`: the measured L2 convention (section 8 of
# `.fb_sel_glmnet_convention()`) is graded against this, not the other way
# round.
.fb_sel_bridge_ridge <- function(A, z, scale) {
  A <- as.matrix(A)
  z <- as.numeric(z)
  if (nrow(A) != length(z)) {
    stop(
      "`A` and `z` must agree on the number of rows: ",
      nrow(A),
      " and ",
      length(z),
      ".",
      call. = FALSE
    )
  }
  if (
    !is.numeric(scale) ||
      length(scale) != 1L ||
      !is.finite(scale) ||
      scale <= 0
  ) {
    stop("`scale` must be a single finite positive number.", call. = FALSE)
  }
  p <- ncol(A)
  if (p == 0L) {
    return(numeric(0))
  }

  lambda_l2 <- 1 / scale^2
  gram <- crossprod(A)
  diag(gram) <- diag(gram) + lambda_l2
  factor_gram <- chol(gram)
  out <- drop(backsolve(
    factor_gram,
    forwardsolve(t(factor_gram), crossprod(A, z))
  ))
  names(out) <- colnames(A)
  out
}


# --- the KKT residual -------------------------------------------------- #

# .fb_sel_kkt_residual() --- how far a point is from the canonical optimum,
# in the max norm of the subdifferential.
#
# For the smooth part g(beta) = 0.5 ||z - A beta||^2 + (l2/2) ||beta||^2 the
# gradient is A'(A beta - z) + l2 beta. Optimality of the full objective is
#
#   beta_j != 0 :  grad_j + l1 sign(beta_j) = 0
#   beta_j == 0 :  |grad_j| <= l1
#
# so the residual is |grad_j + l1 sign(beta_j)| on the active coordinates and
# max(|grad_j| - l1, 0) on the inactive ones. With l1 = 0 it collapses to the
# gradient norm, which is the right thing for ridge.
#
# The measurement is of the point, not of the procedure that produced it. A
# solver that stops early, a provider whose lambda means something else, and
# a closed form with a transcription error all show up here; a convergence
# flag catches none of the three.
.fb_sel_kkt_residual <- function(
  A,
  z,
  beta,
  family = c("ridge", "laplace", "elastic_net"),
  hyper
) {
  family <- match.arg(family)
  A <- as.matrix(A)
  z <- as.numeric(z)
  beta <- as.numeric(beta)
  if (nrow(A) != length(z)) {
    stop(
      "`A` and `z` must agree on the number of rows: ",
      nrow(A),
      " and ",
      length(z),
      ".",
      call. = FALSE
    )
  }
  if (ncol(A) != length(beta)) {
    stop(
      "`beta` must carry one coefficient per column of `A`: ",
      ncol(A),
      " expected, ",
      length(beta),
      " supplied.",
      call. = FALSE
    )
  }
  if (ncol(A) == 0L) {
    return(0)
  }

  penalty <- .fb_sel_canonical_penalty(family, hyper, n = nrow(A))
  .fb_sel_kkt_residual_lambda(
    A,
    z,
    beta,
    penalty$lambda_l1,
    penalty$lambda_l2
  )
}


# .fb_sel_kkt_residual_lambda() --- the same measurement, taken directly in
# the (l1, l2) coordinate.
#
# The path arm walks a lambda grid that belongs to no single prior, so it
# needs the residual without going through a family's hyperparameters. Both
# entry points share this body, which is the point: one subgradient
# condition, measured one way, whatever named it.
.fb_sel_kkt_residual_lambda <- function(A, z, beta, lambda_l1, lambda_l2) {
  beta <- as.numeric(beta)
  if (length(beta) == 0L) {
    return(0)
  }
  gradient <- drop(crossprod(A, drop(A %*% beta) - as.numeric(z))) +
    lambda_l2 * beta

  active <- abs(beta) > 0
  residual <- numeric(length(beta))
  residual[active] <- abs(gradient[active] + lambda_l1 * sign(beta[active]))
  residual[!active] <- pmax(abs(gradient[!active]) - lambda_l1, 0)
  max(residual)
}


# --- the proximal-gradient mode ---------------------------------------- #

# .fb_sel_soft_threshold() --- the proximal operator of l1 ||.||_1.
.fb_sel_soft_threshold <- function(x, threshold) {
  sign(x) * pmax(abs(x) - threshold, 0)
}


# .fb_sel_prox_mode() --- the canonical mode by FISTA, with a KKT
# certificate and an active-set polish.
#
# Accelerated proximal gradient on the smooth part, step 1/L with L the
# largest eigenvalue of A'A plus the ridge coefficient, soft thresholding as
# the prox. Two things are added to the textbook iteration and both earn
# their place.
#
# **Adaptive restart.** FISTA's momentum is not monotone; restarting it when
# the objective rises removes the ripple that otherwise keeps the residual
# an order of magnitude above where the iterate already is.
#
# **Active-set polish.** First-order methods identify the support long before
# they resolve its values, and the last three digits are exactly what a KKT
# tolerance of 1e-10 asks for. Once the support has settled, the restricted
# problem is smooth and has a closed form:
#
#   (A_S' A_S + l2 I) beta_S = A_S' z - l1 sign(beta_S)
#
# solved by Cholesky. The polish is accepted only if the signs it assumed
# survive it and the inactive coordinates stay inside the subdifferential --
# which is not an assertion about the algorithm, it is the KKT residual
# recomputed at the polished point.
#
# **Refusal, not a returned near-miss.** If the residual is still above `tol`
# when the iteration budget runs out, this refuses with
# `selection_prox_mode_not_converged` rather than returning a point that
# looks like an optimum. A mode oracle that hands back its best effort is not
# an oracle.
.fb_sel_prox_mode <- function(
  A,
  z,
  family = c("ridge", "laplace", "elastic_net"),
  hyper,
  tol = 1e-10,
  max_iter = 1e5
) {
  family <- match.arg(family)
  penalty <- .fb_sel_canonical_penalty(family, hyper, n = nrow(as.matrix(A)))
  .fb_sel_prox_mode_lambda(
    A,
    z,
    lambda_l1 = penalty$lambda_l1,
    lambda_l2 = penalty$lambda_l2,
    tol = tol,
    max_iter = max_iter,
    family = family
  )
}


# .fb_sel_prox_mode_lambda() --- the mode in the (l1, l2) coordinate.
#
# The family wrapper above is one line of arithmetic on top of this; the path
# arm calls this directly, because a grid point on a lambda path is not the
# canonical lambda of any prior and inventing hyperparameters to reach it
# would be a fiction in the certificate.
.fb_sel_prox_mode_lambda <- function(
  A,
  z,
  lambda_l1,
  lambda_l2,
  tol = 1e-10,
  max_iter = 1e5,
  family = "canonical"
) {
  A <- as.matrix(A)
  z <- as.numeric(z)
  if (nrow(A) != length(z)) {
    stop(
      "`A` and `z` must agree on the number of rows: ",
      nrow(A),
      " and ",
      length(z),
      ".",
      call. = FALSE
    )
  }
  max_iter <- as.integer(max_iter)
  if (is.na(max_iter) || max_iter < 1L) {
    stop("`max_iter` must be at least one iteration.", call. = FALSE)
  }
  for (part in list(
    list(lambda_l1, "lambda_l1"),
    list(lambda_l2, "lambda_l2")
  )) {
    if (
      !is.numeric(part[[1L]]) ||
        length(part[[1L]]) != 1L ||
        !is.finite(part[[1L]]) ||
        part[[1L]] < 0
    ) {
      stop(
        "`",
        part[[2L]],
        "` must be a single finite non-negative number.",
        call. = FALSE
      )
    }
  }

  p <- ncol(A)
  l1 <- as.numeric(lambda_l1)
  l2 <- as.numeric(lambda_l2)

  if (p == 0L) {
    return(list(
      beta = numeric(0),
      kkt = 0,
      iterations = 0L,
      converged = TRUE,
      family = family,
      lambda_l1 = l1,
      lambda_l2 = l2,
      objective = 0.5 * sum(z^2),
      polished = FALSE
    ))
  }

  gram <- crossprod(A)
  aty <- drop(crossprod(A, z))

  # The exactly smooth corner needs no iteration at all: with no L1 term the
  # canonical objective is a positive-definite quadratic whenever l2 > 0, and
  # a rank-revealing solve covers l2 = 0 as well.
  if (l1 == 0) {
    beta <- .fb_sel_smooth_mode(gram, aty, l2)
    residual <- .fb_sel_kkt_residual_lambda(A, z, beta, l1, l2)
    if (residual > tol) {
      stop(.fb_refusal_condition(
        reason_code = "selection_prox_mode_not_converged",
        message = paste0(
          "fb_select(): the smooth-corner solve for family `",
          family,
          "` left a KKT residual of ",
          format(residual, digits = 3),
          ", above the tolerance ",
          format(tol, digits = 3),
          ". The design is too ill-conditioned for this route to certify a ",
          "mode; no point is returned, because a mode that fails its own ",
          "optimality check is not a mode."
        ),
        kkt = residual,
        tol = tol,
        iterations = 0L
      ))
    }
    return(list(
      beta = stats::setNames(beta, colnames(A)),
      kkt = residual,
      iterations = 0L,
      converged = TRUE,
      family = family,
      lambda_l1 = l1,
      lambda_l2 = l2,
      objective = .fb_sel_canonical_objective(A, z, beta, l1, l2),
      polished = FALSE
    ))
  }

  lipschitz <- max(eigen(gram, symmetric = TRUE, only.values = TRUE)$values) +
    l2
  if (!is.finite(lipschitz) || lipschitz <= 0) {
    lipschitz <- 1
  }
  step <- 1 / lipschitz

  beta <- numeric(p)
  momentum <- beta
  weight <- 1
  objective <- .fb_sel_canonical_objective(A, z, beta, l1, l2)
  best <- .fb_sel_kkt_residual_lambda(A, z, beta, l1, l2)
  polished <- FALSE
  iterations <- 0L
  check_every <- 20L

  for (iteration in seq_len(max_iter)) {
    iterations <- iteration
    gradient <- drop(gram %*% momentum) - aty + l2 * momentum
    proposal <- .fb_sel_soft_threshold(
      momentum - step * gradient,
      step * l1
    )
    proposal_objective <- .fb_sel_canonical_objective(
      A,
      z,
      proposal,
      l1,
      l2
    )
    if (proposal_objective > objective) {
      # Adaptive restart: the momentum has overshot, so drop it and take the
      # plain proximal step from the incumbent.
      weight <- 1
      gradient <- drop(gram %*% beta) - aty + l2 * beta
      proposal <- .fb_sel_soft_threshold(beta - step * gradient, step * l1)
      proposal_objective <- .fb_sel_canonical_objective(
        A,
        z,
        proposal,
        l1,
        l2
      )
    }
    next_weight <- (1 + sqrt(1 + 4 * weight^2)) / 2
    momentum <- proposal + ((weight - 1) / next_weight) * (proposal - beta)
    beta <- proposal
    weight <- next_weight
    objective <- proposal_objective

    if (iteration %% check_every == 0L || iteration == max_iter) {
      best <- .fb_sel_kkt_residual_lambda(A, z, beta, l1, l2)
      if (best <= tol) {
        break
      }
      candidate <- .fb_sel_polish(A, z, gram, aty, beta, l1, l2)
      if (!is.null(candidate)) {
        candidate_kkt <- .fb_sel_kkt_residual_lambda(
          A,
          z,
          candidate,
          l1,
          l2
        )
        if (candidate_kkt < best) {
          beta <- candidate
          momentum <- candidate
          weight <- 1
          objective <- .fb_sel_canonical_objective(A, z, beta, l1, l2)
          best <- candidate_kkt
          polished <- TRUE
        }
        if (best <= tol) {
          break
        }
      }
    }
  }

  if (best > tol) {
    stop(.fb_refusal_condition(
      reason_code = "selection_prox_mode_not_converged",
      message = paste0(
        "fb_select(): the proximal-gradient mode for family `",
        family,
        "` stopped after ",
        iterations,
        " iteration(s) with a KKT residual ",
        "of ",
        format(best, digits = 3),
        ", above the tolerance ",
        format(tol, digits = 3),
        ". No point is returned: a mode that fails its own optimality check ",
        "is not a mode, and a provider that reports convergence anyway is ",
        "the failure this certificate exists to catch."
      ),
      kkt = best,
      tol = tol,
      iterations = iterations
    ))
  }

  list(
    beta = stats::setNames(beta, colnames(A)),
    kkt = best,
    iterations = iterations,
    converged = TRUE,
    family = family,
    lambda_l1 = l1,
    lambda_l2 = l2,
    objective = .fb_sel_canonical_objective(A, z, beta, l1, l2),
    polished = polished
  )
}


# .fb_sel_smooth_mode() --- the L1-free corner: (A'A + l2 I) beta = A'z.
#
# With l2 > 0 the system is positive definite and a Cholesky solve is exact.
# With l2 = 0 -- the no-penalty limit -- it may be singular, and a
# rank-revealing QR returns the estimable-space answer with the aliased
# coordinates at zero, which is the convention `.fb_sel_gls_fit()` already
# uses in this module.
.fb_sel_smooth_mode <- function(gram, aty, lambda_l2) {
  p <- ncol(gram)
  system <- gram
  diag(system) <- diag(system) + lambda_l2
  factor_system <- tryCatch(chol(system), error = function(e) NULL)
  if (!is.null(factor_system)) {
    return(drop(backsolve(
      factor_system,
      forwardsolve(t(factor_system), aty)
    )))
  }
  decomposition <- qr(system)
  beta <- qr.coef(decomposition, aty)
  beta[is.na(beta)] <- 0
  as.numeric(beta)
}


# .fb_sel_polish() --- the exact restricted solve on the current support.
#
# Returns NULL when there is nothing to polish or when the restricted system
# has no Cholesky factor; the caller keeps iterating in that case. It never
# decides on its own that the polished point is better -- the caller
# recomputes the KKT residual and keeps whichever point wins.
.fb_sel_polish <- function(A, z, gram, aty, beta, lambda_l1, lambda_l2) {
  support <- which(abs(beta) > 0)
  if (length(support) == 0L || length(support) > nrow(A)) {
    return(NULL)
  }
  signs <- sign(beta[support])
  system <- gram[support, support, drop = FALSE]
  diag(system) <- diag(system) + lambda_l2
  factor_system <- tryCatch(chol(system), error = function(e) NULL)
  if (is.null(factor_system)) {
    return(NULL)
  }
  rhs <- aty[support] - lambda_l1 * signs
  solved <- drop(backsolve(
    factor_system,
    forwardsolve(t(factor_system), rhs)
  ))
  # A sign that flipped means the support assumption the solve rested on was
  # wrong, so the polished point is not a candidate at all.
  if (any(sign(solved) != signs)) {
    return(NULL)
  }
  out <- numeric(length(beta))
  out[support] <- solved
  out
}


# --- the equivalence level --------------------------------------------- #

# .fb_sel_equivalence_providers() --- the optimiser providers a comparison
# may name, and the families each can express.
#
# `none` is the flat, no-penalty corner: its "mode identity" is the
# generalised least-squares solution, which is the Tier-1 limit section 8.2
# gates alongside ridge and lasso.
.fb_sel_equivalence_providers <- function() {
  list(
    glmnet = c("none", "ridge", "laplace", "elastic_net"),
    closed_form = c("none", "ridge"),
    # The internal proximal solver expresses every family this release
    # ships, grouped ones included: each has a proximal operator written
    # from its own subdifferential and a KKT certificate measured at the
    # point it returns.
    prox = c(
      "none", "ridge", "laplace", "elastic_net",
      .fb_sel_grouped_families()
    ),
    gls_reml_plugin = "none"
  )
}


# .fb_sel_equivalence_level() --- which of section 8.1's three levels a
# comparison carries.
#
# The three levels are ordered by what they permit a reader to conclude, and
# a comparison carries exactly one:
#
#   map_exact_fixed_covariance      the two arms solve the same problem, and
#                                   a numerical oracle gates the claim
#   posterior_predictive_comparable the provider can express the family, but
#                                   the covariance is integrated or sigma is
#                                   free, so the estimands differ and only
#                                   predictive risk on common folds compares
#   no_identity_expected            the provider has no representation of the
#                                   family at all, OR the two arms were never
#                                   scored on one loss at one unit on one set
#                                   of folds; the arms are reported side by
#                                   side and nothing is claimed
#
# The rule that decides is deliberately narrow. `map_exact_fixed_covariance`
# is the only level with a numerical gate behind it, so it is the only one
# that may be reached by default, and it is reached only when every condition
# the identity rests on is stated as true by the caller.
#
# `shared_predictive_table` is the middle level's own gate and it defaults to
# FALSE. `posterior_predictive_comparable` says the arms were scored on the
# same folds, target, loss, unit and candidate design -- and the Bayesian
# penalised route did not write an out-of-fold score at all, so that sentence
# was printed about a comparison nobody had made. It is now printed only when
# the caller states that both arms carry the same unit-level table, and the
# caller states it by having built one.
.fb_sel_equivalence_level <- function(
  prior_family,
  provider,
  covariance_fixed,
  sigma_fixed,
  shared_predictive_table = FALSE
) {
  providers <- .fb_sel_equivalence_providers()
  if (
    !is.character(prior_family) ||
      length(prior_family) != 1L ||
      !nzchar(prior_family)
  ) {
    stop("`prior_family` must be a single non-empty string.", call. = FALSE)
  }
  provider <- match.arg(provider, names(providers))
  flag <- function(value, what) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("`", what, "` must be TRUE or FALSE.", call. = FALSE)
    }
    value
  }
  covariance_fixed <- flag(covariance_fixed, "covariance_fixed")
  sigma_fixed <- flag(sigma_fixed, "sigma_fixed")
  shared_predictive_table <- flag(
    shared_predictive_table, "shared_predictive_table"
  )

  expressible <- providers[[provider]]
  level <- NULL
  reason <- NULL

  if (!prior_family %in% expressible) {
    level <- "no_identity_expected"
    reason <- paste0(
      "`",
      provider,
      "` has no representation of a `",
      prior_family,
      "` prior, so the two arms estimate different objects and are ",
      "reported side by side. It expresses: ",
      paste(expressible, collapse = ", "),
      "."
    )
  } else if (covariance_fixed && sigma_fixed) {
    level <- "map_exact_fixed_covariance"
    reason <- paste0(
      "the covariance and the residual scale are both held fixed, so the ",
      "posterior mode under a `",
      prior_family,
      "` prior is the minimiser of the canonical penalised objective that ",
      "`",
      provider,
      "` solves, up to the stated parameter conversion."
    )
  } else {
    missing_condition <- c(
      if (!covariance_fixed) {
        "the covariance is integrated over rather than plugged in"
      },
      if (!sigma_fixed) "the residual scale is free"
    )
    if (shared_predictive_table) {
      level <- "posterior_predictive_comparable"
      reason <- paste0(
        "`",
        provider,
        "` expresses a `",
        prior_family,
        "` penalty, but ",
        paste(missing_condition, collapse = " and "),
        ", so the two arms estimate different objects. Scored on the same ",
        "folds, target, loss, unit and candidate design, their predictive ",
        "risk is comparable; their coefficients and supports are not."
      )
    } else {
      level <- "no_identity_expected"
      reason <- paste0(
        "`",
        provider,
        "` expresses a `",
        prior_family,
        "` penalty, but ",
        paste(missing_condition, collapse = " and "),
        ", so no mode identity holds; and the two arms carry no common ",
        "out-of-fold score -- one unit-level loss, on one set of folds, ",
        "under one prediction target -- so their predictive risk was never ",
        "compared either. The arms are reported side by side and nothing ",
        "is claimed. `benchmark = TRUE` scores the Bayesian arm on the ",
        "benchmark's own folds and is what makes the middle level ",
        "available."
      )
    }
  }

  list(
    level = level,
    reason = reason,
    prior_family = prior_family,
    provider = provider,
    covariance_fixed = covariance_fixed,
    sigma_fixed = sigma_fixed,
    shared_predictive_table = shared_predictive_table
  )
}


# --- the provider lambda map ------------------------------------------- #

# .fb_sel_response_scale() --- the response scale a provider's penalty is
# measured against, under a named rule.
#
# The rules are the candidates the convention measurement discriminates
# between; the module never picks one from documentation, it picks the one
# that reproduces a closed form.
.fb_sel_response_scale <- function(
  z,
  rule = c("rms", "sd_population", "sd_sample", "none")
) {
  rule <- match.arg(rule)
  z <- as.numeric(z)
  switch(
    rule,
    rms = sqrt(mean(z^2)),
    sd_population = sqrt(mean((z - mean(z))^2)),
    sd_sample = stats::sd(z),
    none = 1
  )
}


# .fb_sel_convention_response_scale() --- the numeric scale a convention
# needs, refused rather than defaulted when it is absent.
.fb_sel_convention_response_scale <- function(convention) {
  if (!is.list(convention)) {
    stop("`convention` must be a measured convention list.", call. = FALSE)
  }
  scale <- convention[["response_scale"]]
  if (
    !is.numeric(scale) ||
      length(scale) != 1L ||
      !is.finite(scale) ||
      scale <= 0
  ) {
    stop(
      "internal contract violation: the lambda map needs the convention's ",
      "`response_scale` for this problem. Fill it with ",
      "`.fb_sel_convention_at(convention, z)` -- defaulting it to one is ",
      "the naive map, and the naive map is measurably wrong for the L2 part.",
      call. = FALSE
    )
  }
  as.numeric(scale)
}


# .fb_sel_canonical_to_glmnet() --- the canonical coordinate as `glmnet`'s.
#
# The canonical penalty is written in the composite coordinate
#
#   pen(beta) = lc [ alpha ||beta||_1 + (1 - alpha)/2 ||beta||_2^2 ]
#
# and `glmnet`, under `standardize = FALSE, intercept = FALSE`, is MEASURED
# to minimise
#
#   1/(2n) ||z - A beta||^2 + lg [ ag ||beta||_1 +
#                                  (1 - ag)/(2 ys) ||beta||_2^2 ]
#
# with `ys` the response scale of the measured convention. Matching the two
# term by term after multiplying through by n gives
#
#   lg = lc (alpha + (1 - alpha) ys) / n
#   ag = alpha / (alpha + (1 - alpha) ys)
#
# The L1 part carries no `ys` and the L2 part carries it once, which is why
# the naive `lambda / n` map is exact for a lasso and off by the response
# scale for a ridge.
.fb_sel_canonical_to_glmnet <- function(
  lambda_canonical,
  alpha,
  n,
  convention
) {
  if (
    !is.numeric(lambda_canonical) ||
      length(lambda_canonical) != 1L ||
      !is.finite(lambda_canonical) ||
      lambda_canonical < 0
  ) {
    stop(
      "`lambda_canonical` must be a single finite non-negative number.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(alpha) ||
      length(alpha) != 1L ||
      !is.finite(alpha) ||
      alpha < 0 ||
      alpha > 1
  ) {
    stop("`alpha` must be a single number in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || n <= 0) {
    stop("`n` must be a single positive number of rows.", call. = FALSE)
  }
  ys <- .fb_sel_convention_response_scale(convention)

  mixture <- alpha + (1 - alpha) * ys
  list(
    lambda = lambda_canonical * mixture / n,
    alpha = alpha / mixture,
    response_scale = ys,
    n = as.numeric(n)
  )
}


# .fb_sel_glmnet_to_canonical() --- the inverse map, exact by construction.
.fb_sel_glmnet_to_canonical <- function(
  lambda_glmnet,
  alpha_glmnet,
  n,
  convention
) {
  if (
    !is.numeric(lambda_glmnet) ||
      length(lambda_glmnet) != 1L ||
      !is.finite(lambda_glmnet) ||
      lambda_glmnet < 0
  ) {
    stop(
      "`lambda_glmnet` must be a single finite non-negative number.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(alpha_glmnet) ||
      length(alpha_glmnet) != 1L ||
      !is.finite(alpha_glmnet) ||
      alpha_glmnet < 0 ||
      alpha_glmnet > 1
  ) {
    stop("`alpha_glmnet` must be a single number in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || n <= 0) {
    stop("`n` must be a single positive number of rows.", call. = FALSE)
  }
  ys <- .fb_sel_convention_response_scale(convention)

  mixture <- alpha_glmnet + (1 - alpha_glmnet) / ys
  list(
    lambda_canonical = n * lambda_glmnet * mixture,
    alpha = alpha_glmnet / mixture,
    response_scale = ys,
    n = as.numeric(n)
  )
}
