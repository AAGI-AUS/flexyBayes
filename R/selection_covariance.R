# Selection covariance --- the whitening substrate for fb_select()
#
# Contract: SPEC_penalised_selection_v3_2026-09-11.md sections 4.3 and 5.1.
# Everything downstream of selection is computed on data whitened by the
# operator built here, so this file carries oracles 1 and 2 of the
# correctness ledger and nothing may be trusted before those are red on
# their own defect first.
#
# Two builders, chosen by structure and never by argument (D3):
#
#   block_cholesky          every random term nests within one observed
#                           grouping factor, so V is block-diagonal and
#                           the Cholesky is taken block by block.
#   random_term_reduction   the random terms are crossed, V is not
#                           block-diagonal, and the candidate block is
#                           constant within the levels of some term. The
#                           problem reduces to weighted regression on the
#                           posterior effects. Built in a later tranche.
#
# The dense reference below is the oracle, not a production path. It
# materialises n x n and refuses above a cap, because the deposit's V at
# 73,341 observations is about 43 GB.


# --- reading the base fit -------------------------------------------- #

# .fb_sel_summary() --- the one place the module reads the fit's summary.
#
# It used to capture output, because summary.flexybayes_inla() printed 67
# lines unconditionally and the substrate reads it once per fold. That is
# fixed at the source now, so this is a plain call. It stays as a named
# function because every fold goes through it and a future change to how
# the base fit is read belongs in one place.
.fb_sel_summary <- function(object) {
  summary(object)
}


# --- variance components ------------------------------------------- #

# .fb_sel_variance_components() --- theta-hat under a named estimand.
#
# `summary(object)$varcomp` reports the marginal posterior of each
# standard deviation: `estimate` is its mean, `std.error` its standard
# deviation. V is linear in the variances, so
#
#   E[V] = sum_k E[sigma_k^2] Z_k Z_k' + E[sigma_e^2] I
#
# holds by linearity alone, with no joint posterior and no independence
# assumption, and E[sigma^2] = mean^2 + sd^2 exactly. That is the whole
# of `plugin_mean_covariance` and it is exact rather than approximate.
#
# Squaring the posterior mean instead gives (E sigma)^2, which is the
# "vector of marginal component medians" the contract forbids. The two
# differ by the posterior variance of the component, so they diverge
# precisely where a component is weakly identified --- which is where a
# wrong covariance does the most damage.
#
# Returns a list with `residual`, a named `terms` vector, both on the
# VARIANCE scale, and the estimand label for the certificate.
.fb_sel_variance_components <- function(
  object,
  estimand = c(
    "plugin_mean_covariance",
    "reml_plugin_covariance",
    "point_estimate_covariance"
  )
) {
  estimand <- match.arg(estimand)
  # Every builder that reads these components writes the residual as
  # sigma^2 I, with no place for a known variance per row
  # (R/known_variance.R).
  ir <- object$fb %||% object$extras$fb_terms
  if (!is.null(.fb_known_variance_term(ir))) {
    .fb_sel_known_variance_refuse_route(paste0(
      "a covariance built from the fitted base's own variance components, ",
      "whose residual is one variance common to every row,"
    ))
  }
  # The fit's variances and the blocks they belong to, through the core's
  # contract (R/marginal_cov.R).
  marginal <- .fb_sel_marginal_cov(object)
  vc <- marginal$varcomp

  structured <- Filter(function(b) !is.null(b$structure), marginal$blocks)
  unhandled <- c(
    if (!is.null(marginal$residual$by)) {
      paste0("sigma_", names(marginal$residual$variance))
    },
    vapply(structured, function(b) b$label, character(1L))
  )
  if (length(unhandled) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_covariance_not_buildable",
      message = paste0(
        "fb_select(): the base fit carries variance components this builder ",
        "does not know how to place in V: ",
        paste(unhandled, collapse = ", "),
        ". Correlation, factor-analytic and known-covariance structures and ",
        "a sectioned residual are not admitted in this release."
      )
    ))
  }

  # Squaring is the only place the estimand enters. Keep it in one
  # expression so the two summaries cannot drift apart.
  to_variance <- switch(estimand,
    plugin_mean_covariance = function(m, s) m^2 + s^2,
    reml_plugin_covariance = function(m, s) m^2,
    point_estimate_covariance = function(m, s) m^2
  )

  # In the order of the variance-component table.
  rows <- vapply(
    marginal$blocks,
    function(b) which(vc$component %in% b$components)[[1L]],
    integer(1L)
  )
  rows <- rows[order(rows)]
  terms <- to_variance(vc$estimate[rows], vc$std.error[rows])
  names(terms) <- names(rows)
  residual <- vc[vc$component == "sigma", , drop = FALSE]

  list(
    residual = to_variance(residual$estimate[[1L]], residual$std.error[[1L]]),
    terms = terms,
    estimand = estimand
  )
}


# --- random-effect designs ----------------------------------------- #

# .fb_sel_random_designs() --- one indicator matrix per random term,
# keyed by the name the variance-component table uses, so Z and sigma^2
# can be matched by name rather than by position. A positional match
# here would be a silent misassignment of variance to structure.
.fb_sel_random_designs <- function(object) {
  rt <- object$fb$random_terms
  d <- object$data

  if (length(rt) == 0L) {
    return(list())
  }

  out <- vector("list", length(rt))
  nms <- character(length(rt))

  for (i in seq_along(rt)) {
    term <- rt[[i]]
    if (identical(term$type, "simple")) {
      nms[[i]] <- term$var
      f <- factor(d[[term$var]])
    } else if (identical(term$type, "nested")) {
      nms[[i]] <- paste0(term$outer, ":", term$inner)
      f <- factor(paste0(
        as.character(d[[term$outer]]), ":", as.character(d[[term$inner]])
      ))
    } else {
      stop(.fb_refusal_condition(
        reason_code = "selection_covariance_not_buildable",
        message = paste0(
          "fb_select(): random term of type '", term$type,
          "' is not admitted by the covariance builder in this release. ",
          "Admitted types: simple, nested."
        )
      ))
    }
    # Indicator, not a contrast: Z Z' must reproduce the block structure.
    out[[i]] <- stats::model.matrix(~ 0 + f)
    colnames(out[[i]]) <- levels(f)
  }

  names(out) <- nms
  out
}


# --- the dense reference (oracle 1) -------------------------------- #

# .fb_sel_V_dense() --- V assembled the obvious way, elementwise, with
# no structural cleverness. This exists to be right rather than fast:
# it is the reference the production builder is checked against, and it
# refuses rather than exhaust memory.
.fb_sel_V_dense <- function(object, estimand = "plugin_mean_covariance",
                            max_n = 4000L) {
  n <- nrow(object$data)
  if (n > max_n) {
    stop(.fb_refusal_condition(
      reason_code = "selection_dense_reference_too_large",
      message = paste0(
        "fb_select(): the dense covariance reference is ", n, " by ", n,
        " (about ", format(round(n^2 * 8 / 2^30, 1), nsmall = 1),
        " GB), above the cap of ", max_n, ". The reference is an oracle for ",
        "small cases, not a production path."
      )
    ))
  }

  vc <- .fb_sel_variance_components(object, estimand)
  Z <- .fb_sel_random_designs(object)

  missing_var <- setdiff(names(Z), names(vc$terms))
  if (length(missing_var) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_covariance_not_buildable",
      message = paste0(
        "fb_select(): no variance component was reported for random term(s) ",
        paste(missing_var, collapse = ", "),
        ". Z and theta must match by name."
      )
    ))
  }

  V <- diag(vc$residual, nrow = n, ncol = n)
  for (nm in names(Z)) {
    V <- V + vc$terms[[nm]] * tcrossprod(Z[[nm]])
  }
  dimnames(V) <- NULL
  V
}


# --- block structure ------------------------------------------------ #

# .fb_sel_block_structure() --- does builder A apply?
#
# It applies when some observed factor g partitions the rows such that
# every random term's levels sit wholly inside one level of g. Then V is
# block-diagonal over g and the Cholesky is taken block by block.
#
# The search is over the random terms' own grouping variables. A crossed
# design (site by year with variety across both) has no such g, and the
# function says so rather than returning a coarse partition that happens
# to be one block.
.fb_sel_block_structure <- function(object) {
  rt <- object$fb$random_terms
  d <- object$data
  n <- nrow(d)

  residual <- object$fb$residual_terms
  if (length(residual) != 1L || !identical(residual[[1L]]$type, "units")) {
    return(list(
      applicable = FALSE,
      reason = "the residual term is not independent units"
    ))
  }

  if (length(rt) == 0L) {
    return(list(
      applicable = TRUE, factor = NULL,
      blocks = list(seq_len(n)), n_blocks = 1L,
      reason = "no random terms: V is diagonal"
    ))
  }

  # Candidate grouping factors: the outermost variable of each term.
  candidates <- unique(vapply(
    rt,
    function(t) if (identical(t$type, "simple")) t$var else t$outer,
    character(1L)
  ))

  # The row-level factor each random term induces.
  term_factors <- lapply(rt, function(t) {
    if (identical(t$type, "simple")) {
      as.character(d[[t$var]])
    } else {
      paste0(as.character(d[[t$outer]]), ":", as.character(d[[t$inner]]))
    }
  })

  for (g in candidates) {
    gv <- as.character(d[[g]])
    # Every term level must map to exactly one level of g, or the term
    # links rows across blocks and V is not block-diagonal.
    ok <- TRUE
    for (tf in term_factors) {
      spans <- tapply(gv, tf, function(x) length(unique(x)))
      if (any(spans > 1L)) {
        ok <- FALSE
        break
      }
    }
    if (ok) {
      idx <- split(seq_len(n), factor(gv))
      return(list(
        applicable = TRUE,
        factor = g,
        blocks = unname(idx),
        block_labels = names(idx),
        n_blocks = length(idx),
        reason = paste0("every random term nests within `", g, "`")
      ))
    }
  }

  list(
    applicable = FALSE,
    reason = paste0(
      "no grouping factor contains every random term: the design is ",
      "crossed, so V is not block-diagonal"
    )
  )
}


# --- the whitener (oracle 2) ---------------------------------------- #

# .fb_sel_whitener() --- W with W V W' = I, held as per-block Cholesky
# factors rather than as a matrix. V^-1 is never formed, and neither is
# W: applying it is a forward solve.
.fb_sel_whitener <- function(object, estimand = "plugin_mean_covariance") {
  bs <- .fb_sel_block_structure(object)
  if (!isTRUE(bs$applicable)) {
    stop(.fb_refusal_condition(
      reason_code = "selection_covariance_not_buildable",
      message = paste0(
        "fb_select(): the block-Cholesky builder does not apply because ",
        bs$reason, ". The random-term reduction is the route for this base, ",
        "and it is not built in this tranche."
      )
    ))
  }

  vc <- .fb_sel_variance_components(object, estimand)
  Z <- .fb_sel_random_designs(object)

  factors <- vector("list", length(bs$blocks))
  for (b in seq_along(bs$blocks)) {
    rows <- bs$blocks[[b]]
    nb <- length(rows)
    Vb <- diag(vc$residual, nrow = nb, ncol = nb)
    for (nm in names(Z)) {
      Zb <- Z[[nm]][rows, , drop = FALSE]
      # Columns empty within this block contribute nothing; dropping them
      # keeps the crossproduct at block scale.
      keep <- colSums(Zb) > 0
      if (any(keep)) {
        Vb <- Vb + vc$terms[[nm]] * tcrossprod(Zb[, keep, drop = FALSE])
      }
    }
    ch <- tryCatch(chol(Vb), error = function(e) NULL)
    if (is.null(ch)) {
      stop(.fb_refusal_condition(
        reason_code = "selection_covariance_not_positive_definite",
        message = paste0(
          "fb_select(): block ", b, " of the covariance is not positive ",
          "definite, so it has no Cholesky factor. No jitter was applied: a ",
          "repair here would be silent misspecification."
        )
      ))
    }
    # V = L L' with L lower triangular; W = L^-1 is a forward solve.
    factors[[b]] <- list(rows = rows, L = t(ch))
  }

  structure(
    list(
      builder = "block_cholesky",
      estimand = vc$estimand,
      n = nrow(object$data),
      block_factor = bs$factor,
      n_blocks = bs$n_blocks,
      factors = factors,
      components = vc
    ),
    class = "fb_sel_whitener"
  )
}


# .fb_sel_whiten() --- apply W to a vector or matrix, block by block.
.fb_sel_whiten <- function(w, M) {
  if (!inherits(w, "fb_sel_whitener")) {
    stop("`w` must be an <fb_sel_whitener>.", call. = FALSE)
  }
  M <- as.matrix(M)
  if (nrow(M) != w$n) {
    stop(.fb_refusal_condition(
      reason_code = "selection_covariance_not_buildable",
      message = paste0(
        "fb_select(): the whitener was built for ", w$n, " rows and was ",
        "applied to ", nrow(M), "."
      )
    ))
  }
  out <- matrix(NA_real_, nrow = nrow(M), ncol = ncol(M))
  colnames(out) <- colnames(M)
  for (f in w$factors) {
    out[f$rows, ] <- forwardsolve(f$L, M[f$rows, , drop = FALSE])
  }
  out
}


# .fb_sel_whitener_as_dense() --- W materialised, for the oracles only.
.fb_sel_whitener_as_dense <- function(w) {
  .fb_sel_whiten(w, diag(1, nrow = w$n, ncol = w$n))
}


#' @export
print.fb_sel_whitener <- function(x, ...) {
  cat("<fb_sel_whitener> builder: ", x$builder, "\n", sep = "")
  cat("  estimand : ", x$estimand, "\n", sep = "")
  cat("  rows     : ", x$n, "\n", sep = "")
  cat(
    "  blocks   : ", x$n_blocks,
    if (!is.null(x$block_factor)) paste0(" over `", x$block_factor, "`") else
      "",
    "\n",
    sep = ""
  )
  invisible(x)
}
