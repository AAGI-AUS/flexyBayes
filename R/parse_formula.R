# Formula parsing for flexyBayes
# Translates ASReml-style formulas into internal term descriptors.
# Not exported.

# Parse a standard R fixed formula into response + intercept + terms list
#
# @param fixed Two-sided formula: response ~ predictors
# @param data data.frame
# @return List with response, intercept (logical), terms (list of
#   term descriptors)
.parse_fixed <- function(fixed, data) {
  if (length(fixed) < 3) {
    stop(.fb_refusal_condition(
      reason_code = "formula_not_two_sided",
      message = "The model formula must be two-sided: response ~ predictors."
    ))
  }
  response <- deparse(fixed[[2]])
  if (!response %in% names(data)) {
    stop(.fb_refusal_condition(
      reason_code = "response_not_in_data",
      message = paste0(
        "Response variable \"",
        response,
        "\" not found in data."
      )
    ))
  }

  tt <- tryCatch(terms(fixed, data = data), error = function(e) terms(fixed))
  intercept <- attr(tt, "intercept") == 1
  term_labels <- attr(tt, "term.labels")

  term_info <- lapply(term_labels, function(lbl) {
    .classify_fixed_term(lbl, data)
  })

  list(response = response, intercept = intercept, terms = term_info)
}

# Classify a single fixed-formula term label
.classify_fixed_term <- function(lbl, data) {
  if (grepl("^I\\(", lbl)) {
    return(list(type = "expression", label = lbl))
  }
  if (grepl(":", lbl, fixed = TRUE)) {
    parts <- strsplit(lbl, ":")[[1]]
    # Per-part classification: factor / continuous / unknown.
    kinds <- vapply(
      parts,
      function(p) {
        if (!p %in% names(data)) {
          return("unknown")
        }
        col <- data[[p]]
        if (is.factor(col) || is.character(col)) "factor" else "continuous"
      },
      character(1)
    )

    # All factors -> factor:factor crossed effect.
    if (all(kinds == "factor")) {
      return(list(type = "factor_interaction", vars = parts, label = lbl))
    }

    # Mixed factor x continuous (v0.2.6). Treatment-
    # coded indexed slopes: the reference level's slope is structurally
    # fixed at zero; the L-1 non-reference levels carry per-level
    # deviation parameters. v0.2 collapsed this to
    # `beta * as_data(factor) * as_data(continuous)`, which silently
    # coerced the factor to its integer codes and produced a posterior
    # disagreeing with lme4::lmer by > 3 sigma; the v0.2.6 path closes
    # that correctness bug. Restricted to two-way mixed terms in v0.2.6
    # (one factor + one continuous); higher-order combinations stay on
    # the legacy `interaction` branch pending v0.3.
    if (
      length(parts) == 2L &&
        sum(kinds == "factor") == 1L &&
        sum(kinds == "continuous") == 1L
    ) {
      fac_idx <- which(kinds == "factor")[[1L]]
      con_idx <- which(kinds == "continuous")[[1L]]
      fac_name <- parts[[fac_idx]]
      con_name <- parts[[con_idx]]
      # Unsupported-contrast detection. Only the treatment-coded path
      # is shipped at v0.2.6; helmert / sdif / sum-to-zero / ordered /
      # user-supplied custom matrices defer to v0.3 representation IR.
      .stop_if_unsupported_contrast(fac_name, data[[fac_name]])
      f <- factor(data[[fac_name]])
      return(list(
        type = "factor_numeric_interaction",
        factor = fac_name,
        continuous = con_name,
        vars = parts,
        levels = levels(f),
        n_levels = nlevels(f),
        label = lbl
      ))
    }

    return(list(type = "interaction", vars = parts, label = lbl))
  }
  if (lbl %in% names(data)) {
    col <- data[[lbl]]
    if (is.factor(col) || is.character(col)) {
      f <- factor(col)
      return(list(
        type = "factor",
        var = lbl,
        n_levels = nlevels(f),
        levels = levels(f),
        label = lbl
      ))
    } else {
      return(list(type = "continuous", var = lbl, label = lbl))
    }
  }
  list(type = "expression", label = lbl)
}

# Parse an ASReml one-sided formula (random or rcov) via AST walking
#
# @param formula One-sided formula: ~ random_terms
# @param data data.frame
# @return List of enriched term descriptors
.parse_formula <- function(formula, data) {
  rhs <- if (length(formula) == 3) formula[[3]] else formula[[2]]
  terms <- .walk(rhs)
  lapply(terms, .enrich, data = data)
}

# Recursive AST walker — converts formula expression tree to term descriptors
.walk <- function(expr) {
  if (is.name(expr)) {
    nm <- as.character(expr)
    if (nm == "units") {
      return(list(list(type = "units")))
    }
    return(list(list(type = "simple", var = nm)))
  }
  if (!is.call(expr)) {
    return(list())
  }

  fn <- as.character(expr[[1]])

  # Binary operators
  if (fn == "+") {
    return(c(.walk(expr[[2]]), .walk(expr[[3]])))
  }

  if (fn == ":") {
    left <- .walk(expr[[2]])
    right <- .walk(expr[[3]])
    l <- if (length(left)) left[[1]] else list(type = "unknown")
    r <- if (length(right)) right[[1]] else list(type = "unknown")
    return(.classify_interaction(l, r))
  }

  # ASReml-specific functions
  switch(
    fn,
    # vm() ingest extends to capture optional
    # named arguments (chol, precision, blocks, low_rank_factor,
    # low_rank_scheme) alongside the historical positional / V=
    # dense-matrix arg. The cov_representation slot is the IR
    # anchor; v0.3.7 emits dense + chol + precision, v0.3.8 emits
    # blocks, low-rank reserves the slot for a future release.
    "vm" = {
      var_name <- .dep(expr, 2)
      # v0.4.0: the fb_cov() constructor front door takes
      # precedence. `cov = fb_cov(...)` is the canonical form; the four
      # legacy keyword carriers (chol / precision / blocks /
      # low_rank_factor) continue to work with a deprecation warning.
      if (!is.null(expr[["cov"]])) {
        cov_rep <- .fb_cov_call_to_representation(expr[["cov"]], "vm", var_name)
        mat_expr <- if (identical(cov_rep$format, "dense")) {
          cov_rep$data
        } else {
          NA_character_
        }
      } else {
        .warn_legacy_vm_kwargs(expr, "vm", var_name)
        mat_expr <- .extract_vm_dense_arg(expr)
        cov_rep <- .cov_representation_from_call(
          fn = "vm",
          mat = mat_expr,
          chol = .dep_named(expr, "chol"),
          precision = .dep_named(expr, "precision"),
          blocks = .dep_named(expr, "blocks"),
          low_rank_factor = .dep_named(expr, "low_rank_factor"),
          low_rank_scheme = .dep_named(expr, "low_rank_scheme")
        )
      }
      list(list(
        type = "vm",
        var = var_name,
        mat = mat_expr,
        cov_representation = cov_rep
      ))
    },
    "ide" = ,
    "id" = list(list(type = "ide", var = .dep(expr, 2))),
    # ped() gains optional use_sparse_precision
    # flag that routes through the sparse-precision format internally.
    "ped" = {
      var_name <- .dep(expr, 2)
      # v0.4.0: ped() also accepts the fb_cov() front
      # door, e.g. ped(animal, cov = fb_cov(A, type = "precision",
      # sparse_precision = TRUE)). The legacy positional A +
      # use_sparse_precision form continues to work unchanged (it is
      # not part of the deprecated keyword set).
      if (!is.null(expr[["cov"]])) {
        cov_rep <- .fb_cov_call_to_representation(
          expr[["cov"]],
          "ped",
          var_name
        )
        mat_expr <- if (identical(cov_rep$format, "dense")) {
          cov_rep$data
        } else {
          NA_character_
        }
      } else {
        mat_expr <- .extract_vm_dense_arg(expr)
        use_sp <- .dep_named(expr, "use_sparse_precision")
        cov_rep <- .cov_representation_from_call(
          fn = "ped",
          mat = mat_expr,
          use_sparse_precision = use_sp
        )
      }
      list(list(
        type = "ped",
        var = var_name,
        mat = mat_expr,
        cov_representation = cov_rep
      ))
    },
    "fa" = list(list(
      type = "fa",
      var = .dep(expr, 2),
      k = as.integer(expr[[3]])
    )),
    "us" = list(list(type = "us", var = .dep(expr, 2))),
    "at" = list(list(
      type = "at",
      var = .dep(expr, 2),
      level = if (length(expr) > 2) .dep(expr, 3) else NULL
    )),
    "ar1" = list(list(type = "ar1", var = .dep(expr, 2))),
    "ar2" = list(list(type = "ar2", var = .dep(expr, 2))),
    "cor" = list(list(type = "cor", var = .dep(expr, 2))),
    "str" = list(list(
      type = "str",
      formula = expr[[2]],
      structure = if (length(expr) > 2) expr[[3]] else NULL
    )),
    "spl" = list(list(type = "spline", var = .dep(expr, 2))),
    # mgcv-style univariate smooth s(x[, k = ..., bs = ...]).
    # v0.1 supports s() only; te() / ti() / t2()
    # tensor-product smooths defer to v0.2. The full original call is
    # captured here; .enrich() evaluates it via mgcv::smoothCon() to
    # build the n x k design matrix.
    # v0.4.0: an optional `representation`
    # argument inside s() requests an approximation scheme for the
    # smooth basis, e.g. s(x, representation = list(scheme =
    # "low_rank_smooth", rank = 5L)) (the list is the de-classed shape
    # fb_approx() returns). It is NOT an mgcv::s() argument,
    # so it is intercepted here and stripped from the call before
    # .enrich() evaluates the cleaned call via mgcv::smoothCon(); the
    # raw spec expression is carried on `approx_expr` and resolved at
    # enrich time (where the basis dimension k and row count n are known
    # for the rank refusal).
    "s" = {
      rep_expr <- if (!is.null(expr[["representation"]])) {
        expr[["representation"]]
      } else {
        NULL
      }
      sm_call <- expr
      if (!is.null(rep_expr)) {
        sm_call[["representation"]] <- NULL
      }
      list(list(
        type = "smooth_mgcv",
        call = sm_call,
        var = .dep(sm_call, 2),
        approx_expr = rep_expr
      ))
    },
    # Tensor-product / multivariate smooths te() / ti() / t2() are not
    # supported (only univariate s() / spl() ship). Refuse explicitly
    # rather than letting them reach the simple-variable default below,
    # which would silently parse e.g. te(x, z) as a term named
    # "te(x, z)" -- a no-silent-fitting charter violation.
    "te" = ,
    "ti" = ,
    "t2" = stop(.fb_refusal_condition(
      reason_code = "tensor_smooth_unsupported",
      message = paste0(
        "Tensor-product / multivariate smooth ",
        deparse(expr),
        " is not supported. flexyBayes fits univariate penalised ",
        "splines only -- use s(x) or spl(x) per smooth dimension. ",
        "Multivariate smooths are planned for a future release."
      )
    )),
    "dsum" = {
      arg <- expr[[2]]
      grp <- tryCatch(
        {
          body_expr <- if (is.call(arg) && as.character(arg[[1]]) == "~") {
            arg[[length(arg)]]
          } else {
            arg
          }
          if (is.call(body_expr) && as.character(body_expr[[1]]) == "|") {
            deparse(body_expr[[3]])
          } else {
            NULL
          }
        },
        error = function(e) NULL
      )
      if (!is.null(grp)) {
        list(list(type = "at_units", var = grp, level = NULL))
      } else {
        list(list(type = "units"))
      }
    },
    "lin" = list(list(type = "continuous", var = .dep(expr, 2))),
    "pol" = list(list(
      type = "polynomial",
      var = .dep(expr, 2),
      degree = if (length(expr) > 2) as.integer(expr[[3]]) else 2L
    )),
    # Default: treat as a simple named variable
    list(list(type = "simple", var = deparse(expr)))
  )
}

# Classify what an a:b interaction means in ASReml context
.classify_interaction <- function(l, r) {
  lt <- l$type
  rt <- r$type
  # FA(k) GxE:  fa(env,k) : id(geno)
  if (lt == "fa" && rt %in% c("id", "ide", "simple")) {
    if (!is.null(l$k) && l$k < 1L) {
      stop(.fb_refusal_condition(
        reason_code = "fa_rank_invalid",
        message = paste0("fa() requires k >= 1; got k = ", l$k, ".")
      ))
    }
    return(list(list(type = "fa_gxe", outer = l$var, k = l$k, inner = r$var)))
  }
  # US GxE:     us(env) : id(geno)
  if (lt == "us" && rt %in% c("id", "ide", "simple")) {
    return(list(list(type = "us_gxe", outer = l$var, inner = r$var)))
  }
  # at(env):units  - heterogeneous residual
  if (lt == "at" && rt == "units") {
    return(list(list(type = "at_units", var = l$var, level = l$level)))
  }
  # at(env):geno   - DIAG(env) x I(geno)
  if (lt == "at" && rt %in% c("simple", "ide", "id")) {
    return(list(list(
      type = "at_simple",
      outer = l$var,
      level = l$level,
      inner = r$var
    )))
  }
  # AR1 spatial
  if (lt == "ar1" && rt %in% c("ar1", "id", "ide", "simple", "units")) {
    return(list(list(
      type = "ar1_spatial",
      row_var = l$var,
      col_var = r$var,
      col_ar1 = (rt == "ar1")
    )))
  }
  # vm(geno,G) : id(env)  - structured GxE with known G
  if (lt == "vm" && rt %in% c("id", "ide", "simple")) {
    return(list(list(
      type = "vm_gxe",
      inner = l$var,
      mat = l$mat,
      outer = r$var
    )))
  }
  # simple:simple  - nested random effect
  if (lt %in% c("simple", "ide", "id") && rt %in% c("simple", "ide", "id")) {
    return(list(list(type = "nested", outer = l$var, inner = r$var)))
  }
  # nested:simple  - three-way combination
  if (lt == "nested" && rt %in% c("simple", "ide", "id")) {
    return(list(list(type = "combo", vars = c(l$outer, l$inner, r$var))))
  }
  # combo:simple   - four-way+ combination
  if (lt == "combo" && rt %in% c("simple", "ide", "id")) {
    return(list(list(type = "combo", vars = c(l$vars, r$var))))
  }
  # fallback
  list(list(type = "interaction_generic", left = l, right = r))
}

# Enrich a term descriptor with data-derived information
.enrich <- function(term, data) {
  .add_var_info <- function(term, var_field = "var") {
    vname <- term[[var_field]]
    if (!is.null(vname) && !is.na(vname) && vname %in% names(data)) {
      col <- data[[vname]]
      f <- factor(col)
      term[[paste0(var_field, "_n")]] <- nlevels(f)
      term[[paste0(var_field, "_levels")]] <- levels(f)
    }
    term
  }

  switch(
    term$type,
    "simple" = ,
    "ide" = ,
    "id" = .add_var_info(term),
    "vm" = .add_var_info(term),
    "ped" = .add_var_info(term),
    "fa" = .add_var_info(term),
    "fa_gxe" = {
      if (term$outer %in% names(data)) {
        f <- factor(data[[term$outer]])
        term$n_outer <- nlevels(f)
        # Data-aware identifiability preflight: a factor-analytic
        # covariance is identifiable only for k < n_outer. The data-free
        # k < 1 floor is enforced earlier in .classify_interaction(); the
        # upper bound needs n_outer, which is known only here.
        if (!is.null(term$k) && term$k >= term$n_outer) {
          stop(.fb_refusal_condition(
            reason_code = "fa_rank_exceeds_dim",
            message = paste0(
              "fa() requires k < the number of levels of \"",
              term$outer, "\"; got k = ", term$k, " with ",
              term$n_outer, " level(s). Use k <= ",
              term$n_outer - 1L, ", or us(", term$outer,
              ") for a full unstructured covariance."
            )
          ))
        }
      }
      if (term$inner %in% names(data)) {
        f <- factor(data[[term$inner]])
        term$n_inner <- nlevels(f)
      }
      term
    },
    "us_gxe" = {
      if (term$outer %in% names(data)) {
        term$n_outer <- nlevels(factor(data[[term$outer]]))
      }
      if (term$inner %in% names(data)) {
        term$n_inner <- nlevels(factor(data[[term$inner]]))
      }
      term
    },
    "at_simple" = {
      if (term$outer %in% names(data)) {
        term$n_outer <- nlevels(factor(data[[term$outer]]))
      }
      if (term$inner %in% names(data)) {
        term$n_inner <- nlevels(factor(data[[term$inner]]))
      }
      term
    },
    "at_units" = {
      if (term$var %in% names(data)) {
        term$n_var <- nlevels(factor(data[[term$var]]))
      }
      term
    },
    "ar1_spatial" = {
      if (!is.null(term$row_var) && term$row_var %in% names(data)) {
        term$n_row <- length(unique(data[[term$row_var]]))
      }
      if (!is.null(term$col_var) && term$col_var %in% names(data)) {
        term$n_col <- length(unique(data[[term$col_var]]))
      }
      term
    },
    "nested" = {
      if (term$outer %in% names(data)) {
        term$n_outer <- nlevels(factor(data[[term$outer]]))
      }
      if (term$inner %in% names(data)) {
        term$n_inner <- nlevels(factor(data[[term$inner]]))
      }
      term
    },
    "factor_numeric_interaction" = {
      # Term descriptor is fully populated by .classify_fixed_term();
      # .enrich() is a no-op here (the random-term enrichment branches
      # do not apply -- this is a fixed term).
      term
    },
    "smooth_mgcv" = {
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop(
          "`mgcv` is required for s() smooth terms. ",
          "Install with: install.packages(\"mgcv\")",
          call. = FALSE
        )
      }
      if (is.null(term$var) || !term$var %in% names(data)) {
        stop(.fb_refusal_condition(
          reason_code = "smooth_variable_not_in_data",
          message = paste0(
            "Smooth s() variable \"",
            term$var,
            "\" not found in data."
          )
        ))
      }
      # Build the design matrix at parse time. The original call is
      # evaluated against mgcv's namespace so user-supplied k / bs /
      # by arguments are honoured. absorb.cons = TRUE folds the sum-
      # to-zero constraint into the basis so we can treat the
      # coefficients as exchangeable normals downstream.
      sm_spec <- eval(term$call, envir = asNamespace("mgcv"))
      sm_obj <- mgcv::smoothCon(
        sm_spec,
        data = data,
        absorb.cons = TRUE,
        scale.penalty = TRUE
      )[[1]]
      term$X <- sm_obj$X
      term$k <- ncol(sm_obj$X)
      term$smooth_label <- sm_obj$label
      # Retain the constructed mgcv Smooth object on the IR.
      # predict.flexybayes(newdata = ...) re-evaluates the basis via
      # mgcv::Predict.matrix(<smooth_obj>, newdata) -- the only correct
      # way to reproduce the training-time knot placement / quantile
      # structure on new data. Without this slot the predict path falls
      # back to stats::model.matrix() which silently produces wrong
      # predictions near or beyond the training range.
      term$smooth_obj <- sm_obj
      # v0.4.0: if the s() carried a
      # `representation` spec, resolve it now (k and n are known) and
      # truncate the basis to its rank-K principal-component
      # approximation. The model then carries B_K = X V_K (n x K) and
      # K coefficients rather than the full n x k basis; V_K and the
      # realised Frobenius capture are recorded on `approx_spec` for
      # the emit path (codegen), the predict-side projection, and
      # validate_approximation().
      if (!is.null(term$approx_expr)) {
        spec <- eval(term$approx_expr)
        if (
          !is.list(spec) ||
            is.null(spec$scheme) ||
            !is.character(spec$scheme) ||
            length(spec$scheme) != 1L
        ) {
          stop(.fb_refusal_condition(
            reason_code = "approximation_spec_invalid",
            message = paste0(
              "s(",
              term$var,
              ", representation = ...): the ",
              "representation spec must be a list (or fb_approx() ",
              "object) carrying a single-string `scheme`; e.g. ",
              "representation = list(scheme = \"low_rank_smooth\", ",
              "rank = 5L)."
            ),
            family_class = "flexybayes_approximation_spec_invalid"
          ))
        }
        # Validate scheme against the locked registry (refuses
        # unknown schemes); only low_rank_smooth has a smooth-basis
        # emit path at v0.4.0.
        .lookup_approximation(spec$scheme)
        if (!identical(spec$scheme, "low_rank_smooth")) {
          stop(.fb_refusal_condition(
            reason_code = "approximation_no_smooth_path",
            message = paste0(
              "s(",
              term$var,
              ", representation = ...): scheme '",
              spec$scheme,
              "' is registered but has no smooth-basis ",
              "emit path at this release; the only smooth ",
              "approximation scheme with an emit path is ",
              "low_rank_smooth."
            ),
            family_class = "flexybayes_approximation_no_smooth_path"
          ))
        }
        tr <- .truncate_smooth_basis(term$X, spec$rank, var = term$var)
        term$X_K <- tr$B_K
        term$approx_spec <- list(
          scheme = spec$scheme,
          rank = tr$rank,
          k = tr$k,
          V_K = tr$V_K,
          singular_values = tr$singular_values,
          frobenius_capture = tr$frobenius_capture
        )
      }
      term
    },
    term
  )
}

# ---------------------------------------------------------------- #
# Contrast-detection guard for factor:continuous interactions       #
# ---------------------------------------------------------------- #
#
# v0.2.6 ships the treatment-coded indexed-slopes emit
# only. Other contrast schemes (`contr.helmert`, `contr.sdif`,
# `contr.sum`, ordered factors, user-supplied custom contrast
# matrices) defer to the v0.3 representation IR. Raise a structured
# refusal at ingest time so the user sees a clean message naming the
# contrast and the deferral target rather than a downstream emit
# failure.
.stop_if_unsupported_contrast <- function(factor_name, col) {
  # Ordered factors -- v0.3 representation IR will carry polynomial
  # contrasts as a first-class IR class; v0.2.6 refuses.
  if (is.ordered(col)) {
    cond <- structure(
      class = c(
        "flexybayes_contrast_unsupported",
        "flexybayes_unsupported_contrast",
        "error",
        "condition"
      ),
      list(
        message = paste0(
          "ordered factor \"",
          factor_name,
          "\" not supported in factor:continuous interaction ",
          "(treatment contrasts only). Polynomial-contrast support for ",
          "ordered factors is deferred to the representation IR."
        ),
        call = NULL,
        contrast = "ordered",
        factor_name = factor_name,
        deferral_target = "the representation IR"
      )
    )
    stop(cond)
  }
  # Custom / non-treatment contrasts -- read the contrasts attribute.
  # Base R stores either a character(1) function name (when set via
  # `options(contrasts = ...)` followed by `factor(...)`) or a
  # numeric L x (L-1) matrix (when set via `contrasts(f) <- ...`).
  # The named-function path is rare in practice; the matrix path is
  # the dominant case for users who explicitly opt into helmert /
  # sdif / sum. Both must refuse.
  ctr <- attr(col, "contrasts")
  if (!is.null(ctr)) {
    ctr_name <- if (is.character(ctr) && length(ctr) == 1L) {
      ctr
    } else if (is.matrix(ctr) || is.numeric(ctr)) {
      .identify_contrast_matrix(ctr, nlevels(factor(col)))
    } else {
      "<custom_matrix>"
    }
    if (!identical(ctr_name, "contr.treatment")) {
      cond <- structure(
        class = c(
          "flexybayes_contrast_unsupported",
          "flexybayes_unsupported_contrast",
          "error",
          "condition"
        ),
        list(
          message = paste0(
            "contrast \"",
            ctr_name,
            "\" on factor \"",
            factor_name,
            "\" not supported in factor:continuous interaction ",
            "(treatment contrasts only). General contrast ",
            "support is deferred to the representation IR."
          ),
          call = NULL,
          contrast = ctr_name,
          factor_name = factor_name,
          deferral_target = "the representation IR"
        )
      )
      stop(cond)
    }
  }
  # Session-level default-contrasts trap: if the user set
  # `options(contrasts = c("contr.helmert", "contr.poly"))` (or any
  # non-treatment default for unordered factors) and the factor
  # column has no explicit contrasts attribute, the contrast would
  # be picked up at model.matrix() time and would silently mismatch
  # the indexed-emit treatment-coded construction. Refuse at ingest
  # so the user sees the structured deferral message.
  opt <- getOption("contrasts")
  if (
    !is.null(opt) &&
      length(opt) >= 1L &&
      !identical(opt[[1L]], "contr.treatment")
  ) {
    cond <- structure(
      class = c(
        "flexybayes_contrast_unsupported",
        "flexybayes_unsupported_contrast",
        "error",
        "condition"
      ),
      list(
        message = paste0(
          "session-level default contrast \"",
          opt[[1L]],
          "\" applies to unordered factor \"",
          factor_name,
          "\" in factor:continuous interaction; only ",
          "the treatment-coded emit is shipped. General contrast support is ",
          "deferred to the representation IR. Workaround: ",
          "`contrasts(",
          factor_name,
          ") <- contr.treatment(",
          "nlevels(",
          factor_name,
          "))` on the input data."
        ),
        call = NULL,
        contrast = opt[[1L]],
        factor_name = factor_name,
        deferral_target = "the representation IR"
      )
    )
    stop(cond)
  }
  invisible(NULL)
}

# Heuristic: compare a user-set contrast matrix against the
# canonical contr.{helmert, sum, treatment, sdif, poly} matrices at
# the same level count. Returns the matching contrast function name
# (length-1 character) or "<custom_matrix>" when none matches. The
# comparison ignores dimnames since `contrasts(f) <- contr.*(L)`
# typically drops them.
.identify_contrast_matrix <- function(mat, L) {
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
  mat <- unname(mat)
  candidates <- list(
    contr.treatment = function(n) stats::contr.treatment(n),
    contr.helmert = function(n) stats::contr.helmert(n),
    contr.sum = function(n) stats::contr.sum(n),
    contr.poly = function(n) stats::contr.poly(n)
  )
  if (requireNamespace("MASS", quietly = TRUE)) {
    candidates[["contr.sdif"]] <- function(n) MASS::contr.sdif(n)
  }
  for (nm in names(candidates)) {
    ref <- tryCatch(
      unname(as.matrix(candidates[[nm]](L))),
      error = function(e) NULL
    )
    if (is.null(ref)) {
      next
    }
    if (
      identical(dim(ref), dim(mat)) &&
        isTRUE(all.equal(ref, mat, check.attributes = FALSE, tolerance = 1e-10))
    ) {
      return(nm)
    }
  }
  "<custom_matrix>"
}
