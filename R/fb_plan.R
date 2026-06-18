# fb_plan() --- plan-only dispatch surface
#
# Returns the dispatch + preflight + representation + memory decision
# the routing layer would make, without firing the backend. The fifth
# verb in the v0.3.8 API spine (flexybayes, fb_brms, fb_greta,
# fb_plan, triangulate); validate_approximation lands at v0.4.0 as
# the sixth.
#
# Design notes:
#
#  - Brms-style formula entry (`y ~ x + s(z) + (1 | g)`); the asreml-
#    style entry (`flexybayes(fixed = y ~ x, random = ~g, plan = TRUE)`)
#    builds the IR via `fb_from_asreml()` and then short-circuits
#    through the same internal `.fb_plan_from_ir()` constructor so the
#    routing-policy decision is identical regardless of ingest path.
#
#  - Preflight is always run (memory ceiling check + per-term estimator),
#    bypassing the 1e5-row dispatcher threshold. The whole point of
#    fb_plan() is to surface the routing decision *for inspection*
#    before any compute is committed --- silently skipping preflight on
#    small data would defeat that contract.
#
#  - cov_validation_policy is the dial set by `flexyBayes.trust_pd`
#    (v0.3.8). "full_pd" when trust_pd is FALSE; "structural_only"
#    when TRUE; "n/a" when the model carries no known-matrix terms.
#    Reported on the plan and on backend_decision()$preflight_summary
#    after fit.
#
#  - The print form is a flight checklist (one line per decision the
#    dispatch layer would make, with a one-line justification). The
#    summary form dumps the full preflight + representation tables.
#    The as.data.frame() coercion is wide single-row for programmatic
#    consumers; stable column ordering by the constant
#    `.FB_PLAN_DF_COLS` below.

#' Plan a flexyBayes fit without firing the backend
#'
#' Returns the dispatch + preflight + representation + memory decision
#' the routing layer would make, without running MCMC or Laplace
#' approximation. Useful for verifying the backend chosen, memory
#' estimate, and any structural refusals before paying the fit cost.
#'
#' @param formula     a brms-style two-sided formula
#'   (`y ~ x + s(z) + (1 | g)`).
#' @param data        a data.frame.
#' @param backend     one of `"greta"`, `"inla"`, `"brms"`, `"auto"`
#'   (default `"auto"`).
#' @param priors      optional `fb_prior()` list; defaults to the v0.2
#'   uniform-on-SD default.
#' @param known_matrices  named list of structured-covariance matrices
#'   referenced by `vm()` or `ped()` terms.
#' @param family,link standard family/link arguments.
#' @param weights     optional observation weights.
#' @param aggregate   `"auto"` / `TRUE` / `FALSE` --- as on `flexybayes()`.
#' @param memory_ceiling_gb  optional override for the preflight memory
#'   ceiling (defaults to `flexyBayes.preflight_ceiling_gb` option, or
#'   `flexyBayes.preflight_ram_fraction` x available RAM).
#' @param predict_plan optional `list(newdata = ..., chunk_size = ...)`
#'   to compute a prediction-shape plan. Plan-only;
#'   does not fire `predict()`.
#' @param ...         currently unused; reserved for future plan inputs.
#'
#' @return an `<fb_plan>` classed list. See `print.fb_plan()` for the
#'   surface; `summary.fb_plan()` for the verbose dump;
#'   `as.data.frame.fb_plan()` for the programmatic-consumer shape.
#'
#' @export
fb_plan <- function(
  formula,
  data,
  backend = c("auto", "greta", "inla", "brms"),
  priors = NULL,
  known_matrices = list(),
  family = "gaussian",
  link = NULL,
  weights = NULL,
  aggregate = "auto",
  memory_ceiling_gb = NULL,
  predict_plan = NULL,
  ...
) {
  .check_approximate_scheme(backend)
  backend <- match.arg(backend)
  aggregate <- .normalise_aggregate(aggregate)

  the_call <- match.call()
  data_name <- deparse(substitute(data))

  fb <- fb_from_brms(
    formula = formula,
    data = data,
    family = family,
    link = link,
    prior = priors,
    weights = weights,
    known_matrices = known_matrices
  )

  .fb_plan_from_ir(
    fb = fb,
    data = data,
    backend = backend,
    known_matrices = known_matrices,
    aggregate = aggregate,
    memory_ceiling_gb = memory_ceiling_gb,
    predict_plan = predict_plan,
    the_call = the_call,
    data_name = data_name
  )
}


# ---------------------------------------------------------------- #
# Internal constructor: builds <fb_plan> from a built IR.            #
# ---------------------------------------------------------------- #
#
# Used by:
#   - fb_plan()                  (brms-style entry; IR built upstream)
#   - flexybayes(plan = TRUE)    (asreml-style; IR built upstream)
#   - fb_brms(plan = TRUE)       (brms-style; IR built upstream)
#
# Every caller hands in the IR + the dispatcher-relevant inputs. This
# function never builds an IR itself; it consumes one. That way the
# ingest-path choice (asreml vs brms) does not appear in the routing
# decision.
.fb_plan_from_ir <- function(
  fb,
  data,
  backend,
  known_matrices,
  aggregate,
  memory_ceiling_gb = NULL,
  predict_plan = NULL,
  the_call = NULL,
  data_name = NULL
) {
  # ---- preflight (always; bypass dispatcher threshold) ---------- #
  pf_dataset <- .fb_dataset(data)
  preflight <- .fb_preflight(
    fb_ir = fb,
    fb_dataset = pf_dataset,
    memory_ceiling_gb = memory_ceiling_gb,
    known_matrices = known_matrices
  )

  preflight_refused <- inherits(preflight$refusal, "fb_preflight_refusal")

  # ---- aggregation plan (read off preflight) -------------------- #
  agg_plan <- preflight$aggregation_plan
  agg <- list(
    eligible = isTRUE(agg_plan$eligible),
    reason_codes = if (is.null(agg_plan)) {
      character(0L)
    } else {
      (agg_plan$reason_codes %||% character(0L))
    },
    K = if (is.null(agg_plan)) {
      NA_integer_
    } else {
      (agg_plan$K %||% NA_integer_)
    },
    N = preflight$n_rows %||% NA_real_,
    ratio = if (
      is.null(agg_plan) ||
        is.null(agg_plan$K) ||
        is.na(agg_plan$K) ||
        agg_plan$K == 0L
    ) {
      NA_real_
    } else {
      preflight$n_rows / agg_plan$K
    }
  )

  # ---- gate (when backend in {inla, auto}) ---------------------- #
  inla_installed <- requireNamespace("INLA", quietly = TRUE)
  gretaR_activated <- isTRUE(getOption("flexyBayes.gretaR_activated", FALSE))

  gate_outcome <- NA_character_
  gate_checks <- NULL
  gate_primary_rule <- NA_character_

  if (backend %in% c("inla", "auto")) {
    gated <- tryCatch(lgm_gate(fb, preflight = preflight), error = function(e) {
      e
    })
    if (inherits(gated, "error")) {
      # Gate itself errored (e.g., a structural input it could not
      # handle); record the failure on the plan, route as refuse_structural.
      gate_outcome <- "refuse_structural"
      gate_checks <- list(list(
        rule_id = "lgm_gate_error",
        reason = conditionMessage(gated)
      ))
      gate_primary_rule <- "lgm_gate_error"
    } else if (is_lgm_refusal(gated)) {
      primary <- gated$failures[[1L]]
      is_memory <- identical(primary$rule_id, "memory_feasibility_inla")
      gate_outcome <- if (is_memory) {
        "refuse_memory"
      } else {
        "refuse_structural"
      }
      gate_checks <- gated$failures
      gate_primary_rule <- primary$rule_id %||% NA_character_
    } else {
      gate_outcome <- "accept"
      gate_checks <- gated$checks %||% NULL
    }
  }

  # ---- resolve routing ----------------------------------------- #
  preflight_outcome <- if (preflight_refused) {
    "refuse_memory"
  } else {
    "clear"
  }

  routing <- .resolve_routing(
    user_request = backend,
    gate_outcome = if (is.na(gate_outcome)) {
      NA_character_
    } else {
      gate_outcome
    },
    preflight_outcome = if (backend %in% c("inla", "auto")) {
      preflight_outcome
    } else {
      NA_character_
    },
    inla_installed = inla_installed,
    gretaR_activated = gretaR_activated
  )

  # ---- aggregation overrides chosen backend when eligible ------ #
  # The dispatcher routes aggregate = "auto" or TRUE silently to the
  # aggregated-INLA path when the plan declares eligibility. Mirror
  # that decision on the <fb_plan> so the plan matches the fit it
  # would produce.
  if (
    (isTRUE(aggregate) || identical(aggregate, "auto")) &&
      isTRUE(agg$eligible) &&
      !preflight_refused &&
      inla_installed &&
      backend != "brms"
  ) {
    chosen_backend <- "inla"
    chosen_path <- "aggregated_inla"
    chosen_reason <- "gaussian_aggregated_eligible"
  } else {
    chosen_backend <- routing$chosen_backend
    chosen_path <- routing$reason_code
    chosen_reason <- routing$reason_code
  }

  # ---- representation plan ------------------------------------- #
  representation_plan <- .representation_plan_from_preflight(preflight)
  representation_plan <- .annotate_block_diagonal_counts(
    representation_plan,
    fb,
    known_matrices
  )

  # ---- cov_validation_policy ----------------------------------- #
  cov_policy <- .cov_validation_policy(known_matrices = known_matrices, fb = fb)

  # ---- known-matrix metadata ------------------------------------ #
  km <- .fb_plan_known_matrix_summary(known_matrices)

  # ---- predict_plan shape (optional) --------------------------- #
  pred_plan <- if (is.null(predict_plan)) {
    NULL
  } else {
    .fb_plan_predict_shape(predict_plan, fb, data)
  }

  # ---- memory estimate (bytes) --------------------------------- #
  mem_bytes <- preflight$total_estimate_bytes
  if (is.null(mem_bytes) || is.na(mem_bytes)) {
    mem_bytes <- NA_real_
  }

  # ---- engine label (human-readable) --------------------------- #
  engine_label <- .engine_label_for(chosen_backend, chosen_path)

  # ---- representation label (Representation: line) ------------- #
  # Detect a low-rank (or any registered)
  # approximation scheme requested on a smooth so the plan's
  # Representation line carries the [APPROX:] badge.
  plan_approx <- .collect_approx(fb$random_terms)
  approx_scheme <- if (length(plan_approx) > 0L) {
    plan_approx[[1L]]$scheme
  } else {
    NULL
  }
  repr_label <- .representation_label_for(
    chosen_path,
    agg,
    representation_plan,
    approx_scheme
  )

  # ---- assemble ------------------------------------------------- #
  structure(
    list(
      formula = if (!is.null(the_call$formula)) {
        the_call$formula
      } else if (!is.null(the_call$fixed)) {
        the_call$fixed
      } else {
        NULL
      },
      data_name = data_name %||% NA_character_,
      backend_requested = backend,
      backend_chosen = chosen_backend,
      path = chosen_path,
      reason_code = chosen_reason,
      will_fit = !preflight_refused,
      preflight_refused = preflight_refused,
      preflight_refusal = preflight$refusal,
      routing_policy_version = .ROUTING_POLICY_VERSION,
      gate_outcome = gate_outcome,
      gate_primary_rule = gate_primary_rule,
      gate_checks = gate_checks,
      aggregation = agg,
      representation_plan = representation_plan,
      representation_label = repr_label,
      engine_label = engine_label,
      rejected_routes = routing$rejected_routes,
      known_matrix_summary = km,
      cov_validation_policy = cov_policy,
      predict_plan = pred_plan,
      memory_estimate_bytes = mem_bytes,
      preflight = preflight,
      call = the_call
    ),
    class = c("fb_plan", "list")
  )
}


# ---------------------------------------------------------------- #
# Helpers                                                            #
# ---------------------------------------------------------------- #

# .cov_validation_policy() --- maps the flexyBayes.trust_pd option +
# the known-matrices argument into one of three policy values shown
# on the plan and recorded on backend_decision()$preflight_summary.
.cov_validation_policy <- function(known_matrices, fb) {
  has_km <- length(known_matrices) > 0L ||
    any(vapply(
      fb$random_terms,
      function(t) isTRUE(t$type %in% c("vm", "ped")),
      logical(1L)
    ))
  if (!has_km) {
    return("n/a")
  }
  if (isTRUE(getOption("flexyBayes.trust_pd", FALSE))) {
    return("structural_only")
  }
  "full_pd"
}


# .fb_plan_known_matrix_summary() --- per-matrix shape + class summary
# for the plan surface. Used by both the print form and the
# as.data.frame() coercion.
.fb_plan_known_matrix_summary <- function(known_matrices) {
  if (length(known_matrices) == 0L) {
    return(list(
      present = FALSE,
      n_matrices = 0L,
      names = character(0L),
      classes = character(0L),
      dims = list()
    ))
  }
  classes <- vapply(known_matrices, function(m) class(m)[[1L]], character(1L))
  dims <- lapply(known_matrices, function(m) {
    if (is.null(dim(m))) {
      c(length(m), 1L)
    } else {
      dim(m)
    }
  })
  list(
    present = TRUE,
    n_matrices = length(known_matrices),
    names = names(known_matrices) %||% rep("", length(known_matrices)),
    classes = classes,
    dims = dims
  )
}


# .fb_plan_predict_shape() --- derived prediction shape (rows, chunk
# count, per-chunk memory ballpark) from a `predict_plan` arg. Reads
# the supplied newdata dimensions; does not call predict.flexybayes().
.fb_plan_predict_shape <- function(predict_plan, fb, data) {
  if (!is.list(predict_plan)) {
    stop(
      "fb_plan(): `predict_plan` must be a list with at least ",
      "`newdata`.",
      call. = FALSE
    )
  }
  nd <- predict_plan$newdata
  if (is.null(nd)) {
    stop("fb_plan(): `predict_plan$newdata` is required.", call. = FALSE)
  }
  if (!is.data.frame(nd)) {
    stop(
      "fb_plan(): `predict_plan$newdata` must be a data.frame; ",
      "got: ",
      paste(class(nd), collapse = "/"),
      call. = FALSE
    )
  }
  n_new <- nrow(nd)
  chunk_size <- predict_plan$chunk_size %||% n_new
  n_chunks <- if (chunk_size <= 0L || chunk_size >= n_new) {
    1L
  } else {
    ceiling(n_new / chunk_size)
  }
  list(
    n_newrows = n_new,
    chunk_size = as.integer(chunk_size),
    n_chunks = as.integer(n_chunks),
    output_file = predict_plan$output_file %||% NA_character_,
    format = predict_plan$format %||% "auto"
  )
}


# .engine_label_for() --- human-readable engine label for the
# Representation:/Engine: print surface. The label disambiguates
# inference engine + approximation regime in one phrase.
.engine_label_for <- function(backend, path) {
  if (identical(backend, "greta")) {
    return("greta MCMC")
  }
  if (identical(backend, "brms")) {
    return("brms / Stan HMC")
  }
  if (identical(backend, "inla")) {
    if (identical(path, "aggregated_inla")) {
      return("INLA Laplace (aggregated)")
    }
    return("INLA Laplace")
  }
  if (identical(backend, "gretaR")) {
    return("gretaR (R-native; dormant)")
  }
  backend
}


# .representation_label_for() --- the Representation: line companion
# to engine_label. "exact" / "aggregated_exact (N:K ratio)" /
# "exact (block-diagonal, K blocks)" / "n/a".
.representation_label_for <- function(
  path,
  agg,
  representation_plan = NULL,
  approx_scheme = NULL
) {
  # A plan whose formula
  # routes a smooth through an approximation scheme surfaces the same
  # [APPROX: <scheme>] badge the fitted object will carry.
  if (!is.null(approx_scheme)) {
    return(paste0("approximate [APPROX: ", approx_scheme, "]"))
  }
  if (identical(path, "aggregated_inla")) {
    if (!is.na(agg$ratio) && agg$ratio >= 2) {
      return(sprintf("aggregated_exact (compression %.0f:1)", agg$ratio))
    }
    return("aggregated_exact")
  }
  if (!is.null(representation_plan)) {
    block_entries <- Filter(
      function(rp) {
        identical(
          rp$representation_class %||% NA_character_,
          .representation_class("block_diagonal")
        )
      },
      representation_plan
    )
    if (length(block_entries)) {
      k_total <- sum(vapply(
        block_entries,
        function(rp) as.integer(rp$block_count %||% NA_integer_),
        integer(1L)
      ))
      if (!is.na(k_total) && k_total > 0L) {
        return(sprintf("exact (block-diagonal, %d blocks)", k_total))
      }
      return("exact (block-diagonal)")
    }
  }
  "exact"
}


# .annotate_block_diagonal_counts() --- enriches representation_plan
# entries whose class is "block_diagonal" with the block count read
# from known_matrices[[<symbol>]]. The preflight estimator runs from
# fb_dataset metadata only, so it does not have access to
# known_matrices --- we annotate post-hoc here. Idempotent on entries
# that are not block_diagonal or whose carrier is not in
# known_matrices (e.g., a parser-bypass fixture).
.annotate_block_diagonal_counts <- function(
  representation_plan,
  fb,
  known_matrices
) {
  if (is.null(representation_plan) || length(representation_plan) == 0L) {
    return(representation_plan)
  }
  block_terms <- Filter(
    function(t) {
      isTRUE(t$type %in% c("vm", "ped")) &&
        !is.null(t$cov_representation) &&
        identical(t$cov_representation$format, "blocks")
    },
    fb$random_terms
  )
  if (length(block_terms) == 0L) {
    return(representation_plan)
  }

  for (i in seq_along(representation_plan)) {
    entry <- representation_plan[[i]]
    if (
      !identical(
        entry$representation_class %||% NA_character_,
        .representation_class("block_diagonal")
      )
    ) {
      next
    }
    matched <- NULL
    for (term in block_terms) {
      candidates <- unique(c(term$label, term$var, term$deparse))
      if (
        any(vapply(
          candidates,
          function(c) {
            !is.null(c) &&
              (identical(c, entry$term_id) ||
                grepl(c, entry$term_id, fixed = TRUE))
          },
          logical(1L)
        ))
      ) {
        matched <- term
        break
      }
    }
    if (is.null(matched)) {
      next
    }
    sym <- matched$cov_representation$data
    val <- known_matrices[[sym]]
    if (is.list(val) && !is.data.frame(val)) {
      representation_plan[[i]]$block_count <- length(val)
    }
  }
  representation_plan
}


# Stable column ordering for as.data.frame.fb_plan() ----------------- #
.FB_PLAN_DF_COLS <- c(
  "backend_requested",
  "backend_chosen",
  "path",
  "reason_code",
  "will_fit",
  "preflight_refused",
  "routing_policy_version",
  "gate_outcome",
  "gate_primary_rule",
  "aggregation_eligible",
  "aggregation_K",
  "aggregation_N",
  "aggregation_ratio",
  "n_representation_terms",
  "n_rejected_routes",
  "cov_validation_policy",
  "known_matrix_count",
  "memory_estimate_mb",
  "predict_planned",
  "predict_chunks",
  "engine_label",
  "representation_label"
)


# ---------------------------------------------------------------- #
# S3 methods                                                         #
# ---------------------------------------------------------------- #

#' Print an `<fb_plan>` --- flight-checklist form
#' @param x   an `<fb_plan>` object.
#' @param ... unused.
#' @export
print.fb_plan <- function(x, ...) {
  cat("== flexyBayes plan ", strrep("=", 38), "\n", sep = "")
  cat(
    "  Will fit:                ",
    if (isTRUE(x$will_fit)) "yes" else "no  (preflight refused)",
    "\n",
    sep = ""
  )
  cat("  Backend requested:       ", x$backend_requested, "\n", sep = "")
  cat("  Backend chosen:          ", x$backend_chosen, "\n", sep = "")
  cat("  Path:                    ", x$path, "\n", sep = "")
  cat("  Routing policy version:  ", x$routing_policy_version, "\n", sep = "")

  if (!is.na(x$gate_outcome)) {
    cat(
      "  Gate outcome:            ",
      x$gate_outcome,
      if (!is.na(x$gate_primary_rule)) {
        paste0("  (", x$gate_primary_rule, ")")
      } else {
        ""
      },
      "\n",
      sep = ""
    )
  }

  cat(
    "  Aggregation:             ",
    if (isTRUE(x$aggregation$eligible)) {
      sprintf(
        "eligible (K = %s, N:K = %.0f:1)",
        format(x$aggregation$K, big.mark = " "),
        x$aggregation$ratio
      )
    } else if (length(x$aggregation$reason_codes) > 0L) {
      paste0(
        "not eligible (",
        paste(utils::head(x$aggregation$reason_codes, 2L), collapse = "; "),
        ")"
      )
    } else {
      "not in scope"
    },
    "\n",
    sep = ""
  )

  cat("  Cov validation policy:   ", x$cov_validation_policy, "\n", sep = "")

  if (isTRUE(x$known_matrix_summary$present)) {
    cat(
      "  Known matrices:          ",
      x$known_matrix_summary$n_matrices,
      " (",
      paste(x$known_matrix_summary$names, collapse = ", "),
      ")\n",
      sep = ""
    )
  }

  if (
    !is.null(x$representation_plan) &&
      length(x$representation_plan) > 0L
  ) {
    cat("  Representation plan:\n", sep = "")
    for (rp in x$representation_plan) {
      cls <- rp$representation_class %||% "unknown"
      if (
        identical(cls, .representation_class("block_diagonal")) &&
          !is.null(rp$block_count)
      ) {
        cls <- paste0(cls, " (", rp$block_count, " blocks)")
      }
      cat(
        "    ",
        format(rp$term_id, width = 24),
        " -> ",
        format(cls, width = 30),
        "  ",
        rp$justification %||% "",
        "\n",
        sep = ""
      )
    }
  }

  if (length(x$rejected_routes) > 0L) {
    cat("  Rejected routes:\n", sep = "")
    for (rr in x$rejected_routes) {
      cat(
        "    ",
        format(rr$backend, width = 8),
        " -> ",
        rr$reason,
        "\n",
        sep = ""
      )
    }
  }

  if (!is.null(x$predict_plan)) {
    cat(
      "  Predict plan:            ",
      sprintf(
        "n_newrows = %d; %d chunk(s) of %d",
        x$predict_plan$n_newrows,
        x$predict_plan$n_chunks,
        x$predict_plan$chunk_size
      ),
      "\n",
      sep = ""
    )
  }

  if (!is.na(x$memory_estimate_bytes)) {
    cat(
      "  Memory estimate:         ",
      sprintf("~ %.1f MB (preflight)", x$memory_estimate_bytes / 1024^2),
      "\n",
      sep = ""
    )
  }

  # v0.3.10: when the per-term
  # INLA memory estimator ran, surface the breakdown rows directly
  # under the headline memory line. Total + overhead make the
  # representation-class contributions auditable at plan time.
  mem <- x$preflight$memory_estimate
  if (
    inherits(mem, "fb_memory_estimate") &&
      !is.null(mem$breakdown) &&
      nrow(mem$breakdown) > 0L
  ) {
    cat(sprintf(
      "    INLA per-term total:  ~ %.1f MB (overhead %.1fx)\n",
      mem$total / 1024^2,
      mem$overhead_factor %||% 2
    ))
    for (i in seq_len(nrow(mem$breakdown))) {
      r <- mem$breakdown[i, ]
      cat(sprintf(
        "      - %-22s %-26s %8.1f MB\n",
        r$term_label,
        r$representation,
        r$bytes / 1024^2
      ))
    }
  }

  cat("  Representation:          ", x$representation_label, "\n", sep = "")
  cat("  Engine:                  ", x$engine_label, "\n", sep = "")
  cat(strrep("=", 58), "\n", sep = "")
  invisible(x)
}


#' Summarise an `<fb_plan>` --- verbose form
#' @param object an `<fb_plan>` object.
#' @param ...    unused.
#' @export
summary.fb_plan <- function(object, ...) {
  print(object)
  cat("\n-- Preflight detail ", strrep("-", 38), "\n", sep = "")
  if (is.null(object$preflight)) {
    cat("  (no preflight; below the dispatch threshold)\n")
  } else {
    pt <- object$preflight$per_term_estimate
    if (length(pt) > 0L) {
      cat(sprintf(
        "  %-28s  %12s  %s\n",
        "term",
        "bytes",
        "representation_class"
      ))
      for (label in names(pt)) {
        e <- pt[[label]]
        cat(sprintf(
          "  %-28s  %12s  %s\n",
          label,
          format(round(e$design_memory_bytes), big.mark = " "),
          e$representation_class %||% "unknown"
        ))
      }
    }
    cat(sprintf(
      "  Total estimate:     %s bytes (ceiling %s)\n",
      format(round(object$preflight$total_estimate_bytes), big.mark = " "),
      format(round(object$preflight$ceiling_bytes), big.mark = " ")
    ))
  }
  invisible(object)
}


#' Coerce an `<fb_plan>` to data.frame --- one row, stable columns
#'
#' Stable column ordering by the internal vector `.FB_PLAN_DF_COLS`;
#' adding new fields appends rather than reorders.
#'
#' @param x         an `<fb_plan>` object.
#' @param row.names unused.
#' @param optional  unused.
#' @param ...       unused.
#' @export
as.data.frame.fb_plan <- function(x, row.names = NULL, optional = FALSE, ...) {
  row <- list(
    backend_requested = x$backend_requested,
    backend_chosen = x$backend_chosen,
    path = x$path,
    reason_code = x$reason_code,
    will_fit = isTRUE(x$will_fit),
    preflight_refused = isTRUE(x$preflight_refused),
    routing_policy_version = x$routing_policy_version,
    gate_outcome = x$gate_outcome %||% NA_character_,
    gate_primary_rule = x$gate_primary_rule %||% NA_character_,
    aggregation_eligible = isTRUE(x$aggregation$eligible),
    aggregation_K = x$aggregation$K %||% NA_integer_,
    aggregation_N = x$aggregation$N %||% NA_real_,
    aggregation_ratio = x$aggregation$ratio %||% NA_real_,
    n_representation_terms = length(x$representation_plan),
    n_rejected_routes = length(x$rejected_routes),
    cov_validation_policy = x$cov_validation_policy,
    known_matrix_count = x$known_matrix_summary$n_matrices %||% 0L,
    memory_estimate_mb = if (is.na(x$memory_estimate_bytes)) {
      NA_real_
    } else {
      x$memory_estimate_bytes / 1024^2
    },
    predict_planned = !is.null(x$predict_plan),
    predict_chunks = x$predict_plan$n_chunks %||% NA_integer_,
    engine_label = x$engine_label,
    representation_label = x$representation_label
  )
  # Stable column ordering (defensive --- list above already follows
  # .FB_PLAN_DF_COLS, but readers may rely on this contract).
  row <- row[.FB_PLAN_DF_COLS]
  out <- as.data.frame(
    row,
    stringsAsFactors = FALSE,
    row.names = row.names %||% "1"
  )
  out
}
