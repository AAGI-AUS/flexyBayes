# fb_penalty -- the shrinkage-prior grammar
#
# Contract: SPEC_bayesian_selection_v5_2026-09-12.md sections 3.2, 4.3 and
# 4.4, and SPEC_penalised_selection_v3_2026-09-11.md section 6 (the
# target-call contract), Unit A.
#
# The grammar is parsed symbolically, exactly as `fb_prior()` parses
# `normal(0, 2.5)`: `term()`, `ridge()` and the rest are not functions, they
# are the vocabulary of a mini-language, and exporting them would put five
# nouns into the package's namespace that nobody can call. Argument matching
# is delegated to base R via `match.call()` against a stub carrying the
# family's formals, so `ridge(1)` and `ridge(scale = 1)` are one
# specification, and a call base R cannot bind is refused by name rather than
# carried unnamed into the compiler.
#
# The two-sided corpus in tests/testthat/test-penalty-grammar.R mutates every
# argument of every constructor by omission, duplication, unknown naming,
# extra position, wrong type and call shape. That corpus is the reason this
# file validates rather than coerces: v3 section 6 exists because the prior
# DSL accepted eight malformed targets before it was written.


# --- constructor ---------------------------------------------------------- #

#' Declare a joint shrinkage prior over candidate fixed effects
#'
#' `fb_penalty()` writes the shrinkage half of a selection model: which
#' candidate terms carry which prior, at what resolution, and on what scale.
#' It is the second argument of [fb_select()], and it validates at
#' construction so a malformed specification is refused where it was written
#' rather than inside a Stan program.
#'
#' Each argument is a two-sided formula, `target ~ family(...)`, and the
#' constructors on both sides are parsed symbolically: they are the grammar's
#' vocabulary, not exported functions.
#'
#' @section Targets:
#'
#' Targets are resolved by the typed term resolver against the model matrix's
#' assignment map, never by string matching on coefficient names. A numeric
#' main effect is one term and one column; a factor is one term and all of
#' its contrast columns; an interaction is one term and every column it
#' generates.
#'
#' * `term(everything())` -- every candidate term
#' * `term(x1, x2)` -- named terms, as bare symbols
#' * `term(c("x1", "x2"))` -- named terms, as a character vector
#' * `term(matches("^SLGA_"))` -- a regular expression on the term labels
#'
#' `starts_with()` and the rest of the tidyselect vocabulary are not
#' supported: a hidden evaluation environment for one selector is not worth
#' the dependency, and `matches()` covers the use.
#'
#' Later specifications override earlier ones for the coordinates they name,
#' so `term(everything()) ~ ridge(scale = 1)` followed by
#' `term(x1) ~ laplace(scale = 0.5)` is the ordinary idiom. Two selectors of
#' equal breadth naming one coordinate leave no order and are refused.
#' `term(...) ~ none()` marks a candidate explicitly unpenalised.
#'
#' @section Families:
#'
#' | Family | Arguments | Hyperparameters |
#' |---|---|---|
#' | `ridge(scale)` | `scale` > 0, prior sd on the standardised scale | fixed |
#' | `laplace(scale)` | `scale` > 0 | fixed |
#' | `elastic_net(lambda, alpha)` | `lambda` > 0, `alpha` in (0, 1) | fixed |
#' | `regularized_horseshoe(p0, slab_scale, slab_df, df, df_global, scope)` | `p0` in \[1, p - 1\] | global and local scales learned |
#' | `r2d2(mean_R2, prec_R2, cons_D2, scope)` | the explained-variance grammar | learned |
#'
#' `horseshoe()` is an alias of `regularized_horseshoe()`. `lasso()` is an
#' alias of `laplace()`, and announces once per session that the posterior
#' mean of a Laplace prior is not a lasso estimate -- the draws contain no
#' exact zeros, so "selected" is always a named decision and never a property
#' of the posterior.
#'
#' `regularized_horseshoe()` and `r2d2()` are class-wide in brms: they apply
#' to every population-level coefficient except the intercept, and brms
#' refuses a coefficient-specific row alongside them. A model with protected
#' fixed effects beyond the intercept therefore refuses, and names both ways
#' forward: a family that targets blocks exactly, or `scope = "all_fixed"`,
#' under which the protected coefficients also sit under the global-local
#' prior and the certificate records the approximation.
#'
#' The group, exclusive and adaptive families of the earlier grammar
#' (`group_lasso()`, `exclusive_lasso()`, `sparse_group_lasso()`,
#' `group_exclusive()`, `adaptive()`) are withdrawn from this release and
#' refuse by name. Their providers are not shipped, and a family without a
#' provider is a promise rather than a method.
#'
#' @param ... One or more two-sided formulas, `target ~ family(...)`,
#'   resolved in the order written.
#' @param resolution `"term"` (the default) resolves targets to whole model
#'   terms; `"coefficient"` resolves them to single design columns. A
#'   coordinate-wise family aimed at one contrast column of a multi-column
#'   term is refused at `"coefficient"`, because the support would move under
#'   a relevel that changes no model.
#' @param scaling `"train_fold_standardize"` (the default) learns each
#'   candidate column's centre and scale inside the training rows of every
#'   fold and states the prior on the standardised column; `"none"` states it
#'   on the column as supplied.
#'
#' @returns An `fb_penalty` object: a list carrying `$specs` (the parsed
#'   target and family of each formula), `$resolution` and `$scaling`.
#'
#' @examples
#' fb_penalty(term(everything()) ~ ridge(scale = 1))
#'
#' fb_penalty(
#'   term(everything()) ~ ridge(scale = 1),
#'   term(x1) ~ laplace(scale = 0.5),
#'   term(matches("^z")) ~ none()
#' )
#'
#' fb_penalty(term(everything()) ~ regularized_horseshoe(p0 = 3))
#' @export
fb_penalty <- function(...,
                       resolution = c("term", "coefficient"),
                       scaling = c("train_fold_standardize", "none")) {
  resolution <- match.arg(resolution)
  scaling <- match.arg(scaling)

  args <- list(...)
  if (length(args) == 0L) {
    stop(.fb_refusal_condition(
      reason_code = "selection_penalty_spec_malformed",
      message = paste0(
        "fb_penalty(): at least one `target ~ family(...)` formula is ",
        "required. A penalty with no specification is `penalty = NULL`, ",
        "which fb_select() already reads as \"every candidate protected\"."
      )
    ))
  }

  specs <- vector("list", length(args))
  for (i in seq_along(args)) {
    a <- args[[i]]
    if (!inherits(a, "formula")) {
      stop(.fb_refusal_condition(
        reason_code = "selection_penalty_spec_malformed",
        message = paste0(
          "each argument must be a two-sided formula ",
          "`target ~ family(...)`. Argument ", i, " is: ",
          .fb_prior_deparse(a), "."
        )
      ))
    }
    if (length(a) != 3L) {
      stop(.fb_refusal_condition(
        reason_code = "selection_penalty_spec_malformed",
        message = paste0(
          "formula ", i, " must be two-sided (target on the left, family on ",
          "the right). Got: ", .fb_prior_deparse(a), "."
        )
      ))
    }
    parsed <- .fb_pen_parse_family(a[[3]], envir = environment(a))
    specs[[i]] <- list(
      target_label = .fb_prior_deparse(a[[2]]),
      selector = .fb_pen_parse_target(a[[2]]),
      family = parsed$family,
      args = parsed$args
    )
  }

  obj <- list(specs = specs, resolution = resolution, scaling = scaling)
  class(obj) <- c("fb_penalty", "list")
  obj
}


#' @export
print.fb_penalty <- function(x, ...) {
  cat(sprintf(
    "<fb_penalty> %d specification(s), resolution = \"%s\", scaling = \"%s\"\n",
    length(x$specs), x$resolution, x$scaling
  ))
  tbl <- .fb_pen_family_table()
  for (i in seq_along(x$specs)) {
    s <- x$specs[[i]]
    cat(sprintf(
      "  %d. %s ~ %s(%s)\n",
      i, s$target_label, s$family, .fb_pen_format_args(s$args)
    ))
    cat(sprintf(
      "     route %s, hyperparameters %s\n",
      tbl[[s$family]]$route, tbl[[s$family]]$hyperparameter_mode
    ))
  }
  cat("  columns resolve against the model matrix at fit time.\n")
  invisible(x)
}


is_fb_penalty <- function(x) inherits(x, "fb_penalty")


# --- target parsing -------------------------------------------------------- #

# .fb_pen_parse_target() --- the left-hand side of one specification.
#
# Four shapes and no others. Anything call-shaped that is not `everything()`,
# `matches()` or `c()` is refused rather than deparsed: deparsing a call
# head-first is what turned `b(c("x", "z"))` into three coefficient names in
# the prior DSL, and the phantom head survived all the way to an R internals
# error about a formula the user wrote.
.fb_pen_parse_target <- function(expr) {
  if (!is.call(expr) || !identical(as.character(expr[[1L]]), "term")) {
    .fb_pen_stop_target(
      expr,
      paste0(
        "the left-hand side must be a `term(...)` selector. The earlier ",
        "grammar's group() and as_group() are not part of this release."
      )
    )
  }

  args_list <- as.list(expr[-1L])
  if (length(args_list) == 0L) {
    .fb_pen_stop_target(expr, "term() names no candidate.")
  }
  nms <- names(args_list) %||% rep("", length(args_list))
  if (any(nzchar(nms))) {
    .fb_pen_stop_target(
      expr,
      paste0(
        "term() takes positional selectors only; it was given the named ",
        "argument(s) ",
        paste0("`", nms[nzchar(nms)], "`", collapse = ", "), "."
      )
    )
  }
  if (any(vapply(args_list, .fb_pen_is_empty_symbol, logical(1L)))) {
    stop(.fb_refusal_condition(
      reason_code = "grammar_argument_empty",
      message = paste0(
        "fb_penalty(): `", .fb_prior_deparse(expr), "` carries an empty ",
        "argument -- a doubled or trailing comma leaves R's empty symbol ",
        "where a selector belongs."
      )
    ))
  }

  first <- args_list[[1L]]
  if (is.call(first)) {
    if (length(args_list) > 1L) {
      .fb_pen_stop_target(
        expr,
        paste0(
          "a `", as.character(first[[1L]]), "()` selector stands alone; it ",
          "cannot be mixed with further arguments."
        )
      )
    }
    return(.fb_pen_parse_selector_call(first, expr))
  }

  # Bare symbols and scalar strings, one term each.
  names_out <- character(length(args_list))
  for (k in seq_along(args_list)) {
    a <- args_list[[k]]
    ok <- is.symbol(a) ||
      (is.character(a) && length(a) == 1L && !is.na(a) && nzchar(a))
    if (!ok) {
      .fb_pen_stop_target(
        expr,
        paste0(
          "argument ", k, " (`", .fb_prior_deparse(a), "`) is neither a bare ",
          "term name nor a single non-empty string."
        )
      )
    }
    names_out[k] <- as.character(a)
  }
  if (anyDuplicated(names_out) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "grammar_argument_duplicated",
      message = paste0(
        "fb_penalty(): `", .fb_prior_deparse(expr), "` names ",
        paste0(
          "`", unique(names_out[duplicated(names_out)]), "`",
          collapse = ", "
        ),
        " more than once. One coordinate carries one prior, so a repeated ",
        "term in one selector is a specification nobody can act on."
      )
    ))
  }
  list(kind = "names", names = names_out)
}


.fb_pen_parse_selector_call <- function(call_expr, expr) {
  # A call whose head is itself a call --- `pkg::helper(x)`, an inline
  # function --- deparses to more than one string, and every branch below
  # would then compare a length-3 vector to a name. That is the shape that
  # used to reach base R as "the condition has length > 1", an R internals
  # message about a formula the user wrote.
  if (!is.name(call_expr[[1L]])) {
    .fb_pen_stop_target(
      expr,
      paste0(
        "`", .fb_prior_deparse(call_expr), "` is not a selector: term() ",
        "takes everything(), matches(), bare term names or c(\"a\", \"b\"), ",
        "written plainly and not through another call."
      )
    )
  }
  fn <- as.character(call_expr[[1L]])
  inner <- as.list(call_expr[-1L])

  if (identical(fn, "everything")) {
    if (length(inner) > 0L) {
      .fb_pen_stop_target(expr, "everything() takes no argument.")
    }
    return(list(kind = "everything"))
  }

  if (identical(fn, "matches")) {
    if (length(inner) != 1L) {
      .fb_pen_stop_target(
        expr,
        "matches() takes exactly one regular expression."
      )
    }
    pattern <- tryCatch(
      eval(inner[[1L]], envir = baseenv()),
      error = function(e) NULL
    )
    ok <- is.character(pattern) && length(pattern) == 1L &&
      !is.na(pattern) && nzchar(pattern)
    if (!ok) {
      .fb_pen_stop_target(
        expr,
        paste0(
          "matches() takes a single non-empty regular expression written as ",
          "a string. Got: ", .fb_prior_deparse(inner[[1L]]), "."
        )
      )
    }
    return(list(kind = "matches", pattern = pattern))
  }

  if (identical(fn, "c")) {
    value <- tryCatch(
      eval(call_expr, envir = baseenv()),
      error = function(e) NULL
    )
    ok <- is.character(value) && length(value) >= 1L &&
      !anyNA(value) && all(nzchar(value))
    if (!ok) {
      .fb_pen_stop_target(
        expr,
        paste0(
          "c() inside term() must build a character vector of term names ",
          "with no missing or empty entry. Got: ",
          .fb_prior_deparse(call_expr), "."
        )
      )
    }
    if (anyDuplicated(value) > 0L) {
      stop(.fb_refusal_condition(
        reason_code = "grammar_argument_duplicated",
        message = paste0(
          "fb_penalty(): `", .fb_prior_deparse(expr), "` names ",
          paste0("`", unique(value[duplicated(value)]), "`", collapse = ", "),
          " more than once."
        )
      ))
    }
    return(list(kind = "names", names = as.character(value)))
  }

  .fb_pen_stop_target(
    expr,
    paste0(
      "`", fn, "()` is not a supported selector. term() takes ",
      "everything(), a regular expression via matches(), bare term names, ",
      "or c(\"a\", \"b\"). The tidyselect vocabulary -- starts_with(), ",
      "ends_with(), all_of() and the rest -- is deliberately not supported."
    )
  )
}


.fb_pen_stop_target <- function(expr, detail) {
  stop(.fb_refusal_condition(
    reason_code = "selection_penalty_target_unsupported",
    message = paste0(
      "`", .fb_prior_deparse(expr), "` does not name a selectable target: ",
      detail
    )
  ))
}


.fb_pen_is_empty_symbol <- function(x) {
  is.symbol(x) && !nzchar(as.character(x))
}


# --- family parsing -------------------------------------------------------- #

# .fb_pen_parse_family() --- the right-hand side of one specification.
#
# The order matters. Withdrawn families are named before unknown ones, so a
# caller who wrote `group_lasso()` is told it was removed and when it
# returns, rather than that it never existed. Aliases resolve before
# validation, so the argument contract is stated once per canonical family.
.fb_pen_parse_family <- function(expr, envir = baseenv()) {
  # The head has to be a plain name as well as a call. `pkg::ridge(1)`
  # deparses head-first to three strings, and every comparison below would
  # then be made against a length-3 vector rather than against a family name.
  if (!is.call(expr) || !is.name(expr[[1L]])) {
    stop(.fb_refusal_condition(
      reason_code = "selection_penalty_family_unknown",
      message = paste0(
        "the right-hand side must be a family call written plainly, for ",
        "example `ridge(scale = 1)`. Got: ", .fb_prior_deparse(expr), "."
      )
    ))
  }

  fn <- as.character(expr[[1L]])
  if (fn %in% .fb_pen_withdrawn_families()) {
    stop(.fb_refusal_condition(
      reason_code = "selection_penalty_family_withdrawn",
      message = paste0(
        "`", fn, "()` is withdrawn from the public grammar in this release. ",
        "Its provider is not shipped, and a family whose provider does not ",
        "exist is a promise rather than a method. The group, exclusive and ",
        "adaptive families return in a later phase, each with its own ",
        "provider-and-prior contract. Available now: ",
        paste(setdiff(.fb_pen_supported_families(), "none"), collapse = ", "),
        "."
      ),
      penalty_family = fn
    ))
  }

  alias <- .fb_pen_family_aliases()
  canonical <- if (fn %in% names(alias)) alias[[fn]] else fn
  if (!canonical %in% .fb_pen_supported_families()) {
    stop(.fb_refusal_condition(
      reason_code = "selection_penalty_family_unknown",
      message = paste0(
        "`", fn, "()` names no shrinkage family. Supported: ",
        paste(.fb_pen_supported_families(), collapse = ", "),
        "; the aliases ", paste(names(alias), collapse = ", "),
        " resolve to ", paste(unname(alias), collapse = ", "), "."
      ),
      penalty_family = fn
    ))
  }
  if (identical(fn, "lasso")) {
    .fb_pen_lasso_alias_note()
  }

  # `p0` is the horseshoe's calibration of the global scale. R2D2 states the
  # same intent on the explained-variance grammar and has no such quantity,
  # so a caller who wrote it has a prior belief the family cannot carry.
  # Refusing it before the generic argument matcher is deliberate: "r2d2()
  # does not have p0" is true and unhelpful, and the contract names this
  # case as a hyperparameter refusal rather than a spelling one.
  raw_names <- names(as.list(expr[-1L])) %||% character(0)
  if (identical(canonical, "r2d2") && "p0" %in% raw_names) {
    stop(.fb_refusal_condition(
      reason_code = "penalty_hyperparameter_invalid",
      message = paste0(
        "`p0` was given to r2d2(). R2D2 calibrates through the explained ",
        "variance -- mean_R2, prec_R2, cons_D2 -- and has no prior guess at ",
        "the number of non-zero coefficients. Write ",
        "regularized_horseshoe(p0 = ...) for that calibration, or state the ",
        "same belief as an R2 prior."
      )
    ))
  }

  entry <- .fb_pen_family_table()[[canonical]]
  expr <- .fb_pen_name_args(expr, canonical, entry)
  args <- as.list(expr[-1L])
  args <- lapply(args, function(a) {
    tryCatch(eval(a, envir = envir), error = function(e) a)
  })

  .fb_pen_check_args(canonical, args, expr, entry)
  args <- utils::modifyList(entry$defaults, args)
  args <- args[intersect(entry$params, names(args))]
  list(family = canonical, args = args)
}


# Bind a family call's arguments to their canonical names using base R's own
# matcher, after refusing the three shapes whose base-R message would not say
# which rule was broken.
.fb_pen_name_args <- function(expr, fn, entry) {
  params <- entry$params
  args_list <- as.list(expr[-1L])
  nms <- names(args_list) %||% rep("", length(args_list))

  if (any(vapply(args_list, .fb_pen_is_empty_symbol, logical(1L)))) {
    stop(.fb_refusal_condition(
      reason_code = "grammar_argument_empty",
      message = paste0(
        "fb_penalty(): `", .fb_prior_deparse(expr), "` carries an empty ",
        "argument -- a doubled or trailing comma leaves R's empty symbol ",
        "where a value belongs."
      )
    ))
  }

  named <- nms[nzchar(nms)]
  dup <- unique(named[duplicated(named)])
  if (length(dup) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "grammar_argument_duplicated",
      message = paste0(
        "fb_penalty(): `", .fb_prior_deparse(expr), "` binds ",
        paste0("`", dup, "`", collapse = ", "),
        " more than once. Supply each argument of ",
        .fb_pen_signature(fn, params), " at most once."
      )
    ))
  }

  unknown <- setdiff(named, params)
  if (length(unknown) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "grammar_named_argument_unknown",
      message = paste0(
        "fb_penalty(): `", .fb_prior_deparse(expr), "` carries the argument",
        if (length(unknown) > 1L) "s" else "", " ",
        paste0("`", unknown, "`", collapse = ", "),
        ", which ", fn, "() does not have. Supported: ",
        .fb_pen_signature(fn, params), "."
      )
    ))
  }

  if (sum(!nzchar(nms)) > length(setdiff(params, named))) {
    stop(.fb_refusal_condition(
      reason_code = "grammar_arity_too_many",
      message = paste0(
        "fb_penalty(): `", .fb_prior_deparse(expr), "` has more positional ",
        "arguments than ", .fb_pen_signature(fn, params), " reads."
      )
    ))
  }

  stub <- function() NULL
  formals(stub) <- stats::setNames(
    rep(list(quote(expr = )), length(params)), params
  )
  matched <- tryCatch(
    match.call(stub, expr),
    error = function(e) conditionMessage(e)
  )
  if (is.character(matched)) {
    stop(.fb_refusal_condition(
      reason_code = "grammar_named_argument_unknown",
      message = paste0(
        "fb_penalty(): `", .fb_prior_deparse(expr), "` cannot be bound to ",
        .fb_pen_signature(fn, params), " (", matched, ")."
      )
    ))
  }
  matched
}


# Required arguments, then the domain of each supplied one.
.fb_pen_check_args <- function(fn, args, expr, entry) {
  nms <- names(args) %||% rep("", length(args))
  missing_args <- setdiff(entry$required, nms)
  if (length(missing_args) > 0L) {
    stop(.fb_refusal_condition(
      reason_code = "penalty_hyperparameter_invalid",
      message = paste0(
        "`", .fb_prior_deparse(expr), "` omits ",
        paste0("`", missing_args, "`", collapse = ", "),
        ", which ", fn, "() requires. Supported: ",
        .fb_pen_signature(fn, entry$params), "."
      )
    ))
  }
  for (nm in names(args)) {
    .fb_pen_check_hyperparameter(
      value = args[[nm]], name = nm,
      domain = entry$domain[[nm]] %||% "real", fn = fn, expr = expr
    )
  }
  invisible(NULL)
}


.fb_pen_check_hyperparameter <- function(value, name, domain, fn, expr) {
  label <- paste0("`", name, "` in ", .fb_prior_deparse(expr))

  if (identical(domain, "scope")) {
    ok <- is.character(value) && length(value) == 1L && !is.na(value) &&
      value %in% c("candidates", "all_fixed")
    if (!ok) {
      stop(.fb_refusal_condition(
        reason_code = "penalty_hyperparameter_invalid",
        message = paste0(
          label, " must be \"candidates\" or \"all_fixed\". Got ",
          .fb_prior_deparse(value), "."
        )
      ))
    }
    return(invisible(NULL))
  }

  if (!is.numeric(value) || length(value) != 1L) {
    stop(.fb_refusal_condition(
      reason_code = "penalty_hyperparameter_invalid",
      message = paste0(
        label, " must be a single number. Got ",
        if (is.numeric(value)) {
          paste0("a numeric vector of length ", length(value))
        } else if (is.language(value)) {
          paste0("the unevaluated expression `", .fb_prior_deparse(value), "`")
        } else {
          paste0("an object of class ", paste(class(value), collapse = "/"))
        },
        "."
      )
    ))
  }
  if (!is.finite(value)) {
    stop(.fb_refusal_condition(
      reason_code = "penalty_hyperparameter_invalid",
      message = paste0(
        label, " must be a finite number. Got ", format(value), "."
      )
    ))
  }
  ok <- switch(
    domain,
    "positive" = value > 0,
    "nonnegative" = value >= 0,
    "unit_open" = value > 0 && value < 1,
    "real" = TRUE,
    TRUE
  )
  if (isTRUE(ok)) {
    return(invisible(NULL))
  }
  requirement <- switch(
    domain,
    "positive" = "be strictly positive",
    "nonnegative" = "be >= 0",
    "unit_open" = "lie strictly inside (0, 1)",
    "lie in its domain"
  )
  stop(.fb_refusal_condition(
    reason_code = "penalty_hyperparameter_invalid",
    message = paste0(
      label, " must ", requirement, ". Got ", format(value), "."
    )
  ))
}


.fb_pen_signature <- function(fn, params) {
  paste0(fn, "(", paste(params, collapse = ", "), ")")
}


# The lasso alias note. Once per session, because it is a statement about
# what the posterior mean of a Laplace prior is not, and a caller who has
# read it once does not need it on every fold.
.fb_pen_lasso_alias_note <- function() {
  if (.emit_state_get("penalty_lasso_alias_note")) {
    return(invisible(NULL))
  }
  .emit_state_set("penalty_lasso_alias_note", TRUE)
  message(
    "fb_penalty(): lasso() is an alias of laplace(). Draws from a Laplace ",
    "prior contain no exact zeros, so the posterior mean is not a lasso ",
    "estimate; a selected set is a named decision taken afterwards."
  )
  invisible(NULL)
}


# Reset the alias latch. Test-only: the note is once per session by design,
# and a suite that runs the alias twice needs to see it twice.
.fb_pen_reset_alias_note <- function() {
  .emit_state_set("penalty_lasso_alias_note", FALSE)
}
